from os.path import exists, abspath, dirname
import time
import logging
import glob
import warnings
import itertools
import networkx as nx
import multiprocessing as mp
from scipy.special import comb

# Open Babel and RDKit libraries
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem import MolFromSmiles

# Local modules
from luna.config.params import ProjectParams
from luna.mol.features import FeatureExtractor
from luna.mol.fingerprint import generate_fp_for_mols
from luna.mol.entry import Entry, MolFileEntry
from luna.mol.groups import AtomGroupPerceiver
from luna.interaction.contact import get_contacts_with
from luna.interaction.calc import InteractionCalculator
from luna.interaction.fp.shell import ShellGenerator
from luna.interaction.fp.type import IFPType
from luna.wrappers.base import MolWrapper
from luna.util import deprecated, LUNAWarning
from luna.util.default_values import *
from luna.util.exceptions import *
from luna.util.file import *
from luna.util.logging import (new_logging_file, load_default_logging_config,
                               VERBOSITY_LEVEL)
from luna.util.multiprocessing_logging import (start_mp_handler,
                                               MultiProcessingHandler)
from luna.util.jobs import ArgsGenerator, ParallelJobs

from luna.MyBio.PDB.PDBParser import PDBParser
from luna.MyBio.PDB.FTMapParser import FTMapParser
from luna.MyBio.util import download_pdb, save_to_file, get_entity_from_entry
from luna.version import __version__, has_version_compatibility

from sys import setrecursionlimit

# Set a recursion limit to avoid RecursionError with the library pickle.
setrecursionlimit(RECURSION_LIMIT)

logger = load_default_logging_config()

MAX_NPROCS = mp.cpu_count() - 1


class StructureCache:

    def __init__(self, compounds, atm_grps_mngr):
        self.compounds = compounds
        self.atm_grps_mngr = atm_grps_mngr

        self._compounds_name = set([c.full_name
                                    for c in self.compounds])

    def is_compound_cached(self, comp):
        return comp.full_name in self._compounds_name


class EntryResults:

    """Store entry results.

    Parameters
    ----------
    entry: :class:`~luna.mol.entry.Entry`
        An :class:`~luna.mol.entry.Entry` object that represents a molecule or
        an entire chain.
    atm_grps_mngr: :class:`~luna.mol.groups.AtomGroupsManager`
        An :class:`~luna.mol.groups.AtomGroupsManager` object that stores the
        perceived atoms and atom groups in the vicinity given by ``entry``.
    interactions_mngr: :class:`~luna.interaction.calc.InteractionsManager`
        An :class:`~luna.interaction.calc.InteractionsManager` object that
        interactions in the vicinity given by ``entry``.
    ifp: :class:`~luna.interaction.fp.fingerprint.Fingerprint`, optional
        An interaction fingerprint (IFP) generated for ``entry``.
    mfp: RDKit :class:`~rdkit.DataStructs.cDataStructs.ExplicitBitVect` or \
            :class:`~rdkit.DataStructs.cDataStructs.SparseBitVect`, optional
        A molecular fingerprint generated for ``entry``.

    Attributes
    ----------
    entry : :class:`~luna.mol.entry.Entry`
    atm_grps_mngr: :class:`~luna.mol.groups.AtomGroupsManager`
    interactions_mngr: :class:`~luna.interaction.calc.InteractionsManager`
    ifp: :class:`~luna.interaction.fp.fingerprint.Fingerprint`
    mfp: RDKit :class:`~rdkit.DataStructs.cDataStructs.ExplicitBitVect` or \
            :class:`~rdkit.DataStructs.cDataStructs.SparseBitVect`
    version : str
        The LUNA's version with which results were generated.
    """

    def __init__(self,
                 entry,
                 atm_grps_mngr,
                 interactions_mngr,
                 ifp=None,
                 mfp=None):

        self.entry = entry
        self.atm_grps_mngr = atm_grps_mngr
        self.interactions_mngr = interactions_mngr
        self.ifp = ifp
        self.mfp = mfp
        self.version = __version__

    def save(self, output_file, compressed=True):
        """Write the pickled representation of this object to the file
        ``output_file``.

        Parameters
        ----------
        output_file : str
            The output file where the pickled representation will be saved.
        compressed : bool, optional
            If True (the default), compress the pickled representation as a
            gzip file (.gz).

        Raises
        -------
        FileNotCreated
            If the file could not be created.
        """
        pickle_data(self, output_file, compressed)

    @staticmethod
    def load(input_file):
        """Read the pickled representation of an `EntryResults` object from
        the file ``input_file`` and return the reconstituted object hierarchy
        specified therein. ``input_file`` can be a gzip-compressed file.

        Raises
        -------
        PKLNotReadError
            If the file could not be loaded.
        """
        return unpickle_data(input_file)


class Project:

    """Define a LUNA project.

    .. note::
        This class is not intended to be used directly because :meth:`run` is
        not implemented by default. Instead, you should use a class that
        inherits from `Project` and implements :meth:`run`.
        An example is the class :class:`LocalProject` that implements a custom
        :meth:`run` that saves results as local files.

    Parameters
    ----------
    entries : iterable of :class:`~luna.mol.entry.Entry`
        Entries determine the target molecule to which interactions and other
        properties will be calculated. They can be ligands, chains, etc, and
        can be defined in a number of ways. Each entry has an associated PDB
        file that may contain macromolecules (protein, RNA, DNA) and other
        small molecules, water, and ions.
        Refer to :class:`~luna.mol.entry.Entry` for more information.
    working_path : str
        Where project results will be saved.
    pdb_path : str
        Path containing local PDB files or to where the PDB files will be
        downloaded. PDB filenames must match that defined for the entries.
        If not provided, the default PDB path will be used.
    overwrite_path : bool
        If True, allow LUNA to overwrite any existing directory, which may
        remove files from a previous project. The default value is False.
    add_h : bool
        Define if you need to add hydrogens or not. The default value is True.

        .. note::
            To be cautious, it does not add hydrogens to NMR-solved structures
            and ligands initialized from molecular files
            (:class:`~luna.mol.entry.MolFileEntry` objects) as they usually
            already contain hydrogens.
    ph : float
        Control the pH and how the hydrogens are going to be added.
        The default value is 7.4.

        .. note::
            To be cautious, it does not modify the protonation of molecular
            files defined by a :class:`~luna.mol.entry.MolFileEntry` object.
    amend_mol : bool
        If True (the default), try to fix atomic charges, valence, and bond
        types for small molecules and residues at PDB files. Only molecules
        at PDB files are validated because they do not contain charge, valence,
        and bond types, which may cause molecules to be incorrectly perceived.
        More information :ref:`here <Ligands in PDB files>`.

        .. note::
            Molecules from external files
            (:class:`~luna.mol.entry.MolFileEntry` objects) will not be
            modified.
    atom_prop_file : str
        A feature definition file (FDef) containing all information needed to
        define a set of chemical or pharmacophoric features.
        The default value is 'LUNA.fdef', which contains default LUNA features
        definition.
    inter_calc : :class:`~luna.interaction.calc.InteractionCalculator`
        Define which and how interactions are calculated.
    binding_mode_filter : :class:`~luna.interaction.filter.BindingModeFilter`
        Define how to filter interactions based on binding modes.
    calc_mfp : bool
        If True, generate ECFP4 fingerprints for each entry in ``entries``.
        The default value is False.
    mfp_output : str
        If ``calc_mfp`` is True, save ECFP4 fingerprints to file
        ``mfp_output``. If not provided, fingerprints are saved at
        <``working_path``>/results/fingerprints/mfp.csv.
    calc_ifp : bool
        If True (the default), generate LUNA interaction fingerprints (IFPs)
        for each entry in ``entries``.
    ifp_num_levels : int
        The maximum number of iterations for fingerprint generation.
        The default value is 2.
    ifp_radius_step : float
        The multiplier used to increase shell size at each iteration.
        At iteration 0, shell radius is 0 * ``radius_step``, at iteration 1,
        radius is 1 * ``radius_step``, etc. The default value is 5.73171.
    ifp_length : int
        The fingerprint length (total number of bits).
        The default value is 4096.
    ifp_count : bool
        If True (the default), create a count fingerprint
        (:class:`~luna.interaction.fp.fingerprint.CountFingerprint`).
        Otherwise, return a bit fingerprint
        (:class:`~luna.interaction.fp.fingerprint.Fingerprint`).
    ifp_diff_comp_classes :
        If True (the default), include differentiation between compound
        classes. That means structural information originated from
        :class:`~luna.mol.groups.AtomGroup` objects belonging to residues,
        nucleotides,  ligands, or water molecules will be considered different
        even if their structural information are the same.
        This is useful for example to differentiate protein-ligand interactions
        from residue-residue ones.
    ifp_type : :class:`~luna.interaction.fp.type.IFPType`
        The fingerprint type (EIFP, FIFP, or HIFP). The default value is EIFP.
    ifp_output : str
        If ``calc_ifp`` is True, save LUNA interaction fingerprints (IFPs) to
        file ``ifp_output``. If not provided, fingerprints are saved at
        <``working_path``>/results/fingerprints/ifp.csv.
    ifp_sim_matrix_output : str, optional
        If provided, compute Tanimoto similarity between interaction
        fingerprints (IFPs) and save the similarity matrix to
        ``ifp_sim_matrix_output``.
    out_pse : bool
        If True, depict interactions save them as Pymol sessions (PSE file).
        The default value is False. PSE files are saved at
        <``working_path``>/results/pse.
    append_mode : bool
        If True, skip entries from processing if a result for them already
        exists in ``working_path``. This can save processing time in case
        additional entries are to be added to an existing project.
    verbosity : int
        Verbosity level. The higher the verbosity level the more information
        is displayed. Valid values are:

            * 4: DEBUG messages;
            * 3: INFO messages (the default);
            * 2: WARNING messages;
            * 1: ERROR messages;
            * 0: CRITICAL messages.
    logging_enabled : bool
        If True (the default), enable the logging system.
    nproc : int
        The number of CPUs to use. The default value is the ``maximum number
        of CPUs - 1``. If ``nproc`` is smaller than 1 or greater than the
        maximum amount of available CPUs at your PC, then ``nproc`` is set to
        its default value. If you set it to None, LUNA will be run serially.


    Attributes
    ----------
    entries : iterable of :class:`~luna.mol.entry.Entry`
    working_path : str
    pdb_path : str
    overwrite_path : bool
    add_h : bool
    ph : float
    amend_mol : bool
    atom_prop_file : str
    inter_calc : :class:`~luna.interaction.calc.InteractionCalculator`
    binding_mode_filter : :class:`~luna.interaction.filter.BindingModeFilter`
    calc_mfp : bool
    mfp_output : str
    calc_ifp : bool
    ifp_num_levels : int
    ifp_radius_step : float
    ifp_length : int
    ifp_count : bool
    ifp_diff_comp_classes : bool
    ifp_type : :class:`~luna.interaction.fp.type.IFPType`
    ifp_output : str
    out_pse : bool
    out_ifp_sim_matrix : bool
    append_mode : bool
    logging_file : str
        The file to where logging messages are saved.
    version : str
        The LUNA's version with which results were generated.
    errors : list of tuple
        Any errors found during the processing of an entry.
        Each tuple contains the input and the exception raised during the
        execution of a task with that input.
    """

    def __init__(self,
                 entries,
                 working_path,
                 pdb_path=PDB_PATH,
                 overwrite_path=False,

                 add_h=True,
                 ph=7.4,
                 amend_mol=True,
                 atom_prop_file=ATOM_PROP_FILE,

                 inter_calc=None,
                 binding_mode_filter=None,

                 calc_mfp=False,
                 mfp_output=None,

                 calc_ifp=False,
                 ifp_num_levels=2,
                 ifp_radius_step=5.73171,
                 ifp_length=IFP_LENGTH,
                 ifp_count=True,
                 ifp_diff_comp_classes=True,
                 ifp_type=IFPType.EIFP,
                 ifp_output=None,
                 ifp_sim_matrix_output=None,

                 out_pse=False,
                 pse_path=None,

                 use_cache=False,
                 append_mode=False,
                 verbosity=3,
                 logging_enabled=True,
                 nproc=MAX_NPROCS):

        # Property required by self._log()
        self.logging_enabled = logging_enabled

        self._log("info", "LUNA version: %s." % __version__)

        if (inter_calc is not None
                and isinstance(inter_calc, InteractionCalculator) is False):
            msg = ".".join([InteractionCalculator.__module__,
                            InteractionCalculator.__name__])
            raise IllegalArgumentError("The informed interaction "
                                       "configuration must be an instance "
                                       "of %s." % msg)
        elif inter_calc is None:
            self._log("info", "No interaction calculator object was defined "
                      "and the default will be used instead.")

        if append_mode:
            self._log("warning", "Append mode set ON, entries with existing "
                      "results will be skipped from the entries processing.")

        self.entries = entries
        self.working_path = working_path
        self.pdb_path = pdb_path
        self.overwrite_path = overwrite_path
        self.atom_prop_file = atom_prop_file or ATOM_PROP_FILE
        self.ph = ph
        self.amend_mol = amend_mol
        self.add_h = add_h

        if inter_calc is None:
            inter_calc = InteractionCalculator(inter_config=INTERACTION_CONFIG)
        self.inter_calc = inter_calc

        self.binding_mode_filter = binding_mode_filter

        # Molecular fingerprint parameters.
        self.calc_mfp = calc_mfp
        self.mfp_output = mfp_output

        # Interaction fingerprint parameters.
        self.calc_ifp = calc_ifp
        self.ifp_num_levels = ifp_num_levels
        self.ifp_radius_step = ifp_radius_step
        self.ifp_length = ifp_length
        self.ifp_count = ifp_count
        self.ifp_diff_comp_classes = ifp_diff_comp_classes
        self.ifp_type = ifp_type
        self.ifp_output = ifp_output
        self.ifp_sim_matrix_output = ifp_sim_matrix_output

        # PSE parameters
        self.out_pse = out_pse
        self.pse_path = pse_path

        # General parameters.
        self.use_cache = use_cache
        self.append_mode = append_mode
        self.nproc = nproc

        self._loaded_logging_file = False
        self.verbosity = verbosity

        self.version = __version__

        self._paths = ["chunks", "configs", "logs", "pdbs",
                       "results/interactions", "results/fingerprints",
                       "results/pse", "results", "tmp"]
        self.errors = []

        self.cache = None

    def __call__(self):
        raise NotImplementedError("This class is not callable. Use a class "
                                  "that implements this method.")

    @property
    def project_file(self):
        """str: Where the pickled representation of the LUNA project is \
        saved."""
        return "%s/project_v%s.pkl.gz" % (self.working_path, __version__)

    @property
    def results(self):
        """iterable of `EntryResults`: LUNA results for each entry."""
        for entry in self.entries:
            results = self.get_entry_results(entry)
            if results:
                yield results

    @property
    def interactions_mngrs(self):
        """iterable of :class:`~luna.interaction.calc.InteractionsManager`: \
            An :class:`~luna.interaction.calc.InteractionsManager` object \
            for each entry."""
        for entry in self.entries:
            results = self.get_entry_results(entry)
            if results:
                yield results.interactions_mngr

    @property
    def atm_grps_mngrs(self):
        """iterable of :class:`~luna.mol.groups.AtomGroupsManager`: \
            An :class:`~luna.mol.groups.AtomGroupsManager` object for \
            each entry."""
        for entry in self.entries:
            results = self.get_entry_results(entry)
            if results:
                yield results.atm_grps_mngr

    @property
    def ifps(self):
        """iterable of :class:`~luna.interaction.fp.fingerprint.Fingerprint`: \
            An interaction fingerprint (IFP) for each entry."""
        for entry in self.entries:
            results = self.get_entry_results(entry)
            if results:
                yield entry, results.ifp

    @property
    def mfps(self):
        """iterable of \
            RDKit :class:`~rdkit.DataStructs.cDataStructs.ExplicitBitVect` \
            or :class:`~rdkit.DataStructs.cDataStructs.SparseBitVect`: \
                A molecular fingerprint for each entry."""
        for entry in self.entries:
            results = self.get_entry_results(entry)
            if results:
                yield entry, results.mfp

    @property
    def nproc(self):
        """int: The number of CPUs to use."""
        return self._nproc

    @nproc.setter
    def nproc(self, nproc):
        if nproc is not None:
            if not isinstance(nproc, int) or isinstance(nproc, bool):
                msg = ("The number of processes must be an integer value, but "
                       "a(n) %s was provided instead. Therefore, the number "
                       "of processes 'nproc' was set to its maximum accepted "
                       "capacity (%d - 1 = %d)."
                       % (nproc.__class__.__name__, mp.cpu_count(),
                          MAX_NPROCS))
                self._log("warning", msg)
                nproc = MAX_NPROCS

            elif nproc < 1:
                msg = ("It was trying to create an invalid number of "
                       "processes (%s). Therefore, the number of processes "
                       "'nproc' was set to its maximum accepted capacity "
                       "(%d - 1 = %d)." % (str(nproc), mp.cpu_count(),
                                           MAX_NPROCS))
                self._log("warning", msg)
                nproc = MAX_NPROCS

            elif nproc >= mp.cpu_count():
                self._log("warning", "It was trying to create %d processes, "
                          "which is equal to or greater than the maximum "
                          "amount of available CPUs (%d). Therefore, the "
                          "number of processes 'nproc' was set to %d "
                          "to leave at least one CPU free."
                          % (nproc, mp.cpu_count(), MAX_NPROCS))
                nproc = MAX_NPROCS
        else:
            self._log("warning", "The number of processes was set to '%s'. "
                      "Therefore, LUNA will run jobs sequentially." % nproc)

        self._nproc = nproc

    @property
    def logging_enabled(self):
        """bool: If the logging system is enable or not."""
        return self._logging_enabled

    @logging_enabled.setter
    def logging_enabled(self, is_enabled):
        if not is_enabled:
            warnings.warn("Logging mode was set OFF. No logging information "
                          "will be saved from now on.", category=LUNAWarning,
                          stacklevel=2)

            logger.disabled = True
        else:
            warnings.warn("Logging mode was set ON. Logging information will "
                          "be saved from now on.", category=LUNAWarning,
                          stacklevel=2)
            logger.disabled = False

        self._logging_enabled = is_enabled

    @property
    def verbosity(self):
        """int: Verbosity level."""
        return self._verbosity

    @verbosity.setter
    def verbosity(self, verbosity):
        if verbosity not in VERBOSITY_LEVEL:
            v_levels = ", ".join(["%d (%s)"
                                  % (k, logging.getLevelName(v))
                                  for k, v in sorted(VERBOSITY_LEVEL.items())])
            error_msg = ("The informed verbosity level '%s' is not valid. "
                         "The valid levels are: %s."
                         % (repr(verbosity), v_levels))
            raise IllegalArgumentError(error_msg)
        else:
            self._log("info", "Verbosity set to: %d (%s)."
                      % (verbosity,
                         logging.getLevelName(VERBOSITY_LEVEL[verbosity])))

        self._verbosity = VERBOSITY_LEVEL[verbosity]

        # If the logging file has already been loaded, it is necessary
        # to update the logging verbosity level.
        if self._loaded_logging_file:
            self._init_logging_file(self.logging_file)

    def get_entry_results(self, entry):
        """Get results for a given entry.

        Parameters
        ----------
        entry : :class:`~luna.mol.entry.Entry`
            An entry from ``entries``.

        Returns
        -------
         : `EntryResults`
        """
        if isinstance(entry, Entry):
            entry = entry.to_string()

        pkl_file = "%s/chunks/%s.pkl.gz" % (self.working_path, entry)
        try:
            return EntryResults.load(pkl_file)
        except Exception as e:
            self._log("exception", e)

    def _log(self, level, message):
        if self.logging_enabled:
            try:
                getattr(logger, level)(message)
            except Exception:
                raise

    def _log_preferences(self):
        self._log("debug", "New project initialized...")
        params = []
        for key in sorted(self.__dict__):
            if key == "entries":
                params.append("\t\t\t-- # %s = %d"
                              % (key, len(self.__dict__[key])))
            else:
                params.append("\t\t\t-- %s = %s"
                              % (key, str(self.__dict__[key])))
        self._log("debug", "Preferences:\n%s" % "\n".join(params))

    def _init_logging_file(self, logging_filename=None, use_mp_handler=True):
        if self.logging_enabled:
            if not logging_filename:
                logging_filename = new_unique_filename(TMP_FILES)

            try:
                new_logging_file(logging_filename,
                                 logging_level=self.verbosity)

                start_mp_handler()

                self._log("info",
                          "Logging file '%s' initialized successfully."
                          % logging_filename)

                # Print preferences at the new logging file.
                self._log_preferences()

                self._loaded_logging_file = True
            except Exception as e:
                self._log("exception", e)
                raise FileNotCreated("Logging file '%s' could not be created."
                                     % logging_filename)

    def _close_logging_file(self):
        try:
            for handler in logger.handlers:
                if isinstance(handler, MultiProcessingHandler):
                    if isinstance(handler.sub_handler, logging.FileHandler):
                        handler.close()
                        logger.removeHandler(handler)
        except Exception:
            pass

    def _prepare_project_path(self, subdirs=None):
        self._log("info",
                  "Preparing project directory '%s'." % self.working_path)

        if self.pdb_path is None or not is_directory_valid(self.pdb_path):
            new_pdb_path = "%s/pdbs/" % self.working_path
            self._log("warning", "The provided PDB path '%s' is not valid "
                      "or does not exist. Therefore, PDBs will be saved at "
                      "the working path: %s" % (self.pdb_path, new_pdb_path))
            self.pdb_path = new_pdb_path

        if subdirs is None:
            subdirs = self._paths

        # Create main project directory.
        create_directory(self.working_path, self.overwrite_path)
        # Create subdirectories.
        for path in subdirs:
            create_directory("%s/%s" % (self.working_path, path))

        self._log("info", "Project directory '%s' created successfully."
                  % self.working_path)

    def _remove_empty_paths(self):
        for path in self._paths:
            remove_directory("%s/%s" % (self.working_path, path),
                             only_empty_paths=True)

    def remove_duplicate_entries(self):
        """Search and remove duplicate entries from ``entries``."""
        entries = {}
        for entry in self.entries:
            if entry.to_string() not in entries:
                entries[entry.to_string()] = entry
            else:
                self._log("debug", "An entry with id '%s' already exists in "
                          "the list of entries, so the entry %s is a "
                          "duplicate and will be removed."
                          % (entry.to_string(), entry))

        self._log("info", "The remotion of duplicate entries was finished. "
                  "%d entrie(s) were removed."
                  % (len(self.entries) - len(entries)))
        self.entries = list(entries.values())

    def _validate_entry_format(self, entry):
        if not entry.is_valid():
            raise InvalidEntry("Entry '%s' does not match a LUNA's entry "
                               "format." % entry.to_string())

    def verify_pdb_files_existence(self):
        """Verify if a local PDB file exists for each entry in ``entries``.
            If it does not find a given PDB file, then LUNA will try to
            download it from RCSB."""
        all_pdb_ids = set()
        to_download = set()
        for entry in self.entries:
            pdb_file = "%s/%s.pdb" % (self.pdb_path, entry.pdb_id)
            if not exists(pdb_file):
                to_download.add(entry.pdb_id)

            all_pdb_ids.add(entry.pdb_id)

        logger.info("%d PDB file(s) found at '%s' from a total of %d PDB(s). "
                    "So, %d PDB(s) need to be downloaded."
                    % ((len(all_pdb_ids) - len(to_download)), self.pdb_path,
                        len(all_pdb_ids), len(to_download)))

        if to_download:
            args = [(pdb_id, self.pdb_path) for pdb_id in to_download]
            pj = ParallelJobs(self.nproc)
            job_results = pj.run_jobs(args=args, consumer_func=download_pdb,
                                      job_name="Download PDBs")
            errors = job_results.errors

            # Warn the users for any errors found during
            # the entries processing.
            if errors:
                self._log("warning", "Number of PDBs with errors: %d. "
                          "Check the log file to see the complete list of "
                          "PDBs that failed." % len(errors))
                self._log("debug", "PDBs that failed: %s."
                          % ", ".join([e[0][0] for e in errors]))

    def _decide_hydrogen_addition(self, pdb_header, entry):
        if self.add_h:
            if "structure_method" in pdb_header:
                method = pdb_header["structure_method"]
                # If the method is not a NMR type does not add hydrogen
                # as it usually already has hydrogens.
                if method.upper() in NMR_METHODS:
                    self._log("debug", "The structure related to the entry "
                              "'%s' was obtained by NMR, so it will "
                              "not add hydrogens to it." % entry.to_string())
                    return False
            return True
        return False

    def _parse_complex(self, entry):

        pdb_file = "%s/%s.pdb" % (self.pdb_path, entry.pdb_id)
        entry.pdb_file = pdb_file

        pdb_parser = entry.parser
        if pdb_parser is None:
            pdb_parser = PDBParser(PERMISSIVE=True, QUIET=True,
                                   FIX_EMPTY_CHAINS=True,
                                   FIX_ATOM_NAME_CONFLICT=True,
                                   FIX_OBABEL_FLAGS=False)

        if isinstance(pdb_parser, FTMapParser):
            only_compounds = [entry.get_biopython_key(full_id=True)]
            structure = pdb_parser.get_structure(entry.pdb_id,
                                                 pdb_file,
                                                 only_compounds=only_compounds)
            pdb_file = "%s/pdbs/%s.pdb" % (self.working_path,
                                           entry.to_string())
            save_to_file(structure, pdb_file)
            entry.pdb_file = pdb_file
        else:
            structure = pdb_parser.get_structure(entry.pdb_id, pdb_file)

        if isinstance(entry, MolFileEntry):
            structure = entry.get_biopython_structure(structure, pdb_parser)

        ligand = get_entity_from_entry(structure, entry)
        ligand.set_as_target(is_target=True)

        return pdb_parser, structure, ligand

    def _get_perceiver(self, add_h, cache=None):
        feats_factory_func = ChemicalFeatures.BuildFeatureFactory
        feature_factory = feats_factory_func(self.atom_prop_file)
        feature_extractor = FeatureExtractor(feature_factory)

        perceiver = AtomGroupPerceiver(feature_extractor, add_h=add_h,
                                       ph=self.ph, amend_mol=self.amend_mol,
                                       cache=cache,
                                       tmp_path="%s/tmp" % self.working_path)

        return perceiver

    def cache_protein_properties(self, entry):
        self._log("info",
                  "Cache memory activated. Residue information will be "
                  "stored and reutilized for all ligands.")

        pdb_parser, structure, ligand = self._parse_complex(entry)
        add_h = self._decide_hydrogen_addition(pdb_parser.get_header(), entry)

        inter_config = self.inter_calc.inter_config
        radius = inter_config.get("cache_cutoff",
                                  BOUNDARY_CONFIG["cache_cutoff"])
        nb_pairs = get_contacts_with(ligand,
                                     entity=structure[0],
                                     radius=radius,
                                     level='R')
        nb_compounds = set([p[1] for p in nb_pairs])
        mol_objs_dict = {}
        if isinstance(entry, MolFileEntry):
            mol_objs_dict[entry.get_biopython_key()] = entry.mol_obj

        perceiver = self._get_perceiver(add_h)
        atm_grps_mngr = \
            perceiver.perceive_atom_groups(nb_compounds,
                                           mol_objs_dict=mol_objs_dict)

        # Remove any edges involving the ligand.
        valid_edges = set()
        for edge in atm_grps_mngr.graph.edges:
            if any([atm.parent.is_hetatm() or atm.parent.is_metal()
                    for atm in edge]) is False:
                valid_edges.add(edge)
        atm_grps_mngr.graph = nx.Graph()
        atm_grps_mngr.graph.add_edges_from(valid_edges)

        # Remove hetatm and metals from the cache.
        invalid_atm_grps = [ag for ag in atm_grps_mngr
                            if (ag.has_hetatm()
                                or ag.has_metal())]
        atm_grps_mngr.remove_atm_grps(invalid_atm_grps)

        # Update the list of valid compounds used to generate the cache.
        invalid_comps = set([c for ag in invalid_atm_grps
                             for c in ag.compounds
                             if c.is_hetatm() or c.is_metal()])
        nb_compounds = nb_compounds - invalid_comps

        # Remove dummy chain if the ligand comes from an
        # external molecular file.
        if isinstance(entry, MolFileEntry):
            chain = ligand.parent
            chain.detach_child(ligand.id)

            if len(chain.child_list) == 0:
                chain.parent.detach_child(chain.id)

        self.cache = StructureCache(nb_compounds, atm_grps_mngr)

    def _perceive_chemical_groups(self, entry, entity, ligand,
                                  add_h=False, cache=None):
        self._log("debug", "Starting pharmacophore perception "
                  "for entry '%s'" % entry.to_string())

        inter_config = self.inter_calc.inter_config
        radius = inter_config.get("bsite_cutoff",
                                  BOUNDARY_CONFIG["bsite_cutoff"])
        nb_pairs = get_contacts_with(ligand,
                                     entity=entity,
                                     radius=radius,
                                     level='R')

        mol_objs_dict = {}
        if isinstance(entry, MolFileEntry):
            mol_objs_dict[entry.get_biopython_key()] = entry.mol_obj

        nb_compounds = set([x[1] for x in nb_pairs])

        perceiver = self._get_perceiver(add_h, cache)
        atm_grps_mngr = \
            perceiver.perceive_atom_groups(nb_compounds,
                                           mol_objs_dict=mol_objs_dict)

        self._log("debug", "Pharmacophore perception for entry '%s' has "
                  "finished." % entry.to_string())

        return atm_grps_mngr

    def _create_mfp(self, entry):
        if isinstance(entry, MolFileEntry):
            rdmol_lig = MolFromSmiles(MolWrapper(entry.mol_obj).to_smiles())
            rdmol_lig.SetProp("_Name", entry.mol_id)

            return generate_fp_for_mols([rdmol_lig], "morgan_fp")[0]["fp"]
        else:
            # TODO: implement support for other entries.
            self._log("warning", "Currently, it cannot generate molecular "
                      "fingerprints for instances of %s."
                      % entry.__class__.__name__)

    def _create_ifp(self, atm_grps_mngr):
        sg = ShellGenerator(self.ifp_num_levels, self.ifp_radius_step,
                            diff_comp_classes=self.ifp_diff_comp_classes,
                            ifp_type=self.ifp_type)
        sm = sg.create_shells(atm_grps_mngr)

        unique_shells = not self.ifp_count
        return sm.to_fingerprint(fold_to_length=self.ifp_length,
                                 unique_shells=unique_shells,
                                 count_fp=self.ifp_count)

    def _create_ifp_file(self):
        ifp_output = self.ifp_output or ("%s/results/fingerprints/ifp.csv"
                                         % self.working_path)
        with open(ifp_output, "w") as OUT:
            if self.ifp_count:
                OUT.write("ligand_id,on_bits,count\n")
            else:
                OUT.write("ligand_id,on_bits\n")

            for entry, ifp in self.ifps:
                if self.ifp_count:
                    fp_bits_str = "\t".join([str(idx)
                                             for idx in ifp.counts.keys()])
                    fp_count_str = "\t".join([str(count) for count
                                              in ifp.counts.values()])
                    OUT.write("%s,%s,%s\n" % (entry.to_string(), fp_bits_str,
                                              fp_count_str))
                else:
                    fp_bits_str = "\t".join([str(x) for x
                                             in ifp.get_on_bits()])
                    OUT.write("%s,%s\n" % (entry.to_string(), fp_bits_str))

    def _create_mfp_file(self):
        mfp_output = (self.mfp_output or "%s/results/fingerprints/mfp.csv"
                      % self.working_path)
        with open(mfp_output, "w") as OUT:
            OUT.write("ligand_id,on_bits\n")
            for entry, mfp in self.mfps:
                try:
                    bits = mfp.GetOnBits()
                except Exception:
                    try:
                        bits = mfp.GetNonzeroElements().keys()
                    except Exception:
                        error_msg = ("Fingerprint bits cannot be recovered "
                                     "for entry '%s'." % entry.to_string())
                        raise InvalidFingerprintType(error_msg)

                fp_str = "\t".join([str(x) for x in bits])
                OUT.write("%s,%s\n" % (entry.to_string(), fp_str))

    def _calc_similarity(self, res1, res2):
        return "%s,%s,%s" % (res1.entry.to_string(),
                             res2.entry.to_string(),
                             str(res1.ifp.calc_similarity(res2.ifp)))

    def _generate_similarity_matrix(self, output_file):
        nargs = int(comb(len(self.entries), 2))
        args = ArgsGenerator(itertools.combinations(self.results, 2), nargs)

        header = "entry1,entry2,similarity"
        pj = ParallelJobs(self.nproc)
        return pj.run_jobs(args=args, consumer_func=self._calc_similarity,
                           output_file=output_file, output_header=header,
                           job_name="Calculate similarities")

    def run(self):
        """Run LUNA. However, this method is not implemented by default.
        Instead, you should use a class that inherits from `Project` and
        implements :meth:`run`. An example is the class :class:`LocalProject`
        that implements a custom :meth:`run` that saves results as local files.
        """
        self()

    def save(self, output_file, compressed=True):
        """Write the pickled representation of this project to the file
        ``output_file``.

        Parameters
        ----------
        output_file : str
            The output file where the pickled representation will be saved.
        compressed : bool, optional
            If True (the default), compress the pickled representation as a
            gzip file (.gz).

        Raises
        -------
        FileNotCreated
            If the file could not be created.
        """
        pickle_data(self, output_file, compressed)

    @staticmethod
    def get_project_file(working_path):
        """Get the pickled project file at ``working_path``."""
        return "%s/project_v%s.pkl.gz" % (working_path, __version__)

    @staticmethod
    def load(pathname, verbosity=3, logging_enabled=True):
        """Read the pickled representation of a `Project` object from a file or
        project path and return the reconstituted object hierarchy specified
        therein. The ``pathname`` can be a gzip-compressed file.

        Parameters
        ----------
        pathname : str
            A file containing the pickled representation of a `Project` object
            or the project path (``working_path``) from where the pickled
            representation will be recovered.
        verbosity : int
            Verbosity level. The higher the verbosity level the more
            information is displayed. Valid values are:

                * 4: DEBUG messages;
                * 3: INFO messages (the default);
                * 2: WARNING messages;
                * 1: ERROR messages;
                * 0: CRITICAL messages.
        logging_enabled : bool
            If True (the default), enable the logging system.

        Raises
        -------
        CompatibilityError
            If the project version is not compatible with the current
            LUNA version.
        PKLNotReadError
            If the file could not be loaded.
        IllegalArgumentError
            If the provided pathname does not exist or is an invalid
            file/directory.
        """

        # Check if the provided input path is a valid file
        # or a directory containing saved projects.
        if is_file_valid(pathname):
            input_file = pathname
        elif is_directory_valid(pathname):
            project_files = glob.glob("%s/project_v*.pkl.gz" % pathname)
            if len(project_files) == 1:
                input_file = project_files[0]
            elif len(project_files) == 0:
                raise PKLNotReadError("In the provided working path '%s', "
                                      "there is no saved project." % pathname)
            else:
                raise PKLNotReadError("In the provided working path '%s', "
                                      "there are multiple saved projects. "
                                      "Please, specify which one you want "
                                      "to load." % pathname)
        else:
            raise IllegalArgumentError("The provided path '%s' does not exist "
                                       "or is an invalid file/directory."
                                       % pathname)

        if not logging_enabled:
            logger.disabled = True

        logger.info("Reloading project saved in '%s'.\n" % input_file)

        proj_obj = unpickle_data(input_file)

        if has_version_compatibility(proj_obj.version):
            proj_obj._loaded_logging_file = False
            proj_obj.verbosity = verbosity
            proj_obj.logging_enabled = logging_enabled

            # Update the working path if the project has been
            # moved to a different path.
            curr_working_path = dirname(abspath(input_file))

            if proj_obj.working_path != curr_working_path:
                proj_obj.working_path = curr_working_path
                proj_obj.logging_file = ("%s/logs/project.log"
                                         % proj_obj.working_path)

            proj_obj._log("info", "Project reloaded successfully.")
            return proj_obj
        else:
            raise CompatibilityError("The project loaded from '%s' has a "
                                     "version (%s) not compatible with the "
                                     "current %s's version (%s)."
                                     % (input_file, proj_obj.version,
                                        __package__.upper(), __version__))

    @classmethod
    def from_config_file(cls, config_file=None):
        if config_file is not None and not exists(config_file):
            raise OSError("File '%s' does not exist." % config_file)

        proj_params = ProjectParams(config_file, fill_defaults=True)

        return cls(**proj_params)

    def save_config_file(self, config_file=None):
        config_file = config_file or f"{self.working_path}/configs/project.cfg"

        params = ProjectParams.from_project_obj(self)
        params.save_config_file(config_file)


class LocalProject(Project):

    """Define a local LUNA project, i.e., results are saved locally and not to
    a database.

    This class inherits from `Project` and implements
        :meth:`~luna.projects.Project.run`.

    Examples
    --------

    In this minimum example, we will calculate protein-ligand interactions for
    dopamine D4 complexes.

    First, we should define the ligand entries and initialize a new
    :class:`~luna.interaction.calc.InteractionCalculator` object.

    >>> from luna.util.default_values import LUNA_PATH
    >>> from luna.interaction.calc import InteractionCalculator
    >>> entries = list(MolFileEntry.from_file(\
input_file=f"{LUNA_PATH}/tutorial/inputs/MolEntries.txt",
    ...                                       pdb_id="D4", \
mol_file=f"{LUNA_PATH}/tutorial/inputs/ligands.mol2"))
    >>> ic = InteractionCalculator(inter_filter=\
InteractionFilter.new_pli_filter())

    Finally, just create the new LUNA project with desired parameters and call
    :meth:`~luna.projects.Project.run`. Here, we opted to define the parameters
    first as a dict, and then we pass it as an argument to `LocalProject`.

    >>> from luna import LocalProject
    >>> opts = {}
    >>> opts["working_path"] = "%s/Results/Test3" % main_path
    >>> opts["pdb_path"] = f"{LUNA_PATH}/tutorial/inputs/"
    >>> opts["entries"] = entries
    >>> opts["inter_calc"] = ic
    >>> proj_obj = LocalProject(**opts)
    >>> proj_obj.run()

    """

    def __init__(self, entries, working_path, **kwargs):
        super().__init__(entries=entries, working_path=working_path, **kwargs)

    def _process_entry(self, entry):

        start = time.time()

        self._log("debug", "Starting entry processing: %s."
                  % entry.to_string())

        try:
            # Check if the entry is in the correct format.
            # It also accepts entries whose pdb_id is defined by the filename.
            if isinstance(entry, MolFileEntry) is False:
                self._validate_entry_format(entry)

            # Entry results will be saved here.
            pkl_file = "%s/chunks/%s.pkl.gz" % (self.working_path,
                                                entry.to_string())
            if self.append_mode and exists(pkl_file):
                self._log("debug", "Since append mode is set ON, it will "
                          "skip entry '%s' because a result for this entry "
                          " already exists in the working path."
                          % entry.to_string())
                return

            pdb_parser, structure, ligand = self._parse_complex(entry)
            add_h = self._decide_hydrogen_addition(pdb_parser.get_header(),
                                                   entry)

            #
            # Perceive pharmacophoric properties
            #
            atm_grps_mngr = self._perceive_chemical_groups(entry,
                                                           structure[0],
                                                           ligand, add_h,
                                                           self.cache)
            atm_grps_mngr.entry = entry

            #
            # Calculate interactions
            #
            calc_func = self.inter_calc.calc_interactions
            interactions_mngr = calc_func(atm_grps_mngr.atm_grps)
            interactions_mngr.entry = entry

            # Create hydrophobic islands.
            atm_grps_mngr.merge_hydrophobic_atoms(interactions_mngr)

            if self.binding_mode_filter is not None:
                bm_filter_func = interactions_mngr.filter_out_by_binding_mode
                bm_filter_func(self.binding_mode_filter)

            # Generate IFP (Interaction fingerprint)
            ifp = None
            if self.calc_ifp:
                ifp = self._create_ifp(atm_grps_mngr)

            # Generate MFP (Molecular fingerprint)
            mfp = None
            if self.calc_mfp:
                mfp = self._create_mfp(entry)

            # Saving entry results.
            entry_results = EntryResults(entry, atm_grps_mngr,
                                         interactions_mngr, ifp, mfp)
            entry_results.save(pkl_file)

            # Saving interactions to CSV file.
            csv_file = ("%s/results/interactions/%s.csv"
                        % (self.working_path, entry.to_string()))
            interactions_mngr.to_csv(csv_file)

            # Saving interactions into a Pymol session.
            if self.out_pse:
                from luna.interaction.view import InteractionViewer
                pse_path = (self.pse_path
                            or "%s/results/pse/" % self.working_path)
                pse_file = "%s/%s.pse" % (pse_path, entry.to_string())
                piv = InteractionViewer(add_directional_arrows=True)
                piv.new_session([(entry, interactions_mngr,
                                  entry.pdb_file)], pse_file)

            self._log("debug",
                      "Processing of entry '%s' finished successfully."
                      % entry.to_string())

        except Exception:
            self._log("debug",
                      "Processing of entry '%s' failed. "
                      "Check the logs for more information."
                      % entry.to_string())
            raise

        proc_time = time.time() - start
        self._log("debug",
                  "Processing of entry '%s' took %.2fs."
                  % (entry.to_string(), proc_time))

    def _process_ifps(self, entry):
        start = time.time()

        self._log("debug",
                  "Starting IFP processing for entry '%s'."
                  % entry.to_string())

        try:
            pkl_file = ("%s/chunks/%s.pkl.gz"
                        % (self.working_path, entry.to_string()))

            if exists(pkl_file):
                # Reload results.
                entry_results = EntryResults.load(pkl_file)
                atm_grps_mngr = entry_results.atm_grps_mngr

                # Generate a new IFP.
                ifp = self._create_ifp(atm_grps_mngr)

                # Substitute old IFP by the new version and save the project.
                entry_results.ifp = ifp
                entry_results.save(pkl_file)
            else:
                error_msg = ("The IFP of the entry '%s' could not be "
                             "generated because its pickled data file "
                             "'%s' was not found." % (entry.to_string(),
                                                      pkl_file))
                raise FileNotFoundError(error_msg)

        except Exception:
            self._log("debug", "IFP processing for entry '%s' failed. Check "
                      "the logs for more information." % entry.to_string())
            raise

        proc_time = time.time() - start
        self._log("debug", "IFP processing for entry '%s' took %.2fs." %
                  (entry.to_string(), proc_time))

    def _process_mfps(self, entry):
        start = time.time()

        self._log("debug",
                  "Starting MFP processing for entry '%s'."
                  % entry.to_string())

        try:
            pkl_file = ("%s/chunks/%s.pkl.gz"
                        % (self.working_path, entry.to_string()))

            if exists(pkl_file):
                # Reload results.
                entry_results = EntryResults.load(pkl_file)

                # Generate a new MFP.
                mfp = self._create_mfp(entry)

                # Substitute old MFP by the new version and save the project.
                entry_results.mfp = mfp
                entry_results.save(pkl_file)
            else:
                error_msg = ("The MFP of the entry '%s' could not be "
                             "generated because its pickled data file "
                             "'%s' was not found." % (entry.to_string(),
                                                      pkl_file))
                raise FileNotFoundError(error_msg)

        except Exception:
            self._log("debug", "MFP processing for entry '%s' failed. Check "
                      "the logs for more information." % entry.to_string())
            raise

        proc_time = time.time() - start
        self._log("debug", "MFP processing for entry '%s' took %.2fs." %
                  (entry.to_string(), proc_time))

    def __call__(self):

        if self.entries is None or len(self.entries) == 0:
            warnings.warn("There is nothing to be done as no "
                          "entry was informed.", category=LUNAWarning,
                          stacklevel=2)
            return

        if self.working_path is None:
            warnings.warn("A working path must be provided.",
                          category=LUNAWarning, stacklevel=2)
            return

        start = time.time()

        self._prepare_project_path()

        self.logging_file = "%s/logs/project.log" % self.working_path
        self._init_logging_file(self.logging_file)

        self.remove_duplicate_entries()

        self.save_config_file()

        self._log("info", "It will verify the existence of PDB files and "
                  "download them as necessary.")
        self.verify_pdb_files_existence()

        self._log("info", "Entries processing will start. Number of entries "
                  "to be processed is: %d." % len(self.entries))
        self._log("info", "The number of processes was set to: %s."
                  % str(self.nproc))

        if self.use_cache:
            entry = self.entries[0]
            self.cache_protein_properties(entry)

        # Run jobs either in Parallel or Sequentially (nproc = None).
        pj = ParallelJobs(self.nproc)
        job_results = pj.run_jobs(args=self.entries,
                                  consumer_func=self._process_entry,
                                  job_name="Entries processing")
        self.errors = job_results.errors

        # Remove failed entries.
        if self.errors:
            entries_with_error = set([e[0].to_string() for e in self.errors])
            self.entries = [e for e in self.entries
                            if e.to_string() not in entries_with_error]

        # If all molecules failed, it won't try to create fingerprints.
        if len(self.entries) == 0:
            self._log("critical", "Entries processing failed.")
        else:
            self._log("info", "Entries processing finished successfully.")

            # Warn the users for any errors found
            # during the entries processing.
            if self.errors:
                self._log("warning",
                          "Number of entries with errors: %d. Check the log "
                          "file to see the complete list of entries that "
                          "failed." % len(entries_with_error))
                self._log("debug", "Entries that failed: %s."
                          % ", ".join([e for e in entries_with_error]))

            # Generate IFP/MFP files
            if self.calc_ifp:
                self._create_ifp_file()
            if self.calc_mfp:
                self._create_mfp_file()

            if self.ifp_sim_matrix_output and len(self.entries) > 1:
                self._log("info", "Calculating the Tanimoto similarity "
                          "between fingerprints.")
                self._generate_similarity_matrix(self.ifp_sim_matrix_output)

        # Save the whole project information.
        self.save(self.project_file)

        # Remove unnecessary paths.
        self._remove_empty_paths()

        end = time.time()
        self._log("info", "Project creation completed!!!")
        self._log("info", "Total processing time: %.2fs." % (end - start))
        self._log("info", "Results were saved at %s." % self.working_path)
        self._log("info", "You can reload your project from %s.\n\n"
                  % self.project_file)

        # Properly close any filehandlers.
        self._close_logging_file()

    def generate_fps(self):
        """Generate LUNA interaction fingerprints (IFPs) or
        molecular fingerprints (MFPs).

        This function can be used to generate new FPs after a project is run.
        Thus, you can reload your project, vary IFP parameters
        (``ifp_num_levels``, ``ifp_radius_step``, ``ifp_length``,
        ``ifp_count``, ``ifp_diff_comp_classes``, ``ifp_type``,
        ``ifp_output``), and call `generate_fps` to create new IFPs without
        having to run the project from the scratch.

        Examples
        --------

        In the below example, we will assume a LUNA project object
        named ``proj_obj`` already exists.

        >>> from luna.interaction.fp.type import IFPType
        >>> proj_obj.ifp_num_levels = 5
        >>> proj_obj.ifp_radius_step = 1
        >>> proj_obj.ifp_length = 4096
        >>> proj_obj.ifp_type = IFPType.EIFP
        >>> proj_obj.ifp_output = "EIFP-4096__length-5__radius-1.csv"
        >>> proj_obj.generate_ifps()
        """

        if self.entries is None or len(self.entries) == 0:
            warnings.warn("There is nothing to be done as no "
                          "entry was informed.", category=LUNAWarning,
                          stacklevel=2)
            return

        if self.working_path is None:
            warnings.warn("A working path must be provided.",
                          category=LUNAWarning, stacklevel=2)
            return

        start = time.time()

        if self.calc_ifp is False and self.calc_mfp is False:
            warnings.warn("There is nothing to be done as both "
                          "'calc_ifp' and 'calc_mfp' is set to False.",
                          category=LUNAWarning, stacklevel=2)
            return

        self.overwrite_path = False

        if self.ifp_output is None or self.mfp_output is None:
            self._prepare_project_path(subdirs=["results",
                                                "results/fingerprints"])

        # Create a new directory for logs.
        if self.logging_enabled:
            if not exists("%s/logs/" % self.working_path):
                self._prepare_project_path(subdirs=["logs"])
            self._init_logging_file(self.logging_file)

        self._log("info", "Fingerprint generation will start. "
                  "Number of entries to be processed is: %d."
                  % len(self.entries))
        self._log("info", "The number of processes was set to: %s."
                  % str(self.nproc))

        def _create_fps(args, proc_func, file_func, job_name):
            # Run jobs either in Parallel or Sequentially (nproc = None).
            pj = ParallelJobs(self.nproc)
            job_results = pj.run_jobs(args=args,
                                      consumer_func=proc_func,
                                      job_name=job_name)
            errors = job_results.errors

            tmp_entries = self.entries
            # Identify failed entries.
            if errors:
                entries_with_error = set([e[0].to_string()
                                          for e in errors])
                tmp_entries = [e for e in self.entries
                               if e.to_string() not in entries_with_error]

            success = False
            # If all molecules failed, it won't try to create the output file.
            if len(tmp_entries) == 0:
                self._log("critical", f"{job_name} failed.")
            else:
                self._log("info", f"{job_name} finished successfully.")

                success = True

                # Warn the users for any errors found during
                # the entries processing.
                if errors:
                    self._log("warning", "Number of entries with errors: %d. "
                              "Check the log file to see the complete list of "
                              "entries that failed." % len(entries_with_error))
                    self._log("debug", "Entries that failed: %s."
                              % ", ".join([e for e in entries_with_error]))

                # Create an output file by calling the provided function.
                file_func()

            return success, errors

        all_errors = []
        if self.calc_ifp:
            success, errors = _create_fps(self.entries, self._process_ifps,
                                          self._create_ifp_file,
                                          "IFPs generation")
            all_errors.extend(errors)

            if success:
                if self.ifp_sim_matrix_output and len(self.entries) > 1:
                    self._log("info", "Calculating the Tanimoto similarity "
                              "between IFPs.")
                    output_file = self.ifp_sim_matrix_output
                    self._generate_similarity_matrix(output_file)

        if self.calc_mfp:
            success, errors = _create_fps(self.entries, self._process_mfps,
                                          self._create_mfp_file,
                                          "MFPs generation")
            all_errors.extend(errors)

        self.errors = all_errors

        # Remove unnecessary paths.
        self._remove_empty_paths()

        end = time.time()
        self._log("info", "Total processing time: %.2fs." % (end - start))
        self._log("info", "Results were saved at %s.\n\n" % self.working_path)

        # Properly close any filehandlers.
        self._close_logging_file()

    @deprecated("0.12.0", msg="Use the general `generate_fps` instead.")
    def generate_ifps(self):
        if self.entries is None or len(self.entries) == 0:
            warnings.warn("There is nothing to be done as no "
                          "entry was informed.", category=LUNAWarning,
                          stacklevel=2)
            return

        if self.working_path is None:
            warnings.warn("A working path must be provided.",
                          category=LUNAWarning, stacklevel=2)
            return

        start = time.time()

        self.calc_ifp = True
        self.overwrite_path = False

        if self.ifp_output is None:
            self._prepare_project_path(subdirs=["results",
                                                "results/fingerprints"])

        # Create a new directory for logs.
        if self.logging_enabled:
            if not exists("%s/logs/" % self.working_path):
                self._prepare_project_path(subdirs=["logs"])
            self._init_logging_file(self.logging_file)

        self._log("info", "Fingerprint generation will start. "
                  "Number of entries to be processed is: %d."
                  % len(self.entries))
        self._log("info", "The number of processes was set to: %s."
                  % str(self.nproc))

        # Run jobs either in Parallel or Sequentially (nproc = None).
        pj = ParallelJobs(self.nproc)
        job_results = pj.run_jobs(args=self.entries,
                                  consumer_func=self._process_ifps,
                                  job_name="Fingerprint generation")
        self.errors = job_results.errors

        tmp_entries = self.entries
        # Remove failed entries.
        if self.errors:
            entries_with_error = set([e[0].to_string() for e in self.errors])
            tmp_entries = [e for e in self.entries
                           if e.to_string() not in entries_with_error]

        # If all molecules failed, it won't try to create the output file.
        if len(tmp_entries) == 0:
            self._log("critical", "Fingerprint generation failed.")
        else:
            self._log("info", "Fingerprint generation finished successfully.")

            # Warn the users for any errors found during
            # the entries processing.
            if self.errors:
                self._log("warning", "Number of entries with errors: %d. "
                          "Check the log file to see the complete list of "
                          "entries that failed." % len(entries_with_error))
                self._log("debug", "Entries that failed: %s."
                          % ", ".join([e for e in entries_with_error]))

            # Generate IFP/MFP files
            if self.calc_ifp:
                self._create_ifp_file()
            if self.calc_mfp:
                self._create_mfp_file()

            if self.ifp_sim_matrix_output and len(self.entries) > 1:
                self._log("info", "Calculating the Tanimoto similarity "
                          "between fingerprints.")
                self._generate_similarity_matrix(self.ifp_sim_matrix_output)

        # Remove unnecessary paths.
        self._remove_empty_paths()

        end = time.time()
        self._log("info", "Total processing time: %.2fs." % (end - start))
        self._log("info", "Results were saved at %s.\n\n" % self.working_path)

        # Properly close any filehandlers.
        self._close_logging_file()
