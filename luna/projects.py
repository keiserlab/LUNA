from math import ceil
from os.path import exists
from collections import defaultdict
from functools import wraps
import time
import multiprocessing as mp
import logging

from pybel import readfile
from pybel import informats as OB_FORMATS
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem import MolFromPDBBlock, MolFromSmiles

# Local modules
from luna.mol.depiction import ligand_pharm_figure
from luna.mol.clustering import cluster_fps_butina
from luna.mol.features import FeatureExtractor
from luna.mol.fingerprint import generate_fp_for_mols
from luna.mol.entry import MolEntry
from luna.mol.groups import AtomGroupPerceiver
from luna.mol.interaction.contact import get_contacts_for_entity
from luna.mol.interaction.calc import InteractionCalculator
from luna.mol.interaction.conf import InteractionConf
from luna.mol.interaction.fp.shell import ShellGenerator
from luna.mol.wrappers.base import MolWrapper
from luna.mol.wrappers.rdkit import RDKIT_FORMATS, read_multimol_file
from luna.mol.amino_features import DEFAULT_AMINO_ATM_FEATURES
from luna.util.default_values import *
from luna.util.exceptions import *
from luna.util.file import pickle_data, unpickle_data, create_directory, get_file_format, get_unique_filename
from luna.util.logging import new_logging_file, load_default_logging_conf
from luna.util import iter_to_chunks
from luna.util.multiprocessing_logging import start_mp_handler

from luna.MyBio.PDB.PDBParser import PDBParser
from luna.MyBio.selector import ResidueSelector
from luna.MyBio.util import download_pdb, entity_to_string, get_entity_from_entry
from luna.version import __version__

logger = logging.getLogger()

PDB_PARSER = PDBParser(PERMISSIVE=True, QUIET=True, FIX_ATOM_NAME_CONFLICT=True, FIX_OBABEL_FLAGS=False)



class Project:

    def __init__(self,
                 entries,
                 working_path,
                 pdb_path=PDB_PATH,
                 has_local_files=False,
                 overwrite_path=False,
                 db_conf_file=None,
                 pdb_template=None,
                 atom_prop_file=ATOM_PROP_FILE,
                 try_h_addition=True,
                 ph=7.4,
                 amend_mol=True,
                 mol_obj_type='rdkit',

                 inter_conf=INTERACTION_CONF,
                 inter_calc=None,

                 mfp_opts=None,
                 mfp_output=None,

                 ifp_num_levels=7,
                 ifp_radius_step=1,
                 ifp_length=IFP_LENGTH,
                 ifp_count=False,
                 ifp_output=None,

                 similarity_func="BulkTanimotoSimilarity",
                 preload_mol_files=False,
                 default_properties=DEFAULT_AMINO_ATM_FEATURES,
                 butina_cutoff=0.2,
                 run_from_step=None,
                 run_until_step=None,
                 nproc=None):

        if mol_obj_type not in ACCEPTED_MOL_OBJ_TYPES:
            raise IllegalArgumentError("Objects of type '%s' are not currently accepted. "
                                       "The available options are: %s." % (mol_obj_type, ", ".join(ACCEPTED_MOL_OBJ_TYPES)))

        if inter_conf is not None and isinstance(inter_conf, InteractionConf) is False:
            raise IllegalArgumentError("The informed interaction configuration must be an instance of '%s'." % InteractionConf)

        if inter_calc is not None and isinstance(inter_calc, InteractionCalculator) is False:
            raise IllegalArgumentError("The informed interaction configuration must be an instance of '%s'." % InteractionCalculator)

        self.entries = entries
        self.working_path = working_path
        self.pdb_path = pdb_path
        self.has_local_files = has_local_files
        self.overwrite_path = overwrite_path
        self.db_conf_file = db_conf_file
        self.pdb_template = pdb_template
        self.atom_prop_file = atom_prop_file
        self.ph = ph
        self.amend_mol = amend_mol
        self.mol_obj_type = mol_obj_type
        self.try_h_addition = try_h_addition

        # Interaction calculator parameters.
        self.inter_conf = inter_conf

        if inter_calc is None:
            inter_calc = InteractionCalculator(inter_conf=self.inter_conf)
        self.inter_calc = inter_calc

        # Fingerprint parameters.
        self.mfp_opts = mfp_opts
        self.ifp_num_levels = ifp_num_levels
        self.ifp_radius_step = ifp_radius_step
        self.ifp_length = ifp_length
        self.ifp_count = ifp_count
        self.mfp_output = mfp_output
        self.ifp_output = ifp_output

        self.similarity_func = similarity_func
        self.butina_cutoff = butina_cutoff
        self.default_properties = default_properties
        self.run_from_step = run_from_step
        self.run_until_step = run_until_step
        self.preload_mol_files = preload_mol_files
        self.step_controls = {}

        if nproc is None:
            nproc = mp.cpu_count() - 1
        elif nproc < 1:
            logger.warning("It was trying to create an invalid number of processes (%s). Therefore, the number of "
                           "processes 'nproc' was set to its maximum accepted capability (%d)." % (str(nproc), (mp.cpu_count() - 1)))
            nproc = mp.cpu_count() - 1
        elif nproc >= mp.cpu_count():
            logger.warning("It was trying to create %d processes, which is equal to or greater than the maximum "
                           "amount of available CPUs (%d). Therefore, the number of processes 'nproc' was set to %d "
                           "to leave at least one CPU free." % (nproc, mp.cpu_count(), (mp.cpu_count() - 1)))
            nproc = mp.cpu_count() - 1
        self.nproc = nproc

        self.version = __version__

        self._references_updated = False

        load_default_logging_conf()

        self.log_preferences()

    def __call__(self):
        raise NotImplementedError("This class is not callable. Use a class that implements this method.")

    def run(self):
        self()

    def log_preferences(self):
        logger.info("New project initialized...")
        params = ["\t\t-- %s = %s" % (key, str(self.__dict__[key])) for key in sorted(self.__dict__)]
        logger.info("Preferences:\n%s" % "\n".join(params))

    def prepare_project_path(self):
        logger.warning("Initializing path '%s'..." % self.working_path)

        create_directory(self.working_path, self.overwrite_path)
        create_directory("%s/pdbs" % self.working_path)
        create_directory("%s/figures" % self.working_path)
        create_directory("%s/results" % self.working_path)
        create_directory("%s/logs" % self.working_path)
        create_directory("%s/tmp" % self.working_path)

        logger.warning("Path '%s' created successfully!!!" % self.working_path)

    def init_logging_file(self, logging_filename=None, use_mp_handler=True):
        if not logging_filename:
            logging_filename = get_unique_filename(TMP_FILES)

        try:
            new_logging_file(logging_filename)

            if use_mp_handler:
                start_mp_handler()

            logger.warning("Logging file initialized successfully.")

            # Print preferences at the new logging file.
            self.log_preferences()
        except Exception as e:
            logger.exception(e)
            raise FileNotCreated("Logging file could not be created.")

    def validate_entry_format(self, target_entry):
        if not target_entry.is_valid():
            raise InvalidEntry("Entry '%s' does not match a LUNA's entry format." % target_entry.to_string())

    def get_pdb_file(self, pdb_id):
        pdb_file = "%s/%s.pdb" % (self.pdb_path, pdb_id)

        if self.has_local_files:
            if not exists(pdb_file):
                raise FileNotFoundError("The PDB file '%s' was not found." % pdb_file)
        elif not exists(pdb_file):
            working_pdb_path = "%s/pdbs" % self.working_path
            pdb_file = "%s/%s.pdb" % (working_pdb_path, pdb_id)

            try:
                download_pdb(pdb_id=pdb_id, output_path=working_pdb_path)
            except Exception as e:
                logger.exception(e)
                raise FileNotCreated("PDB file '%s' was not created." % pdb_file) from e

        return pdb_file

    def decide_hydrogen_addition(self, pdb_header):
        if self.try_h_addition:
            if "structure_method" in pdb_header:
                method = pdb_header["structure_method"]
                # If the method is not a NMR type does not add hydrogen as it usually already has hydrogens.
                if method.upper() in NMR_METHODS:
                    logger.exception("The structure related to the entry '%s' was obtained by NMR, so it will "
                                     "not add hydrogens to it." % self.current_entry)
                    return False
            return True
        return False

    def perceive_chemical_groups(self, entity, ligand, add_h=False):
        feature_factory = ChemicalFeatures.BuildFeatureFactory(self.atom_prop_file)
        feature_extractor = FeatureExtractor(feature_factory)

        perceiver = AtomGroupPerceiver(feature_extractor, add_h=add_h, ph=self.ph, amend_mol=self.amend_mol,
                                       mol_obj_type=self.mol_obj_type, default_properties=self.default_properties,
                                       tmp_path="%s/tmp" % self.working_path)

        radius = self.inter_conf.boundary_cutoff or BOUNDARY_CONF.boundary_cutoff
        nb_compounds = get_contacts_for_entity(entity, ligand, level='R', radius=radius)

        mol_objs_dict = {}
        if isinstance(self.current_entry, MolEntry):
            mol_objs_dict[self.current_entry.get_biopython_key()] = self.current_entry.mol_obj

        atm_grps_mngr = perceiver.perceive_atom_groups(set([x[1] for x in nb_compounds]), mol_objs_dict=mol_objs_dict)

        logger.info("Chemical group perception finished!!!")

        return atm_grps_mngr

    def get_rdkit_mol(self, entity, target, mol_name="Mol0"):
        target_sel = ResidueSelector({target})
        pdb_block = entity_to_string(entity, target_sel, write_conects=False)
        rdmol = MolFromPDBBlock(pdb_block)
        rdmol.SetProp("_Name", mol_name)

        return rdmol

    def add_mol_obj_to_entries(self):
        mol_files = defaultdict(dict)
        for entry in self.entries:
            if not entry.is_mol_obj_loaded():
                entry.mol_obj_type = self.mol_obj_type
                mol_files[(entry.mol_file, entry.mol_file_ext)][entry.mol_id] = entry
            else:
                logger.info("Molecular object in entry %s was manually informed and will not be reloaded." % entry)

        tool = "Open Babel" if self.mol_obj_type == "openbabel" else "RDKit"
        logger.info("It will try to preload the molecular objects using %s. Total of files to be read: %d."
                    % (tool, len(mol_files)))

        try:
            for mol_file, mol_file_ext in mol_files:
                key = (mol_file, mol_file_ext)
                ext = mol_file_ext or get_file_format(mol_file)

                available_formats = OB_FORMATS if self.mol_obj_type == "openbabel" else RDKIT_FORMATS
                if ext not in available_formats:
                    raise IllegalArgumentError("Extension '%s' informed or assumed from the filename is not a format "
                                               "recognized by %s." % (ext, tool))

                if not exists(mol_file):
                    raise FileNotFoundError("The file '%s' was not found." % mol_file)

                logger.info("Reading the file '%s'. The number of target entries located in this file is %d."
                            % (mol_file, len(mol_files[key])))

                try:
                    if self.mol_obj_type == "openbabel":
                        for ob_mol in readfile(ext, mol_file):
                            mol_id = ob_mol.OBMol.GetTitle().strip()
                            if mol_id in mol_files[key]:
                                entry = mol_files[key][mol_id]
                                entry.mol_obj = ob_mol
                                del(mol_files[key][mol_id])
                                logger.info("A structure to the entry '%s' was found in the file '%s' and loaded "
                                            "into the entry." % (entry, mol_file))

                                # If there is no other molecules to search, just stop the loop.
                                if len(mol_files[key]) == 0:
                                    break
                        else:
                            logger.info("All target ligands located in the file '%s' were successfully loaded." % mol_file)
                    else:
                        targets = list(mol_files[key].keys())
                        for (i, rdk_mol) in enumerate(read_multimol_file(mol_file, ext, targets=targets, removeHs=False)):
                            mol_id = targets[i]
                            entry = mol_files[key][mol_id]
                            # It returns None if the molecule parsing generated errors.
                            if rdk_mol:
                                entry.mol_obj = rdk_mol
                                del(mol_files[key][mol_id])
                                logger.info("A structure to the entry '%s' was found in the file '%s' and loaded "
                                            "into the entry." % (entry, mol_file))
                            else:
                                logger.warning("The structure related to the entry '%s' was found in the file '%s', but it could "
                                               "not be loaded as errors were found while parsing it." % (entry, mol_file))
                except Exception as e:
                    logger.exception(e)
                    raise MoleculeObjectError("An error occurred while parsing the molecular file '%s' with %s and the molecule "
                                              "objects could not be created. Check the logs for more information." %
                                              (mol_file, tool))
        except Exception as e:
            logger.exception(e)
            raise

        invalid_entries = [e for m in mol_files for e in mol_files[m].values()]
        if invalid_entries:
            entries = set(self.entries)
            for entry in invalid_entries:
                entries.remove(entry)
                logger.warning("Entry '%s' was not found or generated errors, so it will be removed from the entries list."
                               % entry)
            logger.warning("%d invalid entries removed." % len(invalid_entries))
            self.entries = entries

    def generate_ligand_figure(self, rdmol, group_types):
        atm_types = defaultdict(set)
        atm_map = {}
        for atm in rdmol.GetAtoms():
            atm_map[atm.GetPDBResidueInfo().GetSerialNumber()] = atm.GetIdx()

        for grp in group_types.atm_grps:
            for atm in grp.atoms:
                atm_id = atm_map[atm.serial_number]
                atm_types[atm_id].update(set(grp.chemicalFeatures))

        output = "%s/figures/%s.svg" % (self.working_path,
                                        rdmol.GetProp("_Name"))
        ligand_pharm_figure(rdmol, atm_types, output, ATOM_TYPES_COLOR)

    def clusterize_ligands(self, fingerprints):
        fps_only = [x["fp"] for x in fingerprints]

        try:
            clusters = cluster_fps_butina(fps_only, cutoff=self.butina_cutoff, similarity_func=self.similarity_func)
        except Exception:
            raise ProcessingFailed("Clustering step failed.")

        lig_clusters = {}
        for i, cluster in enumerate(clusters):
            for mol_id in cluster:
                lig_clusters[fingerprints[mol_id]["mol"]] = i

        return lig_clusters

    def save(self, output_file, compressed=True):
        pickle_data(self, output_file, compressed)

    def update_references(self):
        result_pairs = []
        for entry1, atm_grps_mngr in self.neighborhoods:
            for entry2, inter_mngr in self.interactions:
                if entry1.to_string() == entry2.to_string():
                    # Updating references for entries.
                    entry1 = entry2
                    # Create a new pair of managers to update their atom group/interaction references.
                    result_pairs.append((atm_grps_mngr, inter_mngr))

        for atm_grps_mngr, inter_mngr in result_pairs:
            nb_mapping = {}
            for atm_grp in atm_grps_mngr.atm_grps:
                # First, reset the interactions list of an AtomGroup because the reference to
                # their interaction objects will be updated in the next loop.
                atm_grp.interactions = []
                nb_mapping[atm_grp] = atm_grp

            for inter in inter_mngr.interactions:
                # When the project is unpickled, it creates copies of AtomGroup objects and, therefore,
                # their references in the interactions point out to different objects, i.e., objects whose information
                # is equal but stored in different memory addresses. As a consequence, if an atom group is updated
                # in interactions, its corresponding atom group in AtomGroupsManager is not updated.
                #
                # So, let's update the atom groups in interactions with the atom groups in the AtomGroupsManager.
                # During this process, the reference to this interaction will be automatically updated in the new atom group.
                inter.src_grp = nb_mapping[inter.src_grp]
                inter.trgt_grp = nb_mapping[inter.trgt_grp]

        self._references_updated = True

    @staticmethod
    def load(input_file):
        logger.warning("Reloading project saved in '%s'" % input_file)
        proj_obj = unpickle_data(input_file)
        proj_obj.init_logging_file("%s/logs/project.log" % proj_obj.working_path)

        # Update the memory address of atom groups if it wasn't performed before.
        # The goal is to update circular references to atom groups, which may be lost after pickling a object,
        # especially after using multiprocessing.
        if proj_obj._references_updated is False:
            logger.warning("Reapplying interaction references to atom groups.")
            proj_obj.update_references()

        return proj_obj


class LocalProject(Project):

    def __init__(self, entries, working_path, has_local_files=False, **kwargs):
        super().__init__(entries=entries, working_path=working_path, has_local_files=has_local_files, **kwargs)

    def __call__(self):
        start = time.time()

        self.prepare_project_path()
        self.init_logging_file("%s/logs/project.log" % self.working_path)

        if self.preload_mol_files:
            self.add_mol_obj_to_entries()

        manager = mp.Manager()
        self.interactions = manager.list()
        self.neighborhoods = manager.list()

        self.mfps = []
        self.ifps = []

        chunk_size = ceil(len(self.entries) / self.nproc)
        chunks = iter_to_chunks(self.entries, chunk_size)

        processes = []
        for (i, l) in enumerate(chunks):
            p = mp.Process(name="Chunk %d" % i, target=self._process_entries, args=(l,))
            processes.append(p)
            p.start()

        for p in processes:
            p.join()

        self.interactions = list(self.interactions)
        self.neighborhoods = list(self.neighborhoods)

        # Update memory address for atom groups after performing a multiprocessing procedure.
        # The library multiprocessing stores data into an instance of mp.Manager() by using pickle,
        # and by doing so circular references are lost and must be updated.
        self.update_references()

        self._nb_mapping = {x[0].to_string(): x for x in self.neighborhoods}

        end = time.time()
        logger.info("Processing finished!!!")
        logger.info("Check the results at '%s'." % self.working_path)
        logger.info("Processing time: %.2fs." % (end - start))

    def _process_entries(self, entries):

        # Loop over each entry.
        for target_entry in entries:
            try:
                logger.info("Processing entry: %s." % target_entry)

                self.current_entry = target_entry

                # Check if the entry is in the correct format.
                # It also accepts entries whose pdb_id is defined by the filename.
                if isinstance(target_entry, MolEntry) is False:
                    self.validate_entry_format(target_entry)

                # TODO: allow the user to pass a pdb_file into entries.
                pdb_file = self.get_pdb_file(target_entry.pdb_id)
                target_entry.pdb_file = pdb_file

                structure = PDB_PARSER.get_structure(target_entry.pdb_id, pdb_file)
                add_hydrogen = self.decide_hydrogen_addition(PDB_PARSER.get_header())

                if isinstance(target_entry, MolEntry):
                    structure = target_entry.get_biopython_structure(structure, PDB_PARSER)

                ligand = get_entity_from_entry(structure, target_entry)
                ligand.set_as_target(is_target=True)

                atm_grps_mngr = self.perceive_chemical_groups(structure[0], ligand, add_hydrogen)

                #
                # Calculate interactions
                #
                interactions_mngr = self.inter_calc.calc_interactions(atm_grps_mngr.atm_grps)

                # Create hydrophobic islands.
                atm_grps_mngr.merge_hydrophobic_atoms(interactions_mngr)

                self.neighborhoods.append((target_entry, atm_grps_mngr))
                self.interactions.append((target_entry, interactions_mngr))

                logger.warning("Processing of entry %s finished successfully." % target_entry)

            except Exception as e:
                logger.exception(e)
                logger.warning("Processing of entry %s failed. Check the logs for more information." % target_entry)

    def generate_fps(self, calc_ifp=True, calc_mfp=False, save_files=False):

        chunk_size = ceil(len(self.entries) / self.nproc)
        chunks = iter_to_chunks(self.entries, chunk_size)

        manager = mp.Manager()
        self.ifps = manager.list()
        self.mfps = manager.list()

        processes = []
        for (i, l) in enumerate(chunks):
            p = mp.Process(name="Chunk %d" % i, target=self._process_fps, args=(l, calc_mfp, calc_ifp))
            processes.append(p)
            p.start()

        for p in processes:
            p.join()

        self.ifps = list(self.ifps)
        self.mfps = list(self.mfps)

        # Save files.
        if save_files:
            if calc_mfp:
                self.mfp_output = self.mfp_output or "%s/results/mfp.csv" % self.working_path
                with open(self.mfp_output, "w") as OUT:
                    OUT.write("ligand_id,smiles,on_bits\n")
                    for target_entry, mfp in self.mfps:
                        fp_str = "\t".join([str(x) for x in mfp.GetOnBits()])
                        OUT.write("%s,%s,%s\n" % (target_entry.to_string(), "", fp_str))

            if calc_ifp:
                self.ifp_output = self.ifp_output or "%s/results/ifp.csv" % self.working_path
                with open(self.ifp_output, "w") as OUT:
                    if self.ifp_count:
                        OUT.write("ligand_id,smiles,on_bits,count\n")
                    else:
                        OUT.write("ligand_id,smiles,on_bits\n")

                    for target_entry, ifp in self.ifps:
                        if self.ifp_count:
                            fp_bits_str = "\t".join([str(idx) for idx in ifp.counts.keys()])
                            fp_count_str = "\t".join([str(count) for count in ifp.counts.values()])
                            OUT.write("%s,%s,%s,%s\n" % (target_entry.to_string(), "", fp_bits_str, fp_count_str))
                        else:
                            fp_bits_str = "\t".join([str(x) for x in ifp.get_on_bits()])
                            OUT.write("%s,%s,%s\n" % (target_entry.to_string(), "", fp_bits_str))

    def _process_fps(self, entries, calc_mfp, calc_ifp):

        # Loop over each entry.
        for target_entry in entries:
            try:
                logger.info("Generating fingerprints for entry: %s." % target_entry)

                if calc_ifp:
                    # Recover the neighborhood information for an entry.
                    atm_grps_mngr = self._nb_mapping[target_entry.to_string()][1]

                    shells = ShellGenerator(self.ifp_num_levels, self.ifp_radius_step)
                    sm = shells.create_shells(atm_grps_mngr)

                    unique_shells = not self.ifp_count
                    ifp = sm.to_fingerprint(fold_to_size=self.ifp_length, unique_shells=unique_shells, count_fp=self.ifp_count)
                    data = (target_entry, ifp)
                    self.ifps.append(data)

                if calc_mfp:
                    if isinstance(target_entry, MolEntry):
                        rdmol_lig = MolFromSmiles(MolWrapper(target_entry.mol_obj).to_smiles())
                        rdmol_lig.SetProp("_Name", target_entry.mol_id)

                        data = (target_entry, generate_fp_for_mols([rdmol_lig], "morgan_fp")[0]["fp"])
                        self.mfps.append(data)
                    else:
                        logger.warning("Currently, it cannot generate molecular fingerprints for "
                                       "instances of %s." % target_entry.__class__.__name__)
            except Exception as e:
                logger.exception(e)
                logger.warning("Generation of fingerprints for entry %s failed. Check the logs for more information." % target_entry)
