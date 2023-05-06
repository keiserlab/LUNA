import tempfile
import logging
from io import StringIO
from os.path import exists
from shutil import move as rename_pdb_file


from openbabel.pybel import readfile
from openbabel.pybel import Molecule as PybelWrapper

from rdkit.Chem import (MolFromMolBlock, SanitizeFlags,
                        SanitizeMol, MolToMolFile)

from luna.util.file import (is_directory_valid, new_unique_filename,
                            remove_files)
from luna.util.default_values import ENTRY_SEPARATOR, OPENBABEL
from luna.MyBio.selector import AtomSelector
from luna.MyBio.PDB.PDBList import PDBList
from luna.MyBio.PDB.PDBIO import Select
from luna.MyBio.PDB import Selection
from luna.MyBio.PDB.PDBIO import PDBIO
from luna.MyBio.PDB.PDBParser import PDBParser
from luna.mol.validator import MolValidator
from luna.mol.standardiser import Standardizer
from luna.wrappers.base import MolWrapper
from luna.wrappers.obabel import convert_molecule
from luna.wrappers.rdkit import read_mol_from_file

from luna.util.exceptions import (IllegalArgumentError, MoleculeNotFoundError,
                                  ChainNotFoundError, FileNotCreated,
                                  PDBNotReadError, MoleculeSizeError,
                                  MoleculeObjectError)


logger = logging.getLogger()


ENTITY_LEVEL_NAME = {"A": "Atom",
                     "R": "Residue",
                     "C": "Chain",
                     "M": "Model",
                     "S": "Structure"}


def download_pdb(pdb_id, output_path=".", overwrite=False):
    """Download a PDB file from RCSB.org.

    Parameters
    ----------
    pdb_id : str
        4-symbols structure Id from PDB (e.g. 3J92).
    output_path : str
        Put the PDB file in this directory.
    overwrite : bool
        If True, overwrite any existing PDB files.
    """
    logger.debug("It will try to download the PDB '%s' and store it at the "
                 "directory '%s'." % (pdb_id, output_path))

    if pdb_id is not None and pdb_id.strip() != "":
        if is_directory_valid(output_path):
            new_pdb_file = "%s/%s.pdb" % (output_path, pdb_id)
            # If the file already exists and is not to overwrite
            # the file, do nothing.
            if overwrite is False and not exists(new_pdb_file):
                pdbl = PDBList()
                pdbl.retrieve_pdb_file(pdb_id, pdir=output_path,
                                       file_format="pdb", overwrite=overwrite)
                logger.debug("Download of the PDB '%s' completed." % pdb_id)

                # Rename files.
                cur_pdb_file = '%s/pdb%s.ent' % (output_path, pdb_id.lower())
                rename_pdb_file(cur_pdb_file, new_pdb_file)
            else:
                logger.debug("File '%s' already exists. It will not be "
                             "downloaded again." % new_pdb_file)
    else:
        raise IllegalArgumentError("An empty PDB id ('%s') was informed."
                                   % pdb_id)


def parse_from_file(pdb_id, file):
    """Read a PDB file and return a
    :class:`~luna.MyBio.PDB.Structure.Structure` object.

    Parameters
    ----------
    pdb_id : str
        The structure identifier.
    file : str
        Pathname of the PDB file.

    Returns
    -------
    structure : :class:`~luna.MyBio.PDB.Structure.Structure`
        The parsed PDB file as a
        :class:`~luna.MyBio.PDB.Structure.Structure` object.

    Raises
    ------
    PDBNotReadError
        If the PDB file could not be parsed.
    """
    try:
        parser = PDBParser(PERMISSIVE=True, QUIET=True, FIX_EMPTY_CHAINS=True,
                           FIX_ATOM_NAME_CONFLICT=True, FIX_OBABEL_FLAGS=False)
        structure = parser.get_structure(pdb_id, file)
        return structure
    except Exception as e:
        logger.exception(e)
        raise PDBNotReadError("File '%s' could not be parsed." % file)


def save_to_file(entity, output_file, select=Select(), write_conects=True,
                 write_end=True, preserve_atom_numbering=True, sort=False):
    """Write a Structure object (or a subset of a
    :class:`~luna.MyBio.PDB.Structure.Structure` object) into a file.

    Parameters
    ----------
    entity : :class:`~luna.MyBio.PDB.Entity.Entity`
        The PDB object to be saved.
    output_file : str
        Save the selected atoms to this file.
    select : :class:`~luna.MyBio.PDB.PDBIO.Select`
        Decide which atoms will be saved at the PDB output.
        By default, all atoms are accepted.
    write_conects : bool
        If True, write CONECT records.
    write_end : bool
        If True, write the END record.
    preserve_atom_numbering : bool
        If True, preserve the atom numbering.
        Otherwise, re-enumerate the atom serial numbers.
    """
    try:
        io = PDBIO()
        io.set_structure(entity)
        io.save(output_file, select=select,
                write_conects=write_conects,
                write_end=write_end,
                preserve_atom_numbering=preserve_atom_numbering,
                sort=sort)

    except Exception as e:
        logger.exception(e)
        raise FileNotCreated("PDB file '%s' could not be created."
                             % output_file)


def entity_to_string(entity, select=Select(),
                     write_conects=True,
                     write_end=True,
                     preserve_atom_numbering=True):
    """Convert a Structure object (or a subset of a
    :class:`~luna.MyBio.PDB.Structure.Structure` object) to string.

    This function works on a structural level. That means if ``entity`` is not
    a :class:`~luna.MyBio.PDB.Structure.Structure` object, the structure will
    be recovered directly from ``entity``. Therefore, use ``select`` to select
    specific chains, residues, and atoms from the structure object.

    Parameters
    ----------
    entity : :class:`~luna.MyBio.PDB.Entity.Entity`
        The PDB object to be converted.
    select : :class:`~luna.MyBio.PDB.PDBIO.Select`
        Decides which atoms will be saved at the output.
        By default, all atoms are accepted.
    write_conects : bool
        If True, writes CONECT records.
    write_end : bool
        If True, writes the END record.
    preserve_atom_numbering : bool
        If True, preserve the atom numbering.
        Otherwise, re-enumerate the atom serial numbers.

    Returns
    -------
    str
        The converted ``entity``.
    """
    fh = StringIO()
    io = PDBIO()
    io.set_structure(entity)
    io.save(fh, select=select,
            write_conects=write_conects,
            write_end=write_end,
            preserve_atom_numbering=preserve_atom_numbering)
    fh.seek(0)
    return ''.join(fh.readlines())


def get_entity_from_entry(entity, entry, model=0):
    """Get a :class:`~luna.MyBio.PDB.Entity.Chain` or
    :class:`~luna.MyBio.PDB.Residue.Residue` instance based on the
    provided entry ``entry``.

    Parameters
    ----------
    entity : :class:`~luna.MyBio.PDB.Entity.Entity`
        The PDB object to recover the target entry from.
    entry : :class:`~luna.entry.Entry`
        The target entry.
    model : int
        The PDB model where the entry could be found. The default value is 0.

    Returns
    -------
    :class:`~luna.MyBio.PDB.Entity.Entity`
        The target entry.

    Raises
    ------
    MoleculeNotFoundError
        If the entry's molecule was not found in ``entity``.
    ChainNotFoundError
        If the entry's chain was not found in ``entity``.
    """
    structure = entity.get_parent_by_level("S")
    model = structure[model]

    if entry.chain_id in model.child_dict:
        chain = model[entry.chain_id]

        if entry.comp_name is not None and entry.comp_num is not None:
            ligand_key = entry.get_biopython_key()

            if ligand_key in chain.child_dict:
                target_entity = chain[ligand_key]
            else:
                error_msg = ("Ligand '%s' does not exist in PDB '%s'."
                             % (entry.to_string(ENTRY_SEPARATOR),
                                structure.get_id()))
                raise MoleculeNotFoundError(error_msg)
        else:
            target_entity = chain
    else:
        error_msg = ("The informed chain id '%s' for the ligand entry '%s' "
                     "does not exist in the PDB '%s'."
                     % (entry.chain_id, entry.to_string(ENTRY_SEPARATOR),
                        structure.get_id()))
        raise ChainNotFoundError()

    return target_entity


def biopython_entity_to_mol(entity,
                            select=Select(),
                            amend_mol=True,
                            template=None,
                            add_h=False,
                            ph=None,
                            metals_coord=None,
                            wrapped=True,
                            openbabel=OPENBABEL,
                            tmp_path=None,
                            keep_tmp_files=False):
    """Convert an object :class:`~luna.MyBio.PDB.Entity.Entity` to a
    molecular object (:class:`~luna.wrappers.base.MolWrapper`,
    :class:`rdkit.Chem.rdchem.Mol`, or :class:`openbabel.pybel.Molecule`).

    Parameters
    ----------
    entity : :class:`~luna.MyBio.PDB.Entity.Entity`
        The entity to be converted.
    select : :class:`~luna.MyBio.PDB.PDBIO.Select`
        Decide which atoms will be consired.
        By default, all atoms are accepted.
    amend_mol : bool
        If True, validate and standardize molecules.
        The validation is employed with
        :class:`~luna.mol.validator.MolValidator`.
        Currently, only residues are standardized.
    template : :class:`luna.mol.template.Template`, optional
        Fix the converted molecule's bond order based on the bond orders
        in a template molecule. The ``template`` should implement
        assign_bond_order().
    add_h : bool
        If True, add hydrogen to the converted molecule.
    ph : float, optional
        Add hydrogens considering pH ``ph``.
    wrapped : bool
        If True, wrap the molecular object with
        :class:`~luna.wrappers.base.MolWrapper`.
    openbabel : str
        Pathname to Open Babel.
    tmp_path : str
        A temporary directory to where temporary files will be saved.
        If not provided, the system's default temporary directory will
        be used instead.
    keep_tmp_files : bool
        If True, keep all temporary files. Otherwise, removes them in the end.

    Returns
    -------
    mol_obj : :class:`~luna.wrappers.base.MolWrapper`, \
                :class:`rdkit.Chem.rdchem.Mol`, \
                or :class:`openbabel.pybel.Molecule`
        The converted molecule.
    ignored_atoms : list
        List of ignored atoms. Currently, ignored atoms may
        contain only metals.
    """

    tmp_path = tmp_path or tempfile.gettempdir()

    logger.debug("It will try to create a new MOL object from the provided "
                 "entity.")
    logger.debug("Temporary files will be saved at '%s'." % tmp_path)

    # First it saves the selection into a PDB file and then it converts the
    # file to .mol. I had to do it because the OpenBabel 2.4.1 had a
    # problem with some molecules containing aromatic rings. In such cases,
    # the aromatic ring was being wrongly perceived and some atoms received
    # more double bonds than it was expected. The version 2.3.2 works better.
    # Therefore, I recomend using Open Babel 2.3.3 instead.
    filename = new_unique_filename(tmp_path)
    pdb_file = '%s_pdb-file.pdb' % filename

    logger.debug("First: it will try to create a new PDB file (%s) "
                 "from the provided entity." % pdb_file)
    # Apparently, Open Babel creates a bug when it tries to parse a file with
    # CONECTS containing serial numbers with more than 4 digits.
    # E.g.: 1OZH:A:HE3:1406, line CONECT162811627916282.
    # By setting preserve_atom_numbering to False, it solves the problem.
    save_to_file(entity,
                 pdb_file,
                 select,
                 preserve_atom_numbering=False,
                 sort=True)

    ini_input_file = pdb_file
    if template is not None:
        if entity.level == "R" and entity.is_hetatm():
            # Note that the template molecule should have no explicit hydrogens
            # else the algorithm will fail.
            rdmol = read_mol_from_file(pdb_file, mol_format="pdb",
                                       removeHs=True)
            new_rdmol = template.assign_bond_order(rdmol, entity.resname)

            ini_input_file = '%s_tmp-mol-file.mol' % filename
            MolToMolFile(new_rdmol, ini_input_file)

            if not keep_tmp_files:
                remove_files([pdb_file])
        else:
            logger.warning("It cannot apply a template on the provided entity "
                           "because it should be a single compound "
                           "(Residue class).")

    # Convert the PDB file to Mol file with the proper protonation
    # and hydrogen addition if required.
    mol_file = '%s_mol-file.mol' % filename
    ob_opt = {"error-level": 5}
    logger.debug("Next: it will try to convert the PDB file to "
                 " .mol using Open Babel.")
    if add_h:
        logger.debug("Hydrogens will be added to the molecule.")
        if ph is not None:
            ob_opt["p"] = ph
        else:
            ob_opt["h"] = ""
    convert_molecule(ini_input_file, output_file=mol_file,
                     opts=ob_opt, openbabel=openbabel)

    # Currently, ignored atoms are only metals.
    ignored_atoms = []

    mol_obj = None
    if amend_mol:
        logger.debug("Now: a validation will be performed and "
                     "it will try to fix some errors.")

        try:
            mol_obj = next(readfile("mol", mol_file))
        except Exception:
            error_msg = ("An error occurred while parsing the file '%s' with "
                         "Open Babel and the molecule object could not be "
                         "created. Check the logs for more information."
                         % mol_file)
            raise MoleculeObjectError(error_msg)

        # Standardize a specific set of atoms if provided,
        # otherwise use all atoms from entity.
        if isinstance(select, AtomSelector):
            target_atoms = select.entries
        else:
            target_atoms = Selection.unfold_entities(entity, 'A')

        # If the add_h property is set to False, the code will not remove
        # any existing hydrogens from the PDB structure.
        # In these situations, the list of atoms may contain hydrogens.
        # But, we do not need to attribute properties to hydrogens.
        # We just need them to correctly set properties to heavy atoms.
        # So let's just ignore them.
        target_atoms = [atm for atm in target_atoms
                        if select.accept_atom(atm) and atm.element != "H"]

        mol_obj = MolWrapper(mol_obj)
        if mol_obj.get_num_heavy_atoms() != len(target_atoms):
            error_msg = ("The number of heavy atoms in the PDB selection "
                         "and in the MOL file are different.")
            raise MoleculeSizeError(error_msg)

        # Ignore hydrogen atoms.
        atm_obj_list = [atm for atm in mol_obj.get_atoms()
                        if atm.get_atomic_num() != 1]

        # Pairs of AtomWrapper and Atom (Biopython) objects
        # to standardize.
        atom_pairs = []
        pdb_mapping = {}
        for i, atm_obj in enumerate(atm_obj_list):
            atom_pairs.append((atm_obj, target_atoms[i]))

            pdb_mapping[atm_obj.get_idx()] = target_atoms[i]

        updated_metals_coord = {}
        # If there is something to be standardised.
        if atom_pairs:
            rs = Standardizer()
            mol_obj.unwrap().DeleteHydrogens()
            rs.standardize(atom_pairs, metals_coord=metals_coord)

            updated_metals_coord = rs.metals_coord

            # After standardizing residues, we need to recreate the Mol
            # file, otherwise implicit hydrogens will not be included
            # in the MOL object and, therefore, their coordinates could
            # not be accessed. If you try to generate coordinates directly
            # from the object, hydrogens will be incorrectly placed.
            mol_obj = PybelWrapper(mol_obj.unwrap())
            new_mol_file = '%s_tmp-mol-file.mol' % filename
            mol_obj.write("mol", new_mol_file, overwrite=True)

            # Overwrite mol_file by converting the new molecular file using
            # the user specified parameters. Note that right now it will add
            # explicit hydrogens to the molecules according to the provided pH.
            convert_molecule(new_mol_file, output_file=mol_file, 
                             opts=ob_opt, openbabel=OPENBABEL)

            # Let's finally read the correct and standardized molecular file.
            try:
                mol_obj = next(readfile("mol", mol_file))
            except Exception:
                error_msg = ("An error occurred while parsing the file "
                             "'%s' with Open Babel and the molecule object"
                             "could not be created. Check the logs for "
                             "more information." % mol_file)
                raise MoleculeObjectError(error_msg)

            # Remove temporary files.
            if not keep_tmp_files:
                remove_files([new_mol_file])

        mv = MolValidator(metals_coord=updated_metals_coord)
        is_valid = mv.validate_mol(mol_obj, pdb_mapping)
        logger.debug('Validation finished!!!')

        if not is_valid:
            logger.warning("The molecular file '%s' contain invalid atoms. "
                           "Check the logs for more information." % mol_file)

        # Validate molecule using RDKit sanitization methods.
        try:
            aux_mol = (PybelWrapper(mol_obj.unwrap())
                       if isinstance(mol_obj, MolWrapper)
                       else mol_obj)

            # The sanitization is set off. We will apply it in the next step.
            aux_mol = MolFromMolBlock(aux_mol.write('mol'),
                                      sanitize=False, removeHs=False)
            # Sanitize molecule is applied now, so we will be able to catch
            # the exceptions raised by RDKit, otherwise it would
            # not be possible.
            SanitizeMol(aux_mol, SanitizeFlags.SANITIZE_ALL)
        except Exception:
            error_msg = ("An error occurred while parsing the molecular block "
                         "with RDKit. The block was generated by Open Babel "
                         "from the file '%s'. Check the logs for more "
                         "information." % mol_file)
            raise MoleculeObjectError(error_msg)
    else:
        try:
            # Create a new Mol object.
            mol_obj = next(readfile("mol", mol_file))
        except Exception:
            error_msg = ("An error occurred while parsing the file '%s' and "
                         "the molecule object could not be created. "
                         "Check the logs for more information."
                         % mol_file)
            raise MoleculeObjectError(error_msg)

    # Remove temporary files.
    if not keep_tmp_files:
        remove_files([ini_input_file, mol_file])

    if wrapped:
        mol_obj = MolWrapper(mol_obj)
    elif isinstance(mol_obj, MolWrapper):
        mol_obj = mol_obj.unwrap()

    return mol_obj, ignored_atoms
