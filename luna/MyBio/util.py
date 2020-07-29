import tempfile
import logging
from io import StringIO
from os.path import exists
from shutil import move as rename_pdb_file
from itertools import product, combinations

from openbabel import openbabel as ob
from openbabel.pybel import readfile
from openbabel.pybel import Molecule as PybelWrapper

from rdkit.Chem import MolFromMolBlock, MolFromMolFile, SanitizeFlags, SanitizeMol

from luna.util.file import is_directory_valid, get_unique_filename, remove_files
from luna.util.default_values import ENTRY_SEPARATOR, OPENBABEL
from luna.MyBio.selector import AtomSelector
from luna.MyBio.PDB.PDBList import PDBList
from luna.MyBio.PDB.PDBIO import Select
from luna.MyBio.PDB import Selection
from luna.MyBio.PDB.PDBIO import PDBIO
from luna.MyBio.PDB.PDBParser import PDBParser
from luna.mol.validator import MolValidator
from luna.mol.standardiser import ResiduesStandardiser
from luna.mol.wrappers.base import MolWrapper
from luna.mol.wrappers.obabel import convert_molecule

from luna.util.exceptions import (IllegalArgumentError, MoleculeNotFoundError, ChainNotFoundError,
                                  FileNotCreated, PDBNotReadError, MoleculeSizeError, MoleculeObjectError)


logger = logging.getLogger()


def download_pdb(pdb_id, output_path=".", overwrite=False):
    """Download a PDB file from RCSB.org.

        @param pdb_id: 4-symbols structure Id from PDB (e.g. 3J92).
        @type pdb_code: string

        @param output_path: put the PDB file in this directory.
        @type  output_path: string
    """
    logger.debug("Trying to download the PDB '%s' and store it at the directory '%s'." % (pdb_id, output_path))

    try:
        if pdb_id is not None and pdb_id.strip() != "":
            if is_directory_valid(output_path):
                new_pdb_file = "%s/%s.pdb" % (output_path, pdb_id)
                # If the file already exists and is not to overwrite the file, do nothing.
                if overwrite is False and not exists(new_pdb_file):
                    pdbl = PDBList()
                    pdbl.retrieve_pdb_file(pdb_id, pdir=output_path, file_format="pdb", overwrite=overwrite)
                    logger.debug("Download of the PDB '%s' completed." % pdb_id)

                    # Rename files.
                    cur_pdb_file = '%s/pdb%s.ent' % (output_path, pdb_id.lower())
                    rename_pdb_file(cur_pdb_file, new_pdb_file)
                else:
                    logger.debug("File '%s' already exists. It will not be downloaded again." % new_pdb_file)
        else:
            raise IllegalArgumentError("An empty PDB id ('%s') was informed." % pdb_id)
    except Exception as e:
        logger.exception(e)
        raise


def parse_from_file(id, file):
    """Read a PDB file and return a Structure object.

        @param id: the id that will be used for the structure
        @type id: string

        @param file: name of the PDB file
        @type file: string
    """
    try:
        parser = PDBParser(PERMISSIVE=1, QUIET=True)
        structure = parser.get_structure(id, file)
        return structure
    except Exception as e:
        logger.exception(e)
        raise PDBNotReadError("File '%s' not parsed as a PDB file.")


def save_to_file(entity, output_file, select=Select(), write_conects=True, write_end=True,
                 preserve_atom_numbering=True):
    """ Write a Structure object (or a subset of a Structure object) into a file.

        @param entity: the PDB object to be saved
        @type entity: object

        @param output_file: the name of the new PDB file
        @type output_file: string

        @param select: a filtering function. Default: it extracts everything
        @type select: Select

        @param write_conects: decide if it is necessary to write CONECT fields.
        @type write_conects: boolean

        @param write_end: decide if it is necessary to write END fields.
        @type write_end: boolean

        @param preserve_atom_numbering: decide if it is necessary to re-enumerate the atom serial numbers.
        @type preserve_atom_numbering: boolean
    """
    try:
        io = PDBIO()
        io.set_structure(entity)
        io.save(output_file, select=select, write_conects=write_conects, write_end=write_end,
                preserve_atom_numbering=preserve_atom_numbering)
    except Exception as e:
        logger.exception(e)
        raise FileNotCreated("PDB file '%s' could not be created." % output_file)


def entity_to_string(entity, select=Select(), write_conects=True, write_end=True, preserve_atom_numbering=True):
    fh = StringIO()
    io = PDBIO()
    io.set_structure(entity)
    io.save(fh, select=select, write_conects=write_conects, write_end=write_end, preserve_atom_numbering=preserve_atom_numbering)
    fh.seek(0)
    return ''.join(fh.readlines())


def get_entity_level_name():
    return {
        "A": "Atom",
        "R": "Residue",
        "C": "Chain",
        "M": "Model",
        "S": "Structure"
    }


def get_entity_from_entry(entity, entry, model=0):
    structure = entity.get_parent_by_level("S")
    model = structure[model]

    if entry.chain_id in model.child_dict:
        chain = model[entry.chain_id]

        if entry.comp_name is not None and entry.comp_num is not None:
            ligand_key = entry.get_biopython_key()

            if ligand_key in chain.child_dict:
                target_entity = chain[ligand_key]
            else:
                raise MoleculeNotFoundError("Ligand '%s' does not exist in the PDB '%s'."
                                            % (entry.to_string(ENTRY_SEPARATOR), structure.get_id()))
        else:
            target_entity = chain
    else:
        raise ChainNotFoundError("The informed chain id '%s' for the ligand entry '%s' does not exist in the PDB '%s'." %
                                 (entry.chain_id, entry.to_string(ENTRY_SEPARATOR), structure.get_id()))

    return target_entity


def is_covalently_bound(mybio_atm1, mybio_atm2):
    # Distance atom-atom
    dist = mybio_atm1 - mybio_atm2
    # Covalent radius
    cov1 = ob.GetCovalentRad(ob.GetAtomicNum(mybio_atm1.element))
    cov2 = ob.GetCovalentRad(ob.GetAtomicNum(mybio_atm2.element))

    # OpenBabel thresholds.
    if 0.4 <= dist <= cov1 + cov2 + 0.45:
        return True
    return False


def get_residue_cov_bonds(residue, select=Select()):

    # Get valid atoms according to the provided selection function.
    trgt_res_atms = {atm.name: atm for atm in residue.get_atoms() if select.accept_atom(atm)}

    cov_bonds = []
    for atm1, atm2 in combinations(trgt_res_atms.values(), 2):
        if is_covalently_bound(atm1, atm2):
            cov_bonds.append((atm1, atm2))

    return cov_bonds


def get_residue_neighbors(residue, select=Select()):

    # Get valid atoms according to the provided selection function.
    trgt_res_atms = {atm.name: atm for atm in residue.get_atoms() if select.accept_atom(atm)}

    if residue.is_residue():
        if "N" not in trgt_res_atms:
            logger.debug("There is a missing N in the residue %s. It may have been filtered out by the provided "
                         "selection function. So, the predecessor residue cannot be identified." % residue)
        if "C" not in trgt_res_atms:
            logger.debug("There is a missing C in the residue %s. It may have been filtered out by the provided "
                         "selection function. So, the successor residue cannot be identified." % residue)

        # If neither N and C are available in the valid atom list.
        if "N" not in trgt_res_atms and "C" not in trgt_res_atms:
            return {}

        neighbors = {}
        # If the chain has a residue coming before the target residue.
        if residue.idx - 1 >= 0:
            # First residue before the target in the chain list.
            prev_res = residue.parent.child_list[residue.idx - 1]
            # Get valid atoms according to the provided selection function.
            prev_res_atms = {atm.name: atm for atm in prev_res.get_atoms() if select.accept_atom(atm)}

            # A peptide bond exists between the C of one amino acid and the N of another.
            if "C" in prev_res_atms and "N" in trgt_res_atms:
                if is_covalently_bound(trgt_res_atms["N"], prev_res_atms["C"]):
                    neighbors["previous"] = prev_res
                else:
                    logger.debug("The first residue before %s is too distant to fulfill the covalent thresholds. "
                                 "It may be an indication of bad atom positioning or that there are missing residues." % residue)
            elif "C" not in prev_res_atms:
                logger.debug("There is a missing C in %s, the first residue before %s in the chain list. "
                             "It may have been filtered out by the provided selection function. "
                             "So, the predecessor residue cannot be identified." % (prev_res, residue))
        # Otherwise, it could mean that the residue is the first one in the sequence or there are missing residues.
        else:
            logger.debug("The residue %s seems not to have any predecessor residue. "
                         "It may be the first in the chain sequence or there are missing residues." % residue)

        # If the chain has a residue coming after the target residue.
        if residue.idx + 1 < len(residue.parent.child_list):
            # First residue after the target in the chain list.
            next_res = residue.parent.child_list[residue.idx + 1]
            # Get valid atoms according to the provided selection function.
            next_res_atms = {atm.name: atm for atm in next_res.get_atoms() if select.accept_atom(atm)}

            # A peptide bond exists between the C of one amino acid and the N of another.
            if "C" in trgt_res_atms and "N" in next_res_atms:
                if is_covalently_bound(trgt_res_atms["C"], next_res_atms["N"]):
                    neighbors["next"] = next_res
                else:
                    logger.debug("The first residue after %s is too distant to fulfill the covalent thresholds. "
                                 "It may be an indication of bad atom positioning or that there are missing residues." % residue)
            elif "N" in next_res_atms:
                logger.debug("There is a missing N in %s, the first residue after %s in the chain list. "
                             "It may have been filtered out by the provided selection function. "
                             "So, the predecessor residue cannot be identified." % (next_res, residue))
        # Otherwise, it could mean that the residue is the last one in the sequence or there are missing residues.
        else:
            logger.debug("The residue %s seems not to have any successor residue. "
                         "It may be the last in the chain sequence or there are missing residues." % residue)
        return neighbors
    else:
        # First residue before the target in the chain list.
        prev_res = residue.parent.child_list[residue.idx - 1]
        # Get valid atoms according to the provided selection function.
        prev_res_atms = {atm.name: atm for atm in prev_res.get_atoms() if select.accept_atom(atm)}

        neighbors = {}
        # If the chain has a residue coming before the target residue.
        if residue.idx - 1 >= 0:
            # First residue before the target in the chain list.
            prev_res = residue.parent.child_list[residue.idx - 1]
            # Get valid atoms according to the provided selection function.
            prev_res_atms = {atm.name: atm for atm in prev_res.get_atoms() if select.accept_atom(atm)}

            for trgt_atm, prev_atm in product(trgt_res_atms.values(), prev_res_atms.values()):
                if is_covalently_bound(trgt_atm, prev_atm):
                    neighbors["previous"] = prev_res
                    break

            if "previous" not in neighbors:
                logger.debug("The first residue after %s is too distant to fulfill the covalent thresholds. "
                             "It may be an indication of bad atom positioning or that there are missing residues." % residue)
        # Otherwise, it could mean that the residue is the first one in the sequence or there are missing residues.
        else:
            logger.debug("The residue %s seems not to have any predecessor residue. "
                         "It may be the first in the chain sequence or there are missing residues." % residue)

        # If the chain has a residue coming after the target residue.
        if residue.idx + 1 < len(residue.parent.child_list):
            # First residue after the target in the chain list.
            next_res = residue.parent.child_list[residue.idx + 1]
            # Get valid atoms according to the provided selection function.
            next_res_atms = {atm.name: atm for atm in next_res.get_atoms() if select.accept_atom(atm)}

            # Check each pair of atoms for covalently bonded atoms.
            for trgt_atm, next_atm in product(trgt_res_atms.values(), next_res_atms.values()):
                if is_covalently_bound(trgt_atm, next_atm):
                    neighbors["next"] = next_res
                    break

            if "next" not in neighbors:
                logger.debug("The first residue after %s is too distant to fulfill the covalent thresholds. "
                             "It may be an indication of bad atom positioning or that there are missing residues." % residue)
        # Otherwise, it could mean that the residue is the last one in the sequence or there are missing residues.
        else:
            logger.debug("The residue %s seems not to have any successor residue. "
                         "It may be the last in the chain sequence or there are missing residues." % residue)
        return neighbors
    return {}


def biopython_entity_to_mol(entity, select=Select(), validate_mol=True, standardize_mol=True,
                            add_h=False, ph=None, break_metal_bonds=False, mol_obj_type="rdkit", openbabel=OPENBABEL,
                            tmp_path=None, keep_tmp_files=False):

    tmp_path = tmp_path or tempfile.gettempdir()

    logger.debug("It will try to create a new MOL object from the provided entity.")
    logger.debug("Temporary files will be saved at '%s'." % tmp_path)

    # First it saves the selection into a PDB file and then it converts the file to .mol.
    # I had to do it because the OpenBabel 2.4.1 had a problem with some molecules containing aromatic rings.
    # In such cases, the aromatic ring was being wrongly perceived and some atoms received more double bonds than it
    # was expected. The version 2.3.2 works better. Therefore, I defined this version manually (openbabel property).
    filename = get_unique_filename(tmp_path)
    pdb_file = '%s_pdb-file.pdb' % filename

    logger.debug("First: try to create a new PDB file (%s) from the provided entity." % pdb_file)
    # Apparently, Open Babel creates a bug when it tries to parse a file with CONECTS containing serial numbers with more than 4 digits.
    # E.g.: 1OZH:A:HE3:1406, line CONECT162811627916282.
    # By setting preserve_atom_numbering to False, it seems the problem is solved.
    save_to_file(entity, pdb_file, select, preserve_atom_numbering=False)

    # Convert the PDB file to Mol file with the proper protonation and hydrogen addition if required.
    mol_file = '%s_mol-file.mol' % filename
    ob_opt = {"error-level": 5}
    logger.debug("Next: try to convert the PDB file to .mol using Open Babel.")
    if add_h:
        logger.debug("Hydrogens will be added to the molecule.")
        if ph is not None:
            ob_opt["p"] = ph
        else:
            ob_opt["h"] = ""
    convert_molecule(pdb_file, mol_file, opt=ob_opt, openbabel=openbabel)

    # Currently, ignored atoms are only metals.
    ignored_atoms = []

    mol_obj = None
    if validate_mol or standardize_mol:
        logger.debug("Now: a validation will be performed and it will try to fix some errors.")

        try:
            mol_obj = next(readfile("mol", mol_file))
        except Exception:
            raise MoleculeObjectError("An error occurred while parsing the file '%s' with Open Babel and the molecule "
                                      "object could not be created. Check the logs for more information." % mol_file)

        if standardize_mol:
            # Standardize a specific set of atoms if provided, otherwise use all atoms from entity.
            if isinstance(select, AtomSelector):
                target_atoms = select.entries
            else:
                target_atoms = Selection.unfold_entities(entity, 'A')

            # If the add_h property is set to False, the code will not remove any existing hydrogens from the PDB structure.
            # In these situations, the list of atoms may contain hydrogens. But, we do not need to attribute properties to hydrogens.
            # We just need them to correctly set properties to heavy atoms. So let's just ignore them.
            target_atoms = [atm for atm in target_atoms if select.accept_atom(atm) and atm.element != "H"]

            mol_obj = MolWrapper(mol_obj)
            if mol_obj.get_num_heavy_atoms() != len(target_atoms):
                raise MoleculeSizeError("The number of heavy atoms in the PDB selection and in the MOL file are different.")

            # Ignore hydrogen atoms.
            atm_obj_list = [atm for atm in mol_obj.get_atoms() if atm.get_atomic_num() != 1]

            atom_pairs = []
            for i, atm_obj in enumerate(atm_obj_list):
                if target_atoms[i].parent.is_residue():
                    atom_pairs.append((atm_obj, target_atoms[i]))

            rs = ResiduesStandardiser(break_metal_bonds=break_metal_bonds)
            mol_obj.unwrap().DeleteHydrogens()
            rs.standardise(atom_pairs)

            for i, atm_obj in enumerate(atm_obj_list):
                if atm_obj.get_id() in rs.removed_atoms:
                    ignored_atoms.append(target_atoms[i])

            # After standardizing residues, we need to recreate the Mol file, otherwise implicit hydrogens will not be included in the
            # MOL object and, therefore, their coordinates could not be accessed. If you try to generate coordinates directly from the
            # object, hydrogens will be incorrectly placed.
            mol_obj = PybelWrapper(mol_obj.unwrap())
            new_mol_file = '%s_tmp-mol-file.mol' % filename
            mol_obj.write("mol", new_mol_file, overwrite=True)
            # Overwrite mol_file by converting the new molecular file using the user specified parameters.
            # Note that right now it will add explicit hydrogens to the molecules according to the provided pH.
            convert_molecule(new_mol_file, mol_file, opt=ob_opt, openbabel=OPENBABEL)

            # Let's finally read the correct and standardized molecular file.
            try:
                mol_obj = next(readfile("mol", mol_file))
            except Exception:
                raise MoleculeObjectError("An error occurred while parsing the file '%s' with Open Babel and the molecule "
                                          "object could not be created. Check the logs for more information." % mol_file)

            # Remove temporary files.
            if not keep_tmp_files:
                remove_files([new_mol_file])

        mv = MolValidator()
        is_valid = mv.validate_mol(mol_obj)
        logger.debug('Validation finished!!!')

        if not is_valid:
            logger.warning("The molecular file '%s' contain invalid atoms. Check the logs for more information." % mol_file)

        if mol_obj_type == "rdkit":
            try:
                # The sanitization is set off. We will apply it in the next step.
                mol_obj = MolFromMolBlock(mol_obj.write('mol'), sanitize=False, removeHs=False)
                # Sanitize molecule is applied now, so we will be able to catch the exceptions raised by RDKit,
                # otherwise it would not be possible.
                SanitizeMol(mol_obj, SanitizeFlags.SANITIZE_ALL)
            except Exception:
                raise MoleculeObjectError("An error occurred while parsing the molecular block with RDKit. The block was "
                                          "generated by Open Babel from the file '%s'. Check the logs for more information." % mol_file)
    else:
        try:
            # Create a new Mol object.
            if mol_obj_type == "openbabel":
                mol_obj = next(readfile("mol", mol_file))
            else:
                # The sanitization is set off. We will apply it in the next statement.
                mol_obj = MolFromMolFile(mol_file, sanitize=False, removeHs=False)
                # Sanitize molecule is applied now, so we will be able to catch the exceptions raised by RDKit,
                # otherwise it would not be possible.
                SanitizeMol(mol_obj, SanitizeFlags.SANITIZE_ALL)
        except Exception:
            tool = "Open Babel" if mol_obj_type == "openbabel" else "RDKit"
            raise MoleculeObjectError("An error occurred while parsing the file '%s' with %s and the molecule "
                                      "object could not be created. Check the logs for more information." % (mol_file, tool))

    # Remove temporary files.
    if not keep_tmp_files:
        remove_files([pdb_file, mol_file])

    return mol_obj, ignored_atoms
