from util.exceptions import (IllegalArgumentError, MoleculeNotFoundError, ChainNotFoundError)
from util.file import is_directory_valid
from util.default_values import ENTRY_SEPARATOR

from MyBio.PDB.PDBList import PDBList
from MyBio.PDB.PDBIO import Select
from MyBio.PDB.PDBIO import PDBIO
from MyBio.PDB.PDBParser import PDBParser
from util.exceptions import (FileNotCreated, PDBNotReadError)

from io import StringIO
from shutil import move as rename_pdb_file

import logging

logger = logging.getLogger()


def download_pdb(pdb_id, output_path=".", output_file=None, overwrite=False):
    """Download a PDB file from RCSB.org.

        @param pdb_id: 4-symbols structure Id from PDB (e.g. 3J92).
        @type pdb_code: string

        @param output_path: put the PDB file in this directory.
        @type  output_path: string
    """
    logger.info("Trying to download the PDB '%s' and store it at the directory '%s'."
                % (pdb_id, output_path))

    try:
        pdb_id = pdb_id.lower()
        if (pdb_id is not None and pdb_id.strip() != ""):
            if (is_directory_valid(output_path)):
                pdbl = PDBList()
                pdbl.retrieve_pdb_file(pdb_id, pdir=output_path, file_format="pdb", overwrite=overwrite)

                if output_file:
                    pdb_file = '%s/pdb%s.ent' % (output_path, pdb_id)
                    rename_pdb_file(pdb_file, output_file)
        else:
            raise IllegalArgumentError("Inform a non empty PDB id")
    except Exception as e:
        logger.exception(e)
        raise

    logger.info("Download complete!!")


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
    io.save(fh, select=select, write_conects=write_conects, write_end=write_end,
            preserve_atom_numbering=preserve_atom_numbering)
    fh.seek(0)
    return ''.join(fh.readlines())


def get_entity_level_name():
    return {
        "A": "atom",
        "R": "residue",
        "C": "chain",
        "M": "model",
        "S": "structure"
    }


def get_entity_from_entry(entity, entry, model=0):
        structure = entity.get_parent_by_level("S")
        model = structure[model]

        if entry.chain_id in model.child_dict:
            chain = model[entry.chain_id]

            if entry.comp_name and entry.comp_num:
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
