from util.exceptions import IllegalArgumentError
from file.validator import is_directory_valid

from MyBio.PDB.PDBList import PDBList
from MyBio.PDB.PDBIO import Select
from MyBio.PDB.PDBIO import PDBIO
from util.exceptions import FileNotCreated

from io import StringIO
from shutil import move as rename_pdb_file

import logging
logger = logging.getLogger(__name__)


def download_pdb(pdb_id, output_path=".", output_file=None, overwrite=False):
    """Download a PDB file from RCSB.org.

        @param pdb_id: 4-symbols structure Id from PDB (e.g. 3J92).
        @type pdb_code: string

        @param output_path: put the PDB file in this directory.
        @type  output_path: string
    """
    logger.info("Trying to download the PDB '%s' and store it at the "
                "directory '%s'." % (pdb_id, output_path))

    try:
        pdb_id = pdb_id.lower()
        if (pdb_id is not None and pdb_id.strip() != ""):
            if (is_directory_valid(output_path)):
                pdbl = PDBList()
                pdbl.retrieve_pdb_file(pdb_id, pdir=output_path,
                                       file_format="pdb", overwrite=overwrite)

                if output_file:
                    pdb_file = '%s/pdb%s.ent' % (output_path, pdb_id)
                    rename_pdb_file(pdb_file, output_file)
        else:
            raise IllegalArgumentError("Inform a non empty PDB id")
    except Exception as e:
        logger.exception(e)
        raise

    logger.info("Download complete!!")


def try_parse_from_pdb(id, file):
    """Read a PDB file and return a Structure object.

        @param id: the id that will be used for the structure
        @type id: string

        @param file: name of the PDB file
        @type file: string
    """
    from MyBio.PDB.PDBParser import PDBParser
    from util.exceptions import PDBNotReadError

    try:
        parser = PDBParser(PERMISSIVE=1, QUIET=True)
        structure = parser.get_structure(id, file)
        return structure
    except Exception as e:
        logger.exception(e)
        raise PDBNotReadError("File '%s' not parsed as a PDB file.")


def try_save_2pdb(pdbObject, outputFile, selector=Select()):
    """Save a PDB file applying some filtering in it.

        @param pdbObject: the PDB object to be saved
        @type pdbObject: object

        @param outputFile: the name of the new PDB file
        @type outputFile: string

        @param selector: a filtering definition. DEFAULT: extracts everything
        @type selector: Select
    """
    try:
        io = PDBIO()
        io.set_structure(pdbObject)
        io.save(outputFile, selector,
                preserve_atom_numbering=True)
    except Exception as e:
        logger.exception(e)
        raise FileNotCreated("PDB file '%s' could not be created."
                             % outputFile)


def pdb_object_2block(entity, select=Select(),
                      preserve_atom_numbering=True,
                      write_conects=True):

    fh = StringIO()
    io = PDBIO()
    io.set_structure(entity)
    io.save(fh, select=select, preserve_atom_numbering=preserve_atom_numbering,
            write_conects=write_conects)
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
