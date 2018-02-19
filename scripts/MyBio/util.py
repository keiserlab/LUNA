from util.exceptions import IllegalArgumentError
from file.validator import is_directory_valid

from MyBio.PDB.Select import Select

import logging
logger = logging.getLogger(__name__)


def download_pdb(pdbId, outputPath="."):
    """Download a PDB file from RCSB.org.

        @param pdbId: 4-symbols structure Id from PDB (e.g. 3J92).
        @type pdb_code: string

        @param outputPath: put the PDB file in this directory.
        @type  outputPath: string
    """
    from MyBio.PDB.PDBList import PDBList

    logger.info("Trying to download the PDB '%s' and store it at the "
                "directory '%s'." % (pdbId, outputPath))

    try:

        if (pdbId is not None and pdbId.strip() != ""):
            if (is_directory_valid(outputPath)):
                pdbl = PDBList()
                pdbl.retrieve_pdb_file(pdbId, pdir=outputPath,
                                       file_format="pdb", overwrite=True)
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
    from MyBio.PDB.PDBIO import PDBIO
    from util.exceptions import FileNotCreated

    try:
        io = PDBIO()
        io.set_structure(pdbObject)
        io.save(outputFile, selector,
                preserve_atom_numbering=True)
    except Exception as e:
        logger.exception(e)
        raise FileNotCreated("PDB file '%s' could not be created."
                             % outputFile)
