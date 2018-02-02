from util.exceptions import IllegalArgumentError
from file.validator import is_directory_valid

from Bio.PDB.PDBIO import Select

import logging
logger = logging.getLogger(__name__)


def download_pdb(pdbId, outputPath="."):
    """Download a PDB file from RCSB.org.

        @param pdbId: 4-symbols structure Id from PDB (e.g. 3J92).
        @type pdb_code: string

        @param outputPath: put the PDB file in this directory.
        @type  outputPath: string
    """
    from bio.PDBList import PDBList

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
    from Bio.PDB.PDBParser import PDBParser
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
    from Bio.PDB.PDBIO import PDBIO
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


class VersionControl():

    def __init__(self):
        from bio.PDBList import PDBList
        from urllib.error import URLError, HTTPError
        import logging
        logger = logging.getLogger(__name__)

        pdbl = PDBList()
        try:
            new, modified, obsolete = pdbl.get_recent_changes()
        except HTTPError as e:
            logger.exception(e)
            raise HTTPError("PDB server could not be reached.")
        except URLError as e:
            logger.exception(e)
            raise URLError("PDB server could not be reached.")

        self.new = set(new)
        self.modified = set(modified)
        self.obsolete = set(obsolete)

    def is_pdb_outdated(self, pdbId):
        if (pdbId in self.new or pdbId in self.modified):
            return True
        else:
            return False

    def is_pdb_obsolete(self, pdbId):
        if (pdbId in self.obsolete):
            return True
        else:
            return False


class ChainSelector(Select):
    def __init__(self, entries):
        self.entries = entries

    def accept_chain(self, chain):
        return True if (chain in self.entries) else False


class Extractor():
    """Selects everything and creates a new PDB output.
    """

    def __init__(self, model):
        self.model = model

    def extract_chains(self, chains, outputFile):
        selChains = set()
        for chain in chains:
            if (chain in self.model.child_dict):
                selChains.add(self.model[chain])
            else:
                logger.warning("Chain %s does not exist." % chain)

        try_save_2pdb(self.model, outputFile, ChainSelector(selChains))
