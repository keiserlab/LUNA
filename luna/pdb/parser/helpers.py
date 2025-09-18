import logging
from pathlib import Path
from io import StringIO

from luna.pdb.parser.base import PDBParser
from luna.util.exceptions import PDBNotReadError


logger = logging.getLogger()


def load_from_file(pdb_file, pdb_id=None):
    """
    Load a PDB structure from file using the default PDBParser configuration.

    Parameters
    ----------
    pdb_file : str or Path
        Path to the PDB file.
    pdb_id : str or None
        The structure identifier. If None, the filename will be used as id.
 
    Returns
    -------
    Structure
        The parsed Biopython structure object.

    Raises
    ------
    PDBNotReadError
        If the file cannot be parsed as a valid PDB structure.
    """
    
    try:
        pdb_id = pdb_id or Path(pdb_file).stem

        parser = PDBParser(PERMISSIVE=True, 
                           QUIET=True, 
                           FIX_EMPTY_CHAINS=True,
                           FIX_ATOM_NAME_CONFLICT=True, 
                           FIX_OBABEL_FLAGS=False)
        return parser.get_structure(pdb_id, pdb_file)
        
    except Exception as e:
        logger.exception(f"Failed to parse file '{pdb_file}': {str(e)}")
        raise PDBNotReadError(f"File '{pdb_file}' could not be parsed.")


def load_from_string(pdb_id: str, pdb_string: str):
    """
    Load a PDB structure from a string using the default PDBParser configuration.

    Parameters
    ----------
    pdb_id : str
        The structure identifier.
    pdb_string : str
        The full PDB content as a string.

    Returns
    -------
    Structure
        The parsed Biopython structure/

    Raises
    ------
    PDBNotReadError
        If the PDB string could not be parsed.
    """
    try:
        parser = PDBParser(PERMISSIVE=True,
                           QUIET=True,
                           FIX_EMPTY_CHAINS=True,
                           FIX_ATOM_NAME_CONFLICT=True,
                           FIX_OBABEL_FLAGS=False)

        handle = StringIO(pdb_string)
        return parser.get_structure_from_pdb_block(pdb_id, handle.read())

    except Exception as e:
        logger.exception(f"Failed to parse PDB string for '{pdb_id}': {str(e)}")
        raise PDBNotReadError(f"PDB string for '{pdb_id}' could not be parsed.")
