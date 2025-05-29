import logging
from os.path import exists
from pathlib import Path
from shutil import move as rename_pdb_file

from Bio.PDB import PDBList

from luna.util.file import is_directory_valid
from luna.util.exceptions import IllegalArgumentError


logger = logging.getLogger(__name__)


PDB_SERVER = "https://files.wwpdb.org"


def download_pdb(pdb_id, 
                 output_path=".", 
                 overwrite=False,
                 server=PDB_SERVER):
    """Download a PDB file from RCSB.org.

    Parameters
    ----------
    pdb_id : str
        4-symbols structure Id from PDB (e.g. 3J92).
    output_path : str
        Directory where the PDB file should be saved.
    overwrite : bool
        If False, skip download if file already exists.
    """
    
    if not pdb_id or not pdb_id.strip():
        raise IllegalArgumentError(f"An empty PDB ID ('{pdb_id}') was provided.")

    pdb_id = pdb_id.strip()
    output_path = Path(output_path)
    
    if not is_directory_valid(output_path):
        raise IllegalArgumentError(f"Output path '{output_path}' is not valid or writable.")

    output_path.mkdir(parents=True, exist_ok=True)
    target_path = output_path / f"{pdb_id}.pdb"

    if target_path.exists() and not overwrite:
        logger.debug(f"PDB file '{target_path}' already exists. Skipping download.")
        return str(target_path)

    logger.debug(f"Downloading PDB '{pdb_id}' to '{target_path}'...")

    pdbl = PDBList(server=PDB_SERVER)
    pdbl.retrieve_pdb_file(pdb_id, 
                           pdir=output_path,
                           file_format="pdb", 
                           overwrite=overwrite)
    
    downloaded_file = output_path / f"pdb{pdb_id.lower()}.ent"
    if downloaded_file.exists():
        rename_pdb_file(str(downloaded_file), str(target_path))
        logger.debug(f"PDB '{pdb_id}' downloaded and renamed to '{target_path}'.")
    else:
        logger.warning(f"Expected file '{downloaded_file}' not found after download.")

    return str(target_path)
