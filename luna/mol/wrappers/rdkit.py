from rdkit.Chem import MolFromMol2File, MolFromPDBFile, MolFromMolFile, MolFromMolBlock, MolFromMol2Block, SanitizeFlags, SanitizeMol
from xopen import xopen

from luna.util.file import get_file_format
from luna.util.exceptions import IllegalArgumentError

import logging

logger = logging.getLogger()

RDKIT_FORMATS = ("mol2", "mol", "mdl", "sdf", "sd", "pdb")


def read_mol_from_file(mol_file, mol_format, sanitize=True, removeHs=True):
    if mol_format == "mol2":
        # First it creates the molecule without applying the sanitization function.
        rdk_mol = MolFromMol2File(mol_file, sanitize=False, removeHs=removeHs)
    if mol_format == "pdb":
        # First it creates the molecule without applying the sanitization function.
        rdk_mol = MolFromPDBFile(mol_file, sanitize=False, removeHs=removeHs)
    elif mol_format in RDKIT_FORMATS:
        # First it creates the molecule without applying the sanitization function.
        rdk_mol = MolFromMolFile(mol_file, sanitize=False, removeHs=removeHs)
    else:
        raise IllegalArgumentError("Invalid '%s' format. The accepted formats are: %s." % (mol_format, ",".join(RDKIT_FORMATS)))

    # Now it sanitizes the molecule if necessary.
    # We do it here in order to catch an Exception while performing the sanitization.
    if sanitize:
        try:
            SanitizeMol(rdk_mol, SanitizeFlags.SANITIZE_ALL)
        except Exception as e:
            logger.exception(e)
            return None
    return rdk_mol


def new_mol_from_block(block, mol_format, sanitize=True, removeHs=True):
    if mol_format == "mol2":
        # First it creates the molecule without applying the sanitization function.
        rdk_mol = MolFromMol2Block(block, sanitize=False, removeHs=removeHs)
    elif mol_format in RDKIT_FORMATS:
        # First it creates the molecule without applying the sanitization function.
        rdk_mol = MolFromMolBlock(block, sanitize=False, removeHs=removeHs)
    else:
        raise IllegalArgumentError("Invalid '%s' format. The accepted formats are: %s." % (mol_format, ",".join(RDKIT_FORMATS)))

    # Now it sanitizes the molecule if necessary.
    # We do it here in order to catch an Exception while performing the sanitization.
    if sanitize:
        try:
            SanitizeMol(rdk_mol, SanitizeFlags.SANITIZE_ALL)
        except Exception as e:
            logger.exception(e)
            return None
    return rdk_mol


def read_multimol_file(mol_file, targets=None, mol_format=None, sanitize=True, removeHs=True):

    ext = mol_format or get_file_format(mol_file, ignore_compression=True)

    if ext not in RDKIT_FORMATS:
        raise IllegalArgumentError("Format '%s' informed or assumed from the filename is invalid. The accepted formats are: %s."
                                   % (ext, ",".join(RDKIT_FORMATS)))

    with xopen(mol_file, "r") as IN:
        if targets is not None:
            targets = set(targets)

        mol = []
        while True:
            try:
                line = next(IN)
                # Ignore new lines before a molecule block.
                if len(mol) == 0 and line == "\n":
                    continue

                if ext == "mol2":
                    if line.startswith("#"):
                        continue
                    # New molecule block definition is starting...
                    if line.startswith("@<TRIPOS>MOLECULE"):
                        # New molecule identified but an old one already exists.
                        if mol:
                            # If a target list is not informed, create a new molecule object.
                            if targets is None:
                                # Create a new RDKit object
                                yield(new_mol_from_block("".join(mol), ext, sanitize, removeHs))
                            # Otherwise, create a new molecule only if it is in the list.
                            elif mol[1].strip() in targets:
                                targets.remove(mol[1].strip())
                                # Create a new RDKit object
                                yield(new_mol_from_block("".join(mol), ext, sanitize, removeHs))
                        # Restart the molecule block.
                        mol = []
                    mol.append(line)
                else:
                    if line.startswith("M  END"):
                        mol.append(line)
                        # If a target list is informed, create a new molecule only if it is in the list.
                        if targets is None:
                            # Create a new RDKit object
                            yield(new_mol_from_block("".join(mol), ext, sanitize, removeHs))
                        elif mol[0].strip() in targets:
                            # Create a new RDKit object
                            yield(new_mol_from_block("".join(mol), ext, sanitize, removeHs))
                            targets.remove(mol[0].strip())
                        # Restart the molecule block.
                        mol = []
                    elif line.startswith("$$$$") is False:
                        mol.append(line)
            except StopIteration:
                if mol:
                    mol_id = mol[1].strip() if ext == "mol2" else mol[0].strip()
                    # If a target list is informed, create a new molecule only if it is in the list.
                    if targets is None or (targets is not None and mol_id in targets):
                        # Create a new RDKit object
                        yield(new_mol_from_block("".join(mol), ext, sanitize, removeHs))
                break

            # If all target compounds were already found, just break the loop.
            if targets is not None and len(targets) == 0:
                break
