from rdkit.Chem import MolFromMol2File, MolFromPDBFile, MolFromMolFile, MolFromMolBlock, MolFromMol2Block, SanitizeFlags, SanitizeMol
from xopen import xopen

from luna.util.file import get_file_format
from luna.util.exceptions import IllegalArgumentError

import re
import logging

logger = logging.getLogger()

RDKIT_FORMATS = ("mol2", "mol", "mdl", "sdf", "sd", "pdb")

REGEX_MOL_FILE = re.compile(r' (V2000|V3000)$')


def read_mol_from_file(mol_file, mol_format, sanitize=True, removeHs=True):
    if mol_format == "mol2":
        # First it creates the molecule without applying the sanitization function.
        rdk_mol = MolFromMol2File(mol_file, sanitize=False, removeHs=removeHs)
    elif mol_format == "pdb":
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

    def apply_mol_format(lines):
        # Check if the MOL file contains a valid header comprising a Title line, Program/file timestamp line, and a Comment line.
        # Thus, the Counts line must be in line 4 (index 3). If it's not, add blank lines to match the format.
        # This amendment is necessary otherwise RDKit will crash.
        if not REGEX_MOL_FILE.search(lines[3]):
            logger.warning("While parsing the file '%s', we found a molecule starting at line #%d that does not contain a valid header. "
                           "We will add empty lines to its header to match the MOL file format." % (mol_file, mol_starts_at))

            counts_line_pos = None
            for i, line in enumerate(lines):
                if REGEX_MOL_FILE.search(line):
                    counts_line_pos = i
                    break

            missing_lines = ["\n"] * (3 - counts_line_pos)
            lines = missing_lines + lines

        if lines[-1].strip() != "M  END":
            logger.warning("While parsing the file '%s', we found a molecule starting at line #%d that does not contain the END line. "
                           "We will add it to the block end to match the MOL file format." % (mol_file, mol_starts_at))
            lines.append("M  END\n")
        return lines

    ext = mol_format or get_file_format(mol_file, ignore_compression=True)

    if ext not in RDKIT_FORMATS:
        raise IllegalArgumentError("Format '%s' informed or assumed from the filename is invalid. The accepted formats are: %s."
                                   % (ext, ",".join(RDKIT_FORMATS)))

    with xopen(mol_file, "r") as IN:
        if targets is not None:
            targets = set(targets)

        mol = []
        ignore_lines = False
        line_count = 0
        while True:
            try:
                line = IN.readline()
                line_count += 1
                # readline() returns empty strings when EOF is reached.
                if not line:
                    raise StopIteration
                # Ignore new lines before a molecule block.
                if len(mol) == 0 and line == "\n":
                    continue
                # Save the line in which the new molecule starts.
                if len(mol) == 0:
                    mol_starts_at = line_count

                if ext == "mol2":
                    if line.startswith("#"):
                        continue
                    # New molecule block definition is starting...
                    if line.startswith("@<TRIPOS>MOLECULE"):
                        # New molecule identified but an old one already exists.
                        if mol:
                            mol_id = mol[1].strip()
                            # If a target list is not informed, create a new molecule object.
                            if targets is None:
                                # Create a new RDKit object
                                rdk_mol = new_mol_from_block("".join(mol), ext, sanitize, removeHs)
                                yield ((rdk_mol, mol_id))
                            # Otherwise, create a new molecule only if it is in the list.
                            elif mol_id in targets:
                                targets.remove(mol_id)
                                # Create a new RDKit object
                                rdk_mol = new_mol_from_block("".join(mol), ext, sanitize, removeHs)
                                yield ((rdk_mol, mol_id))
                        # Restart the molecule block.
                        mol = []
                    mol.append(line)
                else:
                    if line.startswith("M  END"):
                        mol.append(line)
                        # Fix the MOL block in cases where the header or end lines do not match the MOL file format.
                        mol = apply_mol_format(mol)
                        mol_id = mol[0].strip()

                        # If a target list is informed, create a new molecule only if it is in the list.
                        if targets is None:
                            # Create a new RDKit object
                            rdk_mol = new_mol_from_block("".join(mol), ext, sanitize, removeHs)
                            yield((rdk_mol, mol_id))

                        elif mol[0].strip() in targets:
                            targets.remove(mol_id)
                            # Create a new RDKit object
                            rdk_mol = new_mol_from_block("".join(mol), ext, sanitize, removeHs)
                            yield((rdk_mol, mol_id))
                        # Restart the molecule block.
                        mol = []
                        # After finding a molecule block, ignore any following lines until it finds a line "$$$$".
                        ignore_lines = True
                    elif line.startswith("$$$$") is False:
                        # Ignore lines util it finds a line "$$$$".
                        if ignore_lines:
                            continue
                        mol.append(line)
                    else:
                        ignore_lines = False

            except StopIteration:
                if mol:
                    if ext == "mol":
                        # Fix the MOL block in cases where the header or end lines do not match the MOL file format.
                        mol = apply_mol_format(mol)

                    mol_id = mol[1].strip() if ext == "mol2" else mol[0].strip()

                    # If a target list is informed, create a new molecule only if it is in the list.
                    if targets is None or (targets is not None and mol_id in targets):
                        # Create a new RDKit object
                        rdk_mol = new_mol_from_block("".join(mol), ext, sanitize, removeHs)
                        yield ((rdk_mol, mol_id))
                break

            # If all target compounds were already found, just break the loop.
            if targets is not None and len(targets) == 0:
                break
