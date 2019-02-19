from util.exceptions import IllegalArgumentError
from util.file import get_file_format
from rdkit.Chem import (MolFromMolBlock, MolFromMol2Block)


FORMATS = ("mol2", "mol", "mdl", "sdf", "sd")


def new_mol_from_block(block, mol_format, sanitize=True, removeHs=True):
    if mol_format == "mol2":
        return MolFromMol2Block(block, sanitize=sanitize, removeHs=removeHs)
    elif mol_format in FORMATS:
        return MolFromMolBlock(block, sanitize=sanitize, removeHs=removeHs)
    else:
        raise IllegalArgumentError("Invalid '%s' format. The accepted formats are: %s." % (mol_format, ",".join(FORMATS)))


def read_multimol_file(mol_file, mol_file_format=None, targets=None, sanitize=True, removeHs=True):

    ext = mol_file_format or get_file_format(mol_file)

    if ext not in FORMATS:
        raise IllegalArgumentError("Format '%s' informed or assumed by the filename is invalid. The accepted formats are: %s."
                                   % (ext, ",".join(FORMATS)))

    with open(mol_file, "r") as IN:

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
                            # If a target list is informed, create a new molecule only if it is in the list.
                            if targets is None:
                                # Create a new RDKit object
                                yield(new_mol_from_block("".join(mol), ext, sanitize, removeHs))
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
                mol_id = mol[1].strip() if ext == "mol2" else mol[0].strip()
                # If a target list is informed, create a new molecule only if it is in the list.
                if targets is None or (targets is not None and mol_id in targets):
                    # Create a new RDKit object
                    yield(new_mol_from_block("".join(mol), ext, sanitize, removeHs))
                break

            # If all target compounds were already found, just break the loop.
            if targets is not None and len(targets) == 0:
                break
