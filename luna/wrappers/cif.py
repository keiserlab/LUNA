import re


def get_atom_names_by_id(cif_file):
    """Read a single-molecule CIF file and return the molecule's atom names.

    In the current version, if applied on multi-molecular CIF files,
    only the first molecule's atom names are returned.

    Returns
    -------
     : dict
    """

    regex = re.compile(' {1,}')

    atom_names = {}
    with open(cif_file, "r") as IN:
        loop_read = False
        pdbx_ordinal_read = False

        multi_atom = False
        start_read = False

        last_line = -1

        for i, line in enumerate(IN.readlines()):
            line = line.strip()

            if line.startswith("loop_"):
                loop_read = True
            else:
                if line.startswith("_chem_comp_atom.comp_id"):
                    if loop_read and i == last_line + 1:
                        multi_atom = True
                    else:
                        start_read = True
                elif line.startswith("_chem_comp_atom.pdbx_ordinal"):
                    pdbx_ordinal_read = True
                elif pdbx_ordinal_read and line.startswith("#"):
                    break

                if start_read:
                    if multi_atom:
                        cols = regex.split(line)
                        atom_names[len(atom_names)] = cols[1]
                    elif line.startswith("_chem_comp_atom.atom_id"):
                        atom_names[len(atom_names)] = line[line.rfind(' ') + 1:]
                        break

            if pdbx_ordinal_read and multi_atom:
                start_read = True

            last_line = i

    return atom_names
