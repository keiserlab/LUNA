from input.util import format_selector


def check_molecule_existance(model, chain, molecule):

    if (chain not in model.child_dict):
            return (False, "Chain %s not found in PDB %s" %
                    (chain, model.get_parent().id))
    elif (molecule not in model[chain].child_dict):
        return (False, "Molecule %s not found in chain %s.%s" %
                (format_selector(molecule), chain, model.get_parent().id))
    else:
        return True, None
