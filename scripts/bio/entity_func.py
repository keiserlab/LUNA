def get_parent_by_level(atom, level):

    if (level == "A"):
        return atom
    elif (level == "R"):
        return atom.get_parent()
    elif (level == "C"):
        return get_parent_by_level(atom, 'R').get_parent()
    elif (level == "M"):
        return get_parent_by_level(atom, 'C').get_parent()
    elif (level == "S"):
        return get_parent_by_level(atom, 'M').get_parent()


def get_entity_level_name():
    return {
        "A": "atom",
        "R": "residue",
        "C": "chain",
        "M": "model",
        "S": "structure"
    }
