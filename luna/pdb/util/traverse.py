from luna.util.default_values import ENTRY_SEPARATOR
from luna.util.exceptions import MoleculeNotFoundError, ChainNotFoundError



ENTITY_LEVEL_NAME = {"A": "Atom",
                     "R": "Residue",
                     "C": "Chain",
                     "M": "Model",
                     "S": "Structure"}


def get_entity_from_entry(entity, entry, model_id=0):
    """
    Retrieve a chain or residue from a structure based on an ``entry`` object.

    Parameters
    ----------
    entity : :class:`Bio.PDB.Entity.Entity`
        The PDB object to recover the target entry from.
    entry : :class:`luna.entry.Entry`
        Entry object defining the chain ID or residue ID.
    model_id : int
        Model id to search for the entry.

    Returns
    -------
    :class:`Bio.PDB.Entity.Entity`
        The matching Biopython entity.

    Raises
    ------
    MoleculeNotFoundError
        If the entry's molecule was not found in ``entity``.
    ChainNotFoundError
        If the entry's chain was not found in ``entity``.
    """

    structure = entity.get_parent_by_level("S")
    model = structure[model_id]

    if entry.chain_id not in model.child_dict:
        raise ChainNotFoundError(
            f"Chain '{entry.chain_id}' for entry '{entry.to_string(ENTRY_SEPARATOR)}' "
            f"not found in structure '{structure.get_id()}'."
        )
        
    chain = model[entry.chain_id]

    if entry.comp_name is not None and entry.comp_num is not None:
        ligand_key = entry.get_biopython_key()

        if ligand_key not in chain.child_dict:
            raise MoleculeNotFoundError(
                f"Ligand '{entry.to_string(ENTRY_SEPARATOR)}' not found in PDB '{structure.get_id()}'."
            )
        return chain[ligand_key]

    return chain
