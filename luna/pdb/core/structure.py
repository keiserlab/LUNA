from Bio.PDB.Structure import Structure as BioStructure

from luna.pdb.core.entity import Entity


class Structure(BioStructure, Entity):
    """
    Custom Structure subclass for macromolecular structures.

    Supports SMCRA traversal and stores the original file path.

    Parameters
    ----------
    id : str
        Structure identifier.
    pdb_file : str, optional
        Path to the original PDB file.
    """

    def __init__(self, id, pdb_file=None):
        super().__init__(id)
        
        self.pdb_file = pdb_file

    def __repr__(self):
        return f"<Structure (LUNA subclass) id={self.id}>"
