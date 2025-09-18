from Bio.PDB.Model import Model as BioModel

from luna.pdb.core.entity import Entity


class Model(BioModel, Entity):
    """
    Custom Model subclass with SMCRA traversal support.

    Parameters
    ----------
    id : int
        Model index.
    serial_num : int, optional
        Optional serial number for PDB output.
    """

    def __init__(self, id, serial_num=None):
        super().__init__(id, serial_num)

    def __repr__(self):
        return f"<Model (LUNA subclass) id={self.id}>"

    def __lt__(self, other):
        """Allow sorting of models by ID."""
        return self.id < other.id

    def get_atoms(self):
        """Return atoms."""
        for r in self.get_residues():
            for a in r.get_unpacked_list():
                yield a
