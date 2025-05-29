from Bio.PDB.Chain import Chain as BioChain

from luna.pdb.core.entity import Entity


class Chain(BioChain, Entity):
    """
    Custom Chain subclass with support for SMCRA traversal
    and target flagging for processing.

    Parameters
    ----------
    id : str
        Chain identifier (e.g., "A").
    """

    def __init__(self, id):
        super().__init__(id)
        
        # By default: no chain is a target for processing.
        self._is_target = False

    def __repr__(self):
        return f"<Chain id={self.id}>"

    def __lt__(self, other):
        return self.id < other.id

    def is_target(self):
        """Return True if the chain is marked as a processing target."""
        return self._is_target

    def set_as_target(self, is_target=True):
        """
        Mark or unmark the chain and all its residues as targets.

        Parameters
        ----------
        is_target : bool
            Whether to mark the chain and its residues as a target.
        """
        self._is_target = is_target
        for res in self.get_residues():
            if hasattr(res, "set_as_target"):
                res.set_as_target(is_target)
