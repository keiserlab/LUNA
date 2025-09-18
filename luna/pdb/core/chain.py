from Bio.PDB.Chain import Chain as BioChain

from luna.pdb.core.entity import Entity


class Chain(BioChain, Entity):
    """
    Custom Chain subclass with support for SMCRA traversal
    and reference flagging for processing.

    Parameters
    ----------
    id : str
        Chain identifier (e.g., "A").
    """

    def __init__(self, id):
        super().__init__(id)
        
        # By default: no chain is a reference for processing.
        self._is_reference = False

    def __repr__(self):
        return f"<Chain (LUNA subclass) id={self.id}>"

    def __lt__(self, other):
        return self.id < other.id

    def is_reference(self):
        """
        Return True if this entity is marked as the reference for interaction analysis.

        This flag indicates that the user selected this entity
        as the *main subject* of interaction evaluation — interactions will be computed
        between this entity and others nearby (e.g., ligands, water, ions).

        Returns
        -------
        bool
            True if this entity is marked for interaction reference.
        """
        return self._is_reference

    def set_as_reference(self, is_reference=True):
        """
        Mark or unmark the chain and all its residues as references.

        This method sets a flag indicating that this chain should be 
        treated as the main subject for calculating interactions 
        with nearby molecules.

        Parameters
        ----------
        is_reference : bool, optional
            If True, mark the entity as the interaction reference.
            If False, unmark it. Defaults to True.

        See Also
        --------
        is_reference : Check whether this entity is currently marked.
        """
        self._is_reference = is_reference
        for res in self.get_residues():
            if hasattr(res, "set_as_reference"):
                res.set_as_reference(is_reference)
