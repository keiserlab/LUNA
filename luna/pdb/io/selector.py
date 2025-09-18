from Bio.PDB.PDBIO import Select


DEFAULT_ALTLOC = ("A", "1")


class Selector(Select):
    """
    Base selector for filtering atoms during PDB output.

    Parameters
    ----------
    keep_hydrog : bool
        If True (default), keeps all hydrogens. Otherwise, filter them out.
        Hydrogens are atoms whose element is either "H" or "D" (deuterium).
    keep_altloc : bool
        If True (default), keeps all atoms. Otherwise, keeps only atoms whose
        alternate location is in the list ``altloc``.
    altloc : iterable
        List of valid alternate location identifiers. 
        The default valid values are "A" and "1".
    """

    def __init__(self, 
                 keep_hydrog=True, 
                 keep_altloc=True, 
                 altloc=DEFAULT_ALTLOC):
        self.keep_hydrog = keep_hydrog
        self.keep_altloc = keep_altloc
        self.altloc = altloc

    def accept_atom(self, atom):
        """Return True if the atom should be written."""
        
        # Hydrogen and Deuterium
        if not self.keep_hydrog and atom.element in ["H", "D"]:
            return False
            
        if self.keep_altloc:
            return True

        return not atom.is_disordered() or atom.get_altloc() in self.altloc


class ChainSelector(Selector):
    """
    Select atoms from a specific set of chains.

    Parameters
    ----------
    entries : iterable of :class:`~luna.pdb.core.chain.Chain`
        Chains to include in the output.
    """

    def __init__(self, entries, **kwargs):
        self.entries = entries
        super().__init__(**kwargs)

    def accept_chain(self, chain):
        """Decide if the chain is valid or not."""
        return True if (chain in self.entries) else False

    def accept_residue(self, res):
        """Decide if the residue is valid or not."""
        return self.accept_chain(res.get_parent())

    def accept_atom(self, atom):
        """Decide if the atom is valid or not."""
        return super().accept_atom(atom) and self.accept_residue(atom.get_parent())


class ResidueSelector(Selector):
    """
    Select atoms from a specific set of residues.

    Parameters
    ----------
    entries : iterable of :class:`~luna.pdb.core.residue.Residue`
        Residues to include in the output.
    """
    
    def __init__(self, entries, **kwargs):
        self.entries = entries
        super().__init__(**kwargs)

    def accept_residue(self, res):
        """Decide if the residue is valid or not."""
        return res in self.entries

    def accept_atom(self, atom):
        """Decide if the atom is valid or not."""
        return super().accept_atom(atom) and self.accept_residue(atom.get_parent())


class ResidueSelectorBySeqNum(Selector):
    """
    Select residues by sequence number.

    Parameters
    ----------
    entries : iterable of int
        Residue sequence numbers to include.
    """

    def __init__(self, entries, **kwargs):
        self.entries = entries
        super().__init__(**kwargs)

    def accept_residue(self, res):
        """Decide if the residue is valid or not."""
        return True if (res.get_id()[1] in self.entries) else False

    def accept_atom(self, atom):
        """Decide if the atom is valid or not."""
        return super().accept_atom(atom) and self.accept_residue(atom.get_parent())


class AtomSelector(Selector):
    """
    Select a specific set of atoms.

    Parameters
    ----------
    entries : iterable of :class:`~luna.pdb.core.atom.Atom`
        Atom objects to include in the output.
    """
    
    def __init__(self, entries, **kwargs):
        self.entries = entries
        super().__init__(**kwargs)

    def accept_atom(self, atom):
        return super().accept_atom(atom) and atom in self.entries
