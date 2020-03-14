from luna.MyBio.PDB.PDBIO import Select


DEFAULT_ALTLOC = ("A", "1")


class Selector(Select):

    def __init__(self, keep_hydrog=True, keep_altloc=True, altloc=DEFAULT_ALTLOC):
        self.keep_hydrog = keep_hydrog
        self.keep_altloc = keep_altloc
        self.altloc = altloc

    def accept_atom(self, atom):
        # Hydrogen and Deuterium
        if not self.keep_hydrog and atom.element in ["H", "D"]:
            return False

        if self.keep_altloc:
            return True
        else:
            return not atom.is_disordered() or atom.get_altloc() in self.altloc


class ResidueSelectorByResSeq(Selector):

    def __init__(self, entries, **kwargs):
        self.entries = entries
        super().__init__(**kwargs)

    def accept_residue(self, res):
        return True if (res.get_id()[1] in self.entries) else False


class ChainSelector(Selector):

    def __init__(self, entries, **kwargs):
        self.entries = entries
        super().__init__(**kwargs)

    def accept_chain(self, chain):
        return True if (chain in self.entries) else False


class ResidueSelector(Selector):

    def __init__(self, entries, **kwargs):
        self.entries = entries
        super().__init__(**kwargs)

    def accept_residue(self, res):
        return res in self.entries


class AtomSelector(Selector):
    pass
