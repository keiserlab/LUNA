from MyBio.PDB.PDBIO import Select


class ResidueSelectorByResSeq(Select):
    def __init__(self, entries):
        self.entries = entries

    def accept_residue(self, res):
        return True if (res.get_id()[1] in self.entries) else False


class ResidueSelector(Select):
    def __init__(self, entries):
        self.entries = entries

    def accept_residue(self, res):
        return True if (res in self.entries) else False


class ChainSelector(Select):
    def __init__(self, entries):
        self.entries = entries

    def accept_chain(self, chain):
        return True if (chain in self.entries) else False
