from Bio.PDB.Chain import Chain


class myChain(Chain):

    def __lt__(self, r2):
        return self.id < r2.id
