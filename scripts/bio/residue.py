from Bio.PDB.Residue import Residue


class myResidue(Residue):

    def __lt__(self, r2):
        return self.id < r2.id
