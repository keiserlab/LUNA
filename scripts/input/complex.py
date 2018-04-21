

class Complex():

    def __init__(self, pdb, chain, ligName=None, ligNumber=None):

        self.pdb = pdb
        self.chain = chain
        self.ligName = ligName
        self.ligNumber = ligNumber

    def to_string(self, sep=":"):
        entry = [self.pdb, self.chain]

        if self.ligName:
            entry.append(str(self.ligName))
        if self.ligNumber:
            entry.append(str(self.ligNumber))

        return sep.join(entry)


class DBComplex(Complex):

    def __init__(self, complexId, pdb, chain, ligName=None, ligNumber=None):

        self.id = complexId

        super().__init__(pdb, chain, ligName, ligNumber)
