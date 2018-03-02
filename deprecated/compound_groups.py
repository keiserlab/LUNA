
class CompoundGroups:

    def __init__(self, myBioPDBResidue, atomGroups=[]):

        self.residue = myBioPDBResidue
        self.atomGroups = atomGroups

    def add_group(self, group):
        self.atomGroups += [group]
