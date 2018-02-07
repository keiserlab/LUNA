import openbabel as ob


# THIS FUNCTION IS PRETTY SLOW...


def extract_residue(obmol, targetIdx):

    copy = ob.OBMol(obmol)

    targetIdx = set(targetIdx)

    atoms = list(ob.OBMolAtomIter(copy))

    for atm in reversed(atoms):
        if (atm.GetResidue().GetIdx() not in targetIdx):
            copy.DeleteAtom(atm)

    return copy
