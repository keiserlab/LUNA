
class AtomGroup():

    def __init__(self, atoms, chemicalFeatures):
        import interaction.util as iu

        self.atoms = atoms
        self.chemicalFeatures = chemicalFeatures
        self.coords = iu.get_coordinatesnp(atoms)
        self._centroid = iu.calc_centroidnp(self.coords)
        self._normal = None
        self._isNormalCalculated = False

    def get_serial_numbers(self):
        return [a.get_serial_number() for a in self.atoms]

    @property
    def centroid(self):
        return self._centroid

    @property
    def normal(self):
        if self._isNormalCalculated is False:
            self._calc_normal()

        return self._normal

    def _calc_normal(self):
        import interaction.plane_regression as plane

        if self._isNormalCalculated is False:
            self._normal = plane.calc_normal(self.coords)
            self._isNormalCalculated = True


class CompoundGroups():

    def __init__(self, myBioPDBResidue, atomGroups=None):
        self.residue = myBioPDBResidue

        if (atomGroups is None):
            atomGroups = []

        self.atomGroups = atomGroups

    def add_group(self, group):
        self.atomGroups += [group]


def find_compound_groups(myBioPDBResidue, featureExtractor, ph=None,
                         hasExplicitHydrogen=False):

    from MyBio.PDB.PDBIO import PDBIO
    from MyBio.selector import ResidueSelector
    from interaction.contact import get_contacts_for_entity
    from mol.vicinity_atom import VicinityAtom
    from mol.coordinates import (Coordinate, NeighbourhoodCoordinates)

    from rdkit.Chem import MolFromMolBlock
    from rdkit.Chem.rdDepictor import Compute2DCoords

    from openbabel import (OBMolAtomIter, OBAtomAtomIter)

    from pybel import readstring
    from io import StringIO

    COV_CUTOFF = 1.99999999999

    model = myBioPDBResidue.get_parent_by_level('M')
    proximal = get_contacts_for_entity(model=model, source=myBioPDBResidue,
                                       radius=COV_CUTOFF, level='R')

    covResidues = set()
    for pair in proximal:
        if (pair[1] != myBioPDBResidue):
            covResidues.add(pair[1])

    resSel = ResidueSelector(list(covResidues) + [myBioPDBResidue])
    resSel = ResidueSelector([myBioPDBResidue])

    fh = StringIO()
    io = PDBIO()
    io.set_structure(model)
    io.save(fh, select=resSel, preserve_atom_numbering=True,
            write_conects=True)
    fh.seek(0)
    pdbBlock = ''.join(fh.readlines())

    print(pdbBlock)
    print()

    from rdkit import Chem
    from rdkit.Chem import rdmolfiles    
    from mol.charge_model import perceive_formal_charge
    from rdkit.Chem import rdmolops

    rdMol = Chem.MolFromPDBBlock(pdbBlock, removeHs=False)

    # rdMol = Chem.MolFromSmiles("[H]c1[c]c([H])c([H])c1[H]", sanitize=False)
    

    # rdMol = Chem.MolFromPDBBlock(pdbBlock, removeHs=False, sanitize=False)

    rdMol.UpdatePropertyCache(strict=False)
    rdmolops.AddHs(rdMol, explicitOnly=True)

    # add_formal_charges(rdMol)

    atoms = list(myBioPDBResidue.get_atoms())

    for i, rdAtom in enumerate(rdMol.GetAtoms()):
        if (rdAtom.GetAtomicNum() == 1):
            continue

        print(atoms[i], rdAtom.GetAtomicNum())
        print("Degree:", rdAtom.GetDegree())
        print("ExpValence:", rdAtom.GetExplicitValence(), "ImpValence:", rdAtom.GetImplicitValence())
        print("Formal charge:", rdAtom.GetFormalCharge())
        print("Num Hs:", rdAtom.GetTotalNumHs())

        predictedFormalCharge = perceive_formal_charge(rdAtom)

        if (predictedFormalCharge is not None):
            rdAtom.SetFormalCharge(predictedFormalCharge)
            rdAtom.UpdatePropertyCache()

        print("--------------------------")        
        print("Predicted charge:", predictedFormalCharge)
        print("ExpValence:", rdAtom.GetExplicitValence(), "ImpValence:", rdAtom.GetImplicitValence())
        print("Formal charge:", rdAtom.GetFormalCharge())
        print("Num Hs:", rdAtom.GetTotalNumHs())
        print()
        for neighbor in rdAtom.GetNeighbors():
            print(neighbor.GetAtomicNum())

        print()
        print()



    print(rdmolfiles.MolToSmarts(rdMol))


    return

    exit()






    # from Chem import rdPartialCharges
    # rdPartialCharges.ComputeGasteigerCharges(rdMol)

    # https://gist.github.com/greglandrum/7f546d1e35c2df537c68a64d887793b8

    for i, rdAtom in enumerate(rdMol.GetAtoms()):

        if (rdAtom.GetAtomicNum() != 8):
            continue
        print(atoms[i], rdAtom.GetAtomicNum())
        print(rdAtom.GetFormalCharge())
        print(rdAtom.GetDegree())
        print(rdAtom.GetExplicitValence(), rdAtom.GetImplicitValence())
        print("H", rdAtom.GetNumExplicitHs(), rdAtom.GetNumImplicitHs())
        print(rdAtom.GetNoImplicit())

        rdAtom.SetNoImplicit(True)        
        print(rdAtom.GetNoImplicit())
        print("H", rdAtom.GetNumExplicitHs(), rdAtom.GetNumImplicitHs())    
        print(rdAtom.GetFormalCharge())
        print()        

    print(rdmolfiles.MolToSmarts(rdMol))

    exit()

    obMol = readstring("pdb", pdbBlock)
    if (ph is not None):
        obMol.OBMol.AddHydrogens(True, True, ph)

    atomMap = {}
    hydrogCoords = {}

    filteredOBAtom = [a for a in OBMolAtomIter(obMol.OBMol)
                      if a.GetAtomicNum() != 1]

    for idx, obAtom in enumerate(filteredOBAtom):
        # Ignore hydrogen atoms
        if (obAtom.GetAtomicNum() != 1):
            obResidue = obAtom.GetResidue()
            serialNumber = obResidue.GetSerialNum(obAtom)
            atomMap[idx] = serialNumber

            if (obResidue.GetName()=="LYS"):             
                print(list(myBioPDBResidue.get_atoms())[obAtom.GetIdx()-1])
                print(obAtom.GetIdx()-1, obAtom.GetAtomicNum(), obResidue.GetSerialNum(obAtom))
                obAtom.SetImplicitValence(10)
                print(obAtom.GetValence(), obAtom.GetImplicitValence(), obAtom.GetHvyValence(), obAtom.GetHeteroValence())

                # print("Charge: ", obAtom.GetFormalCharge())
                # Add negative charge to 
                # obAtom.SetFormalCharge(-1)
                # print("Charge: ", implicitHs)
                print()

            nbCoordinates = []
            for nbObAtom in OBAtomAtomIter(obAtom):
                # Only hydrogen coordinates
                if nbObAtom.GetAtomicNum() == 1:
                    coords = Coordinate(nbObAtom.GetX(),
                                        nbObAtom.GetY(),
                                        nbObAtom.GetZ(),
                                        atomicNumber=1)

                    nbCoordinates.append(coords)

            hydrogCoords[serialNumber] \
                = NeighbourhoodCoordinates(nbCoordinates)

    molBlock = obMol.write('mol')

    print()
    print(molBlock)
    print()
    print(obMol.write('smi'))
    print()
    print()

    from rdkit.Chem import rdmolfiles

    rdMol = MolFromMolBlock(molBlock)
    print(rdmolfiles.MolToSmarts(rdMol))

    exit()

    groupFeatures = featureExtractor.get_features_by_groups(rdMol, atomMap)

    targetAtoms = {}
    for atom in myBioPDBResidue.get_atoms():
        # Ignore hydrogen atoms
        if (atom.element != "H"):
            nbCoords = hydrogCoords[atom.get_serial_number()]
            vicinityAtom = VicinityAtom(atom, nbCoords)

            targetAtoms[atom.get_serial_number()] = vicinityAtom

    compoundGroups = CompoundGroups(myBioPDBResidue)

    for gKey in groupFeatures:
        groupObj = groupFeatures[gKey]

        isGroupValid = True
        atoms = []
        for atomIdx in groupObj["atomIds"]:
            if (atomIdx not in targetAtoms):
                isGroupValid = False
                break
            else:
                atoms.append(targetAtoms[atomIdx])

        # This group has atoms that do not belong to the target residue.
        if (isGroupValid is False):
            continue

        atomGroup = AtomGroup(atoms, groupObj["features"])

        print(atoms, groupObj["features"])

        compoundGroups.add_group(atomGroup)

    return compoundGroups
