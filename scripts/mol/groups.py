
class AtomGroup():

    def __init__(self, atoms, chemicalFeatures):
        import interaction.util as iu

        self.atoms = atoms
        self.chemicalFeatures = chemicalFeatures
        self.coords = iu.get_coordinatesnp(atoms)
        self._centroid = iu.calc_centroidnp(self.coords)
        self._normal = None
        self._isNormalCalculated = False

    def has_atom(self, atom):
        return atom in self.atoms

    def get_serial_numbers(self):
        return [a.get_serial_number() for a in self.atoms]

    @property
    def compound(self):
        return self.atoms[0].get_parent()

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

    def __repr__(self):
        return '<AtomGroup: [%s]' % ', '.join([str(x) for x in self.atoms])


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
    proximal = get_contacts_for_entity(entity=model, source=myBioPDBResidue,
                                       radius=COV_CUTOFF, level='R')

    covResidues = set()
    for pair in proximal:
        if (pair[1] != myBioPDBResidue):
            covResidues.add(pair[1])

    resSel = ResidueSelector(list(covResidues) + [myBioPDBResidue])

    fh = StringIO()
    io = PDBIO()
    io.set_structure(model)
    io.save(fh, select=resSel, preserve_atom_numbering=True,
            write_conects=True)
    fh.seek(0)
    pdbBlock = ''.join(fh.readlines())

    obMol = readstring("pdb", pdbBlock)
    if (ph is not None):
        obMol.OBMol.AddHydrogens(True, True, ph)

    atomMap = {}
    hydrogCoords = {}

    filteredOBAtom = [a for a in OBMolAtomIter(obMol.OBMol)
                      if a.GetAtomicNum() != 1]

    isWater = myBioPDBResidue.get_id()[0] == "W"
    for idx, obAtom in enumerate(filteredOBAtom):
        # Ignore hydrogen atoms
        if (obAtom.GetAtomicNum() != 1):
            obResidue = obAtom.GetResidue()
            serialNumber = obResidue.GetSerialNum(obAtom)
            atomMap[idx] = serialNumber

            nbCoordinates = []
            if (isWater is False):
                for nbObAtom in OBAtomAtomIter(obAtom):
                    coords = Coordinate(nbObAtom.GetX(),
                                        nbObAtom.GetY(),
                                        nbObAtom.GetZ(),
                                        atomicNumber=nbObAtom.GetAtomicNum())

                    nbCoordinates.append(coords)

            hydrogCoords[serialNumber] \
                = NeighbourhoodCoordinates(nbCoordinates)

    molBlock = obMol.write('mol')

    rdMol = MolFromMolBlock(molBlock)

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

        compoundGroups.add_group(atomGroup)

    return compoundGroups
