from rdkit import Chem
from rdkit import DataStructs
from mol.fingerprint import FingerprintGenerator

path = '../tmp/rdkit'

infile1 = "%s/1BGA.7DG.ph7.mol" % path
infile2 = "%s/1BGA.7DG.ph14.mol" % path

infile1 = "%s/1AC.ph7.mol" % path
infile2 = "%s/1AC.ph14.mol" % path

m1 = Chem.MolFromMolFile(infile1)
m2 = Chem.MolFromMolFile(infile2)
# m1 = Chem.MolFromSmiles('c1ccccn1')
# m2 = Chem.MolFromSmiles('c1ccco1')

fpGen1 = FingerprintGenerator(m1)
fpGen2 = FingerprintGenerator(m2)

# TOPOLOGICAL FP
fp1 = fpGen1.topological_fp()


# MACCS Keys FP
fp1 = fpGen1.maccs_keys_fp()
fp2 = fpGen2.maccs_keys_fp()

# Atom pairs FP
fp1 = fpGen1.atom_pairs_fp()
fp2 = fpGen2.atom_pairs_fp()

# Atom pairs FP
fp1 = fpGen1.torsion_fp()
fp2 = fpGen2.torsion_fp()

# Morgan FP / Circular FP
fp1 = fpGen1.morgan_fp(2, nBits=2048, features=False, type=3)
fp2 = fpGen2.morgan_fp(2, nBits=2048, features=False, type=3)

# print(fp1.ToBitString())
# print(fp2.ToBitString())

# 2D Pharmacophre FP
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D.SigFactory import SigFactory

fdef = None
if (fdef is None):
    fdefName = 'data/MinimalFeatures.fdef'

featFactory = ChemicalFeatures.BuildFeatureFactory(fdefName)
sigFactory = SigFactory(featFactory, minPointCount=2,
                        maxPointCount=3)
sigFactory.SetBins([(0, 2), (2, 5), (5, 8)])
sigFactory.Init()

fp1 = fpGen1.pharm2d_fp(sigFactory)
fp2 = fpGen2.pharm2d_fp(sigFactory)

# SIMILARITY
print(DataStructs.DiceSimilarity(fp1, fp2))
