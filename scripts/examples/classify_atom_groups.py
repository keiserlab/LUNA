from rdkit import Chem
from rdkit.Chem import ChemicalFeatures

from mol.features import FeatureExtractor
from os.path import (basename, splitext)

import glob
import json

import pprint
pp = pprint.PrettyPrinter(indent=4)

atom_names_map = "../tmp/pharma_rules/atom_name_map.json"
atom_names = None
with open(atom_names_map, "r") as IN:
    atom_names = json.load(IN)

# Amino acids
sdf_path = "../tmp/pharma_rules/amino"
sdf_files = glob.glob('%s/ALA*.sdf' % sdf_path)
# Ligands
sdf_path = "../tmp/pharma_rules/ligand"

lig = "ALA"
sdf_files = glob.glob('%s/%s*.sdf' % (sdf_path, lig))

#mols = [Chem.MolFromMolFile(m) for m in sdf_files]


# Sulfinic (Thiourea dioxide)
# RDKit original: DefineFeature AcidicGroup [C,S](=[O,S,P])-[O;H1,H0&-1]
# LUNA: DefineFeature SulfinicGroup [S](=[O])-[O;H1,H0&-1]
# smiles = 'C(=N)(N)S(=O)O'

# Sulfonic acid (Taurine)
# RDKit original: DefineFeature AcidicGroup [C,S](=[O,S,P])-[O;H1,H0&-1]
#    -> Captura dois grupos de Sulfinic
# LUNA: DefineFeature SulfonicGroup [S](=[O])(=O)[O;H1,H0&-1]
#    -> Captura apenas o grupo sulfonic acid
# smiles = 'NCCS(=O)(=O)O'


# Phosphinic acids (Hypophosphorous acid)
### RDKIT Original: DefineFeature AcidicGroup [C,S](=[O,S,P])-[O;H1,H0&-1]
# Não captura.
### RDKIT DIP2: DefineFeature AcidicGroup [C,S,P](=[O,S])-[O;H1,H0&-1]
# Captura o grupo
# smiles = 'OP=O'
# smiles = 'CC(C)CP(O)(=O)CC(C)C'


# Phosphonic acids (Glyphosate)
### RDKIT Original: DefineFeature AcidicGroup [C,S](=[O,S,P])-[O;H1,H0&-1]
# Não captura.
### RDKIT DIP2: DefineFeature AcidicGroup [C,S,P](=[O,S])-[O;H1,H0&-1]
# Captura dois grupos de Phosphinic acids
# smiles = 'OC(=O)CNCP(=O)(O)O'

# Phosphoric acid
# smiles = '[H]OP(=O)(O[H])O[H]'

# Phosphate examples
# [4-[3-(4-bromophenyl)-3-oxidanylidene-propyl]-6-methyl-5-oxidanyl-pyridin-3-yl]methyl phosphate
# smiles = 'COP([O-])([O-])=O'
# Glycerophosphoric acid
# smiles = 'OCC(O)COP(O)(O)=O'
# MONOISOPROPYLPHOSPHORYLSERINE
# smiles = 'CC(C)O[P@@](O)(=O)OC[C@H](N)C(O)=O'
# ATP
# smiles = 'Nc1ncnc2n(cnc12)[C@@H]1O[C@H](CO[P@](O)(=O)O[P@@](O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]1O'

# Carboxylic acid (Carbonic acid)
# Captura 2 grupos
# smiles = 'OC(=O)O'
# Carboxylate (Acetate)
# smiles = 'CC(=O)[O-]'

# Tetrazole
# 3-(pyrimidin-2-yl)-N-[3-(1H-tetrazol-5-yl)phenyl]benzamide
# smiles = 'O=C(Nc1cccc(c1)-c1nnn[nH]1)c1cccc(c1)-c1ncccn1'

# Hydroxamic acid
# Rhodotorulic acid
# smiles = 'CC(=O)N(CCCC1C(=O)NC(C(=O)N1)CCCN(C(=O)C)O)O'
# smiles = 'CC(=O)N(CCCC1C(=O)NC(C(=O)N1)CCCN(C(=O)C)O)[O-1]'
# Test molecule
# smiles = 'CCCCC(=O)[N-1]O'
# smiles = 'CCCCC(-[O-1])=NO'


# Sulfonamides
# 2,5-Dichlorothiophene-3-Sulfonamide
# smiles = 'C1=C(SC(=C1S(=O)(=O)N)Cl)Cl'
# PDB id: 032
# smiles = 'CCCS(=O)(=O)Nc1ccc(F)c(C(=O)c2c[nH]c3ncc(cc23)-c2ccc(Cl)cc2)c1F'
# PDB id: 03T
# smiles = 'NS(=O)(=O)c1cc2ccccc2s1'
# PDB id: 051
# smiles = 'OC(=O)[C@@H]1CC[C@H](NS(=O)(=O)c2ccc3cc(OCc4ccc(cc4F)C#N)ccc3c2)[C@@H](C1)C(O)=O'

# Should capture this one: it has a tertiary sulfonamide
# smiles = 'CC1=C(N(C(=C1C#N)N(S(=O)(=O)C)S(=O)(=O)C)CC2=CC=CC=C2)C'


# Nitro group
smiles = '[O-][N+](=O)c1cc(ccc1NCc1ccco1)C(F)(F)F'

# smiles = 'NCC'

# Nao pega o N.
# smiles = 'C1CCC(=CC1)N1CCOCC1'

# Piperidine
# smiles = 'C1CCNCC1'

# Cyanamide
# smiles = 'N#CN'

# smiles = 'NC=C'

# smiles = 'NC=C'

# smiles = 'NS=O'
# smiles = 'NC=O'
#smiles = 'Nc1ccccc1'

# Morpholine
#smiles = 'C1COCCN1'

# Amidine
# smiles = 'n1ccccc1N'
# smiles = 'Nc1ncc[nH]1'
#smiles = 'CC(=O)NC=N'
# smiles = 'CC(=C)NC=N'

# Guanidine
# smiles = 'NC(N)=N'
# smiles = 'NC(=N)NC=O'
# smiles = 'OC(=O)N[C@@H]([C@H]1CCNC(=N)N1)C(O)=O'
# Arginine
# smiles = 'N[C@@H](CCCNC(N)=[NH2+])C(O)=O'


# Sulfonium
#smiles = 'O[C@@H]1C[S@+]2CCCC[C@H]2[C@@H]1O'

# 4-aminopyridin
#smiles = 'NC1=CC=NC=C1'
# smiles = 'NCCNc1ccncc1'

# 2-aminopyridin
#smiles = 'NC1=NC=CC=C1'
#smiles = 'Nc1ccccn1'
#smiles = 'Nc1cc(Cn2c(C(O)=O)c(-c3ccc[nH]c3=O)c3cc(ccc23)C(F)(F)F)ccn1'

# Imidazole-like
#smiles = 'Oc1nc2c([nH]1)[nH]c(=O)[nH]c2=O'
#smiles = 'C1CNCC2(C1)CCN(CC2)c1ncnc2[nH]cnc12'
# Imidazole
#smiles = 'Nc1ncc([nH]1)[C@H]1CC[C@@H](CC1)NC(=O)[C@@H]1C=CCn2n1c(=O)n(CCS(=O)(=O)c1ccccc1)c2=O'
#smiles = 'NC(=O)c1ccc(cc1)-c1nc(c([nH]1)-c1ccc2OCOc2c1)-c1ccccn1'
#smiles = 'NC(=O)c1ccc(cc1)-c1nc(c([nH]1)-c1ccc2OCOc2c1)O'
#smiles = 'Cc1c[nH]c(=O)[nH]1'

# Hydroxamic acid
# Rhodotorulic acid
#smiles = 'CC(=O)N(CCCC1C(=O)NC(C(=O)N1)CCCN(C(=O)C)O)O'
# smiles = 'CC(=O)N(CCCC1C(=O)NC(C(=O)N1)CCCN(C(=O)C)O)[O-1]'
# Test molecule
# smiles = 'CCCCC(=O)[N-1]O'
# smiles = 'CCCCC(-[O-1])=NO'

#smiles = 'CCCN(O)C(C)=O'


#################
## Donor test  ##
##             ##
#Pyrazole
# smiles = 'N1C=CC=N1'

# Glycine
# smiles = 'NCC(O)=O'
# smiles = 'NCC([O-])=O'

# Tertiary amine
# smiles = 'CN(C=O)C(C)=C'
# smiles = 'CN(C=N)C(C)=C'
# smiles = 'CN(C=C)C(C)=C'

# Histidine
# smiles = 'N[C@@H](Cc1c[nH]c[nH+]1)C(O)=O'

# Tetrazole
# smiles = 'c1nnn[nH]1'
# 3-(pyrimidin-2-yl)-N-[3-(1H-tetrazol-5-yl)phenyl]benzamide
# smiles = 'O=C(Nc1cccc(c1)-c1nnn[nH]1)c1cccc(c1)-c1ncccn1'

# Sulfonamides
# 2,5-Dichlorothiophene-3-Sulfonamide
#smiles = 'C1=C(SC(=C1S(=O)(=O)N)Cl)Cl'
# PDB id: 032
#smiles = 'CCCS(=O)(=O)Nc1ccc(F)c(C(=O)c2c[nH]c3ncc(cc23)-c2ccc(Cl)cc2)c1F'
# PDB id: 03T
#smiles = 'NS(=O)(=O)c1cc2ccccc2s1'
# PDB id: 051
#smiles = 'OC(=O)[C@@H]1CC[C@H](NS(=O)(=O)c2ccc3cc(OCc4ccc(cc4F)C#N)ccc3c2)[C@@H](C1)C(O)=O'
# Should capture this one: it has a tertiary sulfonamide
#smiles = 'CC1=C(N(C(=C1C#N)N(S(=O)(=O)C)S(=O)(=O)C)CC2=CC=CC=C2)C'


# Acyl sulfonamides
# PDB id: 08F
#smiles = 'CS(=O)(=O)Nc1cccc(c1)S(=O)(=O)NC(=O)c1c(-c2ccc[nH]c2=O)c2cc(Cl)ccc2n1Cc1ccnc(N)c1'
# PDB id: 0C1
#smiles = 'CS(=O)(=O)NC(=O)c1c(-c2ccc[nH]c2=O)c2c3OCCc3ccc2n1Cc1cc(F)ccc1F'
# Aryl sulfonamides
# Phenyl bounding to the S
# smiles = 'O=S(=O)(N)c1c(Cl)cc2c(c1)S(=O)(=O)NCN2'
# Phenyl bounding to the N
# smiles = 'S(=O)(=O)Nc1c(Cl)cc2c(c1)S(=O)(=O)NCN2'

# Barbituric acid
# smiles = 'O=C1NC(=O)C(N2CCN(CC2)c2ncccn2)(C(=O)N1)c1ccc(Oc2ccccc2)cc1'
# 1,3-diazinane-2,4,6-trione
# smiles = 'O=C1CC(=O)NC(=O)N1'

# Thiazolidinediones and rhodanines
#
# Rhodanine
# smiles = 'O=C1NC(=S)SC1'
# Rosiglitazone
# smiles = 'O=C1NC(=O)SC1Cc3ccc(OCCN(c2ncccc2)C)cc3'

# Phosphoramide
# smiles = 'NP(=O)(N)N'
# smiles = 'CN(C)CN(C)P(N)(N)=O'
# Simple tertiary amide
# smiles = 'N(C)(C)C(=O)'

# Allow double bonded S in acidic groups like carboxylic acid, solfonic acid, etc
# smiles = 'OC(=O)CNCP(O)(O)=S'
# smiles = 'OC(=O)CNCP(O)(O)=S'
# smiles = 'NC(=N)S(O)=S'
# smiles = 'OC(=O)CCC(O)=S'
# smiles = 'OP=S'
# smiles = 'OC(=O)CNCP(=S)(O)O'
# smiles = '[H]OP(=S)(O[H])O[H]'
# smiles = 'COP([O-])([O-])=S'
# smiles = 'OCC(O)COP(O)(O)=S'
# smiles = 'CC(C)O[P@@](O)(=S)OC[C@H](N)C(O)=S'
# smiles = 'Nc1ncnc2n(cnc12)[C@@H]1O[C@H](CO[P@](O)(=O)O[P@@](O)(=O)OP(O)(O)=S)[C@@H](O)[C@H]1O'
# smiles = 'C(=N)(N)S(=S)O'
# smiles = 'NCCS(=O)(=S)O'
# smiles = 'O=S(=S)(O)O'

# Doador - excluir O[C,S,P][=O,S]
# smiles = 'OC(=O)CCC(O)=S'
# smiles = 'OC(=O)CCC(O)=P'
# smiles = 'OC(=S)CCC(O)=P'
# smiles = 'OS(=S)CCC(O)=P'
# smiles = 'OS(=O)CCC(O)=P'
# smiles = 'OP(=O)CCC(O)=P'
# smiles = 'OP(=S)CCC(O)=P'


# Sulfate
# smiles = '[O-]S(=O)(=O)[O-]'
# smiles = 'OS(=O)(=O)Oc1ccc(cc1)[N+]([O-])=O'
# smiles = '[O-]P(=O)([O-])OP(=O)([O-])[O-]'

################################################
#
# ACCEPTOR test
#
################################################

# Imidazole
# smiles = 'N1C=CN=C1'

# Indazole
# smiles = 'N1N=CC2=CC=CC=C12'

# Tetrazole
# smiles = 'c1nnn[nH]1'
# smiles = 'O=C(Nc1cccc(c1)-c1nnn[nH]1)c1cccc(c1)-c1ncccn1'

# Purine
# smiles = 'N1C=NC2=NC=NC=C12'

# Pyrrole
# smiles = 'N1C=CC=C1'

# Pyridine
# smiles = 'C1=CC=NC=C1'

# Methazon
# smiles = 'CN=O'

# Aniline
# smiles = 'NC1=CC=CC=C1'

# Guanidine vs Amidine.
# E.g. of a tertiary nitrogen from a guanidine that is captured by the AmidineLikeN rule.
# Such rule not allow the N to be considered a donor atom.
# smiles = 'NC(=N)N1CCC[C@@H](C1)C(O)=O'

# 2-aminopyridin.
# smiles = 'NC1=NC=CC=C1'

# Nao pegou como amidine...
# smiles = 'C1CN=C2CCCCCN2C1'

# PDB id: DH4
# Nao pegou como amidine...
# smiles = "CC1(C)S[C@H](N[C@H]1C(O)=O)[C@@H](C=O)\\N=C\\N1CCCCCC1"

# N,N,N'-trimethylethanimidamide
# Nao pegou como amidine...
# smiles = 'CN=C(C)N(C)C'

# Pyrazine
# smiles = 'C1=CN=CC=N1'
# Pyrazin-2-amine
# It has an amidine-like substructure
# smiles = 'NC1=NC=CN=C1'
# 1H-pyrrol-2-amine
# It has an amidine-like substructure
# smiles = 'NC1=CC=CN1'


##########################
#
#  Dipole-dipole tests
#
##########################

#
# Nucleophiles
#

# Fluorine
# PDB: 1P2A
# Lig: 5BN
# smiles = 'NCCNc1cc(-c2ccc[nH]2)c2C(=O)Nc3ccc(F)c1c23'

# Chlorine
# PDB: 7STD
# Lig: CRP
# smiles = 'CC[C@]1([C@@H](C)C1(Cl)Cl)C(=O)N[C@H](C)c1ccc(Cl)cc1'

# Bromine
# PDB: 4STD
# Lig: BFS
# smiles = 'C[C@@H](NC(=O)c1cc(F)ccc1O)c1ccc(Br)cc1'

# Iodine
# PDB: 1THA
# Lig: T33
# smiles = 'N[C@@H](Cc1ccc(Oc2ccc(O)c(I)c2)c(I)c1)C(O)=O'

# Carbonyl: C2C=O
# PDB: 1LHV
# Lig: NOG
# smiles = 'CC[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@H]34)[C@@H]1CC[C@@]2(O)C#C'

# Carbonyl: R2C=O
# Lig: 0A0
# smiles = 'C[C@](N)(CC(O)=O)C(O)=O'

# Nitrile/Cyano
# Lig: 03U
# smiles = "CC\\C(O)=C(/C#N)C(=O)Nc1ccc(-c2ccc(F)cc2)c(c1)C(=O)OC"

# Nitro
# PDB: 2ROY
# Lig: P28
# smiles = 'CC(=O)N[C@@H](Cc1ccc(Oc2cc(c(O)c(c2)[N+]([O-])=O)[N+]([O-])=O)cc1)C(O)=O'

# R2NO: N-O
# Lig: 0N7
# smiles = 'ON1C(=O)Cc2ccccc2C1=O'
# R2NO: N=O
# Lig: R1A
# smiles = 'CC1(C)C=C(CSSC[C@H](N)C(O)=O)C(C)(C)[N+]1=O'
# R2NO: N-O + N=O
# Pubchem CID: 28467
# smiles = 'C1=CC=C2C(=C1)N(C(=C([N+]2=O)CO)CO)[O-]'
# R2NO: N-O
# Pubchem CID: 15602694
# smiles = 'CC1=C(C(=[N+](C(=C1Cl)C)[O-])Cl)C2=NC(=C3C=C(C(=O)C(=C3)O)[N+](=O)[O-])ON2'
# ligand: 129
# smiles = 'ON(CCP(O)(O)=O)C=O'
# ligand: 0N7
# smiles = 'ON1C(=O)Cc2ccccc2C1=O'
# ligand: AL0
# smiles = 'N[C@@H](CN(O)N=O)C(O)=O'

# Lig: R3N
# smiles = 'COc1cc(cc(OC)c1OC)C(F)(F)C(=O)N1CCCC[C@H]1C(=O)O[C@@H](CCCc1ccccc1)CCCc1cccnc1'
# smiles = 'CN(O)C'

# Sulfonyl-like
# PDB: 1IDB
# Lig: 0DO
# smiles = 'Cc1cccc(C)c1OCC(=O)N[C@@H](Cc1ccccc1)[C@H](O)CN1CC[C@H](C[C@H]1C(=O)NC(C)(C)C)S(=O)(=O)c1ccncc1'

# Sulfinyl-like
# Lig: 3TQ
# smiles = 'CC(C)CN1CCC(CC1)[S@](=O)c1ccc(CNC(=O)c2cc3ccncc3o2)cc1'
# Lig: DSS
# smiles = 'CSC[S@](C)=O'

# DMSO
# smiles = 'CS(C)=O'

# Imine
# PDB: 1NJF, Lig: ADP
# smiles = 'Nc1ncnc2n(cnc12)[C@@H]1O[C@H](CO[P@@](O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]1O'
# methyl(propan-2-ylidene)amine
# smiles = 'CN=C(C)C'

# Ether
# PDB: 1NM9
# Lig: AGL
# smiles = 'C[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1N'

# Furan
# smiles = 'c1ccoc1'

# Thioether
# PDB: 1JFH, Lig: MA2
# smiles = 'CS[C@@H]1[C@@H](CO)O[C@H](O)[C@H](O)[C@H]1O'
# Thiophene
# smiles = 'c1ccsc1'
# Benzenethiophene
# smiles = 's2c1ccccc1cc2'
# Tetrahydrothiophene
# smiles = 'S1CCCC1'

# Alcohol
# PDB: 1MRW
# Lig: K57
# smiles = 'Cc1c(O)cccc1C(=O)N[C@@H](Cc1ccccc1)[C@H](O)C(=O)N1CSC(C)(C)[C@H]1C(=O)NC(C)(C)C'

# Alcohol in Carboxylic acid
# Lig: 0A0
# smiles = 'C[C@](N)(CC(O)=O)C(O)=O'

# Water
# smiles = '[OH2]'

#############################

#
# Nucleophiles
#

# Amidine
# smiles = 'CC(=C)NC=N'
# smiles = 'CC(N)=N'

# Guanidine
# smiles = 'OC(=O)N[C@@H]([C@H]1CCNC(=N)N1)C(O)=O'

# Nitro
# PDB: 2ROY
# Lig: P28
# smiles = 'CC(=O)N[C@@H](Cc1ccc(Oc2cc(c(O)c(c2)[N+]([O-])=O)[N+]([O-])=O)cc1)C(O)=O'


# Fluorine
# PDB: 1P2A
# Lig: 5BN
# smiles = 'NCCNc1cc(-c2ccc[nH]2)c2C(=O)Nc3ccc(F)c1c23'


# Chalcogen donor

# Thiazole
# smiles = 'n1ccsc1'

# 1,3,4-thiadiazole
# smiles = 'S1C=NN=C1'

# 1,2,5-thiadiazole
# No chalcogen donor
# smiles = 'S1N=CC=N1'

# Methionine
# smiles = 'CSCCC(N)C(O)=O'

# Cystheine
# No chalcogen donor
# smiles = 'NC(CS)C(O)=O'

# C-S
# No chalcogen donor
# smiles = 'CS'

# 2-amino-3-[(2-aminoethyl)disulfanyl]propanoic acid
# smiles = 'NCCSSCC(N)C(O)=O'



# Aromatic rings

# TRYPTOPHAN
# smiles = 'N[C@@H](Cc1c[nH]c2ccccc12)C(O)=O'

# Hydrophobic
# smiles = 'Nc1ncnc(n1)-c1cc(ccc1O)C#N'


#
# Ketene acetal
# Target: ligand TEO from PDB.
# This ligand has the substructure 'OC(O)(=*)' which is negatively ionizable.
#
# smiles = 'O[C@H](\C=C(/O)[O-])C([O-])=O'
# smiles = 'OC1=CCC2=C(O1)C=CC=C2'
# smiles = 'OC1=CCC=CO1'
# smiles = '[O-1]C1=CCC=CO1'
# smiles = 'OC(O)=C'
# smiles = 'OC([O-1])=C'
# smiles = 'CCC1=CC(=O)C(C)=C(O)O1'
# smiles = 'CCC1=CC(=O)C(C)=C([O-1])O1'
# smiles = 'OC1=CC=CO1'
#
# Não deve capturar a estrutura nesses exemplos
#
# smiles = 'CCC1=CC(=O)C(C)=C(OC)O1'
# smiles = 'COC(OC)=C[C@@H](O)C([O-])=O'
# smiles = 'COC1=CCC=CO1'
# smiles = 'OC1=CCC2=C(O1)C=CC=C2'


# Charged atoms test

# Nitrate

# Azide
# Azidomethane
# smiles = 'CN=[N+]=[N-]'
# DG9 simplified
# smiles = 'CC(CO[P+]([O-])=O)[C@@H](O)C(N)=O'
# [(Octan-2-yl)oxy](oxo)phosphaniumolate
# smiles = 'CCCCCCC(C)O[P+](=O)[O-]'
# (phosphooxy)methane
# smiles = 'CO[P+]([O-])=O'
# phenyl(phosphooxy)phosphinate;zirconium(4+)
# smiles = 'C1=CC=C(C=C1)P(=O)([O-])O[P+](=O)[O-].C1=CC=C(C=C1)P(=O)([O-])O[P+](=O)[O-].[Zr+4]'


# METALS

# TODO: Some of the metalic complex does not work in RDKit.

# HEM
# smiles = 'CC1=C(CCC(O)=O)C2=Cc3c(CCC(O)=O)c(C)c4C=C5C(C)=C(C=C)C6=[N]5[Fe]5([N]2=C1C=c1c(C=C)c(C)c(=C6)n51)n34'
# smiles = 'CC1=C(C2=CC3=C(C(=C([N-]3)C=C4C(=C(C(=N4)C=C5C(=C(C(=N5)C=C1[N-]2)C=C)C)C=C)C)C)CCC(=O)O)CCC(=O)O.[Fe+2]'
# smiles = 'CC1=C(C2=CC3=NC(=CC4=C(C(=C([N-]4)C=C5C(=C(C(=N5)C=C1[N-]2)C)CCC(=O)O)CCC(=O)O)C)C(=C3C)C=C)C=C.[Fe]'
# smiles = 'CC1=C(C2=CC3=NC(=CC4=C(C(=C([N-]4)C=C5C(=C(C(=N5)C=C1[N-]2)C(CCC=C(C)CCC=C(C)CCC=C(C)C)O)C)C=O)CCC(=O)O)C(=C3C)CCC(=O)O)C=C.[Fe+2]'

# Cobalt
# smiles = 'C[C@H](CNC(=O)CC[C@]1(C)[C@@H](CC(N)=O)[C@H]2[N]3=C1C(C)=C1[C@@H](CCC(N)=O)C(C)(C)C4=CC5=[N]6C(=C(C)C7=[N]([C@]2(C)[C@@](C)(CC(N)=O)[C@@H]7CCC(N)=O)[Co++]36N14)[C@@](C)(CC(N)=O)[C@@H]5CCC(N)=O)O[P@@](O)(=O)O[C@@H]1[C@@H](CO)O[C@@H]([C@@H]1O)n1cnc2cc(C)c(C)cc12'
# smiles = 'CC1=C2[C@@H](CCC(N)=O)[C@](C)(CC(N)=O)[C@]3(C)[C@H]4[C@H](CC(N)=O)[C@@](C)(CCC(=O)NCCO[P@](O)(=O)O[C@@H]5[C@@H](CO)O[C@@H]([C@@H]5O)n5cnc6ccc(O)cc56)C5=[N+]4[Co@@]4(N23)[N+]2=C1[C@@](C)(CC(N)=O)[C@H](CCC(N)=O)C2=CC1=[N+]4C([C@@H](CCC(N)=O)C1(C)C)=C5C'
# smiles = '[Co]123456789[B@]%10%11[B]%12%13%14[B]%15%16([B]%17%18[B]%19%20[B]%12%15%17[B]%10%13%19[C]1%11%20[C]2%18(B3%16)CNS(N)(=O)=O)B4%14.[B]123[B]4%10[B]11[B]44[B]%10%10(B5)B6[C@@]74[C]821B9B3%10'
# smiles = '[B]1234[B@]56[Co]789%10%11([B@@]%12%13[B]77%14[B]%15%16%17[B]%18%19%20B%12[B]%137%15%18)(C[B]8%14%16[B]%17%19B%20)([C]78%12[B]%13%14%15[B]17([B]2%131[B]35([B]%141[B@@]98%15)B%10)[C]46%11%12CCNS(=O)(=O)N)C'
# smiles = 'C[C@@H](O)CNC(=O)CC[C@]1(C)[C@@H](CC(N)=O)C2[N-]3C1=C(C)C1=[N]4C(=CC5=[N]6C(=C(C)C7=[N]([C@]2(C)[C@@](C)(CC(N)=O)[C@@H]7CCC(N)=O)[Co++]346)[C@@](C)(CC(N)=O)[C@@H]5CCC(N)=O)C(C)(C)[C@@H]1CCC(N)=O'
# smiles = '[H][N]([H])([H])[Co+3]([N]([H])([H])[H])([N]([H])([H])[H])[N]([H])([H])[H]'

# Be
# smiles = 'F[Be-](F)F'

# Mg
# smiles = 'CCC1=C(C)C2=[N-]3C1=Cc1c(C)c([C@H](C)O)c4C=C5C(C)=C6C(=O)[C@H](C(=O)OC)C7=C8[C@@H](CCC(=O)OC\C=C(/C)CC\C=C(/C)CCC=C(C)C)[C@H](C)C(=C2)N8[Mg@]3(n14)[N-]5=C67'



# Checking ligands with error:

# 09N
# smiles = 'Cn1cnc(c1)-c1nc(C(O)=O)c(O)c(=O)[nH]1'

# 4PO
# smiles = 'CC(C)(C)[N+](\[O-])=C\c1cc[n+]([O-])cc1'


# ALA
# smiles = 'C[C@H](N)C(O)=O'

# ASP
# smiles = 'N[C@@H](CC(N)=O)C(O)=O'

# smiles = 'C#N'


# Guanidine-like
# smiles = 'NC(N)=N'
# smiles = 'CCNC(N)=N'
# # Guanidine-like tautomer
# smiles = 'CCN=C(N)N'
# smiles = 'CC[NH+]=C(N)N'


# Aromatic diformamide substructures.
# Methyluracil
# smiles = 'Cn1ccc(=O)[nH]c1=O'
# # Uracil
# smiles = 'O=C1NC=CC(=O)N1'



# # Chalcogen donors analysis

# # Ligand id: MET
# smiles = 'CSCC[C@H](N)C(O)=O'

# # Ligand id: 798
# smiles = 'O=C(c1ccc(OC[C@H]2CCCN2)cc1)c1ccc(cc1)-c1ccsc1'

# # Ligand id: MCV
# smiles = 'COc1ccc(OC)c(CCc2csc3nc(N)nc(N)c23)c1'

# # Ligand id: BLT
# smiles = 'OC[C@@H](OS([O-])([O-])[O-])[C@@H](O)C[Se+]1C[C@@H](O)[C@H](O)[C@H]1CO'

# # Ligand id: 30V
# smiles = 'N[C@@H](CS[Se]c1ccccc1C(N)=O)C(O)=O'

# # Ligand id: TPP
# smiles = 'Cc1c(CCO[P@@](O)(=O)OP(O)(O)=O)sc[n+]1Cc1cnc(C)nc1N'

# # Ligand id: MTD
# # smiles = 'C[Te]CC([O-])=O'

# # Not contain chalcogen donor
# # Ligand id: CYS
# # smiles = 'N[C@@H](CS)C(O)=O'

# #
# # Isothiazoles: Chalcogen donors
# #
# smiles = 'S1C=CC=N1'

# # Ligand id: 49J
# smiles = 'COc1cc(ccc1N)-c1cnc2c(snc2c1)N1CCOCC1'

# # Ligand id: NF6
# smiles = 'CN(C)c1ccc2C(=O)N(CCNCc3nnc(C)s3)C(=O)c3cccc1c23'

# # Ligand id: 69V
# smiles = 'O=C(Cc1ccccc1)Nc1nnc(OC2CCN(CC2)c2nnc(NC(=O)Cc3ccccc3)s2)s1'

# # Ligand id: A65
# smiles = 'CSc1nsc(SC)c1CO'

# # Ligand id: LTI
# smiles = 'NC(=O)\\N=c1/s[nH]c(S[C@@H]2CCCCc3ccccc23)c1C(N)=O'

# # Not contain chalcogen donor
# # Ligand id: GBD
# # smiles = 'OC(=O)c1nsnc1O'

# smiles = 'CSc1cc([nH]n1)-c1sc(nc1C)-c1cccs1'

# smiles = 'N[C@@H](CS)C(O)=O'



# Hydrophobic analysis

# TRP
smiles = 'N[C@@H](Cc1c[nH]c2ccccc12)C(O)=O'
# TYR
smiles = 'N[C@@H](Cc1ccc(O)cc1)C(O)=O'
# PHE
smiles = 'N[C@@H](Cc1ccccc1)C(O)=O'

# Fluorine
# PDB: 1P2A, Lig: 5BN
smiles = 'NCCNc1cc(-c2ccc[nH]2)c2C(=O)Nc3ccc(F)c1c23'

# # Chlorine
# # PDB: 7STD, Lig: CRP
# smiles = 'CC[C@]1([C@@H](C)C1(Cl)Cl)C(=O)N[C@H](C)c1ccc(Cl)cc1'

# # Bromine
# # PDB: 4STD, Lig: BFS
# smiles = 'C[C@@H](NC(=O)c1cc(F)ccc1O)c1ccc(Br)cc1'

# # Iodine
# # PDB: 1THA, Lig: T33
# smiles = 'N[C@@H](Cc1ccc(Oc2ccc(O)c(I)c2)c(I)c1)C(O)=O'

smiles = 'Nc1ncnc(n1)-c1cc(ccc1O)C#N'
smiles = 'C1=CC2=CC5=CC=C(C=C4C=CC(C=C3C=CC(=CC1=N2)N3)=N4)N5'
smiles = 'C1CCC2CCCCC2C1'
smiles = 'N1NNN2NNNNN2N1'
smiles = 'C1CCCNCC1'

# 9-membered aromatic ring
# smiles = 'C1=CC=CNC=CC=C1'
# 8-membered aromatic ring
smiles = '[CH-]1\\C=C/[CH-]\\C=C/C=C\\1'

# 7-membered ring
smiles = 'O1C=CC=CC=C1'
smiles = 'C1CCCSCC1'
# smiles = 'C1CCCOCC1'
# smiles = 'C1=CC=COC=C1'
# smiles = 'S1C=CC=CC=C1'
# smiles = 'C1CCCNCC1'
# smiles = 'C1=CC=CNC=C1'

# smiles = 'c1cccccc1'

# smiles = 'C1=CCCC1'
# smiles = 'C1CSCCN1'
# smiles = 'C1CCCCC1'

# smiles = 'N=1\\C=C/C=C\\C=C/C=1'
# smiles = 'c1ccccccc1'
# smiles = 'c1cccnccc1'

# Benzene
# smiles = 'C1=CC=CC=C1'
# Pyridine
# smiles = 'C1=CC=NC=C1'

# Pyrrole
# smiles = 'N1C=CC=C1'
# Pyrazine
# smiles = 'C1=CN=CC=N1'


smiles = 'N=C(n1cccn1)n1cccn1'
smiles = 'NC(=N)n1cccn1'

smiles = 'CCCCCCCCCCCCCC[C@@H](O)[C@@H](O)[C@H](CO[C@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)NC(=O)CCCCCCCCC\\C=C/C\\C=C/CCCCC'

print()
print(smiles)
print()

luna_fdef = '../data/LUNA.fdef'
# luna_fdef = '../data/BaseFeatures_DIP2_NoMicrospecies.fdef'

luna_bff = ChemicalFeatures.BuildFeatureFactory(luna_fdef)
luna_extractor = FeatureExtractor(luna_bff)
feature_extractor = luna_extractor


mols = [Chem.MolFromSmiles(smiles)]

CHEMICAL_FEATURE_IDS = {
    "Aromatic": 1,
    "Acceptor": 2,
    "Donor": 3,
    "Hydrophobe": 4,
    "Hydrophobic": 5,
    "LumpedHydrophobe": 13,
    "Negative": 6,
    "Positive": 7,
    "NegativelyIonizable": 8,
    "PositivelyIonizable": 9,
    "HalogenDonor": 10,
    "HalogenAcceptor": 11,
    "Metal": 12,
    "WeakDonor": 14,
    "WeakAcceptor": 15,
    "Electrophile": 16,
    "Nucleophile": 17,
    "ChalcogenDonor": 19,
    "Amide": 21,
    "Atom": 22
}

for rdmol in mols:
    group_features = feature_extractor.get_features_by_groups(rdmol)
    for grp in group_features:
        atm_ids = group_features[grp]["atm_ids"]
        atom_symbols = [rdmol.GetAtomWithIdx(i).GetSymbol() for i in atm_ids]

        features = []
        for f in group_features[grp]["features"]:
            if f.name in CHEMICAL_FEATURE_IDS:
                features.append(f.name)

        print("%s, %s:\t%s" % (atm_ids, atom_symbols, features))

exit()

for rdmol in mols:
    mol_name = rdmol.GetProp("_Name")
    print(mol_name)
    atom_names_by_id = atom_names[mol_name]

    group_features = feature_extractor.get_features_by_groups(rdmol)
    for grp in group_features:
        atm_ids = group_features[grp]["atm_ids"]
        atm_ids = [atom_names_by_id[str(x)] for x in atm_ids]
        print("%s: %s" % (group_features[grp]["features"], atm_ids))

    exit()