
from rdkit.Chem import ChemicalFeatures
from MyBio.selector import ResidueSelectorByResSeq
from MyBio.PDB.PDBParser import PDBParser
from mol.chemical_feature import FeatureExtractor
from mol.groups import find_compound_groups
from mol.obabel import convert_molecule

from interaction.calc_interactions import calc_all_interactions

from os.path import exists

import pprint

pp = pprint.PrettyPrinter(indent=4)

fdefName = '../data/BaseFeatures.fdef'
featFactory = ChemicalFeatures.BuildFeatureFactory(fdefName)
fExt = FeatureExtractor(featFactory)

pdbPath = "../tmp/pharm"

pdbId = "3QQK"
# pdbId = "1IL3"
pdbFile = "%s/%s.pdb" % (pdbPath, pdbId)

# FETCHED WITH BIOPYTHON
# pdbId = "pdb1eve"
# pdbId = "pdb1p5e"
# pdbId = "pdb3pfp"
# pdbId = "pdb1wbg"
# pdbId = "pdb1r9o"
# pdbFile = "%s/%s.ent" % (pdbPath, pdbId)

newPdbFile = "%s/%s.chainA.ph7.pdb" % (pdbPath, pdbId)
# opt = {"h": None, "error-level": 5}
opt = {"p": 7, "error-level": 5}

if exists(newPdbFile) is False:
    print("It will convert...")
    convert_molecule(pdbFile, newPdbFile, opt=opt)

parser = PDBParser(PERMISSIVE=True)
parser = PDBParser(PERMISSIVE=True,
                   QUIET=True,
                   FIX_ATOM_NAME_CONFLICT=True,
                   FIX_OBABEL_FLAGS=True
                   )

struct = parser.get_structure(pdbId, newPdbFile)
chainA = struct[0]["A"]

# targetMolecules = set([497])
# resSel = ResidueSelectorByResSeq(targetMolecules)

# 3QQK
# lig = ("H_X02", 497, " ")
# nearbyResidues = [chainA[81], chainA[83]]
# nearbyResidues = [chainA[134]]

# 1IL3
# lig = ("H_7DG", 301, " ")
# nearbyResidues = [chainA[80], chainA[193], chainA[123]]
# nearbyResidues = [chainA[1]]

# 1EVE
# lig = ("H_E20", 2001, " ")
# nearbyResidues = [chainA[330]]

# 1P5E
# lig = ("H_TBS", 301, " ")
# nearbyResidues = [chainA[81]]
# nearbyResidues = [chainA[6]]

# 3PFP
# It has ligand with phosphate, MG
# lig = ("H_035", 463, " ")
# nearbyResidues = [chainA[81]]

# 1WBG
# It has a halogenated ligand interacting with an aromatic ring
# chainB = struct[0]["B"]
# lig = [chainB[("H_L03", 1248, " ")]]
# nearbyResidues = [chainB[228]]


from interaction.contact import get_contacts_for_entity

# Target ligand
# 3QQK
targetLigand = chainA[("H_X02", 497, " ")]

# 1R9O
# targetLigand = chainA[("H_FLP", 501, " ")]


# All nearby residues and other molecules
# nearbyResidues = get_contacts_for_entity(chainA, targetLigand, level='R')
# compounds = set([x[1] for x in nearbyResidues])

nearbyResidues = [chainA[81], targetLigand]
compounds = nearbyResidues

targetCompoundGroups = {}
for res in compounds:
    compoundGroups = find_compound_groups(res, fExt)
    targetCompoundGroups[res] = compoundGroups

ligandGroups = targetCompoundGroups[targetLigand]

targetGroups = [targetCompoundGroups[x] for x in targetCompoundGroups
                if x.get_id()[0] != " "]
nearbyGroups = [targetCompoundGroups[x] for x in targetCompoundGroups
                if x != targetLigand]

allInteractions = calc_all_interactions(targetGroups, nearbyGroups)

from util import logging_ini
from util.config_parser import Config
from database.loader import *

from database.util import rows_2dict
from database.luna_model import *
from database.interactions import *

from sqlalchemy import inspect

from data.summary import *
from data.frequency import *


iniFile = "../data/.mysql.ini"
config = Config(iniFile)

dbConfig = config.get_section_map("database")


db = DBLoader(**dbConfig)

# map_interactions_tables(db)

im = InteractionManager(db)

# im.insert_interactions(allInteractions)

# groupTypes = count_group_types(ligandGroups)

# interactionTypes = count_interaction_types(allInteractions)

# frequency = calc_residues_frequency(allInteractions)










# db = QueryExecutor(**dbConfig)
# db.start_session()



# db.new_mapper(InterType, "inter_type")
# db.new_mapper(CompoundType, "compound_type")
# db.new_mapper(Interaction, "interaction")
# db.new_mapper(InterDependInter, "inter_depend_inter")
# db.new_mapper(Group, "group")
# db.new_mapper(Atom, "atom")
# db.new_mapper(GroupAtom, "group_atom")
# db.new_mapper(PDB, "pdb")

# usermapper = mapper(User, users, properties={
    # 'emails': relation(emailmapper),  # Here's where the magic happens
# })

# interTypeRows = db.session.query(InterType).all()

# interactionId = {}
# for row in interTypeRows:
    # interactionId[row.type] = row.id

# insert_interactions(allInteractions, db)

# Access columns
# for column in inspect(GroupCompoundType).columns:
    # print(column.foreign_keys)
