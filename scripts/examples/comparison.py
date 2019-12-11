from rdkit import Chem
from rdkit.Chem import ChemicalFeatures

from mol.chemical_feature import FeatureExtractor
from os.path import (basename, splitext)

import glob
import json

import pprint
pp = pprint.PrettyPrinter(indent=4)

pharm_symbols = {
    "Donor": "D",
    "WeakDonor": "d",
    "HalogenDonor": "x",
    "ChalcogenDonor": "y",
    "Acceptor": "A",
    "WeakAcceptor": "a",
    "PositiveIonizable": "+",
    "NegativeIonizable": "-",
    "Nucleophile": "n",
    "Electrophile": "e",
    "Aromatic": "r",
    "Hydrophobic": "h"
}

atom_names_map = "../tmp/pharma_rules/atom_name_map.json"
atom_names = None
with open(atom_names_map, "r") as IN:
    atom_names = json.load(IN)

luna_fdef = '../data/LUNA.fdef'
luna_bff = ChemicalFeatures.BuildFeatureFactory(luna_fdef)
luna_extractor = FeatureExtractor(luna_bff)

feature_extractor = luna_extractor

working_path = "../tmp/pharma_rules"

ligands_filename = '%s/amino_list' % working_path
ligands_list = set()
with open(ligands_filename, "r") as IN:
    for line in IN.readlines():
        ligands_list.add(line.strip())

sdf_path = '%s/ligand' % working_path

sdf_files = ["%s/%s_ideal.sdf" % (sdf_path, x) for x in sorted(ligands_list)]

mols = [Chem.MolFromMolFile(f) for f in sdf_files]

for rdmol in mols:
    mol_name = rdmol.GetProp("_Name")
    atom_names_by_id = atom_names[mol_name]
    print(mol_name)

    features_by_atom = feature_extractor.get_features_by_atoms(rdmol)
    all_features = []
    for (atm, cf) in sorted(features_by_atom.items()):
        features = [pharm_symbols[c.name] for c in cf]
        all_features.append("/".join(features))

    print(";".join(all_features))
    print()
    output_file = "%s/results/%s.RDKIT.pmapper" % (working_path, mol_name)
    with open(output_file, "w") as OUT:
        OUT.write(";".join(all_features))

print("Done!!!")
