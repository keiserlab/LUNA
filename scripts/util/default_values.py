from mol.interaction.conf import (DefaultInteractionConf, InteractionConf)
from mol.depiction import ColorPallete

ENTRIES_SEPARATOR = ":"

NAPOLI_PATH = "/media/data/Workspace/nAPOLI_v2/tmp/nAPOLI"
PDB_PATH = "%s/public/pdb" % NAPOLI_PATH
TMP_FILES = "%s/public/tmp" % NAPOLI_PATH
DB_CONF_FILE = "../data/.mysql.ini"
ATOM_PROP_FILE = "../data/Napoli.fdef"
INTERACTION_CONF = DefaultInteractionConf()

BOUNDARY_CONF = InteractionConf({"boundary_cutoff": 7})

NMR_METHODS = ["SOLID-STATE NMR", "SOLUTION NMR"]

ATOM_TYPES_COLOR = ColorPallete({
    "Acceptor": (252, 141, 89),
    "Donor": (145, 191, 219),
    "Aromatic": (224, 243, 248),
    "Hydrophobic": (254, 224, 144),
    "Hydrophobe": (254, 224, 144),
    "LumpedHydrophobe": (254, 224, 144),
    "PosIonizable": (69, 117, 180),
    "NegIonizable": (215, 48, 3),
    "HalogenAcceptor": (215, 48, 3),
    "HalogenDonor": (215, 48, 3),
})

CHEMICAL_FEATURES_IDS = {
    "Aromatic": 1,
    "Acceptor": 2,
    "Donor": 3,
    "Hydrophobe": 4,
    "Hydrophobic": 5,
    "Negative": 6,
    "Positive": 7,
    "Negatively ionizable": 8,
    "Positively ionizable": 9,
    "Halogen donor": 10,
    "Halogen acceptor": 11,
    "Metal": 12,
    "Lumped hydrophobe": 13,
    "Weak donor": 14,
    "Weak acceptor": 15,
    "Electrophile": 16,
    "Nucleophile": 17,
    "Chalcogen donor": 19,
    "Amide": 21,
    "Atom": 22
}

INTERACTIONS_IDS = {
    "Proximal": 0,
    "Hydrogen bond": 1,
    "Attractive": 2,
    "Salt bridge": 3,
    "Cation-pi": 4,
    "Edge-to-face pi-stacking": 5,
    "Face-to-face pi-stacking": 6,
    "Parallel-displaced pi-stacking": 7,
    "Hydrophobic": 8,
    "Halogen bond": 9,
    "Repulsive": 10,
    "Water-bridged hydrogen bond": 11,
    "Pi-stacking": 12,
    "Amide-aromatic stacking": 13,
    "Weak hydrogen bond": 14,
    "Covalent bond": 15,
    "Atom overlap": 16,
    "Van der Waals clash": 17,
    "Van der Waals": 18
}

IFP_FINGERPRINT_LENGTH = 1024
