from interaction.calc_interactions import (DefaultInteractionConf,
                                           InteractionConf)
from mol.depiction import ColorPallete

from MyBio.PDB.PDBParser import PDBParser

ENTRIES_SEPARATOR = ":"

DEFAULT_NAPOLI_PATH = "/media/data/Workspace/nAPOLI_v2/tmp/nAPOLI"
DEFAULT_PDB_PATH = "%s/public/pdb" % DEFAULT_NAPOLI_PATH
DEFAULT_TMP_FILES = "%s/public/tmp" % DEFAULT_NAPOLI_PATH

DEFAULT_DB_CONF_FILE = "../data/.mysql.ini"

DEFAULT_ATOM_PROP_FILE = "../data/Napoli.fdef"

DEFAULT_INTERACTION_CONF = DefaultInteractionConf()

BOUNDARY_CONF = InteractionConf({"boundary_cutoff": 7})

PDB_PARSER = PDBParser(PERMISSIVE=True,
                       QUIET=True,
                       FIX_ATOM_NAME_CONFLICT=True,
                       FIX_OBABEL_FLAGS=True)

NMR_METHODS = ["SOLID-STATE NMR", "SOLUTION NMR"]

DEFAULT_ATOM_TYPES_COLOR = ColorPallete({
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

DEFAULT_CHEMICAL_FEATURES_IDS = {
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
    "Chalcogen donor": 19
}

DEFAULT_INTERACTIONS_IDS = {
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
    "Pi-stacking": 12
}
