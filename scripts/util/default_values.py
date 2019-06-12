from mol.interaction.conf import DefaultInteractionConf, InteractionConf
from util import ColorPallete

from os import path

ENTRY_SEPARATOR = ":"

NAPOLI_PATH = path.abspath(path.join(path.realpath(__file__), '../../../', 'tmp/nAPOLI'))
PDB_PATH = "%s/public/pdb" % NAPOLI_PATH
TMP_FILES = "%s/public/tmp" % NAPOLI_PATH

CONF_PATH = path.abspath(path.join(path.realpath(__file__), '../../../', 'data'))
DB_CONF_FILE = "%s/.mysql.ini" % CONF_PATH
ATOM_PROP_FILE = "%s/Napoli.fdef" % CONF_PATH

INTERACTION_CONF = DefaultInteractionConf()
BOUNDARY_CONF = InteractionConf({"boundary_cutoff": 7})

NMR_METHODS = ["SOLID-STATE NMR", "SOLUTION NMR"]

COV_SEARCH_RADIUS = 2.2

IFP_LENGTH = 1024

ACCEPTED_MOL_OBJ_TYPES = ("rdkit", "openbabel")

OPENBABEL = "/usr/bin/obabel"

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
}, (255, 255, 255))

CHEMICAL_FEATURE_IDS = {
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

INTERACTION_IDS = {
    "Proximal": 0,
    "Hydrogen bond": 1,
    "Ionic": 2,
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
    "Van der Waals": 18,
    "Chalcogen bond": 19,
    "Chalcogen-pi": 20,
    "Halogen-pi": 21,
    "Orthogonal multipolar": 22,
    "Parallel multipolar": 23,
    "Antiparallel multipolar": 24,
    "Tilted multipolar": 25,
    "Multipolar": 26,
    "Cation-nucleophile": 27,
    "Anion-electrophile": 28,
    "Unfavorable anion-nucleophile": 29,
    "Unfavorable cation-electrophile": 30,
    "Unfavorable nucleophile-nucleophile": 31,
    "Unfavorable electrophile-electrophile": 32,
}

INTERACTION_SHORT_NAMES = {
    "Proximal": "prox",
    "Hydrogen bond": "hbond",
    "Ionic": "ionic",
    "Salt bridge": "salt_bridge",
    "Cation-pi": "cation-pi",
    "Edge-to-face pi-stacking": "t-shaped_pi-stack",
    "Face-to-face pi-stacking": "sand_pi-stack",
    "Parallel-displaced pi-stacking": "disp_pi-stack",
    "Hydrophobic": "hphobe",
    "Halogen bond": "xbond",
    "Halogen bond": "x-pi",
    "Chalcogen bond": "ybond",
    "Chalcogen-pi": "y-pi",
    "Repulsive": "repuls",
    "Water-bridged hydrogen bond": "water_hbond",
    "Pi-stacking": "pi-stack",
    "Amide-aromatic stacking": "amide_pi-stack",
    "Weak hydrogen bond": "weak_hbond",
    "Orthogonal multipolar": "ort_multipol",
    "Parallel multipolar": "par_multipol",
    "Antiparallel multipolar": "antipar_multipol",
    "Tilted multipolar": "tilted_multipol",
    "Multipolar": "multipol",
    "Unfavorable nucleophile-nucleophile": "unf_nucleop_nucleop",
    "Unfavorable electrophile-electrophile": "unf_electrop_electrop",
    "Cation-nucleophile": "cat_nucleop",
    "Anion-electrophile": "ani_electrop",
    "Unfavorable anion-nucleophile": "unf_ani_nucleop",
    "Unfavorable cation-electrophile": "unf_cat_electrop",
    "Covalent bond": "cov",
    "Atom overlap": "atm_overlap",
    "Van der Waals clash": "vdw_clash",
    "Van der Waals": "vdw",
}

PYMOL_INTERACTION_COLOR = ColorPallete({
    "Proximal": "gray60",
    "Hydrogen bond": "tv_blue",
    "Water-bridged hydrogen bond": "lightblue",
    "Weak hydrogen bond": "lightteal",
    "Ionic": "green",
    "Salt bridge": "forest",
    "Cation-pi": "salmon",
    "Edge-to-face pi-stacking": "tv_red",
    "Face-to-face pi-stacking": "tv_red",
    "Parallel-displaced pi-stacking": "tv_red",
    "Pi-stacking": "tv_red",
    "Amide-aromatic stacking": "raspberry",
    "Hydrophobic": "orange",
    "Halogen bond": "aquamarine",
    "Halogen-pi": "aquamarine",
    "Chalcogen bond": "lightorange",
    "Chalcogen-pi": "lightorange",
    "Repulsive": "violetpurple",
    "Covalent bond": "black",
    "Van der Waals clash": "gray90",
    "Van der Waals": "gray50",
    "Orthogonal multipolar": "paleyellow",
    "Parallel multipolar": "paleyellow",
    "Antiparallel multipolar": "paleyellow",
    "Tilted multipolar": "paleyellow",
    "Multipolar": "paleyellow",
    "Unfavorable nucleophile-nucleophile": "dirtyviolet",
    "Unfavorable electrophile-electrophile": "dirtyviolet",
    "Cation-nucleophile": "palegreen",
    "Anion-electrophile": "palegreen",
    "Unfavorable anion-nucleophile": "dirtyviolet",
    "Unfavorable cation-electrophile": "dirtyviolet"
}, "white")
