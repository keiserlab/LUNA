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
    "Hydrophobic": 5,
    "Halogen bond": 6,
    "Repulsive": 7,
    "Water-bridged hydrogen bond": 8,
    "Amide-aromatic stacking": 9,
    "Weak hydrogen bond": 10,
    "Covalent bond": 11,
    "Atom overlap": 12,
    "Van der Waals clash": 13,
    "Van der Waals": 14,
    "Chalcogen bond": 15,
    "Chalcogen-pi": 16,
    "Halogen-pi": 17,
    "Orthogonal multipolar": 18,
    "Parallel multipolar": 19,
    "Antiparallel multipolar": 20,
    "Tilted multipolar": 21,
    "Multipolar": 22,
    "Cation-nucleophile": 23,
    "Anion-electrophile": 24,
    "Unfavorable anion-nucleophile": 25,
    "Unfavorable cation-electrophile": 26,
    "Unfavorable nucleophile-nucleophile": 27,
    "Unfavorable electrophile-electrophile": 28,
    "Pi-stacking": 29,
    "Face-to-face pi-stacking": 30,
    "Face-to-edge pi-stacking": 31,
    "Face-to-slope pi-stacking": 32,
    "Edge-to-edge pi-stacking": 33,
    "Edge-to-face pi-stacking": 34,
    "Edge-to-slope pi-stacking": 35,
    "Displaced face-to-face pi-stacking": 36,
    "Displaced face-to-edge pi-stacking": 37,
    "Displaced face-to-slope pi-stacking": 38
}

INTERACTION_SHORT_NAMES = {
    "Proximal": "prox",
    "Hydrogen bond": "hbond",
    "Ionic": "ionic",
    "Salt bridge": "salt_bridge",
    "Cation-pi": "cation-pi",
    "Hydrophobic": "hphobe",
    "Halogen bond": "xbond",
    "Halogen bond": "x-pi",
    "Chalcogen bond": "ybond",
    "Chalcogen-pi": "y-pi",
    "Repulsive": "repuls",
    "Water-bridged hydrogen bond": "water_hbond",
    "Amide-aromatic stacking": "amide_stack",
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
    "Pi-stacking": "pi-stack",
    "Face-to-face pi-stacking": "face-to-face_stack",
    "Face-to-edge pi-stacking": "face-to-edge_stack",
    "Face-to-slope pi-stacking": "face-to-slope_stack",
    "Edge-to-edge pi-stacking": "edge-to-edge_stack",
    "Edge-to-face pi-stacking": "edge-to-face_stack",
    "Edge-to-slope pi-stacking": "edge-to-slope_stack",
    "Displaced face-to-face pi-stacking": "disp_face-to-face_stack",
    "Displaced face-to-edge pi-stacking": "disp_face-to-edge_stack",
    "Displaced face-to-slope pi-stacking": "disp_face-to-slope_stack"
}

PYMOL_INTERACTION_COLOR = ColorPallete({
    "Proximal": "gray60",
    "Hydrogen bond": "tv_blue",
    "Water-bridged hydrogen bond": "lightblue",
    "Weak hydrogen bond": "lightteal",
    "Ionic": "green",
    "Salt bridge": "forest",
    "Cation-pi": "salmon",
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
    "Unfavorable cation-electrophile": "dirtyviolet",
    "Pi-stacking": "tv_red",
    "Face-to-face pi-stacking": "tv_red",
    "Face-to-edge pi-stacking": "tv_red",
    "Face-to-slope pi-stacking": "tv_red",
    "Edge-to-edge pi-stacking": "tv_red",
    "Edge-to-face pi-stacking": "tv_red",
    "Edge-to-slope pi-stacking": "tv_red",
    "Displaced face-to-face pi-stacking": "tv_red",
    "Displaced face-to-edge pi-stacking": "tv_red",
    "Displaced face-to-slope pi-stacking": "tv_red"
}, "white")
