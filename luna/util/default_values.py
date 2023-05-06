from os import path

from luna.interaction.config import DefaultInteractionConfig, InteractionConfig
from luna.util import ColorPallete


ENTRY_SEPARATOR = ":"

LUNA_PATH = path.abspath(path.join(path.realpath(__file__), '../../../'))
OUTPUT_PATH = "%s/output/public" % LUNA_PATH
PDB_PATH = "%s/pdb" % OUTPUT_PATH
TMP_FILES = "%s/tmp" % OUTPUT_PATH

CONFIG_PATH = path.abspath(path.join(path.realpath(__file__), '../../', 'data'))
DB_CONFIG_FILE = "%s/mysql.ini" % CONFIG_PATH
ATOM_PROP_FILE = "%s/LUNA.fdef" % CONFIG_PATH
LIGAND_EXPO_FILE = "%s/ligand_expo.tsv" % CONFIG_PATH
MIN_FDEF_FILE = '%s/MinimalFeatures.fdef' % CONFIG_PATH

INTERACTION_CONFIG = DefaultInteractionConfig()
BOUNDARY_CONFIG = InteractionConfig({"bsite_cutoff": 6.2,
                                     "cache_cutoff": 10})

NMR_METHODS = ["SOLID-STATE NMR", "SOLUTION NMR"]

COV_SEARCH_RADIUS = 2.2
METAL_COMPLEX_DIST = 2.8

IFP_LENGTH = 4096

ACCEPTED_MOL_OBJ_TYPES = ("rdkit", "openbabel")

OPENBABEL = "obabel"
ANTECHAMBER = "antechamber"
DOCK6_PATH = "/opt/dock6"
DOCK6 = f"{DOCK6_PATH}/bin/dock6"

RECURSION_LIMIT = 600000

ATOM_TYPES_COLOR = ColorPallete({
    "Aromatic": (255, 187, 204),
    "Acceptor": (158, 154, 200),
    "Donor": (17, 170, 119),
    "Hydrophobe": (254, 224, 144),
    "Hydrophobic": (254, 224, 144),
    "Negative": (251, 128, 114),
    "Positive": (31, 120, 180),
    "NegIonizable": (251, 128, 114),
    "PosIonizable": (31, 120, 180),
    "NegativelyIonizable": (251, 128, 114),
    "PositivelyIonizable": (31, 120, 180),
    "HalogenDonor": (196, 156, 148),
    "Metal": (96, 125, 139),
    "WeakDonor": (161, 217, 155),
    "WeakAcceptor": (188, 194, 220),
    "Electrophile": (251, 138, 124),
    "Nucleophile": (158, 218, 229),
    "ChalcogenDonor": (255, 193, 7),
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
    "Displaced face-to-slope pi-stacking": 38,
    "Single bond": 39,
    "Double bond": 40,
    "Triple bond": 41,
    "Aromatic bond": 42,
    "Other bond": 43,
    "Metal coordination": 44,
}

INTERACTION_SHORT_NAMES = {
    "Proximal": "prox",
    "Hydrogen bond": "hbond",
    "Ionic": "ionic",
    "Salt bridge": "salt_bridge",
    "Cation-pi": "cation-pi",
    "Hydrophobic": "hphobe",
    "Halogen bond": "xbond",
    "Halogen-pi": "x-pi",
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
    "Displaced face-to-slope pi-stacking": "disp_face-to-slope_stack",
    "Single bond": "single-bond",
    "Double bond": "double-bond",
    "Triple bond": "triple-bond",
    "Aromatic bond": "aromatic-bond",
    "Other bond": "other-bond",
    "Metal coordination": "metal-coord",
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
    "Atom overlap": "gray40",
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
    "Displaced face-to-slope pi-stacking": "tv_red",
    "Single bond": "black",
    "Double bond": "black",
    "Triple bond": "black",
    "Aromatic bond": "black",
    "Other bond": "black",
    "Metal coordination": "olive"
}, "white")

PYMOL_INTERACTION_COLOR_AS_RGB = ColorPallete({
    "Proximal": (0.6, 0.6, 0.6),
    "Hydrogen bond": (0.3, 0.3, 1.0),
    "Water-bridged hydrogen bond": (0.75, 0.75, 1.0),
    "Weak hydrogen bond": (0.4, 0.7, 0.7),
    "Ionic": (0.0, 1.0, 0.0),
    "Salt bridge": (0.2, 0.6, 0.2),
    "Cation-pi": (1.0, 0.6, 0.6),
    "Amide-aromatic stacking": (0.70, 0.30, 0.40),
    "Hydrophobic": (1.0, 0.5, 0.0),
    "Halogen bond": (0.5, 1, 1),
    "Halogen-pi": (0.5, 1, 1),
    "Chalcogen bond": (1.0, 0.8, 0.5),
    "Chalcogen-pi": (1.0, 0.8, 0.5),
    "Repulsive": (0.55, 0.25, 0.60),
    "Covalent bond": (0.0, 0.0, 0.0),
    "Atom overlap": (0.4, 0.4, 0.4),
    "Van der Waals clash": (0.9, 0.9, 0.9),
    "Van der Waals": (0.5, 0.5, 0.5),
    "Orthogonal multipolar": (1.0, 1.0, 0.5),
    "Parallel multipolar": (1.0, 1.0, 0.5),
    "Antiparallel multipolar": (1.0, 1.0, 0.5),
    "Tilted multipolar": (1.0, 1.0, 0.5),
    "Multipolar": (1.0, 1.0, 0.5),
    "Unfavorable nucleophile-nucleophile": (0.70, 0.50, 0.50),
    "Unfavorable electrophile-electrophile": (0.70, 0.50, 0.50),
    "Cation-nucleophile": (0.65, 0.9, 0.65),
    "Anion-electrophile": (0.65, 0.9, 0.65),
    "Unfavorable anion-nucleophile": (0.70, 0.50, 0.50),
    "Unfavorable cation-electrophile": (0.70, 0.50, 0.50),
    "Pi-stacking": (1.0, 0.2, 0.2),
    "Face-to-face pi-stacking": (1.0, 0.2, 0.2),
    "Face-to-edge pi-stacking": (1.0, 0.2, 0.2),
    "Face-to-slope pi-stacking": (1.0, 0.2, 0.2),
    "Edge-to-edge pi-stacking": (1.0, 0.2, 0.2),
    "Edge-to-face pi-stacking": (1.0, 0.2, 0.2),
    "Edge-to-slope pi-stacking": (1.0, 0.2, 0.2),
    "Displaced face-to-face pi-stacking": (1.0, 0.2, 0.2),
    "Displaced face-to-edge pi-stacking": (1.0, 0.2, 0.2),
    "Displaced face-to-slope pi-stacking": (1.0, 0.2, 0.2),
    "Single bond": (0.0, 0.0, 0.0),
    "Double bond": (0.0, 0.0, 0.0),
    "Triple bond": (0.0, 0.0, 0.0),
    "Aromatic bond": (0.0, 0.0, 0.0),
    "Other bond": (0.0, 0.0, 0.0),
    "Metal coordination": (0.0, 0.0, 0.0)
}, (1.0, 1.0, 1.0))


ARTIFACTS_LIST = ["ACE", "ACT", "BME", "CSD", "CSW", "EDO", "FMT", "GOL", "MSE", "NAG", "NO3", "PO4", "SGM", "SO4", "TPO"]
