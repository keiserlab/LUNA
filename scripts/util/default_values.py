from interaction.calc_interactions import (DefaultInteractionConf,
                                           InteractionConf)
from mol.depiction import ColorPallete

from MyBio.PDB.PDBParser import PDBParser

ENTRIES_SEPARATOR = ":"

DEFAULT_NAPOLI_PATH = "/media/shared/UFMG/Workspace/nAPOLI_v2/tmp/nAPOLI"
DEFAULT_PDB_PATH = "%s/public/pdb" % DEFAULT_NAPOLI_PATH
DEFAULT_TMP_FILES = "%s/public/tmp" % DEFAULT_NAPOLI_PATH

DEFAULT_DB_CONF_FILE = "../data/.mysql.ini"

DEFAULT_ATOM_PROP_FILE = "../data/BaseFeatures.fdef"

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
