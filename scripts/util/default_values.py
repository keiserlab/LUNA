from interaction.calc_interactions import DefaultInteractionConf

ENTRIES_SEPARATOR = ":"

NAPOLI_PATH = "/media/shared/UFMG/Workspace/nAPOLI_v2/tmp/nAPOLI"
DEFAULT_PDB_PATH = "%s/public/pdb" % NAPOLI_PATH

DEFAULT_DB_CONF_FILE = "../data/.mysql.ini"

DEFAULT_ATOM_PROP_FILE = "../data/BaseFeatures.fdef"

DEFAULT_INTERACTION_CONF = DefaultInteractionConf()

# PLI_DEFAULT_VALUES = {
#     "pdbPath": None,
#     "dbConf": "../data/.mysql.ini",
#     "overwritePath": False,

#     "pdbTemplate": None,
#     "chainTemplate": None,

#     "atomPropFile": "../data/BaseFeatures.fdef",
#     "interCriteria": DefaultInteractionConf(),
#     "overviewFirstThenFilter": True,
#     "ph": 7,

#     "fpFunction": "pharm2d_fp",

#     "runFromStep": 1,
#     "runUntilStep": -1
# }
