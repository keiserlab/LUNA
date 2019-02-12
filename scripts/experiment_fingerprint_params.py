from mol.entry import *

from napoli import *
from util.default_values import *

import numpy as np
import util.stringcase as case

from util.file import create_directory


output_path = "%s/tmp/fp_experiments" % NAPOLI_PATH

mol_file = "%s/data/small_diverse_stratified.mol2" % output_path

entries = []
# input_file = "%s/data/zinc_test_list" % output_path
input_file = "%s/tmp/local_proj/mols/small_diverse_stratified_INPUT" % NAPOLI_PATH
with open(input_file, "r") as IN:
    for l in IN.readlines():
        entries.append(MolEntry("d4_receptor", l.strip(), mol_file))

print("Number of entries to be processed: %d." % len(entries))


experiment_settings = [
        [("experiment", "Interactions importance - no interactions"), ("set", [("add_non_cov", False), ("add_atom_atom", False), ("add_proximal", True)])],
        [("experiment", "Interactions importance - with interactions"), ("set", [("add_non_cov", True), ("add_atom_atom", True), ("add_proximal", False)])],
        [("experiment", "Interactions importance - all"), ("set", [("add_non_cov", True), ("add_atom_atom", True), ("add_proximal", True)])],

        [("experiment", "Hydrogen importance - no hydrogen"), ("set", [("add_hydrogen", False)])],
        [("experiment", "Hydrogen importance - with hydrogen"), ("set", [("add_hydrogen", True)])],

        [("experiment", "Shell parameters - fingerprint size"), ("list", [32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384]), ("set_values_to", "ifp_length")],
        [("experiment", "Shell parameters - number of levels"), ("range", (1, 16)), ("step", 1), ("set_values_to", "ifp_num_levels")],
        [("experiment", "Shell parameters - radius steps"), ("range", (1, 8)), ("step", 0.5), ("set_values_to", "ifp_radius_step")],
]

default_opts = {}
default_opts["overwrite_path"] = False
default_opts["has_local_files"] = True
default_opts["add_hydrogen"] = False
default_opts["preload_mol_files"] = False
default_opts["calc_mfp"] = False
default_opts["pdb_path"] = "%s/data" % output_path

processes = []
for s in experiment_settings:
    opt = dict(default_opts)

    settings = dict(s)

    project_name = None
    values = []
    for (key, val) in settings.items():
        if key == "experiment":
            project_name = case.snakecase(val.lower())
        elif key == "range":
            step = settings.get("step", 1)
            values = list(np.arange(val[0], val[1], step))
        elif key == "list":
            values = list(val)
        elif key == "set":
            for (param_name, param_val) in val:
                opt[param_name] = param_val

    project_path = "%s/%s" % (output_path, project_name)
    create_directory(project_path)

    if values:
        for val in values:
            opt = dict(opt)
            opt["working_path"] = "%s/%s" % (project_path, "value_%s" % str(val))
            opt["ifp_output"] = "%s/results/IFP__%s__val--%s.csv" % (opt["working_path"], project_name, str(val))

            if "set_values_to" in settings:
                opt[settings["set_values_to"]] = val

            processes.append(opt)
    else:
        opt["working_path"] = project_path
        opt["ifp_output"] = "%s/results/IFP__%s.csv" % (project_path, project_name)
        processes.append(opt)

print()
for opt in processes:
    from util import logging_ini
    print("-------------------------------------------------")
    print()
    print("Starting new process...")
    print("Project: %s" % opt["working_path"])
    print("Parameters:")
    for (key, val) in opt.items():
        print("\t-- %s:\t%s" % (key, str(val)))

    opt["entries"] = entries

    try:
        pli_obj = Fingerprint_PLI_Project(**opt)
        pli_obj()
        print(">>>> DONE!!!!")
    except Exception as e:
        print(e)
        print(">>>> ERROR FOUND")
    print("#################################################")
    print()
    print()
