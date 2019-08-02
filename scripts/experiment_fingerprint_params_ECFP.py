from mol.entry import *

from napoli import *
from util.default_values import *

import numpy as np
import util.stringcase as case

from util.file import create_directory, get_filename

import copy


#
#
# Data load
#
#

# Default
output_path = "%s/tmp/fp_experiments" % NAPOLI_PATH
pdb_path = "%s/data" % output_path

#
# SMALL dataset
#
mol_file = "%s/data/small_diverse_stratified.mol2" % output_path
input_file = "%s/data/small_diverse_stratified_INPUT" % output_path

#
# MEDIUM dataset
#
mol_file = "%s/data/medium_diverse_stratified.mol2" % output_path
input_file = "%s/data/medium_diverse_stratified_INPUT" % output_path


# Create the entries
entries = []
with open(input_file, "r") as IN:
    for l in IN.readlines():
        entries.append(MolEntry("d4_receptor", l.strip(), mol_file))

print("Number of entries to be processed: %d." % len(entries))


#
#
# Project params and Experiments
#
#

mfp_length = 4096

experiment_settings = [
        # [("experiment", "MorganFP parameters - fingerprint size"), ("list", [1024, 2048, 4096, 8192, 16384]), ("set_values_to_dict", ("mfp_opts", "nBits"))],
        [("experiment", "%s__MorganFP parameters - radius size" % get_filename(input_file).upper()), ("range", (1, 9)), ("step", 1),
         ("set_values_to_dict", ("mfp_opts", "radius"))],
]


#
#
# Set default options
#
#

default_opts = {}
default_opts["overwrite_path"] = False
default_opts["has_local_files"] = True
default_opts["preload_mol_files"] = True
default_opts["try_h_addition"] = False
default_opts["mol_obj_type"] = 'rdkit'

default_opts["calc_ifp"] = False
default_opts["calc_mfp"] = True
default_opts["mfp_opts"] = {"fp_function": "morgan_fp", "nBits": mfp_length}

default_opts["pdb_path"] = pdb_path


#
#
# Loading the experiments to be processed
#
#

processes = []
for s in experiment_settings:
    opt = copy.deepcopy(default_opts)

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
            opt = copy.deepcopy(opt)
            opt["working_path"] = "%s/%s" % (project_path, "value_%s" % str(val))
            opt["mfp_output"] = "%s/results/MFP__%s__val--%s.csv" % (opt["working_path"], project_name, str(val))

            if "set_values_to_dict" in settings:
                dict_key, opt_key = settings["set_values_to_dict"]
                opt[dict_key][opt_key] = val
            processes.append(opt)
    else:
        opt["working_path"] = project_path
        opt["mfp_output"] = "%s/results/MFP__%s.csv" % (project_path, project_name)
        processes.append(opt)


#
#
# Running
#
#

print()
for opt in processes:
    print("-------------------------------------------------")
    print()
    print("Starting new process...")
    print("Project: %s" % opt["working_path"])
    print("Parameters:")
    for (key, val) in opt.items():
        if key == "entries":
            continue
        print("\t-- %s:\t%s" % (key, str(val)))

    opt["entries"] = entries

    try:
        pli_obj = FingerprintProject(**opt)
        pli_obj()
        print(">>>> DONE!!!!")
    except Exception as e:
        print(e)
        print(">>>> ERROR FOUND")
    print("#################################################")
    print()
    print()
