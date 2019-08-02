from mol.entry import *

from napoli import *

from mol.interaction.filter import InteractionFilter
from mol.interaction.calc import InteractionCalculator

from util.default_values import *
from util.file import get_filename


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

#
#
# Project params
#
#

ifp_num_levels = 7
ifp_radius_step = 1
ifp_length = 4096
add_proximal = False

project_name = "%s___Run__levels-%d__radius-%.1f__length-%d__proximal-%s" % (get_filename(input_file).upper(), ifp_num_levels,
                                                                             ifp_radius_step, ifp_length, add_proximal)
working_path = "%s/%s" % (output_path, project_name)
ifp_file = "%s/results/IFP__%s.csv" % (working_path, project_name)

ic = InteractionCalculator(inter_filter=InteractionFilter.new_pli_filter(), add_proximal=add_proximal)

#
#
# Set running options
#
#

print("Number of entries to be processed: %d." % len(entries))

opt = {}
opt["overwrite_path"] = True
opt["has_local_files"] = True
opt["preload_mol_files"] = True
opt["try_h_addition"] = False
opt["calc_mfp"] = False
opt["mol_obj_type"] = 'rdkit'

opt["ifp_num_levels"] = ifp_num_levels
opt["ifp_radius_step"] = ifp_radius_step
opt["ifp_length"] = ifp_length

opt["inter_calc"] = ic

opt["pdb_path"] = pdb_path
opt["working_path"] = working_path
opt["ifp_output"] = ifp_file

opt["entries"] = entries


#
#
# Running
#
#

print("-------------------------------------------------")
print()
print("Starting new process...")
print("Project: %s" % opt["working_path"])
print("Parameters:")
for (key, val) in opt.items():
    print("\t-- %s:\t%s" % (key, str(val)))

try:
    pli_obj = FingerprintProject(**opt)
    pli_obj()
    print(">>>> DONE!!!!")
except Exception as e:
    import traceback
    print(traceback.format_exc())
    print(">>>>", e)
    print(">>>> ERROR FOUND")
print()
print("#################################################")
print()
print()
