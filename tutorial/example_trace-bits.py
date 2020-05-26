from mol.entry import *

from luna import *
from util.default_values import *

from mol.interaction.filter import InteractionFilter
from mol.interaction.calc import InteractionCalculator
from mol.interaction.view import InteractionViewer


working_path = "./tmp/Jessica_data_evaluation"


# Create the entries
entry_str = "trajectory_ligand_14:A:UNK:522"


entry_str = "1KNE:P:M3L:9"

c1 = CompoundEntry.from_string(entry_str)

entries = [c1]


ic = InteractionCalculator(inter_filter=InteractionFilter.new_pli_filter())


#
#
# Set running options
#
#

print("Number of entries to be processed: %d." % len(entries))

opt = {}
opt["overwrite_path"] = True
opt["has_local_files"] = False
opt["preload_mol_files"] = False
opt["try_h_addition"] = False
opt["calc_mfp"] = False
opt["mol_obj_type"] = 'rdkit'
opt["amend_mol"] = True

opt["inter_calc"] = ic

opt["pdb_path"] = working_path
opt["working_path"] = "%s/luna_results" % working_path

opt["entries"] = entries

pli_obj = LocalProject(**opt)
pli_obj.run()


for inter_tuple in pli_obj.interactions:
    output_file = "%s/%s.pse" % (opt["pdb_path"], inter_tuple[0].to_string())
    piv = InteractionViewer(add_directional_arrows=True)
    piv.new_session([inter_tuple], output_file)
    print()
    print()

    # for i in inter_tuple[1].interactions:
    #     if i.type == "Hydrogen bond":
    #         print("-------------------")
    #         print(i)
    #         print()
    #         print(i.src_grp, "\t", i.src_grp.features)
    #         print(i.trgt_grp, "\t", i.trgt_grp.features)
    #         print()
    #         print(i.params)
    #         print()
    # print()
    # print()

print("DONE!!!")
print()
