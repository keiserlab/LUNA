from mol.entry import *

from napoli import *
from util.default_values import *

from mol.interaction.view import InteractionViewer
from mol.interaction.fp.view import ShellViewer

from mol.interaction.filter import InteractionFilter
from mol.interaction.calc import InteractionCalculator


entry_str = "3QQK:A:X02:497"

print(entry_str)

c1 = CompoundEntry.from_string(entry_str)
entries = [c1]

ic = InteractionCalculator(inter_filter=InteractionFilter.new_pli_filter())

working_path = "%s/tmp/local_proj" % NAPOLI_PATH

print("Number of entries to be processed: %d." % len(entries))

opt = {}
opt["working_path"] = working_path
opt["overwrite_path"] = False
opt["entries"] = entries
opt["inter_calc"] = ic
opt["mol_obj_type"] = 'rdkit'
opt["try_h_addition"] = True
opt["amend_mol"] = True

pli_obj = LocalProject(**opt)
pli_obj.run()

print("Creating a new Pymol session...")
output_file = "%s/%s.pse" % (working_path, c1.to_string())
piv = InteractionViewer(add_directional_arrows=True)
piv.new_session(pli_obj.interactions, output_file)

print("Creating a new CSV file with all interactions...")

interactions_set = set()
for pdb_file, inter_mngr in pli_obj.interactions:
    for i in inter_mngr.interactions:
        ag1 = ";".join(sorted(["/".join(a.full_atom_name.split("/")) for a in i.src_grp.atoms]))
        ag2 = ";".join(sorted(["/".join(a.full_atom_name.split("/")) for a in i.trgt_grp.atoms]))
        ag1, ag2 = sorted([ag1, ag2])
        interactions_set.add((ag1, ag2, i.type))

block = "\n".join([",".join(k) for k in sorted(interactions_set)])

inter_file = "%s/%s.csv" % (working_path, c1.to_string())
with open(inter_file, "w") as OUT:
    OUT.write(block)

print("DONE!!!")
