from collections import Counter, defaultdict
import pandas as pd
from luna.mol.entry import MolEntry


class InteractingResidues:

    def __init__(self, level1=None, level2=None, level3=None):
        level1 = level1 or []
        level2 = level2 or []
        level3 = level3 or []

        self.level1 = Counter(level1)
        self.level2 = Counter(level2)
        self.level3 = Counter(level3)

    def by_residue(self):
        return self.level

    def by_interaction(self):
        return self.level2

    def by_atom(self):
        return self.level3

    def update(self, interacting_residues):
        self.level1 += interacting_residues.level1
        self.level2 += interacting_residues.level2
        self.level3 += interacting_residues.level3


def get_interacting_residues(interactions, targets, key_map={}):
    res = set()
    res_inter = set()
    res_inter_atm = set()

    for i in interactions:
        if (i.src_grp.compound in targets and i.trgt_grp.compound.is_residue()):

            # TODO: fix this code to accept the new format of AtomGroups as the atoms can be from different compounds
            # if (any([c in targets for c in i.src_grp.compounds]) and i.trgt_grp.is_residue())

            res.add(i.trgt_grp.compound)
            res_inter.add((i.trgt_grp.compound, i.type))
            res_inter_atm.add((i.trgt_grp.compound, i.type, i.trgt_grp))
        elif (i.trgt_grp.compound in targets and i.src_grp.compound.is_residue()):
            res.add(i.src_grp.compound)
            res_inter.add((i.src_grp.compound, i.type))
            res_inter_atm.add((i.src_grp.compound, i.type, i.src_grp))

    return InteractingResidues(res, res_inter, res_inter_atm)


class ResidueFrequencies(object):

    def __init__(self):
        self.entries = []
        self.interactions = []
        self.freq = 0

    def add(self, entry, interactions):
        self.entries.append(entry)
        self.interactions = list(set(self.interactions + list(interactions)))
        self.freq += 1


def calculate_residues_frequency(interaction_tuples):

    global_freq = defaultdict(ResidueFrequencies)
    for entry, inter_mngr in interaction_tuples:
        interacting_res_mapping = defaultdict(list)
        for inter in inter_mngr:
            comps1 = tuple(sorted([(r.parent.parent.id, r.parent.id, r.resname, r.id[1:]) for r in inter.src_grp.compounds]))
            comps2 = tuple(sorted([(r.parent.parent.id, r.parent.id, r.resname, r.id[1:]) for r in inter.trgt_grp.compounds]))
            comps1, comps2 = tuple(sorted([comps1, comps2]))

            interacting_res_mapping[(comps1, comps2, inter.type)].append(inter)

        for k in interacting_res_mapping:
            global_freq[k].add(entry, interacting_res_mapping[k])

    return global_freq


def generate_residue_matrix(interactions_mngrs, by_interaction=True):

    data_by_entry = defaultdict(lambda: defaultdict(int))
    residues = set()

    for inter_mngr in interactions_mngrs:
        entry = inter_mngr.entry

        for inter in inter_mngr:
            # Continue if no target is in the interaction.
            if not inter.src_grp.has_target() and not inter.trgt_grp.has_target():
                continue
            # Ignore interactions involving the same compounds.
            if inter.src_grp.compounds == inter.trgt_grp.compounds:
                continue

            if inter.src_grp.has_hetatm():
                comp1 = sorted(inter.src_grp.compounds)
                comp2 = sorted(inter.trgt_grp.compounds)
            elif inter.trgt_grp.has_hetatm():
                comp1 = sorted(inter.trgt_grp.compounds)
                comp2 = sorted(inter.src_grp.compounds)
            else:
                comp1 = sorted(inter.src_grp.compounds)
                comp2 = sorted(inter.trgt_grp.compounds)
                comp1, comp2 = sorted([comp1, comp2])

            comp1 = ";".join(["%s/%s/%d%s" % (r.parent.id, r.resname, r.id[1], r.id[2].strip()) for r in comp1])
            comp2 = ";".join(["%s/%s/%d%s" % (r.parent.id, r.resname, r.id[1], r.id[2].strip()) for r in comp2])

            entry_id = entry.mol_id if isinstance(entry, MolEntry) else entry.to_string()
            if by_interaction:
                key = (entry_id, inter.type)
            else:
                key = entry_id

            data_by_entry[key][comp2] += 1

            residues.add(comp2)

    heatmap_data = defaultdict(list)

    if by_interaction:
        entries = set([k[0] for k in data_by_entry.keys()])
        interactions = set([k[1] for k in data_by_entry.keys()])

        for e in entries:
            for i in interactions:
                for res in residues:
                    heatmap_data["entry"].append(e)
                    heatmap_data["interaction"].append(i)
                    heatmap_data["residues"].append(res)
                    heatmap_data["frequency"].append(data_by_entry[(e, i)][res])

    else:
        for key in data_by_entry:
            for res in data_by_entry[key]:
                heatmap_data["entry"].append(key)

                heatmap_data["residues"].append(res)
                heatmap_data["frequency"].append(data_by_entry[key][res])

    df = pd.DataFrame.from_dict(heatmap_data)

    if by_interaction:
        return pd.pivot_table(df, index=['entry', 'interaction'], columns='residues', values='frequency', fill_value=0)
    else:
        return pd.pivot_table(df, index='entry', columns='residues', values='frequency', fill_value=0)
