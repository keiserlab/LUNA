from collections import Counter, defaultdict


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
