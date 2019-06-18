from collections import defaultdict
import logging


logger = logging.getLogger()


# TODO: fix this code to accept the new format of AtomGroups as the atoms can be from different compounds

def count_group_types(compound, key_map={}):
    grp_types = defaultdict(int)

    for group in compound.atomGroups:
        for feature in group.chemicalFeatures:
            key = feature.format_name()
            if key_map:
                if key in key_map:
                    key = key_map[key]
                else:
                    logger.warning("Does not exist a corresponding mapping to the key '%s'. It will be ignored." % key)
                    continue

            grp_types[key] += 1

    return grp_types


def count_interaction_types(interactions, targets=None, key_map={}):

    interaction_types = defaultdict(int)
    seen_pairs = set()
    for i in interactions:
        contain_trgts = True
        if targets is not None:
            if (i.src_grp.compound not in targets and i.trgt_grp.compound not in targets):
                contain_trgts = False

        if contain_trgts:
            pair_key1 = (i.type, i.src_grp, i.trgt_grp)
            pair_key2 = (i.type, i.trgt_grp, i.src_grp)

            if pair_key1 in seen_pairs or pair_key2 in seen_pairs:
                continue

            seen_pairs.add(pair_key1)
            key = i.type
            if key_map:
                if i.type in key_map:
                    key = key_map[i.type]
                else:
                    continue

            interaction_types[key] += 1

    return interaction_types
