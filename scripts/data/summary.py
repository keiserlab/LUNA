from collections import defaultdict


def count_group_types(compound, key_map={}):
    grp_types = defaultdict(int)

    for group in compound.atomGroups:
        for feature in group.chemicalFeatures:
            key = feature.name
            if key_map:
                if feature.name in key_map:
                    key = key_map[feature.name]
                else:
                    continue

            grp_types[key] += 1

    return grp_types


def count_interaction_types(interactions, key_map={}):

    interaction_types = defaultdict(int)
    seen_pairs = set()
    for i in interactions:
        pair_key1 = (i.type, i.comp1, i.comp2)
        pair_key2 = (i.type, i.comp2, i.comp1)

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
