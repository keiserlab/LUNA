from collections import defaultdict


def count_group_types(compound, keyMap={}):
    groupTypes = defaultdict(int)

    for group in compound.atomGroups:
        for feature in group.chemicalFeatures:

            key = feature.name
            if keyMap:
                if feature.name in keyMap:
                    key = keyMap[feature.name]
                else:
                    continue

            groupTypes[key] += 1

    return groupTypes


def count_interaction_types(interactions, keyMap={}):

    interactionTypes = defaultdict(int)

    for interaction in interactions:

        key = interaction.type
        if keyMap:
            if interaction.type in keyMap:
                key = keyMap[interaction.type]
            else:
                continue

        interactionTypes[key] += 1

    return interactionTypes
