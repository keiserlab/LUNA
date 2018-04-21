from util.exceptions import EntityLevelError

import logging
logger = logging.getLogger(__name__)


def all_contacts_nh_search(model, radius=7, level='A'):

    try:
        from MyBio.entity_func import get_entity_level_name
        from MyBio.PDB.NeighborSearch import NeighborSearch

        if (level == 'S' or level == 'M' or level == 'C'):
            raise EntityLevelError("Minimum entity level to be chosen is: "
                                   "R (residues) or A (atoms)")

        entitiesName = get_entity_level_name()
        if (level not in entitiesName):
            raise EntityLevelError("The defined level %s does not exist"
                                   % (level))

        logger.info("Trying to select all contacts in the PDB file %s."
                    % (model.get_parent().id))

        allAtoms = list(model.get_atoms())
        ns = NeighborSearch(allAtoms)

        pairs = ns.search_all(radius, level)

        logger.info("Number of nearby %s(s) found: %d."
                    % (entitiesName[level], len(pairs)))

        return pairs
    except Exception as e:
        logger.exception(e)
        raise


def get_contacts_for_entity(model, source, target=None, radius=7, level='A'):

    try:
        from MyBio.entity_func import get_entity_level_name
        from MyBio.PDB import (Selection, NeighborSearch)
        from itertools import product

        if (level == 'S' or level == 'M' or level == 'C'):
            raise EntityLevelError("Minimum entity level to be chosen is: "
                                   "R (residues) or A (atoms)")

        entitiesName = get_entity_level_name()
        if (level not in entitiesName):
            raise EntityLevelError("The defined level %s does not exist"
                                   % (level))

        sourceAtoms = Selection.unfold_entities([source], 'A')

        targetAtoms = []
        if (target is None):
            targetAtoms = list(model.get_atoms())
        else:
            targetAtoms = Selection.unfold_entities(target, 'A')

        ns = NeighborSearch(targetAtoms)

        entities = set()
        for atom in sourceAtoms:
            entity = atom.get_parent_by_level(level)

            nearbyEntities = ns.search(atom.coord, radius, level)
            pairs = set(product([entity], nearbyEntities))

            entities.update(pairs)

        logger.info("Number of nearby %s(s) found: %d."
                    % (entitiesName[level], len(entities)))

        return entities
    except Exception as e:
        logger.exception(e)
        raise
