from util.exceptions import (EntityLevelError,
                             ChainNotFoundError,
                             ResidueNotFoundError)

import logging
logger = logging.getLogger(__name__)


def all_contacts_nh_search(model, radius=7, level='A'):

    try:
        from bio.entity_func import get_entity_level_name
        from Bio.PDB import NeighborSearch

        if (level == 'S' or level == 'M'):
            raise EntityLevelError("Minimum entity level to be chosen is: "
                                   "R (residues) or A (atoms)")

        entitiesName = get_entity_level_name()
        if (level not in entitiesName):
            raise EntityLevelError("The defined level %s does not exist"
                                   % (level))

        logger.info("Trying to select all contacts in the PDB file %s"
                    % (model.get_parent().id))

        allAtoms = list(model.get_atoms())
        ns = NeighborSearch(allAtoms)

        pairs = ns.search_all(radius, level)

        logger.info("Number of nearby %s(s) found: %d"
                    % (entitiesName[level], len(pairs)))

        return pairs
    except Exception as e:
        logger.exception(e)
        raise


def get_contacts_for_entity(model, source, target=None, radius=7, level='A'):

    try:
        from bio.entity_func import (get_entity_level_name,
                                     get_parent_by_level)
        from Bio.PDB import (Selection, NeighborSearch)
        from itertools import product

        if (level == 'S' or level == 'M'):
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
            entity = get_parent_by_level(atom, level)

            nearbyEntities = ns.search(atom.coord, radius, level)
            pairs = set(product([entity], nearbyEntities))

            entities.update(pairs)

        logger.info("Number of nearby %s(s) found: %d"
                    % (entitiesName[level], len(entities)))

        return entities
    except Exception as e:
        logger.exception(e)
        raise


def nearby_entities_nb_search(model, target,
                              radius=7, level='A'):

    try:
        from util.default_values import get_entity_level_name
        from Bio.PDB import (Selection, NeighborSearch)

        entitiesName = get_entity_level_name()
        if (level not in entitiesName):
            raise EntityLevelError("The defined level %s does not exist"
                                   % (level))

        if targetChain not in model.child_dict:
            raise ChainNotFoundError("The defined chain %s does not exist"
                                     "in the PDB file %s"
                                     % (targetChain, model.get_parent().id))

        chain = model[targetChain]
        targetResidues = None
        if (resTuple is not None):
            logger.info("Trying to select %s(s) in contact with "
                        "%s from chain %s..." % (entitiesName[level],
                                                 str(resTuple), targetChain))

            if chain.has_id(resTuple) is False:
                raise ResidueNotFoundError("The defined ligand %s "
                                           "does not exist in the chain "
                                           "%s.%s" % (str(resTuple),
                                                      model.get_parent().id,
                                                      targetChain))
            else:
                targetResidues = [chain[resTuple]]
        else:
            logger.info("Trying to select %s(s) in contact with the chain %s"
                        % (entitiesName[level], targetChain))

            targetResidues = list(chain.get_residues())

        targetAtoms = Selection.unfold_entities(targetResidues, 'A')

        allAtoms = list(model.get_atoms())
        ns = NeighborSearch(allAtoms)

        entities = set()
        for atom in targetAtoms:
            entities.update(ns.search(atom.coord, radius, level))

        logger.info("Number of nearby %s(s) found: %d"
                    % (entitiesName[level], len(entities)))

        return entities

    except Exception as e:
        logger.exception(e)
        raise


# def filter_interacting_entities(model, targetChain, resTuple, exclusion):

