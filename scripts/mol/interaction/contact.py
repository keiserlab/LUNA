from util.exceptions import EntityLevelError

from MyBio.util import get_entity_level_name
from MyBio.PDB.NeighborSearch import NeighborSearch
from MyBio.PDB import Selection

from itertools import product

import logging
logger = logging.getLogger(__name__)


def all_contacts_nh_search(entity, radius=7, level='A'):
    try:
        if level == 'S' or level == 'M' or level == 'C':
            raise EntityLevelError('Minimum entity level to be chosen is: R (residues) or A (atoms)')

        entity_names = get_entity_level_name()
        if level not in entity_names:
            raise EntityLevelError('The defined level %s does not exist' % level)

        logger.info('Trying to select all contacts in the PDB file %s.'
                    % entity.get_parent_by_level('S').id)

        all_atoms = list(entity.get_atoms())
        ns = NeighborSearch(all_atoms)
        pairs = ns.search_all(radius, level)

        logger.info('Number of nearby %s(s) found: %d.' % (entity_names[level], len(pairs)))

        return pairs
    except Exception as e:
        logger.exception(e)
        raise


def get_contacts_for_entity(entity, source, target=None, radius=7, level='A'):
    try:
        if level == 'S' or level == 'M' or level == 'C':
            raise EntityLevelError('Minimum entity level to be chosen is: R (residues) or A (atoms)')

        entity_names = get_entity_level_name()
        if level not in entity_names:
            raise EntityLevelError('The defined level %s does not exist' % level)

        source_atoms = Selection.unfold_entities([source], 'A')
        target_atoms = []
        if target is None:
            target_atoms = list(entity.get_atoms())
        else:
            target_atoms = Selection.unfold_entities(target, 'A')

        ns = NeighborSearch(target_atoms)
        entities = set()
        for atom in source_atoms:
            entity = atom.get_parent_by_level(level)
            nb_entities = ns.search(atom.coord, radius, level)
            pairs = set(product([entity], nb_entities))
            entities.update(pairs)

        logger.info('Number of nearby %s(s) found: %d.' % (entity_names[level], len(entities)))

        return entities
    except Exception as e:
        logger.exception(e)
        raise
