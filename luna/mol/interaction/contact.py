from itertools import product


from luna.util.default_values import COV_SEARCH_RADIUS, BOUNDARY_CONF
from luna.util.exceptions import EntityLevelError
from luna.MyBio.util import get_entity_level_name, is_covalently_bound
from luna.MyBio.PDB.NeighborSearch import NeighborSearch
from luna.MyBio.PDB import Selection


import logging

logger = logging.getLogger()


def all_contacts_search(entity, radius=BOUNDARY_CONF.boundary_cutoff, level='A'):
    try:
        if level == 'S' or level == 'M' or level == 'C':
            raise EntityLevelError("Maximum entity level to be chosen is: R (residues) or A (atoms)")

        entity_names = get_entity_level_name()
        if level not in entity_names:
            raise EntityLevelError("The defined level '%s' does not exist" % level)

        logger.debug("Trying to select all contacts in the PDB file %s." % entity.get_parent_by_level('S').id)

        all_atoms = list(entity.get_atoms())
        ns = NeighborSearch(all_atoms)
        pairs = ns.search_all(radius, level)

        logger.debug("Number of nearby %s(s) found: %d." % (entity_names[level].lower(), len(pairs)))

        return pairs
    except Exception as e:
        logger.exception(e)
        raise


def get_contacts_for_entity(entity, source, target=None, radius=BOUNDARY_CONF.boundary_cutoff, level='A'):
    try:
        if level == 'S' or level == 'M' or level == 'C':
            raise EntityLevelError("Maximum entity level to be chosen is: R (residues) or A (atoms)")

        entity_names = get_entity_level_name()
        if level not in entity_names:
            raise EntityLevelError("The defined level '%s' does not exist" % level)

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

        logger.debug("Number of nearby %s(s) found: %d." % (entity_names[level].lower(), len(entities)))
        return entities
    except Exception as e:
        logger.exception(e)
        raise


def get_cov_contacts_for_entity(entity, source, target=None):
    entities = get_contacts_for_entity(entity, source, target, radius=COV_SEARCH_RADIUS, level='A')

    cov_bonds = set()
    for atm1, atm2 in entities:
        if is_covalently_bound(atm1, atm2):
            cov_bonds.add(tuple(sorted([atm1, atm2], key=lambda x: x.serial_number)))

    return cov_bonds
