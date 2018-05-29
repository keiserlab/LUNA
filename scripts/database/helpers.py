from database.napoli_model import *

from sqlalchemy.orm import relationship

import logging
logger = logging.getLogger(__name__)


class MapperConfiguration():

    def __init__(self, table_class, table_name, properties=None):
        self.table_class = table_class
        self.table_name = table_name
        self.properties = properties


class Manager:

    def __init__(self, db, mapper_list):
        self._db = db

        for conf in mapper_list:
            self.db.new_mapper(conf.table_class,
                               conf.table_name,
                               conf.properties)

    @property
    def db(self):
        return self._db


class LigandEntryManager(Manager):

    def __init__(self, db):
        mapper_conf = MapperConfiguration(LigandEntry, "ligand_entry")
        super().__init__(db, [mapper_conf])

    def get_entries(self, projectId):
        ligand_entries = (self.db.session
                          .query(LigandEntry)
                          .filter(LigandEntry.inter_proj_id == projectId)
                          .all())

        return ligand_entries


class InteractionManager(Manager):

    def __init__(self, db, mapper_list):
        mapper_list.append(MapperConfiguration(InterType, "inter_type"))
        mapper_list.append(MapperConfiguration(InterDependInter, "inter_depend_inter"))
        mapper_list.append(MapperConfiguration(GroupCompoundType, "group_compound_type"))
        mapper_list.append(MapperConfiguration(CompoundType, "compound_type"))
        mapper_list.append(MapperConfiguration(GroupAtom, "group_atom"))
        mapper_list.append(MapperConfiguration(Atom, "atom"))

        grp_comp_type_tbl = db.get_table("group_compound_type")
        grp_atom_tbl = db.get_table("group_atom")
        grp_properties = {
            "compound_types": relationship(CompoundType,
                                           secondary=grp_comp_type_tbl,
                                           backref='group'),
            "atoms": relationship(Atom,
                                  secondary=grp_atom_tbl,
                                  backref='group')
        }
        mapper_list.append(MapperConfiguration(Group, "group",
                                               properties=grp_properties))

        super().__init__(db, mapper_list)

    def insert_interactions(self, interactions, ligand_entity):
        # Create a copy of the interactions to modify the list.
        auxList = list(interactions)

        interInDB = {}
        # Control if an interaction with dependents failed once.
        retries = set()
        while auxList:
            interaction = auxList.pop(0)

            # If an interaction has dependecies, i.e., it depends on
            # any other interaction to exist.
            # It occurs, for example, in Water-bridged hydrogen bonds
            # as this interaction depends on two or more hydrogen bonds
            # to exist.
            if "dependOf" in interaction.params:
                dependOfInter = interaction.params["dependOf"]
                existDependencies = True
                invalidDependencies = False
                for requiredInter in dependOfInter:
                    if requiredInter not in interInDB:
                        existDependencies = False
                    else:
                        if interInDB[requiredInter] is None:
                            invalidDependencies = True

                # If any required interaction is invalid,
                # it does not add the current interaction
                if invalidDependencies:
                    logger.info("The interaction '%s' depend on some invalid "
                                "interactions: %s." %
                                (interaction, str(dependOfInter)))

                # If all required interactions has already
                # been added to the DB.
                if existDependencies:
                    dbRequiredInteractions = ([interInDB[x]
                                               for x in dependOfInter])

                    dbInteraction = self.new_interaction(interaction,
                                                         ligand_entity)

                    if dbInteraction:
                        dbInteraction.dependencies = dbRequiredInteractions
                    else:
                        logger.info("The interaction '%s' does "
                                    "not have a valid type." % interaction)
                else:
                    # If the interaction is in the 'retries' list
                    # it means that it is the second try to the
                    # same interaction. So, ignore it.
                    if interaction in retries:
                        logger.info("The interaction '%s' depend on not "
                                    "existing interactions: %s." %
                                    (interaction, str(dependOfInter)))
                    else:
                        retries.add(interaction)
                        auxList.append(interaction)
            else:
                dbInteraction = self.new_interaction(interaction,
                                                     ligand_entity)
                if dbInteraction is None:
                    logger.info("The interaction '%s' does "
                                "not have a valid type." % interaction)

                interInDB[interaction] = dbInteraction

        # Approve all insertions.
        self.db.approve_session()

    def new_interaction(self, interaction, ligand_entity):
        g1 = self.new_atom_group(interaction.comp1)
        g2 = self.new_atom_group(interaction.comp2)

        interType = interaction.type
        dbInterType = (self.db.session.query(InterType).
                       filter(InterType.type == interType).
                       all())

        params = dict(interaction.params)
        if ("dependOf" in params):
            del params["dependOf"]

        dbInteraction = Interaction(group_id1=g1.id,
                                    group_id2=g2.id,
                                    inter_type_id=dbInterType[0].id,
                                    **params)

        if ligand_entity:
            dbInteraction.ligand_entities.append(ligand_entity)

        self.db.session.add(dbInteraction)

        return dbInteraction

    def new_atom_group(self, group):
        oriComp = group.atoms[0].get_parent()

        db_group = Group(comp_pdb_id=oriComp.get_parent_by_level('S').id,
                         comp_chain_id=oriComp.get_parent_by_level('C').id,
                         comp_name=oriComp.get_parent_by_level('R').resname,
                         comp_num=oriComp.get_parent_by_level('R').get_id()[1],
                         comp_icode=oriComp.get_parent_by_level('R').get_id()[2],
                         name="",
                         coord_x=group.centroid[0],
                         coord_y=group.centroid[1],
                         coord_z=group.centroid[2])
        self.db.session.add(db_group)

        # TODO: o nome dos grupos esta ficando como?
        featureNames = [x.name for x in group.chemicalFeatures]

        # If the 'featureNames' list is empty or has non-existing values,
        # the DB selection will return an empty list.
        # The chemical features can be empty if the group does not match
        # any chemical definition.
        # Despite it, two atoms can establish an interaction of some type,
        # e.g., a Covalent bond.
        dbChemicalFeatures = (self.db.session.query(CompoundType).
                              filter(CompoundType.type.in_(featureNames)).
                              all())

        for dbObj in dbChemicalFeatures:
            db_group.compound_types.append(dbObj)

        for atom in group.atoms:
            db_atom = self.new_atom(atom)
            db_group.atoms.append(db_atom)

        return db_group

    def new_atom(self, atom):
        db_atom = Atom(name=atom.name,
                       serial_number=atom.serial_number,
                       altloc=atom.altloc,
                       element=atom.element,
                       coord_x=atom.coord[0],
                       coord_y=atom.coord[1],
                       coord_z=atom.coord[2])

        self.db.session.add(db_atom)

        return db_atom


class RCSBInteractionManager(InteractionManager):

    def __init__(self, db):
        mapper_list = []
        mapper_list.append(MapperConfiguration(LigandEstablishInteraction,
                                               "ligand_establish_interaction"))

        ligand_inter_tbl = db.get_table("ligand_establish_interaction")
        inter_inter_tbl = db.get_table("inter_depend_inter")
        inter_tbl = db.get_table("interaction")
        group_tbl = db.get_table("group")

        # TODO: colocar join com tipo de interações
        # TODO: colocar join com tipo de átomos
        inter_properties = {
            "dependencies": relationship(Interaction,
                                         primaryjoin=(inter_tbl.c.id ==
                                                      inter_inter_tbl.c.interaction_id),
                                         secondaryjoin=(inter_tbl.c.id ==
                                                        inter_inter_tbl.c.interaction_id1),
                                         secondary=inter_inter_tbl,
                                         backref='interaction'),
            "ligand_entities": relationship(Ligand,
                                            secondary=ligand_inter_tbl,
                                            backref='interaction'),
            "group1": relationship(Group,
                                   backref='interactions1',
                                   primaryjoin=(inter_tbl.c.group_id1 ==
                                                group_tbl.c.id)),
            "group2": relationship(Group,
                                   backref='interactions2',
                                   primaryjoin=(inter_tbl.c.group_id2 ==
                                                group_tbl.c.id)),
            "inter_type": relationship(InterType)
        }
        mapper_list.append(MapperConfiguration(Interaction, inter_tbl,
                                               properties=inter_properties))

        ligand_tbl = db.get_table("ligand")
        ligand_properties = {
            "interactions": relationship(Interaction,
                                         secondary=ligand_inter_tbl,
                                         backref='ligand')
        }
        mapper_list.append(MapperConfiguration(Ligand, ligand_tbl,
                                               properties=ligand_properties))

        super().__init__(db, mapper_list)

    def select_interactions(self, join_filter, inter_filters):

        query = (self.db.session.query(Interaction)
                 .join(Interaction.ligand_entities)
                 .filter(*join_filter.rules))

        interactions = set()
        for inter_filter in inter_filters:
            inter = query.filter(*inter_filter.rules).all()
            interactions.update(inter)

        aux_set = set(interactions)
        keep_inter = set()
        remove_inter = set()
        for db_inter in interactions:
            if db_inter.dependencies:
                remove = False
                for db_req_inter in db_inter.dependencies:
                    if db_req_inter not in interactions:
                        remove = True
                        break

                if remove:
                    remove_inter.update(set(db_inter.dependencies))
                    aux_set.remove(db_inter)
                else:
                    keep_inter.update(set(db_inter.dependencies))

        for db_inter in remove_inter:
            if db_inter in aux_set and db_inter not in keep_inter:
                aux_set.remove(db_inter)

        return list(aux_set)


class ProjectInteractionManager(InteractionManager):

    def __init__(self, db):
        mapper_list = []
        mapper_list.append(MapperConfiguration(LigandEntry, "ligand_entry"))
        mapper_list.append(MapperConfiguration(LigandEntryHasInteractions,
                                               "ligand_entry_has_interactions"))

        ligand_inter_tbl = db.get_table("ligand_entry_has_interactions")
        inter_inter_tbl = db.get_table("inter_depend_inter")
        inter_tbl = db.get_table("interaction")
        inter_properties = {
            "dependencies": relationship(Interaction,
                                         primaryjoin=(inter_tbl.c.id ==
                                                      inter_inter_tbl.c.interaction_id),
                                         secondaryjoin=(inter_tbl.c.id ==
                                                        inter_inter_tbl.c.interaction_id1),
                                         secondary=inter_inter_tbl,
                                         backref='interaction'),
            "ligand_entities": relationship(LigandEntry,
                                            secondary=ligand_inter_tbl,
                                            backref='interaction')
        }
        mapper_list.append(MapperConfiguration(Interaction, inter_tbl,
                                               properties=inter_properties))

        super().__init__(db, mapper_list)


    # def select_interactions(self, filterRules):
    #     query = self.db.session.query(Interaction)

    #     interactions = set()

    #     for filterRule in filterRules:
    #         inter = query.filter(*filterRule.rules).all()
    #         interactions.update(inter)

    #     invalidInteractions = set()
    #     for dbInter in interactions:
    #         if dbInter.dependencies:
    #             hasAllDependencies = True
    #             for dbRequiredInter in dbInter.dependencies:
    #                 if dbRequiredInter not in interactions:
    #                     hasAllDependencies = False
    #                     break

    #             if not hasAllDependencies:
    #                 invalidInteractions.add(dbInter)

    #     interactions.difference_update(invalidInteractions)

    #     return interactions



