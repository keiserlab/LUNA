from database.napoli_model import *

from sqlalchemy.orm import relationship

import logging
logger = logging.getLogger(__name__)


class Manager:

    def __init__(self, db, kwargs):
        self._db = db

        for targetClass, targetTable in kwargs.items():
            self._db.new_mapper(targetClass, targetTable)


class ComplexManager(Manager):

    def __init__(self, db):
        kwargs = {Complex: "complex"}

        super().__init__(db, kwargs)

    def get_complexes(self, projectId):

        dbComplexes = (self._db.session
                           .query(Complex)
                           .filter(Complex.inter_proj_id == projectId)
                           .all())

        return dbComplexes


class InteractionManager():

    def __init__(self, db):
        self._db = db
        self._map_tables()

    def _map_tables(self):
        self._db.new_mapper(InterType, "inter_type")

        self._db.new_mapper(Complex, "complex")

        self._db.new_mapper(ComplexInter, "complex_inter")
        relComplexInter = self._db._metadata.tables['complex_inter']

        self._db.new_mapper(InterDependInter, "inter_depend_inter")
        relInterInter = self._db._metadata.tables['inter_depend_inter']

        interTable = self._db.get_table("interaction")
        self._db.new_mapper(Interaction, interTable, properties={
                            "dependencies": relationship(Interaction,
                                                         primaryjoin=interTable.c.id==relInterInter.c.interaction_id,
                                                         secondaryjoin=interTable.c.id==relInterInter.c.interaction_id1,
                                                         secondary=relInterInter,
                                                         backref='interaction'),
                            "complexes": relationship(Complex,
                                                      secondary=relComplexInter,
                                                      backref='interaction')
                            })

        self._db.new_mapper(GroupCompoundType, "group_compound_type")
        self._db.new_mapper(CompoundType, "compound_type")
        self._db.new_mapper(GroupAtom, "group_atom")
        self._db.new_mapper(Atom, "atom")

        relGroupComp = self._db._metadata.tables['group_compound_type']
        relGroupAtom = self._db._metadata.tables['group_atom']
        self._db.new_mapper(Group, "group", properties={
                            "compound_types": relationship(CompoundType,
                                                           secondary=relGroupComp,
                                                           backref='group'),
                            "atoms": relationship(Atom,
                                                  secondary=relGroupAtom,
                                                  backref='group')
                            })
        self._db.new_mapper(PDB, "pdb")

    def insert_interactions(self, interactions, dbComplex=None):
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

                    dbInteraction = self._new_interaction(interaction, dbComplex)

                    if dbInteraction is None:
                        logger.info("The interaction '%s' does "
                                    "not have a valid type." % interaction)
                    else:
                        dbInteraction.dependencies = dbRequiredInteractions
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
                dbInteraction = self._new_interaction(interaction, dbComplex)

                if dbInteraction is None:
                    logger.info("The interaction '%s' does "
                                "not have a valid type." % interaction)

                interInDB[interaction] = dbInteraction

        # Approve all insertions.
        self._db.approve_session()

    def _new_interaction(self, interaction, dbComplex=None):
        if not interaction.type:
            return

        g1 = self._new_atom_group(interaction.comp1)
        g2 = self._new_atom_group(interaction.comp2)

        interType = interaction.type
        dbInterType = (self._db.session.query(InterType).
                       filter(InterType.type == interType).
                       all())

        pdbId = None
        params = dict(interaction.params)
        if ("dependOf" in params):
            del params["dependOf"]

        dbInteraction = Interaction(group_id1=g1.id,
                                    group_id2=g2.id,
                                    inter_type_id=dbInterType[0].id,
                                    pdb_id=pdbId,
                                    **params
                                    )

        if dbComplex:
            dbInteraction.complexes.append(dbComplex)

        self._db.session.add(dbInteraction)

        return dbInteraction

    def _new_atom_group(self, group):
        oriComp = group.atoms[0].get_parent()

        dbGroup = Group(orig_comp_pdb=oriComp.get_parent_by_level('S').id,
                        orig_comp_chain=oriComp.get_parent_by_level('C').id,
                        orig_comp_name=oriComp.get_parent_by_level('R').resname,
                        orig_comp_num=oriComp.get_parent_by_level('R').get_id()[1],
                        orig_comp_icode=oriComp.get_parent_by_level('R').get_id()[2],
                        name="",
                        coord_x=group.centroid[0],
                        coord_y=group.centroid[1],
                        coord_z=group.centroid[2])

        self._db.session.add(dbGroup)

        # TODO: o nome dos grupos esta ficando como?
        featureNames = [x.name for x in group.chemicalFeatures]

        # If the 'featureNames' list is empty or has non-existing values,
        # the DB selection will return an empty list.
        # The chemical features can be empty if the group does not match
        # any chemical definition.
        # Despite it, two atoms can establish an interaction of some type, e.g., a Covalent bond.
        dbChemicalFeatures = (self._db.session.query(CompoundType).
                              filter(CompoundType.type.in_(featureNames)).
                              all())

        for dbObj in dbChemicalFeatures:
            dbGroup.compound_types.append(dbObj)

        for atom in group.atoms:
            dbAtom = self._new_atom(atom)
            dbGroup.atoms.append(dbAtom)

        return dbGroup

    def _new_atom(self, atom):
        dbAtom = Atom(name=atom.name,
                      serial_number=atom.serial_number,
                      altloc=atom.altloc,
                      element=atom.element,
                      coord_x=atom.coord[0],
                      coord_y=atom.coord[1],
                      coord_z=atom.coord[2]
                      )
        self._db.session.add(dbAtom)
        return dbAtom

    # def select_interactions(self, interConf):
    #     self._db.new_session()

    #     query = self._db.session.query(Interaction)

    #     a = Interaction.ha_dist_hb_inter <= interConf.max_ha_dist_hb_inter
    #     b = Interaction.da_dist_hb_inter <= interConf.max_da_dist_hb_inter

    #     hbond = query.filter(a, b)

    #     print(len(query.all()))
    #     print(len(hbond.all()))

    def select_interactions(self, filterRules):
        query = self._db.session.query(Interaction)

        interactions = set()

        for filterRule in filterRules:
            inter = query.filter(*filterRule.rules).all()
            interactions.update(inter)

        invalidInteractions = set()
        for dbInter in interactions:
            if dbInter.dependencies:
                hasAllDependencies = True
                for dbRequiredInter in dbInter.dependencies:
                    if dbRequiredInter not in interactions:
                        hasAllDependencies = False
                        break

                if not hasAllDependencies:
                    invalidInteractions.add(dbInter)

        interactions.difference_update(invalidInteractions)

        return interactions
