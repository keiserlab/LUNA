from Bio.PDB.Polypeptide import is_aa

from sqlalchemy import inspect
from sqlalchemy.orm import relationship

from database.luna_model import *
from database.filters import FilterRules
from database.helpers import MapperConfiguration

from mol.entry import DBEntry
from mol.interaction.type import InteractionType
from mol.groups import AtomGroup
from mol.features import ChemicalFeature


def object_as_dict(obj):
    return {c.key: getattr(obj, c.key)
            for c in inspect(obj).mapper.column_attrs}


def rows_2dict(rows):
    return [object_as_dict(r) for r in rows]


def get_ligand_tbl_join_filter(ligand_entry, target_class):
    rules = []
    rules.append(target_class.pdb_id == ligand_entry.pdb_id)
    rules.append(target_class.chain_id == ligand_entry.chain_id)
    rules.append(target_class.lig_name == ligand_entry.lig_name)
    rules.append(target_class.lig_num == ligand_entry.lig_num)
    rules.append(target_class.lig_icode == ligand_entry.lig_icode)
    join_filter = FilterRules("Join to %s" % ligand_entry, rules)

    return join_filter


def default_interaction_filters(interIdByType, interConf):
    filters = []

    ##################################################
    interType = "Hydrogen bond"
    if interType in interIdByType:
        rules = []
        rules.append(Interaction.inter_type_id == interIdByType[interType])
        rules.append(Interaction.da_dist_hb_inter <= interConf.max_da_dist_hb_inter)
        rules.append(Interaction.ha_dist_hb_inter <= interConf.max_ha_dist_hb_inter)
        rules.append(Interaction.dha_ang_hb_inter >= interConf.min_dha_ang_hb_inter)
        filters.append(FilterRules(interType, rules))

        # Hydrogen bonds in which the Hydrogen atom of the Donor atom could not be defined.
        # Commonly, it would just occur for Water molecules.
        rules = []
        rules.append(Interaction.inter_type_id == interIdByType[interType])
        rules.append(Interaction.da_dist_hb_inter <= interConf.max_da_dist_hb_inter)
        rules.append(Interaction.da_dist_hb_inter - 1 <= interConf.max_ha_dist_hb_inter)
        rules.append(Interaction.dha_ang_hb_inter == -1)
        filters.append(FilterRules(interType, rules))

    ##################################################
    interType = "Repulsive"
    if interType in interIdByType:
        rules = []
        rules.append(Interaction.inter_type_id == interIdByType[interType])
        rules.append(Interaction.dist_repuls_inter <= interConf.max_dist_repuls_inter)
        filters.append(FilterRules(interType, rules))

    ##################################################
    interType = "Attractive"
    if interType in interIdByType:
        rules = []
        rules.append(Interaction.inter_type_id == interIdByType[interType])
        rules.append(Interaction.dist_attract_inter <= interConf.max_dist_attract_inter)
        filters.append(FilterRules(interType, rules))

    ##################################################
    interType = "Edge-to-face pi-stacking"
    if interType in interIdByType:
        rules = []
        rules.append(Interaction.inter_type_id == interIdByType[interType])
        rules.append(Interaction.cc_dist_pi_pi_inter <= interConf.max_cc_dist_pi_pi_inter)
        rules.append(Interaction.dihed_ang_pi_pi_inter >= interConf.min_dihed_ang_pi_pi_inter)
        filters.append(FilterRules(interType, rules))

    ##################################################
    interType = "Face-to-face pi-stacking"
    if interType in interIdByType:
        rules = []
        rules.append(Interaction.inter_type_id == interIdByType[interType])
        rules.append(Interaction.cc_dist_pi_pi_inter <= interConf.max_cc_dist_pi_pi_inter)
        rules.append(Interaction.dihed_ang_pi_pi_inter < interConf.min_dihed_ang_pi_pi_inter)
        rules.append(Interaction.disp_ang_pi_pi_inter <= interConf.max_disp_ang_pi_pi_inter)
        filters.append(FilterRules(interType, rules))

    ##################################################
    interType = "Parallel-displaced pi-stacking"
    if interType in interIdByType:
        rules = []
        rules.append(Interaction.inter_type_id == interIdByType[interType])
        rules.append(Interaction.cc_dist_pi_pi_inter <= interConf.max_cc_dist_pi_pi_inter)
        rules.append(Interaction.dihed_ang_pi_pi_inter < interConf.min_dihed_ang_pi_pi_inter)
        rules.append(Interaction.disp_ang_pi_pi_inter > interConf.max_disp_ang_pi_pi_inter)
        filters.append(FilterRules(interType, rules))

    ##################################################
    interType = "Hydrophobic"
    if interType in interIdByType:
        rules = []
        rules.append(Interaction.inter_type_id == interIdByType[interType])
        rules.append(Interaction.dist_hydrop_inter <= interConf.max_dist_hydrop_inter)
        filters.append(FilterRules(interType, rules))

    ##################################################
    interType = "Cation-pi"
    if interType in interIdByType:
        rules = []
        rules.append(Interaction.inter_type_id == interIdByType[interType])
        rules.append(Interaction.dist_cation_pi_inter <= interConf.max_dist_cation_pi_inter)
        filters.append(FilterRules(interType, rules))

    ##################################################
    interType = "Halogen bond"
    if interType in interIdByType:
        # Halogen bond involving a PI system (aromatic ring)
        rules = []
        rules.append(Interaction.inter_type_id == interIdByType[interType])
        rules.append(Interaction.xc_dist_xbond_inter <= interConf.max_xc_dist_xbond_inter)
        rules.append(Interaction.disp_ang_xbond_inter <= interConf.max_disp_ang_xbond_inter)
        rules.append(Interaction.cxa_ang_xbond_inter >= interConf.min_cxa_ang_xbond_inter)
        filters.append(FilterRules(interType, rules))

        # Halogen bond
        rules = []
        rules.append(Interaction.inter_type_id == interIdByType[interType])
        rules.append(Interaction.xa_dist_xbond_inter <= interConf.max_xa_dist_xbond_inter)
        rules.append(Interaction.cxa_ang_xbond_inter >= interConf.min_cxa_ang_xbond_inter)
        rules.append(Interaction.xar_ang_xbond_inter >= interConf.min_xar_ang_xbond_inter)
        filters.append(FilterRules(interType, rules))

    ##################################################
    interType = "Water-bridged hydrogen bond"
    if interType in interIdByType:
        rules = []
        rules.append(Interaction.inter_type_id == interIdByType[interType])
        filters.append(FilterRules(interType, rules))

    ##################################################
    interType = "Salt bridge"
    if interType in interIdByType:
        rules = []
        rules.append(Interaction.inter_type_id == interIdByType[interType])
        filters.append(FilterRules(interType, rules))

    return filters


def format_db_ligand_entries(db_ligand_entries):
    entries = []
    for db_entry in db_ligand_entries:
        new_entry = DBEntry(db_entry.id, db_entry.pdb_id, db_entry.chain_id, db_entry.lig_name, db_entry.lig_num)
        entries.append(new_entry)

    return entries


def format_db_atoms(mybio_residue, db_atoms):
    atoms = []
    for db_atom in db_atoms:
        mybio_atom = mybio_residue[db_atom.name]
        atoms.append(NbAtom(mybio_atom, []))

    return atoms


def format_db_chemical_features(db_chem_features):
    chem_features = []
    for chem_feature in db_chem_features:
        chem_features.append(ChemicalFeature(chem_feature.type))

    return chem_features


def format_db_group(mybio_entity, db_group):
    chain_id = db_group.comp_chain_id
    comp_name = db_group.comp_name
    comp_num = db_group.comp_num
    comp_icode = db_group.comp_icode

    # TODO: Pode ser um acido nucleico tb!!!!
    if is_aa(comp_name):
        hetfield = " "
    elif comp_name == "HOH":
        hetfield = "W"
    else:
        hetfield = "H_%s" % comp_name

    if not comp_icode:
        comp_icode = " "

    comp_key = (hetfield, comp_num, comp_icode)
    mybio_residue = mybio_entity[chain_id][comp_key]

    atoms = format_db_atoms(mybio_residue, db_group.atoms)
    chem_features = format_db_chemical_features(db_group.compound_types)

    atom_group = AtomGroup(atoms, chem_features)

    return atom_group


def format_db_interactions(mybio_entity, db_interactions):
    db_sorted_inter = sorted(db_interactions, key=lambda x: x.id)

    if mybio_entity.level != "S":
        mybio_entity = mybio_entity.get_parent_by_level("S")
    mybio_entity = mybio_entity[0]

    interactions = []
    for db_inter in db_sorted_inter:
        atom_group1 = format_db_group(mybio_entity, db_inter.group1)
        atom_group2 = format_db_group(mybio_entity, db_inter.group2)

        params = {}
        cols = list(inspect(db_inter).mapper.column_attrs)[4:]
        for col in cols:
            value = getattr(db_inter, col.key)
            if value:
                params[col.key] = value

        interactions.append(InteractionType(atom_group1,
                                            atom_group2,
                                            db_inter.inter_type.type,
                                            params))

    return interactions


def get_default_mappers_list(db):
        mapper_list = []
        mapper_list.append(MapperConfiguration(CompTypeCount, "comp_type_count"))
        mapper_list.append(MapperConfiguration(InterTypeCount, "inter_type_count"))
        mapper_list.append(MapperConfiguration(Status, "status"))

        ##################################################################
        # Project related tables
        ##################################################################
        mapper_list.append(MapperConfiguration(Project, "project"))
        mapper_list.append(MapperConfiguration(ProjectStep, "proj_step"))
        mapper_list.append(MapperConfiguration(ProjectType, "project_type",
                                               properties={"steps": relationship(ProjectStep)}))
        mapper_list.append(MapperConfiguration(ProjectStepDetail,
                                               "proj_step_detail"))
        mapper_list.append(MapperConfiguration(ProjectStepMessage,
                                               "proj_step_message"))
        mapper_list.append(MapperConfiguration(InterTypeFreqByCluster,
                                               "inter_type_freq_by_cluster"))
        mapper_list.append(MapperConfiguration(CompTypeFreqByCluster,
                                               "comp_type_freq_by_cluster"))

        ##################################################################
        # Interaction related tables
        ##################################################################
        mapper_list.append(MapperConfiguration(InterType, "inter_type"))
        mapper_list.append(MapperConfiguration(InterDependInter,
                                               "inter_depend_inter"))
        mapper_list.append(MapperConfiguration(GroupCompoundType,
                                               "group_compound_type"))
        mapper_list.append(MapperConfiguration(CompoundType, "compound_type"))
        mapper_list.append(MapperConfiguration(GroupAtom, "group_atom"))
        mapper_list.append(MapperConfiguration(Atom, "atom"))

        grp_comp_type_tbl = db.get_table("group_compound_type")
        grp_atom_tbl = db.get_table("group_atom")
        grp_props = {
            "compound_types": relationship(CompoundType,
                                           secondary=grp_comp_type_tbl,
                                           backref='group'),
            "atoms": relationship(Atom,
                                  secondary=grp_atom_tbl,
                                  backref='group')
        }
        mapper_list.append(MapperConfiguration(Group, "group",
                                               properties=grp_props))

        inter_inter_tbl = db.get_table("inter_depend_inter")
        inter_tbl = db.get_table("interaction")

        ##################################################################
        # RCSB related tables
        ##################################################################
        mapper_list.append(MapperConfiguration(ExperimentalTechnique, "experimental_tech"))
        struct_props = {"experimental_tech": relationship(ExperimentalTechnique)}
        mapper_list.append(MapperConfiguration(Structure, "structure",
                                               properties=struct_props))

        mapper_list.append(MapperConfiguration(LigandEstablishInteraction,
                                               "ligand_establish_interaction"))
        ligand_inter_tbl = db.get_table("ligand_establish_interaction")
        group_tbl = db.get_table("group")
        inter_props = {
            "dependencies": relationship(RCSBInteraction,
                                         primaryjoin=(inter_tbl.c.id ==
                                                      inter_inter_tbl.c.interaction_id1),
                                         secondaryjoin=(inter_tbl.c.id ==
                                                        inter_inter_tbl.c.interaction_id2),
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
        mapper_list.append(MapperConfiguration(RCSBInteraction, inter_tbl,
                                               properties=inter_props))

        ligand_tbl = db.get_table("ligand")
        ligand_props = {
            "interactions": relationship(RCSBInteraction,
                                         secondary=ligand_inter_tbl,
                                         backref='ligand')
        }
        mapper_list.append(MapperConfiguration(Ligand, ligand_tbl,
                                               properties=ligand_props))

        # ##################################################################
        # User project related tables
        # ##################################################################
        mapper_list.append(MapperConfiguration(LigandEntryHasStepMessage, "lig_entry_has_step_msg"))

        lig_entry_msg_tbl = db.get_table("lig_entry_has_step_msg")
        ligand_entry_tbl = db.get_table("ligand_entry")
        lig_entry_props = {
            "step_messages": relationship(ProjectStepMessage,
                                          secondary=lig_entry_msg_tbl,
                                          backref="ligand_entry")
        }
        mapper_list.append(MapperConfiguration(LigandEntry, ligand_entry_tbl,
                                               properties=lig_entry_props))
        mapper_list.append(MapperConfiguration(LigandEntryHasInteractions,
                                               "ligand_entry_has_interactions"))

        ligand_inter_tbl = db.get_table("ligand_entry_has_interactions")
        inter_props = {
            "dependencies": relationship(ProjectInteraction,
                                         primaryjoin=(inter_tbl.c.id ==
                                                      inter_inter_tbl.c.interaction_id1),
                                         secondaryjoin=(inter_tbl.c.id ==
                                                        inter_inter_tbl.c.interaction_id2),
                                         secondary=inter_inter_tbl,
                                         backref='interaction'),
            "ligand_entities": relationship(LigandEntry,
                                            secondary=ligand_inter_tbl,
                                            backref='interaction')
        }
        mapper_list.append(MapperConfiguration(ProjectInteraction, inter_tbl,
                                               properties=inter_props))

        return mapper_list
