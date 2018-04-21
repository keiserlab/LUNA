from sqlalchemy import inspect
from database.napoli_model import Interaction
from database.filters import FilterRules

from input.complex import DBComplex


def object_as_dict(obj):

    return {c.key: getattr(obj, c.key)
            for c in inspect(obj).mapper.column_attrs}


def rows_2dict(rows):
    return [object_as_dict(r) for r in rows]


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
        rules.append(Interaction.ha_dist_hb_inter == -1)
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


def prepare_complex_entries(dbComplexes):
    complexes = []
    for dbComplex in dbComplexes:
        newComplex = DBComplex(dbComplex.id,
                               dbComplex.pdb,
                               dbComplex.chain,
                               dbComplex.lig_name,
                               dbComplex.lig_num)
        complexes.append(newComplex)

    return complexes
