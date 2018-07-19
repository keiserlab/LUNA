class Ligand(object):
    pass


class LigandEstablishInteraction(object):
    pass


class Project(object):
    pass


class LigandEntry(object):

    def __init__(self, pdb_id, chain_id, lig_name=None, lig_num=None,
                 lig_icode=None, cluster=None, status_id=None,
                 process_error_id=None, inter_proj_id=None):

        self.pdb_id = pdb_id
        self.chain_id = chain_id
        self.lig_name = lig_name
        self.lig_num = lig_num
        self.lig_icode = lig_icode
        self.cluster = cluster
        self.status_id = status_id
        self.process_error_ir = process_error_id
        self.inter_proj_id = inter_proj_id


class LigandEntryHasInteractions(object):
    pass


class CompTypeCount(object):

    def __init__(self, ligand_entry_id, inter_proj_id,
                 compound_type_id, total):

        self.ligand_entry_id = ligand_entry_id
        self.inter_proj_id = inter_proj_id
        self.compound_type_id = compound_type_id
        self.total = total


class InterTypeCount(object):

    def __init__(self, ligand_entry_id, inter_proj_id,
                 inter_type_id, total):

        self.ligand_entry_id = ligand_entry_id
        self.inter_proj_id = inter_proj_id
        self.inter_type_id = inter_type_id
        self.total = total


class InterType(object):
    pass


class Interaction(object):

    def __init__(self, group_id1, group_id2, inter_type_id,
                 da_dist_hb_inter=None,
                 ha_dist_hb_inter=None, dha_ang_hb_inter=None,
                 dist_repuls_inter=None, dist_attract_inter=None,
                 cc_dist_pi_pi_inter=None, dihed_ang_pi_pi_inter=None,
                 disp_ang_pi_pi_inter=None, dist_hydrop_inter=None,
                 xa_dist_xbond_inter=None, xc_dist_xbond_inter=None,
                 cxa_ang_xbond_inter=None, xar_ang_xbond_inter=None,
                 disp_ang_xbond_inter=None, dist_cation_pi_inter=None,
                 cov_dist=None):

        self.group_id1 = group_id1
        self.group_id2 = group_id2
        self.inter_type_id = inter_type_id

        self.da_dist_hb_inter = da_dist_hb_inter
        self.ha_dist_hb_inter = ha_dist_hb_inter
        self.dha_ang_hb_inter = dha_ang_hb_inter

        self.dist_repuls_inter = dist_repuls_inter
        self.dist_attract_inter = dist_attract_inter

        self.cc_dist_pi_pi_inter = cc_dist_pi_pi_inter
        self.dihed_ang_pi_pi_inter = dihed_ang_pi_pi_inter
        self.disp_ang_pi_pi_inter = disp_ang_pi_pi_inter

        self.dist_hydrop_inter = dist_hydrop_inter

        self.xa_dist_xbond_inter = xa_dist_xbond_inter
        self.xc_dist_xbond_inter = xc_dist_xbond_inter
        self.cxa_ang_xbond_inter = cxa_ang_xbond_inter
        self.xar_ang_xbond_inter = xar_ang_xbond_inter
        self.disp_ang_xbond_inter = disp_ang_xbond_inter

        self.dist_cation_pi_inter = dist_cation_pi_inter

        self.cov_dist = cov_dist


class RCSBInteraction(Interaction):
    pass


class ProjectInteraction(Interaction):
    pass


class InterDependInter(object):
    pass


class CompoundType(object):
    pass


class Group(object):

    def __init__(self, comp_pdb_id, comp_chain_id, comp_name,
                 comp_num, comp_icode, name,
                 coord_x, coord_y, coord_z):

        self.comp_pdb_id = comp_pdb_id
        self.comp_chain_id = comp_chain_id
        self.comp_name = comp_name
        self.comp_num = comp_num
        self.comp_icode = comp_icode
        self.name = name
        self.coord_x = coord_x
        self.coord_y = coord_y
        self.coord_z = coord_z


class GroupCompoundType(object):
    pass


class GroupAtom(object):
    pass


class Atom(object):

    def __init__(self, name, serial_number, altloc, element,
                 coord_x, coord_y, coord_z):

        self.name = name
        self.serial_number = serial_number
        self.altloc = altloc
        self.element = element
        self.coord_x = coord_x
        self.coord_y = coord_y
        self.coord_z = coord_z


class Structure(object):
    pass


class ExperimentalTechnique(object):
    pass


class Status(object):
    pass


class ProjectStep(object):
    pass


class ProjectType(object):
    pass


class ProjectStepDetail(object):

    def __init__(self, project_id, proj_step_id, status_id=None,
                 warning=None, progress_perc=None):

        self.project_id = project_id
        self.proj_step_id = proj_step_id
        self.status_id = status_id
        self.warning = warning
        self.progress_perc = progress_perc


class ProjectStepMessage(object):

    def __init__(self, project_id, proj_step_id,
                 status_id=None, message=None, is_critical=None):
        self.project_id = project_id
        self.proj_step_id = proj_step_id
        self.status_id = status_id
        self.message = message

        # TODO: Verificar se vou remover esse parametro.
        self.is_critical = is_critical


class LigandEntryHasStepMessage(object):
    pass


class InterTypeFreqByCluster(object):

    def __init__(self, inter_proj_id, inter_type_id, cluster,
                 inter_count, ligand_count, cumulative_num):
        self.inter_proj_id = inter_proj_id
        self.inter_type_id = inter_type_id
        self.cluster = cluster
        self.inter_count = inter_count
        self.ligand_count = ligand_count
        self.cumulative_num = cumulative_num


class CompTypeFreqByCluster(object):

    def __init__(self, inter_proj_id, compound_type_id, cluster,
                 comp_count, ligand_count, cumulative_num):
        self.inter_proj_id = inter_proj_id
        self.compound_type_id = compound_type_id
        self.cluster = cluster
        self.comp_count = comp_count
        self.ligand_count = ligand_count
        self.cumulative_num = cumulative_num
