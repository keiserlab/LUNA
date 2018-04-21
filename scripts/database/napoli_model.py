class Project(object):
    pass


class Complex(object):
    pass


class ComplexInter(object):
    pass


class CompTypeCount(object):

    def __init__(self, complex_id, complex_inter_proj_id,
                 compound_type_id, total):

        self.complex_id = complex_id
        self.complex_inter_proj_id = complex_inter_proj_id
        self.compound_type_id = compound_type_id
        self.total = total


class InterTypeCount(object):

    def __init__(self, complex_id, complex_inter_proj_id,
                 inter_type_id, total):

        self.complex_id = complex_id
        self.complex_inter_proj_id = complex_inter_proj_id
        self.inter_type_id = inter_type_id
        self.total = total


class InterType(object):
    pass


class Interaction(object):

    def __init__(self, group_id1, group_id2, inter_type_id,
                 pdb_id=None, da_dist_hb_inter=None,
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
        self.pdb_id = pdb_id

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


class InterDependInter(object):
    pass


class CompoundType(object):
    pass


class Group(object):

    def __init__(self, orig_comp_pdb, orig_comp_chain, orig_comp_name,
                 orig_comp_num, orig_comp_icode, name,
                 coord_x, coord_y, coord_z):

        self.orig_comp_pdb = orig_comp_pdb
        self.orig_comp_chain = orig_comp_chain
        self.orig_comp_name = orig_comp_name
        self.orig_comp_num = orig_comp_num
        self.orig_comp_icode = orig_comp_icode
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


class PDB(object):
    pass
