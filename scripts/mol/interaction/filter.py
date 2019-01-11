from mol.interaction.conf import DefaultInteractionConf
from operator import (le, ge)


class InteractionFilter:

    def __init__(self, funcs=None, inter_conf=DefaultInteractionConf(),
                 ignore_self_inter=True, ignore_prot_prot=True, ignore_prot_nucl=True,
                 ignore_prot_lig=True, ignore_nucl_nucl=True, ignore_nucl_lig=True,
                 ignore_lig_lig=True, ignore_h2o_h2o=True, ignore_any_h2o=False,
                 ignore_multi_comps=False, ignore_mixed_class=False):

        if funcs is None:
            funcs = self._default_functions()
        self._funcs = funcs

        self.inter_conf = inter_conf

        self.ignore_self_inter = ignore_self_inter
        self.ignore_prot_prot = ignore_prot_prot
        self.ignore_prot_nucl = ignore_prot_nucl
        self.ignore_prot_lig = ignore_prot_lig
        self.ignore_nucl_nucl = ignore_nucl_nucl
        self.ignore_nucl_lig = ignore_nucl_lig
        self.ignore_lig_lig = ignore_lig_lig
        self.ignore_h2o_h2o = ignore_h2o_h2o
        self.ignore_any_h2o = ignore_any_h2o
        self.ignore_multi_comps = ignore_multi_comps
        self.ignore_mixed_class = ignore_mixed_class

    @classmethod
    def new_pli_filter(cls, ignore_prot_lig=False, ignore_any_h2o=False, **kwargs):
        return cls(ignore_prot_lig=ignore_prot_lig, ignore_any_h2o=ignore_any_h2o, **kwargs)

    @classmethod
    def new_ppi_filter(cls, ignore_prot_prot=False, ignore_any_h2o=False, **kwargs):
        return cls(ignore_prot_prot=ignore_prot_prot, ignore_any_h2o=ignore_any_h2o, **kwargs)

    @classmethod
    def new_pni_filter(cls, ignore_prot_nucl=False, ignore_any_h2o=False, **kwargs):
        return cls(ignore_prot_nucl=ignore_prot_nucl, ignore_any_h2o=ignore_any_h2o, **kwargs)

    @classmethod
    def new_nni_filter(cls, ignore_nucl_nucl=False, ignore_any_h2o=False, **kwargs):
        return cls(ignore_nucl_nucl=ignore_nucl_nucl, ignore_any_h2o=ignore_any_h2o, **kwargs)

    @classmethod
    def new_nli_filter(cls, ignore_nucl_lig=False, ignore_any_h2o=False, **kwargs):
        return cls(ignore_nucl_lig=ignore_nucl_lig, ignore_any_h2o=ignore_any_h2o, **kwargs)

    @property
    def funcs(self):
        return self._funcs

    def _default_functions(self):
        return {"Hydrogen bond": self.filter_hbond,
                "Attractive": self.filter_attractive,
                "Cation-pi": self.filter_cation_pi,
                "Edge-to-face pi-stacking": self.filter_pi_pi,
                "Face-to-face pi-stacking": self.filter_pi_pi,
                "Parallel-displaced pi-stacking": self.filter_pi_pi,
                "Pi-stacking": self.filter_pi_pi,
                "Hydrophobic": self.filter_hydrop,
                "Halogen bond": self.filter_xbond,
                "Repulsive": self.filter_repulsive}

    def is_valid_pair(self, atm_grp1, atm_grp2):

        # It will always ignore interactions involving the same atom groups.
        # Loops in the graph is not permitted and does not make any sense.
        if atm_grp1 == atm_grp2:
            return False

        # It will always ignore interactions involving atoms and the group to which they belong to.
        # For example, the centroid of an aromatic ring cannot interact with an atom
        # that belongs to the ring. It is a type of Loop.
        if atm_grp1.contain_group(atm_grp2) or atm_grp2.contain_group(atm_grp1):
            return False

        # If one of the groups contain atoms from different compounds.
        has_multi_comps = (len(atm_grp1.compounds) > 1 or len(atm_grp2.compounds) > 1)
        if self.ignore_multi_comps and has_multi_comps:
            return False

        # If one of the groups contain compounds from different classes as, for instance, residue and ligand.
        # It means that compounds from different classes are covalently bonded to each other.
        has_any_mixed = (atm_grp1.is_mixed() or atm_grp2.is_mixed())
        if self.ignore_mixed_class and has_any_mixed:
            return False

        # It ignores interactions involving the same compounds if it is required.
        # As each group may have atoms from different compounds, we can check if there is at least
        # one common compound between the two groups. Remember that if two or more compounds exist in a group,
        # it means that these molecules are covalently bonded and should be considered the same molecule.
        # For example: a carbohydrate can be expressed in a PDB as its subparts:
        # GLC + GLC = CBI
        is_same_compounds = len(atm_grp1.compounds.intersection(atm_grp2.compounds)) >= 1
        if self.ignore_self_inter and is_same_compounds:
            return False
        # It accepts all pairs composed by atom groups from a same compound when
        # ignore_self_inter is set off.
        elif not self.ignore_self_inter and is_same_compounds:
            return True

        # It ignores protein-protein interactions if it is required.
        is_prot_prot = (atm_grp1.is_aminoacid() and atm_grp2.is_aminoacid())
        if self.ignore_prot_prot and is_prot_prot:
            return False

        # It ignores protein-nucleic acid interactions if it is required.
        is_prot_nucl = ((atm_grp1.is_aminoacid() and atm_grp2.is_nucleotide()) or
                        (atm_grp1.is_nucleotide() and atm_grp2.is_aminoacid()))
        if self.ignore_prot_nucl and is_prot_nucl:
            return False

        # It ignores protein-ligand interactions if it is required.
        is_prot_lig = ((atm_grp1.is_aminoacid() and atm_grp2.is_hetatm()) or
                       (atm_grp1.is_hetatm() and atm_grp2.is_aminoacid()))
        if self.ignore_prot_lig and is_prot_lig:
            return False

        # It ignores nucleic acid-nucleic acid interactions if it is required.
        is_nucl_nucl = (atm_grp1.is_nucleotide() and atm_grp2.is_nucleotide())
        if self.ignore_nucl_nucl and is_nucl_nucl:
            return False

        # It ignores nucleic acid-ligand interactions if it is required.
        is_nucl_lig = ((atm_grp1.is_nucleotide() and atm_grp2.is_hetatm()) or
                       (atm_grp1.is_hetatm() and atm_grp2.is_nucleotide()))
        if self.ignore_nucl_lig and is_nucl_lig:
            return False

        # It ignores ligand-ligand interactions if it is required.
        is_lig_lig = (atm_grp1.is_hetatm() and atm_grp2.is_hetatm())
        if self.ignore_lig_lig and is_lig_lig:
            return False

        # It ignores interactions of other molecule types with water.
        # It enables the possibility of identifying water-bridged interaction.
        # Eg: residue-water, ligand-water = residue -- water -- ligand.
        is_any_h2o = (atm_grp1.is_water() or atm_grp2.is_water())
        if self.ignore_any_h2o and is_any_h2o:
            return False

        # It ignores interactions involving two waters if it is required.
        # If it is on, it will produce water-bridged interactions of multiple levels
        # Eg: residue -- h2o -- h2o -- ligand, residue -- residue -- h2o -- h2o -- ligand.
        is_h2o_h2o = (atm_grp1.is_water() and atm_grp2.is_water())
        if self.ignore_h2o_h2o and is_h2o_h2o:
            return False

        return True

    def is_valid_interaction(self, interaction, validate_pair=True):
        if validate_pair:
            if not self.is_valid_pair(interaction.atm_grp1, interaction.atm_grp2):
                return False

        if (interaction.type not in self.funcs):
            return True
        else:
            func = self.funcs[interaction.type]
            return func(interaction)

    def filter_hbond(self, interaction):
        is_valid = False

        if self.is_within_boundary(interaction.da_dist_hb_inter, "max_da_dist_hb_inter", le):

            if (interaction.ha_dist_hb_inter == -1 and interaction.dha_ang_hb_inter == -1):
                # If the angle DHA and the distance HA is equal -1, it is
                # expected that the molecules to be water molecules. "
                if (interaction.atm_grp1.compound.is_water() or interaction.atm_grp2.compound.is_water()):
                    ha_dist = interaction.da_dist_hb_inter - 1
                    if (self.is_within_boundary(ha_dist, "max_ha_dist_hb_inter", le)):
                        is_valid = True
            else:
                if (self.is_within_boundary(interaction.ha_dist_hb_inter, "max_ha_dist_hb_inter", le) and
                        self.is_within_boundary(interaction.dha_ang_hb_inter, "min_dha_ang_hb_inter", ge)):
                        is_valid = True
        return is_valid

    def filter_attractive(self, interaction):
        return self.is_within_boundary(interaction.dist_attract_inter, "max_dist_attract_inter", le)

    def filter_cation_pi(self, interaction):
        return self.is_within_boundary(interaction.dist_cation_pi_inter, "max_dist_cation_pi_inter", le)

    def filter_pi_pi(self, interaction):
        is_valid = False
        if self.is_within_boundary(interaction.cc_dist_pi_pi_inter, "max_cc_dist_pi_pi_inter", le):
            # If the angle criteria were defined, test if the interaction fits the requirements.
            if ("min_dihed_ang_pi_pi_inter" in self.inter_conf.conf and
                    "max_disp_ang_pi_pi_inter" in self.inter_conf.conf):
                if self.is_within_boundary(interaction.dihed_ang_pi_pi_inter, "min_dihed_ang_pi_pi_inter", ge):
                    # It overwrites a general Pi-stacking interaction.
                    interaction.type = "Edge-to-face pi-stacking"
                elif self.is_within_boundary(interaction.disp_ang_pi_pi_inter, "max_disp_ang_pi_pi_inter", le):
                    # It overwrites a general Pi-stacking interaction.
                    interaction.type = "Face-to-face pi-stacking"
                else:
                    # It overwrites a general Pi-stacking interaction.
                    interaction.type = "Parallel-displaced pi-stacking"
            else:
                interaction.type = "Pi-stacking"

            # It requires that at least the distance fits the criteria to an interaction to be considered valid.
            # For example, if the angle criteria were not defined, it is not possible to attribute an explicit type
            # of aromatic stacking, however the interaction should be considered valid.
            is_valid = True

        return is_valid

    def filter_hydrop(self, interaction):
        return self.is_within_boundary(interaction.dist_hydrop_inter, "max_dist_hydrop_inter", le)

    def filter_xbond(self, interaction):
        # Halogen bond involving a PI system (aromatic ring)
        if (interaction.atm_grp1 == "Aromatic" or interaction.atm_grp2 == "Aromatic"):
            return (self.is_within_boundary(interaction.xc_dist_xbond_inter, "max_xc_dist_xbond_inter", le) and
                    self.is_within_boundary(interaction.disp_ang_xbond_inter, "max_disp_ang_xbond_inter", le) and
                    self.is_within_boundary(interaction.cxa_ang_xbond_inter, "min_cxa_ang_xbond_inter", ge))
        # Classical halogen bond, i.e., does not envolve a PI system
        else:
            return (self.is_within_boundary(interaction.xa_dist_xbond_inter, "max_xa_dist_xbond_inter", le) and
                    self.is_within_boundary(interaction.cxa_ang_xbond_inter, "min_cxa_ang_xbond_inter", ge) and
                    self.is_within_boundary(interaction.xar_ang_xbond_inter, "min_xar_ang_xbond_inter", ge))

    def filter_repulsive(self, interaction):
        return self.is_within_boundary(interaction.dist_repuls_inter, "max_dist_repuls_inter", le)

    def is_within_boundary(self, value, key, func):
        if key not in self.inter_conf.conf:
            return True
        else:
            return func(value, self.inter_conf.conf[key])
