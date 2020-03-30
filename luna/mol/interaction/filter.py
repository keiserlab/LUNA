from operator import le, ge

from luna.mol.interaction.conf import DefaultInteractionConf


class InteractionFilter:

    def __init__(self, funcs=None, inter_conf=DefaultInteractionConf(),
                 ignore_self_inter=True, ignore_intra_chain=True, ignore_inter_chain=True,
                 ignore_res_res=True, ignore_res_nucl=True, ignore_res_hetatm=True,
                 ignore_nucl_nucl=True, ignore_nucl_hetatm=True, ignore_hetatm_hetatm=True,
                 ignore_h2o_h2o=True, ignore_any_h2o=False, ignore_multi_comps=False, ignore_mixed_class=False):

        if funcs is None:
            funcs = self._default_functions()

        self._funcs = funcs

        self.inter_conf = inter_conf

        self.ignore_self_inter = ignore_self_inter
        self.ignore_intra_chain = ignore_intra_chain
        self.ignore_inter_chain = ignore_inter_chain
        self.ignore_res_res = ignore_res_res
        self.ignore_res_nucl = ignore_res_nucl
        self.ignore_res_hetatm = ignore_res_hetatm
        self.ignore_nucl_nucl = ignore_nucl_nucl
        self.ignore_nucl_hetatm = ignore_nucl_hetatm
        self.ignore_hetatm_hetatm = ignore_hetatm_hetatm
        self.ignore_h2o_h2o = ignore_h2o_h2o
        self.ignore_any_h2o = ignore_any_h2o
        self.ignore_multi_comps = ignore_multi_comps
        self.ignore_mixed_class = ignore_mixed_class

    @classmethod
    def new_pli_filter(cls, ignore_res_hetatm=False, ignore_any_h2o=False, ignore_self_inter=False, **kwargs):
        return cls(ignore_res_hetatm=ignore_res_hetatm, ignore_any_h2o=ignore_any_h2o, ignore_self_inter=ignore_self_inter, **kwargs)

    @classmethod
    def new_ppi_filter(cls, ignore_res_res=False, ignore_inter_chain=False, ignore_any_h2o=False, ignore_self_inter=False, **kwargs):
        return cls(ignore_res_res=ignore_res_res, ignore_any_h2o=ignore_any_h2o, ignore_self_inter=ignore_self_inter, **kwargs)

    @classmethod
    def new_pni_filter(cls, ignore_res_nucl=False, ignore_any_h2o=False, ignore_self_inter=False, **kwargs):
        return cls(ignore_res_nucl=ignore_res_nucl, ignore_any_h2o=ignore_any_h2o, ignore_self_inter=ignore_self_inter, **kwargs)

    @classmethod
    def new_nni_filter(cls, ignore_nucl_nucl=False, ignore_any_h2o=False, ignore_self_inter=False, **kwargs):
        return cls(ignore_nucl_nucl=ignore_nucl_nucl, ignore_any_h2o=ignore_any_h2o, ignore_self_inter=ignore_self_inter, **kwargs)

    @classmethod
    def new_nli_filter(cls, ignore_nucl_hetatm=False, ignore_any_h2o=False, ignore_self_inter=False, **kwargs):
        return cls(ignore_nucl_hetatm=ignore_nucl_hetatm, ignore_any_h2o=ignore_any_h2o, ignore_self_inter=ignore_self_inter, **kwargs)

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

    def is_valid_pair(self, src_grp, trgt_grp):

        # It will always ignore interactions involving the same atom groups.
        # Loops in the graph is not permitted and does not make any sense.
        if src_grp == trgt_grp:
            return False

        # It will always ignore interactions involving atoms and the group to which they belong to.
        # For example, the centroid of an aromatic ring cannot interact with an atom
        # that belongs to the ring. It is a type of Loop.
        if src_grp.contain_group(trgt_grp) or trgt_grp.contain_group(src_grp):
            return False

        # If one of the groups contain atoms from different compounds.
        has_multi_comps = (len(src_grp.compounds) > 1 or len(trgt_grp.compounds) > 1)
        if self.ignore_multi_comps and has_multi_comps:
            return False

        # If one of the groups contain compounds from different classes as, for instance, residue and ligand.
        # It means that compounds from different classes are covalently bonded to each other.
        has_any_mixed = (src_grp.is_mixed() or trgt_grp.is_mixed())
        if self.ignore_mixed_class and has_any_mixed:
            return False

        # It ignores interactions involving the same compounds if required.
        # As each group may have atoms from different compounds, we can check if there is at least
        # one common compound between the two groups. Remember that if two or more compounds exist in a group,
        # it means that these molecules are covalently bonded and should be considered the same molecule.
        # For example: a carbohydrate can be expressed in a PDB as its subparts:
        #      E.g.: GLC + GLC = CBI
        # The same applies to any group formed after covalently bonding a residue to a hetatm (ligand or non-standard amino acid
        # represented as hetatm)
        is_same_compounds = len(src_grp.compounds.intersection(trgt_grp.compounds)) >= 1
        if self.ignore_self_inter and is_same_compounds:
            return False

        # Check if two groups contain the same chains and if both of them contain only one chain. The second condition removes
        # groups containing residues of different chains as may occur due to disulfide bonds.
        # Note, however, that this flag will be used only as a filter for intra-interactions in protein/RNA/DNA chains.
        is_same_chain = src_grp.get_chains() == trgt_grp.get_chains() and len(src_grp.get_chains()) == 1

        # Filters for residue-residue interactions if required.
        is_res_res = (src_grp.is_residue() and trgt_grp.is_residue())
        if is_res_res:
            # Ignore all residue-residue interactions.
            if self.ignore_res_res:
                return False
            # Ignore all intra-chain interactions involving two residues.
            elif self.ignore_intra_chain and is_same_chain:
                return False
            elif self.ignore_inter_chain and not is_same_chain:
                return False

        # It ignores residue-nucleic acid interactions if required.
        is_res_nucl = ((src_grp.is_residue() and trgt_grp.is_nucleotide()) or
                       (src_grp.is_nucleotide() and trgt_grp.is_residue()))
        if self.ignore_res_nucl and is_res_nucl:
            return False

        # It ignores residue-ligand interactions if required.
        is_res_hetatm = ((src_grp.is_residue() and trgt_grp.is_hetatm()) or
                         (src_grp.is_hetatm() and trgt_grp.is_residue()))
        if self.ignore_res_hetatm and is_res_hetatm:
            return False

        # Filters for nucleic acid-nucleic acid interactions if required.
        is_nucl_nucl = (src_grp.is_nucleotide() and trgt_grp.is_nucleotide())
        if is_nucl_nucl:
            # Ignore all nucleic acid-nucleic acid interactions
            if self.ignore_nucl_nucl:
                return False
            # Ignore all intra-chain interactions involving two nucleic acids (RNA/DNA chains).
            elif self.ignore_intra_chain and is_same_chain:
                return False
            # Ignore all inter-chain interactions involving two nucleic acids (RNA/DNA chains).
            elif self.ignore_inter_chain and not is_same_chain:
                return False

        # It ignores nucleic acid-ligand interactions if required.
        is_nucl_hetatm = ((src_grp.is_nucleotide() and trgt_grp.is_hetatm()) or
                          (src_grp.is_hetatm() and trgt_grp.is_nucleotide()))
        if self.ignore_nucl_hetatm and is_nucl_hetatm:
            return False

        # It ignores ligand-ligand interactions if required.
        is_hetatm_hetatm = (src_grp.is_hetatm() and trgt_grp.is_hetatm())
        if self.ignore_hetatm_hetatm and is_hetatm_hetatm:
            return False

        # It ignores interactions of other molecule types with water.
        # It enables the possibility of identifying water-bridged interaction.
        # Eg: residue-water, ligand-water = residue -- water -- ligand.
        is_any_h2o = (src_grp.is_water() or trgt_grp.is_water())
        if self.ignore_any_h2o and is_any_h2o:
            return False

        # It ignores interactions involving two waters if required.
        # if on, it will produce water-bridged interactions of multiple levels
        # Eg: residue -- h2o -- h2o -- ligand, residue -- residue -- h2o -- h2o -- ligand.
        is_h2o_h2o = (src_grp.is_water() and trgt_grp.is_water())
        if self.ignore_h2o_h2o and is_h2o_h2o:
            return False

        return True

    def is_valid_interaction(self, interaction, validate_pair=True):
        if validate_pair:
            if not self.is_valid_pair(interaction.src_grp, interaction.trgt_grp):
                return False

        if interaction.type not in self.funcs:
            return True
        else:
            func = self.funcs[interaction.type]
            return func(interaction)

    def filter_hbond(self, interaction):
        pass

    def filter_attractive(self, interaction):
        pass

    def filter_cation_pi(self, interaction):
        pass

    def filter_pi_pi(self, interaction):
        pass

    def filter_hydrop(self, interaction):
        pass

    def filter_xbond(self, interaction):
        pass

    def filter_repulsive(self, interaction):
        pass

    def is_within_boundary(self, value, key, func):
        if key not in self.inter_conf.conf:
            return True
        else:
            return func(value, self.inter_conf.conf[key])
