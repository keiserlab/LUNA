import re
from ast import literal_eval

from luna.mol.entry import REGEX_RESNUM_ICODE
from luna.mol.interaction.conf import DefaultInteractionConf
from luna.util.config_parser import Config
from luna.util.exceptions import IllegalArgumentError


AROMATIC_STACKINGS = ["pi-stacking", "face-to-face pi-stacking", "face-to-edge pi-stacking", "face-to-slope pi-stacking",
                      "edge-to-edge pi-stacking", "edge-to-face pi-stacking", "edge-to-slope pi-stacking", "displaced face-to-face pi-stacking",
                      "displaced face-to-edge pi-stacking", "displaced face-to-slope pi-stacking"]


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
    def new_pli_filter(cls, ignore_res_hetatm=False, ignore_hetatm_hetatm=False, ignore_any_h2o=False, ignore_self_inter=False, **kwargs):
        return cls(ignore_res_hetatm=ignore_res_hetatm, ignore_hetatm_hetatm=ignore_hetatm_hetatm, ignore_any_h2o=ignore_any_h2o,
                   ignore_self_inter=ignore_self_inter, **kwargs)

    @classmethod
    def new_ppi_filter(cls, ignore_res_res=False, ignore_inter_chain=False, ignore_intra_chain=False, ignore_any_h2o=False,
                       ignore_self_inter=False, **kwargs):

        return cls(ignore_res_res=ignore_res_res, ignore_inter_chain=ignore_inter_chain, ignore_intra_chain=ignore_intra_chain,
                   ignore_any_h2o=ignore_any_h2o, ignore_self_inter=ignore_self_inter, **kwargs)

    @classmethod
    def new_pni_filter(cls, ignore_res_nucl=False, ignore_inter_chain=False, ignore_intra_chain=False, ignore_any_h2o=False,
                       ignore_self_inter=False, **kwargs):

        return cls(ignore_res_nucl=ignore_res_nucl, ignore_inter_chain=ignore_inter_chain, ignore_intra_chain=ignore_intra_chain,
                   ignore_any_h2o=ignore_any_h2o, ignore_self_inter=ignore_self_inter, **kwargs)

    @classmethod
    def new_nni_filter(cls, ignore_nucl_nucl=False, ignore_inter_chain=False, ignore_intra_chain=False,
                       ignore_any_h2o=False, ignore_self_inter=False, **kwargs):

        return cls(ignore_nucl_nucl=ignore_nucl_nucl, ignore_inter_chain=ignore_inter_chain, ignore_intra_chain=ignore_intra_chain,
                   ignore_any_h2o=ignore_any_h2o, ignore_self_inter=ignore_self_inter, **kwargs)

    @classmethod
    def new_nli_filter(cls, ignore_nucl_hetatm=False, ignore_hetatm_hetatm=False, ignore_any_h2o=False, ignore_self_inter=False, **kwargs):
        return cls(ignore_nucl_hetatm=ignore_nucl_hetatm, ignore_hetatm_hetatm=ignore_hetatm_hetatm, ignore_any_h2o=ignore_any_h2o,
                   ignore_self_inter=ignore_self_inter, **kwargs)

    @property
    def funcs(self):
        return self._funcs

    def _default_functions(self):
        return {}

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


class BindingModeCondition:

    def __init__(self, condition_str):
        self.accept_all = False
        self.accept_all_chains = False
        self.accept_all_comps = False
        self.accept_all_atoms = False
        self.accept_all_comp_nums = False

        self.chain = None
        self.comp_name = None
        self.comp_num = None
        self.comp_icode = None
        self.atom = None

        self._condition_repr = condition_str

        self._parse_condition(condition_str.upper())

    def _parse_condition(self, condition_str):
        # Accept everything.
        if condition_str == "*":
            self.accept_all = True
        else:
            chain, comp_name, comp_num, atom = condition_str.split("/")

            if chain == "*":
                self.accept_all_chains = True
            else:
                self.chain = chain

            if comp_name == "*":
                self.accept_all_comps = True
            else:
                self.comp_name = comp_name

            if comp_num == "*":
                self.accept_all_comp_nums = True
            else:
                # Separate ligand number from insertion code.
                matched = REGEX_RESNUM_ICODE.match(comp_num)
                if matched:
                    comp_num = matched.group(1)
                    try:
                        assert float(comp_num).is_integer()
                        comp_num = int(comp_num)
                    except (ValueError, AssertionError):
                        raise IllegalArgumentError("The informed compound number '%s' is invalid. It must be an integer." % str(comp_num))

                    icode = None if matched.group(2) == "" else matched.group(2)
                else:
                    raise IllegalArgumentError("The compound number and its insertion code (if applicable) '%s' is invalid. "
                                               "It must be an integer followed by one insertion code character when applicable."
                                               % comp_num)

                self.comp_num = comp_num
                self.comp_icode = icode

            if atom == "*":
                self.accept_all_atoms = True
            else:
                self.atom = atom

    def is_valid(self, atm_grp):

        # Accept everything.
        if self.accept_all:
            return True

        # Accept everything.
        if self.accept_all_chains and self.accept_all_comps and self.accept_all_comp_nums and self.accept_all_atoms:
            return True

        # Tries to identify the first valid compound in the atom group.
        for atm in atm_grp.atoms:
            comp = atm.parent

            is_chain_valid, is_comp_valid, is_comp_num_valid, is_atom_valid = (self.accept_all_chains, self.accept_all_comps,
                                                                               self.accept_all_comp_nums, self.accept_all_atoms)

            # print("\t => ", comp, atm.name)
            # print("\t     - Initial: ", is_chain_valid, is_comp_valid, is_comp_num_valid, is_atom_valid)

            if self.chain is not None and self.chain == comp.parent.id:
                is_chain_valid = True

            if self.comp_name is not None and self.comp_name == comp.resname:
                is_comp_valid = True

            if self.comp_num is not None and self.comp_num == comp.id[1]:
                icode = None if comp.id[2].strip() == "" else comp.id[2]

                if self.comp_icode == icode:
                    is_comp_num_valid = True

            if self.atom is not None:

                # Verify element equality.
                if self.atom.endswith("*"):
                    elem = self.atom.rstrip("*")

                    if elem == atm.element:
                        is_atom_valid = True

                # Verify atom name equality.
                elif self.atom == atm.name:
                    is_atom_valid = True

            # print("\t     -   Final: ", is_chain_valid, is_comp_valid, is_comp_num_valid, is_atom_valid)
            # print()

            if is_chain_valid and is_comp_valid and is_comp_num_valid and is_atom_valid:
                return True

        return False

    def __repr__(self):
        return "<BindingModeCondition: %s" % self._condition_repr


class BindingModeFilter:

    def __init__(self, config):
        self.config = config

    @classmethod
    def from_config_file(cls, config_file):
        filtering_config = {}

        config = Config(config_file)

        for inter_type in config.sections():
            params = config.get_section_map(inter_type)

            inter_type = inter_type.lower()

            accept_all = False

            values = []
            if "accept_all" in params:
                accept_all = literal_eval(params["accept_all"])
                if accept_all is True:
                    values = ["*"]

            elif "accept_only" in params and not accept_all:
                values = literal_eval(params["accept_only"])

                if "*" in values:
                    values = ["*"]

            conditions = [BindingModeCondition(condition) for condition in values]
            filtering_config[inter_type] = conditions

        return cls(filtering_config)

    def is_valid(self, inter):

        if inter.type.lower() in self.config:
            inter_type = inter.type.lower()

        elif "aromatic stacking" in self.config and inter.type.lower() in AROMATIC_STACKINGS:
            inter_type = "aromatic stacking"

        elif "*" in self.config:
            inter_type = "*"

        else:
            return False

        for condition in self.config[inter_type]:
            is_src_grp_valid = condition.is_valid(inter.src_grp)
            is_trgt_grp_valid = condition.is_valid(inter.trgt_grp)

            if is_src_grp_valid or is_trgt_grp_valid:
                return True

        return False
