from os import path
from ast import literal_eval
from collections import defaultdict

from luna.mol.entry import REGEX_RESNUM_ICODE
from luna.util.config import Config
from luna.util.exceptions import IllegalArgumentError


AROMATIC_STACKINGS = ["pi-stacking",
                      "face-to-face pi-stacking",
                      "face-to-edge pi-stacking",
                      "face-to-slope pi-stacking",
                      "edge-to-edge pi-stacking",
                      "edge-to-face pi-stacking",
                      "edge-to-slope pi-stacking",
                      "displaced face-to-face pi-stacking",
                      "displaced face-to-edge pi-stacking",
                      "displaced face-to-slope pi-stacking"]


class InteractionFilter:
    """ Filter interactions based on their components.

    Parameters
    ----------
    ignore_self_inter : bool
        If True, ignore interactions involving atoms of the same compound.
    ignore_intra_chain : bool
        If True, ignore intra-chain interactions (e.g., interactions
        between residues in the same protein chain).
    ignore_inter_chain : bool
        If True, ignore interactions between different chains.
    ignore_res_res : bool
        If True, ignore residue-residue interactions.
    ignore_res_nucl : bool
        If True, ignore residue-nucleotide interactions.
    ignore_res_hetatm : bool
        If True, ignore residue-ligand interactions.
    ignore_nucl_nucl : bool
        If True, ignore nucleotide-nucleotide interactions.
    ignore_nucl_hetatm : bool
        If True, ignore nucleotide-ligand interactions.
    ignore_hetatm_hetatm : bool
        If True, ignore ligand-ligand interactions.
    ignore_h2o_h2o : bool
        If True, ignore water-water interactions.
    ignore_any_h2o : bool
        If True, ignore all interactions involving water.
    ignore_multi_comps : bool
        If True, ignore interactions established by atom groups composed of
        multiple compounds (e.g.: amides formed by peptide bonds involve
        two residues).
    ignore_mixed_class : bool
        If True, ignore interactions established by atom groups comprising
        mixed compound classes (e.g. residues and ligands bound by a
        covalent bond).
    """

    def __init__(self, ignore_self_inter=True, ignore_intra_chain=True,
                 ignore_inter_chain=True, ignore_res_res=True,
                 ignore_res_nucl=True, ignore_res_hetatm=True,
                 ignore_nucl_nucl=True, ignore_nucl_hetatm=True,
                 ignore_hetatm_hetatm=True, ignore_h2o_h2o=True,
                 ignore_any_h2o=False, ignore_multi_comps=False,
                 ignore_mixed_class=False):

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
    def new_pli_filter(cls, ignore_res_hetatm=False,
                       ignore_hetatm_hetatm=False, ignore_any_h2o=False,
                       ignore_self_inter=False, **kwargs):
        """Initialize the default filter for protein-ligand interactions.

        Returns
        -------
         : `InteractionFilter`
        """
        return cls(ignore_res_hetatm=ignore_res_hetatm,
                   ignore_hetatm_hetatm=ignore_hetatm_hetatm,
                   ignore_any_h2o=ignore_any_h2o,
                   ignore_self_inter=ignore_self_inter, **kwargs)

    @classmethod
    def new_ppi_filter(cls, ignore_res_res=False, ignore_inter_chain=False,
                       ignore_intra_chain=False, ignore_any_h2o=False,
                       ignore_self_inter=False, **kwargs):
        """Initialize the default filter for protein-protein interactions.

        Returns
        -------
         : `InteractionFilter`
        """
        return cls(ignore_res_res=ignore_res_res,
                   ignore_inter_chain=ignore_inter_chain,
                   ignore_intra_chain=ignore_intra_chain,
                   ignore_any_h2o=ignore_any_h2o,
                   ignore_self_inter=ignore_self_inter, **kwargs)

    @classmethod
    def new_pni_filter(cls, ignore_res_nucl=False, ignore_inter_chain=False,
                       ignore_intra_chain=False, ignore_any_h2o=False,
                       ignore_self_inter=False, **kwargs):
        """Initialize the default filter for protein-nucleotide interactions.

        Returns
        -------
         : `InteractionFilter`
        """
        return cls(ignore_res_nucl=ignore_res_nucl,
                   ignore_inter_chain=ignore_inter_chain,
                   ignore_intra_chain=ignore_intra_chain,
                   ignore_any_h2o=ignore_any_h2o,
                   ignore_self_inter=ignore_self_inter, **kwargs)

    @classmethod
    def new_nni_filter(cls, ignore_nucl_nucl=False, ignore_inter_chain=False,
                       ignore_intra_chain=False, ignore_any_h2o=False,
                       ignore_self_inter=False, **kwargs):
        """Initialize the default filter for nucleotide-nucleotide
        interactions.

        Returns
        -------
         : `InteractionFilter`
        """
        return cls(ignore_nucl_nucl=ignore_nucl_nucl,
                   ignore_inter_chain=ignore_inter_chain,
                   ignore_intra_chain=ignore_intra_chain,
                   ignore_any_h2o=ignore_any_h2o,
                   ignore_self_inter=ignore_self_inter, **kwargs)

    @classmethod
    def new_nli_filter(cls, ignore_nucl_hetatm=False,
                       ignore_hetatm_hetatm=False, ignore_any_h2o=False,
                       ignore_self_inter=False, **kwargs):
        """Initialize the default filter for nucleotide-ligand interactions.

        Returns
        -------
         : `InteractionFilter`
        """
        return cls(ignore_nucl_hetatm=ignore_nucl_hetatm,
                   ignore_hetatm_hetatm=ignore_hetatm_hetatm,
                   ignore_any_h2o=ignore_any_h2o,
                   ignore_self_inter=ignore_self_inter, **kwargs)

    @classmethod
    def from_config_file(cls, config_file):
        """Initialize from a configuration file.

        Parameters
        ----------
        config_file : str
            The configuration file pathname.

        Returns
        -------
         : `InteractionFilter`
        """

        if not path.exists(config_file):
            raise OSError("File '%s' does not exist." % config_file)

        params = Config(config_file)

        inter_filter = None
        if "default" in params.sections():
            params_dict = params.get_section_map("default")
            filter_type = params_dict.pop("type", None)

            if filter_type is not None:
                filter_type = filter_type.upper()

                available_filters = ["PLI", "PPI", "PNI", "NNI", "NLI"]
                if filter_type not in available_filters:
                    raise KeyError("The accepted default filter types are: "
                                   "%s." % ", ".join(available_filters))

                if filter_type == "PLI":
                    inter_filter = cls.new_pli_filter()

                elif filter_type == "PPI":
                    inter_filter = cls.new_ppi_filter()

                elif filter_type == "PNI":
                    inter_filter = cls.new_pni_filter()

                elif filter_type == "NNI":
                    inter_filter = cls.new_nni_filter()

                elif filter_type == "NLI":
                    inter_filter = cls.new_nli_filter()

        if "ignore" in params.sections():
            params_dict = params.get_section_map("ignore")
            params_dict = {"ignore_" + k: v
                           for k, v in params_dict.items()}

            if inter_filter:
                for prop, val in params_dict.items():
                    setattr(inter_filter, prop, val)
            else:
                inter_filter = cls(**params_dict)

        return inter_filter

    def is_valid_pair(self, src_grp, trgt_grp):
        """Evaluate if a pair of atom groups are valid according to the flags
        defined in this class.

        src_grp, trgt_grp : :class:`luna.mol.groups.AtomGroup`
        """

        # It will always ignore interactions involving the same atom groups.
        # Loops in the graph is not permitted and does not make any sense.
        if src_grp == trgt_grp:
            return False

        # It will always ignore interactions involving atoms and the group
        # to which they belong to. For example, the centroid of an aromatic
        # ring cannot interact with an atom that belongs to the ring.
        # It is a type of Loop.
        if src_grp.contain_group(trgt_grp) or trgt_grp.contain_group(src_grp):
            return False

        # If one of the groups contain atoms from different compounds.
        has_multi_comps = (len(src_grp.compounds) > 1
                           or len(trgt_grp.compounds) > 1)
        if self.ignore_multi_comps and has_multi_comps:
            return False

        # If one of the groups contain compounds from different classes as,
        # for instance, residue and ligand. It means that compounds from
        # different classes are covalently bonded to each other.
        has_any_mixed = (src_grp.is_mixed() or trgt_grp.is_mixed())
        if self.ignore_mixed_class and has_any_mixed:
            return False

        # It ignores interactions involving the same compounds if required.
        # As each group may have atoms from different compounds, we can check
        # if there is at least one common compound between the two groups.
        # Remember that if two or more compounds exist in a group, it means
        # that these compounds are covalently bonded and should be considered
        # the same compound.
        # For example: a carbohydrate can be expressed in a PDB as its
        # subparts:
        #      E.g.: GLC + GLC = CBI
        # The same applies to any group formed after covalently bonding a
        # residue to a hetatm (ligand or non-standard amino acid
        # represented as hetatm)
        is_same_compounds = \
            len(src_grp.compounds.intersection(trgt_grp.compounds)) >= 1
        if self.ignore_self_inter and is_same_compounds:
            return False

        # Check if two groups contain the same chains and if both of them
        # contain only one chain. The second condition removes groups
        # containing residues of different chains as may occur due to
        # disulfide bonds. Note, however, that this flag will be used only
        # as a filter for intra-interactions in protein/RNA/DNA chains.
        is_same_chain = (src_grp.get_chains() == trgt_grp.get_chains()
                         and len(src_grp.get_chains()) == 1)

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
        is_res_nucl = ((src_grp.is_residue() and trgt_grp.is_nucleotide())
                       or (src_grp.is_nucleotide() and trgt_grp.is_residue()))
        if self.ignore_res_nucl and is_res_nucl:
            return False

        # It ignores residue-ligand interactions if required.
        is_res_hetatm = ((src_grp.is_residue() and trgt_grp.is_hetatm())
                         or (src_grp.is_hetatm() and trgt_grp.is_residue()))
        if self.ignore_res_hetatm and is_res_hetatm:
            return False

        # Filters for nucleic acid-nucleic acid interactions if required.
        is_nucl_nucl = (src_grp.is_nucleotide() and trgt_grp.is_nucleotide())
        if is_nucl_nucl:
            # Ignore all nucleic acid-nucleic acid interactions
            if self.ignore_nucl_nucl:
                return False
            # Ignore all intra-chain interactions involving
            # two nucleic acids (RNA/DNA chains).
            elif self.ignore_intra_chain and is_same_chain:
                return False
            # Ignore all inter-chain interactions involving
            # two nucleic acids (RNA/DNA chains).
            elif self.ignore_inter_chain and not is_same_chain:
                return False

        # It ignores nucleic acid-ligand interactions if required.
        is_nucl_hetatm = ((src_grp.is_nucleotide() and trgt_grp.is_hetatm())
                          or (src_grp.is_hetatm()
                              and trgt_grp.is_nucleotide()))
        if self.ignore_nucl_hetatm and is_nucl_hetatm:
            return False

        # It ignores ligand-ligand interactions if required.
        is_hetatm_hetatm = (src_grp.is_hetatm() and trgt_grp.is_hetatm())
        if self.ignore_hetatm_hetatm and is_hetatm_hetatm:
            return False

        # It ignores interactions of other compound types with water.
        # It enables the possibility of identifying water-bridged interaction.
        # Eg: residue-water, ligand-water = residue -- water -- ligand.
        is_any_h2o = (src_grp.is_water() or trgt_grp.is_water())
        if self.ignore_any_h2o and is_any_h2o:
            return False

        # It ignores interactions involving two waters if required.
        # If ON, it will produce water-bridged interactions of multiple levels
        # E.g.: residue -- h2o -- h2o -- ligand
        #       residue -- residue -- h2o -- h2o -- ligand.
        is_h2o_h2o = (src_grp.is_water() and trgt_grp.is_water())
        if self.ignore_h2o_h2o and is_h2o_h2o:
            return False

        return True

    def save_config_file(self, config_file):
        """Save the interaction filter parameters into a configuration file.

        Parameters
        ----------
        config_file : str
            The output configuration file.
        """
        with open(config_file, "w") as OUT:
            OUT.write("[ignore]\n")
            OUT.write("self_inter = %s\n" % self.ignore_self_inter)
            OUT.write("intra_chain = %s\n" % self.ignore_intra_chain)
            OUT.write("inter_chain = %s\n" % self.ignore_inter_chain)
            OUT.write("res_res = %s\n" % self.ignore_res_res)
            OUT.write("res_nucl = %s\n" % self.ignore_res_nucl)
            OUT.write("res_hetatm = %s\n" % self.ignore_res_hetatm)
            OUT.write("nucl_nucl = %s\n" % self.ignore_nucl_nucl)
            OUT.write("nucl_hetatm = %s\n" % self.ignore_nucl_hetatm)
            OUT.write("hetatm_hetatm = %s\n" % self.ignore_hetatm_hetatm)
            OUT.write("h2o_h2o = %s\n" % self.ignore_h2o_h2o)
            OUT.write("any_h2o = %s\n" % self.ignore_any_h2o)
            OUT.write("multi_comps = %s\n" % self.ignore_multi_comps)
            OUT.write("mixed_class = %s\n" % self.ignore_mixed_class)


class BindingModeCondition:

    """Define binding mode conditions to filter interactions.

    Parameters
    ----------

    condition : str
        A string defining which chains, compounds, or atoms should be accepted.
        If ``condition`` is the wildcard '*', then all chains, compounds, and
        atoms will be considered valid. Otherwise, ``condition`` should have
        the format '<CHAIN ID>/<COMPOUND NAME>/<COMPOUND NUMBER>/<ATOM>'.
        Wildcards are accepted for each one of these fields.
        For example:

            * '\\*/HIS/\\*/\\*': represents all histidines' \
                                 atoms from all chains.
            * 'A/CBL/\\*/\\*' represents all ligands named CBL \
                            from chain A.
            * 'B/HIS/\\*/N\\*' represents all histidines' nitrogens \
                             from chain B.

    Attributes
    ----------
    accept_all : bool
        If True, all chains, compounds, and atoms will be considered valid.
    accept_all_chains : bool
        If True, all chains will be considered valid.
    accept_all_comps : bool
        If True, all compound names will be considered valid.
    accept_all_comp_nums : bool
        If True, all compound numbers (residue sequence number in the PDB
        format) will be considered valid.
    accept_all_atoms : bool
        If True, all atoms will be considered valid.
    chain_id : str or None
        If provided, accept only chains whose id matches ``chain_id``.
    comp_name : str or None
        If provided, accept only compounds whose name matches ``comp_name``.
    comp_num : int or None
        If provided, accept only compounds whose sequence number matches
        ``comp_num``.
    comp_icode : str or None
        If provided, accept only compounds whose insertion code matches
        ``comp_icode``.
    atom : str or None
        If provided, accept only atoms whose name matches ``atom``.
    """

    def __init__(self, condition):
        self.accept_all = False
        self.accept_all_chains = False
        self.accept_all_comps = False
        self.accept_all_comp_nums = False
        self.accept_all_atoms = False

        self.chain_id = None
        self.comp_name = None
        self.comp_num = None
        self.comp_icode = None
        self.atom = None

        self._condition_repr = condition

        if condition is not None:
            self._parse_condition(condition.upper())

    def _parse_condition(self, condition):
        # Accept everything.
        if condition == "*":
            self.accept_all = True
        else:
            chain, comp_name, comp_num, atom = condition.split("/")

            if chain == "*":
                self.accept_all_chains = True
            else:
                self.chain_id = chain

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
                        error_msg = ("The informed compound number '%s' is "
                                     "invalid. It must be an integer."
                                     % str(comp_num))
                        raise IllegalArgumentError(error_msg)

                    icode = (None if matched.group(2) == ""
                             else matched.group(2))
                else:
                    error_msg = ("The compound number and its insertion code "
                                 "(if applicable) '%s' is invalid. It must be "
                                 "an integer followed by one insertion code "
                                 "character when applicable." % comp_num)
                    raise IllegalArgumentError(error_msg)

                self.comp_num = comp_num
                self.comp_icode = icode

            if atom == "*":
                self.accept_all_atoms = True
            else:
                self.atom = atom

    def is_valid(self, atm_grp):
        """Check if an atom group is valid or not based on this condition.

        atm_grp : :class:`luna.mol.groups.AtomGroup`
        """

        # Accept everything.
        if self.accept_all:
            return True

        # Accept everything.
        if (self.accept_all_chains and self.accept_all_comps
                and self.accept_all_comp_nums and self.accept_all_atoms):
            return True

        # Tries to identify the first valid compound in the atom group.
        for atm in atm_grp.atoms:
            comp = atm.parent

            is_chain_valid = self.accept_all_chains
            is_comp_valid = self.accept_all_comps
            is_comp_num_valid = self.accept_all_comp_nums
            is_atom_valid = self.accept_all_atoms

            if self.chain_id is not None and self.chain_id == comp.parent.id:
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

            if (is_chain_valid and is_comp_valid
                    and is_comp_num_valid and is_atom_valid):
                return True

        return False

    def __repr__(self):
        return "<BindingModeCondition: %s>" % self._condition_repr


class BindingModeFilter:

    """Filter interactions based on a set of binding mode conditions.

    Parameters
    ----------
    config : dict of {str : iterable}
        A dict defining binding modes and how interactions should be validated.
        Each key represents an interaction type and values are an iterable of
        `BindingModeCondition` instances.
    """

    def __init__(self, config, names_map=None):
        self.config = config
        self.names_map = names_map or {}

    @classmethod
    def from_config_file(cls, config_file):
        """Initialize from a configuration file.

        Parameters
        ----------
        config_file : str
            The configuration file pathname.

        Returns
        -------
         : `BindingModeFilter`

        Examples
        --------

        It follows an example of a configuration file::

            ; To configurate an interaction type, create a new line and define
            ; the interaction: [New interaction]. Then you can define whether
            ; or not all interactions must be accepted by setting 'accept_only'
            ; to True or False.

            ; If you want to specify binding modes, use the variable
            ; 'accept_only', which expects a list of strings in the format:
            ; <CHAIN ID>/<COMPOUND NAME>/<COMPOUND NUMBER>/<ATOM>.
            ;
            ; Wildcards are accepted for the expected fields.
            ; For example, "*/HIS/*/*" represents all histidines' atoms from \
all chains.
            ;               "A/CBL/*/*" represents all ligands named CBL from \
chain A.
            ;               "B/HIS/*/N*" represents all histidines' nitrogens \
from chain B.

            [Hydrogen bond]
            accept_only = ["A/LYS/245/*", "*/HIS/*/*"]

            [Hydrophobic]
            accept_all = True

            [Cation-pi]
            accept_only = ["*"]
            accept_all = True

            [Weak hydrogen bond]
            accept_all = False
            accept_only = ["*/THR/434/O*"]

            [Face-to-edge pi-stacking]
            accept_only = ["*"]

            [Aromatic stacking]
            accept_all = True

            [*]
            accept_all = False
        """

        if not path.exists(config_file):
            raise OSError("File '%s' does not exist." % config_file)

        filtering_config = {}
        names_mapping = {}

        config = Config(config_file)
        for inter_type in config.sections():
            params = config.get_section_map(inter_type)

            ori_name = inter_type
            inter_type = inter_type.lower()

            accept_all = False
            values = []
            if "accept_all" in params:
                accept_all = literal_eval(params["accept_all"])
                if accept_all is True:
                    values = ["*"]
                else:
                    values = [None]

            elif "accept_only" in params and not accept_all:
                values = literal_eval(params["accept_only"])

                if "*" in values:
                    values = ["*"]

            conditions = [BindingModeCondition(condition)
                          for condition in values]

            filtering_config[inter_type] = conditions
            names_mapping[inter_type] = ori_name

        return cls(filtering_config, names_mapping)

    def is_valid(self, inter):
        """Check if an interaction is valid or not based on this binding mode
        configuration.

        inter : :class:`luna.interaction.type.InteractionType`
        """
        if inter.type.lower() in self.config:
            inter_type = inter.type.lower()

        elif ("aromatic stacking" in self.config
                and inter.type.lower() in AROMATIC_STACKINGS):
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

    def save_config_file(self, config_file):
        """Save the binding mode filter parameters into a
        configuration file.

        Parameters
        ----------
        config_file : str
            The output configuration file.
        """
        args_by_section = defaultdict(list)

        for inter, conditions in self.config.items():
            section_name = self.names_map.get(inter, inter)

            multi_conditions = []
            for condition in conditions:
                if (condition.accept_all is False
                        and condition._condition_repr is None):
                    args_by_section[section_name].append(("accept_all", False))

                elif condition.accept_all is True:
                    args_by_section[section_name].append(("accept_all", True))

                else:
                    multi_conditions.append(condition._condition_repr)

            if multi_conditions:
                args_by_section[section_name].append(("accept_only",
                                                      str(multi_conditions)))

        with open(config_file, "w") as OUT:
            for k in args_by_section:
                OUT.write("[%s]\n" % k)
                for arg, val in args_by_section[k]:
                    OUT.write("%s = %s\n" % (arg, val))
                OUT.write("\n")
