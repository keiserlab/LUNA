import ast
from os import path

from luna.util.config import Config
from luna.mol.entry import Entry
from luna.interaction.calc import InteractionCalculator, DEFAULT_SOLVENTS
from luna.interaction.config import (InteractionConfig,
                                     DefaultInteractionConfig)
from luna.interaction.filter import InteractionFilter, BindingModeFilter
from luna.interaction.fp.type import IFPType
from luna.util.default_values import ATOM_PROP_FILE, ACCEPTED_MOL_OBJ_TYPES


class ProjectParams:

    def __init__(self, config_file=None):

        if config_file is None:
            config_file = path.abspath(path.join(path.realpath(__file__),
                                                 "../default.cfg"))

        self.config_file = config_file

        self._parse()

    def _parse_value(self, params, param_name, fallback=None):

        try:
            value = params.get(param_name)
        except ValueError:
            return fallback

        try:
            return ast.literal_eval(value)
        except (ValueError, SyntaxError):
            return value

    def _validate_type(self, value, expected_type,
                       section, param_name, is_none_valid=False):
        if is_none_valid and value is None:
            return value

        if type(value) != expected_type:
            # Ok, if it's an int and the expected value is float.
            if expected_type == float and type(value) == int:
                return value
            # Ok, if it's an integer represented as float.
            if (expected_type == int and type(value) == float
                    and value.is_integer()):
                return int(value)

            raise TypeError("The property '%s' from section '%s' "
                            "must be a(n) %s."
                            % (param_name, section, expected_type.__name__))
        return value

    def _parse_entries_params(self, params):
        section = 'entries'

        entries_file = self._validate_type(params.pop("input", None),
                                           str, section, 'input', True)

        pdb_id = self._validate_type(params.pop("pdb_id", None),
                                     str, section, 'pdb_id', True)

        mol_file = self._validate_type(params.pop("mol_file", None),
                                       str, section, 'mol_file', True)

        entries_sep = self._validate_type(params.pop("entries_sep", ":"),
                                          str, section, 'entries_sep', True)

        fields_sep = self._validate_type(params.pop("fields_sep", ","),
                                         str, section, 'fields_sep', True)

        entries = None
        if entries_file:
            entries = list(Entry.from_file(entries_file, pdb_id, mol_file,
                                           entries_sep, fields_sep))

        return {"entries": entries}

    def _parse_paths_params(self, params):
        section = 'paths'

        working_path = self._validate_type(params.pop("working_path", None),
                                           str, section, 'working_path', True)

        pdb_path = self._validate_type(params.pop("pdb_path", None),
                                       str, section, 'pdb_path', True)

        overwrite_path = self._validate_type(params.pop("overwrite_path",
                                                        False),
                                             bool, section, 'overwrite_path')

        params = {"pdb_path": pdb_path, "overwrite_path": overwrite_path}

        if working_path is not None:
            params["working_path"] = working_path

        return params

    def _parse_hydrogens_params(self, params):
        section = 'hydrogens'

        add_h = self._validate_type(params.pop("add", True),
                                    bool, section, 'add')

        ph = self._validate_type(params.pop("ph", 7.4),
                                 float, section, 'ph')

        return {"add_h": add_h, "ph": ph}

    def _parse_standardization_params(self, params):
        section = 'standardization'
        amend_mol = self._validate_type(params.pop("amend_mol", True),
                                        bool, section, 'amend_mol')
        return {"amend_mol": amend_mol}

    def _parse_features_params(self, params):
        section = 'features'

        config_file = self._validate_type(params.pop("config_file", None),
                                          str, section, 'config_file', True)
        if config_file is None:
            config_file = ATOM_PROP_FILE

        return {"atom_prop_file": config_file}

    def _new_inter_calc(self, params):
        section = 'interaction'

        config_file = self._validate_type(params.pop("config_file", None),
                                          str, section, 'config_file', True)
        if config_file is None:
            inter_config = DefaultInteractionConfig()
        else:
            inter_config = InteractionConfig.from_config_file(config_file)

        filter_file = self._validate_type(params.pop("filter_file", None),
                                          str, section, 'filter_file', True)
        inter_filter = None
        if filter_file is not None:
            inter_filter = InteractionFilter.from_config_file(filter_file)

        add_non_cov = self._validate_type(params.pop("add_non_cov", True),
                                          bool, section, 'add_non_cov')

        add_cov = self._validate_type(params.pop("add_cov", True),
                                      bool, section, 'add_cov')

        add_proximal = self._validate_type(params.pop("add_proximal", False),
                                           bool, section, 'add_proximal')

        add_atom_atom = self._validate_type(params.pop("add_atom_atom", True),
                                            bool, section, 'add_atom_atom')

        add_dep_inter = self._validate_type(params.pop("add_dependent_inter",
                                                       False),
                                            bool, section,
                                            'add_dependent_inter')

        prop_val = params.pop("add_h2o_pairs_with_no_target", False)
        add_h2o_pairs = self._validate_type(prop_val, bool, section,
                                            'add_h2o_pairs_with_no_target')

        sdonor_rules = self._validate_type(params.pop("strict_donor_rules",
                                                      True),
                                           bool, section, 'strict_donor_rules')

        prop_name = "strict_weak_donor_rules"
        swdonor_rules = self._validate_type(params.pop(prop_name, True),
                                            bool, section, prop_name)

        prop_name = "lazy_comps_list"
        lazy_comps_list = self._validate_type(params.pop(prop_name, None),
                                              list, section, prop_name, True)
        if lazy_comps_list is None:
            lazy_comps_list = DEFAULT_SOLVENTS

        params = {
            "inter_config": inter_config,
            "inter_filter": inter_filter,
            "add_non_cov": add_non_cov,
            "add_cov": add_cov,
            "add_proximal": add_proximal,
            "add_atom_atom": add_atom_atom,
            "add_dependent_inter": add_dep_inter,
            "add_h2o_pairs_with_no_target": add_h2o_pairs,
            "strict_donor_rules": sdonor_rules,
            "strict_weak_donor_rules": swdonor_rules,
            "lazy_comps_list": lazy_comps_list,
        }

        return {"inter_calc": InteractionCalculator(**params)}

    def _new_binding_mode_filter(self, params):
        section = 'binding_mode_filter'

        config_file = self._validate_type(params.pop("config_file", None),
                                          str, section, 'config_file', True)

        bmf = None
        if config_file:
            bmf = BindingModeFilter.from_config_file(config_file)
        return {"binding_mode_filter": bmf}

    def _parse_mfp_params(self, params):
        section = 'mfp'

        calc_mfp = self._validate_type(params.pop("calculate", False),
                                       bool, section, 'calculate')

        mfp_output = self._validate_type(params.pop("output", None), str,
                                         section, 'output', True)

        return {"calc_mfp": calc_mfp, "mfp_output": mfp_output}

    def _parse_ifp_params(self, params):
        section = 'ifp'

        calc_ifp = self._validate_type(params.pop("calculate", False),
                                       bool, section, 'calculate')

        num_levels = self._validate_type(params.pop("num_levels", 2),
                                         int, section, 'num_levels')

        radius_step = self._validate_type(params.pop("radius_step", 5.73171),
                                          float, section, 'radius_step')

        ifp_length = self._validate_type(params.pop("length", 4096), int,
                                         section, 'length')

        ifp_count = self._validate_type(params.pop("count", True), bool,
                                        section, 'count')

        diff_comp_classes = params.pop("diff_comp_classes", True)
        diff_comp_classes = self._validate_type(diff_comp_classes, bool,
                                                section, 'diff_comp_classes')

        ifp_type = params.pop("type", "EIFP")
        try:
            ifp_type = IFPType[ifp_type]
        except KeyError:
            valid_ifps = ", ".join([e.name for e in IFPType])
            raise KeyError("The accepted IFP types are: "
                           "%s." % valid_ifps)

        ifp_output = self._validate_type(params.pop("output", None), str,
                                         section, 'output', True)

        sim_matrix_output = params.pop("sim_matrix_output", IFPType.EIFP)
        sim_matrix_output = self._validate_type(sim_matrix_output, str,
                                                section, 'sim_matrix_output',
                                                True)

        return {
            "calc_ifp": calc_ifp,
            "ifp_num_levels": num_levels,
            "ifp_radius_step": radius_step,
            "ifp_length": ifp_length,
            "ifp_count": ifp_count,
            "ifp_diff_comp_classes": diff_comp_classes,
            "ifp_type": ifp_type,
            "ifp_output": ifp_output,
            "ifp_sim_matrix_output": sim_matrix_output,
        }

    def _parse_pse_params(self, params):
        section = 'pse'
        out_pse = self._validate_type(params.pop("generate", False), bool,
                                      section, 'generate')
        return {"out_pse": out_pse}

    def _parse_general_params(self, params):
        section = 'general'

        mol_obj_type = self._validate_type(params.pop("mol_obj_type", 'rdkit'),
                                           str, section, 'mol_obj_type')
        if mol_obj_type not in ACCEPTED_MOL_OBJ_TYPES:
            raise KeyError("Invalid choise '%s' for property 'mol_obj_type' "
                           "from section '%s'. Choose from: %s."
                           % (mol_obj_type, section,
                              ", ".join(ACCEPTED_MOL_OBJ_TYPES)))

        append_mode = self._validate_type(params.pop("append_mode", False),
                                          bool, section, 'append_mode')

        verbosity = self._validate_type(params.pop("verbosity", 3),
                                        int, section, 'verbosity')

        logging_enabled = self._validate_type(params.pop("logging_enabled",
                                                         True),
                                              bool, section, 'logging_enabled')

        nproc = self._validate_type(params.pop("nproc", -1),
                                    int, section, 'nproc', True)

        return {
            "mol_obj_type": mol_obj_type,
            "append_mode": append_mode,
            "verbosity": verbosity,
            "logging_enabled": logging_enabled,
            "nproc": nproc,
        }

    def _parse(self):
        config = Config(self.config_file)

        proj_params = {}
        for section in config.sections():
            params_dict = config.get_section_map(section)

            params_dict = {param_name: self._parse_value(params_dict,
                                                         param_name)
                           for param_name in params_dict}

            if section == "entries":
                params_dict = self._parse_entries_params(params_dict)

            if section == "paths":
                params_dict = self._parse_paths_params(params_dict)

            if section == "hydrogens":
                params_dict = self._parse_hydrogens_params(params_dict)

            if section == "standardization":
                params_dict = self._parse_standardization_params(params_dict)

            if section == "features":
                params_dict = self._parse_features_params(params_dict)

            if section == "interactions":
                params_dict = self._new_inter_calc(params_dict)

            if section == "binding_mode_filter":
                params_dict = self._new_binding_mode_filter(params_dict)

            if section == "mfp":
                params_dict = self._parse_mfp_params(params_dict)

            if section == "ifp":
                params_dict = self._parse_ifp_params(params_dict)

            if section == "pse":
                params_dict = self._parse_pse_params(params_dict)

            if section == "general":
                params_dict = self._parse_general_params(params_dict)

            proj_params.update(params_dict)

        self.params = proj_params
