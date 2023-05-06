import ast
from os import path
from pathlib import Path
import logging

from luna.util.config import Config
from luna.mol.entry import *
from luna.interaction.calc import InteractionCalculator
from luna.interaction.config import InteractionConfig
from luna.interaction.filter import InteractionFilter, BindingModeFilter
from luna.interaction.fp.type import IFPType
from luna.util.logging import VERBOSITY_LEVEL
from luna.util.exceptions import IllegalArgumentError


class ProjectParams(dict):

    """Define LUNA project parameters.

    Parameters
    ----------
    params : dict
        Set LUNA parameters from ``params``.

    fill_defaults : bool
        If True, initialize `ProjectParams` with default LUNA values.

        .. note::
            Any parameter provided in ``params`` will overwrite the
            default values.
    """

    def __init__(self, params=None, fill_defaults=False):

        params_to_parse = {}

        # Loads in default LUNA values.
        if fill_defaults:
            luna_path = path.abspath(path.join(path.realpath(__file__),
                                               '../..'))

            default_config_file = f"{luna_path}/config/default.cfg"
            config = Config(default_config_file)

            params_to_parse = config.as_dict(only_props=True)

            if params_to_parse["feat_cfg"] == "None":
                atom_prop_file = f"{luna_path}/data/LUNA.fdef"
                params_to_parse["feat_cfg"] = atom_prop_file

            if params_to_parse["inter_cfg"] == "None":
                inter_config_file = f"{luna_path}/interaction/config.cfg"
                params_to_parse["inter_cfg"] = inter_config_file

        # Overwrites default values with provided params.
        if params:
            params_to_parse.update(params)

        self._validate(params_to_parse)

    @classmethod
    def from_project_obj(cls, proj_obj):
        """Initialize a `ProjectParams` from a `~luna.projects.Project` object.

        Returns
        -------
         : `ProjectParams`
        """
        params = {"add_h": proj_obj.add_h,
                  "amend_mol": proj_obj.amend_mol,
                  "append_mode": proj_obj.append_mode,
                  "feat_cfg": proj_obj.atom_prop_file,
                  "binding_mode_filter": proj_obj.binding_mode_filter,
                  "ifp": proj_obj.calc_ifp,
                  "mfp": proj_obj.calc_mfp,
                  "entries": proj_obj.entries,
                  "ifp_count": proj_obj.ifp_count,
                  "ifp_diff_comp_classes": proj_obj.ifp_diff_comp_classes,
                  "ifp_length": proj_obj.ifp_length,
                  "ifp_num_levels": proj_obj.ifp_num_levels,
                  "ifp_output": proj_obj.ifp_output,
                  "ifp_radius_step": proj_obj.ifp_radius_step,
                  "ifp_sim_matrix_output": proj_obj.ifp_sim_matrix_output,
                  "ifp_type": proj_obj.ifp_type,
                  "inter_calc": proj_obj.inter_calc,
                  "logging_enabled": proj_obj.logging_enabled,
                  "mfp_output": proj_obj.mfp_output,
                  "nproc": proj_obj.nproc,
                  "pse": proj_obj.out_pse,
                  "overwrite_path": proj_obj.overwrite_path,
                  "pdb_path": proj_obj.pdb_path,
                  "ph": proj_obj.ph,
                  "pse_path": proj_obj.pse_path,
                  "use_cache": proj_obj.use_cache,
                  "verbosity": proj_obj.verbosity,
                  "working_path": proj_obj.working_path}

        return cls(params)

    def save_config_file(self, config_file):
        """Save the project parameters into a configuration file.

        Parameters
        ----------
        config_file : str
            The output configuration file.
        """
        with open(config_file, "w") as OUT:
            OUT.write("[entries]\n")
            for k, v in self._get_entry_params(config_file).items():
                OUT.write(f"{k} = {v}\n")
            OUT.write("\n")

            OUT.write("[paths]\n")
            OUT.write("working_path = %s\n" % self["working_path"])
            OUT.write("pdb_path = %s\n" % self["pdb_path"])
            OUT.write("overwrite_path = %s\n" % self["overwrite_path"])
            OUT.write("\n")

            OUT.write("[hydrogens]\n")
            OUT.write("add_h = %s\n" % self["add_h"])
            OUT.write("ph = %s\n" % self["ph"])
            OUT.write("\n")

            OUT.write("[standardization]\n")
            OUT.write("amend_mol = %s\n" % self["amend_mol"])
            OUT.write("\n")

            OUT.write("[features]\n")
            OUT.write("feat_cfg = %s\n" % self["atom_prop_file"])
            OUT.write("\n")

            OUT.write("[interactions]\n")
            for k, v in self._get_inter_params(config_file).items():
                OUT.write(f"{k} = {v}\n")
            OUT.write("\n")

            OUT.write("[binding_mode_filter]\n")
            for k, v in self._get_binding_filter_params(config_file).items():
                OUT.write(f"{k} = {v}\n")
            OUT.write("\n")

            OUT.write("[mfp]\n")
            OUT.write("mfp = %s\n" % self["calc_mfp"])
            OUT.write("mfp_output = %s\n" % self["mfp_output"])
            OUT.write("\n")

            OUT.write("[ifp]\n")
            OUT.write("ifp = %s\n" % self["calc_ifp"])
            OUT.write("ifp_num_levels = %s\n" % self["ifp_num_levels"])
            OUT.write("ifp_radius_step = %s\n" % self["ifp_radius_step"])
            OUT.write("ifp_length = %s\n" % self["ifp_length"])
            OUT.write("ifp_count = %s\n" % self["ifp_count"])
            OUT.write("ifp_type = %s\n" % self["ifp_type"].name)
            OUT.write("ifp_diff_comp_classes = %s\n"
                      % self["ifp_diff_comp_classes"])
            OUT.write("ifp_output = %s\n" % self["ifp_output"])
            OUT.write("ifp_sim_matrix_output = %s\n"
                      % self["ifp_sim_matrix_output"])
            OUT.write("\n")

            OUT.write("[pse]\n")
            OUT.write("pse = %s\n" % self["out_pse"])
            OUT.write("pse_path = %s\n" % self["pse_path"])
            OUT.write("\n")

            verbosity = [k for k, v in VERBOSITY_LEVEL.items()
                         if v == self["verbosity"]].pop()

            OUT.write("[general]\n")
            OUT.write("append_mode = %s\n" % self["append_mode"])
            OUT.write("use_cache = %s\n" % self["use_cache"])
            OUT.write("verbosity = %s\n" % verbosity)
            OUT.write("logging_enabled = %s\n" % self["logging_enabled"])
            OUT.write("nproc = %s\n" % self["nproc"])
            OUT.write("\n")

    def _get_entry_params(self, config_file):
        sep = ":"
        entries = []
        for e in self["entries"]:
            if isinstance(e, MolEntry):
                entries.append(e.to_string(sep))

            elif isinstance(e, MolFileEntry):
                fields = [str(e.pdb_id), str(e.mol_id), str(e.mol_file),
                          str(e.is_multimol_file)]
                entries.append(",".join(fields))

        path = Path(config_file).parent.resolve()
        filename = Path(config_file).stem

        entry_file = f"{path}/{filename}.entries.txt"
        with open(entry_file, "w") as OUT:
            for e in entries:
                OUT.write(f"{e}\n")

        return {"entries": entry_file,
                "pdb_id": None,
                "mol_file": None,
                "entries_sep": ":",
                "fields_sep": ","}

    def _get_inter_params(self, config_file):
        path = Path(config_file).parent.resolve()
        filename = Path(config_file).stem

        ic = self["inter_calc"]

        inter_config = ic.inter_config
        inter_cfg_file = None
        if isinstance(inter_config, InteractionConfig):
            inter_cfg_file = f"{path}/{filename}.inter.txt"
            inter_config.save_config_file(inter_cfg_file)

        inter_filter = ic.inter_filter
        filter_cfg_file = None
        if isinstance(inter_filter, InteractionFilter):
            filter_cfg_file = f"{path}/{filename}.filter.txt"
            inter_filter.save_config_file(filter_cfg_file)

        inter_filter = ic.inter_filter

        return {"inter_cfg": inter_cfg_file,
                "filter_cfg": filter_cfg_file,
                "add_non_cov": ic.add_non_cov,
                "add_cov": ic.add_cov,
                "add_proximal": ic.add_proximal,
                "add_atom_atom": ic.add_atom_atom,
                "add_dependent_inter": ic.add_dependent_inter,
                "add_h2o_pairs_with_no_target":
                    ic.add_h2o_pairs_with_no_target,
                "strict_donor_rules": ic.strict_donor_rules,
                "strict_weak_donor_rules": ic.strict_weak_donor_rules,
                "lazy_comps_list": ",".join(ic.lazy_comps_list)}

    def _get_binding_filter_params(self, config_file):
        path = Path(config_file).parent.resolve()
        filename = Path(config_file).stem

        bmf = self["binding_mode_filter"]
        bmf_cfg_file = None
        if isinstance(bmf, BindingModeFilter):
            bmf_cfg_file = f"{path}/{filename}.bind.txt"
            bmf.save_config_file(bmf_cfg_file)

        return {"bind_cfg": bmf_cfg_file}

    def _get_value(self, params, param_name, dtype,
                   fallback=None, is_none_valid=True):

        value = params.get(param_name, fallback)

        try:
            value = ast.literal_eval(value)
        except (ValueError, SyntaxError):
            pass

        if is_none_valid and value is None:
            return value

        if type(value) != dtype:
            # Ok, if it's an int and the expected value is float.
            if dtype == float and type(value) == int:
                return value
            # Ok, if it's an integer represented as float.
            if (dtype == int and type(value) == float
                    and value.is_integer()):
                return int(value)

            raise TypeError("The parameter '%s' must be a(n) %s."
                            % (param_name, dtype.__name__))

        return value

    def _parse_entries_params(self, params):
        try:
            entries_file = self._get_value(params, "entries", str)

            pdb_id = self._get_value(params, "pdb_id", str)
            mol_file = self._get_value(params, "mol_file", str)
            entries_sep = self._get_value(params, "entries_sep", str)
            fields_sep = self._get_value(params, "fields_sep", str)

            mol_obj_type = params["mol_obj_type"]

            entries = None
            if entries_file:
                entries = list(Entry.from_file(entries_file, pdb_id, mol_file,
                                               entries_sep, fields_sep,
                                               mol_obj_type=mol_obj_type))
        except IllegalArgumentError:
            raise
        except Exception:
            entries = self._get_value(params, "entries", list)

        return {"entries": entries}

    def _parse_paths_params(self, params):
        working_path = self._get_value(params, "working_path", str)
        pdb_path = self._get_value(params, "pdb_path", str)
        overwrite_path = self._get_value(params, "overwrite_path", bool)

        return {"working_path": working_path,
                "pdb_path": pdb_path,
                "overwrite_path": overwrite_path}

    def _parse_hydrogens_params(self, params):
        add_h = self._get_value(params, "add_h", bool)
        ph = self._get_value(params, "ph", float)

        return {"add_h": add_h, "ph": ph}

    def _parse_standardization_params(self, params):
        amend_mol = self._get_value(params, "amend_mol", bool)

        return {"amend_mol": amend_mol}

    def _parse_features_params(self, params):
        config_file = self._get_value(params, "feat_cfg", str)

        return {"atom_prop_file": config_file}

    def _parse_filter_params(self, params):
        filter_types = ["pli", "ppi", "pni", "nni", "nli"]

        if "default_filter" in params:
            default_filter = params["default_filter"]

            if default_filter not in filter_types:
                raise IllegalArgumentError("The informed default filter '%s' "
                                           "is not valid. The valid default "
                                           "filters are: %s."
                                           % (default_filter,
                                              ", ".join(filter_types)))

            func_name = f"new_{default_filter}_filter"
            func = getattr(InteractionFilter, func_name)
            return func()

        else:
            filter_file = self._get_value(params, "filter_cfg", str)
            inter_filter = None
            if filter_file is not None:
                inter_filter = InteractionFilter.from_config_file(filter_file)

            return inter_filter

    def _parse_inter_params(self, params):

        try:
            ic = self._get_value(params, "inter_calc",
                                 InteractionCalculator,
                                 is_none_valid=False)

        except Exception:
            config_file = self._get_value(params, "inter_cfg", str)

            inter_config = None
            if config_file is not None:
                inter_config = InteractionConfig.from_config_file(config_file)

            inter_filter = self._parse_filter_params(params)

            add_non_cov = self._get_value(params, "add_non_cov", bool)

            add_cov = self._get_value(params, "add_cov", bool)

            add_proximal = self._get_value(params, "add_proximal", bool)

            add_atom_atom = self._get_value(params, "add_atom_atom", bool)

            add_dep_inter = self._get_value(params, "add_dependent_inter",
                                            bool)

            add_h2o_pairs = self._get_value(params,
                                            "add_h2o_pairs_with_no_target",
                                            bool)

            sdonor_rules = self._get_value(params, "strict_donor_rules", bool)

            swdonor_rules = self._get_value(params, "strict_weak_donor_rules",
                                            bool)

            lazy_comps_list = self._get_value(params, "lazy_comps_list", str)

            if lazy_comps_list is not None:
                lazy_comps_list = lazy_comps_list.split(",")

            ic = None
            if inter_config is not None:
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

                ic = InteractionCalculator(**params)

        return {"inter_calc": ic}

    def _parse_binding_mode_filter(self, params):
        try:
            bmf = self._get_value(params, "binding_mode_filter",
                                  BindingModeFilter,
                                  is_none_valid=False)
        except Exception:
            config_file = self._get_value(params, "bind_cfg", str)

            bmf = None
            if config_file is not None:
                bmf = BindingModeFilter.from_config_file(config_file)

        return {"binding_mode_filter": bmf}

    def _parse_mfp_params(self, params):
        calc_mfp = self._get_value(params, "mfp", bool)
        mfp_output = self._get_value(params, "mfp_output", str)

        return {"calc_mfp": calc_mfp, "mfp_output": mfp_output}

    def _parse_ifp_params(self, params):

        calc_ifp = self._get_value(params, "ifp", bool)

        num_levels = self._get_value(params, "ifp_num_levels", int)

        radius_step = self._get_value(params, "ifp_radius_step", float)

        ifp_length = self._get_value(params, "ifp_length", int)

        ifp_count = self._get_value(params, "ifp_count", bool)

        ifp_diff_comp_classes = self._get_value(params,
                                                "ifp_diff_comp_classes",
                                                bool)

        try:
            ifp_type = self._get_value(params, "ifp_type", str)

            if ifp_type is not None:
                try:
                    ifp_type = IFPType[ifp_type]
                except KeyError:
                    valid_ifps = ", ".join([e.name for e in IFPType])
                    raise KeyError("The accepted IFP types are: "
                                   "%s." % valid_ifps)
        except Exception:
            ifp_type = self._get_value(params, "ifp_type", IFPType)

        ifp_output = self._get_value(params, "ifp_output", str)

        ifp_sim_matrix_output = self._get_value(params,
                                                "ifp_sim_matrix_output", str)

        return {
            "calc_ifp": calc_ifp,
            "ifp_num_levels": num_levels,
            "ifp_radius_step": radius_step,
            "ifp_length": ifp_length,
            "ifp_count": ifp_count,
            "ifp_diff_comp_classes": ifp_diff_comp_classes,
            "ifp_type": ifp_type,
            "ifp_output": ifp_output,
            "ifp_sim_matrix_output": ifp_sim_matrix_output,
        }

    def _parse_pse_params(self, params):
        out_pse = self._get_value(params, "pse", bool)
        pse_path = self._get_value(params, "pse_path", str)

        return {"out_pse": out_pse, "pse_path": pse_path}

    def _parse_general_params(self, params):

        append_mode = self._get_value(params, "append_mode", bool)

        use_cache = self._get_value(params, "use_cache", bool)

        verbosity = self._get_value(params, "verbosity", int)

        if verbosity is not None:
            try:
                verbosity = [verbosity for k, v in VERBOSITY_LEVEL.items()
                             if v == verbosity].pop()
            except Exception:
                if verbosity not in VERBOSITY_LEVEL:
                    lvls = sorted(VERBOSITY_LEVEL.items())
                    msg = ("The informed verbosity level '%s' is not valid. "
                           "The valid levels are: %s."
                           % (repr(verbosity),
                              ", ".join(["%d (%s)"
                                         % (k, logging.getLevelName(v))
                                         for k, v in lvls])))
                    raise IllegalArgumentError(msg)

        logging_enabled = self._get_value(params, "logging_enabled", bool)

        nproc = self._get_value(params, "nproc", int)

        return {
            "append_mode": append_mode,
            "use_cache": use_cache,
            "verbosity": verbosity,
            "logging_enabled": logging_enabled,
            "nproc": nproc,
        }

    def _validate(self, params):
        validated_params = {}

        # Entries
        validated_params.update(self._parse_entries_params(params))

        # Paths
        validated_params.update(self._parse_paths_params(params))

        # Hydrogens
        validated_params.update(self._parse_hydrogens_params(params))

        # Standardization
        validated_params.update(self._parse_standardization_params(params))

        # Features
        validated_params.update(self._parse_features_params(params))

        # Interactions
        validated_params.update(self._parse_inter_params(params))

        # Binding mode filter
        validated_params.update(self._parse_binding_mode_filter(params))

        # MFPs
        validated_params.update(self._parse_mfp_params(params))

        # IFPs
        validated_params.update(self._parse_ifp_params(params))

        # PSE
        validated_params.update(self._parse_pse_params(params))

        # General
        validated_params.update(self._parse_general_params(params))

        super().__init__(validated_params)
