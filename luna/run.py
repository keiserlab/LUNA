import argparse
import os
import shutil
import configparser

from luna.projects import LocalProject, EntryResults
from luna.config.params import ProjectParams
from luna.interaction.view import InteractionViewer
from luna.interaction.fp.type import IFPType
from luna.util.file import create_directory
from luna.interaction.calc import InteractionsManager
from luna.util.logging import VERBOSITY_LEVEL

IFP_TYPES = {
    "EIFP": IFPType.EIFP,
    "FIFP": IFPType.FIFP,
    "HIFP": IFPType.HIFP,
}


class NegateAction(argparse.Action):
    def __call__(self, parser, ns, values, option):
        setattr(ns, self.dest, not option.startswith('--no-'))


class CountIFPAction(argparse.Action):
    def __call__(self, parser, ns, values, option):
        setattr(ns, self.dest, option == "-C")


def get_parser():
    parser = argparse.ArgumentParser()

    group = parser.add_argument_group('Entries')
    group.add_argument('-e', "--entries", type=str,
                       help="an input file containing a list of "
                            "entries to process.")
    group.add_argument("-p", "--prot", dest="pdb_id",
                       type=str,
                       help="A local PDB filename. "
                            "Mandatory if there is at least one entry "
                            "containing only the ligand name. "
                            "This only applies to ligands to be read "
                            "from a multimolecular file.")
    group.add_argument('-l', "--lig", dest="mol_file", type=str,
                       help="a molecular file containing 1 or more ligands. "
                            "Mandatory if there is at least one entry "
                            "containing only the ligand name. "
                            "This only applies to ligands to be read "
                            "from a multimolecular file.")
    group.add_argument('--esep', dest="entries_sep", type=str,
                       help="the entries separator. Default: ':'. "
                            "E.g.: 3QQK:A:X02:497.")
    group.add_argument('--fsep', dest="fields_sep", type=str,
                       help="the separator used to separate fields "
                            "(protein id, ligand name, molecular file, "
                            "multimol flag) in <ENTRIES_FILE>. "
                            "Default: ','.")

    group = parser.add_argument_group('Paths')
    group.add_argument('-w', '--workdir', dest="working_path", type=str,
                       help="the path where the project and "
                            "its results will be saved.")
    group.add_argument('--pdbdir', dest="pdb_path", type=str,
                       help="path containing local PDB files or to where the "
                            "PDB files will be downloaded. Default: None.")
    group.add_argument('--overwrite', '--no-overwrite', dest="overwrite_path",
                       action=NegateAction, nargs=0,
                       help="overwrite or not an existing project. "
                            " Default: None.")

    group = parser.add_argument_group('Hydrogens')
    group.add_argument('--addh', '--no-addh', dest="add_h",
                       action=NegateAction, nargs=0,
                       help="turn off hydrogen addition.")
    group.add_argument('--ph', type=float,
                       help="add hydrogens appropriate for pH. Default: 7.4.")

    group = parser.add_argument_group('Standardization')
    group.add_argument('--amend', '--no-amend', dest="amend_mol",
                       action=NegateAction, nargs=0,
                       help="fix atomic charges, valence, and bond "
                            "types for small molecules and residues at "
                            "PDB files.")

    group = parser.add_argument_group('Features')
    group.add_argument('-F', '--feat', dest="feat_cfg", type=str,
                       help="a feature definition file (FDef) containing "
                            "chemical and pharmacophoric rules. "
                            "If not provided, the default LUNA features "
                            "definition will be used instead.")

    group = parser.add_argument_group('Interactions')
    group.add_argument('-I', '--inter', dest="inter_cfg", type=str,
                       help="a config file defining interaction parameters. "
                            "If not provided, the default LUNA config "
                            "file will be used instead.")
    group.add_argument('--noncov', '--no-noncov', dest="add_non_cov",
                       action=NegateAction, nargs=0,
                       help="turn on/off non-covalent interactions.")
    group.add_argument('--cov', '--no-cov', dest="add_cov",
                       action=NegateAction, nargs=0,
                       help="turn on/off covalent interactions.")
    group.add_argument('--prox', '--no-prox', dest="add_proximal",
                       action=NegateAction, nargs=0,
                       help="turn on/off proximal contacts.")
    group.add_argument('--aa', '--no-aa', dest="add_atom_atom",
                       action=NegateAction, nargs=0,
                       help="turn on/off atom-atom interactions.")
    group.add_argument('--adi', '--no-adi', dest="add_dependent_inter",
                       action=NegateAction, nargs=0,
                       help="turn on/off dependent interactions.")
    group.add_argument('--awp', '--no-awp',
                       dest="add_h2o_pairs_with_no_target",
                       action=NegateAction, nargs=0,
                       help="accept or not water that doesn't interact "
                            "with the target ligand.")
    group.add_argument('--sdr', '--no-sdr', dest="strict_donor_rules",
                       action=NegateAction, nargs=0,
                       help="turn on/off strict donor rules.")
    group.add_argument('--swdr', '--no-swdr', dest="strict_weak_donor_rules",
                       action=NegateAction, nargs=0,
                       help="turn on/off strict weak donor rules.")
    group.add_argument('--lazy', dest="lazy_comps_list", type=str,
                       help="a comma-separated list of molecule names "
                            "(PDB ids) to ignore explicit hydrogen "
                            "position during the calculation of (weak) "
                            "hydrogen bonds. "
                            "Default: HOH,DOD,WAT,H2O,OH2,NH3,NH4.")

    group = parser.add_argument_group('Interactions filter')
    action = group.add_mutually_exclusive_group()
    action.add_argument('--filter', dest="filter_cfg", type=str,
                        help="a config file defining interaction filters. "
                             "If not provided, all interaction pairs will be "
                             "accepted by default.")

    action.add_argument('--df', dest="default_filter", type=str,
                        choices=['pli', 'ppi', 'pni', 'nni', 'nli'],
                        help="use default filters for interactions, where "
                             "pli, ppi, pni, nni, and nli stand for "
                             "protein-ligand, protein-protein, "
                             "protein-nucleotide, nucleotide-nucleotide, and "
                             "nucleotide-ligand interactions respectively")

    group = parser.add_argument_group('Binding mode filters')
    group.add_argument('--bind', dest="bind_cfg", type=str,
                       help="a config file defining binding mode rules "
                            "to filter interactions.")

    group = parser.add_argument_group('LUNA interaction fingerprints')
    group.add_argument('--ifp', '--no-ifp', dest="ifp",
                       action=NegateAction, nargs=0,
                       help="generate or not LUNA interaction fingerprints.")
    group.add_argument('-L', dest="ifp_num_levels", type=int,
                       help="the number of level defines the number of "
                            "iterations to construct the fingerprint. "
                            "Default: 2.")
    group.add_argument('-R', dest="ifp_radius_step",
                       type=float,
                       help="the radius growth rate defines the multiplier "
                            "to increase the sphere size at each level. "
                            "Default: 5.73171.")
    group.add_argument('-S', dest="ifp_length", type=int,
                       help="the fingerprint length. Default: 4096.")
    group.add_argument('-B', '-C', dest="ifp_count",
                       action=CountIFPAction, nargs=0,
                       help="use bit fingerprints. The default value "
                            "is False, which imply that count fingerprints "
                            "are used instead.")
    group.add_argument('-T', dest="ifp_type", type=str,
                       choices=['EIFP', 'HIFP', 'FIFP'],
                       help="the fingerprint type. Default: EIFP.")
    group.add_argument('--diff', '--no-diff', dest="ifp_diff_comp_classes",
                       action=NegateAction, nargs=0,
                       help="turn on/off differentiation between "
                            "compound classes.")
    group.add_argument('--ifp-out', dest="ifp_output", type=str,
                       help="the fingerprint output file. "
                            "Default: <WORKING_PATH>/results/"
                            "fingerprints/ifp.csv.")
    group.add_argument('--ifp-matrix', dest="ifp_sim_matrix_output",
                       type=str, nargs='?', const=True,
                       help="if provided, save a similarity "
                            "matrix to <IFP_MATRIX>. "
                            "Optionally, use the flag without "
                            "parameters to save the matrix at the "
                            "default path: <WORKING_PATH>/results/"
                            "fingerprints/ifp_sim_matrx.csv")

    group = parser.add_argument_group('Molecular fingerprints')
    group.add_argument('--mfp', '--no-mfp', dest="mfp",
                       action=NegateAction, nargs=0,
                       help="generate or not molecular fingerprints.")
    group.add_argument('--mfp-out', dest="mfp_output", type=str,
                       help="the molecular fingerprint output file. "
                            "Default: <WORKING_PATH>/results/"
                            "fingerprints/mfp.csv.")

    group = parser.add_argument_group('Pymol sessions')
    group.add_argument('--pse', '--no-pse', dest="pse",
                       action=NegateAction, nargs=0,
                       help="export or not interactions to Pymol.")
    group.add_argument('--psedir', dest="pse_path", type=str,
                       help="save PSE files to <PSE_PATH>. "
                            "Default: <WORKING_PATH>/results/pse/.")

    # General options
    parser.add_argument('--lib', dest="mol_obj_type", type=str,
                        choices=['rdkit', 'openbabel'],
                        help="which library (RDKit or Open Babel) to use "
                             "to parse molecules. Default: rdkit.")
    parser.add_argument('--append', '--no-append', dest="append_mode",
                        action=NegateAction, nargs=0,
                        help="skip or not entries when a result for them "
                             "already exists in <WORKING_PATH>.")
    parser.add_argument('--cache', '--no-cache', dest="use_cache",
                        action=NegateAction, nargs=0,
                        help="cache or not protein information"
                             "to save processing time.")
    parser.add_argument("-v", dest="verbosity", type=int,
                        choices=sorted(VERBOSITY_LEVEL.keys()),
                        help="verbosity level. Default: 3.")
    parser.add_argument('--log', '--no-log', dest="logging_enabled",
                        action=NegateAction, nargs=0,
                        help="turn on/off the logging system.")
    parser.add_argument("--nproc", type=int,
                        help="the number of processors to use. "
                             "Default: -1 (use all available CPUs).")
    parser.add_argument('--fork', dest="fork_project", type=str,
                        help="if provided, copy an existing "
                             "project to <WORKING_PATH>.")
    parser.add_argument('-c', '--conf', dest="luna_config", type=str,
                        help="load LUNA parameters from the "
                             "provided config file.")

    return parser


def main():

    conf_parser = argparse.ArgumentParser(add_help=False)
    conf_parser.add_argument("-c", "--conf", dest="luna_config", type=str,
                             help="Load LUNA parameters from the "
                                  "provided config file.")
    args, remaining_argv = conf_parser.parse_known_args()

    # Read the provided config file and set its values to 'config_values'.
    # This dict will then be used to overwrite the default valuies.
    config_values = {}
    if args.luna_config:
        config = configparser.ConfigParser()
        config.read([args.luna_config])

        for section in config.sections():
            config_values.update(dict(config.items(section)))

    # Get the list of expected arguments.
    parser = get_parser()
    # Set default values with parameters provided in the config file.
    parser.set_defaults(**config_values)
    # Overwrite all parameters using values provided through command line.
    args = parser.parse_args(remaining_argv)

    params = dict(args._get_kwargs())

    # Remove all, but params provided through a config file or
    # through command line.
    params = {k: v for k, v in params.items()
              if k in config_values or v is not None}

    # If provided without parameter, this argument will be set to True.
    # That means the default directory should be used instead.
    # However, ProjectParams expects a string. To avoid raising an error,
    # we remove it from 'params' to deal with this argument separatelly.
    ifp_sim_matrix_output = params.pop("ifp_sim_matrix_output", None)

    # Parse params.
    # It will load in default values for every parameter not provided.
    params = ProjectParams(params, fill_defaults=True)

    required = []
    if params.get("entries", None) is None:
        required.append("-e/--entries")

    if params.get("working_path", None) is None:
        required.append("-w/--workdir")

    if required:
        parser.error("the following arguments are required: %s. "
                     "It can also be provided in a config file (-c/--conf)."
                     % ", ".join(required))

    # Fork a project
    if args.fork_project:
        fork_proj_pkl_file = LocalProject.get_project_file(args.fork_project)
        if not os.path.exists(fork_proj_pkl_file):
            raise OSError("No valid project was found at '%s'"
                          % args.fork_project)
        shutil.copytree(args.fork_project, params["working_path"])

    # As the same project can be reutilized using different IFP parameters
    # and binding mode filters, we remove the below args from 'params'
    # so that we first process or load a project, and only then we apply the
    # binding mode filter, calculate IFPs, and PSE sessions.
    binding_mode_filter = params.get("binding_mode_filter", None)
    calc_mfp = params.get("calc_mfp", False)
    calc_ifp = params.get("calc_ifp", False)
    out_pse = params.get("out_pse", False)
    out_pse = params.get("out_pse", False)

    # Set those parameters again using its turned off values.
    params["binding_mode_filter"] = None
    params["calc_mfp"] = False
    params["calc_ifp"] = False
    params["out_pse"] = False

    proj_pkl_file = LocalProject.get_project_file(params["working_path"])
    if not os.path.exists(proj_pkl_file) or params["overwrite_path"]:
        print(u"\u25a8 Creating a new LUNA project: '%s'...\n"
              % params["working_path"])

        pli_obj = LocalProject(**params)
        pli_obj.run()
    else:
        print(u"\u25a8 Reloading existing project: '%s'...\n"
              % params["working_path"])
        pli_obj = LocalProject.load(proj_pkl_file,
                                    logging_enabled=params["logging_enabled"])

    if binding_mode_filter:
        pli_obj.binding_mode_filter = binding_mode_filter

        for agm in pli_obj.atm_grps_mngrs:
            im = InteractionsManager(agm.get_all_interactions(), agm.entry)
            im.filter_out_by_binding_mode(binding_mode_filter)

            entry_results = EntryResults(agm.entry, agm, im, None, None)
            pkl_file = "%s/chunks/%s.pkl.gz" % (pli_obj.working_path,
                                                agm.entry.to_string())
            entry_results.save(pkl_file)

            csv_file = ("%s/results/interactions/%s.csv"
                        % (pli_obj.working_path, agm.entry.to_string()))
            im.to_csv(csv_file)

    if calc_ifp or calc_mfp:
        if calc_ifp:
            ifp_output = params["ifp_output"]
            if ifp_output is None:
                create_directory("%s/results/fingerprints"
                                 % params["working_path"])
                ifp_output = ("%s/results/fingerprints/ifp.csv"
                              % params["working_path"])

            pli_obj.calc_ifp = True
            pli_obj.ifp_num_levels = params["ifp_num_levels"]
            pli_obj.ifp_radius_step = params["ifp_radius_step"]
            pli_obj.ifp_length = params["ifp_length"]
            pli_obj.ifp_type = params["ifp_type"]
            pli_obj.ifp_diff_comp_classes = params["ifp_diff_comp_classes"]
            pli_obj.ifp_count = params["ifp_count"]
            pli_obj.ifp_output = ifp_output

            if ifp_sim_matrix_output:
                if ifp_sim_matrix_output is True:
                    default_output = ("%s/results/fingerprints/"
                                      "ifp_sim_matrix.csv"
                                      % pli_obj.working_path)
                    ifp_sim_matrix_output = default_output

                pli_obj.ifp_sim_matrix_output = ifp_sim_matrix_output

        if calc_mfp:
            mfp_output = params["mfp_output"]
            if mfp_output is None:
                create_directory("%s/results/fingerprints"
                                 % params["working_path"])
                mfp_output = ("%s/results/fingerprints/mfp.csv"
                              % params["working_path"])
            pli_obj.calc_mfp = True

        pli_obj.generate_fps()

    if out_pse:
        pli_obj.out_pse = True

        pli_obj.pse_path = (params["pse_path"]
                            or f"{pli_obj.working_path}/results/pse")
        create_directory(pli_obj.pse_path)

        for ic in pli_obj.interactions_mngrs:
            pse_file = f"{pli_obj.pse_path}/%s.pse" % ic.entry.to_string()
            viewer = InteractionViewer(add_directional_arrows=False)
            viewer.new_session([(ic.entry, ic, pli_obj.pdb_path)], pse_file)

    pli_obj.save_config_file()


if __name__ == '__main__':
    print()
    main()
