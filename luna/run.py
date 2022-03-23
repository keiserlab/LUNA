import argparse
import os
import shutil

from luna.projects import LocalProject
from luna.version import __version__ as version

from luna.mol.entry import MolFileEntry
from luna.interaction.filter import InteractionFilter, BindingModeFilter
from luna.interaction.calc import InteractionCalculator
from luna.interaction.view import InteractionViewer
from luna.interaction.fp.type import IFPType
from luna.util.file import get_filename, create_directory
from luna.interaction.calc import InteractionsManager
from luna.projects import EntryResults

IFP_TYPES = {
    "EIFP": IFPType.EIFP,
    "FIFP": IFPType.FIFP,
    "HIFP": IFPType.HIFP,
}


def get_entries(args):
    with open(args.entries_file, "r") as IN:
        entry_ids = set()
        for row in IN:
            row = row.strip()
            if row == "":
                continue
            entry_ids.add(row)

        pdb_id = get_filename(args.pdb_file)

        for ligand_id in entry_ids:
            yield MolFileEntry.from_mol_file(pdb_id, ligand_id, mol_file=args.ligand_file, is_multimol_file=True)


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('-p', "--prot", dest="pdb_file", type=str, required=True,
                        help="the protein PDB file")

    parser.add_argument('-e', "--entries", dest="entries_file", type=str, required=True,
                        help="an input file containing a list of ligand ids to process")

    parser.add_argument('-l', "--lig", dest="ligand_file", type=str, required=True,
                        help="a molecular file containing 1 or more ligands")

    parser.add_argument('-w', dest="working_path", type=str, required=True,
                        help="the path where the project and its results will be saved")

    parser.add_argument('--out_ifp', dest="out_ifp", action="store_true",
                        help="defines whether it should generate LUNA interaction fingerprints")

    parser.add_argument('-L', dest="ifp_num_levels", type=int, default=2,
                        help="the number of level defines the number of iterations to construct the fingerprint. "
                             "Default: 2")

    parser.add_argument('-R', dest="ifp_radius_step", type=float, default=5.73171,
                        help="the radius growth rate defines the multiplier to increase the sphere size at each level. "
                             "Default: 5.73171")

    parser.add_argument('-S', dest="ifp_length", type=int, default=4096,
                        help="the fingerprint length. Default: 4096")

    parser.add_argument('-T', dest="ifp_type", type=str, default="EIFP", choices=['EIFP', 'HIFP', 'FIFP'],
                        help="the fingerprint type. Default: EIFP")

    parser.add_argument('-B', dest="bit_ifp", action="store_true",
                        help="defines whether it should use bit fingerprints. The default value is False, "
                             "which imply that count fingerprints are used instead")

    parser.add_argument('-O', dest="ifp_output", type=str,
                        help="the fingerprint output file. Default: <WORKING_PATH>/results/fingerprints/ifp.csv")

    parser.add_argument('--sim_matrix_output', type=str, required=False,
                        help="the path where the similarity matrix will be saved. "
                             "If not provided, it won't be generated")

    parser.add_argument('--filter_binding_modes', dest="binding_modes_file", type=str,
                        help="the path of a file containing binding modes to filter")

    parser.add_argument('--out_pse', dest="out_pse", action="store_true",
                        help="defines whether it should export interactions to Pymol")

    parser.add_argument('--pse_path', dest="pse_path", type=str, required=False,
                        help="the path where Pymol sessions (PSE files) will be saved. "
                             "Default: <WORKING_PATH>/results/pse/")

    parser.add_argument('--overwrite', dest="overwrite", action="store_true",
                        help="defines whether it should overwrite an existing project")

    parser.add_argument('-f', '--fork_project', dest="fork_project", type=str, required=False,
                        help="If provided, copy an existing project to <WORKING_PATH>")

    parser.add_argument("--nproc", dest="nproc", type=int,
                        help="the number of processors to use")

    args = parser.parse_args()

    # Fork a project
    if args.fork_project:
        fork_proj_pkl_file = "%s/project_v%s.pkl.gz" % (args.fork_project, version)
        if not os.path.exists(fork_proj_pkl_file):
            raise OSError("No valid project was found at '%s'" % args.fork_project)
        shutil.copytree(args.fork_project, args.working_path)

    proj_pkl_file = "%s/project_v%s.pkl.gz" % (args.working_path, version)
    if not os.path.exists(proj_pkl_file) or args.overwrite:
        print(u"\u25a8 Creating a new LUNA project: '%s'...\n" % args.working_path)

        entries = list(get_entries(args))

        pdb_path = os.path.dirname(os.path.realpath(args.pdb_file))
        ic = InteractionCalculator(inter_filter=InteractionFilter.new_pli_filter(), strict_donor_rules=True, strict_weak_donor_rules=True)

        opts = {}
        opts["working_path"] = args.working_path
        opts["pdb_path"] = pdb_path
        opts["entries"] = entries
        opts["overwrite_path"] = args.overwrite
        opts["inter_calc"] = ic
        opts["mol_obj_type"] = 'rdkit'
        opts["add_h"] = True
        opts["amend_mol"] = True
        opts["calc_ifp"] = False
        opts["out_pse"] = False

        if args.nproc:
            opts["nproc"] = args.nproc

        pli_obj = LocalProject(**opts)
        pli_obj.run()
    else:
        print(u"\u25a8 Reloading existing project: '%s'...\n" % args.working_path)
        pli_obj = LocalProject.load(proj_pkl_file)
    print()

    if args.binding_modes_file:
        binding_mode_filter = BindingModeFilter.from_config_file(args.binding_modes_file)

        for agm in pli_obj.atm_grps_mngrs:
            im = InteractionsManager(agm.get_all_interactions(), agm.entry)
            im.filter_out_by_binding_mode(binding_mode_filter)

            entry_results = EntryResults(agm.entry, agm, im, None, None)
            pkl_file = "%s/chunks/%s.pkl.gz" % (pli_obj.working_path, agm.entry.to_string())
            entry_results.save(pkl_file)

            csv_file = "%s/results/interactions/%s.csv" % (pli_obj.working_path, agm.entry.to_string())
            im.to_csv(csv_file)

    if args.out_ifp:
        ifp_output = args.ifp_output
        if ifp_output is None:
            create_directory("%s/results/fingerprints" % args.working_path)
            ifp_output = "%s/results/fingerprints/ifp.csv" % args.working_path

        pli_obj.calc_ifp = True
        pli_obj.ifp_num_levels = args.ifp_num_levels
        pli_obj.ifp_radius_step = args.ifp_radius_step
        pli_obj.ifp_length = args.ifp_length
        pli_obj.ifp_type = IFP_TYPES[args.ifp_type]

        pli_obj.ifp_count = True
        if args.bit_ifp:
            pli_obj.ifp_count = False

        pli_obj.ifp_output = ifp_output

        if args.sim_matrix_output:
            pli_obj.ifp_sim_matrix_output = args.sim_matrix_output

        pli_obj.generate_ifps()

    if args.out_pse:
        pse_path = args.pse_path or f"{pli_obj.working_path}/results/pse"
        create_directory(pse_path)

        for ic in pli_obj.interactions_mngrs:
            pse_file = f"{pse_path}/%s.pse" % ic.entry.to_string()
            viewer = InteractionViewer(add_directional_arrows=False)
            viewer.new_session([(ic.entry, ic, pli_obj.pdb_path)], pse_file)


if __name__ == '__main__':
    print()
    main()
