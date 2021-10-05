import argparse
import os
import sys
from os import path


scripts = path.abspath(path.join(path.realpath(__file__), '../..'))
sys.path.append(scripts)

from luna import LocalProject
from luna.version import __version__ as version

from luna.mol.entry import MolEntry
from luna.mol.interaction.filter import InteractionFilter, BindingModeFilter
from luna.mol.interaction.calc import InteractionCalculator
from luna.mol.interaction.fp.type import IFPType
from luna.util.file import get_filename, create_directory


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

        target_id = get_filename(args.pdb_file)

        for ligand_id in entry_ids:
            yield MolEntry.from_mol_file(target_id, ligand_id, mol_file=args.ligand_file, is_multimol_file=True)


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

    parser.add_argument('-O', dest="ifp_output", type=str,
                        help="the fingerprint output file. Default: <WORKING_PATH>/results/fingerprints/ifp.csv")

    parser.add_argument('--filter_binding_modes', dest="binding_modes_file", type=str,
                        help="the path of a file containing binding modes to filter")

    parser.add_argument('--overwrite', dest="overwrite", action="store_true",
                        help="defines whether it should overwrite an existing project")

    parser.add_argument('--out_pse', dest="out_pse", action="store_true",
                        help="defines whether it should export interactions to Pymol")

    parser.add_argument('--out_sim', dest="out_sim", action="store_true",
                        help="defines whether it should generate a similarity matrix")

    parser.add_argument("--nproc", dest="nproc", type=int,
                        help="the number of processors to use")

    args = parser.parse_args()

    proj_pkl_file = "%s/project_v%s.pkl.gz" % (args.working_path, version)

    if not os.path.exists(proj_pkl_file) or args.overwrite:
        print(u"\u25a8 Creating a new LUNA project: '%s'...\n" % args.working_path)

        entries = list(get_entries(args))

        pdb_path = os.path.dirname(os.path.realpath(args.pdb_file))
        ic = InteractionCalculator(inter_filter=InteractionFilter.new_pli_filter(), strict_donor_rules=True, strict_weak_donor_rules=True)

        opt = {}
        opt["working_path"] = args.working_path
        opt["pdb_path"] = pdb_path
        opt["entries"] = entries
        opt["overwrite_path"] = args.overwrite
        opt["inter_calc"] = ic
        opt["mol_obj_type"] = 'rdkit'
        opt["try_h_addition"] = True
        opt["amend_mol"] = True
        opt["calc_ifp"] = False
        opt["out_pse"] = args.out_pse

        if args.nproc:
            opt["nproc"] = args.nproc

        if args.binding_modes_file:
            opt["binding_mode_filter"] = BindingModeFilter.from_config_file(args.binding_modes_file)

        pli_obj = LocalProject(**opt)
        pli_obj.run()
    else:
        print(u"\u25a8 Reloading existing project: '%s'...\n" % args.working_path)
        pli_obj = LocalProject.load(proj_pkl_file)
    print()

    ifp_output = args.ifp_output
    if ifp_output is None:
        create_directory("%s/results/fingerprints" % args.working_path)
        ifp_output = "%s/results/fingerprints/ifp.csv" % args.working_path

    pli_obj.calc_ifp = True
    pli_obj.ifp_num_levels = args.ifp_num_levels
    pli_obj.ifp_radius_step = args.ifp_radius_step
    pli_obj.ifp_length = args.ifp_length
    pli_obj.ifp_type = IFP_TYPES[args.ifp_type]
    pli_obj.ifp_output = ifp_output
    pli_obj.out_ifp_sim_matrix = args.out_sim

    pli_obj.generate_ifps()


if __name__ == '__main__':
    print()
    proj = main()
