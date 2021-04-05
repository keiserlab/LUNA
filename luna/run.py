import argparse
import os

from luna import LocalProject
from luna.version import __version__ as version

from luna.mol.entry import MolEntry
from luna.mol.interaction.filter import InteractionFilter
from luna.mol.interaction.calc import InteractionCalculator
from luna.mol.interaction.fp.type import IFPType

from luna.util.file import get_filename


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
                        help="the number of level defines the number of iterations to construct the fingerprint")

    parser.add_argument('-R', dest="ifp_radius_step", type=float, default=5.73171,
                        help="the radius growth rate defines the multiplier to increase the sphere size at each level")

    parser.add_argument('-S', dest="ifp_length", type=int, default=4096,
                        help="the fingerprint length")

    parser.add_argument('-T', dest="ifp_type", type=str, default="EIFP", choices=['EIFP', 'HIFP', 'FIFP'],
                        help="the fingerprint type")

    args = parser.parse_args()

    proj_pkl_file = "%s/project_v%s.pkl.gz" % (args.working_path, version)
    if not os.path.exists(proj_pkl_file):
        print(u"\u25a8 Creating a new LUNA project: '%s'...\n" % args.working_path)

        entries = list(get_entries(args))

        pdb_path = os.path.dirname(os.path.realpath(args.pdb_file))
        ic = InteractionCalculator(inter_filter=InteractionFilter.new_pli_filter(), strict_donor_rules=True, strict_weak_donor_rules=True)

        opt = {}
        opt["working_path"] = args.working_path
        opt["pdb_path"] = pdb_path
        opt["entries"] = entries
        opt["overwrite_path"] = False
        opt["inter_calc"] = ic
        opt["mol_obj_type"] = 'rdkit'
        opt["try_h_addition"] = True
        opt["amend_mol"] = True

        pli_obj = LocalProject(**opt)
        pli_obj.run()
    else:
        print(u"\u25a8 Reloading existing project: '%s'...\n" % args.working_path)
        pli_obj = LocalProject.load(proj_pkl_file)
    print()

    pli_obj.calc_ifp = True
    pli_obj.ifp_num_levels = args.ifp_num_levels
    pli_obj.ifp_radius_step = args.ifp_radius_step
    pli_obj.ifp_length = args.ifp_length
    pli_obj.ifp_type = IFP_TYPES[args.ifp_type]

    pli_obj.generate_ifps()


if __name__ == '__main__':
    print()
    proj = main()
