import unittest
import sys
from os.path import dirname, abspath, exists
from io import StringIO
import pandas as pd

luna_path = dirname(dirname(abspath(__file__)))
sys.path.append(luna_path)

from luna import LocalProject
from luna.mol.entry import *
from luna.interaction.filter import InteractionFilter
from luna.interaction.calc import InteractionCalculator
from luna.interaction.fp.shell import ShellGenerator
from luna.interaction.fp.type import IFPType
from luna.util.file import create_directory, remove_files
from luna.util.default_values import LUNA_PATH
from luna.version import __version__ as version


example_path = f"{LUNA_PATH}/example"
inputs_path = f"{example_path}/inputs"
lig_file = f"{inputs_path}/ligands.mol2"

results_path = f"{LUNA_PATH}/tests/project"
ifp_path = f"{results_path}/results/fingerprints/"

test_ligs = ["ZINC000012442563", "ZINC000343043015"]
expected_entries = set(["protein:%s" % lig_id for lig_id in test_ligs])


class FingerprintTest(unittest.TestCase):

    def _get_project_results(self):

        proj_pkl_file = "%s/project_v%s.pkl.gz" % (results_path, version)
        if not exists(proj_pkl_file):
            print(u"\u25a8 Project not found locally. It will be recreated at: '%s'...\n" % results_path)

            entries = [MolFileEntry.from_mol_file("protein", lig_id, mol_file=lig_file,
                                                  is_multimol_file=True, autoload=True) for lig_id in test_ligs]

            ic = InteractionCalculator(inter_filter=InteractionFilter.new_pli_filter())

            opt = {}
            opt["working_path"] = results_path
            opt["pdb_path"] = inputs_path
            opt["entries"] = entries
            opt["overwrite_path"] = True
            opt["inter_calc"] = ic
            opt["mol_obj_type"] = 'rdkit'
            opt["add_h"] = True
            opt["amend_mol"] = True
            opt["calc_ifp"] = False
            opt["verbosity"] = 1
            opt["logging_enabled"] = False

            pli_obj = LocalProject(**opt)
            pli_obj.run()
        else:
            print(u"\u25a8 Reloading existing project: '%s'...\n" % results_path)
            pli_obj = LocalProject.load(proj_pkl_file, verbosity=1, logging_enabled=False)

        return pli_obj

    def _generate_manual_ifp(self, agm, ifp_num_levels=2, ifp_radius_step=5.73171):
        shells = ShellGenerator(ifp_num_levels, ifp_radius_step, diff_comp_classes=True, ifp_type=IFPType.EIFP)
        sm = shells.create_shells(agm)
        ifp = sm.to_fingerprint(count_fp=True, unique_shells=False, fold_to_length=4096)
        ifp_bits_str = "\t".join([str(idx) for idx in ifp.counts.keys()])
        ifp_count_str = "\t".join([str(count) for count in ifp.counts.values()])

        return "%s,%s,%s" % (agm.entry.to_string(), ifp_bits_str, ifp_count_str)

    def test_automatic_ifps_creation(self):

        pli_obj = self._get_project_results()

        self.assertEqual(expected_entries,
                         set([e.to_string() for e in pli_obj.entries]))

        create_directory(ifp_path)

        ifp_out1 = f"{ifp_path}/ifp_out1.csv"
        if exists(ifp_out1):
            remove_files([ifp_out1])

        self.assertRaises(FileNotFoundError, open, ifp_out1, "r")

        pli_obj.calc_ifp = True
        pli_obj.ifp_output = ifp_out1
        pli_obj.generate_fps()
        self.assertTrue(exists(ifp_out1),
                        "File '%s' was not found." % ifp_out1)

        expected_output = StringIO('ligand_id,on_bits,count\nprotein:ZINC000012442563,83\t84\t87\t135\t157\t180\t187\t211\t273\t306\t349\t414\t520\t542\t592\t598\t604\t615\t668\t705\t764\t807\t820\t882\t911\t926\t931\t1020\t1058\t1095\t1164\t1189\t1210\t1211\t1215\t1218\t1226\t1261\t1327\t1337\t1392\t1394\t1412\t1423\t1471\t1501\t1543\t1561\t1590\t1629\t1653\t1668\t1694\t1719\t1740\t1773\t1837\t1901\t1904\t1912\t1928\t2042\t2062\t2067\t2077\t2093\t2135\t2141\t2210\t2230\t2254\t2352\t2422\t2433\t2441\t2453\t2455\t2487\t2535\t2537\t2593\t2609\t2635\t2657\t2664\t2696\t2733\t2758\t2763\t2769\t2788\t2793\t2798\t2879\t2894\t2914\t2937\t3015\t3021\t3072\t3107\t3109\t3111\t3117\t3126\t3148\t3155\t3163\t3194\t3210\t3257\t3266\t3276\t3280\t3299\t3395\t3431\t3465\t3482\t3491\t3513\t3582\t3642\t3658\t3682\t3683\t3702\t3707\t3722\t3762\t3820\t3823\t3871\t3887\t3896\t3908\t3937\t4003\t4007\t4021\t4035\t4047\t4048\t4056\t4062,2\t4\t1\t1\t27\t1\t1\t1\t1\t1\t1\t1\t1\t1\t45\t2\t1\t2\t3\t2\t4\t5\t1\t53\t20\t2\t1\t5\t1\t2\t38\t1\t37\t9\t3\t2\t1\t3\t2\t1\t1\t1\t1\t53\t1\t1\t1\t1\t1\t4\t2\t1\t1\t2\t1\t1\t42\t1\t2\t1\t2\t2\t1\t18\t2\t2\t2\t42\t3\t9\t41\t1\t13\t5\t1\t1\t1\t7\t44\t1\t1\t3\t4\t1\t1\t1\t1\t2\t1\t1\t1\t3\t2\t37\t1\t1\t2\t1\t2\t1\t1\t32\t5\t5\t3\t1\t1\t1\t3\t4\t1\t18\t1\t1\t2\t9\t1\t2\t1\t38\t1\t3\t1\t1\t1\t1\t1\t2\t5\t1\t1\t1\t1\t1\t39\t5\t3\t1\t1\t1\t1\t4\t3\t3\t39\nprotein:ZINC000343043015,3\t36\t83\t84\t157\t180\t211\t243\t349\t389\t423\t542\t547\t573\t592\t596\t598\t604\t615\t629\t668\t688\t753\t764\t797\t807\t820\t882\t911\t919\t926\t931\t962\t988\t1012\t1020\t1095\t1164\t1210\t1211\t1213\t1215\t1218\t1261\t1300\t1327\t1380\t1394\t1412\t1423\t1424\t1492\t1543\t1561\t1629\t1653\t1719\t1773\t1837\t1899\t1912\t2042\t2052\t2067\t2093\t2115\t2138\t2141\t2149\t2210\t2230\t2254\t2311\t2319\t2346\t2352\t2422\t2427\t2433\t2453\t2487\t2535\t2537\t2625\t2635\t2654\t2664\t2755\t2758\t2763\t2788\t2793\t2879\t2914\t2924\t2937\t3015\t3021\t3095\t3109\t3111\t3117\t3126\t3185\t3210\t3266\t3299\t3318\t3395\t3431\t3465\t3482\t3491\t3513\t3516\t3582\t3707\t3722\t3887\t3896\t3908\t3937\t3999\t4047\t4048\t4056\t4062,1\t1\t2\t3\t25\t2\t1\t1\t2\t1\t1\t2\t1\t1\t35\t1\t2\t2\t2\t1\t3\t1\t1\t3\t1\t5\t1\t49\t19\t1\t2\t2\t1\t1\t1\t4\t1\t33\t34\t8\t1\t2\t1\t2\t1\t2\t1\t1\t1\t49\t1\t1\t1\t1\t3\t2\t3\t1\t37\t2\t1\t3\t1\t14\t1\t1\t1\t37\t1\t3\t9\t36\t1\t1\t1\t1\t12\t2\t5\t1\t7\t40\t1\t1\t3\t1\t1\t1\t2\t1\t1\t2\t33\t1\t1\t1\t1\t2\t1\t30\t4\t4\t3\t1\t3\t14\t3\t1\t5\t1\t3\t2\t33\t1\t1\t2\t1\t4\t1\t34\t4\t2\t1\t3\t3\t2\t33\n')
        expected_df = pd.read_csv(expected_output)
        curr_df = pd.read_csv(ifp_out1)
        self.assertTrue(expected_df.equals(curr_df))

        pli_obj.calc_ifp = True
        pli_obj.ifp_length = 4
        pli_obj.generate_fps()
        curr_df = pd.read_csv(ifp_out1)
        self.assertFalse(expected_df.equals(curr_df))

    def test_manual_ifps_creation(self):

        pli_obj = self._get_project_results()

        self.assertEqual(expected_entries,
                         set([e.to_string() for e in pli_obj.entries]))

        expected_results = {"atm_grps_len": [("protein:ZINC000012442563",
                                              485),
                                             ("protein:ZINC000343043015",
                                              422)],
                            "inter_len": [("protein:ZINC000012442563", 73), ("protein:ZINC000343043015", 52)],
                            "ifps": ['protein:ZINC000012442563,83\t84\t87\t135\t157\t180\t187\t211\t273\t306\t349\t414\t520\t542\t592\t598\t604\t615\t668\t705\t764\t807\t820\t882\t911\t926\t931\t1020\t1058\t1095\t1164\t1189\t1210\t1211\t1215\t1218\t1226\t1261\t1327\t1337\t1392\t1394\t1412\t1423\t1471\t1501\t1543\t1561\t1590\t1629\t1653\t1668\t1694\t1719\t1740\t1773\t1837\t1901\t1904\t1912\t1928\t2042\t2062\t2067\t2077\t2093\t2135\t2141\t2210\t2230\t2254\t2352\t2422\t2433\t2441\t2453\t2455\t2487\t2535\t2537\t2593\t2609\t2635\t2657\t2664\t2696\t2733\t2758\t2763\t2769\t2788\t2793\t2798\t2879\t2894\t2914\t2937\t3015\t3021\t3072\t3107\t3109\t3111\t3117\t3126\t3148\t3155\t3163\t3194\t3210\t3257\t3266\t3276\t3280\t3299\t3395\t3431\t3465\t3482\t3491\t3513\t3582\t3642\t3658\t3682\t3683\t3702\t3707\t3722\t3762\t3820\t3823\t3871\t3887\t3896\t3908\t3937\t4003\t4007\t4021\t4035\t4047\t4048\t4056\t4062,2\t4\t1\t1\t27\t1\t1\t1\t1\t1\t1\t1\t1\t1\t45\t2\t1\t2\t3\t2\t4\t5\t1\t53\t20\t2\t1\t5\t1\t2\t38\t1\t37\t9\t3\t2\t1\t3\t2\t1\t1\t1\t1\t53\t1\t1\t1\t1\t1\t4\t2\t1\t1\t2\t1\t1\t42\t1\t2\t1\t2\t2\t1\t18\t2\t2\t2\t42\t3\t9\t41\t1\t13\t5\t1\t1\t1\t7\t44\t1\t1\t3\t4\t1\t1\t1\t1\t2\t1\t1\t1\t3\t2\t37\t1\t1\t2\t1\t2\t1\t1\t32\t5\t5\t3\t1\t1\t1\t3\t4\t1\t18\t1\t1\t2\t9\t1\t2\t1\t38\t1\t3\t1\t1\t1\t1\t1\t2\t5\t1\t1\t1\t1\t1\t39\t5\t3\t1\t1\t1\t1\t4\t3\t3\t39',
                                     'protein:ZINC000343043015,3\t36\t83\t84\t157\t180\t211\t243\t349\t389\t423\t542\t547\t573\t592\t596\t598\t604\t615\t629\t668\t688\t753\t764\t797\t807\t820\t882\t911\t919\t926\t931\t962\t988\t1012\t1020\t1095\t1164\t1210\t1211\t1213\t1215\t1218\t1261\t1300\t1327\t1380\t1394\t1412\t1423\t1424\t1492\t1543\t1561\t1629\t1653\t1719\t1773\t1837\t1899\t1912\t2042\t2052\t2067\t2093\t2115\t2138\t2141\t2149\t2210\t2230\t2254\t2311\t2319\t2346\t2352\t2422\t2427\t2433\t2453\t2487\t2535\t2537\t2625\t2635\t2654\t2664\t2755\t2758\t2763\t2788\t2793\t2879\t2914\t2924\t2937\t3015\t3021\t3095\t3109\t3111\t3117\t3126\t3185\t3210\t3266\t3299\t3318\t3395\t3431\t3465\t3482\t3491\t3513\t3516\t3582\t3707\t3722\t3887\t3896\t3908\t3937\t3999\t4047\t4048\t4056\t4062,1\t1\t2\t3\t25\t2\t1\t1\t2\t1\t1\t2\t1\t1\t35\t1\t2\t2\t2\t1\t3\t1\t1\t3\t1\t5\t1\t49\t19\t1\t2\t2\t1\t1\t1\t4\t1\t33\t34\t8\t1\t2\t1\t2\t1\t2\t1\t1\t1\t49\t1\t1\t1\t1\t3\t2\t3\t1\t37\t2\t1\t3\t1\t14\t1\t1\t1\t37\t1\t3\t9\t36\t1\t1\t1\t1\t12\t2\t5\t1\t7\t40\t1\t1\t3\t1\t1\t1\t2\t1\t1\t2\t33\t1\t1\t1\t1\t2\t1\t30\t4\t4\t3\t1\t3\t14\t3\t1\t5\t1\t3\t2\t33\t1\t1\t2\t1\t4\t1\t34\t4\t2\t1\t3\t3\t2\t33']}

        results = list(pli_obj.results)

        atm_grps_len = [(r.entry.to_string(), len(r.atm_grps_mngr)) for r in results]
        self.assertEqual(expected_results["atm_grps_len"], atm_grps_len)

        inter_len = [(r.entry.to_string(), len(r.interactions_mngr)) for r in results]
        self.assertEqual(expected_results["inter_len"], inter_len)

        idx = 0
        ifp_as_str1 = self._generate_manual_ifp(results[idx].atm_grps_mngr)
        self.assertEqual(expected_results["ifps"][idx], ifp_as_str1)

        idx = 1
        ifp_as_str2 = self._generate_manual_ifp(results[idx].atm_grps_mngr)
        self.assertEqual(expected_results["ifps"][idx], ifp_as_str2)

        self.assertNotEqual(ifp_as_str1, ifp_as_str2)

        #
        # Must fail as it's using a different number of levels.
        idx = 0
        ifp_as_str = self._generate_manual_ifp(results[idx].atm_grps_mngr, ifp_num_levels=3)
        self.assertNotEqual(expected_results["ifps"][idx], ifp_as_str)

        idx = 1
        ifp_as_str = self._generate_manual_ifp(results[idx].atm_grps_mngr, ifp_num_levels=3)
        self.assertNotEqual(expected_results["ifps"][idx], ifp_as_str)

        idx = 0
        ifp_as_str = self._generate_manual_ifp(results[idx].atm_grps_mngr, ifp_num_levels=5, ifp_radius_step=1)
        self.assertNotEqual(expected_results["ifps"][idx], ifp_as_str)

        idx = 1
        ifp_as_str = self._generate_manual_ifp(results[idx].atm_grps_mngr, ifp_num_levels=5, ifp_radius_step=1)
        self.assertNotEqual(expected_results["ifps"][idx], ifp_as_str)


if __name__ == '__main__':
    unittest.main()
