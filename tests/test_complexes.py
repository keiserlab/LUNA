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

project_path = f"{LUNA_PATH}/tests/projects/test_complexes"
ifp_file = f"{project_path}/results/fingerprints/ifp.csv"

complexes = ["6JWU:A:CA:806", "1I2C:A:UPG:402", "1KVT:A:NAD:340",
             "3MJM:A:DOR:1410", "3QL8:A:X01:300", "6GST:A:GSH:218",
             "1USN:A:ZN:258", "4B3U:A:NWL:400", "3I4Y:A:35C:1",
             "9icc:A:CR:339", "1THA:B:HOH:137"]


class ComplexesTest(unittest.TestCase):

    def _get_project_results(self):

        proj_pkl_file = "%s/project_v%s.pkl.gz" % (project_path, version)
        if not exists(proj_pkl_file):
            print(u"\u25a8 Project not found locally. "
                  "It will be recreated at: '%s'...\n" % project_path)

            entries = [MolEntry.from_string(entry)
                       for entry in complexes]

            inter_filter = InteractionFilter.new_pli_filter()
            ic = InteractionCalculator(inter_filter=inter_filter)

            opt = {}
            opt["working_path"] = project_path
            opt["pdb_path"] = inputs_path
            opt["entries"] = entries
            opt["overwrite_path"] = True
            opt["inter_calc"] = ic
            opt["add_h"] = True
            opt["amend_mol"] = True
            opt["calc_ifp"] = True
            opt["verbosity"] = 1
            opt["logging_enabled"] = False

            pli_obj = LocalProject(**opt)
            pli_obj.run()
        else:
            print(u"\u25a8 Reloading existing project: '%s'...\n"
                  % project_path)
            pli_obj = LocalProject.load(proj_pkl_file,
                                        verbosity=1,
                                        logging_enabled=False)

        return pli_obj

    def test_complex_results(self):

        pli_obj = self._get_project_results()

        self.assertEqual(set(complexes),
                         set([e.to_string() for e in pli_obj.entries]))

        self.assertTrue(exists(ifp_file),
                        "File '%s' was not found." % ifp_file)

        expected_output = StringIO("ligand_id,on_bits,count\n6JWU:A:CA:806,84\t157\t172\t209\t411\t585\t592\t598\t615\t809\t882\t911\t1020\t1026\t1058\t1164\t1210\t1211\t1215\t1218\t1305\t1327\t1423\t1452\t1528\t1629\t1653\t1837\t1934\t1992\t2067\t2070\t2093\t2141\t2210\t2215\t2230\t2254\t2422\t2487\t2535\t2635\t2788\t2879\t2898\t3021\t3109\t3111\t3117\t3126\t3201\t3210\t3266\t3491\t3534\t3707\t3722\t3896\t4047\t4048\t4056\t4062,3\t15\t1\t1\t1\t1\t9\t1\t1\t2\t15\t6\t1\t2\t1\t13\t12\t3\t5\t1\t2\t1\t15\t1\t1\t1\t2\t16\t1\t1\t3\t2\t1\t15\t3\t1\t1\t16\t6\t1\t14\t1\t1\t12\t1\t1\t20\t1\t1\t3\t1\t3\t3\t13\t1\t1\t1\t12\t1\t2\t3\t9\n1I2C:A:UPG:402,61\t83\t84\t90\t109\t111\t138\t157\t167\t172\t176\t187\t211\t212\t223\t225\t273\t284\t349\t353\t388\t401\t414\t435\t470\t503\t532\t536\t542\t575\t584\t592\t598\t604\t615\t618\t658\t668\t677\t692\t700\t717\t727\t750\t764\t807\t820\t825\t861\t868\t882\t896\t911\t926\t931\t941\t962\t964\t987\t1013\t1020\t1026\t1042\t1051\t1057\t1058\t1068\t1090\t1095\t1108\t1156\t1164\t1185\t1210\t1211\t1215\t1218\t1233\t1257\t1261\t1270\t1305\t1312\t1324\t1327\t1340\t1354\t1412\t1423\t1424\t1459\t1492\t1527\t1528\t1540\t1543\t1561\t1576\t1577\t1599\t1606\t1629\t1644\t1653\t1672\t1700\t1708\t1719\t1736\t1747\t1759\t1784\t1837\t1865\t1868\t1888\t1899\t1904\t1930\t1938\t1968\t2042\t2055\t2062\t2067\t2070\t2072\t2077\t2085\t2093\t2134\t2138\t2141\t2159\t2193\t2210\t2230\t2237\t2245\t2254\t2259\t2287\t2314\t2319\t2326\t2344\t2352\t2422\t2433\t2442\t2460\t2465\t2482\t2483\t2487\t2492\t2528\t2534\t2535\t2537\t2551\t2585\t2593\t2635\t2659\t2673\t2674\t2742\t2752\t2774\t2788\t2793\t2842\t2879\t2894\t2909\t2921\t2943\t2952\t2961\t2980\t3015\t3021\t3031\t3070\t3072\t3109\t3111\t3115\t3117\t3126\t3127\t3132\t3137\t3194\t3196\t3207\t3210\t3216\t3218\t3236\t3246\t3250\t3266\t3269\t3299\t3395\t3414\t3430\t3431\t3465\t3469\t3478\t3487\t3491\t3497\t3499\t3508\t3510\t3540\t3551\t3562\t3582\t3591\t3609\t3642\t3668\t3705\t3707\t3722\t3742\t3743\t3770\t3798\t3813\t3821\t3838\t3859\t3896\t3908\t3915\t3937\t3982\t4014\t4020\t4047\t4048\t4056\t4062\t4063\t4092,1\t3\t9\t1\t1\t1\t1\t24\t1\t2\t1\t3\t1\t1\t1\t1\t2\t1\t1\t1\t1\t1\t2\t1\t1\t1\t2\t1\t1\t1\t1\t29\t3\t1\t1\t1\t2\t2\t1\t1\t1\t4\t1\t1\t2\t1\t1\t1\t2\t1\t56\t1\t20\t3\t1\t1\t1\t1\t1\t1\t4\t2\t1\t9\t1\t4\t1\t1\t3\t1\t1\t43\t2\t44\t7\t19\t2\t1\t8\t4\t1\t1\t1\t1\t3\t2\t4\t1\t56\t1\t1\t1\t1\t3\t2\t1\t4\t1\t1\t1\t1\t2\t1\t2\t1\t1\t1\t1\t1\t1\t1\t3\t52\t1\t1\t1\t1\t4\t3\t1\t1\t1\t1\t2\t12\t1\t1\t2\t2\t8\t3\t1\t51\t1\t1\t2\t7\t1\t4\t48\t2\t6\t4\t3\t1\t1\t7\t10\t1\t1\t1\t2\t1\t1\t5\t4\t1\t1\t42\t8\t1\t2\t1\t2\t2\t1\t17\t1\t1\t1\t1\t1\t1\t35\t4\t1\t2\t1\t1\t1\t1\t1\t3\t1\t1\t2\t46\t4\t1\t4\t2\t1\t6\t1\t6\t1\t1\t9\t2\t1\t1\t1\t1\t9\t1\t1\t7\t1\t1\t1\t1\t2\t1\t2\t46\t1\t1\t1\t10\t1\t1\t1\t1\t1\t1\t2\t1\t1\t7\t3\t1\t1\t1\t5\t2\t1\t1\t2\t40\t6\t1\t1\t1\t1\t1\t2\t1\t6\t27\t1\t1\n1KVT:A:NAD:340,7\t60\t83\t84\t92\t93\t138\t142\t147\t154\t157\t172\t187\t191\t211\t216\t225\t263\t273\t286\t353\t370\t379\t414\t481\t508\t518\t527\t532\t534\t590\t592\t596\t598\t615\t618\t658\t668\t676\t677\t700\t705\t717\t755\t764\t794\t796\t802\t807\t820\t836\t882\t900\t911\t926\t964\t987\t1018\t1020\t1023\t1026\t1051\t1055\t1058\t1060\t1068\t1090\t1095\t1099\t1154\t1164\t1210\t1211\t1215\t1218\t1228\t1257\t1261\t1263\t1286\t1305\t1327\t1340\t1349\t1376\t1393\t1412\t1423\t1488\t1527\t1528\t1540\t1543\t1552\t1561\t1580\t1606\t1629\t1653\t1669\t1704\t1708\t1736\t1747\t1773\t1784\t1792\t1837\t1861\t1870\t1899\t1904\t1912\t1930\t1972\t1992\t1994\t2062\t2067\t2070\t2077\t2085\t2093\t2138\t2141\t2146\t2159\t2210\t2230\t2245\t2251\t2254\t2259\t2272\t2287\t2314\t2319\t2330\t2334\t2352\t2377\t2391\t2422\t2433\t2442\t2446\t2485\t2487\t2492\t2520\t2535\t2537\t2585\t2591\t2635\t2643\t2648\t2664\t2674\t2720\t2740\t2781\t2788\t2789\t2798\t2805\t2840\t2845\t2846\t2879\t2885\t2893\t2894\t2898\t2914\t2921\t2935\t2953\t2957\t2981\t3015\t3021\t3072\t3082\t3109\t3111\t3115\t3117\t3121\t3126\t3132\t3156\t3161\t3188\t3194\t3203\t3210\t3212\t3216\t3221\t3227\t3231\t3246\t3253\t3261\t3266\t3269\t3275\t3299\t3342\t3375\t3395\t3416\t3431\t3469\t3478\t3491\t3508\t3510\t3534\t3538\t3562\t3573\t3597\t3631\t3632\t3642\t3654\t3669\t3687\t3695\t3705\t3707\t3722\t3785\t3798\t3813\t3880\t3893\t3896\t3908\t3910\t4013\t4038\t4047\t4048\t4056\t4062\t4084\t4092,1\t1\t1\t2\t1\t1\t2\t1\t1\t1\t19\t3\t1\t1\t1\t1\t1\t1\t1\t1\t1\t2\t1\t2\t1\t1\t1\t1\t1\t1\t1\t29\t1\t3\t2\t1\t5\t1\t1\t1\t3\t1\t4\t1\t1\t2\t1\t1\t3\t1\t1\t52\t2\t15\t1\t2\t2\t1\t1\t1\t2\t9\t1\t6\t1\t2\t1\t3\t1\t1\t41\t41\t7\t21\t3\t1\t8\t4\t1\t1\t2\t3\t1\t1\t1\t1\t1\t52\t2\t1\t1\t2\t1\t1\t1\t1\t1\t3\t5\t1\t1\t1\t1\t1\t1\t1\t1\t51\t1\t1\t1\t2\t1\t2\t1\t1\t1\t1\t10\t2\t1\t2\t7\t1\t48\t1\t1\t3\t8\t2\t1\t50\t2\t1\t4\t4\t3\t1\t1\t7\t1\t1\t10\t4\t1\t1\t1\t6\t4\t1\t40\t3\t2\t1\t3\t1\t1\t1\t17\t1\t1\t1\t1\t2\t1\t1\t2\t1\t1\t33\t1\t2\t4\t1\t1\t2\t2\t1\t1\t1\t2\t3\t1\t1\t45\t1\t3\t3\t1\t3\t6\t1\t1\t1\t3\t1\t2\t1\t1\t1\t1\t1\t1\t1\t1\t10\t1\t2\t1\t1\t1\t7\t1\t1\t1\t4\t48\t1\t10\t2\t1\t3\t1\t1\t1\t1\t1\t2\t1\t1\t1\t1\t3\t3\t1\t4\t1\t1\t1\t32\t6\t1\t1\t2\t3\t2\t10\t24\t1\t1\n3MJM:A:DOR:1410,84\t126\t157\t178\t187\t196\t211\t273\t361\t412\t414\t471\t487\t512\t592\t598\t615\t668\t755\t764\t778\t814\t820\t832\t841\t847\t864\t882\t911\t1058\t1065\t1085\t1096\t1103\t1128\t1164\t1210\t1211\t1215\t1261\t1301\t1327\t1354\t1364\t1404\t1412\t1423\t1424\t1435\t1466\t1508\t1543\t1561\t1580\t1602\t1636\t1653\t1732\t1752\t1833\t1835\t1837\t1872\t1899\t1904\t1934\t1976\t1988\t2018\t2025\t2055\t2062\t2067\t2077\t2093\t2114\t2141\t2145\t2210\t2228\t2230\t2254\t2259\t2314\t2352\t2416\t2422\t2465\t2476\t2487\t2519\t2535\t2537\t2611\t2615\t2664\t2674\t2702\t2724\t2753\t2772\t2787\t2788\t2879\t2886\t2898\t2912\t2983\t3009\t3015\t3021\t3072\t3085\t3095\t3109\t3116\t3117\t3120\t3126\t3143\t3194\t3210\t3212\t3246\t3266\t3298\t3382\t3491\t3556\t3591\t3642\t3645\t3722\t3740\t3746\t3896\t3902\t3908\t3993\t4000\t4048\t4056\t4062,2\t1\t7\t1\t2\t1\t1\t2\t1\t1\t2\t1\t1\t1\t18\t2\t1\t4\t1\t5\t2\t1\t5\t1\t1\t1\t1\t28\t4\t1\t1\t1\t1\t1\t1\t19\t18\t1\t12\t4\t1\t2\t1\t1\t1\t1\t28\t1\t1\t1\t1\t5\t1\t2\t3\t1\t1\t2\t1\t1\t1\t23\t2\t1\t4\t1\t1\t1\t1\t1\t1\t2\t9\t3\t2\t2\t21\t1\t2\t1\t2\t22\t2\t2\t5\t1\t2\t1\t1\t1\t1\t17\t1\t3\t2\t2\t1\t1\t1\t3\t1\t1\t1\t17\t1\t1\t1\t1\t1\t3\t2\t2\t1\t2\t18\t1\t2\t2\t2\t1\t6\t2\t1\t2\t8\t2\t2\t20\t1\t1\t3\t1\t1\t2\t1\t19\t1\t2\t1\t1\t1\t2\t13\n3QL8:A:X01:300,54\t55\t84\t109\t138\t145\t157\t212\t328\t384\t481\t532\t543\t592\t609\t659\t668\t764\t788\t820\t830\t878\t882\t911\t920\t1020\t1026\t1058\t1095\t1164\t1205\t1210\t1211\t1215\t1218\t1305\t1412\t1423\t1435\t1446\t1528\t1543\t1560\t1564\t1602\t1628\t1629\t1644\t1653\t1671\t1681\t1721\t1745\t1837\t1894\t1899\t2067\t2070\t2093\t2141\t2210\t2245\t2254\t2422\t2485\t2490\t2506\t2535\t2635\t2788\t2840\t2879\t2921\t2931\t2950\t3000\t3109\t3111\t3126\t3127\t3208\t3210\t3266\t3309\t3393\t3395\t3467\t3478\t3491\t3510\t3520\t3534\t3550\t3570\t3634\t3707\t3785\t3832\t3864\t3896\t3900\t3908\t4047\t4048\t4056\t4062,1\t1\t3\t1\t1\t1\t16\t1\t1\t1\t1\t1\t1\t12\t1\t1\t1\t1\t1\t1\t1\t1\t29\t12\t1\t2\t2\t3\t3\t20\t1\t19\t5\t7\t1\t2\t1\t29\t1\t1\t2\t1\t1\t1\t1\t3\t2\t1\t3\t1\t1\t1\t1\t25\t1\t1\t3\t2\t3\t25\t3\t1\t25\t8\t1\t1\t1\t23\t2\t2\t1\t16\t1\t2\t1\t1\t24\t2\t3\t1\t1\t3\t3\t1\t1\t4\t1\t1\t22\t1\t1\t2\t1\t1\t1\t3\t1\t1\t1\t20\t1\t5\t2\t3\t6\t9\n6GST:A:GSH:218,49\t83\t84\t157\t172\t187\t211\t223\t261\t273\t349\t414\t446\t466\t530\t532\t542\t592\t598\t604\t615\t672\t742\t764\t807\t827\t882\t911\t926\t931\t982\t1026\t1054\t1058\t1065\t1122\t1164\t1180\t1208\t1210\t1211\t1215\t1261\t1305\t1327\t1354\t1423\t1424\t1428\t1508\t1528\t1534\t1561\t1581\t1583\t1653\t1665\t1673\t1719\t1732\t1821\t1833\t1837\t1846\t1848\t1866\t1881\t1904\t1912\t1934\t1979\t1989\t1991\t2021\t2042\t2062\t2067\t2070\t2077\t2093\t2141\t2148\t2210\t2230\t2254\t2314\t2352\t2397\t2422\t2433\t2487\t2535\t2537\t2542\t2553\t2566\t2669\t2733\t2757\t2786\t2792\t2793\t2839\t2847\t2879\t3015\t3021\t3072\t3085\t3095\t3109\t3115\t3117\t3120\t3125\t3126\t3194\t3210\t3266\t3272\t3299\t3349\t3384\t3431\t3465\t3478\t3491\t3494\t3582\t3642\t3654\t3705\t3707\t3722\t3844\t3872\t3896\t3902\t3937\t3988\t4007\t4050\t4056\t4062\t4071,1\t1\t2\t17\t3\t2\t2\t1\t1\t2\t2\t2\t1\t1\t1\t3\t2\t14\t1\t2\t1\t1\t1\t2\t1\t1\t29\t9\t1\t1\t1\t1\t1\t4\t2\t1\t20\t1\t1\t19\t2\t23\t4\t1\t1\t1\t29\t1\t1\t1\t1\t1\t1\t1\t1\t3\t1\t1\t1\t2\t1\t1\t26\t1\t1\t1\t1\t4\t1\t1\t1\t1\t1\t1\t2\t2\t8\t1\t2\t5\t25\t1\t1\t3\t25\t2\t4\t1\t8\t1\t1\t21\t1\t1\t1\t1\t1\t1\t1\t1\t1\t2\t1\t1\t15\t4\t1\t2\t2\t2\t25\t1\t3\t1\t1\t1\t6\t2\t8\t1\t2\t1\t1\t1\t2\t3\t24\t1\t2\t2\t1\t1\t5\t3\t1\t1\t17\t2\t2\t1\t1\t1\t7\t14\t1\n1USN:A:ZN:258,157\t544\t592\t598\t615\t668\t764\t820\t882\t911\t1164\t1210\t1211\t1215\t1218\t1327\t1423\t1543\t1602\t1629\t1653\t1837\t2067\t2077\t2141\t2197\t2210\t2230\t2254\t2422\t2487\t2535\t2635\t2879\t3021\t3109\t3117\t3126\t3266\t3382\t3491\t3642\t3722\t3740\t3896\t4047\t4048\t4056\t4062,4\t1\t19\t2\t2\t3\t3\t3\t15\t2\t9\t9\t1\t3\t1\t2\t15\t3\t3\t1\t1\t10\t8\t1\t10\t1\t1\t2\t10\t2\t2\t9\t1\t9\t2\t9\t2\t1\t8\t3\t9\t1\t2\t1\t10\t1\t1\t1\t19\n4B3U:A:NWL:400,6\t84\t94\t106\t109\t111\t126\t157\t263\t308\t349\t360\t532\t542\t586\t592\t595\t598\t604\t615\t668\t689\t710\t712\t764\t807\t820\t882\t910\t911\t931\t945\t962\t1020\t1026\t1058\t1116\t1164\t1210\t1211\t1215\t1218\t1261\t1265\t1305\t1327\t1354\t1364\t1394\t1412\t1423\t1437\t1441\t1454\t1543\t1561\t1575\t1602\t1629\t1653\t1719\t1740\t1837\t1902\t1912\t1924\t2042\t2067\t2070\t2077\t2093\t2141\t2159\t2210\t2230\t2254\t2259\t2319\t2352\t2353\t2371\t2405\t2422\t2433\t2437\t2487\t2522\t2535\t2537\t2601\t2635\t2664\t2758\t2772\t2788\t2840\t2879\t2914\t2921\t2931\t2937\t3021\t3095\t3109\t3111\t3117\t3126\t3159\t3210\t3216\t3246\t3266\t3299\t3393\t3395\t3455\t3465\t3478\t3488\t3491\t3496\t3510\t3552\t3562\t3642\t3654\t3655\t3707\t3722\t3748\t3770\t3887\t3896\t3902\t3908\t3940\t4018\t4019\t4047\t4048\t4051\t4055\t4056\t4062,1\t3\t1\t1\t1\t1\t1\t17\t1\t1\t1\t1\t1\t1\t1\t21\t1\t3\t1\t3\t2\t1\t1\t1\t2\t1\t2\t26\t1\t9\t1\t1\t1\t1\t1\t1\t1\t19\t21\t3\t6\t2\t4\t1\t1\t3\t1\t1\t1\t2\t26\t1\t1\t1\t2\t1\t1\t1\t1\t2\t1\t1\t24\t1\t1\t1\t1\t9\t1\t1\t2\t23\t2\t2\t4\t23\t2\t1\t2\t1\t1\t1\t6\t1\t1\t5\t1\t22\t1\t1\t1\t1\t2\t1\t2\t1\t19\t1\t1\t1\t1\t3\t1\t25\t1\t3\t2\t1\t3\t1\t1\t9\t1\t1\t5\t1\t1\t2\t1\t21\t1\t3\t1\t2\t1\t1\t1\t2\t3\t1\t1\t1\t18\t1\t5\t1\t1\t1\t1\t1\t1\t1\t3\t19\n3I4Y:A:35C:1,31\t84\t135\t157\t172\t187\t273\t414\t457\t532\t592\t598\t601\t615\t668\t677\t764\t820\t882\t911\t1020\t1058\t1164\t1188\t1194\t1210\t1211\t1215\t1265\t1327\t1423\t1543\t1561\t1602\t1653\t1799\t1837\t1838\t1904\t1912\t2062\t2067\t2077\t2093\t2141\t2181\t2182\t2210\t2228\t2230\t2250\t2254\t2314\t2328\t2355\t2422\t2426\t2487\t2535\t2537\t2604\t2668\t2735\t2737\t2781\t2788\t2793\t2879\t2898\t3021\t3072\t3076\t3086\t3109\t3111\t3115\t3117\t3126\t3194\t3210\t3218\t3266\t3395\t3466\t3478\t3491\t3582\t3642\t3647\t3677\t3707\t3722\t3736\t3896\t3908\t3937\t3974\t4048\t4056\t4062,1\t5\t1\t16\t2\t3\t3\t3\t1\t1\t21\t3\t2\t2\t2\t2\t3\t2\t35\t13\t2\t2\t24\t1\t1\t22\t5\t5\t1\t3\t35\t2\t1\t1\t3\t1\t30\t1\t6\t1\t3\t11\t4\t3\t29\t1\t1\t4\t1\t2\t2\t29\t2\t1\t1\t11\t1\t2\t25\t1\t1\t1\t1\t1\t1\t3\t1\t21\t2\t2\t3\t2\t1\t22\t2\t1\t4\t3\t9\t5\t1\t11\t2\t1\t2\t26\t1\t4\t1\t1\t2\t2\t1\t28\t4\t1\t2\t3\t5\t17\n9icc:A:CR:339,157\t199\t313\t314\t349\t441\t446\t450\t463\t470\t486\t506\t542\t550\t583\t604\t631\t653\t658\t717\t817\t857\t882\t890\t896\t911\t916\t931\t957\t1005\t1051\t1113\t1115\t1140\t1164\t1210\t1211\t1215\t1250\t1257\t1261\t1332\t1340\t1367\t1371\t1374\t1418\t1423\t1517\t1719\t1724\t1740\t1742\t1771\t1837\t1855\t1859\t1899\t1909\t2042\t2064\t2085\t2103\t2141\t2144\t2210\t2224\t2245\t2254\t2256\t2287\t2314\t2334\t2352\t2422\t2428\t2468\t2492\t2535\t2537\t2556\t2591\t2641\t2665\t2674\t2719\t2731\t2788\t2828\t2832\t2879\t2894\t2898\t2916\t2935\t2972\t3029\t3109\t3126\t3132\t3133\t3213\t3243\t3299\t3335\t3416\t3465\t3478\t3485\t3491\t3510\t3655\t3669\t3884\t3896\t3933\t4043\t4048\t4082,2\t1\t1\t1\t1\t1\t1\t1\t4\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t2\t2\t4\t9\t1\t1\t1\t1\t1\t2\t1\t1\t2\t2\t1\t5\t6\t1\t5\t2\t1\t1\t1\t1\t1\t1\t1\t1\t9\t1\t1\t1\t1\t1\t2\t9\t1\t1\t2\t2\t1\t1\t1\t2\t7\t1\t3\t1\t1\t9\t1\t1\t4\t1\t3\t1\t1\t1\t3\t4\t2\t1\t1\t1\t3\t2\t1\t2\t1\t1\t1\t4\t1\t2\t1\t1\t1\t1\t7\t3\t3\t1\t2\t1\t1\t1\t1\t1\t1\t1\t6\t1\t1\t1\t1\t6\t1\t1\t1\t1\n1THA:B:HOH:137,83\t84\t592\t911\t1058\t1210\t1215\t1423\t1629\t1653\t1837\t2067\t2093\t2230\t2254\t2535\t3109\t3491\t4047,2\t1\t5\t5\t1\t6\t1\t9\t1\t1\t7\t1\t1\t2\t7\t9\t4\t7\t1\n")
        expected_df = pd.read_csv(expected_output)
        curr_df = pd.read_csv(ifp_file)
        self.assertTrue(expected_df.equals(curr_df))

        #
        # Atomic groups.
        expected_output = set([("6JWU:A:CA:806", 154),
                               ("1I2C:A:UPG:402", 614),
                               ("1KVT:A:NAD:340", 584),
                               ("3MJM:A:DOR:1410", 279),
                               ("3QL8:A:X01:300", 253),
                               ("6GST:A:GSH:218", 306),
                               ("1USN:A:ZN:258", 127),
                               ("4B3U:A:NWL:400", 286),
                               ("3I4Y:A:35C:1", 319),
                               ("9icc:A:CR:339", 122),
                               ("1THA:B:HOH:137", 71)])
        n_atm_grps = set([(agm.entry.to_string(), len(agm))
                          for agm in pli_obj.atm_grps_mngrs])
        self.assertTrue(expected_output == n_atm_grps)

        #
        # Interactions.
        expected_output = set([("6JWU:A:CA:806", 20),
                               ("1I2C:A:UPG:402", 258),
                               ("1KVT:A:NAD:340", 291),
                               ("3MJM:A:DOR:1410", 106),
                               ("3QL8:A:X01:300", 72),
                               ("6GST:A:GSH:218", 121),
                               ("1USN:A:ZN:258", 4),
                               ("4B3U:A:NWL:400", 83),
                               ("3I4Y:A:35C:1", 57),
                               ("9icc:A:CR:339", 57),
                               ("1THA:B:HOH:137", 0)])
        n_interactionss = set([(im.entry.to_string(), len(im))
                               for im in pli_obj.interactions_mngrs])
        self.assertTrue(expected_output == n_interactionss)
