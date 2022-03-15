import unittest

import sys
from os.path import dirname, abspath

sys.path.append(dirname(dirname(abspath(__file__))))

from luna.mol.entry import *
from luna.util.exceptions import InvalidEntry, IllegalArgumentError, MoleculeObjectError, MoleculeObjectTypeError, MoleculeNotFoundError


class EntryTest(unittest.TestCase):

    def test_init(self):
        # Missing obligatory parameter.
        self.assertRaises(TypeError, Entry, "3QL8")

        # Invalid chain id format.
        self.assertRaises(InvalidEntry, Entry, "3QL8", "AA")

        # Invalid chain id format.
        self.assertRaises(InvalidEntry, Entry, "3QL8", "+1")

        # Missing obligatory parameter when a compound information is provided.
        self.assertRaises(IllegalArgumentError, Entry, "3QL8", "A", "NAG")
        self.assertRaises(IllegalArgumentError, Entry, "3QL8", "A", comp_name="NAG")
        self.assertRaises(IllegalArgumentError, Entry, "3QL8", "A", comp_num=100)

        # Valid PDB and chain ID formats.
        self.assertEqual(Entry("3QL8", "A").to_string(), "3QL8:A")
        self.assertEqual(Entry("3QL8", "1").to_string(), "3QL8:1")
        self.assertEqual(Entry("3QL8", "A", None, None).to_string(), "3QL8:A")

        # Valid filenames.
        long_filename = "Q411SM6ZMPI6CBUKZSX10IA9M2923PUIHY88URUWC0R1ZOLQLJKXK9VATXRI3CM9W7RU409M19CG6S55TQ5L8L03HIKKL2HSEBTJ2VT7Y4G2JTUB558JOG2F3H9O524MMSYM0M3A347Q2JN77XG4MLAR3XI3EULR6WMVZ6U6NMASHPSUBE28ETOB0S1M0PTDX4KW45KLYC0JL5E77HKLAQK97NNOVB15B4SGVEVSZ90HZ3IMWZLF361AM0M0NSS"
        self.assertEqual(Entry(long_filename, "B").to_string(), "%s:B" % long_filename)
        self.assertEqual(Entry('file1.my_test=1a', "B").to_string(), "file1.my_test=1a:B")

        # Long filename.
        self.assertRaises(InvalidEntry, Entry, "4Q411SM6ZMPI6CBUKZSX10IA9M2923PUIHY88URUWC0R1ZOLQLJKXK9VATXRI3CM9W7RU409M19CG6S55TQ5L8L03HIKKL2HSEBTJ2VT7Y4G2JTUB558JOG2F3H9O524MMSYM0M3A347Q2JN77XG4MLAR3XI3EULR6WMVZ6U6NMASHPSUBE28ETOB0S1M0PTDX4KW45KLYC0JL5E77HKLAQK97NNOVB15B4SGVEVSZ90HZ3IMWZLF361AM0M0NSS", "B")

        # Valid compound entries.
        self.assertEqual(Entry("3QL8", "A", comp_name="X02", comp_num=100).to_string(), "3QL8:A:X02:100")
        self.assertEqual(Entry("3QL8", "A", "X02", 100).to_string(), "3QL8:A:X02:100")
        self.assertEqual(Entry("3QL8", "A", "X02", "100").to_string(), "3QL8:A:X02:100")
        self.assertEqual(Entry("3QL8", "A", 150, "100").to_string(), "3QL8:A:150:100")
        self.assertEqual(Entry("3QL8", "A", "X02", -100).to_string(), "3QL8:A:X02:-100")
        self.assertEqual(Entry("3QL8", "A", "X02", "+100").to_string(), "3QL8:A:X02:100")
        self.assertEqual(Entry("3NS9", "A", "NS9", "0").to_string(), "3NS9:A:NS9:0")

        # Valid compound entries with different compound name forms.
        self.assertEqual(Entry("3QL8", "A", "X", 533).to_string(), "3QL8:A:X:533")
        self.assertEqual(Entry("3QL8", "A", "XA", 533).to_string(), "3QL8:A:XA:533")
        self.assertEqual(Entry("3QL8", "A", "XAA", 533).to_string(), "3QL8:A:XAA:533")
        self.assertEqual(Entry("3QL8", "A", "Cl-", 533).to_string(), "3QL8:A:Cl-:533")
        self.assertEqual(Entry("3QL8", "A", "Na+", 533).to_string(), "3QL8:A:Na+:533")
        self.assertEqual(Entry("3QL8", "A", "NA+", 533).to_string(), "3QL8:A:NA+:533")
        self.assertEqual(Entry("3QL8", "A", "NA-", 533).to_string(), "3QL8:A:NA-:533")
        self.assertEqual(Entry("3QL8", "A", "N++", 533).to_string(), "3QL8:A:N++:533")
        self.assertEqual(Entry("3QL8", "A", "O--", 533).to_string(), "3QL8:A:O--:533")

        # Valid compound entries with Icode
        self.assertEqual(Entry("3QL8", "A", "X01", "100", "A").to_string(), "3QL8:A:X01:100A")

        # Invalid compound entries.
        self.assertRaises(InvalidEntry, Entry, "3QL8", "A", "X043", 10)
        self.assertRaises(IllegalArgumentError, Entry, "3QL8", "A", "X01", "XAA")
        self.assertRaises(IllegalArgumentError, Entry, "3QL8", "A", "X01", "X")
        self.assertRaises(IllegalArgumentError, Entry, "3QL8", "A", None, "X")
        self.assertRaises(IllegalArgumentError, Entry, "3QL8", "A", "X01", None)

        # Invalid compound entries with Icode
        self.assertRaises(IllegalArgumentError, Entry, "3QL8", "A", "X01", 100, "1")
        self.assertRaises(IllegalArgumentError, Entry, "3QL8", "A", "X01", 100, 1)

    def test_property_update(self):
        entry = Entry("3QL8", "A")

        with self.assertRaises(AttributeError):
            entry.pdb_id = "3QQK"

        with self.assertRaises(AttributeError):
            entry.chain_id = "C"

        with self.assertRaises(AttributeError):
            entry.comp_name = "X05"

        with self.assertRaises(AttributeError):
            entry.comp_num = 100

        with self.assertRaises(AttributeError):
            entry.comp_icode = "B"

        entry.is_hetatm = False
        self.assertEqual(entry.is_hetatm, False)

        entry.sep = "_"
        self.assertEqual(entry.to_string(), "3QL8_A")

    def test_from_string(self):
        # Missing obligatory parameter.
        self.assertRaises(IllegalArgumentError, Entry.from_string, "3QL8")

        # Invalid chain id format.
        self.assertRaises(InvalidEntry, Entry, "3QL8", "AA")

        # Invalid chain id format.
        self.assertRaises(InvalidEntry, Entry.from_string, "3QL8:+1")

        # Missing obligatory parameter when a compound information is provided.
        self.assertRaises(IllegalArgumentError, Entry.from_string, "3QL8:A:NAG")
        self.assertRaises(IllegalArgumentError, Entry.from_string, "3QL8:A:100")

        # Valid PDB and chain ID formats.
        self.assertEqual(Entry.from_string("3QL8:A").to_string(), "3QL8:A")
        self.assertEqual(Entry.from_string("3QL8:1").to_string(), "3QL8:1")
        self.assertEqual(Entry.from_string("3QL8-A", sep="-").to_string(), "3QL8-A")
        self.assertEqual(Entry.from_string("3QL8_A", sep="_").to_string(), "3QL8_A")

        # Valid filenames.
        long_filename = "Q411SM6ZMPI6CBUKZSX10IA9M2923PUIHY88URUWC0R1ZOLQLJKXK9VATXRI3CM9W7RU409M19CG6S55TQ5L8L03HIKKL2HSEBTJ2VT7Y4G2JTUB558JOG2F3H9O524MMSYM0M3A347Q2JN77XG4MLAR3XI3EULR6WMVZ6U6NMASHPSUBE28ETOB0S1M0PTDX4KW45KLYC0JL5E77HKLAQK97NNOVB15B4SGVEVSZ90HZ3IMWZLF361AM0M0NSS"
        self.assertEqual(Entry.from_string("%s:B" % long_filename).to_string(), "%s:B" % long_filename)
        self.assertEqual(Entry.from_string('file1.my_test=1a:B').to_string(), "file1.my_test=1a:B")

        # Long filename.
        self.assertRaises(InvalidEntry, Entry.from_string, "4Q411SM6ZMPI6CBUKZSX10IA9M2923PUIHY88URUWC0R1ZOLQLJKXK9VATXRI3CM9W7RU409M19CG6S55TQ5L8L03HIKKL2HSEBTJ2VT7Y4G2JTUB558JOG2F3H9O524MMSYM0M3A347Q2JN77XG4MLAR3XI3EULR6WMVZ6U6NMASHPSUBE28ETOB0S1M0PTDX4KW45KLYC0JL5E77HKLAQK97NNOVB15B4SGVEVSZ90HZ3IMWZLF361AM0M0NSS:B")

        # Valid compound entries.
        self.assertEqual(Entry.from_string("3QL8:A:X02:100").to_string(), "3QL8:A:X02:100")
        self.assertEqual(Entry.from_string("3QL8:A:150:100").to_string(), "3QL8:A:150:100")
        self.assertEqual(Entry.from_string("3QL8:A:X02:-100").to_string(), "3QL8:A:X02:-100")
        self.assertEqual(Entry.from_string("3QL8:A:X02:+100").to_string(), "3QL8:A:X02:100")

        # Valid compound entries with different compound name forms.
        self.assertEqual(Entry.from_string("3QL8:A:X:533").to_string(), "3QL8:A:X:533")
        self.assertEqual(Entry.from_string("3QL8:A:XA:533").to_string(), "3QL8:A:XA:533")
        self.assertEqual(Entry.from_string("3QL8:A:XAA:533").to_string(), "3QL8:A:XAA:533")
        self.assertEqual(Entry.from_string("3QL8:A:Cl-:533").to_string(), "3QL8:A:Cl-:533")
        self.assertEqual(Entry.from_string("3QL8:A:Na+:533").to_string(), "3QL8:A:Na+:533")
        self.assertEqual(Entry.from_string("3QL8:A:NA+:533").to_string(), "3QL8:A:NA+:533")
        self.assertEqual(Entry.from_string("3QL8:A:NA-:533").to_string(), "3QL8:A:NA-:533")
        self.assertEqual(Entry.from_string("3QL8:A:N++:533").to_string(), "3QL8:A:N++:533")
        self.assertEqual(Entry.from_string("3QL8:A:O--:533").to_string(), "3QL8:A:O--:533")
        self.assertEqual(Entry.from_string("3QL8:A:X02:9999").to_string(), "3QL8:A:X02:9999")

        # Valid compound entries with Icode
        self.assertEqual(Entry.from_string("3QL8:A:X01:100A").to_string(), "3QL8:A:X01:100A")

        # Invalid compound entries.
        self.assertRaises(InvalidEntry, Entry.from_string, "3QL8:A:X043:10")
        self.assertRaises(InvalidEntry, Entry.from_string, "3QL8:A:X02:99999")
        self.assertRaises(IllegalArgumentError, Entry.from_string, "3QL8:A:X04:10:A")
        self.assertRaises(IllegalArgumentError, Entry.from_string, "3QL8:A:X04:10:A:10")
        self.assertRaises(IllegalArgumentError, Entry.from_string, "3QL8:A:X01:XAA")
        self.assertRaises(IllegalArgumentError, Entry.from_string, "3QL8:A:X01:X")
        self.assertRaises(IllegalArgumentError, Entry.from_string, "3QL8:A:None:X")
        self.assertRaises(IllegalArgumentError, Entry.from_string, "3QL8:A:X01")

    def test_get_biopython_key(self):
        entry = Entry("3QL8", "A")
        self.assertEqual(entry.get_biopython_key(), "A")

        entry = Entry("3QL8", "A", "X02", 104)
        self.assertEqual(entry.get_biopython_key(), ("H_X02", 104, " "))

        entry = Entry("3QL8", "A", "X02", 104, "A")
        self.assertEqual(entry.get_biopython_key(), ("H_X02", 104, "A"))

        entry = Entry("3QL8", "A", "X02", 104, " ")
        self.assertEqual(entry.get_biopython_key(), ("H_X02", 104, " "))

        entry = Entry("3QL8", "A", "X02", 104, None)
        self.assertEqual(entry.get_biopython_key(), ("H_X02", 104, " "))

        entry = Entry("3QL8", "A", "X02", 104, is_hetatm=False)
        self.assertEqual(entry.get_biopython_key(), (" ", 104, " "))

        entry = Entry("3QL8", "A", "HOH", 104, is_hetatm=False)
        self.assertEqual(entry.get_biopython_key(), ("W", 104, " "))

        entry = Entry("3QL8", "A", "HOH", 104, is_hetatm=True)
        self.assertEqual(entry.get_biopython_key(), ("W", 104, " "))

        entry = Entry("3QL8", "A", "WAT", 104, is_hetatm=False)
        self.assertEqual(entry.get_biopython_key(), ("W", 104, " "))


if __name__ == '__main__':
    unittest.main()
