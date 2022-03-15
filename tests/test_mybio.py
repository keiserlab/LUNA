import unittest

import sys
from os.path import dirname, abspath

sys.path.append(dirname(dirname(abspath(__file__))))

from luna.MyBio.PDB.Model import Model
from luna.MyBio.PDB.Chain import Chain
from luna.MyBio.extractor import Extractor


class ExtractorTest(unittest.TestCase):

    def test_init(self):
        # Wrong type.
        self.assertRaises(TypeError, Extractor, "3QL8")
        self.assertRaises(TypeError, Extractor, 1)
        self.assertRaises(TypeError, Extractor, int)
        self.assertRaises(TypeError, Extractor, self)
        self.assertRaises(TypeError, Extractor, Chain)
        self.assertRaises(TypeError, Extractor, None)

        # Wrong initialization.
        with self.assertRaises(TypeError):
            Extractor(Chain())
        # Wrong initialization.
        with self.assertRaises(TypeError):
            Extractor(Model())

        # Extractor correctly initialized.
        model = Model(0)
        extractor = Extractor(model)
        self.assertEqual(extractor.entity.id, model.id)
        self.assertNotEqual(extractor.entity.id, "A")
        self.assertNotEqual(extractor.entity.id, 1)

        extractor.extract_chains(["A"], "test.pdb")

        # Extractor correctly initialized.
        chain = Chain("A")
        extractor = Extractor(chain)
        self.assertEqual(extractor.entity.id, chain.id)
        self.assertNotEqual(extractor.entity.id, "B")
        self.assertNotEqual(extractor.entity.id, 1)


if __name__ == '__main__':
    unittest.main()
