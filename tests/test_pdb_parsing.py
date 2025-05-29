import unittest
import tempfile
from pathlib import Path
from io import StringIO

from luna.pdb.parser.base import PDBParser
from luna.pdb.parser.helpers import load_from_string, load_from_file
from luna.pdb.io.helpers import save_to_file

from luna.pdb.core.structure import Structure
from luna.pdb.core.model import Model
from luna.pdb.core.chain import Chain
from luna.pdb.core.residue import Residue
from luna.pdb.core.atom import Atom


PDB_STRING = """
HEADER    TEST PDB
ATOM      1  N   MET A   1      11.104  13.207  -5.219  1.00 13.79           N
ATOM      2  CA  MET A   1      10.749  11.821  -5.624  1.00 14.13           C
ATOM      3  C   MET A   1       9.258  11.611  -5.967  1.00 13.66           C
ATOM      4  O   MET A   1       8.646  10.571  -6.307  1.00 13.46           O
ATOM      5  CB  MET A   1      11.504  11.289  -6.826  1.00 15.48           C
TER
END
"""


class TestParsingFromString(unittest.TestCase):

    def setUp(self):
        self.structure = load_from_string("TEST", PDB_STRING)

    def test_structure_type(self):
        self.assertIsInstance(self.structure, Structure)

    def test_model_type(self):
        for model in self.structure:
            self.assertIsInstance(model, Model)

    def test_chain_type(self):
        for model in self.structure:
            for chain in model:
                self.assertIsInstance(chain, Chain)

    def test_residue_type(self):
        for chain in self.structure.get_chains():
            for residue in chain:
                self.assertIsInstance(residue, Residue)

    def test_atom_type(self):
        for atom in self.structure.get_atoms():
            self.assertIsInstance(atom, Atom)

    def test_hierarchy_name(self):
        first_atom = next(self.structure.get_atoms())
        self.assertIsInstance(first_atom.hierarchy_name, str)
        self.assertTrue("TEST" in first_atom.hierarchy_name)
        self.assertTrue(first_atom.hierarchy_name == 'TEST/0/A/MET/1/N')

    def test_hierarchy_traversal(self):
        first_atom = next(self.structure.get_atoms())
        
        entity = first_atom.get_parent_by_level('R')
        self.assertIsInstance(entity, Residue)

        entity = first_atom.get_parent_by_level('C')
        self.assertIsInstance(entity, Chain)

        entity = first_atom.get_parent_by_level('M')
        self.assertIsInstance(entity, Model)

        entity = first_atom.get_parent_by_level('S')
        self.assertIsInstance(entity, Structure)

        res = first_atom.parent
        entity = res.get_parent_by_level('M')
        self.assertIsInstance(entity, Model)

        with self.assertRaises(ValueError):
            res.get_parent_by_level('A') 


class TestLoadSavePDB(unittest.TestCase):

    def test_load_save_reload_consistency(self):
        # 1. Load from string
        struct1 = load_from_string("TEST", PDB_STRING)
        self.assertIsInstance(struct1, Structure)

        # 2. Save to temp file
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "test.pdb"
            save_to_file(struct1, out_path)

            # 3. Load from file
            struct2 = load_from_file("TEST", out_path)
            self.assertIsInstance(struct2, Structure)

            # 4. Compare atoms
            atoms1 = list(struct1.get_atoms())
            atoms2 = list(struct2.get_atoms())

            self.assertEqual(len(atoms1), len(atoms2))
            for a1, a2 in zip(atoms1, atoms2):
                self.assertEqual(a1.get_name(), a2.get_name())
                self.assertTrue((a1.coord == a2.coord).all())


if __name__ == '__main__':
    unittest.main()