import unittest

import sys
from os.path import dirname, abspath

sys.path.append(dirname(dirname(abspath(__file__))))

from luna.util.progress_tracker import *


class ProgressResultTest(unittest.TestCase):

    def test_init(self):
        # Wrong type.
        self.assertRaises(TypeError, ProgressResult, "3QL8")
        self.assertRaises(TypeError, ProgressResult, 1)
        self.assertRaises(TypeError, ProgressResult, ProgressData(1, 20))

    def test_attribute_access(self):
        # Access 'inputs' when there is invalid objects in 'results'.
        with self.assertRaises(TypeError):
            pr = ProgressResult(["a", "b", "c"])
            pr.inputs

        # Access 'outputs' when there is invalid objects in 'results'.
        with self.assertRaises(TypeError):
            pr = ProgressResult(["a", "b", "c"])
            pr.outputs

        # Access 'errors' when there is invalid objects in 'results'.
        with self.assertRaises(TypeError):
            pr = ProgressResult(["a", "b", "c"])
            pr.errors

    def test_append(self):

        # Appending invalid data to 'results'
        with self.assertRaises(TypeError):
            pr = ProgressResult()
            pr.append(1)

        # Appending invalid data to 'results'
        with self.assertRaises(TypeError):
            pr = ProgressResult()
            pr.append("string")

        # Appending invalid data to 'results'
        with self.assertRaises(TypeError):
            pr = ProgressResult()
            pr.append(1)

        # Appending invalid data to 'results'
        with self.assertRaises(TypeError):
            pr = ProgressResult([ProgressData(1, 10), ProgressData(2, 11)])
            pr.append(1)

        # Check valid appending.
        pr = ProgressResult()
        pr.append(ProgressData(1, 20))
        self.assertEqual(pr.inputs, [1])

        # Check valid appending.
        pr = ProgressResult([ProgressData(1, 10), ProgressData(2, 11)])
        pr.append(ProgressData(3, 12))
        self.assertEqual(pr.inputs, [1, 2, 3])

        # Check outputs for valid appending.
        pr = ProgressResult()
        pr.append(ProgressData(3, 12, output_data="a"))
        self.assertEqual(pr.outputs, [(3, "a")])
        self.assertNotEqual(pr.outputs, [(3, "b")])

        # Check outputs for valid initialization.
        pr = ProgressResult([ProgressData(3, 12, output_data="a")])
        self.assertEqual(pr.outputs, [(3, "a")])
        self.assertNotEqual(pr.outputs, [(3, "b")])
        self.assertNotEqual(pr.outputs, [(1, "a")])


if __name__ == '__main__':
    unittest.main()
