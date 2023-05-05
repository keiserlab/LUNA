from os import path
from collections import defaultdict
from pdbecif.mmcif_io import CifFileReader

from luna.util.file import parse_json_file

import logging
logger = logging.getLogger()


class MoleculesData(dict):
    """Generic class to define precomputed molecular data such as
    default bonds, atomic invariants, and physicochemical
    properties.

    Parameters
    ----------
    data : dict, optional
        A dict containing molecules' data.
    """

    def __init__(self, data=None):
        data = data or {}

        super().__init__(data)

    @classmethod
    def from_json_file(cls, json_file):

        data = parse_json_file(json_file)

        return cls(data)


class DefaultResidueData(MoleculesData):
    """Default parameters for calculating interactions in LUNA."""

    def __init__(self, **kwargs):

        file_path = path.abspath(path.join(path.realpath(__file__),
                                           '../../data'))
        filename = "precomputed_residue_data.json"
        default_file = f"{file_path}/{filename}"

        data = MoleculesData.from_json_file(default_file)

        super().__init__(data)
