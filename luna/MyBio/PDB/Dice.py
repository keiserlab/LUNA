# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

###################################################################

# Modifications included by Alexandre Fassio (alexandrefassio@dcc.ufmg.br).
# Date: 19/02/2018.

# 1) Inherit inhouse modifications. Package: MyBio.

# Each line or block with modifications contain a MODBY tag.

###################################################################

"""Code for chopping up (dicing) a structure.

This module is used internally by the Bio.PDB.extract() function.
"""

import re
import warnings

# MODBY: Alexandre Fassio
# Inherit inhouse modifications. Package: MyBio.
from luna.MyBio.PDB.PDBIO import PDBIO

from Bio import BiopythonWarning

_hydrogen = re.compile("[123 ]*H.*")


class ChainSelector(object):
    """Only accepts residues with right chainid, between start and end.

    Remove hydrogens, waters and ligands. Only use model 0 by default.
    """

    def __init__(self, chain_id, start, end, model_id=0):
        self.chain_id = chain_id
        self.start = start
        self.end = end
        self.model_id = model_id

    def accept_model(self, model):
        # model - only keep model 0
        if model.get_id() == self.model_id:
            return 1
        return 0

    def accept_chain(self, chain):
        if chain.get_id() == self.chain_id:
            return 1
        return 0

    def accept_residue(self, residue):
        # residue - between start and end
        hetatm_flag, resseq, icode = residue.get_id()
        if hetatm_flag != " ":
            # skip HETATMS
            return 0
        if icode != " ":
            warnings.warn("WARNING: Icode %s at position %s"
                          % (icode, resseq), BiopythonWarning)
        if self.start <= resseq <= self.end:
            return 1
        return 0

    def accept_atom(self, atom):
        # atoms - get rid of hydrogens
        name = atom.get_id()
        if _hydrogen.match(name):
            return 0
        else:
            return 1


def extract(structure, chain_id, start, end, filename):
    """Write out selected portion to filename."""
    sel = ChainSelector(chain_id, start, end)
    io = PDBIO()
    io.set_structure(structure)
    io.save(filename, sel)


if __name__ == "__main__":

    # MODBY: Alexandre Fassio
    # Inherit inhouse modifications. Package: MyBio.
    from luna.MyBio.PDB.PDBParser import PDBParser

    import sys

    p = PDBParser()
    s = p.get_structure("scr", sys.argv[1])

    extract(s, " ", 1, 100, "out.pdb")
