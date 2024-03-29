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



"""Superimpose two structures."""

from __future__ import print_function

import numpy

from Bio.SVDSuperimposer import SVDSuperimposer

# MODBY: Alexandre Fassio
# Inherit inhouse modifications. Package: MyBio.
from luna.MyBio.PDB.PDBExceptions import PDBException


class Superimposer(object):
    """Rotate/translate one set of atoms on top of another,
    thereby minimizing the RMSD.
    """

    def __init__(self):
        self.rotran = None
        self.rms = None

    def set_atoms(self, fixed, moving):
        """Put (translate/rotate) the atoms in fixed on the atoms in
        moving, in such a way that the RMSD is minimized.

        @param fixed: list of (fixed) atoms
        @param moving: list of (moving) atoms
        @type fixed,moving: [L{Atom}, L{Atom},...]
        """
        if not (len(fixed) == len(moving)):
            raise PDBException("Fixed and moving atom lists differ in size")
        l = len(fixed)
        fixed_coord = numpy.zeros((l, 3))
        moving_coord = numpy.zeros((l, 3))
        for i in range(0, len(fixed)):
            fixed_coord[i] = fixed[i].get_coord()
            moving_coord[i] = moving[i].get_coord()
        sup = SVDSuperimposer()
        sup.set(fixed_coord, moving_coord)
        sup.run()
        self.rms = sup.get_rms()
        self.rotran = sup.get_rotran()

    def apply(self, atom_list):
        """Rotate/translate a list of atoms."""
        if self.rotran is None:
            raise PDBException("No transformation has been calculated yet")
        rot, tran = self.rotran
        rot = rot.astype('f')
        tran = tran.astype('f')
        for atom in atom_list:
            atom.transform(rot, tran)
