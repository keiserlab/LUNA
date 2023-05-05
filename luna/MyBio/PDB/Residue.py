# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


###################################################################

# Modifications included by Alexandre Fassio (alexandrefassio@dcc.ufmg.br).

# Date: 02/19/2018.
# 1) The module __lt__ was overwritten to allow sorting in Python 3
# 2) Inherit inhouse modifications. Package: MyBio.

# Date: 02/21/2018.
# 1) Added function to check if a residue is a water molecule.
# 2) Added function to check if a residue is an hetero group.
# 3) Added function to check if a residue is an amino acid.
# 4) Added function to check if a residue is a nucleic acid.

# Date: 12/17/2018
# 1) Added property "_is_target" to control residues that will be targets for some processing.
# 2) Added function "is_target() to verify the status of the is_target variable.
# 3) Added function "set_as_target()" to allow the definition if a residue is a target or not.

# Date: 01/10/2019
# 1) Added function "get_class()" to get the compound class: amino acid, hetatm, water, nucleotide, or unknown.

# Data: 04/08/2019
# 1) New property idx added to Residue. This property stores the position of a residue in the list stored by its parent chain.

# Each line or block with modifications contain a MODBY tag.

###################################################################


"""Residue class, used by Structure objects."""

# My Stuff
import warnings
from Bio import BiopythonDeprecationWarning

# MODBY: Alexandre Fassio
# Inherit inhouse modifications. Package: MyBio.
from luna.MyBio.PDB.PDBExceptions import PDBConstructionException
from luna.MyBio.PDB.Entity import Entity, DisorderedEntityWrapper
from Bio.PDB.Polypeptide import is_aa

_atom_name_dict = {}
_atom_name_dict["N"] = 1
_atom_name_dict["CA"] = 2
_atom_name_dict["C"] = 3
_atom_name_dict["O"] = 4


# PDB ids for the following atoms were not found: Fr, Ra, Sc, Nb, Tc, Sn, Hf
METALS = ["LI", "NA", "K", "RB", "CS", "MG", "CA", "SR", "BA", "V", "CR",
          "MN", "MN3", "FE2", "FE", "CO", "3CO", "NI", "3NI", "CU1", "CU",
          "CU3", "ZN", "AL", "GA", "ZR", "MO", "4MO", "6MO", "RU", "RH",
          "RH3", "PD", "AG", "CD", "IN", "W", "RE", "OS", "OS4", "IR",
          "IR3", "PT", "PT4", "AU", "AU3", "HG", "TL", "PB", "BS3", "0BE",
          "4TI", "Y1", "YT3", "TA0"]


class Residue(Entity):
    """Represents a residue. A Residue object stores atoms."""

    def __init__(self, id, resname, segid, idx, at_line=None):
        self.level = "R"
        self.disordered = 0
        self.resname = resname
        self.segid = segid

        # MODBY: Alexandre Fassio
        # This property stores the position of a residue in the list stored by its chain parent.
        self.idx = idx

        # MODBY: Alexandre Fassio
        # This optional property stores the line number in the PDB file where this residue starts.
        self.at_line = at_line

        # MODBY: Alexandre Fassio
        # By default: no residue is a target for any calculations.
        self._is_target = False

        self.cluster_id = None

        Entity.__init__(self, id)

    # Special methods

    def __repr__(self):
        resname = self.get_resname()
        hetflag, resseq, icode = self.get_id()
        full_id = (resname, hetflag, resseq, icode)
        return "<Residue %s het=%s resseq=%s icode=%s>" % full_id

    # MODBY: Alexandre Fassio.
    # __lt__ method overwritten.
    def __lt__(self, r2):
        return self.idx < r2.idx

    # Private methods

    def _sort(self, a1, a2):
        """Sort the Atom objects.

        Atoms are sorted alphabetically according to their name,
        but N, CA, C, O always come first.

        Arguments:
         - a1, a2 - Atom objects

        """
        name1 = a1.get_name()
        name2 = a2.get_name()
        if name1 == name2:
            return(cmp(a1.get_altloc(), a2.get_altloc()))
        if name1 in _atom_name_dict:
            index1 = _atom_name_dict[name1]
        else:
            index1 = None
        if name2 in _atom_name_dict:
            index2 = _atom_name_dict[name2]
        else:
            index2 = None
        if index1 and index2:
            return cmp(index1, index2)
        if index1:
            return -1
        if index2:
            return 1

        return cmp(name1, name2)

    # Public methods

    @property
    def full_name(self):
        full_name = "%s/%s/%s" % self.get_full_id()[0:3]
        res_name = "%s/%d%s" % (self.resname, self.id[1], self.id[2].strip())
        full_name += "/%s" % res_name
        return full_name

    @property
    def metal_coordination(self):
        metal_coordination = {}
        for a in self.get_atoms():
            if a.has_metal_coordination():
                metal_coordination[a] = a.metal_coordination
        return metal_coordination

    def add(self, atom):
        """Add an Atom object.

        Checks for adding duplicate atoms, and raises a
        PDBConstructionException if so.
        """
        atom_id = atom.get_id()
        if self.has_id(atom_id):
            raise PDBConstructionException(
                "Atom %s defined twice in residue %s" % (atom_id, self))
        Entity.add(self, atom)

    def sort(self):
        self.child_list.sort(self._sort)

    def flag_disordered(self):
        """Set the disordered flag."""
        self.disordered = 1

    def is_disordered(self):
        """Return 1 if the residue contains disordered atoms."""
        return self.disordered

    # MODBY: Alexandre Fassio
    # Check if a residue is a water.
    def is_water(self):
        """Return 1 if the residue is a water molecule."""
        return self.get_id()[0] == "W"

    # MODBY: Alexandre Fassio
    # Check if a compound is a metal.
    def is_metal(self):
        """Return True if the residue is a metal."""
        return (self.get_id()[0].startswith("H_")
                and self.resname in METALS)

    # MODBY: Alexandre Fassio
    # Check if a compound is a hetero group.
    def is_hetatm(self):
        """Return True if the residue is an hetero group."""
        return (self.get_id()[0].startswith("H_")
                and not self.is_metal())

    # MODBY: Alexandre Fassio
    # Check if a residue is an amino acid.
    def is_residue(self):
        """Return True if the residue is an amino acid."""
        return self.get_id()[0] == " " and is_aa(self.resname)

    # MODBY: Alexandre Fassio
    # Check if a residue is a nucleotide.
    def is_nucleotide(self):
        """Return True if the residue is a nucleic acid."""
        return self.get_id()[0] == " " and not self.is_residue()

    # MODBY: Alexandre Fassio
    # Check if a residue is a target, i.e., if it will be
    # used for any calculations.
    def is_target(self):
        return self._is_target

    # MODBY: Alexandre Fassio
    # Get the compound class.
    def get_class(self):
        if self.is_water():
            return "Water"
        if self.is_hetatm():
            return "Hetatm"
        if self.is_residue():
            return "Residue"
        if self.is_nucleotide():
            return "Nucleotide"
        return "Unknown"

    def get_resname(self):
        return self.resname

    def get_unpacked_list(self):
        """Returns the list of all atoms, unpack DisorderedAtoms."""
        atom_list = self.get_list()
        undisordered_atom_list = []
        for atom in atom_list:
            if atom.is_disordered():
                undisordered_atom_list = (undisordered_atom_list + atom.disordered_get_list())
            else:
                undisordered_atom_list.append(atom)
        return undisordered_atom_list

    def get_segid(self):
        return self.segid

    def get_atoms(self):
        for a in self:
            yield a

    def get_atom(self):
        warnings.warn("`get_atom` has been deprecated and we intend to remove it"
                      " in a future release of Biopython. Please use `get_atoms` instead.",
                     BiopythonDeprecationWarning)
        for a in self:
            yield a

    def get_flag(self):
        if self.get_id()[0].startswith("H_"):
            return "H"
        return self.get_id()[0]

    # MODBY: Alexandre Fassio
    # Define if the residue is a target or not.
    def set_as_target(self, is_target=True):
        self._is_target = is_target

    def as_json(self):
        if self.is_water():
            comp_repr = "sphere"
        elif self.is_hetatm():
            if len(self.child_list) == 1 or len([atm for atm in self.child_list if atm.element != "H"]) == 1:
                comp_repr = "sphere"
            else:
                comp_repr = "stick"
        else:
            comp_repr = "stick"

        return {"chain": self.parent.id,
                "name": self.resname,
                "number": self.id[1],
                "icode": self.id[2],
                "repr": comp_repr}


class DisorderedResidue(DisorderedEntityWrapper):
    """DisorderedResidue is a wrapper around two or more Residue objects.

    It is used to represent point mutations (e.g. there is a Ser 60 and a Cys 60
    residue, each with 50 % occupancy).
    """

    def __init__(self, id):
        DisorderedEntityWrapper.__init__(self, id)

    def __repr__(self):
        resname = self.get_resname()
        hetflag, resseq, icode = self.get_id()
        full_id = (resname, hetflag, resseq, icode)
        return "<DisorderedResidue %s het=%s resseq=%i icode=%s>" % full_id

    def add(self, atom):
        residue = self.disordered_get()
        if not atom.is_disordered() == 2:
            # Atoms in disordered residues should have non-blank
            # altlocs, and are thus represented by DisorderedAtom objects.
            resname = residue.get_resname()
            het, resseq, icode = residue.get_id()
            # add atom anyway, if PDBParser ignores exception the atom will be part of the residue
            residue.add(atom)
            raise PDBConstructionException(
                "Blank altlocs in duplicate residue %s (%s, %i, %s)"
                % (resname, het, resseq, icode))
        residue.add(atom)

    def sort(self):
        """Sort the atoms in the child Residue objects."""
        for residue in self.disordered_get_list():
            residue.sort()

    def disordered_add(self, residue):
        """Add a residue object and use its resname as key.

        Arguments:
         - residue - Residue object

        """
        resname = residue.get_resname()
        # add chain parent to residue
        chain = self.get_parent()
        residue.set_parent(chain)
        assert(not self.disordered_has_id(resname))
        self[resname] = residue
        self.disordered_select(resname)
