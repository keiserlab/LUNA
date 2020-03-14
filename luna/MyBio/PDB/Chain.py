# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

###################################################################

# Modifications included by Alexandre Fassio (alexandrefassio@dcc.ufmg.br).
# Date: 19/02/2018.

# 1) The module __lt__ was overwritten to allow sorting in Python 3
# 2) Inherit inhouse modifications. Package: MyBio.

# Date: 02/22/2019
# 1) Added property "_is_target" to control chains that will be targets for some processing.
# 2) Added function "is_target() to verify the status of the is_target variable.
# 3) Added function "set_as_target()" to allow the definition if a chain is a target or not.

# Each line or block with modifications contain a MODBY tag.

###################################################################


"""Chain class, used in Structure objects."""

# MODBY: Alexandre Fassio
# Inherit inhouse modifications. Package: MyBio.
from luna.MyBio.PDB.Entity import Entity


class Chain(Entity):
    def __init__(self, id):
        self.level = "C"

        # MODBY: Alexandre Fassio
        # By default: no chain is a target for any calculations.
        self._is_target = False

        Entity.__init__(self, id)

    # Special methods

    # MODBY: Alexandre Fassio.
    # __lt__ method overwritten.
    def __lt__(self, c2):
        return self.id < c2.id

    # Private methods

    def _sort(self, r1, r2):
        """Sort function for residues in a chain (PRIVATE).

        Residues are first sorted according to their hetatm records.
        Protein and nucleic acid residues first, hetatm residues next,
        and waters last. Within each group, the residues are sorted according
        to their resseq's (sequence identifiers). Finally, residues with the
        same resseq's are sorted according to icode.

        Arguments:

        - r1, r2 - Residue objects

        """
        hetflag1, resseq1, icode1 = r1.id
        hetflag2, resseq2, icode2 = r2.id
        if hetflag1 != hetflag2:
            return cmp(hetflag1[0], hetflag2[0])
        elif resseq1 != resseq2:
            return cmp(resseq1, resseq2)
        return cmp(icode1, icode2)

    def _translate_id(self, id):
        """Translate sequence identifer to tuple form (PRIVATE).

        A residue id is normally a tuple (hetero flag, sequence identifier,
        insertion code). Since for most residues the hetero flag and the
        insertion code are blank (i.e. " "), you can just use the sequence
        identifier to index a residue in a chain. The _translate_id method
        translates the sequence identifier to the (" ", sequence identifier,
        " ") tuple.

        Arguments:

        - id - int, residue resseq

        """
        if isinstance(id, int):
            id = (' ', id, ' ')
        return id

    # Special methods

    def __getitem__(self, id):
        """Return the residue with given id.

        The id of a residue is (hetero flag, sequence identifier, insertion code).
        If id is an int, it is translated to (" ", id, " ") by the _translate_id
        method.

        Arguments:

        - id - (string, int, string) or int

        """
        id = self._translate_id(id)
        return Entity.__getitem__(self, id)

    def __contains__(self, id):
        """True if a residue with given id is present in this chain.

        Arguments:

        - id - (string, int, string) or int

        """
        id = self._translate_id(id)
        return Entity.__contains__(self, id)

    def __delitem__(self, id):
        """Delete item.

        Arguments:

        - id - (string, int, string) or int

        """
        id = self._translate_id(id)
        return Entity.__delitem__(self, id)

    def __repr__(self):
        return "<Chain id=%s>" % self.get_id()

    # Public methods

    def get_unpacked_list(self):
        """Return a list of undisordered residues.

        Some Residue objects hide several disordered residues
        (DisorderedResidue objects). This method unpacks them,
        ie. it returns a list of simple Residue objects.
        """
        unpacked_list = []
        for residue in self.get_list():
            if residue.is_disordered() == 2:
                for dresidue in residue.disordered_get_list():
                    unpacked_list.append(dresidue)
            else:
                unpacked_list.append(residue)
        return unpacked_list

    def has_id(self, id):
        """Return 1 if a residue with given id is present.

        The id of a residue is (hetero flag, sequence identifier, insertion code).

        If id is an int, it is translated to (" ", id, " ") by the _translate_id
        method.

        Arguments:

        - id - (string, int, string) or int

        """
        id = self._translate_id(id)
        return Entity.has_id(self, id)

    # MODBY: Alexandre Fassio
    # Check if a chain is a target, i.e., if it will be used for any calculations.
    def is_target(self):
        return self._is_target

    # Public

    def get_residues(self):
        for r in self:
            yield r

    def get_atoms(self):
        for r in self.get_residues():
            for a in r:
                yield a

    # MODBY: Alexandre Fassio
    # Define if all the residues in the chain are a target or not.
    def set_as_target(self, is_target=True):
        self._is_target = is_target

        for r in self.get_residues():
            r.set_as_target(is_target)
