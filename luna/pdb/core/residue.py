from Bio.PDB.Residue import Residue as BioResidue
from Bio.PDB.Residue import DisorderedResidue as BioDisorderedResidue
from Bio.PDB.Polypeptide import is_aa

from luna.pdb.core.entity import Entity


# List of metal residue names (based on PDB naming conventions)
METALS = ["LI", "NA", "K", "RB", "CS", "MG", "CA", "SR", "BA", "V", "CR",
          "MN", "MN3", "FE2", "FE", "CO", "3CO", "NI", "3NI", "CU1", "CU",
          "CU3", "ZN", "AL", "GA", "ZR", "MO", "4MO", "6MO", "RU", "RH",
          "RH3", "PD", "AG", "CD", "IN", "W", "RE", "OS", "OS4", "IR",
          "IR3", "PT", "PT4", "AU", "AU3", "HG", "TL", "PB", "BS3", "0BE",
          "4TI", "Y1", "YT3", "TA0"]


class Residue(BioResidue, Entity):
    """
    Custom Residue subclass with support for indexing, classification,
    hierarchical naming, and processing flags.

    Parameters
    ----------
    id : tuple
        Biopython residue ID (hetflag, resseq, icode).
    resname : str
        Residue name (e.g., "ALA", "HOH").
    segid : str
        Segment ID.
    idx : int, optional
        Index in parent chain.
    at_line : int, optional
        Line number in the original PDB.
    """

    def __init__(self, id, resname, segid, idx=None, at_line=None):
        super().__init__(id, resname, segid)

        # Store position within parent chain
        self.idx = idx
        # Optional PDB line number
        self.at_line = at_line 

        # Used by FTMapParser, optional.
        self.cluster_id = None
        
        # By default: no residue is a reference for processing.
        self._is_reference = False
    
    def __repr__(self):
        """Return the residue full id."""
        resname = self.get_resname()
        hetflag, resseq, icode = self.get_id()
        full_id = (resname, hetflag, resseq, icode)
        return "<Residue (LUNA subclass) %s het=%s resseq=%s icode=%s>" % full_id

    def __lt__(self, other):
        return self.idx < other.idx

    @property
    def metal_coordination(self):
        """
        Return atoms in this residue that participate in metal coordination.

        Returns
        -------
        dict
            Atom → set of coordination partners.
        """
        return {
            atom: atom.metal_coordination
            for atom in self.get_atoms()
            if getattr(atom, "has_metal_coordination", lambda: False)()
        }

    def is_reference(self):
        """
        Return True if this entity is marked as the reference for interaction analysis.

        This flag indicates that the user selected this entity
        as the *main subject* of interaction evaluation — interactions will be computed
        between this entity and others nearby (e.g., ligands, water, ions).

        Returns
        -------
        bool
            True if this entity is marked for interaction reference.
        """
        return self._is_reference

    def set_as_reference(self, is_reference=True):
        """
        Mark this entity as the reference for interaction analysis.

        This method sets a flag indicating that this residue should be
        treated as the main subject for calculating interactions 
        with nearby molecules.

        Parameters
        ----------
        is_reference : bool, optional
            If True, mark the entity as the interaction reference.
            If False, unmark it. Defaults to True.

        See Also
        --------
        is_reference : Check whether this entity is currently marked.
        """
        self._is_reference = is_reference

    def is_water(self):
        """Return True if the residue is a water molecule."""
        return self.get_id()[0] == "W"

    def is_metal(self):
        """Return True if the residue is a metal."""
        return (self.get_id()[0].startswith("H_")
                and self.resname in METALS)

    def is_hetatm(self):
        """Return True if the residue is an hetero group."""
        return (self.get_id()[0].startswith("H_")
                and not self.is_metal())

    def is_residue(self):
        """Return True if the residue is an amino acid."""
        return self.get_id()[0] == " " and is_aa(self.resname)

    def is_nucleotide(self):
        """Return True if the residue is a nucleic acid."""
        return self.get_id()[0] == " " and not self.is_residue()

    def is_monoatom(self):
        """Return True if the residue contains only one atom (or one non-hydrogen)."""
        atoms = list(self.get_atoms())
        return len(atoms) == 1 or len([a for a in atoms if getattr(a, "element", None) != "H"]) == 1

    def get_class(self):
        """Get the compound class."""
        if self.is_water():
            return "Water"
        if self.is_hetatm():
            return "Hetatm"
        if self.is_residue():
            return "Residue"
        if self.is_nucleotide():
            return "Nucleotide"
        return "Unknown"

    def get_flag(self):
        """
        Return a shorthand flag for residue classification.

        Returns
        -------
        str
            "H" for HETATM or metal, otherwise the original hetflag.
        """
        return "H" if self.get_id()[0].startswith("H_") else self.get_id()[0]

    def as_json(self):
        """
        Serialize this residue to a compact dictionary representation.

        Fields:
            chain, name, number, icode, repr

        'repr' will be 'sphere' for waters or monoatomic hetatms, 'stick' otherwise.
        """
        comp_repr = (
            "sphere" if self.is_water() or (self.is_hetatm() and self.is_monoatom()) else "stick"
        )
        return {
            "chain": self.get_parent().id,
            "name": self.resname,
            "number": self.id[1],
            "icode": self.id[2],
            "repr": comp_repr
        }


class DisorderedResidue(BioDisorderedResidue):
    """
    LUNA's DisorderedResidue placeholder.

    Currently identical to Biopython's, but allows future extensions and consistent imports.
    """
    def __init__(self, id):
        super().__init__(id)

    def __repr__(self):
        """Return disordered residue full identifier."""
        if self.child_dict:
            resname = self.get_resname()
            hetflag, resseq, icode = self.get_id()
            full_id = (resname, hetflag, resseq, icode)
            return "<DisorderedResidue (LUNA subclass) %s het=%s resseq=%i icode=%s>" % full_id
        else:
            return "<Empty DisorderedResidue (LUNA subclass)>"