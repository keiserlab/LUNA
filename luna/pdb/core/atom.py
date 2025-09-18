from Bio.PDB.Atom import Atom as BioAtom
from Bio.PDB.Atom import DisorderedAtom as BioDisorderedAtom


class Atom(BioAtom):
    """
    A custom Atom class that supports hierarchical traversal and
    additional coordination metadata.

    Inherits from Biopython's Atom and adds:
    - `metal_coordination` set
    - `hierarchy_name` property
    - `get_parent_by_level()` for SMCRA traversal
    - Custom sorting via __lt__

    Parameters
    ----------
    name : str
        Atom name (e.g., "CA").
    coord : array_like
        Cartesian coordinates.
    bfactor : float
        B-factor value.
    occupancy : float
        Occupancy.
    altloc : str
        Alternate location indicator.
    fullname : str
        Full atom name with padding (e.g., " CA ").
    serial_number : int
        Atom serial number.
    element : str, optional
        Atomic element.
    """
    def __init__(self, name, coord, bfactor, occupancy, altloc,
                 fullname, serial_number, element=None):
        
        super().__init__(name, coord, bfactor, occupancy, altloc,
                         fullname, serial_number, element)
        self.metal_coordination = set()

    def __repr__(self):
        return f"<Atom (LUNA subclass) {self.get_id()}>"

    def __lt__(self, other):
        """
        Atom sorting based on structure hierarchy and residue index.

        Ensures original PDB order is preserved.
        """
        id1 = self.get_full_id()[0:3] + (self.get_parent().idx,) + self.get_full_id()[4:]
        id2 = other.get_full_id()[0:3] + (other.get_parent().idx,) + other.get_full_id()[4:]
        return id1 < id2

    @property
    def hierarchy_name(self):
        """
        Return full hierarchical name of the atom.
    
        Format: structure/model/chain/resname/resid/atom[-altloc]
        """
        full_id = self.get_full_id()
        parent = self.get_parent()
    
        parent_path = '/'.join(str(x) for x in full_id[0:3])
        resname = f"{parent.resname}/{parent.id[1]}{parent.id[2].strip()}"
        atom_name = self.name + (f"-{self.altloc}" if self.altloc != " " else "")
    
        return f"{parent_path}/{resname}/{atom_name}"

    def has_metal_coordination(self) -> bool:
        """Return True if the atom is metal-coordinated."""
        return bool(self.metal_coordination)

    # As the Atom class only partly implements the Entity interface, 
    # this function was included in both classes.
    def get_parent_by_level(self, level):
        """
        Return the parent entity at the specified SMCRA level.

        Parameters
        ----------
        level : str
            One of "A", "R", "C", "M", "S".

        Returns
        -------
        Entity
            The parent entity at the requested level.
        """
        if self.level == level:
            return self
        else:
            return self.parent.get_parent_by_level(level)


class DisorderedAtom(BioDisorderedAtom):
    """
    LUNA's DisorderedAtom placeholder.

    Currently identical to Biopython's, but allows future extensions and consistent imports.
    """
    def __init__(self, id):
        super().__init__(id)

    def __repr__(self):
        """Return disordered atom identifier."""
        if self.child_dict:
            return f"<DisorderedAtom (LUNA subclass) {self.get_id()}>"
        else:
            return f"<Empty DisorderedAtom (LUNA subclass) {self.get_id()}>"
