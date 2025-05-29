
from Bio.PDB.Entity import Entity as BioEntity


class Entity(BioEntity):
    """
    Base class for LUNA's Structure, Model, Chain, and Residue classes.

    Wraps child management, identifier logic, and hierarchical traversal.
    """

    @property
    def hierarchy_name(self):
        """
        Return a full hierarchical name for this entity.

        For example:
            - Structure:       "1ABC"
            - Model:           "1ABC/0"
            - Chain:           "1ABC/0/A"
            - Residue:         "1ABC/0/A/GLY/42"
            - Atom:            override in Atom class

        Returns
        -------
        str
            Hierarchical identifier.
        """
        full_id = self.get_full_id()

        if self.level == "S":
            return str(full_id[0])

        if self.level == "M":
            return f"{full_id[0]}/{full_id[1]}"

        if self.level == "C":
            return f"{full_id[0]}/{full_id[1]}/{full_id[2]}"

        if self.level == "R":
            parent_path = '/'.join(str(x) for x in full_id[0:3])
            resname = f"{self.resname}/{self.id[1]}{self.id[2].strip()}"
            return f"{parent_path}/{resname}"

        # Default fallback
        return "/".join(str(x) for x in full_id)

    def get_parent_by_level(self, level):
        """
        Traverse up the SMCRA hierarchy to retrieve a parent at a specified level.

        Parameters
        ----------
        level : str
            One of: "A", "R", "C", "M", "S"

        Returns
        -------
        Entity
            The parent entity at the requested level.

        Raises
        ------
        ValueError
            If the level is invalid or not higher than the current level.
        """
        levels = ("A", "R", "C", "M", "S")

        if level not in levels:
            raise ValueError(f"Level must be one of {levels}.")
        if levels.index(level) < levels.index(self.level):
            raise ValueError(
                f"Cannot traverse up to '{level}' from level '{self.level}'. "
                "Valid hierarchy: A < R < C < M < S."
            )
        if self.level == level:
            return self
        return self.parent.get_parent_by_level(level)