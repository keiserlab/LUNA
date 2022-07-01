from enum import Enum, auto

from luna.wrappers.base import BondType


import logging

logger = logging.getLogger()


METALS = ["Li", "Na", "K", "Rb", "Cs", "Fr", "Be", "Mg", "Ca", "Sr", "Ba",
          "Ra", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
          "Al", "Ga", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
          "Cd", "In", "Sn", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au",
          "Hg", "Tl", "Pb", "Bi"]
METAL_ATOM = "[%s]" % ",".join(METALS)


class ResidueStandard(Enum):
    """An enumeration of protonation states for standard residues."""

    # TODO: add other residues to the standard

    # HID as in Amber: histidine with hydrogen on the delta nitrogen.
    HID = auto()
    # HIE as in Amber: histidine with hydrogen on the epsilon nitrogen.
    HIE = auto()
    # HIP as in Amber: histidine with hydrogens on both nitrogens; this is positively charged.
    HIP = auto()

    # CYS as in Amber: protonated cysteine
    CYS = auto()
    # CYM as in Amber: deprotonated cysteine.
    CYM = auto()


# TODO: create a new MetalStandardiser or LigandStandardiser for common ligands.

class ResiduesStandardiser:
    """Standardize residues.

    Parameters
    ----------
    break_metal_bonds : bool
        If True, break covalent bonds with metals and correct the topology of the involved atoms.
    his_type : {:class:`~ResidueStandard.HID`, :class:`~ResidueStandard.HIE`, :class:`~ResidueStandard.HIP`}
        Define which histidine protonation state to use. Currently, this option is still not been used.

    """

    # TODO: create filters for specific patterns like HIS tautomers or CYS:S- vs CYS:SH

    # TODO: some residues still have no WARNING for unexpected atoms bound to the atoms N, O, S.

    # TODO: I need to verify for HIS tautomers and metals around it, which will influenciate on the Hydrogen placement.

    # TODO: Metals will influentiate where the Hydrogen should be placed in TYR.

    def __init__(self, break_metal_bonds=False, his_type=ResidueStandard.HIE):

        self.break_metal_bonds = break_metal_bonds

        # TODO: not being used
        self.his_type = his_type

    def standardise(self, atom_pairs):
        """Standardize residues.

        Parameters
        ----------
        atom_pairs : iterable of tuple of (:class:`~luna.wrappers.base.MolWrapper`, :class:`~luna.MyBio.PDB.Atom.Atom`)
            The atoms to standardise.
        """
        pdb_map = {atm_obj.get_idx(): pdb_atm for atm_obj, pdb_atm in atom_pairs}

        self.found_metals = {}
        self.removed_atoms = []

        for atm_obj, pdb_atm in atom_pairs:

            # Atom: N
            # Sanity check for all N nitrogens with invalid bond types.
            #
            # E.g.: 3QQL:A:GLY:11.
            if pdb_atm.name == "N":

                # Non-N-terminal N where the carbonyl O is involved in a covalent bond with an atom not comprised
                # by the metal rules. So, it's better not to update anything here.
                if atm_obj.matches_smarts(f"N-,=[C;X4,X3]([#6])[OX2;$(O([C;X4,X3])[!#1]);!$(O{METAL_ATOM})]"):
                    logger.debug("While checking for inconsistencies in the atom N of the residue %s, it was found an unexpected "
                                 "atom covalently bound to the neighboring atom O. So, N will not be amended." % pdb_atm.parent)
                    continue
                # Non-N-terminal N where the carbonyl O is involved in a covalent bond with a metal.
                # However, 'break_metal_bonds' was set to False and there isn't anything to be done so.
                elif atm_obj.matches_smarts(f"N-,=[C;X4,X3]([#6])[OX2]{METAL_ATOM}") and self.break_metal_bonds is False:
                    logger.debug("While checking for inconsistencies in the atom N of the residue %s, it was found a metal "
                                 "covalently bound to its atom O. However, nothing will be done because 'break_metal_bonds' "
                                 "was set to False." % pdb_atm.parent)
                    continue

                # Any non-N-terminal N not from PRO.
                if (pdb_atm.parent.resname != "PRO" and atm_obj.get_neighbors_number(True) == 2
                    and (atm_obj.get_valence() != 3 or atm_obj.get_charge() != 0
                         or atm_obj.has_only_bond_type(BondType.SINGLE) is False
                         or atm_obj.get_h_count() != 1 or atm_obj.is_in_ring() or atm_obj.is_aromatic())):

                    self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=1)

                # Any N-terminal N not from PRO.
                elif (pdb_atm.parent.resname != "PRO" and atm_obj.get_neighbors_number(True) == 1
                        and (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 1
                             or atm_obj.has_only_bond_type(BondType.SINGLE) is False
                             or atm_obj.get_h_count() != 3 or atm_obj.is_in_ring() or atm_obj.is_aromatic())):

                    self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=1, implicit_h_count=3)

                # Non-N-terminal N from PRO.
                elif (pdb_atm.parent.resname == "PRO" and atm_obj.get_neighbors_number(True) == 3
                        and (atm_obj.get_valence() != 3 or atm_obj.get_charge() != 0
                             or atm_obj.has_only_bond_type(BondType.SINGLE) is False
                             or atm_obj.get_h_count() != 0 or atm_obj.is_in_ring() is False or atm_obj.is_aromatic())):

                    self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=0, in_ring=True)

                # N-terminal N from PRO.
                elif (pdb_atm.parent.resname == "PRO" and atm_obj.get_neighbors_number(True) == 2
                        and (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 1
                             or atm_obj.has_only_bond_type(BondType.SINGLE) is False
                             or atm_obj.get_h_count() != 2 or atm_obj.is_in_ring() is False or atm_obj.is_aromatic())):

                    self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=1, implicit_h_count=2, in_ring=True)

            # Atom: C.
            # Sanity check for all C carbons with invalid bond types.
            elif pdb_atm.name == "C":

                fix_atom = False

                # If it is bonded to 3 heavy atoms.
                if atm_obj.get_neighbors_number(True) == 3:
                    # Main chain carbonyl O involved in a metallic coordination.
                    #
                    #   It identifies only the form generated by Open Babel, where the double bond between O and C becomes
                    #       a single bond, and it adds an additional single bond with a Metal.
                    #
                    #   Note that if the list of atoms is sorted in a way that O comes first than C in this loop, then the
                    #       double bond would have already been amended and the bond with the metal would have already been removed.
                    #
                    #   E.g.: 6JWU:A:GLU:42.
                    #         1DKA:A:PRO:99.
                    #
                    if atm_obj.matches_smarts(f"[C;X4,X3]([$([OX2]{METAL_ATOM})])([#6])(-,=[N,O])"):
                        # If it is necessary to break covalent bonds with metals.
                        if self.break_metal_bonds:
                            fix_atom = True
                        else:
                            logger.debug("While checking for inconsistencies in the atom C of the residue %s, it was found a metal "
                                         "covalently bound to its atom O. However, nothing will be done because 'break_metal_bonds' "
                                         "was set to False." % pdb_atm.parent)

                    # If the neighboring oxygen is bound to something else not comprised in the previous rule, it is better not to
                    # update anything. Otherwise, fix the C.
                    elif (atm_obj.matches_smarts(f"[C;X4,X3]([OX2;$(O([C;X4,X3])[!#1]);!$(O{METAL_ATOM})])([#6])(-,=[N,O])") is False
                            and atm_obj.matches_smarts("[CX3](=[OX1])([#6])[NX3]") is False
                            and atm_obj.matches_smarts("[CX3](=[OX1])([#6])[O;H1,H0&-1]") is False):
                        fix_atom = True

                    # Alert for unexpected atoms bound to the oxygen O.
                    elif atm_obj.matches_smarts(f"[C;X4,X3]([OX2;$(O([C;X4,X3])[!#1]);!$(O{METAL_ATOM})])([#6])(-,=[N,O])"):
                        logger.debug("While checking for inconsistencies in the atom C of the residue %s, it was found an unexpected "
                                     "atom covalently bound to its atom O. So, C will not be amended." % pdb_atm.parent)

                    if fix_atom:
                        bond_types = []
                        for bond_obj in atm_obj.get_bonds():
                            partner_obj = bond_obj.get_partner_atom(atm_obj)
                            if partner_obj.get_idx() in pdb_map and pdb_map[partner_obj.get_idx()].name == "O":
                                bond_types.append((bond_obj, BondType.DOUBLE))
                            else:
                                bond_types.append((bond_obj, BondType.SINGLE))

                        self._fix_atom(atm_obj, bond_types=bond_types, charge=0, implicit_h_count=0)

                # If it is bonded only to 2 heavy atoms, then it is a C-terminal with a missing OXT or there are missing residues.
                elif atm_obj.get_neighbors_number(True) == 2:

                    logger.debug("The residue %s seems not to have any successor residue. "
                                 "It may be the last in the chain sequence or there are missing residues." % pdb_atm.parent)

                    # Main chain carbonyl O involved in a metallic coordination.
                    #
                    #   It identifies only the form generated by Open Babel, where the double bond between O and C becomes
                    #       a single bond, and it adds an additional single bond with a Metal.
                    #
                    #   Note that if the list of atoms is sorted in a way that O comes first than C in this loop, then the
                    #       double bond would have already been amended and the bond with the metal would have already been removed.
                    #
                    #   E.g.: 6JWU:A:GLU:42.
                    #         1DKA:A:THR:98.
                    #
                    if atm_obj.matches_smarts(f"[CX4]([$([OX2]{METAL_ATOM})])([#6])"):
                        # If it is necessary to break covalent bonds with metals.
                        if self.break_metal_bonds:
                            fix_atom = True
                        else:
                            logger.debug("While checking for inconsistencies in the atom C of the residue %s, it was found a metal "
                                         "covalently bound to its atom O. However, nothing will be done because 'break_metal_bonds' "
                                         "was set to False." % pdb_atm.parent)

                    # If the neighboring oxygen is bound to something else not comprised in the previous rule, it is better not to
                    # update anything. Otherwise, fix the C.
                    elif (atm_obj.matches_smarts(f"[CX4]([OX2;$(O([CX4])[!#1]);!$(O{METAL_ATOM})])([#6])") is False and
                            atm_obj.matches_smarts("[CX3](=[OX1])([#6])") is False):
                        fix_atom = True

                    # Alert for unexpected atoms bound to the oxygen O.
                    elif atm_obj.matches_smarts(f"[CX4]([OX2;$(O([CX4])[!#1]);!$(O{METAL_ATOM})])([#6])"):
                        logger.debug("While checking for inconsistencies in the atom C of the residue %s, it was found an unexpected "
                                     "atom covalently bound to its atom O. So, C will not be amended." % pdb_atm.parent)

                    if fix_atom:
                        bond_types = []
                        for bond_obj in atm_obj.get_bonds():
                            partner_obj = bond_obj.get_partner_atom(atm_obj)
                            if partner_obj.get_idx() in pdb_map and pdb_map[partner_obj.get_idx()].name == "O":
                                bond_types.append((bond_obj, BondType.DOUBLE))
                            else:
                                bond_types.append((bond_obj, BondType.SINGLE))

                        self._fix_atom(atm_obj, bond_types=bond_types, charge=0, implicit_h_count=1)

            # Atom: CA
            # Sanity check for all CA with invalid bond types.
            # E.g.: 3QQL:A:GLY:11.
            elif pdb_atm.name == "CA":

                # Any CA not from GLY/PRO.
                if (pdb_atm.parent.resname not in ["GLY", "PRO"] and atm_obj.get_neighbors_number(True) == 3 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.has_only_bond_type(BondType.SINGLE) is False or
                            atm_obj.get_h_count() != 1 or atm_obj.is_in_ring() or atm_obj.is_aromatic())):

                    self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=1)

                elif (pdb_atm.parent.resname == "PRO" and atm_obj.get_neighbors_number(True) == 3 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.has_only_bond_type(BondType.SINGLE) is False or
                            atm_obj.get_h_count() != 1 or atm_obj.is_in_ring() is False or atm_obj.is_aromatic())):

                    self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=1, in_ring=True)

                # CA from GLY.
                elif (pdb_atm.parent.resname == "GLY" and atm_obj.get_neighbors_number(True) == 2 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.has_only_bond_type(BondType.SINGLE) is False or
                            atm_obj.get_h_count() != 2 or atm_obj.is_in_ring() or atm_obj.is_aromatic())):

                    self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=2)

            # Atom: CB
            # Sanity check for all CB carbons with invalid bond types.
            # E.g.: 1THA:B:GLN:63.
            elif pdb_atm.name == "CB":

                # Any CB not from ALA, ILE, PRO, THR, VAL.
                if (pdb_atm.parent.resname not in ["ALA", "ILE", "PRO", "THR", "VAL"] and atm_obj.get_neighbors_number(True) == 2 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.has_only_bond_type(BondType.SINGLE) is False or
                            atm_obj.get_h_count() != 2 or atm_obj.is_in_ring() or atm_obj.is_aromatic())):

                    self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=2)

                # CB from ILE, THR, VAL.
                elif (pdb_atm.parent.resname in ["ILE", "THR", "VAL"] and atm_obj.get_neighbors_number(True) == 3 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.has_only_bond_type(BondType.SINGLE) is False or
                            atm_obj.get_h_count() != 1 or atm_obj.is_in_ring() or atm_obj.is_aromatic())):

                    self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=1)

                # CB from ALA.
                elif (pdb_atm.parent.resname == "ALA" and atm_obj.get_neighbors_number(True) == 1 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.has_only_bond_type(BondType.SINGLE) is False or
                            atm_obj.get_h_count() != 3 or atm_obj.is_in_ring() or atm_obj.is_aromatic())):

                    self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=3)

                # CB from PRO.
                elif (pdb_atm.parent.resname == "PRO" and atm_obj.get_neighbors_number(True) == 2 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.has_only_bond_type(BondType.SINGLE) is False or
                            atm_obj.get_h_count() != 2 or atm_obj.is_in_ring() is False or atm_obj.is_aromatic())):

                    self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=2, in_ring=True)

            # Atom: CG
            # Sanity check for all CG carbons with invalid bond types.
            elif pdb_atm.name == "CG":

                # E.g.: 1THA:B:GLN:63, 6COD:A:LYS:121.
                if (pdb_atm.parent.resname in ["ARG", "GLN", "GLU", "LYS", "MET"] and atm_obj.get_neighbors_number(True) == 2 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.has_only_bond_type(BondType.SINGLE) is False or
                            atm_obj.get_h_count() != 2 or atm_obj.is_in_ring() or atm_obj.is_aromatic())):

                    self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=2)

                elif (pdb_atm.parent.resname == "PRO" and atm_obj.get_neighbors_number(True) == 2 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.has_only_bond_type(BondType.SINGLE) is False or
                            atm_obj.get_h_count() != 2 or atm_obj.is_in_ring() is False or atm_obj.is_aromatic())):

                    self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=2, in_ring=True)

                elif (pdb_atm.parent.resname in ["ASN", "ASP"] and atm_obj.get_neighbors_number(True) == 3 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.get_degree() != 3 or
                            atm_obj.get_h_count() != 0 or atm_obj.is_in_ring() or atm_obj.is_aromatic())):

                    fix_atom = False
                    # ASN/ASP with metallic coordination perceived as covalent bond.
                    #
                    #   It identifies only the form generated by Open Babel, where the double bond between OD1 and CG becomes
                    #       a single bond, and it adds an additional single bond with a Metal.
                    #
                    #   Note that if the list of atoms is sorted in a way that OD1 comes first than CG in this loop, then the
                    #       double bond would have already been amended and the bond with the metal would have already been removed.
                    #
                    #   E.g.: 6JWU:B:ASP:210, 4FVR:A:ASN:678
                    #
                    if (atm_obj.matches_smarts(f"[CX4](N)([#6])[OX2]{METAL_ATOM}") or
                            atm_obj.matches_smarts(f"[C;X4,X3](-,=[O])([#6])[OX2]({METAL_ATOM})")):
                        # If it is necessary to break covalent bonds with metals.
                        if self.break_metal_bonds:
                            fix_atom = True
                        else:
                            logger.debug("While checking for inconsistencies in the atom CG of the residue %s, it was found a metal "
                                         "covalently bound to its atom OD1. However, nothing will be done because 'break_metal_bonds' "
                                         "was set to False." % pdb_atm.parent)

                    # If the neighboring oxygen is bound to something else not comprised in the previous rule, it is better not to
                    # update anything. Otherwise, fix the CG.
                    elif (atm_obj.matches_smarts("[CX4](N)([#6])[OX2][!#1]") is False and
                            atm_obj.matches_smarts("[C;X4,X3](-,=[O])([#6])[OX2]([!#1])") is False):
                        fix_atom = True

                    # Alert for unexpected atoms bound to the neighboring nitrogen/oxygen.
                    elif atm_obj.matches_smarts("[C;X4,X3](-,=[N,O])([#6])[OX2][!#1]"):
                        logger.debug("While checking for inconsistencies in the atom CG of the residue %s, it was found an unexpected "
                                     "atom covalently bound to a neighboring nitrogen/oxygen. So, CG will not be amended."
                                     % pdb_atm.parent)

                    if fix_atom:
                        bond_types = []
                        for bond_obj in atm_obj.get_bonds():
                            partner_obj = bond_obj.get_partner_atom(atm_obj)

                            if partner_obj.get_idx() in pdb_map:
                                # Redundant rules: OD1 and OD2 section already fixes the bonds with CG.
                                if partner_obj.get_idx() in pdb_map and pdb_map[partner_obj.get_idx()].name == "OD1":
                                    # If the OD1 atom is bound to the below metals, update the CG-OD1 bond to CG=OD1.
                                    if partner_obj.matches_smarts(f"[OX2]({METAL_ATOM})"):
                                        bond_types.append((bond_obj, BondType.DOUBLE))
                                    else:
                                        bond_types.append((bond_obj, BondType.SINGLE))

                                # Redundant rules: OD1 and OD2 section already fixes the bonds with CG.
                                elif partner_obj.get_idx() in pdb_map and pdb_map[partner_obj.get_idx()].name == "OD2":
                                    # If the second oxygen (OD1) is bound to some atom not comprised by the expected metals, but OD2 is,
                                    # then update the CG-OD2 bond to CG=OD2.
                                    if partner_obj.matches_smarts(f"[OX2]({METAL_ATOM})[CX4]([#6])"
                                                                  f"[OX2;$(O([CX4])[!#1]);!$(O{METAL_ATOM})]"):
                                        bond_types.append((bond_obj, BondType.DOUBLE))
                                    else:
                                        bond_types.append((bond_obj, BondType.SINGLE))

                                # Any other atom bound to CG
                                else:
                                    bond_types.append((bond_obj, BondType.SINGLE))

                        self._fix_atom(atm_obj, bond_types=bond_types, charge=0, implicit_h_count=0)

                elif (pdb_atm.parent.resname == "LEU" and atm_obj.get_neighbors_number(True) == 3 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.has_only_bond_type(BondType.SINGLE) is False or
                            atm_obj.get_h_count() != 1 or atm_obj.is_in_ring() or atm_obj.is_aromatic())):

                    self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=1)

                elif (pdb_atm.parent.resname in ["HIS", "PHE", "TRP", "TYR"] and atm_obj.get_neighbors_number(True) == 3 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.get_degree() != 3 or
                            atm_obj.get_h_count() != 0 or atm_obj.is_aromatic() is False)):

                    bond_types = []
                    for bond_obj in atm_obj.get_bonds():
                        partner_obj = bond_obj.get_partner_atom(atm_obj)
                        if partner_obj.get_idx() in pdb_map:
                            # Single bonds with CB (all)
                            if pdb_map[partner_obj.get_idx()].name == "CB":
                                bond_types.append((bond_obj, BondType.SINGLE))
                            # Single bond with ND1 (only in HIS)
                            elif pdb_map[partner_obj.get_idx()].name == "ND1":
                                bond_types.append((bond_obj, BondType.SINGLE, True))
                            # Double bonds with CD1 (all, except HIS)
                            elif pdb_map[partner_obj.get_idx()].name == "CD1":
                                bond_types.append((bond_obj, BondType.DOUBLE, True))
                            # Bonds with CD2
                            elif pdb_map[partner_obj.get_idx()].name == "CD2":
                                # Double bond as in HIS.
                                if pdb_atm.parent.resname == "HIS":
                                    bond_types.append((bond_obj, BondType.DOUBLE, True))
                                # Single bonds (all, except HIS)
                                else:
                                    bond_types.append((bond_obj, BondType.SINGLE, True))

                    self._fix_atom(atm_obj, bond_types=bond_types, charge=0, implicit_h_count=0, is_aromatic=True)

            # Atom: CG1
            # Sanity check for all CG1 carbons with invalid bond types.
            elif pdb_atm.name == "CG1":
                # CG1 from ILE.
                if (pdb_atm.parent.resname == "ILE" and atm_obj.get_neighbors_number(True) == 2 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.has_only_bond_type(BondType.SINGLE) is False or
                            atm_obj.get_h_count() != 2 or atm_obj.is_in_ring() or atm_obj.is_aromatic())):

                    self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=2)

                # CG1 from VAL.
                elif (pdb_atm.parent.resname == "VAL" and atm_obj.get_neighbors_number(True) == 1 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.has_only_bond_type(BondType.SINGLE) is False or
                            atm_obj.get_h_count() != 3 or atm_obj.is_in_ring() or atm_obj.is_aromatic())):

                    self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=3)

            # Atom: CG2
            # Sanity check for all CG2 carbons with invalid bond types.
            elif (pdb_atm.name == "CG2" and atm_obj.get_neighbors_number(True) == 1 and
                    (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.has_only_bond_type(BondType.SINGLE) is False or
                        atm_obj.get_h_count() != 3 or atm_obj.is_in_ring() or atm_obj.is_aromatic())):

                self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=3)

            # Atom: CD
            # Sanity check for all CD carbons with invalid bond types.
            elif pdb_atm.name == "CD":

                if (pdb_atm.parent.resname in ["ARG", "LYS"] and atm_obj.get_neighbors_number(True) == 2 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.has_only_bond_type(BondType.SINGLE) is False or
                            atm_obj.get_h_count() != 2 or atm_obj.is_in_ring() or atm_obj.is_aromatic())):

                    self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=2)

                # E.g.: 4FVR:A:PRO:792
                elif (pdb_atm.parent.resname == "PRO" and atm_obj.get_neighbors_number(True) == 2 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.has_only_bond_type(BondType.SINGLE) is False or
                            atm_obj.get_h_count() != 2 or atm_obj.is_in_ring() is False or atm_obj.is_aromatic())):

                    self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=2, in_ring=True)

                # E.g.: 1THA:B:GLN:63
                # E.g.: 6JWU:A:GLU:42 containing metallic cordination.
                elif (pdb_atm.parent.resname in ["GLN", "GLU"] and atm_obj.get_neighbors_number(True) == 3 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.get_degree() != 3 or
                            atm_obj.get_h_count() != 0 or atm_obj.is_in_ring() or atm_obj.is_aromatic())):

                    fix_atom = False
                    # GLN/GLU with metallic coordination perceived as covalent bond.
                    #
                    #   It identifies only the form generated by Open Babel, where the double bond between OE1 and CD becomes
                    #       a single bond, and it adds an additional single bond with a Metal.
                    #
                    #   Note that if the list of atoms is sorted in a way that OE1 comes first than CD in this loop, then the
                    #       double bond would have already been amended and the bond with the metal would have already been removed.
                    #
                    if (atm_obj.matches_smarts(f"[CX4](N)([#6])[OX2]{METAL_ATOM}") or
                            atm_obj.matches_smarts(f"[C;X4,X3](-,=[O])([#6])[OX2]({METAL_ATOM})")):
                        # If it is necessary to break covalent bonds with metals.
                        if self.break_metal_bonds:
                            fix_atom = True
                        else:
                            logger.debug("While checking for inconsistencies in the atom CD of the residue %s, it was found a metal "
                                         "covalently bound to its atom OE1. However, nothing will be done because 'break_metal_bonds' "
                                         "was set to False." % pdb_atm.parent)

                    # If the neighboring oxygen is bound to something else not comprised in the previous rule, it is better not to
                    # update anything. Otherwise, fix the CD.
                    elif (atm_obj.matches_smarts("[CX4](N)([#6])[OX2][!#1]") is False and
                            atm_obj.matches_smarts("[C;X4,X3](-,=[O])([#6])[OX2][!#1]") is False):
                        fix_atom = True

                    # Alert for unexpected atoms bound to the neighboring nitrogen/oxygen.
                    elif atm_obj.matches_smarts("[C;X4,X3](-,=[N,O])([#6])[OX2][!#1]"):
                        logger.debug("While checking for inconsistencies in the atom CD of the residue %s, it was found an unexpected "
                                     "atom covalently bound to a neighboring nitrogen/oxygen. So, CD will not be amended."
                                     % pdb_atm.parent)

                    if fix_atom:
                        bond_types = []
                        for bond_obj in atm_obj.get_bonds():
                            partner_obj = bond_obj.get_partner_atom(atm_obj)

                            if partner_obj.get_idx() in pdb_map:
                                # Redundant rules: OE1 and OE2 section already fixes the bonds with CD.
                                if partner_obj.get_idx() in pdb_map and pdb_map[partner_obj.get_idx()].name == "OE1":
                                    # If the OE1 atom is bound to the below metals, update the CD-OE1 bond to CD=OE1.
                                    if partner_obj.matches_smarts(f"[OX2]({METAL_ATOM})"):
                                        bond_types.append((bond_obj, BondType.DOUBLE))
                                    else:
                                        bond_types.append((bond_obj, BondType.SINGLE))

                                # Redundant rules: OE1 and OE2 section already fixes the bonds with CD.
                                elif partner_obj.get_idx() in pdb_map and pdb_map[partner_obj.get_idx()].name == "OE2":
                                    # If the second oxygen (OE1) is bound to some atom not comprised by the expected metals, but OE2 is,
                                    # then update the CD-OE2 bond to CD=OE2.
                                    if partner_obj.matches_smarts(f"[OX2]({METAL_ATOM})[CX4]([#6])"
                                                                  f"[OX2;$(O([CX4])[!#1]);!$(O{METAL_ATOM})]"):
                                        bond_types.append((bond_obj, BondType.DOUBLE))
                                    else:
                                        bond_types.append((bond_obj, BondType.SINGLE))

                                # Any other atom bound to CD
                                else:
                                    bond_types.append((bond_obj, BondType.SINGLE))

                        self._fix_atom(atm_obj, bond_types=bond_types, charge=0, implicit_h_count=0)

            # Atom: CD1
            # Sanity check for all CD1 carbons with invalid bond types.
            elif pdb_atm.name == "CD1":

                if (pdb_atm.parent.resname in ["ILE", "LEU"] and atm_obj.get_neighbors_number(True) == 1 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.has_only_bond_type(BondType.SINGLE) is False or
                            atm_obj.get_h_count() != 3 or atm_obj.is_in_ring() or atm_obj.is_aromatic())):

                    self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=3)

                elif (pdb_atm.parent.resname in ["PHE", "TRP", "TYR"] and atm_obj.get_neighbors_number(True) == 2 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.get_degree() != 3 or
                            atm_obj.get_h_count() != 1 or atm_obj.is_aromatic() is False)):

                    bond_types = []
                    for bond_obj in atm_obj.get_bonds():
                        partner_obj = bond_obj.get_partner_atom(atm_obj)
                        if partner_obj.get_idx() in pdb_map:
                            # Double bonds with CG
                            if pdb_map[partner_obj.get_idx()].name == "CG":
                                bond_types.append((bond_obj, BondType.DOUBLE, True))
                            # Single bonds with CE1 and NE1
                            else:
                                bond_types.append((bond_obj, BondType.SINGLE, True))

                    self._fix_atom(atm_obj, bond_types=bond_types, charge=0, implicit_h_count=1, is_aromatic=True)

            # Atom: CD2
            # Sanity check for all CD2 carbons with invalid bond types.
            elif pdb_atm.name == "CD2":

                if (pdb_atm.parent.resname == "LEU" and atm_obj.get_neighbors_number(True) == 1 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.has_only_bond_type(BondType.SINGLE) is False or
                            atm_obj.get_h_count() != 3 or atm_obj.is_in_ring() or atm_obj.is_aromatic())):

                    self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=3)

                elif (pdb_atm.parent.resname in ["HIS", "PHE", "TYR"] and atm_obj.get_neighbors_number(True) == 2 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.get_degree() != 3 or
                            atm_obj.get_h_count() != 1 or atm_obj.is_aromatic() is False)):

                    bond_types = []
                    for bond_obj in atm_obj.get_bonds():
                        partner_obj = bond_obj.get_partner_atom(atm_obj)
                        if partner_obj.get_idx() in pdb_map:
                            # Bonds with CG
                            if pdb_map[partner_obj.get_idx()].name == "CG":
                                # Double bond with CG (only in HIS)
                                if pdb_atm.parent.resname == "HIS":
                                    bond_types.append((bond_obj, BondType.DOUBLE, True))
                                # Single bonds with CG (all, except HIS)
                                else:
                                    bond_types.append((bond_obj, BondType.SINGLE, True))
                            # Double bonds with CE2 (all, except HIS)
                            elif pdb_map[partner_obj.get_idx()].name == "CE2":
                                bond_types.append((bond_obj, BondType.DOUBLE, True))
                            # Single bonds with NE2 (only HIS)
                            elif pdb_map[partner_obj.get_idx()].name == "NE2":
                                bond_types.append((bond_obj, BondType.SINGLE, True))

                    self._fix_atom(atm_obj, bond_types=bond_types, charge=0, implicit_h_count=1, is_aromatic=True)

                elif (pdb_atm.parent.resname == "TRP" and atm_obj.get_neighbors_number(True) == 3 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.get_degree() != 3 or
                            atm_obj.get_h_count() != 0 or atm_obj.is_aromatic() is False)):

                    bond_types = []
                    for bond_obj in atm_obj.get_bonds():
                        partner_obj = bond_obj.get_partner_atom(atm_obj)
                        if partner_obj.get_idx() in pdb_map:
                            # Double bonds with CE2
                            if pdb_map[partner_obj.get_idx()].name == "CE2":
                                bond_types.append((bond_obj, BondType.DOUBLE, True))
                            # Single bonds with CG and CE3
                            else:
                                bond_types.append((bond_obj, BondType.SINGLE, True))

                    self._fix_atom(atm_obj, bond_types=bond_types, charge=0, implicit_h_count=0, is_aromatic=True)

            # Atom: CE
            # Sanity check for all CE carbons with invalid bond types.
            elif pdb_atm.name == "CE":

                # CE from LYS.
                if (pdb_atm.parent.resname == "LYS" and atm_obj.get_neighbors_number(True) == 2 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.has_only_bond_type(BondType.SINGLE) is False or
                            atm_obj.get_h_count() != 2 or atm_obj.is_in_ring() or atm_obj.is_aromatic())):

                    self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=2)

                # CE from MET.
                elif (pdb_atm.parent.resname == "MET" and atm_obj.get_neighbors_number(True) == 1 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.has_only_bond_type(BondType.SINGLE) is False or
                            atm_obj.get_h_count() != 3 or atm_obj.is_in_ring() or atm_obj.is_aromatic())):

                    self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=3)

            # Atom: CE1
            # Sanity check for all CE1 carbons with invalid bond types.
            elif pdb_atm.name == "CE1":

                if (pdb_atm.parent.resname in ["HIS", "PHE", "TYR"] and atm_obj.get_neighbors_number(True) == 2 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.get_degree() != 3 or
                            atm_obj.get_h_count() != 1 or atm_obj.is_aromatic() is False)):

                    bond_types = []
                    for bond_obj in atm_obj.get_bonds():
                        partner_obj = bond_obj.get_partner_atom(atm_obj)
                        if partner_obj.get_idx() in pdb_map:
                            # Single bonds with CD1 (all, except HYS)
                            if pdb_map[partner_obj.get_idx()].name == "CD1":
                                bond_types.append((bond_obj, BondType.SINGLE, True))
                            # Double bonds with CZ (all, except HYS)
                            elif pdb_map[partner_obj.get_idx()].name == "CZ":
                                bond_types.append((bond_obj, BondType.DOUBLE, True))
                            # Double bonds with ND1 (only HYS)
                            elif pdb_map[partner_obj.get_idx()].name == "ND1":
                                bond_types.append((bond_obj, BondType.DOUBLE, True))
                            # Single bonds with NE2 (only HYS)
                            elif pdb_map[partner_obj.get_idx()].name == "NE2":
                                # CE1 - NE2
                                bond_types.append((bond_obj, BondType.SINGLE, True))

                    self._fix_atom(atm_obj, bond_types=bond_types, charge=0, implicit_h_count=1, is_aromatic=True)

            # Atom: CE2
            # Sanity check for all CE2 carbons with invalid bond types.
            elif pdb_atm.name == "CE2":

                if (pdb_atm.parent.resname in ["PHE", "TYR"] and atm_obj.get_neighbors_number(True) == 2 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.get_degree() != 3 or
                            atm_obj.get_h_count() != 1 or atm_obj.is_aromatic() is False)):

                    bond_types = []
                    for bond_obj in atm_obj.get_bonds():
                        partner_obj = bond_obj.get_partner_atom(atm_obj)
                        if partner_obj.get_idx() in pdb_map:
                            # Double bonds with CD2
                            if pdb_map[partner_obj.get_idx()].name == "CD2":
                                bond_types.append((bond_obj, BondType.DOUBLE, True))
                            # Single bonds with CZ
                            elif pdb_map[partner_obj.get_idx()].name == "CZ":
                                bond_types.append((bond_obj, BondType.SINGLE, True))

                    self._fix_atom(atm_obj, bond_types=bond_types, charge=0, implicit_h_count=1, is_aromatic=True)

                elif (pdb_atm.parent.resname == "TRP" and atm_obj.get_neighbors_number(True) == 3 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.get_degree() != 3 or
                            atm_obj.get_h_count() != 0 or atm_obj.is_aromatic() is False)):

                    bond_types = []
                    for bond_obj in atm_obj.get_bonds():
                        partner_obj = bond_obj.get_partner_atom(atm_obj)
                        if partner_obj.get_idx() in pdb_map:
                            # Double bond with CD2
                            if pdb_map[partner_obj.get_idx()].name == "CD2":
                                bond_types.append((bond_obj, BondType.DOUBLE, True))
                            # Single bonds with NE1 and CZ2
                            else:
                                bond_types.append((bond_obj, BondType.SINGLE, True))

                    self._fix_atom(atm_obj, bond_types=bond_types, charge=0, implicit_h_count=0, is_aromatic=True)

            # Atom: CZ
            # Sanity check for all CZ carbons with invalid bond types.
            elif pdb_atm.name == "CZ":

                if (pdb_atm.parent.resname == "PHE" and atm_obj.get_neighbors_number(True) == 2 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.get_degree() != 3 or
                            atm_obj.get_h_count() != 1 or atm_obj.is_aromatic() is False)):

                    bond_types = []
                    for bond_obj in atm_obj.get_bonds():
                        partner_obj = bond_obj.get_partner_atom(atm_obj)
                        if partner_obj.get_idx() in pdb_map:
                            # Double bond with CE1
                            if pdb_map[partner_obj.get_idx()].name == "CE1":
                                bond_types.append((bond_obj, BondType.DOUBLE, True))
                            # Single bond with CE2
                            else:
                                bond_types.append((bond_obj, BondType.SINGLE, True))

                    self._fix_atom(atm_obj, bond_types=bond_types, charge=0, implicit_h_count=1, is_aromatic=True)

                elif (pdb_atm.parent.resname == "TYR" and atm_obj.get_neighbors_number(True) == 3 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.get_degree() != 3 or
                            atm_obj.get_h_count() != 0 or atm_obj.is_aromatic() is False)):

                    bond_types = []
                    for bond_obj in atm_obj.get_bonds():
                        partner_obj = bond_obj.get_partner_atom(atm_obj)
                        if partner_obj.get_idx() in pdb_map:
                            # Double bond with CE1
                            if pdb_map[partner_obj.get_idx()].name == "CE1":
                                bond_types.append((bond_obj, BondType.DOUBLE, True))
                            # Single bond with CE2
                            elif pdb_map[partner_obj.get_idx()].name == "CE2":
                                bond_types.append((bond_obj, BondType.SINGLE, True))
                            # Single bond with OH
                            else:
                                bond_types.append((bond_obj, BondType.SINGLE))

                    self._fix_atom(atm_obj, bond_types=bond_types, charge=0, implicit_h_count=0, is_aromatic=True)

                elif (pdb_atm.parent.resname == "ARG" and atm_obj.get_neighbors_number(True) == 3 and
                        (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.get_degree() != 3 or
                            atm_obj.get_h_count() != 0 or atm_obj.is_in_ring() or atm_obj.is_aromatic())):

                    bond_types = []
                    for bond_obj in atm_obj.get_bonds():
                        partner_obj = bond_obj.get_partner_atom(atm_obj)
                        if partner_obj.get_idx() in pdb_map:
                            # Double bond with NH2
                            if pdb_map[partner_obj.get_idx()].name == "NH2":
                                bond_types.append((bond_obj, BondType.DOUBLE))
                            # Single bonds with NE and NH1
                            else:
                                bond_types.append((bond_obj, BondType.SINGLE))

                    self._fix_atom(atm_obj, bond_types=bond_types, charge=0, implicit_h_count=0)

            # Atom: CE3, CZ2, CZ3, CH2
            # Sanity check for TRP carbons with invalid bond types.
            elif (pdb_atm.name in ["CE3", "CZ2", "CZ3", "CH2"] and atm_obj.get_neighbors_number(True) == 2 and
                    (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 0 or atm_obj.get_degree() != 3 or
                        atm_obj.get_h_count() != 1 or atm_obj.is_aromatic() is False)):

                bond_types = []
                for bond_obj in atm_obj.get_bonds():
                    partner_obj = bond_obj.get_partner_atom(atm_obj)
                    if partner_obj.get_idx() in pdb_map:

                        # Double bonds
                        if ((pdb_atm.name == "CE3" and pdb_map[partner_obj.get_idx()].name == "CZ3") or
                                (pdb_atm.name == "CZ3" and pdb_map[partner_obj.get_idx()].name == "CE3") or
                                (pdb_atm.name == "CZ2" and pdb_map[partner_obj.get_idx()].name == "CH2") or
                                (pdb_atm.name == "CH2" and pdb_map[partner_obj.get_idx()].name == "CZ2")):

                            bond_types.append((bond_obj, BondType.DOUBLE, True))
                        # Single bonds
                        else:
                            bond_types.append((bond_obj, BondType.SINGLE, True))

                self._fix_atom(atm_obj, bond_types=bond_types, charge=0, implicit_h_count=1, is_aromatic=True)

            # Atom: ND1
            # Sanity check for HIS:ND1 with invalid bond types.
            elif pdb_atm.name == "ND1":

                fix_atom = False
                # HIS with metallic coordination perceived as covalent bond.
                #
                # Sometimes, the aromatic ring in HIS becomes a simple ring, therefore we should check for any
                # type of nitrogen in a ring, i.e, aromatic or aliphatic. Note that the aromatic ring
                # will be fixed after removing the bond with the metal.
                #
                # E.g.: 1USN:A:HIS:179
                #       1IUZ:A:HIS:37
                #
                if atm_obj.matches_smarts(f"[#7;R]-{METAL_ATOM}"):
                    # If it is necessary to break covalent bonds with metals.
                    if self.break_metal_bonds:
                        self._remove_metallic_bond(atm_obj)
                        fix_atom = True
                    else:
                        logger.debug("While checking for inconsistencies in the atom ND1 of the residue %s, it was found a metal "
                                     "covalently bound to it. However, nothing will be done because 'break_metal_bonds' "
                                     "was set to False." % pdb_atm.parent)

                # If the nitrogen is bound to something else not comprised in the previous rule, it is better not to
                # update anything. Otherwise, fix ND1.
                elif (atm_obj.matches_smarts("[#7;R]([#6])([#6])[!#1]") is False and
                      atm_obj.get_neighbors_number(True) == 2 and (atm_obj.get_valence() != 3 or atm_obj.get_charge() != 0 or
                                                                   atm_obj.get_degree() != 2 or atm_obj.get_h_count() != 0 or
                                                                   atm_obj.is_aromatic() is False)):
                    fix_atom = True

                # Alert for unexpected atoms bound to the nitrogen ND1.
                elif atm_obj.matches_smarts("[#7;R]([#6])([#6])[!#1]"):
                    logger.debug("While checking for inconsistencies in the atom ND1 of the residue %s, it was found an unexpected "
                                 "atom covalently bound to it. So, ND1 will not be amended." % pdb_atm.parent)

                if fix_atom:
                    bond_types = []
                    for bond_obj in atm_obj.get_bonds():
                        partner_obj = bond_obj.get_partner_atom(atm_obj)
                        if partner_obj.get_idx() in pdb_map:
                            # Double bond with CE1
                            if pdb_map[partner_obj.get_idx()].name == "CE1":
                                bond_types.append((bond_obj, BondType.DOUBLE, True))
                            # Single bonds with CG
                            else:
                                bond_types.append((bond_obj, BondType.SINGLE, True))

                    self._fix_atom(atm_obj, bond_types=bond_types, charge=0, implicit_h_count=0, is_aromatic=True)

            elif (pdb_atm.name == "ND1" and atm_obj.get_neighbors_number(True) == 2 and
                    (atm_obj.get_valence() != 3 or atm_obj.get_charge() != 0 or atm_obj.get_degree() != 2 or
                        atm_obj.get_h_count() != 0 or atm_obj.is_aromatic() is False)):

                bond_types = []
                for bond_obj in atm_obj.get_bonds():
                    partner_obj = bond_obj.get_partner_atom(atm_obj)
                    if partner_obj.get_idx() in pdb_map:
                        # Double bond with CE1
                        if pdb_map[partner_obj.get_idx()].name == "CE1":
                            bond_types.append((bond_obj, BondType.DOUBLE, True))
                        # Single bonds with CG
                        else:
                            bond_types.append((bond_obj, BondType.SINGLE, True))

                self._fix_atom(atm_obj, bond_types=bond_types, charge=0, implicit_h_count=0, is_aromatic=True)

            # Atom: ND2
            # Sanity check for ASN:ND2 with invalid bond types.
            elif (pdb_atm.name == "ND2" and atm_obj.get_neighbors_number(True) == 1 and
                    (atm_obj.get_valence() != 3 or atm_obj.get_charge() != 0 or atm_obj.has_only_bond_type(BondType.SINGLE) is False or
                        atm_obj.get_h_count() != 2 or atm_obj.is_in_ring() or atm_obj.is_aromatic())):

                # ASN with metallic coordination perceived as covalent bond.
                #
                #   It identifies two forms:
                #       - the firt consists of the usual form generated by Open Babel, where the double bond between
                #           OD1 and CG becomes a single bond, and it adds an additional single bond with a Metal.
                #
                #       - the second form may appear after correcting CG first, so the missing double bond would have already been
                #           fixed.
                #
                #   E.g.: 4FVR:A:ASN:678
                #
                if atm_obj.matches_smarts(f"N[C;X4,X3]([#6])-,=[OX2]{METAL_ATOM}"):
                    # If it is necessary to break covalent bonds with metals.
                    if self.break_metal_bonds:
                        self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=2)
                    else:
                        logger.debug("While checking for inconsistencies in the atom ND2 of the residue %s, it was found a metal "
                                     "covalently bound to the atom OD1. However, nothing will be done because 'break_metal_bonds' "
                                     "was set to False." % pdb_atm.parent)

                # If the oxygen is bound to something else not comprised in the previous rules, it is better not to update anything.
                # Otherwise, fix the ND2.
                elif atm_obj.matches_smarts("N[CX4]([#6])[OX2][!#1]") is False:

                    self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=2)

                # Alert for unexpected atoms bound to the oxygen OD1.
                elif atm_obj.matches_smarts("N[CX4]([#6])[OX2][!#1]"):
                    logger.debug("While checking for inconsistencies in the atom ND2 of the residue %s, it was found an unexpected "
                                 "atom covalently bound to the oxygen OD1. So, ND2 will not be amended." % pdb_atm.parent)

            # Atom: NE
            # Sanity check for ARG:NE with invalid bond types.
            elif (pdb_atm.name == "NE" and atm_obj.get_neighbors_number(True) == 2 and
                    (atm_obj.get_valence() != 3 or atm_obj.get_charge() != 0 or atm_obj.has_only_bond_type(BondType.SINGLE) is False or
                        atm_obj.get_h_count() != 1 or atm_obj.is_in_ring() or atm_obj.is_aromatic())):

                self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=1)

            # Atom: NE1
            # Sanity check for TRP:NE1 with invalid bond types.
            elif (pdb_atm.name == "NE1" and atm_obj.get_neighbors_number(True) == 2 and
                    (atm_obj.get_valence() != 3 or atm_obj.get_charge() != 0 or atm_obj.get_degree() != 3 or
                        atm_obj.get_h_count() != 1 or atm_obj.is_aromatic() is False)):

                bond_types = [(bond_obj, BondType.SINGLE, True) for bond_obj in atm_obj.get_bonds()]
                self._fix_atom(atm_obj, bond_types=bond_types, charge=0, implicit_h_count=1, is_aromatic=True)

            # Atom: NE2
            # Sanity check for NE2 nitrogens with invalid bond types.
            elif pdb_atm.name == "NE2":

                if pdb_atm.parent.resname == "HIS":
                    # HIS with metallic coordination perceived as covalent bond.
                    #
                    # Sometimes, the aromatic ring in HIS becomes a simple ring, therefore we check for any
                    # type of nitrogen in a ring, i.e, aromatic or aliphatic. Note that the aromatic ring
                    # will be fixed after removing the bond with the metal.
                    #
                    # E.g.: 1USN:A:HIS:151
                    #       1USN:A:HIS:166
                    if atm_obj.matches_smarts(f"[#7;R]{METAL_ATOM}"):
                        # If it is necessary to break covalent bonds with metals.
                        if self.break_metal_bonds:
                            self._remove_metallic_bond(atm_obj)

                            bond_types = [(bond_obj, BondType.SINGLE, True) for bond_obj in atm_obj.get_bonds()]
                            self._fix_atom(atm_obj, bond_types=bond_types, charge=0, implicit_h_count=1, is_aromatic=True)
                        else:
                            logger.debug("While checking for inconsistencies in the atom NE2 of the residue %s, it was found a metal "
                                         "covalently bound to it. However, nothing will be done because 'break_metal_bonds' "
                                         "was set to False." % pdb_atm.parent)

                    # If the nitrogen is bound to something else not comprised in the previous rule, it is better not to
                    # update anything. Otherwise, fix NE2.
                    elif (atm_obj.matches_smarts("[#7;R]([#6])([#6])[!#1]") is False and
                          atm_obj.get_neighbors_number(True) == 2 and (atm_obj.get_valence() != 3 or atm_obj.get_charge() != 0 or
                                                                       atm_obj.get_degree() != 3 or atm_obj.get_h_count() != 1 or
                                                                       atm_obj.is_aromatic() is False)):

                        bond_types = [(bond_obj, BondType.SINGLE, True) for bond_obj in atm_obj.get_bonds()]
                        self._fix_atom(atm_obj, bond_types=bond_types, charge=0, implicit_h_count=1, is_aromatic=True)

                    # Alert for unexpected atoms bound to the nitrogen NE2.
                    elif atm_obj.matches_smarts("[#7;R]([#6])([#6])[!#1]"):
                        logger.debug("While checking for inconsistencies in the atom NE2 of the residue %s, it was found an unexpected "
                                     "atom covalently bound to it. So, NE2 will not be amended." % pdb_atm.parent)

                elif (pdb_atm.parent.resname == "GLN" and atm_obj.get_neighbors_number(True) == 1 and
                        (atm_obj.get_valence() != 3 or atm_obj.get_charge() != 0 or atm_obj.has_only_bond_type(BondType.SINGLE) is False or
                            atm_obj.get_h_count() != 2 or atm_obj.is_in_ring() or atm_obj.is_aromatic())):

                    # GLN with metallic coordination perceived as covalent bond.
                    #
                    #   It identifies two forms:
                    #       - the firt consists of the usual form generated by Open Babel, where the double bond between
                    #           OE1 and CD becomes a single bond, and it adds an additional single bond with a Metal.
                    #
                    #       - the second form may appear after correcting CD first, so the missing double bond would have already been
                    #           fixed.
                    #
                    if atm_obj.matches_smarts(f"N[C;X4,X3]([#6])-,=[OX2]{METAL_ATOM}"):
                        # If it is necessary to break covalent bonds with metals.
                        if self.break_metal_bonds:
                            self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=2)
                        else:
                            logger.debug("While checking for inconsistencies in the atom NE2 of the residue %s, it was found a metal "
                                         "covalently bound to the atom OE1. However, nothing will be done because 'break_metal_bonds' "
                                         "was set to False." % pdb_atm.parent)

                    # If the oxygen is bound to something else not comprised in the previous rules, it is better not to update anything.
                    # Otherwise, fix the NE2.
                    elif atm_obj.matches_smarts("N[CX4]([#6])[OX2][!#1]") is False:
                        self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=2)

                    # Alert for unexpected atoms bound to the oxygen OE1.
                    elif atm_obj.matches_smarts("N[CX4]([#6])[OX2][!#1]"):
                        logger.debug("While checking for inconsistencies in the atom NE2 of the residue %s, it was found an unexpected "
                                     "atom covalently bound to the oxygen OE1. So, NE2 will not be amended." % pdb_atm.parent)

            # Atom: NZ
            # Sanity check for LYS:NZ with invalid bond types.
            elif (pdb_atm.name == "NZ" and atm_obj.get_neighbors_number(True) == 1 and
                    (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 1 or atm_obj.get_degree() != 4 or
                        atm_obj.get_h_count() != 3 or atm_obj.is_in_ring() or atm_obj.is_aromatic())):

                self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=1, implicit_h_count=3)

            # Atom: NH1
            # Sanity check for ARG:NH1 with invalid bond types.
            elif (pdb_atm.name == "NH1" and atm_obj.get_neighbors_number(True) == 1 and
                    (atm_obj.get_valence() != 3 or atm_obj.get_charge() != 0 or atm_obj.get_degree() != 3 or
                        atm_obj.get_h_count() != 2 or atm_obj.is_in_ring() or atm_obj.is_aromatic())):

                self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=2)

            # Atom: NH2
            # Sanity check for ARG:NH2 with invalid bond types.
            elif (pdb_atm.name == "NH2" and atm_obj.get_neighbors_number(True) == 1 and
                    (atm_obj.get_valence() != 4 or atm_obj.get_charge() != 1 or atm_obj.get_degree() != 3 or
                        atm_obj.get_h_count() != 2 or atm_obj.is_in_ring() or atm_obj.is_aromatic())):

                self._fix_atom(atm_obj, bond_types=[BondType.DOUBLE], charge=1, implicit_h_count=2)

            # Atom: O, OD1, and OE1
            # Sanity check for OD1/OE1 oxygens with invalid bond types.
            elif pdb_atm.name in ["O", "OD1", "OE1"]:

                # Any O oxygen or ASN/ASP/GLN/GLU OD1/OE1 with metallic coordination perceived as covalent bond.
                #
                #   It identifies two forms:
                #       - the firt consists of the usual form generated by Open Babel where the double bonds between O and C,
                #           OD1 and CG, or OE1 and CD become single bonds, while it adds an additional single bond with a Metal.
                #
                #       - the second form may appear after correcting C/CG/CD first, so the missing double bond would have already been
                #           fixed, but it makes the oxygens to have incorrect valence and degree, so we still need to fix them.
                #           And of course, we still need to remove the bond with the metal.
                #
                # E.g.: 6JWU:A:GLU:42, 6JWU:B:ASP:210, 4FVR:A:ASN:678.
                #
                if (atm_obj.matches_smarts(f"[OX2]({METAL_ATOM})-,=[C;X4,X3]([#6])(-,=[N])") or
                        atm_obj.matches_smarts(f"[OX2]({METAL_ATOM})-,=[C;X4,X3]([#6])(-,=[O])")):
                    # If it is necessary to break covalent bonds with metals.
                    if self.break_metal_bonds:
                        self._remove_metallic_bond(atm_obj)
                        self._fix_atom(atm_obj, bond_types=[BondType.DOUBLE], charge=0, implicit_h_count=0)
                    else:
                        logger.debug("While checking for inconsistencies in the atom %s of the residue %s, it was found a metal "
                                     "covalently bound to it. However, nothing will be done because 'break_metal_bonds' "
                                     "was set to False." % (pdb_atm.name, pdb_atm.parent))

                # If the oxygen is bound to something else not comprised in the previous rules, it is better not to update anything.
                # Otherwise, fix the oxygen.
                #
                # Note that it doesn't check out the situation of the other oxygen (OXT or OD2/OE2 in ASP/GLU), because O/OD1/OE1 will
                # always have a double bond no matter the other oxygen has or not a bond with metals or other atoms.
                #
                # E.g.: 3QQL:A:GLY:11 contain incorrectly perceived bonds.
                #
                elif (atm_obj.get_neighbors_number(True) == 1 and (atm_obj.get_valence() != 2 or atm_obj.get_charge() != 0 or
                                                                   atm_obj.has_only_bond_type(BondType.DOUBLE) is False or
                                                                   atm_obj.get_h_count() != 0 or atm_obj.is_in_ring() or
                                                                   atm_obj.is_aromatic())):

                    self._fix_atom(atm_obj, bond_types=[BondType.DOUBLE], charge=0, implicit_h_count=0)

                # Alert for unexpected atoms bound to the oxygens.
                elif atm_obj.matches_smarts(f"[OX2;$(O([C;X4,X3])[!#1]);!$(O{METAL_ATOM})]-,=[C;X4,X3]"):
                    logger.debug("While checking for inconsistencies in the atom %s of the residue %s, it was found an unexpected "
                                 "atom covalently bound to it. So, %s will not be amended." % (pdb_atm.name, pdb_atm.parent,
                                                                                                 pdb_atm.name))

            # Atom: OXT, OD2, and OE2
            # Sanity check for OXT/OD2/OE2 oxygens with invalid bond types.
            elif pdb_atm.name in ["OXT", "OD2", "OE2"]:

                # ASP/GLU or any OXT with metallic coordination perceived as covalent bond.
                #
                #   It considers monodentate (only OXT/OD2/E2) and bidentate (both oxygens in ASP/GLU or main chain O and OXT)
                #   interactions with metals.
                #
                #   For the second oxygen (O/OD1/OE1), we also consider that it may appear after correcting C/CG/CD first, so the
                #   missing double bond would have already been fixed.
                #
                #   E.g.: 6JWU:B:ASP:210.
                #
                if (atm_obj.matches_smarts(f"[OX2]({METAL_ATOM})[C;X4,X3]([#6])-,=[OX2;$(O{METAL_ATOM})]") or
                        atm_obj.matches_smarts(f"[OX2]({METAL_ATOM})[CX3]([#6])=[OX1]")):
                    # If it is necessary to break covalent bonds with metals.
                    if self.break_metal_bonds:
                        self._remove_metallic_bond(atm_obj)
                        self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=-1, implicit_h_count=0)
                    else:
                        logger.debug("While checking for inconsistencies in the atom %s of the residue %s, it was found a metal "
                                     "covalently bound to it. However, nothing will be done because 'break_metal_bonds' "
                                     "was set to False." % (pdb_atm.name, pdb_atm.parent))

                # If the second oxygen (OD1/OE1) is bound to something else not comprised in the previous rules,
                # then it won't be fixed and, therefore, we should update OD2/OE2 as follows: remove the bond with
                # the metal and then substitute the single bond with CG/CD for a double bond. Note that OD2/OE2
                # will play the role of OD1/OE1 after the update because that oxygen won't be amended.
                #
                # Note, we also consider that it may appear after correcting C/CG/CD first, so the missing double
                # bond would have already been fixed.
                elif atm_obj.matches_smarts(f"[OX2]({METAL_ATOM})-,=[C;X4,X3]([#6])[OX2;$(O([C;X4,X3])[!#1])]"):
                    # If it is necessary to break covalent bonds with metals.
                    if self.break_metal_bonds:
                        self._remove_metallic_bond(atm_obj)
                        self._fix_atom(atm_obj, bond_types=[BondType.DOUBLE], charge=0, implicit_h_count=0)
                    else:
                        logger.debug("While checking for inconsistencies in the atom %s of the residue %s, it was found a metal "
                                     "covalently bound to it. However, nothing will be done because 'break_metal_bonds' "
                                     "was set to False." % (pdb_atm.name, pdb_atm.parent))

                # Capture cases where the second oxygen is bound to a metal atom, but not OXT/OD2/OE2.
                elif atm_obj.matches_smarts(f"[OX1]=[CX3]([#6])[OX2;$(O{METAL_ATOM})]"):
                    # If it is necessary to break covalent bonds with metals.
                    if self.break_metal_bonds:
                        self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=-1, implicit_h_count=0)
                    else:
                        logger.debug("While checking for inconsistencies in the atom %s of the residue %s, it was found a metal "
                                     "covalently bound to the other carboxyl oxygen. However, nothing will be done because "
                                     "'break_metal_bonds' was set to False." % (pdb_atm.name, pdb_atm.parent))

                # It fixes invalid oxygens, excluing cases where OD2/OE2 is bound to an atom not comprised in the
                # previous rules and cases where OD2/OE2 has a double bond with CG/CD and the second oxygen (OD1/OE1)
                # is bound to two atoms, one of which is an unexpected atom not considered in our standardization function
                elif (atm_obj.matches_smarts("[OX2]([!#1])[C;X4,X3]") is False and
                        atm_obj.matches_smarts(f"[OX1]=[CX3]([#6])[OX2;$(O([CX3])[!#1]);!$(O{METAL_ATOM})]") is False and
                        atm_obj.get_neighbors_number(True) == 1 and (atm_obj.get_valence() != 1 or atm_obj.get_charge() != -1 or
                                                                     atm_obj.get_degree() != 1 or atm_obj.get_h_count() != 0 or
                                                                     atm_obj.is_in_ring() or atm_obj.is_aromatic())):
                    self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=-1, implicit_h_count=0)

                # Alert for unexpected atoms bound to the oxygens.
                elif atm_obj.matches_smarts("[OX2]([!#1])[C;X4,X3]"):
                    logger.debug("While checking for inconsistencies in the atom %s of the residue %s, it was found an unexpected "
                                 "atom covalently bound to it. So, %s will not be amended." % (pdb_atm.name, pdb_atm.parent,
                                                                                                 pdb_atm.name))

                # Alert for unexpected atoms bound to the oxygens.
                elif atm_obj.matches_smarts(f"[OX1]=[CX3]([#6])[OX2;$(O([CX3])[!#1]);!$(O{METAL_ATOM})]"):
                    logger.debug("While checking for inconsistencies in the atom %s of the residue %s, it was found an unexpected "
                                 "atom covalently bound to the other carboxyl oxygen. So, %s was treated as a carbonyl oxygen." %
                                 (pdb_atm.name, pdb_atm.parent, pdb_atm.name))

            # Atom: OG, OG1, and OH
            # Sanity check for OG/OG1/OH oxygens with invalid bond types.
            elif pdb_atm.name in ["OG", "OG1", "OH"]:

                # TYR/SER/THR with metallic coordination perceived as covalent bond.
                #
                # Bertini et al. 2007. Biological Inorganic Chemistry: Structure and Reactivity.
                #
                # E.g: 1TFD:A:TYR:188
                #
                if atm_obj.matches_smarts(f"[OX2]{METAL_ATOM}"):
                    # If it is necessary to break covalent bonds with metals.
                    if self.break_metal_bonds:
                        self._remove_metallic_bond(atm_obj)
                        self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=1)
                    else:
                        logger.debug("While checking for inconsistencies in the atom %s of the residue %s, it was found a metal "
                                     "covalently bound to it. However, nothing will be done because 'break_metal_bonds' "
                                     "was set to False." % (pdb_atm.name, pdb_atm.parent))

                # If the oxygen is bound to something else not comprised in the previous rule, it is better not to
                # update anything. Otherwise, fix OG/OG1/OH.
                elif atm_obj.get_neighbors_number(True) == 1 and (atm_obj.get_valence() != 2 or atm_obj.get_charge() != 0 or
                                                                  atm_obj.has_only_bond_type(BondType.SINGLE) is False or
                                                                  atm_obj.get_h_count() != 1 or atm_obj.is_in_ring() or
                                                                  atm_obj.is_aromatic()):

                    self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=1)

                # Alert for unexpected atoms bound to the hydroxyl oxygens.
                elif atm_obj.matches_smarts(f"[OX2;$(O([#6])[!#1]);!$(O{METAL_ATOM})]"):
                    logger.debug("While checking for inconsistencies in the atom %s of the residue %s, it was found an unexpected "
                                 "atom covalently bound to it. So, %s will not be amended." % (pdb_atm.name, pdb_atm.parent,
                                                                                                 pdb_atm.name))

            # Atom: SD
            # Sanity check for MET:SD with invalid bond types.
            #
            # Although MET can be involved in the coordination of metals, due to its two single bonds, it will never be incorrectly
            # perceived as being covalently bound to metals. Therefore, we do not need to check for it.
            #
            elif (pdb_atm.name == "SD" and atm_obj.get_neighbors_number(True) == 2 and
                    (atm_obj.get_valence() != 2 or atm_obj.get_charge() != 0 or atm_obj.has_only_bond_type(BondType.SINGLE) is False or
                        atm_obj.get_h_count() != 0 or atm_obj.is_in_ring() or atm_obj.is_aromatic())):

                self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=0)

            # Atom: SG
            # Sanity check for CYS:SG with invalid bond types.
            elif pdb_atm.name == "SG":

                # CYS with metallic coordination perceived as covalent bond.
                #
                # Harding et al. 2010. Metals in protein structures: a review of their principal features.
                # Bertini et al. 2007. Biological Inorganic Chemistry: Structure and Reactivity.
                #
                # E.g: 1IUZ:A:CYS:84
                #
                if atm_obj.matches_smarts(f"[SX2]{METAL_ATOM}"):
                    # If it is necessary to break covalent bonds with metals.
                    if self.break_metal_bonds:
                        self._remove_metallic_bond(atm_obj)
                        self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=-1, implicit_h_count=0)
                    else:
                        logger.debug("While checking for inconsistencies in the atom SG of the residue %s, it was found a metal "
                                     "covalently bound to it. However, nothing will be done because 'break_metal_bonds' "
                                     "was set to False." % pdb_atm.parent)

                # If the sulfur is bound to something else not comprised in the previous rule, it is better not to
                # update anything. Otherwise, fix SG.
                elif (atm_obj.matches_smarts("[SX2](C)[!#1]") is False and
                        (atm_obj.get_valence() != 2 or atm_obj.get_charge() != 0 or atm_obj.get_degree() != 2 or
                         atm_obj.has_only_bond_type(BondType.SINGLE) is False or atm_obj.is_in_ring() or
                         atm_obj.is_aromatic())):

                    # Get the number of implicit hydrogen count based on the bonds SG has, i.e., if it is the usual cysteine, then
                    # it has 1 implicit hydrogen, but if it establishes a disulfide bond, then it has no hydrogen.
                    implicit_h_count = 0 if atm_obj.get_neighbors_number(True) == 2 else 1

                    self._fix_atom(atm_obj, bond_types=[BondType.SINGLE], charge=0, implicit_h_count=implicit_h_count)

                # Alert for unexpected atoms bound to the hydroxyl oxygens.
                elif atm_obj.matches_smarts("[SX2](C)[!#1]"):
                    logger.debug("While checking for inconsistencies in the atom SG of the residue %s, it was found an unexpected "
                                 "atom covalently bound to it. So, SG will not be amended." % pdb_atm.parent)

        # Only remove metals with no bond to any residue.
        for metal_obj in self.found_metals.values():
            if len(metal_obj.get_bonds()) == 0:
                self.removed_atoms.append(metal_obj.get_id())
                atm_obj.parent.unwrap().DeleteAtom(metal_obj)

    def _fix_atom(self, atm_obj, charge=None, remove_explict_h=True, implicit_h_count=None,
                  bond_types=None, is_aromatic=False, in_ring=False):

        # Remove all current explicit hydrogens.
        if remove_explict_h:
            self._remove_explict_hydrogens(atm_obj)

        # Convert bond types.
        if bond_types is not None:
            self._fix_bonds(atm_obj, bond_types)

        # All aromatic atoms will have their in_ring property set to True.
        if is_aromatic:
            in_ring = True

        # Set if atom belongs to a ring or not.
        atm_obj.set_in_ring(in_ring)

        # Set if an atom is aromatic or not.
        atm_obj.set_as_aromatic(is_aromatic)

        if implicit_h_count is not None:
            atm_obj.unwrap().SetImplicitHCount(implicit_h_count)

        # Set charge.
        if charge is not None:
            atm_obj.set_charge(charge)

    def _remove_explict_hydrogens(self, atm_obj):
        delete_hs = []
        for b in atm_obj.get_bonds():
            if b.get_partner_atom(atm_obj).get_symbol() == "H":
                delete_hs.append(b.get_partner_atom(atm_obj))

        for hs_obj in delete_hs:
            atm_obj.parent.unwrap().DeleteAtom(hs_obj.unwrap())

    def _fix_bonds(self, atm_obj, bond_types):
        # bond_types becomes an empty list if bond_types is None.
        bond_types = bond_types or []

        # If bond_types is an empty list, do nothing.
        if len(bond_types) == 0:
            return

        # All bonds will be converted to the same type.
        elif len(bond_types) == 1 and isinstance(bond_types[0], BondType):
            new_bond_type = bond_types[0]
            for bond_obj in atm_obj.get_bonds():
                bond_obj.set_bond_type(new_bond_type)

        # It expects a list of tuples where the first element is the bond to be updated, the second element should be the new type,
        # and the third optional element is a flag to indicate if the bond is aromatic or not.
        else:
            for bond_info in bond_types:
                bond_obj, new_bond_type = bond_info[0:2]
                # Update bond type
                bond_obj.set_bond_type(new_bond_type)

                # If the tuple contains three elements, it must be a boolean to define if it is an aromatic bond or not.
                if len(bond_info) == 3:
                    is_aromatic = bond_info[2]
                    bond_obj.set_as_aromatic(is_aromatic)

    def _remove_metallic_bond(self, atm_obj):

        bonds_to_remove = []
        for bond_obj in atm_obj.get_bonds():
            partner_obj = bond_obj.get_partner_atom(atm_obj)

            if partner_obj.get_symbol() in METALS:
                bonds_to_remove.append((partner_obj, bond_obj))

        # Remove metallic bond and mark metal to remotion.
        for metal_obj, bond_obj in bonds_to_remove:

            self.found_metals[metal_obj.get_id()] = metal_obj

            atm_obj.parent.unwrap().DeleteBond(bond_obj.unwrap())
