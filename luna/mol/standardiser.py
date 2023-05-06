from operator import xor
from enum import Enum, auto
from collections import defaultdict, Counter

from Bio.PDB.Polypeptide import is_aa
from luna.mol.precomp_data import DefaultResidueData
from luna.wrappers.base import MolWrapper, BondType, OBBondType
from luna.MyBio.neighbors import get_residue_neighbors

from openbabel.openbabel import OBSmartsPattern

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

    # Aspartate with double bond on OD1.
    ASP = auto()
    # Same as ASP.
    ASP_OD1 = auto()
    # Aspartate with double bond on OD2.
    ASP_OD2 = auto()

    # Glutamate with double bond on OE1.
    GLU = auto()
    # Same as GLU.
    GLU_OE1 = auto()
    # Glutamate with double bond on OE2.
    GLU_OE2 = auto()

    # HIS: same as HIE.
    HIS = auto()
    # HID as in Amber: histidine with hydrogen on the delta nitrogen.
    HID = auto()
    # HIE as in Amber: histidine with hydrogen on the epsilon nitrogen.
    HIE = auto()
    # HIP as in Amber: histidine with hydrogens on both nitrogens, i.e.,
    # this is positively charged.
    HIP = auto()
    # Negatively charged HIS: both ND1 and NE2 are deprotonated.
    # This is usually used for HIS coordinating two metals.
    HIS_D = auto()

    # CYS as in Amber: neutral cysteine
    CYS = auto()
    # CYM as in Amber: deprotonated cysteine.
    CYM = auto()

    # LYS as in Amber: protonated lysine
    LYS = auto()
    # LYN as in Amber: neutral lysine.
    LYN = auto()

    # Default SER
    SER = auto()
    # SER with the hydroxyl group deprotonated.
    SER_D = auto()

    # Default THR
    THR = auto()
    # THR with the hydroxyl group deprotonated.
    THR_D = auto()

    # Default TYR
    TYR = auto()
    # TYR with the hydroxyl group deprotonated.
    TYR_D = auto()


class Standardizer:

    def __init__(self,
                 asp_type=ResidueStandard.ASP_OD1,
                 cys_type=ResidueStandard.CYS,
                 glu_type=ResidueStandard.GLU_OE1,
                 his_type=ResidueStandard.HIE,
                 lys_type=ResidueStandard.LYS,
                 res_data=DefaultResidueData()):

        self.asp_type = asp_type
        self.cys_type = cys_type
        self.glu_type = glu_type
        self.his_type = his_type
        self.lys_type = lys_type

        self.res_data = res_data

    def _validate_atoms_list(self, res, atom_names, precomp_data=None):
        # Current list of atoms found in this residue.
        atom_names = set(atom_names)

        precomp_data = precomp_data or self.res_data

        # Set of expected atom names.
        data = precomp_data.get(res.resname, {})
        props = data.get("atoms", {})
        expected_atms = set([a for a in props if "," not in a])

        if len(expected_atms) == 0:
            return

        # Identify missing atoms, except OXT.
        missing_atoms = [a for a in expected_atms - atom_names
                         if a != "OXT"]

        # Warn about the missing atoms.
        if missing_atoms:
            logger.warning("The following atoms were not found for "
                           "residue %s: %s."
                           % (res.full_name,
                              ", ".join(sorted(missing_atoms))))

        # Unexpected atoms found in this residue.
        other_atoms = [a for a in atom_names - expected_atms]

        # Warn about the unexpected atom names identified.
        if other_atoms:
            logger.warning("The following unexpected atoms were found "
                           "for residue %s: %s."
                           % (res.full_name,
                              ", ".join(sorted(other_atoms))))

    def _resolve_resname(self, res):

        atom_pairs = self._comps[res]

        # Check if ASP is bound to metals.
        if res.resname == "ASP":
            try:
                OD2_atm_obj = [atm_obj for atm_obj, pdb_atm in atom_pairs
                               if pdb_atm.name == "OD2"][0]
            except IndexError:
                # In case OD2 cannot be recovered, return the
                # ASP type defined by the user.
                return self.asp_type.name

            # OD1 -- M.
            smarts1 = f"[OX1]=[CX3]([#6])[OX2;$(O{METAL_ATOM})]"
            # OD2 -- M.
            smarts2 = f"[OX2]({METAL_ATOM})[CX3]([#6])=[OX1]"

            coords = self.metals_coord.get(res, {})
            O_atms = [a.name for a in coords if a.name in ["OD1", "OD2"]]

            # If OD1 is bound to a metal, move the double bond to OD2.
            if (OD2_atm_obj.matches_smarts(smarts1)
                    or ("OD1" in O_atms and "OD2" not in O_atms)):
                return "ASP_OD2"
            # If OD2 is bound to a metal, move the double bond to OD1.
            elif (OD2_atm_obj.matches_smarts(smarts2)
                    or ("OD1" not in O_atms and "OD2" in O_atms)):
                return "ASP_OD1"
            # Otherwise, return the ASP type defined by the user.
            return self.asp_type.name

        # Check if CYS is bound to metals.
        elif res.resname == "CYS":
            try:
                SG_atm_obj = [atm_obj for atm_obj, pdb_atm in atom_pairs
                              if pdb_atm.name == "SG"][0]
            except IndexError:
                # In case SG cannot be recovered, return the
                # CYS type defined by the user.
                return self.cys_type.name

            coords = self.metals_coord.get(res, {})
            trgt_atms = [a.name for a in coords if a.name in ["SG"]]

            # If SG is bound to a metal, deprotonate it.
            if (SG_atm_obj.matches_smarts(f"[S]{METAL_ATOM}")
                    or "SG" in trgt_atms):
                return "CYM"
            # Otherwise, return the CYS type defined by the user.
            return self.cys_type.name

        # Check if GLU is bound to metals.
        elif res.resname == "GLU":
            try:
                OE2_atm_obj = [atm_obj for atm_obj, pdb_atm in atom_pairs
                               if pdb_atm.name == "OE2"][0]
            except IndexError:
                # In case OE2 cannot be recovered, return the
                # GLU type defined by the user.
                return self.glu_type.name

            # OE1 -- M.
            smarts1 = f"[OX1]=[CX3]([#6])[OX2;$(O{METAL_ATOM})]"
            # OE2 -- M.
            smarts2 = f"[OX2]({METAL_ATOM})[CX3]([#6])=[OX1]"

            coords = self.metals_coord.get(res, {})
            O_atms = [a.name for a in coords if a.name in ["OE1", "OE2"]]

            # If OE1 is bound to a metal, move the double bond to OE2.
            if (OE2_atm_obj.matches_smarts(smarts1)
                    or "OE1" in O_atms):
                return "GLU_OE2"
            # If OE2 is bound to a metal, move the double bond to OE1.
            elif (OE2_atm_obj.matches_smarts(smarts2)
                    or "OE2" in O_atms):
                return "GLU_OE1"
            # Otherwise, return the GLU type defined by the user.
            return self.glu_type.name

        elif res.resname == "HIS":
            try:
                ND1_atm_obj = [atm_obj for atm_obj, pdb_atm in atom_pairs
                               if pdb_atm.name == "ND1"][0]
            except IndexError:
                # In case ND1 cannot be recovered, return the
                # HIS type defined by the user.
                return self.his_type.name

            # HIS bound to two different metals:  M1 <--- ND1 ... NE2 ---> M2.
            smarts1 = f"[#7;R]({METAL_ATOM})~[#6;R]~[#7;R]{METAL_ATOM}"
            # ND1 ---> M1
            smarts2 = f"[#7;R]{METAL_ATOM}"
            # NE2 ---> M1
            smarts3 = f"[#7;R]~[#6;R]~[#7;R]{METAL_ATOM}"

            coords = self.metals_coord.get(res, {})
            N_atms = [a.name for a in coords if a.name in ["ND1", "NE2"]]

            # If both Ns are bound to the metal, deprotonate the HIS.
            if (ND1_atm_obj.matches_smarts(smarts1)
                    or len(N_atms) == 2):
                return "HIS_D"
            # If ND1 is bound to the metal, the H is added to NE2.
            elif (ND1_atm_obj.matches_smarts(smarts2)
                    or "ND1" in N_atms):
                return "HIE"
            # If NE2 is bound to the metal, the H is added to ND1.
            elif (ND1_atm_obj.matches_smarts(smarts3)
                    or "NE2" in N_atms):
                return "HID"
            # Otherwise, return the HIS type defined by the user.
            return self.his_type.name

        elif res.resname == "LYS":
            try:
                NZ_atm_obj = [atm_obj for atm_obj, pdb_atm in atom_pairs
                              if pdb_atm.name == "NZ"][0]
            except IndexError:
                # In case NZ cannot be recovered, return the
                # LYS type defined by the user.
                return self.lys_type.name

            coords = self.metals_coord.get(res, {})
            trgt_atms = [a.name for a in coords if a.name in ["NZ"]]

            # If NZ is bound to a metal, deprotonate it.
            if (NZ_atm_obj.matches_smarts(f"[N]{METAL_ATOM}")
                    or "NZ" in trgt_atms):
                return "LYN"
            # Otherwise, return the LYS type defined by the user.
            return self.lys_type.name

        elif res.resname == "TYR":
            try:
                OH_atm_obj = [atm_obj for atm_obj, pdb_atm in atom_pairs
                              if pdb_atm.name == "OH"][0]
            except IndexError:
                # In case OG/OG1/OH cannot be recovered, return the
                # residue type defined by the user.
                return res.resname

            coords = self.metals_coord.get(res, {})
            trgt_atms = [a.name for a in coords
                         if a.name == "OH"]

            # If the oxygen is bound to a metal, deprotonate it.
            if (OH_atm_obj.matches_smarts(f"[O]{METAL_ATOM}")
                    or "OH" in trgt_atms):
                return res.resname + "_D"
            # Otherwise, return the residue type defined by the user.
            return res.resname

        # If any rule applies to this residue, return its current name.
        return res.resname

    def _check_metal_coordination(self, atm_obj, pdb_atm):

        res = pdb_atm.parent
        n_heavy_atoms = \
            atm_obj.get_neighbors_number(only_heavy_atoms=True)

        # If this residue is not bound to any atom, then its valence
        # is correct.
        if n_heavy_atoms == 0:
            return res.resname

        # Metal properties.
        data = self.res_data.get(res.resname, {})
        props = data.get("atoms", {})

        # Coordination number.
        coords = props[pdb_atm.name].get("coords", [])
        coords = [int(i) for i in coords.split(",")]

        # Evaluate if the coordination number is correct.
        if n_heavy_atoms not in coords:
            msg = ("An unexpected coordination number (%d) has been "
                   "found for metal %s.")
            logger.warning(msg % (n_heavy_atoms, pdb_atm.full_name))

    def _get_partner_atm(self, atm_obj, pdb_atm, bond_obj):
        partner_obj = \
            bond_obj.get_partner_atom(atm_obj)
        partner_idx = partner_obj.get_idx()

        # If this atom's partner cannot be identified, it is better
        # not to update this atom.
        if partner_idx not in self._pdb_mapping:
            self._events["ignore"].add(pdb_atm)

            msg = ("The atom %s from %s will not be amended "
                   "because an unidentified atom is covalently bound "
                   "to it." % (pdb_atm.name, pdb_atm.parent.full_name))
            logger.warning(msg)

            return None, None

        partner_atm = self._pdb_mapping[partner_idx]

        return partner_obj, partner_atm

    def _resolve_bonds(self, resname, atm_obj, pdb_atm, precomp_data=None):

        precomp_data = precomp_data or self.res_data

        # All expected residue bonds.
        res_dict = precomp_data.get(resname, {})
        res_bonds = res_dict.get("bonds", {})

        # Add disulfide bridge as a standard bond to CYS.
        if pdb_atm.parent.resname == "CYS":
            res_bonds["SG,SG"] = {"type": "SINGLE", "is_aromatic": False}

        # Bonds found for atom 'atom_obj', identified by their key.
        found_bonds = set()
        # Set of atom pairs.
        found_atm_pairs = set()

        # Loop through the bonds formed by this atom to check if there
        # is any unexpected bond (e.g., metal bonds).
        for bond_obj in atm_obj.get_bonds():

            # Get this atom's partner.
            partner_obj, partner_atm = \
                self._get_partner_atm(atm_obj, pdb_atm, bond_obj)

            if partner_obj is None or partner_atm is None:
                return

            bond_key = ",".join(sorted([pdb_atm.name, partner_atm.name]))

            # Break any bonds with metals.
            if (atm_obj.get_symbol() in METALS
                    or partner_obj.get_symbol() in METALS):
                self._events["break"].add(((atm_obj, pdb_atm),
                                           (partner_obj, partner_atm),
                                           bond_obj))
                continue

            # Break any bonds with water.
            if pdb_atm.parent.is_water() or partner_atm.parent.is_water():
                self._events["break"].add(((atm_obj, pdb_atm),
                                           (partner_obj, partner_atm),
                                           bond_obj))
                continue

            # If an unexpected bond between residue atoms is found, break it.
            if (is_aa(pdb_atm.parent) and is_aa(partner_atm.parent)
                    and bond_key not in res_bonds):

                # In case there is a peptide bond between a residue and a
                # non-standard residue, whose N has a name different than 'N',
                # we can ignore it if the atom C formed an amide with the N.
                if pdb_atm.name == "C" and partner_obj.get_symbol() == "N":
                    # Correct amide group.
                    if atm_obj.matches_smarts("[$([CX3](=O)N)]"):
                        continue

                    # Fix C-N bond from an amide group incorrectly parsed by
                    # Open Babel.
                    elif atm_obj.matches_smarts("[$([CX3](O)=N)]"):
                        bond_type = BondType["SINGLE"]
                        # Add a new modification event for this bond.
                        self._events["modbonds"].add(((atm_obj,
                                                       pdb_atm),
                                                      (partner_obj,
                                                       partner_atm),
                                                      bond_obj, bond_type,
                                                      False))
                        continue

                self._events["break"].add(((atm_obj, pdb_atm),
                                           (partner_obj, partner_atm),
                                           bond_obj))

                if pdb_atm.parent == partner_atm.parent:
                    msg = ("An unexpected bond between atoms %s and %s from "
                           "%s was found. It will break the bond and "
                           "standardize the atoms."
                           % (pdb_atm, partner_atm, pdb_atm.parent))
                else:
                    msg = ("An unexpected bond between atoms %s and %s from "
                           "residues %s and %s was found. It will break the "
                           "bond and standardize the atoms."
                           % (pdb_atm, partner_atm,
                              pdb_atm.parent, partner_atm.parent))
                logger.warning(msg)
                continue

            # If this atom's partner is not a residue, it is better not to
            # update this atom.
            elif ((is_aa(pdb_atm.parent) and not is_aa(partner_atm.parent))
                    or bond_key not in res_bonds):

                self._events["ignore"].add(pdb_atm)

                msg = ("The atom %s from %s will not be amended "
                       "because an unexpected atom (%s) from %s is covalently "
                       "bound to it." % (pdb_atm.name,
                                         pdb_atm.parent.full_name,
                                         partner_atm.name,
                                         partner_atm.parent.full_name))
                logger.warning(msg)
                return

            # Update the sets with the new identified pair.
            found_bonds.add(bond_key)
            found_atm_pairs.add(tuple(sorted([pdb_atm, partner_atm])))

            # Get bond type from the file definition.
            bond_type = BondType[res_bonds[bond_key]["type"]]
            # Get bond aromaticity from the file definition.
            is_aromatic = res_bonds[bond_key]["is_aromatic"]

            # If the current bond is as expected, don't modify it.
            if (bond_type == bond_obj.get_bond_type()
                    and is_aromatic == bond_obj.is_aromatic()):
                continue

            # Add a new modification event for this bond.
            self._events["modbonds"].add(((atm_obj, pdb_atm),
                                          (partner_obj, partner_atm),
                                          bond_obj, bond_type,
                                          is_aromatic))

        # Set of expected bonds.
        expected_bonds = set([b for b in res_bonds
                              for an in b.split(",") if pdb_atm.name == an])

        # Set of missing bonds.
        # It ignores missing disulfide bonds as they may not exist.
        missing_bonds = [b for b in expected_bonds - set(found_bonds)
                         if b != "SG,SG"]

        # Only verify missing bonds for residues.
        if missing_bonds and pdb_atm.parent.is_residue():
            nb_residues = get_residue_neighbors(pdb_atm.parent, verbose=False)

            # Loop through each missing bond and add them in case
            # both atoms exist.
            for bond_key in missing_bonds:
                bond_pair = bond_key.split(",")

                partner_atm_name = (bond_pair[1]
                                    if bond_pair[0] == pdb_atm.name
                                    else bond_pair[0])

                partner_res = pdb_atm.parent
                if pdb_atm.name == "N" and partner_atm_name == "C":
                    # There's nothing to be done because the previous
                    # residue was not selected or is missing.
                    if not nb_residues.has_predecessor():
                        continue
                    partner_res = nb_residues.predecessor

                elif pdb_atm.name == "C" and partner_atm_name == "N":
                    # There's nothing to be done because the next
                    # residue was not selected or is missing.
                    if not nb_residues.has_successor():
                        continue
                    partner_res = nb_residues.successor

                trg_atm_pairs = []
                # If the partner residue exists, search for the partner
                # atom by name.
                if partner_res in self._comps:
                    atom_pairs = self._comps[partner_res]
                    trg_atm_pairs = [(o, p) for o, p in atom_pairs
                                     if (p.parent == partner_res
                                         and p.name == partner_atm_name)]

                # Skip missing bonds related to C- and N-terminal residues.
                if (bond_key in ["C,N", "C,OXT"]
                        and len(trg_atm_pairs) == 0):
                    continue

                # If the partner atom could not be identified, skip it.
                elif len(trg_atm_pairs) == 0:
                    self._events["ignore"].add(pdb_atm)

                    msg = ("The bond '%s' was not found for residue %s. "
                           "However, it cannot be added because the atom %s "
                           "seems not to exist."
                           % (bond_key, pdb_atm.parent.full_name,
                              partner_atm_name))
                    logger.warning(msg)
                    continue

                # Partner atom's object
                partner_obj, partner_atm = trg_atm_pairs[0]
                # Bond type
                bond_name = res_bonds[bond_key]["type"]
                bond_type = OBBondType[bond_name]
                # Bond aromaticity
                is_aromatic = res_bonds[bond_key]["is_aromatic"]

                # Create a new bond modification event.
                self._events["modbonds"].add(((atm_obj, pdb_atm),
                                              (partner_obj, partner_atm),
                                              None, bond_type,
                                              is_aromatic))

    def _resolve_carbox_invariants(self, trgt_atm_obj, trgt_pdb_atm):

        smarts1 = "[OX1]=[CX3][OX1H0-,OX2H1]"
        smarts2 = "[OX1H0-,OX2H1][CX3]=[OX1]"
        if (not trgt_atm_obj.matches_smarts(smarts1)
                and not trgt_atm_obj.matches_smarts(smarts2)):
            return

        bond_obj = trgt_atm_obj.get_bonds()[0]

        C_atm_obj = bond_obj.get_partner_atom(trgt_atm_obj)

        other_O_atm_obj = None
        other_O_pdb_atm = None

        for C_bond_obj in C_atm_obj.get_bonds():
            C_partner_obj = \
                C_bond_obj.get_partner_atom(C_atm_obj)

            if (C_partner_obj.get_symbol() == "O"
                    and C_partner_obj.get_id() != trgt_atm_obj.get_id()):

                other_O_atm_obj = C_partner_obj
                other_O_pdb_atm = self._pdb_mapping[C_partner_obj.get_idx()]

        # If atomic data from the other carbox O were not found, return None.
        if other_O_atm_obj is None or other_O_pdb_atm is None:
            return

        coords = self.metals_coord.get(trgt_pdb_atm.parent, {})

        # M1 <--- O-C=O ---> M2
        #
        #   if both Os coordinate a metal, it's not necessary to
        #   update their invariants.
        #
        if trgt_pdb_atm in coords and other_O_pdb_atm in coords:
            return

        #     V
        # O-C=O ---> M
        #
        #   if the O with the double bond coordinates a metal,
        #   update its invariant.
        if (trgt_atm_obj.matches_smarts(smarts1)
                and trgt_pdb_atm in coords):

            # Fix the bond: double to single bond.
            C_pdb_atm = self._pdb_mapping[C_atm_obj.get_idx()]
            self._events["modbonds"].add(((trgt_atm_obj, trgt_pdb_atm),
                                          (C_atm_obj, C_pdb_atm),
                                          bond_obj, BondType["SINGLE"],
                                          False))
            return [1, 1, 8, 16, -1, 0, 0]

        # V
        # O-C=O ---> M
        #
        #   if the O with the single bond has a partner O coordinating
        #   a metal, update its invariant.
        elif (trgt_atm_obj.matches_smarts(smarts2)
                and other_O_pdb_atm in coords):

            # Fix the bond: single to double bond.
            C_pdb_atm = self._pdb_mapping[C_atm_obj.get_idx()]
            self._events["modbonds"].add(((trgt_atm_obj, trgt_pdb_atm),
                                          (C_atm_obj, C_pdb_atm),
                                          bond_obj, BondType["DOUBLE"],
                                          False))
            return [1, 2, 8, 16, 0, 0, 0]

    def _resolve_imidazole_invariants(self, trgt_atm_obj, trgt_pdb_atm,
                                      precomp_data=None):

        # Imidazole-like, but not tetrazoles.
        strict_smarts_rule = ("[$(nc[n;H1]),$([n;H1]cn);R1r5;"
                              "!$([nR1r5;$(n:n:n:n:c),$(n:n:n:c:n)])]")

        generic_smarts_rule = ("[$([#7]~[#6]~[#7]);R1r5;"
                               "!$([nR1r5;$(n:n:n:n:c),$(n:n:n:c:n)])]")

        if not trgt_atm_obj.matches_smarts(generic_smarts_rule):
            return

        C_atm_obj = None
        C_bond_obj = None
        for bond_obj in trgt_atm_obj.get_bonds():
            partner_obj = bond_obj.get_partner_atom(trgt_atm_obj)

            if partner_obj.matches_smarts("[$([#6](~[#7])~[#7])]"):
                C_atm_obj = partner_obj
                C_bond_obj = bond_obj
                break

        other_N_pdb_atm = None

        for bond_obj in C_atm_obj.get_bonds():
            C_partner_obj = \
                bond_obj.get_partner_atom(C_atm_obj)

            if (C_partner_obj.get_symbol() == "N"
                    and C_partner_obj.get_id() != trgt_atm_obj.get_id()):

                other_N_pdb_atm = self._pdb_mapping[C_partner_obj.get_idx()]

        is_N_single = False
        is_N_double = False
        if trgt_atm_obj.matches_smarts(strict_smarts_rule):
            if trgt_atm_obj.matches_smarts("[n+0X3H1]"):
                is_N_single = True
            elif trgt_atm_obj.matches_smarts("[n+1X3H1,n+0X2H0]"):
                is_N_double = True

        elif precomp_data is not None:
            res_dict = precomp_data.get(trgt_pdb_atm.parent.resname, {})
            bonds_count = Counter([res_dict["bonds"][b]["type"]
                                   for b in res_dict["bonds"]
                                   for a in b.split(",")
                                   if a == trgt_pdb_atm.name])

            if bonds_count["SINGLE"] == 2:
                is_N_single = True
            elif bonds_count["SINGLE"] == 1 and bonds_count["DOUBLE"] == 1:
                is_N_double = True

        # Both False or both True
        if is_N_single == is_N_double:
            return

        coords = self.metals_coord.get(trgt_pdb_atm.parent, {})

        # M1 <--- n1 ... n2 ---> M2.
        #
        #   if both Ns coordinate a metal, remove the Hs from both Ns.
        if trgt_pdb_atm in coords and other_N_pdb_atm in coords:
            if is_N_single:
                return [2, 2, 7, 14, -1, 0, 1]
            return [2, 3, 7, 14, 0, 0, 1]

        #       v
        #  n1-c=n2 ---> M1.
        #
        #   if the N with the double bond coordinates a metal and have a
        #       H, fix its charge.
        if trgt_pdb_atm in coords and is_N_double:
            return [2, 3, 7, 14, 0, 0, 1]

        #  v
        #  n1-c=n2 ---> M1.
        #
        #   if the N with the single bond has a partner N coordinating
        #       a metal, fix it if necessary.
        if other_N_pdb_atm in coords and is_N_single:
            return [2, 2, 7, 14, 0, 1, 1]

        #       v
        #  n2=c-n1 ---> M1.
        #
        #   if the N with the single bond coordinates a metal, move the
        # double bond to it.
        if trgt_pdb_atm in coords and is_N_single:
            # Fix the bond: single to double bond.
            C_pdb_atm = self._pdb_mapping[C_atm_obj.get_idx()]
            self._events["modbonds"].add(((trgt_atm_obj, trgt_pdb_atm),
                                          (C_atm_obj, C_pdb_atm),
                                          C_bond_obj, BondType["DOUBLE"],
                                          True))
            return [2, 3, 7, 14, 0, 0, 1]

        #  v
        #  n2=c-n1 ---> M1.
        #
        #   if the N with the double bond has a partner N coordinating
        #       a metal, move the single bond to it.
        elif other_N_pdb_atm in coords and is_N_double:

            # Fix the bond: double to single bond.
            C_pdb_atm = self._pdb_mapping[C_atm_obj.get_idx()]
            self._events["modbonds"].add(((trgt_atm_obj, trgt_pdb_atm),
                                          (C_atm_obj, C_pdb_atm),
                                          C_bond_obj, BondType["SINGLE"],
                                          True))
            return [2, 2, 7, 14, 0, 1, 1]

    def _resolve_tetrazole_invariants(self, trgt_atm_obj, trgt_pdb_atm,
                                      precomp_data=None):

        precomp_data = precomp_data or {}

        # Tetrazole-like.
        strict_smarts_rule = ("[nR1r5;$(n:n:n:n:c),$(n:n:n:c:n)]")
        generic_smarts_rule = ("[#7R1r5;$([#7]~[#7]~[#7]~[#7]~[#6]),"
                               "$([#7]~[#7]~[#7]~[#6]~[#7])]")

        if not trgt_atm_obj.matches_smarts(generic_smarts_rule):
            return

        # The next lines will verify if the current N contains
        # a double bond or only single bonds.
        is_N_single = False
        is_N_double = False
        if trgt_atm_obj.matches_smarts(strict_smarts_rule):
            if trgt_atm_obj.matches_smarts("[n;X3H1+0,X2H0-1]"):
                is_N_single = True
            elif trgt_atm_obj.matches_smarts("[nX2H0+0]"):
                is_N_double = True
        elif precomp_data is not None:
            res_dict = precomp_data.get(trgt_pdb_atm.parent.resname, {})
            bonds_count = Counter([res_dict["bonds"][b]["type"]
                                   for b in res_dict["bonds"]
                                   for a in b.split(",")
                                   if a == trgt_pdb_atm.name])

            if bonds_count["SINGLE"] == 2:
                is_N_single = True
            elif bonds_count["SINGLE"] == 1 and bonds_count["DOUBLE"] == 1:
                is_N_double = True

        # Error: return if both False or both True.
        if is_N_single == is_N_double:
            return

        # In the next line, it will capture the tetrazole group that contains
        # the current atom.
        ob_smart = OBSmartsPattern()
        ob_smart.Init("[#6]1~[#7]~[#7]~[#7]~[#7H,#7-1]1")
        ob_smart.Match(trgt_atm_obj.parent.unwrap())

        atoms = []
        matches = [x for x in ob_smart.GetMapList()]
        if matches:
            for match in matches:
                trgt_grp_found = False
                for idx in match:
                    pdb_atm = self._pdb_mapping[idx]
                    atoms.append(pdb_atm)

                    if pdb_atm == trgt_pdb_atm:
                        trgt_grp_found = True

                if trgt_grp_found:
                    break

        # Error: if no atom group has been found by the previous
        # SMARTS definition.
        if not atoms:
            return

        # Count how many metals this tetrazole coordinates.
        n_coords = len([a for a in atoms if a.has_metal_coordination()])

        # If the tetrazole coordinates only one metal, move the negative
        # charge to the N coordinating the metal.
        if n_coords == 1:
            # If the current atom is the one coordinating the metal.
            if trgt_pdb_atm.has_metal_coordination():

                # If the current atom contains a double bond, it will be
                # necessary to standardize all bonds in the tetrazole.
                if is_N_double:
                    # Atoms in the opposite end to the current atom.
                    opp_atm_obj, opp_pdb_atm = None, None
                    # Set of imediate neighbors of the current atom.
                    first_bonded_atms = set()
                    for bond_obj1 in trgt_atm_obj.get_bonds():
                        partner_obj1 = \
                            bond_obj1.get_partner_atom(trgt_atm_obj)
                        pdb_atm1 = self._pdb_mapping[partner_obj1.get_idx()]

                        mod_bond = ((trgt_atm_obj, trgt_pdb_atm),
                                    (partner_obj1, pdb_atm1),
                                    bond_obj1,
                                    BondType["SINGLE"], True)
                        self._events["modbonds"].add(mod_bond)

                        first_bonded_atms.add(pdb_atm1)

                        for bond_obj2 in partner_obj1.get_bonds():
                            partner_obj2 = \
                                bond_obj2.get_partner_atom(partner_obj1)
                            pdb_atm2 = \
                                self._pdb_mapping[partner_obj2.get_idx()]

                            if pdb_atm2 == trgt_pdb_atm:
                                continue
                            if pdb_atm2 not in atoms:
                                continue

                            if opp_atm_obj is None:
                                opp_atm_obj = partner_obj2
                                opp_pdb_atm = pdb_atm2

                            mod_bond = ((partner_obj1, pdb_atm1),
                                        (partner_obj2, pdb_atm2),
                                        bond_obj2,
                                        BondType["DOUBLE"], True)
                            self._events["modbonds"].add(mod_bond)

                    # Standardize the bonds in the opposite end.
                    for bond_obj in opp_atm_obj.get_bonds():
                        partner_obj = \
                            bond_obj.get_partner_atom(opp_atm_obj)
                        pdb_atm = self._pdb_mapping[partner_obj.get_idx()]

                        if pdb_atm in first_bonded_atms:
                            continue

                        mod_bond = ((opp_atm_obj, opp_pdb_atm),
                                    (partner_obj, pdb_atm),
                                    bond_obj,
                                    BondType["SINGLE"], True)
                        self._events["modbonds"].add(mod_bond)
                return [2, 2, 7, 14, -1, 0, 1]
            return [2, 3, 7, 14, 0, 0, 1]

        # If the tetrazole coordinates more than 1 metal, standardize it
        # according to the PDB.
        else:
            self._resolve_bonds(trgt_pdb_atm.parent.resname,
                                trgt_atm_obj,
                                trgt_pdb_atm,
                                precomp_data)

            if is_N_double:
                return [2, 3, 7, 14, 0, 0, 1]
            return [2, 2, 7, 14, -1, 0, 1]

    def _resolve_sulfur_groups_invariants(self, trgt_atm_obj, trgt_pdb_atm,
                                          precomp_data=None):

        precomp_data = precomp_data or {}
        res_dict = precomp_data.get(trgt_pdb_atm.parent.resname, {})
        res_bonds = res_dict.get("bonds", {})

        if not res_bonds:
            return

        bad_sulf_acid_O_smarts = "[$([OH1]S([OH1])([OH1]))]"
        if not trgt_atm_obj.matches_smarts(bad_sulf_acid_O_smarts):
            return

        bond_obj = trgt_atm_obj.get_bonds()[0]

        S_atm_obj = bond_obj.get_partner_atom(trgt_atm_obj)
        S_pdb_atm = self._pdb_mapping[S_atm_obj.get_idx()]

        bond_key = ",".join(sorted([trgt_pdb_atm.name,
                                    S_pdb_atm.name]))
        bond_dict = res_bonds[bond_key]

        if bond_dict["type"] == "SINGLE":
            return [1, 1, 8, 16, -1, 0, 0]

        if bond_dict["type"] == "DOUBLE":
            return [1, 2, 8, 16, 0, 0, 0]

    def _resolve_invariants(self, resname, atm_obj, pdb_atm,
                            atm_names=None, precomp_data=None):

        if pdb_atm in self._events["ignore"]:
            return

        res = pdb_atm.parent
        atm_names = atm_names or []
        precomp_data = precomp_data or self.res_data

        # Expected invariant.
        invariants = None
        # Define if this atom is aromatic.
        is_aromatic = None

        if res.is_hetatm() or res.is_nucleotide():
            coords = self.metals_coord.get(res, {})

            res_dict = precomp_data.get(resname, {})
            res_props = res_dict.get("props", {})

            # Check if the atom aromaticity is incorrect.
            fix_atm_aromaticity = False
            if (pdb_atm.name in res_props
                    and "AROMATIC" in res_props[pdb_atm.name]["features"]
                    and not atm_obj.is_aromatic()):
                fix_atm_aromaticity = True

            # Currently, only standardization rules for O, N, and S were
            # defined. For that reason, only these atoms pass the following
            # condition or atoms whose aromaticity was incorrectly perceived.
            if (atm_obj.get_symbol() not in ["O", "N", "S"]
                    and not fix_atm_aromaticity):
                return

            # Fix O atoms
            if atm_obj.get_symbol() == "O":

                # Carbonate groups.
                carbonate_smarts = \
                    ("[$([OX1]=[CX3]([OX1H0-,OX2H1])[OX1H0-,OX2H1])"
                     ",$([OX1H0-,OX2H1][CX3]([OX1H0-,OX2H1])=[OX1])]")
                # Carboxylic acid, but not carbonate.
                carbox_smarts = ("[$([OX1]=[CX3][OX1H0-,OX2H1])"
                                 ",$([OX1H0-,OX2H1][CX3]=[OX1])"
                                 ";!$(%s)]" % carbonate_smarts)
                phosphate_smarts = "[$([OX2H1]P([OX2H1])([OX2H1])=[OX1])]"
                # O from sulfuric or sulfonic acid perceived as OH.
                bad_sulf_acid_O_smarts = "[$([OH1]S([OH1])([OH1]))]"

                # Carboxylic group O coordinating a metal.
                if atm_obj.matches_smarts(carbox_smarts):
                    invariants = \
                        self._resolve_carbox_invariants(atm_obj, pdb_atm)

                # Force all OH from phosphates to be deprotonated.
                elif atm_obj.matches_smarts(phosphate_smarts):
                    invariants = [1, 1, 8, 16, -1, 0, 0]

                # Phosphoric, phosphonic, or phosphinic acid O
                # coordinating a metal.
                elif (atm_obj.matches_smarts("[$([OX2H1]P(=[OX1]))]")
                        and pdb_atm in coords):
                    invariants = [1, 1, 8, 16, -1, 0, 0]

                # Sulfuric acid O perceived as OH.
                elif atm_obj.matches_smarts(bad_sulf_acid_O_smarts):
                    invariants = \
                        self._resolve_sulfur_groups_invariants(atm_obj,
                                                               pdb_atm,
                                                               precomp_data)
                # Phenol O coordinating a metal as in TYR.
                elif (atm_obj.matches_smarts("[OH+0]c")
                        and pdb_atm in coords):
                    invariants = [1, 1, 8, 16, -1, 0, 0]

            # Fix N atoms.
            elif atm_obj.get_symbol() == "N":

                # Amines.
                amine_smarts = "[$([NX4+;H3,H2,H1;!$(NC=O)][C])]"

                # Tetrazole.
                #
                # This rule doesn't specify aromaticity because some aromatic
                # rings may be perceived a simple cyclic ring by Open Babel.
                generic_tetrazole_N_smarts = \
                    ("[#7R1r5;$([#7]~[#7]~[#7]~[#7]~[#6]),"
                     "$([#7]~[#7]~[#7]~[#6]~[#7])]")

                # Tetrazole.
                #
                # With aromaticity correctly perceived.
                tetrazole_N_smarts = "[nR1r5;$(n:n:n:n:c),$(n:n:n:c:n)]"

                # Imidazole-like, but not tetrazoles.
                #
                # This rule doesn't specify aromaticity because some aromatic
                # rings may be perceived a simple cyclic ring by Open Babel.
                imidazole_smarts = ("[$([#7]~[#6]~[#7]~[#6]~[#6]);R1r5;"
                                    f"!$({generic_tetrazole_N_smarts})]")

                if atm_obj.matches_smarts(generic_tetrazole_N_smarts):
                    if (atm_obj.matches_smarts(tetrazole_N_smarts)
                            or fix_atm_aromaticity):
                        invariants = \
                            self._resolve_tetrazole_invariants(atm_obj,
                                                               pdb_atm,
                                                               precomp_data)
                    if fix_atm_aromaticity:
                        is_aromatic = True

                # Imidazole N coordinating a metal.
                elif atm_obj.matches_smarts(imidazole_smarts):
                    if (atm_obj.matches_smarts('[$(nc[n;H1]),$([n;H1]cn)]')
                            or fix_atm_aromaticity):

                        invariants = \
                            self._resolve_imidazole_invariants(atm_obj,
                                                               pdb_atm,
                                                               precomp_data)

                    if fix_atm_aromaticity:
                        is_aromatic = True

                # Amine N coordinating a metal.
                elif (atm_obj.matches_smarts(amine_smarts)
                        and pdb_atm in coords):
                    invariants = atm_obj.get_atomic_invariants()
                    invariants[4] = 0
                    invariants[5] = invariants[5] - 1

                # Amide N coordinating a metal.
                #
                # Note that only amide groups as in amino acid backbones
                # are considered by this rule.
                elif (atm_obj.matches_smarts("[$([#7;X3H1]C=O)]")
                        and pdb_atm in coords):
                    invariants = [2, 2, 7, 14, -1, 0, 0]

                # Broken amide from residues (including non-standard ones).
                elif (is_aa(pdb_atm.parent)
                        and atm_obj.matches_smarts("[$(N=[CX3]O)]")):
                    invariants = [2, 2, 7, 14, 0, 1, 0]

                # Sulfonamide N coordinating a metal.
                elif (atm_obj.matches_smarts("[$([N;H1,H2]S(=O)(=O))]")
                        and pdb_atm in coords):
                    invariants = atm_obj.get_atomic_invariants()
                    invariants[4] = -1
                    invariants[5] = invariants[5] - 1

            # Fix S atoms
            elif atm_obj.get_symbol() == "S":
                # S from a sulfuric or sulfonic acid where all Os were
                # perceived as OH.
                bad_sulf_acid_S_smarts = "[$(S([OH1])([OH1])([OH1]))]"
                if atm_obj.matches_smarts(bad_sulf_acid_S_smarts):
                    self._resolve_bonds(resname, atm_obj, pdb_atm,
                                        precomp_data)

                # Thiol S coordinating a metal.
                elif (atm_obj.matches_smarts("[SH+0]")
                        and pdb_atm in coords):
                    invariants = [1, 1, 16, 32, -1, 0, 0]

            # Fix C atoms with incorrect aromaticity.
            elif (atm_obj.get_symbol() == "C" and fix_atm_aromaticity
                    and pdb_atm.name in res_props):

                invariants = atm_obj.get_atomic_invariants()

                # Set the correct degree.
                invariants[1] = (res_props[pdb_atm.name]["degree"]
                                 - res_props[pdb_atm.name]["n_Hs"])
                # Set the proper number of Hs.
                invariants[5] = res_props[pdb_atm.name]["n_Hs"]
                # Set the aromaticity flag ON.
                is_aromatic = True

                # Don't update bonds from tetrazole carbons,
                # because their bonds are updated when
                # _resolve_tetrazole_invariants() are called.
                generic_tetrazole_C_smarts = "[$([#6]~[#7]~[#7]~[#7]~[#7])]"
                if not atm_obj.matches_smarts(generic_tetrazole_C_smarts):
                    self._resolve_bonds(resname, atm_obj, pdb_atm,
                                        precomp_data)

        elif res.is_residue():

            if pdb_atm.name in ["C", "N"]:
                nb_residues = \
                    get_residue_neighbors(res, verbose=False)

                coords = self.metals_coord.get(res, {})
                trgt_atms = [a.name for a in coords if a.name in ["N"]]

                if (pdb_atm.name == "N"
                        and nb_residues.predecessor not in self._comps):
                    invariants = [1, 1, 7, 14, 1, 3, 0]
                    if res.resname == "PRO":
                        invariants = [2, 2, 7, 14, 1, 2, 1]

                elif (pdb_atm.name == "N"
                        and atm_obj.matches_smarts(f"[N]{METAL_ATOM}")
                        or "N" in trgt_atms):
                    invariants = [2, 2, 7, 14, -1, 0, 0]

                elif (pdb_atm.name == "C"
                        and nb_residues.successor not in self._comps):

                    n_heavy_atms = \
                        atm_obj.get_neighbors_number(only_heavy_atoms=True)

                    # The following rule accounts for amide bonds between
                    # residues and non-standard residues. That's because the
                    # function get_residue_neighbors() only considers peptide
                    # bonds formed between an atom 'C' and an atom 'N'.
                    # That means the function would ignore a non-standard
                    # residue forming the peptide if its N has a different
                    # name (e.g., N3 as in SNN).
                    amide_like = "[$([CX3](=,-O)=,-N)]"
                    if (nb_residues.successor is None and n_heavy_atms == 3
                            and atm_obj.matches_smarts(amide_like)):
                        invariants = [3, 4, 6, 12, 0, 0, 0]

                    # If it's a C-terminal residue and OXT doesn't exist,
                    # then C has one Hydrogen.
                    elif "OXT" not in atm_names or n_heavy_atms != 3:
                        invariants = [2, 3, 6, 12, 0, 1, 0]

            # CYS:SG in disulfide bridges.
            elif pdb_atm.name == "SG" and atm_obj.matches_smarts("[$(S(C)S)]"):
                invariants = [2, 2, 16, 32, 0, 0, 0]

        # Skip any hetatm and nucleotide if 'invariants' is None.
        #
        # Residues/Water/Metal can pass because their invariants are
        # predefined in the default definition file.
        if (res.is_hetatm() or res.is_nucleotide()) and invariants is None:
            return

        # Load the invariants from the definition file.
        if invariants is None:
            props = precomp_data[resname]["atoms"]
            invariants = props[pdb_atm.name]["invariants"]
            invariants = [int(i) for i in invariants.split(",")]
            is_aromatic = props[pdb_atm.name]["is_aromatic"]

        # Add a new event to modify the atom in case the invariants
        # don't match.
        if (atm_obj.get_atomic_invariants() != invariants
                or is_aromatic is not None):
            self._events["modatms"].add(((atm_obj, pdb_atm),
                                         tuple(invariants),
                                         is_aromatic))

    def _fix_bonds(self, bond_types):
        # bond_types becomes an empty list if bond_types is None.
        bond_types = bond_types or []

        # If bond_types is an empty list, do nothing.
        if len(bond_types) == 0:
            return

        # It expects a list of tuples where the first element is the bond to be
        # updated, the second element should be the new type, and the third
        # optional element is a flag to indicate if the bond is aromatic
        # or not.
        for bond_info in bond_types:
            bond_obj, new_bond_type = bond_info[0:2]
            # Update bond type
            bond_obj.set_bond_type(new_bond_type)

            # If the tuple contains three elements, it must be a boolean to
            # define if it is an aromatic bond or not.
            if len(bond_info) == 3:
                is_aromatic = bond_info[2]
                bond_obj.set_as_aromatic(is_aromatic)

    def _add_bonds(self, mol_obj, new_bonds):
        # new_bonds becomes an empty list if new_bonds is None.
        new_bonds = new_bonds or []

        # If new_bonds is an empty list, do nothing.
        if len(new_bonds) == 0:
            return

        obmol = mol_obj
        if isinstance(mol_obj, MolWrapper):
            mol_obj = mol_obj.unwrap()

        # It expects a list of tuples where the first element is the bond to be
        # updated, the second element should be the new type, and the third
        # optional element is a flag to indicate if the bond is aromatic
        # or not.
        for bond_info in new_bonds:
            begin_idx, end_idx = bond_info[0:2]
            bond_type = bond_info[2]

            # If 'flags' is:
            #   - 0: the bond is not aromatic.
            #   - 2: the bond is aromatic.
            flags = 0

            # If the tuple contains four elements, it must be a boolean to
            # define if it is an aromatic bond or not.
            if len(bond_info) == 4:
                flags = 2 if bond_info[3] is True else 0

            obmol.AddBond(begin_idx, end_idx,
                          bond_type, flags)

    def _remove_explicit_hydrogens(self, atm_obj):
        delete_hs = []
        for b in atm_obj.get_bonds():
            if b.get_partner_atom(atm_obj).get_symbol() == "H":
                delete_hs.append(b.get_partner_atom(atm_obj))

        for hs_obj in delete_hs:
            atm_obj.parent.unwrap().DeleteAtom(hs_obj.unwrap())

    def _fix_atom(self, atm_obj, charge=None, remove_explict_h=True,
                  implicit_h_count=None, is_aromatic=None,
                  in_ring=None):

        # Remove all current explicit hydrogens.
        if remove_explict_h:
            self._remove_explicit_hydrogens(atm_obj)

        # All aromatic atoms will have their in_ring property set to True.
        if is_aromatic:
            in_ring = True

        # Set if atom belongs to a ring or not.
        if in_ring is not None:
            atm_obj.set_in_ring(in_ring)

        # Set if an atom is aromatic or not.
        if is_aromatic is not None:
            atm_obj.set_as_aromatic(is_aromatic)

        if implicit_h_count is not None:
            atm_obj.unwrap().SetImplicitHCount(implicit_h_count)

        # Set charge.
        if charge is not None:
            atm_obj.set_charge(charge)

    def _apply_modifications(self):

        processed_bonds = set()
        for data in self._events["break"]:
            atm_obj, pdb_atm = data[0]
            partner_obj, partner_atm = data[1]
            bond_to_break = data[2]

            bond_key = tuple(sorted([pdb_atm, partner_atm]))

            if bond_key in processed_bonds:
                continue
            if pdb_atm in self._events["ignore"]:
                continue

            if xor(pdb_atm.parent.is_metal(),
                   partner_atm.parent.is_metal()):

                if pdb_atm.parent.is_metal():
                    key, idx = partner_atm, partner_obj.get_idx()
                    val = pdb_atm.parent
                else:
                    key, idx = pdb_atm, atm_obj.get_idx()
                    val = partner_atm.parent

                coords_dict = self.metals_coord[key.parent][key]
                coords_dict["atm_idx"] = idx
                coords_dict["metals"].add(val)

            obmol = atm_obj.parent.unwrap()
            obmol.DeleteBond(bond_to_break.unwrap())

            processed_bonds.add(bond_key)

        for data in self._events["modbonds"]:
            (atm_obj, pdb_atm) = data[0]
            (partner_obj, partner_atm) = data[1]

            bond_to_mod, new_bond_type, is_aromatic = data[2:]

            if (pdb_atm in self._events["ignore"]
                    or partner_atm in self._events["ignore"]):
                continue

            if bond_to_mod is not None:
                self._fix_bonds([(bond_to_mod, new_bond_type, is_aromatic)])
            else:
                obmol = atm_obj.parent
                new_bond = [(atm_obj.get_idx(), partner_obj.get_idx(),
                             new_bond_type.value, is_aromatic)]
                self._add_bonds(obmol, new_bond)

        for data in self._events["modatms"]:
            (atm_obj, pdb_atm) = data[0]
            invariants = data[1]
            is_aromatic = data[2]

            if pdb_atm in self._events["ignore"]:
                continue

            charge, implicit_h_count, in_ring = None, None, None
            if invariants is not None:
                charge = invariants[4]
                implicit_h_count = invariants[5]
                in_ring = bool(invariants[6])

            self._fix_atom(atm_obj,
                           charge=charge,
                           implicit_h_count=implicit_h_count,
                           is_aromatic=is_aromatic,
                           in_ring=in_ring)

    def _get_precomp_data_from_cif_file(self, res):
        # TODO
        return {}

    def _check_mol(self, res):
        atom_pairs = self._comps[res]
        atm_names = [pdb_atm.name
                     for atm_obj, pdb_atm in atom_pairs]

        precomp_data = None
        if res.is_hetatm() or res.is_nucleotide():
            precomp_data = self._get_precomp_data_from_cif_file(res)

        self._validate_atoms_list(res, atm_names, precomp_data)

        resname = res.resname
        if res.is_residue():
            resname = self._resolve_resname(res)
        elif res.is_water():
            resname = "HOH"

        for atm_obj, pdb_atm in atom_pairs:
            if res.is_metal():
                self._check_metal_coordination(atm_obj, pdb_atm)

            if not res.is_hetatm() and not res.is_nucleotide():
                self._resolve_bonds(resname, atm_obj, pdb_atm)

            self._resolve_invariants(resname, atm_obj, pdb_atm,
                                     atm_names, precomp_data)

    def standardize(self, atom_pairs, metals_coord=None):

        if not metals_coord:
            custom_dict = lambda: {"atm_idx": None, "metals": set()}
            metals_coord = defaultdict(lambda: defaultdict(custom_dict))
        self.metals_coord = metals_coord

        self._comps = defaultdict(list)
        self._events = defaultdict(set)
        self._pdb_mapping = {}

        for atm_obj, pdb_atm in atom_pairs:
            self._pdb_mapping[atm_obj.get_idx()] = pdb_atm
            self._comps[pdb_atm.parent].append((atm_obj, pdb_atm))

            coords = self.metals_coord.get(pdb_atm.parent, {})
            if pdb_atm in coords:
                coords[pdb_atm]["atm_idx"] = atm_obj.get_idx()

        for res in self._comps:
            self._check_mol(res)

        self._apply_modifications()
