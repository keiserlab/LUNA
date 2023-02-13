from operator import xor
from enum import Enum, auto
from collections import defaultdict

from luna.mol.precomp_data import DefaultResidueData
from luna.wrappers.base import MolWrapper, BondType, OBBondType
from luna.MyBio.neighbors import get_residue_neighbors


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

    def _validate_atoms_list(self, res, atom_names):
        # Current list of atoms found in this residue.
        atom_names = set(atom_names)

        # Set of expected atom names.
        data = self.res_data.get(res.resname, {})
        props = data.get("atoms", {})
        expected_atms = set([a for a in props if "," not in a])

        # Identify missing atoms, except OXT.
        missing_atoms = [a for a in expected_atms - atom_names
                         if a != "OXT"]

        # Warn about the missing atoms.
        if missing_atoms:
            logger.error("The following atoms were not found for "
                         "residue %s: %s."
                         % (res.full_name,
                            ", ".join(sorted(missing_atoms))))

        # Unexpected atoms found in this residue.
        other_atoms = [a for a in atom_names - expected_atms]

        # Warn about the unexpected atom names identified.
        if other_atoms:
            logger.error("The following unexpected atoms were found "
                         "for residue %s: %s."
                         % (res.full_name,
                            ", ".join(sorted(other_atoms))))

    def _resolve_resname(self, res):

        atom_pairs = self.residues[res]

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

            # If OD1 is bound to a metal, move the double bond to OD2.
            if OD2_atm_obj.matches_smarts(smarts1):
                return "ASP_OD2"
            # If OD2 is bound to a metal, move the double bond to OD1.
            elif OD2_atm_obj.matches_smarts(smarts2):
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

            # If SG is bound to a metal, deprotonate it.
            if SG_atm_obj.matches_smarts(f"[S]{METAL_ATOM}"):
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

            # If OE1 is bound to a metal, move the double bond to OE2.
            if OE2_atm_obj.matches_smarts(smarts1):
                return "GLU_OE2"
            # If OE2 is bound to a metal, move the double bond to OE1.
            elif OE2_atm_obj.matches_smarts(smarts2):
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

            # If both Ns are bound to the metal, deprotonate the HIS.
            if ND1_atm_obj.matches_smarts(smarts1):
                return "HIS_D"
            # If ND1 is bound to the metal, the H is added to NE2.
            elif ND1_atm_obj.matches_smarts(smarts2):
                return "HIE"
            # If NE2 is bound to the metal, the H is added to ND1.
            elif ND1_atm_obj.matches_smarts(smarts3):
                return "HID"
            # Otherwise, return the HIS type defined by the user.
            else:
                return self.his_type.name

        elif res.resname == "LYS":
            try:
                NZ_atm_obj = [atm_obj for atm_obj, pdb_atm in atom_pairs
                              if pdb_atm.name == "NZ"][0]
            except IndexError:
                # In case NZ cannot be recovered, return the
                # LYS type defined by the user.
                return self.lys_type.name

            # If NZ is bound to a metal, deprotonate it.
            if NZ_atm_obj.matches_smarts(f"[N]{METAL_ATOM}"):
                return "LYN"
            # Otherwise, return the LYS type defined by the user.
            return self.lys_type.name

        elif res.resname in ["SER", "THR", "TYR"]:
            try:
                OH_atm_obj = [atm_obj for atm_obj, pdb_atm in atom_pairs
                              if pdb_atm.name in ["OG", "OG1", "OH"]][0]
            except IndexError:
                # In case OG/OG1/OH cannot be recovered, return the
                # residue type defined by the user.
                return res.resname

            # If the oxygen is bound to a metal, deprotonate it.
            if OH_atm_obj.matches_smarts(f"[O]{METAL_ATOM}"):
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
            logger.error(msg % (n_heavy_atoms, pdb_atm.full_name))

    def _get_partner_atm(self, atm_obj, pdb_atm, bond_obj):
        partner_obj = \
            bond_obj.get_partner_atom(atm_obj)
        partner_idx = partner_obj.get_idx()

        # If this atom's partner cannot be identified, it is better
        # not to update this atom.
        if partner_idx not in self._pdb_map:
            self._events["ignore"].add(pdb_atm)

            msg = ("Atom %s of residue %s will not be amended "
                   "because an unidentified atom is covalently bound "
                   "to it." % (pdb_atm, pdb_atm.parent))
            logger.error(msg)

            return None, None

        partner_atm = self._pdb_map[partner_idx]

        return partner_obj, partner_atm

    def _resolve_bonds(self, resname, atm_obj, pdb_atm):

        # All expected residue bonds.
        res_dict = self.res_data.get(resname, {})
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
            if pdb_atm.parent.is_water():
                self._events["break"].add(((atm_obj, pdb_atm),
                                           (partner_obj, partner_atm),
                                           bond_obj))
                continue

            # If an unexpected bond between residue atoms
            # is found, break it.
            if (pdb_atm.parent.is_residue()
                    and partner_atm.parent.is_residue(standard=False)
                    and bond_key not in res_bonds):

                self._events["break"].add(((atm_obj, pdb_atm),
                                           (partner_obj, partner_atm),
                                           bond_obj))

                if pdb_atm.parent == partner_atm.parent:
                    msg = ("An unexpected bond between atoms %s and %s from "
                           "%s was found. It will break the bond and "
                           "standardise the atoms."
                           % (pdb_atm, partner_atm, pdb_atm.parent))
                else:
                    msg = ("An unexpected bond between atoms %s and %s from "
                           "residues %s and %s was found. It will break the "
                           "bond and standardise the atoms."
                           % (pdb_atm, partner_atm,
                              pdb_atm.parent, partner_atm.parent))
                logger.error(msg)
                continue

            # If this atom's partner is not a residue, it is better
            # not to update this atom.
            elif ((pdb_atm.parent.is_residue()
                    and not partner_atm.parent.is_residue(standard=False))
                    or bond_key not in res_bonds):

                self._events["ignore"].add(pdb_atm)

                msg = ("Atom %s of molecule %s will not be amended "
                       "because an unexpected atom (%s) from %s is covalently "
                       "bound to it." % (pdb_atm, pdb_atm.parent,
                                         partner_atm, partner_atm.parent))
                logger.error(msg)
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
                if partner_res in self.residues:
                    atom_pairs = self.residues[partner_res]
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
                           % (bond_key, pdb_atm.parent, partner_atm_name))
                    logger.error(msg)
                    continue

                # Partner atom's object
                partner_obj, partner_atm = trg_atm_pairs[0]
                # Bond type
                bond_name = res_bonds[bond_key]["type"]
                bond_type = OBBondType[bond_name].value
                # Bond aromaticity
                is_aromatic = res_bonds[bond_key]["is_aromatic"]

                # Create a new bond modification event.
                self._events["modbonds"].add(((atm_obj, pdb_atm),
                                              (partner_obj, partner_atm),
                                              None, bond_type,
                                              is_aromatic))

    def _resolve_invariants(self, resname, atm_obj, pdb_atm, atm_names=None):

        if pdb_atm in self._events["ignore"]:
            return

        res = pdb_atm.parent
        atm_names = atm_names or []

        # Expected invariant.
        invariants = None
        # Define if this atom is aromatic.
        is_aromatic = None

        if res.is_water():
            smarts1 = f"[O]({METAL_ATOM}){METAL_ATOM}"
            smarts2 = f"[O]{METAL_ATOM}"
            if (pdb_atm.name == "O"
                    and atm_obj.matches_smarts(smarts1)):
                invariants = [0, 0, 8, 16, -2, 0, 0]

            elif (pdb_atm.name == "O"
                    and atm_obj.matches_smarts(smarts2)):
                invariants = [0, 0, 8, 16, -1, 1, 0]

        elif pdb_atm.name in ["C", "N"]:
            nb_residues = get_residue_neighbors(res,
                                                verbose=False)

            if (pdb_atm.name == "N"
                    and nb_residues.predecessor not in self.residues):

                invariants = [1, 1, 7, 14, 1, 3, 0]
                if res.resname == "PRO":
                    invariants = [2, 2, 7, 14, 1, 2, 1]

            elif (pdb_atm.name == "N"
                    and atm_obj.matches_smarts(f"[N]{METAL_ATOM}")):
                invariants = [2, 2, 7, 14, -1, 0, 0]

            elif (pdb_atm.name == "C"
                    and nb_residues.successor not in self.residues):

                n_heavy_atms = \
                    atm_obj.get_neighbors_number(only_heavy_atoms=True)

                # If it's a C-terminal residue and OXT exists,
                # then C has no Hydrogen.
                if "OXT" not in atm_names or n_heavy_atms != 3:
                    invariants = [2, 3, 6, 12, 0, 1, 0]

        # Load the invariants from the definition file.
        if invariants is None:
            props = self.res_data[resname]["atoms"]
            invariants = props[pdb_atm.name]["invariants"]
            invariants = [int(i) for i in invariants.split(",")]
            is_aromatic = props[pdb_atm.name]["is_aromatic"]

        # Add a new event to modify the atom in case the invariants
        # don't match.
        if atm_obj.get_atomic_invariants() != invariants:
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
                    key = (partner_obj.get_idx(), partner_atm)
                    val = atm_obj
                else:
                    key = (atm_obj.get_idx(), pdb_atm)
                    val = partner_obj

                self.bound_to_metals[key].append(val)

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
                             new_bond_type, is_aromatic)]
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

    def _standardize_residue(self, res):
        atom_pairs = self.residues[res]
        atm_names = [pdb_atm.name
                     for atm_obj, pdb_atm in atom_pairs]
        self._validate_atoms_list(res, atm_names)

        resname = self._resolve_resname(res)
        for atm_obj, pdb_atm in atom_pairs:
            self._resolve_bonds(resname, atm_obj, pdb_atm)
            self._resolve_invariants(resname, atm_obj, pdb_atm, atm_names)

    def _standardize_water(self, res):
        atom_pairs = self.residues[res]

        if len(atom_pairs) != 1:
            logger.error("The water molecule %s contains %d "
                         "heavy atom(s), while one is expected."
                         % (res, len(atom_pairs)))
            return

        atm_obj, pdb_atm = atom_pairs[0]

        self._resolve_bonds("HOH", atm_obj, pdb_atm)
        self._resolve_invariants("HOH", atm_obj, pdb_atm)

    def _standardize_metal(self, res):
        atom_pairs = self.residues[res]

        if len(atom_pairs) != 1:
            logger.error("The metal %s contains %d "
                         "heavy atom(s), while one is expected."
                         % (res, len(atom_pairs)))
            return

        atm_obj, pdb_atm = atom_pairs[0]

        resname = self._resolve_metal_name(pdb_atm, atm_obj)

        self._resolve_bonds(res.resname, atm_obj, pdb_atm)
        self._resolve_invariants(resname, atm_obj, pdb_atm)

    def _check_mol(self, res):
        atom_pairs = self.residues[res]
        atm_names = [pdb_atm.name
                     for atm_obj, pdb_atm in atom_pairs]
        self._validate_atoms_list(res, atm_names)

        resname = res.resname
        if res.is_residue():
            resname = self._resolve_resname(res)
        elif res.is_water():
            resname = "HOH"

        for atm_obj, pdb_atm in atom_pairs:
            if res.is_metal():
                self._check_metal_coordination(atm_obj, pdb_atm)

            self._resolve_bonds(resname, atm_obj, pdb_atm)
            self._resolve_invariants(resname, atm_obj, pdb_atm, atm_names)

    def standardize(self, atom_pairs):

        self.bound_to_metals = defaultdict(list)

        self.metals = {}
        self.events = []

        self.residues = defaultdict(list)
        self._pdb_map = {}

        self._events = defaultdict(set)

        for atm_obj, pdb_atm in atom_pairs:
            self._pdb_map[atm_obj.get_idx()] = pdb_atm
            self.residues[pdb_atm.parent].append((atm_obj, pdb_atm))

        for res in self.residues:
            if res.is_hetatm():
                continue
            self._check_mol(res)

        self._apply_modifications()
