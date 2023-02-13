from rdkit.Chem import SanitizeFlags, SanitizeMol
from openbabel.openbabel import OBSmartsPattern

from luna.wrappers.base import BondType, AtomWrapper, MolWrapper
from luna.mol.charge_model import OpenEyeModel

import logging

logger = logging.getLogger()


class MolValidator:
    """Validate and fix molecules with the errors most commonly
    found when parsing PDB files with Open Babel.

    Parameters
    ----------
    charge_model : :class:`~luna.mol.charge_model.ChargeModel`
        A charge model object. By default, the implementation of OpenEye charge
        model is used.
    fix_nitro : bool
        If True, fix nitro groups whose nitrogens are perceived as having a
        valence of 5 and that are double-bonded to the oxygens.
        SMILES representation of the error: '[$([NX3v5]([!#8])(=O)=O)]'.
    fix_amidine_and_guanidine : bool
        If True, fix amidine and guanidine-like groups.
        When the molecule is ionized, it may happen that the charge is
        incorrectly assigned to the central carbon, which ends up with a +1
        charge and the nitrogen double-bonded to it ends up with a +0 charge.
    fix_valence : bool
        If True, fix the valence of atoms.
        Currently, we only detect and fix errors of quaternary ammonium
        nitrogens. These charged atoms are sometimes perceived as having a
        valence equal to 5 and no charge, i.e., Open Babel considers those
        nitrogens as hypervalent.
    fix_charges : bool
        If True, fix charges on basis of the charge model ``charge_model``.
    """

    def __init__(self, charge_model=OpenEyeModel(), fix_nitro=True,
                 fix_amidine_and_guanidine=True, fix_valence=True,
                 fix_charges=True, bound_to_metals=None):
        self.charge_model = charge_model
        self.fix_nitro = fix_nitro
        self.fix_amidine_and_guanidine = fix_amidine_and_guanidine
        self.fix_valence = fix_valence
        self.fix_charges = fix_charges
        self.bound_to_metals = bound_to_metals or []

    def validate_mol(self, mol_obj):
        """Validate molecule ``mol_obj``.

        It will first try to fix the provided molecule according to the
        initialization flags. If any errors are detected and if it fails
        to fix them, then this function may return False, i.e., the molecule
        is invalid.

        Parameters
        ----------
        mol_obj : :class:`~luna.wrappers.base.MolWrapper`, \
                :class:`rdkit.Chem.rdchem.Mol`, \
                or :class:`openbabel.pybel.Molecule`
            The molecule.

        Returns
        -------
        is_valid : bool
            If the molecule is valid or not.
        """

        if not isinstance(mol_obj, MolWrapper):
            mol_obj = MolWrapper(mol_obj)

        # Validations to be made...
        #       TODO: check aromaticity.
        #       TODO: return list of errors

        if self.bound_to_metals:
            self._fix_atoms_bound_to_metals(mol_obj)

        if self.fix_nitro:
            self._fix_nitro_substructure_and_charge(mol_obj)

        if self.fix_amidine_and_guanidine:
            self._fix_amidine_and_guanidine_charges(mol_obj)

        # Check if the molecule has errors...
        is_valid = True
        for atm_obj in mol_obj.get_atoms():
            if not self._is_valence_valid(atm_obj):
                is_valid = False

            if not self._is_charge_valid(atm_obj):
                is_valid = False

        if not is_valid:
            logger.debug("Invalid molecule: check the logs for more "
                         "information.")

        return is_valid

    def _fix_atoms_bound_to_metals(self, mol_obj):

        # Iterate over each atom that has a dative bond with a metal
        # to evaluate if its charge is correct.
        for idx, pdb_atm in self.bound_to_metals:

            atm_obj = mol_obj.get_atom_by_idx(idx)
            metals = self.bound_to_metals[(idx, pdb_atm)]

            if not atm_obj:
                raise KeyError("Atom #%d does not exist." % idx)

            # Fix HIS:ND1 and HIS:NE2 that are bound to metals.
            if (pdb_atm.name in ["ND1", "NE2"]
                    and pdb_atm.parent.resname == "HIS"):

                if atm_obj.get_charge() == 1:
                    logger.debug("Atom #%d has incorrect charge. It will "
                                 "update its charge from 1 to 0." % idx)

                    atm_obj.set_charge(0)
                    self._remove_explicit_hydrogens(atm_obj)

            # Fix main chain N or LYS:NZ that are bound to metals and
            # that has a charge = +1.
            elif (pdb_atm.name == "N"
                    or (pdb_atm.name == "NZ"
                        and pdb_atm.parent.resname == "LYS")):

                if atm_obj.get_charge() == 1:
                    logger.debug("Atom #%d has incorrect charge. It will "
                                 "update its charge from 1 to 0 and fix "
                                 "its number of bound hydrogens." % idx)

                    # Set residue charge to 0.
                    atm_obj.set_charge(0)

                    # Now, loop over this atom's bonds to identify
                    # all bonds with hydrogens. Then, it will identify
                    # the closest H to the metal and remove it.
                    #
                    # That's because residue can only coordinate
                    # one metal at once.
                    h_to_remove = None
                    min_dist = float('inf')
                    for b in atm_obj.get_bonds():
                        if b.get_partner_atom(atm_obj).get_symbol() == "H":
                            h_obj = b.get_partner_atom(atm_obj)

                            for metal_obj in metals:
                                # Distance between the H and the metal.
                                d = h_obj - metal_obj
                                if d < min_dist:
                                    h_to_remove = h_obj
                                    min_dist = d

                    # Remove the closest H to the metal.
                    if h_to_remove:
                        ob_mol = atm_obj.parent.unwrap()
                        ob_mol.DeleteAtom(h_to_remove.unwrap())

    def _fix_nitro_substructure_and_charge(self, mol_obj):
        ob_smart = OBSmartsPattern()
        # Invalid nitro pattern.
        ob_smart.Init("[$([NX3v5]([!#8])(=O)=O)]")
        if ob_smart.Match(mol_obj.unwrap()):
            logger.debug("One or more invalid nitro substructures "
                         "('*-N(=O)=O') were found. It will try to "
                         "substitute them to '*-[N+]([O-])=O'.")

            # Iterate over each Nitro group in the molecule.
            for ids in ob_smart.GetUMapList():
                # Get the N atom.
                atm_obj = [AtomWrapper(mol_obj.GetAtom(i)) for i in ids][0]
                for bond in atm_obj.get_bonds():
                    partner = bond.get_partner_atom(atm_obj)

                    if (partner.get_symbol() == "O"
                            and bond.get_bond_type() == BondType.DOUBLE):
                        # Change double bond to single bond.
                        bond.set_bond_type(BondType.SINGLE)
                        # Attributes a +1 charge to the N.
                        atm_obj.set_charge(1)
                        # Attributes a -1 charge to the O.
                        partner.set_charge(-1)

                        # It needs to update only one of the oxygen bonds.
                        break

            ob_smart = OBSmartsPattern()
            # Valid nitro pattern.
            ob_smart.Init("[$([NX3v4+](=O)[O-])][!#8]")
            if ob_smart.Match(mol_obj.unwrap()):
                logger.debug("Invalid nitro substructures ('*-N(=O)=O') "
                             "successfully substituted to '*-[N+]([O-])=O'.")

    def _fix_amidine_and_guanidine_charges(self, mol_obj):
        # These errors occur with guanidine-like substructures when the
        # molecule is ionized. It happens that the charge is incorrectly
        # assigned to the central carbon, so the guanidine-like C ends up with
        # a +1 charge and the N with a double bond to the central C ends up
        # with a +0 charge. To fix it, we assign the correct charges to
        # the N (+1) and C (0).

        ob_smart = OBSmartsPattern()
        # Invalid amidine and guanidine pattern.
        ob_smart.Init("[$([NH1X2v3+0](=[CH0X3+1](N)))]")
        if ob_smart.Match(mol_obj.unwrap()):
            logger.debug("One or more amidine/guanidine substructures with no "
                         "charge were found. It will try to attribute a +1 "
                         "charge to the N bound to the central carbon with a "
                         "double bond.")

            # Iterate over each Amidine/Guanidine group in the molecule.
            for ids in ob_smart.GetUMapList():
                # Get the N atom.
                atm_obj = [AtomWrapper(mol_obj.GetAtom(i)) for i in ids][0]

                for bond in atm_obj.get_bonds():
                    partner = bond.get_partner_atom(atm_obj)

                    if partner.get_symbol() == "C":
                        # Attributes a +1 charge to the N and corrects its
                        # degree to three, i.e., the N will become '=NH2+'.
                        atm_obj.set_charge(1)

                        # Set the number of implicit Hydrogens to 1 (=NH2+).
                        #
                        # Explanation: the above initialized Smarts pattern
                        #      checks for nitrogens containing 1 explicit
                        #      hydrogen ([NH1X2v3+0]). Thus, to correctly
                        #      update its valence to 4 (=NH2+), it needs to
                        #      add a new implicit hydrogen.
                        atm_obj.unwrap().SetImplicitHCount(1)

                        # Remove any charges in the C.
                        partner.set_charge(0)

            ob_smart = OBSmartsPattern()
            # Valid amidine and guanidine pattern.
            ob_smart.Init("[$([NH2X3v4+1](=[CH0X3+0](N)))]")
            if ob_smart.Match(mol_obj.unwrap()):
                logger.debug("Invalid amidine/guanidine substructures were "
                             "correctly charged.")

    def _is_valence_valid(self, atm_obj):
        if not isinstance(atm_obj, AtomWrapper):
            atm_obj = AtomWrapper(atm_obj)

        # Currently, atoms other than N are not evaluated as no valence error
        # has been identified for them.
        if atm_obj.get_atomic_num() == 7:

            # Molecules containing quaternary ammonium N.
            #
            #    - Open Babel 2.3.2 incorrectly sets N with:
            #         . Valence = 5
            #         . Charge  = 0
            #         . # Implicit Hs = 0
            #
            #    - Solution: set the number of implicit Hs to 0 and set the
            #                charge to +1.
            #
            if atm_obj.get_valence() == 5 and atm_obj.get_charge() == 0:
                logger.debug("Atom #%d has incorrect valence and charge."
                             % atm_obj.get_idx())

                if self.fix_valence:
                    logger.debug("'Fix valence' option is set on. It will "
                                 "update the valence of atom #%d from %d "
                                 "to 4 and correct its charge."
                                 % (atm_obj.get_idx(), atm_obj.get_valence()))

                    # Fix the number of implicit Hs and charge.
                    atm_obj.unwrap().SetImplicitHCount(0)
                    atm_obj.set_charge(1)
                    return True
                return False

            # Molecules containing quaternary ammonium N.
            #
            #    - Open Babel 3.x incorrectly sets N with:
            #         . Valence = 5
            #         . Charge  = 1
            #         . Num. Implicit Hs = 0
            #         . Num. Explicit Hs = 1
            #
            #    - Solution: remove the H and its bond with the N.
            #
            elif (atm_obj.get_valence() == 5 and atm_obj.get_charge() == 1
                    and atm_obj.get_h_count() == 1):

                logger.debug("Atom #%d has incorrect valence and number of "
                             "hydrogens."
                             % atm_obj.get_idx())

                if self.fix_valence:
                    logger.debug("'Fix valence' option is set on. It will "
                                 "update the valence of atom #%d from %d "
                                 "to 4 and correct the number of hydrogens "
                                 "bound to it."
                                 % (atm_obj.get_idx(), atm_obj.get_valence()))

                    self._remove_explicit_hydrogens(atm_obj)

                    return True
                return False
        return True

    def _is_charge_valid(self, atm_obj):
        if not isinstance(atm_obj, AtomWrapper):
            atm_obj = AtomWrapper(atm_obj)

        expected_charge = self._get_expected_charge(atm_obj)

        if (expected_charge is not None
                and expected_charge != atm_obj.get_charge()):
            logger.debug("Atom #%d has incorrect charges defined."
                         % atm_obj.get_idx())

            if self.fix_charges:
                logger.debug("'Fix charges' option is set on. It will update "
                             "the charge of atom #%d from %d to %d."
                             % (atm_obj.get_idx(), atm_obj.get_charge(),
                                expected_charge))
                atm_obj.set_charge(expected_charge)
                return True
            return False
        return True

    def _get_expected_charge(self, atm_obj):
        if not isinstance(atm_obj, AtomWrapper):
            atm_obj = AtomWrapper(atm_obj)

        return self.charge_model.get_charge(atm_obj)

    def _remove_explicit_hydrogens(self, atm_obj):
        if not isinstance(atm_obj, AtomWrapper):
            atm_obj = AtomWrapper(atm_obj)

        delete_hs = []
        for b in atm_obj.get_bonds():
            if b.get_partner_atom(atm_obj).get_symbol() == "H":
                delete_hs.append(b.get_partner_atom(atm_obj))

        for hs_obj in delete_hs:
            atm_obj.parent.unwrap().DeleteAtom(hs_obj.unwrap())


class RDKitValidator:
    """Check if RDKit molecular objects are valid or not.

    Parameters
    ----------
    sanitize_opts : :class:`rdkit.Chem.rdchem.SanitizeFlags`
        Sanitization operations to be carried out.
    """

    def __init__(self, sanitize_opts=SanitizeFlags.SANITIZE_ALL):
        self.sanitize_opts = sanitize_opts

    def is_valid(self, rdk_mol):
        """Try to sanitize the molecule ``rdk_mol``.
        If it succeeds, returns True. Otherwise, returns False.

        Parameters
        ----------
        rdk_mol : :class:`rdkit.Chem.rdchem.Mol`

        Returns
        -------
        is_valid : bool
            If the molecule is valid or not.
        """
        try:
            SanitizeMol(rdk_mol, sanitizeOps=self.sanitize_opts)
            return True
        except Exception as e:
            logger.exception(e)
            return False
