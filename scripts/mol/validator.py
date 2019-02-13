from rdkit.Chem import (SanitizeFlags, SanitizeMol)

from mol.wrappers.base import (AtomWrapper, MolWrapper)
from mol.charge_model import OpenEyeModel

import logging

logger = logging.getLogger()


class MolValidator:

    def __init__(self, charge_model=OpenEyeModel(), references=None, fix_implicit_valence=False, fix_charges=False):
        self.charge_model = charge_model
        self.references = references
        self.fix_implicit_valence = fix_implicit_valence
        self.fix_charges = fix_charges

    def is_mol_valid(self, mol_obj):
        if not isinstance(mol_obj, MolWrapper):
            mol_obj = MolWrapper(mol_obj)

        # Check if the molecule have any errors...
        # Validations to be made...
        # TODO: check aromaticity.
        # TODO: Return list of errors
        is_mol_valid = True
        for atm in mol_obj.get_atoms():
            atm = AtomWrapper(atm)

            if not self.is_implicit_valence_valid(atm):
                is_mol_valid = False

            if not self.is_charge_valid(atm):
                is_mol_valid = False

        if not is_mol_valid:
            logger.warning("Invalid molecule: check the logs for more information.")

        return is_mol_valid

    def is_implicit_valence_valid(self, atm_obj):
        if not isinstance(atm_obj, AtomWrapper):
            atm_obj = AtomWrapper(atm_obj)

        # Atoms other than N are not currently evaluated because we did not find any similar errors with other atoms.
        if atm_obj.get_atomic_num() == 7:
            # It corrects quaternary ammonium Nitrogen errors.
            # While reading from PDBs with no pH correction, it happens that quaternary ammonium N is perceived as
            # having a valence equal to 5 (v5, hypervalent).
            # Thus, as a consequence, one additional (implicit) hydrogen is added to this N.
            if atm_obj.get_degree() == 4 and atm_obj.get_charge() == 0 and atm_obj.get_implicit_h_count() != 0:
                logger.warning("Atom # %d has incorrect number of implicit hydrogens." % atm_obj.get_id())

                if self.fix_implicit_valence:
                    logger.warning("'Fix implicit valence' option is set on. It will update the implicit "
                                   "valence of atom # %d from %d to 4." % (atm_obj.get_id(), atm_obj.get_valence()))
                    atm_obj.set_implicit_valence(4)
                    return True
                else:
                    return False
        return True

    def is_charge_valid(self, atm_obj):
        if not isinstance(atm_obj, AtomWrapper):
            atm_obj = AtomWrapper(atm_obj)

        expected_charge = self.get_expected_charge(atm_obj)

        if expected_charge is not None and expected_charge != atm_obj.get_charge():
            logger.warning("Atom # %d has incorrect charges defined." % atm_obj.get_id())

            if self.fix_charges:
                logger.warning("'Fix charges' option is set on. It will update the charge of atom # %d from %d to %d."
                               % (atm_obj.get_id(), atm_obj.get_charge(), expected_charge))
                atm_obj.set_charge(expected_charge)
                return True
            else:
                return False
        return True

    def get_expected_charge(self, atm_obj):
        if not isinstance(atm_obj, AtomWrapper):
            atm_obj = AtomWrapper(atm_obj)

        return self.charge_model.get_charge(atm_obj)

    def compare_to_ref(self, mol_obj, ref):
        # to be implemented
        pass


class RDKitValidator:

    def __init__(self, sanitize_opts=SanitizeFlags.SANITIZE_ALL):
        self.sanitize_opts = sanitize_opts

    def is_mol_valid(self, rdk_mol):
        try:
            SanitizeMol(rdk_mol, sanitizeOps=self.sanitize_opts)
            return True
        except Exception as e:
            logger.warning(e)
            return False
