from luna.wrappers.base import AtomWrapper

import logging
logger = logging.getLogger()


# TODO: implement this other model: https://github.com/openbabel/openbabel/blob/master/src/formats/mdlvalence.h


class ChargeModel:
    """Implementation of a charge model."""

    def get_charge(self):
        raise NotImplementedError("Subclasses should implement this.")


class OpenEyeModel(ChargeModel):
    """Implementation of OpenEye charge model."""

    def get_charge(self, atm_obj):
        """Get the formal charge for atom ``atom_obj``.

        Currently, only formal charges for the elements Hydrogen, Carbon,
        Nitrogen, Oxygen, Phosphorus, Sulfur, Chlorine, Fluorine, Bromine,
        Iodine, Magnesium, Calcium, Zinc, Lithium, Sodium, Potassium,
        and Boron can be recovered.

        Parameters
        ----------
        atm_obj : :class:`~luna.wrappers.base.AtomWrapper`
            The target atom.

        Examples
        --------
        >>> from luna.mol.charge_model import OpenEyeModel
        >>> from luna.wrappers.base import MolWrapper
        >>> cm = OpenEyeModel()
        >>> mol_obj = MolWrapper.from_smiles("[NH3+]CC([O-])=O")
        >>> for atm_obj in mol_obj.get_atoms():
        >>>     formal_charge = cm.get_charge(atm_obj)
        >>>     print("Charge for atom #%d (%s): %d" % (atm_obj.get_idx(),
        ...                                             atm_obj.get_symbol(),
        ...                                             formal_charge))
        Charge for atom #0 (N): 1.
        Charge for atom #1 (C): 0.
        Charge for atom #2 (C): 0.
        Charge for atom #3 (O): -1.
        Charge for atom #4 (O): 0.

        """
        formal_charge = None

        atm_num = atm_obj.get_atomic_num()
        valence = atm_obj.get_valence()

        # Hydrogen
        if atm_num == 1:
            if valence != 1:
                formal_charge = 1
            elif valence == 1:
                formal_charge = 0

        # Carbon
        elif atm_num == 6:
            if valence == 3:
                has_polar_nb = False
                nbs = atm_obj.get_neighbors()
                for nb_atm_obj in nbs:
                    nb_atm_obj = AtomWrapper(nb_atm_obj)
                    # Polar neighbor: N, O or S
                    if nb_atm_obj.get_atomic_num() in (7, 8, 16):
                        has_polar_nb = True
                        break
                if has_polar_nb is True:
                    formal_charge = 1
                else:
                    formal_charge = -1
            elif valence == 4:
                formal_charge = 0

        # Nitrogen
        elif atm_num == 7:
            if valence == 2:
                formal_charge = -1
            elif valence == 4:
                formal_charge = 1
            elif valence == 3:
                formal_charge = 0

        # Oxygen
        elif atm_num == 8:
            if valence == 1:
                formal_charge = -1
            elif valence == 3:
                formal_charge = 1
            elif valence == 2:
                formal_charge = 0

        # Phosphorus
        elif atm_num == 15:
            if valence == 4:
                formal_charge = 1

        # Sulfur
        elif atm_num == 16:
            if valence == 1:
                formal_charge = -1
            elif valence == 3:
                formal_charge = 1
            elif valence == 5:
                formal_charge = -1
            elif valence == 4:
                if atm_obj.get_degree() == 4:
                    formal_charge = 2
            elif valence == 2 or valence == 6:
                formal_charge == 0

        # Chlorine
        elif atm_num == 17:
            if valence == 0:
                formal_charge = -1
            elif valence == 4:
                formal_charge = 3
            elif valence == 1:
                formal_charge = 0

        # Fluorine, Bromine, Iodine
        elif atm_num == 9 or atm_num == 35 or atm_num == 53:
            if valence == 0:
                formal_charge = -1
            elif valence == 1:
                formal_charge = 0

        # Magnesium, Calcium, Zinc
        elif atm_num == 12 or atm_num == 20 or atm_num == 30:
            if valence == 0:
                formal_charge = 2
            elif valence == 2:
                formal_charge = 0

        # Lithium, Sodium, Potassium
        elif atm_num == 3 or atm_num == 11 or atm_num == 19:
            if valence == 0:
                formal_charge = 1
            elif valence == 1:
                formal_charge = 0

        # Boron
        elif atm_num == 5:
            # If the valence is four, the formal charge is -1.
            # OBS: there is an error in OpenEye chargel model text,
            # they said that Boron should have a charge of +1 when
            # the valence is 4.
            # However, the MDL Valence Model says it should be -1.
            if valence == 4:
                formal_charge = -1
            elif valence == 3:
                formal_charge = 0

        # Iron: 26
        # If the valence is zero, the formal charge is +3 if the partial charge is 3.0, and +2 otherwise.
        # Copper: 29
        # If the valence is zero, the formal charge is +2 if the partial charge is 2.0, and +1 otherwise.
        # else:
        # For the remaining elements, if the valence of an atom is zero, its formal charge is set from its partial charge.

        return formal_charge
