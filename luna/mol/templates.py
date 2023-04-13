import pandas as pd

from luna.util.default_values import LIGAND_EXPO_FILE
from luna.util.exceptions import MoleculeNotFoundError, MoleculeObjectTypeError
from luna.wrappers.base import MolWrapper

from rdkit.Chem.AllChem import AssignBondOrdersFromTemplate

import logging

logger = logging.getLogger()


class Template:
    """Standardize small molecules based on templates."""

    def assign_bond_order(self):
        """Assign bond order to a molecular object. However, this method
        is not implemented by default. Instead, you should use a class
        that inherits from `Template` and implements :meth:`assign_bond_order`.
        An example is the class `LigandExpoTemplate` that assigns bonds
        order based on ligands SMILES from
        `LigandExpo <http://ligand-expo.rcsb.org/>`_. Therefore, you should
        define your own logic beyond :meth:`assign_bond_order` that meets
        your goals."""
        raise NotImplementedError("Subclasses should implement this.")


class LigandExpoTemplate(Template):
    """Standardize small molecules based on templates (SMILES)
    from `LigandExpo <http://ligand-expo.rcsb.org/>`_.

    Parameters
    ----------
    lig_expo_file : str
        The `LigandExpo <http://ligand-expo.rcsb.org/>`_ file containing
        the SMILES and ligand ids.
    """

    def __init__(self, lig_expo_file=LIGAND_EXPO_FILE):

        self.lig_expo_file = lig_expo_file
        self._data = None

    @property
    def data(self):
        """:py:class:`pandas.DataFrame` : \
                The `LigandExpo <http://ligand-expo.rcsb.org/>`_ data."""
        if self._data is None:
            self._data = pd.read_csv(self.lig_expo_file,
                                     sep="\t+",
                                     na_filter=False,
                                     names=["smiles", "ligand_id"],
                                     usecols=[0, 1],
                                     engine='python').dropna()
        return self._data

    def get_ligand_smiles(self, lig_id):
        """Get SMILES for the ligand ``lig_id``.

        Parameters
        ----------
        lig_id : str
            The ligand identifier (PDB id) in
            `LigandExpo <http://ligand-expo.rcsb.org/>`_.

        Returns
        ----------
        smiles : str or None
            The ligand SMILES.
        """

        data = self.data[self.data["ligand_id"] == lig_id]["smiles"]
        if data.shape[0] == 0:
            return None
        return data.values[0]

    def assign_bond_order(self, mol_obj, lig_id):
        """Assign bond order to a molecular object based on its
        `LigandExpo <http://ligand-expo.rcsb.org/>`_ SMILES.

        Parameters
        ----------
        mol_obj : :class:`~luna.wrappers.base.MolWrapper`, \
                    :class:`rdkit.Chem.rdchem.Mol`, or \
                    :class:`openbabel.pybel.Molecule`
            A molecule to standardise.
        lig_id : str
            The ligand identifier (PDB id) in \
                `LigandExpo <http://ligand-expo.rcsb.org/>`_.

        Returns
        -------
        new_mol : :class:`~luna.wrappers.base.MolWrapper`, \
                    :class:`rdkit.Chem.rdchem.Mol`, or \
                    :class:`openbabel.pybel.Molecule`
            A standardized molecular object of the same type as ``mol_obj``.
        """
        tmp_mol_obj = MolWrapper(mol_obj)
        if tmp_mol_obj.is_rdkit_obj():
            smiles = self.get_ligand_smiles(lig_id)
            # Raise an exception when the template is not found.
            if smiles is None:
                error_msg = ("It is not possible to assign the bond orders to "
                             "the ligand %s because its corresponding "
                             "template was not found at Ligand Expo." % lig_id)
                raise MoleculeNotFoundError(error_msg)

            template = MolWrapper.from_smiles(smiles,
                                              mol_obj_type="rdkit").unwrap()

            # Note that the template molecule should have no explicit
            # hydrogens else the algorithm will fail.
            new_mol = AssignBondOrdersFromTemplate(template,
                                                   tmp_mol_obj.unwrap())

            if isinstance(mol_obj, MolWrapper):
                return MolWrapper(mol_obj)
            return new_mol
        else:
            logger.exception("Objects of type '%s' are not currently accepted."
                             % mol_obj.__class__)
            raise MoleculeObjectTypeError("Objects of type '%s' are not "
                                          "currently accepted."
                                          % mol_obj.__class__)
