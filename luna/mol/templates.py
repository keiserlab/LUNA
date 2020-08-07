import pandas as pd

from luna.util.default_values import LIGAND_EXPO_FILE
from luna.util.exceptions import MoleculeObjectTypeError
from luna.mol.wrappers.base import MolWrapper

from rdkit.Chem.AllChem import AssignBondOrdersFromTemplate

import logging

logger = logging.getLogger()


class LigandExpoTemplate:

    def __init__(self, lig_expo_file=LIGAND_EXPO_FILE):

        self.lig_expo_file = lig_expo_file
        self._data = None

    @property
    def data(self):
        if self._data is None:
            self._data = pd.read_csv(self.lig_expo_file, sep="\t+", na_filter=False,
                                     names=["smiles", "ligand_id"], usecols=[0, 1], engine='python').dropna()
        return self._data

    def as_json(lig_expo_file, json_file):
        ligands = {}
        return ligands

    def get_ligand_smiles(self, lig_id):
        return self.data[self.data["ligand_id"] == lig_id]["smiles"].values[0]

    def assign_bond_order(self, mol_obj, lig_id):
        tmp_mol_obj = MolWrapper(mol_obj)
        if tmp_mol_obj.is_rdkit_obj():
            smiles = self.get_ligand_smiles(lig_id)
            template = MolWrapper.from_smiles(smiles, mol_obj_type="rdkit").unwrap()

            # Note that the template molecule should have no explicit hydrogens else the algorithm will fail.
            new_mol = AssignBondOrdersFromTemplate(template, tmp_mol_obj.unwrap())

            if isinstance(mol_obj, MolWrapper):
                return MolWrapper(mol_obj)
            return new_mol
        else:
            logger.exception("Objects of type '%s' are not currently accepted." % mol_obj.__class__)
            raise MoleculeObjectTypeError("Objects of type '%s' are not currently accepted." % mol_obj.__class__)
