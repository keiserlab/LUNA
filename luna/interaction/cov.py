from itertools import combinations
from openbabel import openbabel as ob

from luna.MyBio.PDB.PDBIO import Select

import logging

logger = logging.getLogger()


def is_covalently_bound(atm1, atm2):
    """Verifies if atoms ``atm1`` and ``atm2`` are covalently bound."""

    # Distance atom-atom
    dist = atm1 - atm2
    # Covalent radius.
    #   Note that we call title() to format atoms' symbol as in Open Babel.
    #   E.g.: ZN becomes Zn.
    cov1 = ob.GetCovalentRad(ob.GetAtomicNum(atm1.element.title()))
    cov2 = ob.GetCovalentRad(ob.GetAtomicNum(atm2.element.title()))

    # OpenBabel thresholds.
    if 0.4 <= dist <= cov1 + cov2 + 0.45:
        return True
    return False


def get_residue_cov_bonds(residue, select=Select()):
    """Get covalently bound atoms from residues or other molecules.

    Parameters
    ----------
    residue : :class:`~luna.MyBio.PDB.Residue.Residue`
        The residue or other molecule from which covalently
        bound atoms will be recovered.
    select : :class:`~luna.MyBio.PDB.PDBIO.Select`
        Decides which atoms will be consired.
        By default, all atoms are accepted.

    Returns
    -------
    list of tuple(:class:`~luna.MyBio.PDB.Atom.Atom`, \
            :class:`~luna.MyBio.PDB.Atom.Atom`)
        List of pairs of covalently bound atoms.
    """

    # Get valid atoms according to the provided selection function.
    trgt_res_atms = {atm.name: atm for atm in residue.get_atoms()
                     if select.accept_atom(atm)}

    cov_bonds = []
    pairs = combinations(trgt_res_atms.values(), 2)
    for atm1, atm2 in pairs:
        if is_covalently_bound(atm1, atm2):
            cov_bonds.append((atm1, atm2))

    return cov_bonds
