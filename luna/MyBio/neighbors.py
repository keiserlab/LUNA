from itertools import product

from luna.MyBio.PDB.PDBIO import Select
from luna.interaction.cov import is_covalently_bound

import logging

logger = logging.getLogger()


class Neighbors:

    def __init__(self, prev_res=None, next_res=None):

        self.predecessor = prev_res
        self.successor = next_res

    def has_predecessor(self):
        return self.predecessor is not None

    def has_successor(self):
        return self.successor is not None


def get_residue_neighbors(residue, select=Select(), verbose=True):
    """Get all neighbors from a residue.

    In the case of an amino acid that is part of a peptide bond, its neighbors
    are any predecessor or successor residues.
    The same idea applies to nucleic acids.

    Parameters
    ----------
    residue : :class:`~luna.MyBio.PDB.Residue.Residue`
        The residue or other molecule from which covalently bound molecules
        will be recovered.
    select : :class:`~luna.MyBio.PDB.PDBIO.Select`
        Decides which atoms will be consired. By default,
        all atoms are accepted.

    Returns
    -------
    neighbors : dict of {str : :class:`~luna.MyBio.PDB.Residue.Residue`}
        A dictionary containing the predecessor (``previous``) or
        successor (``next``) molecules.
    """

    # Get valid atoms according to the provided selection function.
    trgt_res_atms = {atm.name: atm for atm in residue.get_atoms()
                     if select.accept_atom(atm)}

    if residue.is_residue():
        if "N" not in trgt_res_atms and verbose:
            logger.debug("There is a missing N in the residue %s. It may have "
                         "been filtered out by the provided selection "
                         "function. So, the predecessor residue cannot be "
                         "identified." % residue)
        if "C" not in trgt_res_atms and verbose:
            logger.debug("There is a missing C in the residue %s. It may have "
                         "been filtered out by the provided selection "
                         "function. So, the successor residue cannot be "
                         "identified." % residue)

        # If neither N and C are available in the valid atom list.
        if "N" not in trgt_res_atms and "C" not in trgt_res_atms:
            return {}

        neighbors = Neighbors()
        # If the chain has a residue coming before the target residue.
        if residue.idx - 1 >= 0:
            # First residue before the target in the chain list.
            prev_res = residue.parent.child_list[residue.idx - 1]
            # Get valid atoms according to the provided selection function.
            prev_res_atms = {atm.name: atm for atm in prev_res.get_atoms()
                             if select.accept_atom(atm)}

            # A peptide bond exists between the C of one amino acid
            # and the N of another.
            if "C" in prev_res_atms and "N" in trgt_res_atms:
                if is_covalently_bound(trgt_res_atms["N"], prev_res_atms["C"]):
                    neighbors.predecessor = prev_res
                elif verbose:
                    logger.debug("The first residue before %s is too distant "
                                 "to fulfill the covalent bond threshold. It "
                                 "may be an indication of bad atom "
                                 "positioning or that there are missing "
                                 "residues." % residue)
            elif "C" not in prev_res_atms and verbose:
                logger.debug("There is a missing C in %s, the first residue "
                             "before %s in the chain list. It may have been "
                             "filtered out by the provided selection "
                             "function. So, the predecessor residue cannot "
                             "be identified." % (prev_res, residue))
        # Otherwise, it could mean that the residue is the first one in
        # the sequence or there are missing residues.
        elif verbose:
            logger.debug("The residue %s seems not to have any predecessor "
                         "residue. It may be the first in the chain sequence "
                         "or there are missing residues." % residue)

        # If the chain has a residue coming after the target residue.
        if residue.idx + 1 < len(residue.parent.child_list):
            # First residue after the target in the chain list.
            next_res = residue.parent.child_list[residue.idx + 1]
            # Get valid atoms according to the provided selection function.
            next_res_atms = {atm.name: atm for atm in next_res.get_atoms()
                             if select.accept_atom(atm)}

            # A peptide bond exists between the C of one amino acid and
            # the N of another.
            if "C" in trgt_res_atms and "N" in next_res_atms:
                if is_covalently_bound(trgt_res_atms["C"], next_res_atms["N"]):
                    neighbors.successor = next_res
                elif verbose:
                    logger.debug("The first residue after %s is too distant "
                                 "to fulfill the covalent thresholds. It may "
                                 "be an indication of bad atom positioning or "
                                 "that there are missing residues." % residue)
            elif "N" in next_res_atms and verbose:
                logger.debug("There is a missing N in %s, the first residue "
                             "after %s in the chain list. It may have been "
                             "filtered out by the provided selection "
                             "function. So, the predecessor residue cannot be "
                             "identified." % (next_res, residue))
        # Otherwise, it could mean that the residue is the last one in the
        # sequence or there are missing residues.
        elif verbose:
            logger.debug("The residue %s seems not to have any successor "
                         "residue. It may be the last in the chain sequence "
                         "or there are missing residues." % residue)
        return neighbors
    else:
        # First residue before the target in the chain list.
        prev_res = residue.parent.child_list[residue.idx - 1]
        # Get valid atoms according to the provided selection function.
        prev_res_atms = {atm.name: atm for atm in prev_res.get_atoms()
                         if select.accept_atom(atm)}

        neighbors = Neighbors()
        # If the chain has a residue coming before the target residue.
        if residue.idx - 1 >= 0:
            # First residue before the target in the chain list.
            prev_res = residue.parent.child_list[residue.idx - 1]
            # Get valid atoms according to the provided selection function.
            prev_res_atms = {atm.name: atm for atm in prev_res.get_atoms()
                             if select.accept_atom(atm)}

            for trgt_atm, prev_atm in product(trgt_res_atms.values(),
                                              prev_res_atms.values()):
                if is_covalently_bound(trgt_atm, prev_atm):
                    neighbors.predecessor = prev_res
                    break

            if not neighbors.has_predecessor() and verbose:
                logger.debug("The first residue after %s is too distant to "
                             "fulfill the covalent thresholds. It may be an "
                             "indication of bad atom positioning or that "
                             "there are missing residues." % residue)
        # Otherwise, it could mean that the residue is the first one in the
        # sequence or there are missing residues.
        elif verbose:
            logger.debug("The residue %s seems not to have any predecessor "
                         "residue. It may be the first in the chain sequence "
                         "or there are missing residues." % residue)

        # If the chain has a residue coming after the target residue.
        if residue.idx + 1 < len(residue.parent.child_list):
            # First residue after the target in the chain list.
            next_res = residue.parent.child_list[residue.idx + 1]
            # Get valid atoms according to the provided selection function.
            next_res_atms = {atm.name: atm for atm in next_res.get_atoms()
                             if select.accept_atom(atm)}

            # Check each pair of atoms for covalently bonded atoms.
            for trgt_atm, next_atm in product(trgt_res_atms.values(),
                                              next_res_atms.values()):
                if is_covalently_bound(trgt_atm, next_atm):
                    neighbors.successor = next_res
                    break

            if not neighbors.has_successor() and verbose:
                logger.debug("The first residue after %s is too distant to "
                             "fulfill the covalent thresholds. It may be an "
                             "indication of bad atom positioning or that "
                             "there are missing residues." % residue)
        # Otherwise, it could mean that the residue is the last one in the
        # sequence or there are missing residues.
        elif verbose:
            logger.debug("The residue %s seems not to have any successor "
                         "residue. It may be the last in the chain sequence "
                         "or there are missing residues." % residue)
        return neighbors
    return {}
