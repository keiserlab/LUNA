from luna.MyBio.util import save_to_file
from luna.MyBio.selector import ChainSelector, ResidueSelector
from luna.MyBio.PDB.Model import Model
from luna.MyBio.PDB.Chain import Chain

import logging

logger = logging.getLogger()


class Extractor():
    """Extract chains or residues from ``entity``.

    Parameters
    ----------
    entity : :class:`~luna.MyBio.PDB.Entity.Model` or \
                :class:`~luna.MyBio.PDB.Entity.Chain`
        A model or chain object from where chains and residues
        will be extracted, respectively.
    """

    def __init__(self, entity):
        self.entity = entity

    @property
    def entity(self):
        return self._entity

    @entity.setter
    def entity(self, entity):
        if not isinstance(entity, (Model, Chain)):
            raise TypeError("Only objects of type %s or %s are accepted." %
                            (Model.__module__, Chain.__module__))
        self._entity = entity

    def extract_chains(self, chains, output_file):
        """Extract chains from ``entity`` and save it to ``output_file``.

        Parameters
        ----------
        chains : iterable of str
            A sequence of chains to extract from model ``entity``.
        output_file : str
            Save extracted chains to this file.
        """
        if self.entity.level != "M":
            logger.warning("Entity must be a Model object.")
        else:
            chain_sel = set()
            for chain in chains:
                if chain in self.entity.child_dict:
                    chain_sel.add(self.entity[chain])
                else:
                    logger.warning("Chain %s does not exist." % chain)

            if len(chain_sel) > 0:
                save_to_file(self.entity,
                             output_file,
                             ChainSelector(chain_sel))
            else:
                logger.warning("No valid chain to extract.")

    def extract_residues(self, residues, output_file):
        """Extract residues from ``entity`` and save it to ``output_file``.

        Parameters
        ----------
        residues : iterable
            A sequence of residues to extract from chain ``entity``.
        output_file : str
            Save extracted residues to to this file.
        """
        if self.entity.level != "C":
            logger.warning("Entity must be a Chain object.")
        else:
            res_sel = set()
            for res in residues:
                if res in self.entity.child_dict:
                    res_sel.add(self.entity[res])
                else:
                    logger.warning("Residue %s does not exist." % str(res))

            if len(res_sel) > 0:
                save_to_file(self.entity, 
                             output_file,
                             ResidueSelector(res_sel))
            else:
                logger.warning("No valid residue to extract.")
