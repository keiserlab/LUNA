import logging

from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain

from luna.pdb.io.helpers import save_to_file
from luna.pdb.io.selector import ChainSelector, ResidueSelector


logger = logging.getLogger()


class Extractor():
    """
    Extract chains or residues from a model or chain and save them to a PDB file.

    Parameters
    ----------
    entity : :class:`~luna.pdb.core.model.Model` or :class:`~luna.pdb.core.chain.Chain`
        A model or chain from which entities will be extracted.
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
        """
        Extract specific chains from a model and save to a file.

        Parameters
        ----------
        chains : iterable of str
            Chain IDs to extract.
        output_file : str
            Filepath to write the selected chains.
        """
        
        if self.entity.level != "M":
            logger.warning("Entity must be a Model object.")
            return
        
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
        """
        Extract specific residues from a chain and save to a file.

        Parameters
        ----------
        residues : iterable of tuple
            Residue IDs to extract. Each should be a (hetflag, resseq, icode) tuple.
        output_file : str
            Filepath to write the selected residues.
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
