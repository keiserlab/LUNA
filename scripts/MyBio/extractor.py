from MyBio.util import try_save_2pdb
from MyBio.selector import (ChainSelector, ResidueSelector)

import logging
logger = logging.getLogger(__name__)


class Extractor():
    """Selects everything and creates a new PDB output.
    """

    def __init__(self, entity):
        self.entity = entity

    def extract_chains(self, chains, outputFile):
        if (self.entity.level != "M"):
            logger.warning("Target must be a Model object.")
        else:
            selChains = set()
            for chain in chains:
                if (chain in self.target.child_dict):
                    selChains.add(self.entity[chain])
                else:
                    logger.warning("Chain %s does not exist." % chain)

            if (len(selChains) > 0):
                try_save_2pdb(self.entity, outputFile,
                              ChainSelector(selChains))
            else:
                logger.warning("No valid chain to extract.")

    def extract_residues(self, residues, outputFile):
        if (self.entity.level != "C"):
            logger.warning("Target must be a Chain object.")
        else:
            selResidues = set()

            for res in residues:
                if (res in self.entity.child_dict):
                    selResidues.add(self.entity[res])
                else:
                    logger.warning("Residue %s does not exist." % str(res))

            if (len(selResidues) > 0):
                try_save_2pdb(self.entity, outputFile,
                              ResidueSelector(selResidues))
            else:
                logger.warning("No valid residue to extract.")
