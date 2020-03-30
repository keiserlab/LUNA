from luna.MyBio.util import save_to_file
from luna.MyBio.selector import ChainSelector, ResidueSelector

import logging

logger = logging.getLogger()


class Extractor():
    """Selects everything and creates a new PDB output.
    """

    def __init__(self, entity):
        self.entity = entity

    def extract_chains(self, chains, output_file):
        if self.entity.level != "M":
            logger.warning("Target must be a Model object.")
        else:
            chain_sel = set()
            for chain in chains:
                if chain in self.target.child_dict:
                    chain_sel.add(self.entity[chain])
                else:
                    logger.warning("Chain %s does not exist." % chain)

            if len(chain_sel) > 0:
                save_to_file(self.entity, output_file, ChainSelector(chain_sel))
            else:
                logger.warning("No valid chain to extract.")

    def extract_residues(self, residues, output_file):
        if (self.entity.level != "C"):
            logger.warning("Target must be a Chain object.")
        else:
            res_sel = set()
            for res in residues:
                if (res in self.entity.child_dict):
                    res_sel.add(self.entity[res])
                else:
                    logger.warning("Residue %s does not exist." % str(res))

            if len(res_sel) > 0:
                save_to_file(self.entity, output_file, ResidueSelector(res_sel))
            else:
                logger.warning("No valid residue to extract.")
