from subprocess import Popen, PIPE

from luna.wrappers.chemaxon.util import check_license
from luna.wrappers.base import MolWrapper
from luna.util.exceptions import ProcessingFailed

import logging
logger = logging.getLogger()


class StereoisomerGenerator:

    def __init__(self,
                 max_stereoisomers=1000,
                 protect_double_bond_stereo=False,
                 protect_tetrahedral_stereo=False,
                 verify_3d=False,
                 mol_obj_type="rdkit",
                 cxcalc="cxcalc"):

        self.max_stereoisomers = max_stereoisomers
        self.protect_double_bond_stereo = protect_double_bond_stereo
        self.protect_tetrahedral_stereo = protect_tetrahedral_stereo
        self.verify_3d = verify_3d

        self.mol_obj_type = mol_obj_type
        self.cxcalc = cxcalc

    def _prep_opts(self):

        opt_list = []

        # Maximum number of stereoisomers to be generated
        opt_list += ["-m", str(self.max_stereoisomers)]

        # Protect aromaticity
        opt_list += ["-D", str(self.protect_double_bond_stereo)]

        # Protect charge
        opt_list += ["-T", str(self.protect_tetrahedral_stereo)]

        # Exclude antiaromatic compounds
        opt_list += ["-v", str(self.verify_3d)]

        # SMILES output.
        opt_list += ["-f", "smiles"]

        return opt_list

    def _run_cxcalc(self, mol_input):
        opt_list = self._prep_opts()

        self.mol_input = mol_input

        plugin = "stereoisomers"
        args = [self.cxcalc] + [mol_input] + [plugin] + opt_list

        p = Popen(args, stdout=PIPE, stderr=PIPE)
        try:
            stdout, stderr = p.communicate()
        except Exception:
            p.kill()
            raise

        error_lines = stderr.decode().strip().split("\n")
        if len(error_lines) == 1 and error_lines[0] == "":
            logger.debug("Stereoisomers for input '%s' created with success."
                         % mol_input)
        else:
            logger.error(stderr.decode())
            raise ProcessingFailed("cxcalc could not generate stereoisomers. "
                                   "Check the logs for more information.")

        self.stereoisomers = []
        for smi in stdout.decode().strip().split("\n"):
            smi = smi.strip()

            if smi == "":
                continue

            mol_obj = MolWrapper.from_smiles(smi)

            mol_id = "Stereoisomer %d" % len(self.stereoisomers)
            mol_obj.set_name(mol_id)

            self.stereoisomers.append(mol_obj)

    def generate(self, mol_input):

        self.stereoisomers = []

        check_license()

        self._run_cxcalc(mol_input)

        return self.stereoisomers
