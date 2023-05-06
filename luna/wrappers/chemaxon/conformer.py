from subprocess import Popen, PIPE
import tempfile
import os

from luna.wrappers.base import MolWrapper
from luna.wrappers.chemaxon.util import check_license
from luna.util.file import (new_unique_filename, is_file_valid,
                            get_file_format)
from luna.util.exceptions import (InvalidFileFormat, IllegalArgumentError,
                                  ProcessingFailed)

import logging
logger = logging.getLogger()


valid_formats = ["mol", "sdf", "sybyl", "mol2", "pdb", "xyz"]


class ConformerGenerator:

    def __init__(self, cleaner=3, skip_3d=False,
                 add_h=True, use_mmff94=True,
                 optimize=True, optim_criteria=3,
                 mol_obj_type="rdkit",
                 molconvert="molconvert"):

        self.cleaner = cleaner
        self.skip_3d = skip_3d
        self.add_h = add_h
        self.use_mmff94 = use_mmff94
        self.optimize = optimize
        self.optim_criteria = optim_criteria
        self.mol_obj_type = mol_obj_type
        self.molconvert = molconvert

    def _prep_opts(self):
        opt_list = []

        if self.cleaner in [0, 1, 2, 3]:
            if not self.skip_3d:
                opt_list += ["[c]{%d}{0}" % self.cleaner]
            else:
                opt_list += ["[c]{%d}{1}" % self.cleaner]
        else:
            raise IllegalArgumentError("The informed cleaner '%s' is invalid. "
                                       "Available options are 0, 1, 2, or 3."
                                       % str(self.cleaner))

        if self.add_h:
            opt_list += ["[prehydrogenize]"]

        if self.use_mmff94:
            opt_list += ["[mmff94]"]

        if self.optimize:
            if self.optim_criteria in [0, 1, 2, 3]:
                opt_list += ["[o]{1}{0}{%d}" % self.optim_criteria]
            else:
                raise IllegalArgumentError("The informed optimization "
                                           "criteria '%s' is invalid. "
                                           "Available options are "
                                           "0, 1, 2, or 3."
                                           % str(self.optim_criteria))
        else:
            opt_list += ["[o]{0}"]

        return opt_list

    def _run_molconvert(self, mol_input, is_smiles,
                        output_file, output_format):
        opt_list = self._prep_opts()

        input_list = [mol_input]
        if is_smiles:
            input_list = ["-s", mol_input]

        output_list = ["-o", output_file]
        args = ([self.molconvert, "-3:%s" % "".join(opt_list), output_format]
                + input_list + output_list)

        p = Popen(args, stdout=PIPE, stderr=PIPE)
        try:
            stdout, stderr = p.communicate()
        except Exception:
            p.kill()
            raise

        error_lines = stderr.decode().strip().split("\n")

        if len(error_lines) == 1 and error_lines[0] == "":
            error_lines = []

        # The list of errors are ok if they are one of the following.
        info_msgs = ["chemaxon.license.RemoteLicenseClient",
                     "INFO: Failed to retrieve available "
                     "licenses from the server:",
                     "INFO: Failed to retrieve license file from server:"]
        critical_errors_found = False
        for e in error_lines:
            if any([m in e for m in info_msgs]):
                continue
            critical_errors_found = True
            break

        # Raise ProcessingFailed error if any critical errors were found
        # or if the output file is empty.
        if (not critical_errors_found
                and os.path.getsize(self.output_file) > 0):

            logger.debug("Structure for input '%s' created with success."
                         % mol_input)

            # Try to parse the created molecular file.
            self.mol = MolWrapper.from_mol_file(self.output_file,
                                                output_format,
                                                self.mol_obj_type)
        else:
            logger.error(stderr.decode())
            raise ProcessingFailed("molconvert could not generate the 3D "
                                   "structure. Check the logs for more "
                                   "information.")

    def generate(self, mol_input, is_smiles=False,
                 output_file=None, output_format=None):

        check_license()

        if is_smiles:
            if not isinstance(mol_input, str):
                raise TypeError("The parameter 'is_smiles' was set on. "
                                "However, an invalid input type was provided. "
                                "A string was expected instead.")

        elif isinstance(mol_input, str) and not is_file_valid(mol_input):
            raise IllegalArgumentError("'%s' does not exist or "
                                       "is not a file." % mol_input)

        elif not isinstance(mol_input, str):
            raise IllegalArgumentError("'mol_input' must be a "
                                       "string representing a molecular "
                                       "file or SMILES.")

        if isinstance(output_file, str):
            if output_format is None:
                logger.debug("The output file format was not defined. "
                             "It will assume the format by the "
                             "file extension.")
                output_format = get_file_format(output_file, 1)

            if output_format not in valid_formats:
                raise InvalidFileFormat("Output format '%s' is not valid."
                                        % output_format)

            logger.debug("Output format: %s" % output_format)

        elif output_file is not None:
            raise IllegalArgumentError("The 'output_file' should be "
                                       "a string or None.")

        if output_file is None:
            output_file = new_unique_filename(tempfile.gettempdir()) + ".mol"
            output_format = "mol"

        self.mol = None
        self.mol_input = mol_input
        self.output_file = output_file

        self._run_molconvert(mol_input, is_smiles,
                             output_file, output_format)
