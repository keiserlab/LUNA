from subprocess import Popen, PIPE

from luna.util.exceptions import (InvalidFileFormat, IllegalArgumentError,
                                  ProcessingFailed)
from luna.util.default_values import ANTECHAMBER
from luna.util.file import is_file_valid, get_file_format

import logging
logger = logging.getLogger()


INFORMATS = ["ac", "mol2", "pdb", "mpdb", "prepi", "prepc", "gzmat",
             "gcrt", "mopint", "mopcrt", "gout", "mopout", "alc",
             "csd", "mdl", "hin", "rst"]
OUTFORMATS = INFORMATS

CHARGE_METHODS = ["resp", "cm2", "mul", "rc", "bcc", "esp", "gas", "wc"]

ATOM_TYPES = ["gaff", "amber", "bcc", "sybyl"]


def _prep_opts(opts, prefix=""):
    opt_list = []
    if opts is not None:
        for key in opts:
            opt_list.append("-%s%s" % (prefix, key))

            if opts[key] is not None:
                opt_list.append(str(opts[key]))
    return opt_list


def add_charges(input_file, output_file, total_charge, input_format=None,
                output_format=None, charge_method="bcc", atom_type="sybyl",
                opts=None, antechamber=ANTECHAMBER):

    logger.debug("Charges will be added to the molecule at '%s'." % input_file)

    if is_file_valid(input_file):
        if input_format is None:
            logger.debug("Input file format not defined. "
                         "It will assume the format from the file extension.")
            input_format = get_file_format(input_file)

        if input_format not in INFORMATS:
            raise InvalidFileFormat("Input format '%s' is not accepted "
                                    "by Antechamber." % input_format)

        logger.debug("Input format: %s" % input_format)
    else:
        raise OSError("File '%s' does not exist or is not a valid file."
                      % input_format)

    if isinstance(output_file, str):
        if output_format is None:
            logger.debug("Output file format not defined. "
                         "It will assume the format by the file extension.")
            output_format = get_file_format(output_file, 1)

        if output_format not in OUTFORMATS:
            raise InvalidFileFormat("Input format '%s' does not exist."
                                    % output_format)

        logger.debug("Output format: %s" % output_format)
    else:
        raise IllegalArgumentError("'output_file' must be a string.")

    if charge_method not in CHARGE_METHODS:
        raise IllegalArgumentError("The charge method '%s' does not exist "
                                   "in Antechamber." % charge_method)

    if atom_type not in ATOM_TYPES:
        raise IllegalArgumentError("The atom type '%s' does not exist "
                                   "in Antechamber." % atom_type)

    opts = opts or {}
    if "nc" in opts:
        logger.warning("The parameter for the total charge found "
                       "in 'opts' will be overwritten.")
    if "c" in opts:
        logger.warning("The parameter for the charge method found "
                       "in 'opts' will be overwritten.")
    if "at" in opts:
        logger.warning("The parameter for the atom type found "
                       "in 'opts' will be overwritten.")

    opts["nc"] = total_charge
    opts["c"] = charge_method
    opts["at"] = atom_type
    opts["pf"] = "y"
    opt_list = _prep_opts(opts)

    input_list = ['-i', input_file, "-fi", input_format]
    output_list = ['-o', output_file, "-fo", output_format]

    args = [antechamber] + input_list + output_list + opt_list

    p = Popen(args, stdout=PIPE, stderr=PIPE)
    try:
        stdout, stderr = p.communicate()
    except Exception:
        p.kill()
        raise

    error_lines = stderr.decode().strip().split("\n")

    if len(error_lines) == 1 and error_lines[0] == "":
        logger.debug("File '%s' created with success." % output_file)
    else:
        raise ProcessingFailed("Antechamber could not add charges "
                               "to the provided molecule.")
