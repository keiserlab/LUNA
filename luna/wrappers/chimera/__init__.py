from pathlib import Path
from subprocess import Popen, PIPE
from os.path import exists

from luna.util.exceptions import IllegalArgumentError, ProcessingFailed
from luna.wrappers.chimera.add_charges import CHARGE_METHODS
from luna.util.file import is_file_valid


import logging
logger = logging.getLogger()


def add_charges(input_file, output_file, total_charge,
                charge_method="am1-bcc"):

    logger.debug("Charges will be added to the molecule at '%s'." % input_file)

    if not is_file_valid(input_file):
        raise OSError("File '%s' does not exist.")

    if not isinstance(output_file, str):
        raise IllegalArgumentError("'output_file' must be a string.")

    if charge_method not in CHARGE_METHODS:
        raise IllegalArgumentError("The charge method '%s' is not "
                                   "accepted by Chimera." % charge_method)

    script = str(Path(__file__).parent) + "/add_charges.py"
    command = (f"chimera --nogui --nostatus --script '{script} "
               f"-i \"{input_file}\" -o \"{output_file}\" "
               f"-nc {total_charge} -c {charge_method}'")

    p = Popen(command, shell=True, stdout=PIPE, stderr=PIPE)
    try:
        stdout, stderr = p.communicate()
    except Exception:
        p.kill()
        raise

    if exists(output_file):
        logger.debug("File '%s' created with success." % output_file)
    else:
        logger.error(stderr.decode())
        raise ProcessingFailed("Chimera could not add charges "
                               "to the provided molecule.")
