from pybel import informats, outformats
from subprocess import Popen, PIPE

from util.exceptions import FileNotCreated, InvalidFileFormat
from util.default_values import OPENBABEL
from util.file import is_file_valid, try_validate_file, get_file_format


import logging
logger = logging.getLogger()


def mol_to_svg(infile, output, opt=None):
    """Depict a molecule as SVG using OpenBabel.

        @param infile: a file to be converted.
        @type pdb_code: string

        @param output: a path to the converted file.
        @type output: string

        @param opt: a set of depiction options. Check OpenBabel options.
        @type opt: dictionary
        @example opt: {"C": None, "e": None, "P": 500, "d": None}
    """
    logger.info("Trying to generate a 2D diagram to the file '%s'." % infile)

    try:
        try_validate_file(infile)

        opt_list = _get_options(opt, "x")
        args = ["obabel", infile, "-O", output] + opt_list

        p = Popen(args, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()

        svg_created = stderr.decode().strip() == "1 molecule converted"
        if (not is_file_valid(output) or not svg_created):
            raise FileNotCreated("SVG diagram for '%s' not created." % infile)
        else:
            logger.info("SVG diagram for '%s' created." % infile)
    except Exception as e:
        logger.exception(e)
        raise


def _get_options(opt, prefix=""):
    opt_list = []

    if (opt is not None):
        for key in opt:
            if len(key) > 1:
                opt_list.append("--%s%s" % (prefix, key))
            else:
                opt_list.append("-%s%s" % (prefix, key))

            if (opt[key] is not None):
                opt_list.append(str(opt[key]))

    return opt_list


def convert_molecule(infile, output, infile_format=None,
                     output_format=None, opt=None, openbabel=OPENBABEL):
    """Convert a molecule file to other format using OpenBabel.

        @param infile: a file to be converted.
        @type pdb_code: string

        @param output: a path to the converted file.
        @type output: string

        @param infile_format: the format of the infile file.
        @type infile_format: string

        @param output_format: the format of the infile file.
        @type output_format: string

        @param opt: a set of convertion options. Check OpenBabel options.
        @type opt: dictionary
        @example opt: {"p": 7, "AddPolarH": None, c": None, "error-level": 5}
    """
    logger.info("Trying to convert the file '%s' to '%s'." % (infile, output))

    try:
        # It raises an error if it is not valid.
        try_validate_file(infile)

        if infile_format is None:
            logger.warning("Input file format not defined. It will assume the format by the file extension.")
            infile_format = get_file_format(infile)

        if infile_format not in informats:
            raise InvalidFileFormat("Infile format '%s' does not exist." % infile_format)

        logger.info("Input format: %s" % infile_format)

        if output_format is None:
            logger.warning("Output file format not defined. It will assume the format by the file extension.")
            output_format = get_file_format(output, 1)

        if output_format not in outformats:
            raise InvalidFileFormat("Infile format '%s' does not exist." % output_format)

        logger.info("Output format: %s" % output_format)

        opt_list = _get_options(opt)
        args = [openbabel, "-i", infile_format, infile, "-o", output_format, "-O", output] + opt_list

        p = Popen(args, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()

        output_lines = stderr.decode().strip().split("\n")
        is_mol_converted = output_lines[-1] != "0 molecule converted"
        if not is_file_valid(output) or not is_mol_converted:
            raise FileNotCreated("File '%s' not converted to '%s'." % (infile, output))
        else:
            logger.info("File '%s' converted to '%s'." % (infile, output))
    except Exception as e:
        logger.exception(e)
        raise
