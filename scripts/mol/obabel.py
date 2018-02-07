from file.validator import (is_file_valid,
                            try_validate_file)

import logging
logger = logging.getLogger(__name__)


def mol_2svg_obabel(infile, output, opt=None):
    """Depict a molecule as SVG using OpenBabel.

        @param infile: a file to be converted.
        @type pdb_code: string

        @param output: a path to the converted file.
        @type output: string

        @param opt: a set of depiction options. Check OpenBabel options.
        @type opt: dictionary
        @example opt: {"C": None, "e": None, "P": 500, "d": None}
    """
    from subprocess import Popen, PIPE
    from util.exceptions import FileNotCreated

    logger.info("Trying to generate a 2D diagram to the file '%s'." % infile)

    try:
        try_validate_file(infile)

        optList = _get_options(opt, "x")
        args = ["obabel", infile, "-O", output] + optList

        p = Popen(args, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()

        isSVGCreated = stderr.decode().strip() == "1 molecule converted"
        if (not is_file_valid(output) or not isSVGCreated):
            raise FileNotCreated("SVG diagram for '%s' not created." % infile)
        else:
            logger.info("SVG diagram for '%s' created." % infile)
    except Exception as e:
        logger.exception(e)
        raise


def _get_options(opt, prefix=""):
    optList = []

    if (opt is not None):
        for key in opt:
            if (len(key) > 1):
                optList.append("--%s%s" % (prefix, key))
            else:
                optList.append("-%s%s" % (prefix, key))

            if (opt[key] is not None):
                optList.append(str(opt[key]))

    return optList


def convert_molecule(infile, output, infileFormat=None,
                     outputFormat=None, opt=None):
    """Convert a molecule file to other format using OpenBabel.

        @param infile: a file to be converted.
        @type pdb_code: string

        @param output: a path to the converted file.
        @type output: string

        @param infileFormat: the format of the infile file.
        @type infileFormat: string

        @param outputFormat: the format of the infile file.
        @type outputFormat: string

        @param opt: a set of convertion options. Check OpenBabel options.
        @type opt: dictionary
        @example opt: {"p": 7, "c": None, "error-level": 5}
    """
    from subprocess import Popen, PIPE
    from file.util import get_file_format
    from util.exceptions import FileNotCreated

    logger.info("Trying to convert the file '%s' to '%s'." % (infile, output))

    try:
        try_validate_file(infile)

        if (infileFormat is None):
            logger.info("Input file format not defined.")
            logger.info("It will try to figure out the format.")
            infileFormat = get_file_format(infile)
        logger.info("Input format: %s" % infileFormat)

        if (outputFormat is None):
            logger.info("Output file format not defined.")
            logger.info("It will try to figure out the format.")
            outputFormat = get_file_format(output)
        logger.info("Output format: %s" % outputFormat)

        optList = _get_options(opt)
        args = ["obabel", "-i", infileFormat, infile,
                "-o", outputFormat, "-O", output] + optList

        p = Popen(args, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()

        isMolConverted = stderr.decode().strip() == "1 molecule converted"
        if (not is_file_valid(output) or not isMolConverted):
            raise FileNotCreated("File '%s' not converted to '%s'."
                                 % (infile, output))
        else:
            logger.info("File '%s' converted to '%s'." % (infile, output))
    except Exception as e:
        logger.exception(e)
        raise
