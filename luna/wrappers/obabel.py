from subprocess import Popen, PIPE, TimeoutExpired

from openbabel.pybel import informats, outformats

from luna.util.exceptions import FileNotCreated, InvalidFileFormat
from luna.util.default_values import OPENBABEL
from luna.util.file import is_file_valid, validate_file, get_file_format

import logging
logger = logging.getLogger()


def _prep_opts(opts, prefix=""):
    opt_list = []
    if opts is not None:
        for key in opts:
            if len(key) > 1:
                opt_list.append("--%s%s" % (prefix, key))
            else:
                opt_list.append("-%s%s" % (prefix, key))

            if opts[key] is not None:
                opt_list.append(str(opts[key]))
    return opt_list


def mol_to_svg(infile, output, opts=None):
    """Depict a molecule as SVG using Open Babel.

    Parameters
    ----------
    infile : str
        The pathname of a molecular file.
    output : str
        Save the SVG to this file.
    opts : dict
        A set of depiction options. Check `Open Babel <http://openbabel.org/docs/dev/FileFormats/SVG_depiction.html>`_ \
        to discover which options are available.

    Examples
    --------

    In this example, we will depict the molecule ZINC000007786517 as SVG.
    The following options will be used:
        * C: do not draw terminal C (and attached H) explicitly;
        * d: do not display the molecule name;
        * P: image size in pixels.

    >>> from luna.util.default_values import LUNA_PATH
    >>> from luna.wrappers.obabel import mol_to_svg
    >>> mol_to_svg(infile=f"{LUNA_PATH}/tutorial/inputs/ZINC000007786517.mol",
    ...            output="example.svg",
    ...            opts={"C": None, "d": None, "P": 500})

    """

    logger.debug("Trying to generate a 2D diagram to the file '%s'." % infile)

    try:
        validate_file(infile)

        opt_list = _prep_opts(opts, "x")
        args = ["obabel", infile, "-O", output] + opt_list

        p = Popen(args, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()

        svg_created = stderr.decode().strip() == "1 molecule converted"
        if not is_file_valid(output) or not svg_created:
            raise FileNotCreated("SVG diagram for '%s' not created." % infile)
        else:
            logger.debug("SVG diagram for '%s' created." % infile)
    except Exception as e:
        logger.exception(e)
        raise


def convert_molecule(infile, output, infile_format=None,
                     output_format=None, opts=None, openbabel=OPENBABEL):
    """Convert a molecular file to another format using Open Babel.

    Parameters
    ----------

    infile : str
        The pathname of the molecular file to be converted.
    output : str
        Save the converted molecule to this file.
    infile_format : str, optional
        The molecular format of ``infile``.
        If not provided, the format will be defined by the file extension.
    output_format : str, optional
        The molecular format of ``output``.
        If not provided, the format will be defined by the file extension.
    opts : dict
        A set of convertion options. Check `Open Babel <https://openbabel.org/docs/dev/Command-line_tools/babel.html>`_ to discover which options are available.
    openbabel : str, optional
        The Open Babel binary location. If not provided, the default binary ('obabel') will be used.

    Raises
    ------
    InvalidFileFormat
        If the provided molecular formats are not accepted by Open Babel.

    Examples
    --------

    In this example, we will convert the molecule ZINC000007786517 from the format MOL to MOL2
    and add hydrogens to it considering a pH of 7 (option "p").

    >>> from luna.util.default_values import LUNA_PATH
    >>> from luna.wrappers.obabel import convert_molecule
    >>> convert_molecule(infile=f"{LUNA_PATH}/tutorial/inputs/ZINC000007786517.mol",
    ...                  output="example.mol2",
    ...                  opts={"p": 7})
    """
    logger.debug("Trying to convert the file '%s' to '%s'." % (infile, output))

    try:
        # It raises an error if it is not valid.
        validate_file(infile)

        if infile_format is None:
            logger.debug("Input file format not defined. It will assume the format from the file extension.")
            infile_format = get_file_format(infile)

        if infile_format not in informats:
            raise InvalidFileFormat("Infile format '%s' does not exist." % infile_format)

        logger.debug("Input format: %s" % infile_format)

        if output_format is None:
            logger.debug("Output file format not defined. It will assume the format by the file extension.")
            output_format = get_file_format(output, 1)

        if output_format not in outformats:
            raise InvalidFileFormat("Infile format '%s' does not exist." % output_format)

        logger.debug("Output format: %s" % output_format)

        opt_list = _prep_opts(opts)
        args = [openbabel, "-i", infile_format, infile, "-o", output_format, "-O", output] + opt_list

        p = Popen(args, stdout=PIPE, stderr=PIPE)
        try:
            stdout, stderr = p.communicate(timeout=30)
        except TimeoutExpired:
            p.kill()
            raise

        output_lines = stderr.decode().strip().split("\n")
        is_mol_converted = output_lines[-1] != "0 molecule converted"
        if not is_file_valid(output) or not is_mol_converted:
            raise FileNotCreated("File '%s' not converted to '%s'." % (infile, output))
        else:
            logger.debug("File '%s' converted to '%s'." % (infile, output))
    except Exception:
        raise
