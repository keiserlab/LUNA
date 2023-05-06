from subprocess import Popen, PIPE, TimeoutExpired

from openbabel.pybel import informats, outformats

from luna.util.exceptions import (FileNotCreated, InvalidFileFormat,
                                  IllegalArgumentError, ProcessingFailed)
from luna.util.default_values import OPENBABEL
from luna.util.file import is_file_valid, validate_file, get_file_format

import logging
logger = logging.getLogger()


def _prep_opts(opts, prefix=""):
    opt_list = []
    if opts is not None:
        for key in opts:
            if key == "x" and prefix == "":
                opt_list.append("-%s%s" % (key, str(opts[key])))
            else:
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
        A set of depiction options. Check `Open Babel <http://openbabel.org/docs/dev/FileFormats/SVG_depiction.html>`_
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


def convert_molecule(mol_input, input_format=None,
                     output_file=None, output_format=None,
                     opts=None, openbabel=OPENBABEL):
    """Convert a molecular file to another format using Open Babel.

    Parameters
    ----------

    mol_input : str
        The pathname of a molecular file or a SMILES string.
    input_format : str, optional
        The molecular format of ``mol_input``.
        If not provided, the format will be defined by the file ``mol_input``
        extension. If the extension could not be identified, an
        `InvalidFileFormat` exception will be raised.
    output_file : str or None
        Save the converted molecule to this file.
        If None, return the output molecule as string.
    output_format : str, optional
        The molecular format of ``output_file``.
        If not provided, the format will be defined by the ``output_file`` when
        it is not None.
        Otherwise, it will raise an `InvalidFileFormat` exception.
    opts : dict
        A set of convertion options.
        Check `Open Babel <https://openbabel.org/docs/dev/Command-line_tools/babel.html>`_
        to discover which options are available.
    openbabel : str, optional
        The Open Babel binary location.
        If not provided, the default binary ('obabel') will be used.


    Returns
    -------
     : str or None
        Return the converted molecule as a string if ``output_file`` is None.
        Otherwise, return None.

    Raises
    ------
    InvalidFileFormat
        If the provided molecular formats are not accepted by Open Babel.

    Examples
    --------

    In this example, we will convert the molecule ZINC000007786517 from
    the format MOL to MOL2 and add hydrogens to it considering a pH of 7
    (option "p").

    >>> from luna.util.default_values import LUNA_PATH
    >>> from luna.wrappers.obabel import convert_molecule
    >>> convert_molecule(mol_input=f"{LUNA_PATH}/tutorial/inputs/\
ZINC000007786517.mol",
    ...                  output_file="example.mol2",
    ...                  opts={"p": 7})
    """
    logger.info("Molecule or file '%s' will be converted to format '%s'."
                % (mol_input, input_format))

    if is_file_valid(mol_input):
        if input_format is None:
            msg = ("Input file format not defined. "
                   "It will assume the format from the file extension.")
            logger.debug(msg)
            input_format = get_file_format(mol_input)

        if input_format not in informats:
            msg = "Input format '%s' does not exist." % input_format
            raise InvalidFileFormat(msg)

        logger.debug("Input format: %s" % input_format)

    else:
        if input_format != "smi":
            msg = ("File '%s' does not exist or is not a valid file."
                   % mol_input)
            raise OSError(msg)

    if isinstance(output_file, str):
        if output_format is None:
            msg = ("Output file format not defined. It will assume the format "
                   "by the file extension.")
            logger.debug(msg)
            output_format = get_file_format(output_file, 1)

        if output_format not in outformats:
            msg = "Output format '%s' does not exist." % output_format
            raise InvalidFileFormat(msg)

        logger.debug("Output format: %s" % output_format)

    elif output_file is not None:
        msg = "The 'output_file' should be a string or None."
        raise IllegalArgumentError(msg)

    if output_format is None:
        msg = "The output format could not be identified."
        raise IllegalArgumentError(msg)

    if input_format == "smi":
        input_list = [f'-:{mol_input}']
    else:
        input_list = ['-i', input_format, mol_input]

    if output_file is None:
        output_list = ['-o', output_format]
    else:
        output_list = ['-o', output_format, '-O', output_file]

    opts = opts or {}

    timeout = 30
    if "gen3d" in opts or "minimize" in opts:
        timeout = None

    opt_list = _prep_opts(opts)
    args = [openbabel] + input_list + output_list + opt_list

    p = Popen(args, stdout=PIPE, stderr=PIPE)
    try:
        stdout, stderr = p.communicate(timeout=timeout)
    except TimeoutExpired:
        p.kill()
        raise

    error_lines = stderr.decode().strip().split("\n")
    if (error_lines[-1] == "0 molecules converted"
            or error_lines[-1] == "0 molecule converted"):
        raise ProcessingFailed("The provided molecule could not be converted.")

    if output_file is None:
        return stdout.decode().strip()

    logger.debug("File '%s' created with success." % output_file)
