from util.file import (is_file_valid, try_validate_file, get_file_format)
from subprocess import (Popen, PIPE)
from util.exceptions import (FileNotCreated, InvalidFileFormat)

import logging
logger = logging.getLogger(__name__)


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
            if (len(key) > 1):
                opt_list.append("--%s%s" % (prefix, key))
            else:
                opt_list.append("-%s%s" % (prefix, key))

            if (opt[key] is not None):
                opt_list.append(str(opt[key]))

    return opt_list


def convert_molecule(infile, output, infile_format=None,
                     output_format=None, opt=None, openbabel='obabel'):
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

        if (infile_format is None):
            logger.info("Input file format not defined.")
            logger.info("It will try to figure out the format.")
            infile_format = get_file_format(infile)

        if (infile_format not in read_formats()):
            raise InvalidFileFormat("Infile format '%s' does not exist." % infile_format)

        logger.info("Input format: %s" % infile_format)

        if (output_format is None):
            logger.info("Output file format not defined.")
            logger.info("It will try to figure out the format.")
            output_format = get_file_format(output, 1)

        if (output_format not in write_formats()):
            raise InvalidFileFormat("Infile format '%s' does not exist." % output_format)

        logger.info("Output format: %s" % output_format)

        opt_list = _get_options(opt)
        args = [openbabel, "-i", infile_format, infile, "-o", output_format, "-O", output] + opt_list

        p = Popen(args, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()

        output_lines = stderr.decode().strip().split("\n")
        is_mol_converted = output_lines[-1] != "0 molecule converted"
        if (not is_file_valid(output) or not is_mol_converted):
            raise FileNotCreated("File '%s' not converted to '%s'." % (infile, output))
        else:
            logger.info("File '%s' converted to '%s'." % (infile, output))
    except Exception as e:
        logger.exception(e)
        raise


def write_formats():
    return {"acesin": "ACES input format",
            "adf": "ADF cartesian input format",
            "alc": "Alchemy format",
            "ascii": "ASCII format",
            "bgf": "MSI BGF format",
            "box": "Dock 3.5 Box format",
            "bs": "Ball and Stick format",
            "c3d1": "Chem3D Cartesian 1 format",
            "c3d2": "Chem3D Cartesian 2 format",
            "cac": "CAChe MolStruct format",
            "caccrt": "Cacao Cartesian format",
            "cache": "CAChe MolStruct format",
            "cacint": "Cacao Internal format",
            "can": "Canonical SMILES format",
            "cdxml": "ChemDraw CDXML format",
            "cht": "Chemtool format",
            "cif": "Crystallographic Information File",
            "ck": "ChemKin format",
            "cml": "Chemical Markup Language",
            "cmlr": "CML Reaction format",
            "com": "Gaussian 98/03 Input",
            "CONFIG": "DL-POLY CONFIG",
            "CONTCAR": "VASP format",
            "copy": "Copy raw text",
            "crk2d": "Chemical Resource Kit diagram(2D)",
            "crk3d": "Chemical Resource Kit 3D format",
            "csr": "Accelrys/MSI Quanta CSR format",
            "cssr": "CSD CSSR format",
            "ct": "ChemDraw Connection Table format",
            "cub": "Gaussian cube format",
            "cube": "Gaussian cube format",
            "dmol": "DMol3 coordinates format",
            "dx": "OpenDX cube format for APBS",
            "ent": "Protein Data Bank format",
            "fa": "FASTA format",
            "fasta": "FASTA format",
            "feat": "Feature format",
            "fh": "Fenske-Hall Z-Matrix format",
            "fhiaims": "FHIaims XYZ format",
            "fix": "SMILES FIX format",
            "fps": "FPS text fingerprint format (Dalke)",
            "fpt": "Fingerprint format",
            "fract": "Free Form Fractional format",
            "fs": "Fastsearch format",
            "fsa": "FASTA format",
            "gamin": "GAMESS Input",
            "gau": "Gaussian 98/03 Input",
            "gjc": "Gaussian 98/03 Input",
            "gjf": "Gaussian 98/03 Input",
            "gpr": "Ghemical format",
            "gr96": "GROMOS96 format",
            "gro": "GRO format",
            "gukin": "GAMESS-UK Input",
            "gukout": "GAMESS-UK Output",
            "gzmat": "Gaussian Z-Matrix Input",
            "hin": "HyperChem HIN format",
            "inchi": "InChI format",
            "inchikey": "InChIKey",
            "inp": "GAMESS Input",
            "jin": "Jaguar input format",
            "k": "Compare molecules using InChI",
            "lmpdat": "The LAMMPS data format",
            "mcdl": "MCDL format",
            "mcif": "Macromolecular Crystallographic Info",
            "mdl": "MDL MOL format",
            "ml2": "Sybyl Mol2 format",
            "mmcif": "Macromolecular Crystallographic Info",
            "mmd": "MacroModel format",
            "mmod": "MacroModel format",
            "mna": "Multilevel Neighborhoods of Atoms (MNA)",
            "mol": "MDL MOL format",
            "mol2": "Sybyl Mol2 format",
            "mold": "Molden format",
            "molden": "Molden format",
            "molf": "Molden format",
            "molreport": "Open Babel molecule report",
            "mop": "MOPAC Cartesian format",
            "mopcrt": "MOPAC Cartesian format",
            "mopin": "MOPAC Internal",
            "mp": "Molpro input format",
            "mpc": "MOPAC Cartesian format",
            "mpd": "MolPrint2D format",
            "mpqcin": "MPQC simplified input format",
            "mrv": "Chemical Markup Language",
            "msms": "M.F. Sanner's MSMS input format",
            "nul": "Outputs nothing",
            "nw": "NWChem input format",
            "outmol": "DMol3 coordinates format",
            "pcm": "PCModel Format",
            "pdb": "Protein Data Bank format",
            "pdbqt": "AutoDock PDQBT format",
            "png": "PNG 2D depiction",
            "POSCAR": "VASP format",
            "pov": "POV-Ray input format",
            "pqr": "PQR format",
            "pqs": "Parallel Quantum Solutions format",
            "qcin": "Q-Chem input format",
            "report": "Open Babel report format",
            "rsmi": "Reaction SMILES format",
            "rxn": "MDL RXN format",
            "sd": "MDL MOL format",
            "sdf": "MDL MOL format",
            "smi": "SMILES format",
            "smiles": "SMILES format",
            "svg": "SVG 2D depiction",
            "sy2": "Sybyl Mol2 format",
            "tdd": "Thermo format",
            "text": "Read and write raw text",
            "therm": "Thermo format",
            "tmol": "TurboMole Coordinate format",
            "txt": "Title format",
            "txyz": "Tinker XYZ format",
            "unixyz": "UniChem XYZ format",
            "VASP": "VASP format",
            "vmol": "ViewMol format",
            "xed": "XED format",
            "xyz": "XYZ cartesian coordinates format",
            "yob": "YASARA.org YOB format",
            "zin": "ZINDO input format"}


def read_formats():
    return {"abinit": "ABINIT Output Format",
            "acesout": "ACES output format",
            "acr": "ACR format",
            "adfout": "ADF output format",
            "alc": "Alchemy format",
            "arc": "Accelrys/MSI Biosym/Insight II CAR format",
            "axsf": "XCrySDen Structure Format",
            "bgf": "MSI BGF format",
            "box": "Dock 3.5 Box format",
            "bs": "Ball and Stick format",
            "c09out": "Crystal 09 output format",
            "c3d1": "Chem3D Cartesian 1 format",
            "c3d2": "Chem3D Cartesian 2 format",
            "caccrt": "Cacao Cartesian format",
            "can": "Canonical SMILES format",
            "car": "Accelrys/MSI Biosym/Insight II CAR format",
            "castep": "CASTEP format",
            "ccc": "CCC format",
            "cdx": "ChemDraw binary format",
            "cdxml": "ChemDraw CDXML format",
            "cif": "Crystallographic Information File",
            "ck": "ChemKin format",
            "cml": "Chemical Markup Language",
            "cmlr": "CML Reaction format",
            "CONFIG": "DL-POLY CONFIG",
            "CONTCAR": "VASP format",
            "crk2d": "Chemical Resource Kit diagram(2D)",
            "crk3d": "Chemical Resource Kit 3D format",
            "ct": "ChemDraw Connection Table format",
            "cub": "Gaussian cube format",
            "cube": "Gaussian cube format",
            "dat": "Generic Output file format",
            "dmol": "DMol3 coordinates format",
            "dx": "OpenDX cube format for APBS",
            "ent": "Protein Data Bank format",
            "fa": "FASTA format",
            "fasta": "FASTA format",
            "fch": "Gaussian formatted checkpoint file format",
            "fchk": "Gaussian formatted checkpoint file format",
            "fck": "Gaussian formatted checkpoint file format",
            "feat": "Feature format",
            "fhiaims": "FHIaims XYZ format",
            "fract": "Free Form Fractional format",
            "fs": "Fastsearch format",
            "fsa": "FASTA format",
            "g03": "Gaussian Output",
            "g09": "Gaussian Output",
            "g92": "Gaussian Output",
            "g94": "Gaussian Output",
            "g98": "Gaussian Output",
            "gal": "Gaussian Output",
            "gam": "GAMESS Output",
            "gamess": "GAMESS Output",
            "gamin": "GAMESS Input",
            "gamout": "GAMESS Output",
            "got": "GULP format",
            "gpr": "Ghemical format",
            "gro": "GRO format",
            "gukin": "GAMESS-UK Input",
            "gukout": "GAMESS-UK Output",
            "gzmat": "Gaussian Z-Matrix Input",
            "hin": "HyperChem HIN format",
            "HISTORY": "DL-POLY HISTORY",
            "inchi": "InChI format",
            "inp": "GAMESS Input",
            "ins": "ShelX format",
            "jout": "Jaguar output format",
            "log": "Generic Output file format",
            "mcdl": "MCDL format",
            "mcif": "Macromolecular Crystallographic Info",
            "mdl": "MDL MOL format",
            "ml2": "Sybyl Mol2 format",
            "mmcif": "Macromolecular Crystallographic Info",
            "mmd": "MacroModel format",
            "mmod": "MacroModel format",
            "mol": "MDL MOL format",
            "mol2": "Sybyl Mol2 format",
            "mold": "Molden format",
            "molden": "Molden format",
            "molf": "Molden format",
            "moo": "MOPAC Output format",
            "mop": "MOPAC Cartesian format",
            "mopcrt": "MOPAC Cartesian format",
            "mopin": "MOPAC Internal",
            "mopout": "MOPAC Output format",
            "mpc": "MOPAC Cartesian format",
            "mpo": "Molpro output format",
            "mpqc": "MPQC output format",
            "mrv": "Chemical Markup Language",
            "msi": "Accelrys/MSI Cerius II MSI format",
            "nwo": "NWChem output format",
            "out": "Generic Output file format",
            "outmol": "DMol3 coordinates format",
            "output": "Generic Output file format",
            "pc": "PubChem format",
            "pcm": "PCModel Format",
            "pdb": "Protein Data Bank format",
            "pdbqt": "AutoDock PDQBT format",
            "png": "PNG 2D depiction",
            "pos": "POS cartesian coordinates format",
            "POSCAR": "VASP format",
            "pqr": "PQR format",
            "pqs": "Parallel Quantum Solutions format",
            "prep": "Amber Prep format",
            "pwscf": "PWscf format",
            "qcout": "Q-Chem output format",
            "res": "ShelX format",
            "rsmi": "Reaction SMILES format",
            "rxn": "MDL RXN format",
            "sd": "MDL MOL format",
            "sdf": "MDL MOL format",
            "smi": "SMILES format",
            "smiles": "SMILES format",
            "sy2": "Sybyl Mol2 format",
            "t41": "ADF TAPE41 format",
            "tdd": "Thermo format",
            "text": "Read and write raw text",
            "therm": "Thermo format",
            "tmol": "TurboMole Coordinate format",
            "txt": "Title format",
            "txyz": "Tinker XYZ format",
            "unixyz": "UniChem XYZ format",
            "VASP": "VASP format",
            "vmol": "ViewMol format",
            "xml": "General XML format",
            "xsf": "XCrySDen Structure Format",
            "xtc": "XTC format",
            "xyz": "XYZ cartesian coordinates format",
            "yob": "YASARA.org YOB format"}
