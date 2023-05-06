from subprocess import Popen, PIPE
import tempfile
import re

from luna.wrappers.chemaxon.util import check_license
from luna.wrappers.base import MolWrapper
from luna.util.exceptions import ProcessingFailed
from luna.util.file import new_unique_filename

import logging
logger = logging.getLogger()


REGEX_MOL_FILE = re.compile(r' (V2000|V3000)$')


class TautomerGenerator:

    def __init__(self, ph=None, path_length=4,
                 protect_aromaticity=True,
                 protect_charge=True,
                 exclude_antiaromatic=True,
                 protect_double_bond_stereo=False,
                 protect_all_tetra_centers=False,
                 protect_labeled_tetra_centers=False,
                 protect_ester_groups=True,
                 min_percentage=20,
                 only_non_zeros=True,
                 mol_obj_type="rdkit",
                 skip_errors=False,
                 cxcalc="cxcalc"):

        self.ph = ph
        self.path_length = path_length
        self.protect_aromaticity = protect_aromaticity
        self.protect_charge = protect_charge
        self.exclude_antiaromatic = exclude_antiaromatic
        self.protect_double_bond_stereo = protect_double_bond_stereo
        self.protect_all_tetra_centers = protect_all_tetra_centers
        self.protect_labeled_tetra_centers = protect_labeled_tetra_centers
        self.protect_ester_groups = protect_ester_groups
        self.min_percentage = min_percentage
        self.only_non_zeros = only_non_zeros
        self.mol_obj_type = mol_obj_type
        self.skip_errors = skip_errors
        self.cxcalc = cxcalc

        self.percentages = []

    def _prep_opts(self):

        opt_list = []

        if self.ph is not None:
            opt_list += ["-H", str(self.ph)]

        # Path length
        opt_list += ["-l", str(self.path_length)]

        # Protect aromaticity
        opt_list += ["-a", str(self.protect_aromaticity)]

        # Protect charge
        opt_list += ["-C", str(self.protect_charge)]

        # Exclude antiaromatic compounds
        opt_list += ["-e", str(self.exclude_antiaromatic)]

        # Protect double bond stereo
        opt_list += ["-P", str(self.protect_double_bond_stereo)]

        # Protect all tetrahedral stereo
        opt_list += ["-T", str(self.protect_all_tetra_centers)]

        # Protect labeled tetrahedral stereo centers
        opt_list += ["-L", str(self.protect_labeled_tetra_centers)]

        # Protect ester groups
        opt_list += ["-E", str(self.protect_ester_groups)]

        return opt_list

    def _parse(self):
        self.tautomers = []
        self.percentages = []

        def apply_mol_format(lines):
            # Check if the MOL file contains a valid header comprising
            # a Title line, Program/file timestamp line, and a Comment line.
            # Thus, the Counts line must be in line 4 (index 3).
            # If it's not, add blank lines to match the format.
            # This amendment is necessary otherwise RDKit will crash.
            if not REGEX_MOL_FILE.search(lines[3]):
                logger.info("While parsing the file '%s', we found a "
                            "molecule starting at line #%d that does "
                            "not contain a valid header. "
                            "We will add empty lines to its header to "
                            "match the MOL file format."
                            % (self.output_file, mol_starts_at))

                counts_line_pos = None
                for i, line in enumerate(lines):
                    if REGEX_MOL_FILE.search(line):
                        counts_line_pos = i
                        break

                missing_lines = ["\n"] * (3 - counts_line_pos)
                lines = missing_lines + lines

            if lines[-1].strip() != "M  END":
                logger.info("While parsing the file '%s', we found a "
                            "molecule starting at line #%d that does not "
                            "contain the END line. We will add it to the "
                            "block end to match the MOL file format."
                            % (self.output_file, mol_starts_at))
                lines.append("M  END\n")
            return lines

        with open(self.output_file) as IN:
            n_mols = 0
            mol = []
            line_count = 0
            ignore_lines = False
            for line in IN:
                line_count += 1

                # Ignore new lines before a molecule block.
                if len(mol) == 0 and line == "\n":
                    continue

                # Save the line in which the new molecule starts.
                if len(mol) == 0:
                    mol_starts_at = line_count

                if line.startswith("M  END"):
                    mol.append(line)
                    # Fix the MOL block in cases where the header or end
                    # lines do not match the MOL file format.
                    mol = apply_mol_format(mol)

                    # After finding a molecule block, ignore any following
                    # lines until it finds a line "$$$$".
                    ignore_lines = True

                elif line.startswith(">  <TAUTOMER_DISTRIBUTION>"):
                    line = IN.readline()
                    line_count += 1
                    # readline() returns empty strings when EOF is reached.
                    if not line:
                        raise StopIteration

                    perc = float(line.strip())

                    # If the tautomer distribution is less than
                    # the minimum percentage or if it is zero,
                    # ignore the molecule.
                    ignore_mol = False
                    if ((self.min_percentage is not None
                            and perc < self.min_percentage)
                            or (self.only_non_zeros and perc == 0)):
                        ignore_mol = True

                    if not ignore_mol:
                        mol_obj = None
                        try:
                            # Create a new MolWrapper object
                            mol_obj = \
                                MolWrapper.from_mol_block("".join(mol), "sdf",
                                                          self.mol_obj_type)
                        except Exception:
                            if not self.skip_errors:
                                raise

                        if mol_obj is not None:
                            mol_id = "Tautomer %d" % (n_mols + 1)
                            mol_obj.set_name(mol_id)

                            self.tautomers.append(mol_obj)
                            self.percentages.append((mol_id, perc))

                    # Restart the molecule block.
                    mol = []

                elif line.startswith("$$$$") is False:
                    # Ignore lines util it finds a line "$$$$".
                    if ignore_lines:
                        continue
                    mol.append(line)
                else:
                    ignore_lines = False
                    n_mols += 1

    def _run_cxcalc(self, mol_input, output_file=None):
        opt_list = self._prep_opts()

        self.mol_input = mol_input

        if output_file is None:
            output_file = new_unique_filename(tempfile.gettempdir()) + ".sdf"
        self.output_file = output_file
        output = ["-o", output_file]

        plugin = "dominanttautomerdistribution"
        args = [self.cxcalc] + [mol_input] + output + [plugin] + opt_list

        p = Popen(args, stdout=PIPE, stderr=PIPE)
        try:
            stdout, stderr = p.communicate()
        except Exception:
            p.kill()
            raise

        error_lines = stderr.decode().strip().split("\n")

        if len(error_lines) == 1 and error_lines[0] == "":
            logger.debug("Tautomers for input '%s' created with success."
                         % mol_input)
        else:
            logger.error(stderr.decode())
            raise ProcessingFailed("cxcalc could not generate tautomers. "
                                   "Check the logs for more information.")

        self._parse()

    def generate(self, mol_input, output_file=None):
        self.tautomers = []

        check_license()

        self._run_cxcalc(mol_input, output_file)

        return self.tautomers
