# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
###################################################################

Modifications included by Alexandre Fassio (alexandrefassio@dcc.ufmg.br).
    Date: 05/03/2018.

1) Inherit inhouse modifications. Package: MyBio.
2) Now it parses CONECT records.
3) Added an option to PDBParser to correct conflicts in atom names.
    These conflicts can occur mainly when hydrogen atoms are added by OpenBabel.
    In this situation, all hydrogens have the name "H".
4) Added an option to PDBParser to correct flags that were incorrectly set by OpenBabel.
    When OpenBabel adds Hydrogen atoms, it does not set HETATM for hydrogen atoms from heteroatoms.
5) Added an option to correct all OpenBabel errors.

Modifications included by Alexandre Fassio (alexandrefassio@dcc.ufmg.br).
    Date: 01/02/19

1) Now it recognizes different labels for waters: 'HOH', 'DOD', 'WAT', 'H2O', 'OH2'

Modifications included by Alexandre Fassio (alexandrefassio@dcc.ufmg.br).
    Date: 04/09/19

1) Modified the line where the residue name is captured. There, I applied a strip function to remove any whitespace characters.
    E.g.: Sodium id is in the format ' NA' and we want something like 'NA'.

Each line or block with modifications contain a MODBY tag.

###################################################################
"""

# Parser for PDB files.

from __future__ import print_function

import warnings

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use the PDB parser.")

from Bio.File import as_handle

# MODBY: Alexandre Fassio
# Inherit inhouse modifications. Package: MyBio.
from luna.MyBio.PDB.PDBExceptions import PDBConstructionException
from luna.MyBio.PDB.PDBExceptions import PDBConstructionWarning
from luna.MyBio.PDB.StructureBuilder import StructureBuilder
from luna.MyBio.PDB.parse_pdb_header import _parse_pdb_header_list


# MODBY: Alexandre Fassio
# Possible labels for water molecules.
# Ref: http://prody.csb.pitt.edu/manual/reference/atomic/flags.html.
# DOD: deutered water.
WATER_NAMES = ['HOH', 'DOD', 'WAT', 'H2O', 'OH', 'OH2', "O"]

# MODBY: Alexandre Fassio
# Replace empty chains for a default value.
DEFAULT_CHAIN_ID = "z"


# If PDB spec says "COLUMNS 18-20" this means line[17:20]

class PDBParser(object):
    """Parse a PDB file and return a Structure object."""

    def __init__(self, PERMISSIVE=True, get_header=False,
                 structure_builder=None, QUIET=False,
                 FIX_ATOM_NAME_CONFLICT=False,
                 FIX_EMPTY_CHAINS=False,
                 FIX_OBABEL_FLAGS=False,
                 FIX_ALL_OBABEL_ERRORS=False):
        """Create a PDBParser object.

        The PDB parser call a number of standard methods in an aggregated
        StructureBuilder object. Normally this object is instanciated by the
        PDBParser object itself, but if the user provides his/her own
        StructureBuilder object, the latter is used instead.

        Arguments:
         - PERMISSIVE - Evaluated as a Boolean. If false, exceptions in
           constructing the SMCRA data structure are fatal. If true (DEFAULT),
           the exceptions are caught, but some residues or atoms will be missing.
           THESE EXCEPTIONS ARE DUE TO PROBLEMS IN THE PDB FILE!.
         - structure_builder - an optional user implemented StructureBuilder class.
         - QUIET - Evaluated as a Boolean. If true, warnings issued in constructing
           the SMCRA data will be suppressed. If false (DEFAULT), they will be shown.
           These warnings might be indicative of problems in the PDB file!

        """
        if structure_builder is not None:
            self.structure_builder = structure_builder
        else:
            self.structure_builder = StructureBuilder()
        self.header = None
        # MODBY: Alexandre Fassio
        # Now it parses CONECT records.
        self.conects = None
        self.trailer = None
        self.line_counter = 0
        self.PERMISSIVE = bool(PERMISSIVE)
        self.QUIET = bool(QUIET)

        if FIX_ALL_OBABEL_ERRORS:
            FIX_ATOM_NAME_CONFLICT = True
            FIX_OBABEL_FLAGS = True
        # MODBY: Alexandre Fassio
        # Correct conflicts in atom names.
        self.FIX_ATOM_NAME_CONFLICT = FIX_ATOM_NAME_CONFLICT
        # MODBY: Alexandre Fassio
        # Fix empty chains.
        self.FIX_EMPTY_CHAINS = FIX_EMPTY_CHAINS
        # MODBY: Alexandre Fassio
        # Correct flags that were incorrectly set by OpenBabel.
        self.FIX_OBABEL_FLAGS = FIX_OBABEL_FLAGS
        # MODBY: Alexandre Fassio
        # Correct all OpenBabel errors.
        self.FIX_ALL_OBABEL_ERRORS = FIX_ALL_OBABEL_ERRORS

    # Public methods

    def get_structure(self, id, file):
        """Return the structure.

        Arguments:
         - id - string, the id that will be used for the structure
         - file - name of the PDB file OR an open filehandle

        """
        with warnings.catch_warnings():
            if self.QUIET:
                warnings.filterwarnings("ignore", category=PDBConstructionWarning)

            self.header = None
            self.trailer = None

            # Make a StructureBuilder instance (pass id of structure as parameter)
            # MODBY: Alexandre Fassio
            # Now it pass the PDB file as a parameter
            self.structure_builder.init_structure(id, file)

            with as_handle(file, mode='rU') as handle:
                self._parse(handle.readlines())

            self.structure_builder.set_header(self.header)

            # MODBY: Alexandre Fassio
            # Now it parses CONECT records.
            self.structure_builder.set_conects(self.conects)

            # Return the Structure instance
            structure = self.structure_builder.get_structure()

        return structure

    def get_structure_from_pdb_block(self, id, pdb_block):
        with warnings.catch_warnings():
            if self.QUIET:
                warnings.filterwarnings("ignore", category=PDBConstructionWarning)

            self.header = None
            self.trailer = None
            # Make a StructureBuilder instance (pass id of structure as parameter)
            self.structure_builder.init_structure(id)

            lines = pdb_block.split("\n")
            self._parse(lines)

            self.structure_builder.set_header(self.header)

            # MODBY: Alexandre Fassio
            # Now it parses CONECT records.
            self.structure_builder.set_conects(self.conects)

            # Return the Structure instance
            structure = self.structure_builder.get_structure()

        return structure

    def get_header(self):
        """Return the header."""
        return self.header

    def get_trailer(self):
        """Return the trailer."""
        return self.trailer

    # MODBY: Alexandre Fassio
    # Now it parses CONECT records.
    def get_conects(self):
        """Return conects."""
        return self.conects

    # Private methods

    def _parse(self, header_coords_trailer):
        """Parse the PDB file (PRIVATE)."""
        # Extract the header; return the rest of the file
        self.header, coords_trailer = self._get_header(header_coords_trailer)

        # Parse the atomic data; return the PDB file trailer
        self.trailer = self._parse_coordinates(coords_trailer)

        # MODBY: Alexandre Fassio
        # Now it parses CONECT records.
        self.conects = self._parse_conects(self.trailer)

    def _get_header(self, header_coords_trailer):
        """Get the header of the PDB file, return the rest (PRIVATE)."""
        structure_builder = self.structure_builder
        i = 0
        for i in range(0, len(header_coords_trailer)):
            structure_builder.set_line_counter(i + 1)
            line = header_coords_trailer[i]
            record_type = line[0:6]
            if record_type == "ATOM  " or record_type == "HETATM" or record_type == "MODEL ":
                break
        header = header_coords_trailer[0:i]
        # Return the rest of the coords+trailer for further processing
        self.line_counter = i
        coords_trailer = header_coords_trailer[i:]
        header_dict = _parse_pdb_header_list(header)
        return header_dict, coords_trailer

    def _parse_coordinates(self, coords_trailer):
        """Parse the atomic data in the PDB file (PRIVATE)."""
        local_line_counter = 0
        structure_builder = self.structure_builder
        current_model_id = 0
        # Flag we have an open model
        model_open = 0
        current_chain_id = None
        current_segid = None
        current_residue_id = None
        current_resname = None

        conflicts = {}
        for i in range(0, len(coords_trailer)):
            line = coords_trailer[i].rstrip('\n')
            record_type = line[0:6]
            global_line_counter = self.line_counter + local_line_counter + 1
            structure_builder.set_line_counter(global_line_counter)
            if record_type == "ATOM  " or record_type == "HETATM":
                # Initialize the Model - there was no explicit MODEL record
                if not model_open:
                    structure_builder.init_model(current_model_id)
                    current_model_id += 1
                    model_open = 1
                fullname = line[12:16]
                # get rid of whitespace in atom names
                split_list = fullname.split()
                if len(split_list) != 1:
                    # atom name has internal spaces, e.g. " N B ", so
                    # we do not strip spaces
                    name = fullname
                else:
                    # atom name is like " CA ", so we can strip spaces
                    name = split_list[0]
                altloc = line[16]
                # MODBY: Alexandre Fassio
                # Applying a strip() on the residue name because some of them may have external whitespace characters.
                resname = line[17:20].strip()
                chainid = line[21]

                # MODBY: Alexandre Fassio
                # Replace empty chains for a default value.
                if chainid.strip() == "" and self.FIX_EMPTY_CHAINS:
                    chainid = DEFAULT_CHAIN_ID

                try:
                    serial_number = int(line[6:11])
                except Exception:
                    serial_number = 0
                resseq = int(line[22:26].split()[0])  # sequence identifier
                icode = line[26]  # insertion code

                # MODBY: Alexandre Fassio
                # Correct flags that were incorrectly set by OpenBabel.
                # Recognize other water labels.
                if self.FIX_OBABEL_FLAGS:
                    # Incorrect flags are set only for added hydrogen atoms
                    if name == "H":
                        child_dict = structure_builder.model[chainid].child_dict
                        # If this hydrogen belongs to water molecule or to a ligand
                        if resname in WATER_NAMES or ("H_%s" % resname, resseq, icode) in child_dict:
                            record_type = "HETATM"

                if record_type == "HETATM":  # hetero atom flag
                    if resname in WATER_NAMES:
                        hetero_flag = "W"
                    else:
                        hetero_flag = "H"
                else:
                    hetero_flag = " "

                residue_id = (hetero_flag, resseq, icode)
                # atomic coordinates
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except Exception:
                    # Should we allow parsing to continue in permissive mode?
                    # If so, what coordinates should we default to?  Easier to abort!
                    raise PDBConstructionException("Invalid or missing coordinate(s) at line %i."
                                                   % global_line_counter)
                coord = numpy.array((x, y, z), "f")
                # occupancy & B factor
                try:
                    occupancy = float(line[54:60])
                except Exception:
                    self._handle_PDB_exception("Invalid or missing occupancy",
                                               global_line_counter)
                    occupancy = None  # Rather than arbitrary zero or one
                if occupancy is not None and occupancy < 0:
                    # TODO - Should this be an error in strict mode?
                    # self._handle_PDB_exception("Negative occupancy",
                    #                            global_line_counter)
                    # This uses fixed text so the warning occurs once only:
                    warnings.warn("Negative occupancy in one or more atoms", PDBConstructionWarning)
                try:
                    bfactor = float(line[60:66])
                except Exception:
                    self._handle_PDB_exception("Invalid or missing B factor",
                                               global_line_counter)
                    bfactor = 0.0  # The PDB use a default of zero if the data is missing
                segid = line[72:76]
                element = line[76:78].strip().upper()
                if current_segid != segid:
                    current_segid = segid
                    structure_builder.init_seg(current_segid)
                if current_chain_id != chainid:
                    current_chain_id = chainid
                    structure_builder.init_chain(current_chain_id)
                    current_residue_id = residue_id
                    current_resname = resname
                    try:
                        structure_builder.init_residue(resname, hetero_flag, resseq, icode)
                    except PDBConstructionException as message:
                        self._handle_PDB_exception(message, global_line_counter)
                elif current_residue_id != residue_id or current_resname != resname:
                    current_residue_id = residue_id
                    current_resname = resname
                    try:
                        structure_builder.init_residue(resname, hetero_flag, resseq, icode)
                    except PDBConstructionException as message:
                        self._handle_PDB_exception(message, global_line_counter)
                # init atom
                try:
                    # MODBY: Alexandre Fassio
                    # Correct conflicts in atom names.
                    if self.FIX_ATOM_NAME_CONFLICT:
                        if structure_builder.residue.has_id(name):
                            conflict_key = (structure_builder.residue, name)
                            if (conflict_key in conflicts):
                                atual_error = conflicts[conflict_key] + 1
                            else:
                                atual_error = 1

                            # Replace atom name.
                            name += str(atual_error)
                            # Replace fullname (atom name + spaces).
                            parts = [' '] * 4
                            replace_ini = 0 if len(name) == 4 else 1
                            replace_end = len(name) + 1
                            parts[replace_ini:replace_end] = name
                            fullname = ''.join(parts)
                            conflicts[conflict_key] = atual_error

                    structure_builder.init_atom(name, coord, bfactor, occupancy, altloc,
                                                fullname, serial_number, element)
                except PDBConstructionException as message:
                    self._handle_PDB_exception(message, global_line_counter)
            elif record_type == "ANISOU":
                anisou = [float(x) for x in (line[28:35], line[35:42], line[43:49],
                                             line[49:56], line[56:63], line[63:70])]
                # U's are scaled by 10^4
                anisou_array = (numpy.array(anisou, "f") / 10000.0).astype("f")
                structure_builder.set_anisou(anisou_array)
            elif record_type == "MODEL ":
                try:
                    serial_num = int(line[10:14])
                except Exception:
                    self._handle_PDB_exception("Invalid or missing model serial number",
                                               global_line_counter)
                    serial_num = 0
                structure_builder.init_model(current_model_id, serial_num)
                current_model_id += 1
                model_open = 1
                current_chain_id = None
                current_residue_id = None
            elif record_type == "END   " or record_type == "CONECT":
                # End of atomic data, return the trailer
                self.line_counter += local_line_counter
                return coords_trailer[local_line_counter:]
            elif record_type == "ENDMDL":
                model_open = 0
                current_chain_id = None
                current_residue_id = None
            elif record_type == "SIGUIJ":
                # standard deviation of anisotropic B factor
                siguij = [float(x) for x in (line[28:35], line[35:42], line[42:49],
                                             line[49:56], line[56:63], line[63:70])]
                # U sigma's are scaled by 10^4
                siguij_array = (numpy.array(siguij, "f") / 10000.0).astype("f")
                structure_builder.set_siguij(siguij_array)
            elif record_type == "SIGATM":
                # standard deviation of atomic positions
                sigatm = [float(x) for x in (line[30:38], line[38:45], line[46:54],
                                             line[54:60], line[60:66])]
                sigatm_array = numpy.array(sigatm, "f")
                structure_builder.set_sigatm(sigatm_array)
            local_line_counter += 1
        # EOF (does not end in END or CONECT)
        self.line_counter = self.line_counter + local_line_counter
        return []

    # MODBY: Alexandre Fassio
    # Now it parses CONECT records.
    def _parse_conects(self, trailer):
        """Parse the CONECT record data in the PDB file (PRIVATE)."""
        local_line_counter = 0
        structure_builder = self.structure_builder

        structure = structure_builder.structure
        all_serial_numbers = set([a.get_serial_number()
                                  for a in structure.get_atoms()])

        conects = {}
        local_line_counter = 0

        for i in range(0, len(trailer)):
            line = trailer[i].rstrip('\n')
            record_type = line[0:6]

            global_line_counter = self.line_counter + local_line_counter + 1
            structure_builder.set_line_counter(global_line_counter)

            local_line_counter += 1

            if record_type == "CONECT":
                serial_numbers = [line[6:11], line[11:16], line[16:21],
                                  line[21:26], line[26:31]]

                try:
                    orig_serial_number = int(serial_numbers[0].strip())
                except Exception:
                    warnings.warn("Serial number '%s' in CONECT record "
                                  "is invalid. Line %d."
                                  % (orig_serial_number, global_line_counter),
                                  PDBConstructionWarning)
                    continue
                if (orig_serial_number not in all_serial_numbers):
                    warnings.warn("Serial number '%s' in CONECT record "
                                  "does not exist. Line %d."
                                  % (orig_serial_number, global_line_counter),
                                  PDBConstructionWarning)
                    continue

                bonded_atoms = []
                for serial_number in serial_numbers[1:]:
                    serial_number = serial_number.strip()

                    if (serial_number is ''):
                        continue

                    try:
                        serial_number = int(serial_number)
                    except Exception:
                        warnings.warn("Serial number '%s' in CONECT record "
                                      "is invalid. Line %d."
                                      % (serial_number, global_line_counter),
                                      PDBConstructionWarning)
                        continue

                    if (serial_number not in all_serial_numbers):
                        warnings.warn("Serial number '%s' in CONECT record "
                                      "does not exist. Line %d."
                                      % (serial_number, global_line_counter),
                                      PDBConstructionWarning)
                        continue

                    bonded_atoms.append(serial_number)

                if (len(bonded_atoms) == 0):
                    warnings.warn("No bonded atom were found in CONECT record."
                                  " Line %d." % global_line_counter,
                                  PDBConstructionWarning)
                    continue
                else:
                    orig_bonded_atoms = (conects[orig_serial_number]
                                         if orig_serial_number in conects
                                         else [])

                    orig_bonded_atoms += bonded_atoms
                    conects[orig_serial_number] = orig_bonded_atoms

        return conects

    def _handle_PDB_exception(self, message, line_counter):
        """Handle exception (PRIVATE).

        This method catches an exception that occurs in the StructureBuilder
        object (if PERMISSIVE), or raises it again, this time adding the
        PDB line number to the error message.
        """
        message = "%s at line %i." % (message, line_counter)
        if self.PERMISSIVE:
            # just print a warning - some residues/atoms may be missing
            warnings.warn("PDBConstructionException: %s\n"
                          "Exception ignored.\n"
                          "Some atoms or residues may be missing in the data structure."
                          % message, PDBConstructionWarning)
        else:
            # exceptions are fatal - raise again with new message (including line nr)
            raise PDBConstructionException(message)
