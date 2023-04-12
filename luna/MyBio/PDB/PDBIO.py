# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

###################################################################

# Modifications included by Alexandre Fassio (alexandrefassio@dcc.ufmg.br).
# Date: 08/03/2018.

# 1) Inherit inhouse modifications. Package: MyBio.
# 2) StructureBuilder now parses CONECT records.
# 3) set_structure now creates a new Structure object with the same Id as the original structure.

# Each line or block with modifications contain a MODBY tag.

###################################################################

"""Output of PDB files."""

from Bio._py3k import basestring

# MODBY: Alexandre Fassio
# Inherit inhouse modifications. Package: MyBio.
from luna.MyBio.PDB.StructureBuilder import StructureBuilder  # To allow saving chains, residues, etc..
from Bio.Data.IUPACData import atom_weights  # Allowed Elements

_ATOM_FORMAT_STRING = "%s%5i %-4s%c%3s %c%4i%c   %8.3f%8.3f%8.3f%s%6.2f      %4s%2s%2s\n"

_CONECT_FORMAT_STRING = "CONECT%5s%5s%5s%5s%5s\n"


class Select(object):
    """Select everything fo PDB output (for use as a bas class).

    Default selection (everything) during writing - can be used as base class
    to implement selective output. This selects which entities will be written out.
    """

    def __repr__(self):
        return "<Select all>"

    def accept_model(self, model):
        """Overload this to reject models for output."""
        return 1

    def accept_chain(self, chain):
        """Overload this to reject chains for output."""
        return 1

    def accept_residue(self, residue):
        """Overload this to reject residues for output."""
        return 1

    def accept_atom(self, atom):
        """Overload this to reject atoms for output."""
        return 1


class PDBIO(object):
    """Write a Structure object (or a subset of a Structure object) as a string.

    Example:

        >>> p=PDBParser()
        >>> s=p.get_structure("1fat", "1fat.pdb")
        >>> io=PDBIO()
        >>> io.set_structure(s)
        >>> io.save("out.pdb")

    """

    def __init__(self, use_model_flag=0):
        """Creat the PDBIO object.

        @param use_model_flag: if 1, force use of the MODEL record in output.
        @type use_model_flag: int
        """
        self.use_model_flag = use_model_flag

    # private mathods

    def _get_atom_line(self, atom, hetfield, segid, atom_number, resname,
                       resseq, icode, chain_id, charge="  "):
        """Returns an ATOM PDB string (PRIVATE)."""
        if hetfield != " ":
            record_type = "HETATM"
        else:
            record_type = "ATOM  "

        if atom.element:
            element = atom.element.strip().upper()
            if element.capitalize() not in atom_weights:
                raise ValueError("Unrecognised element %r" % atom.element)
            element = element.rjust(2)
        else:
            element = "  "

        name = atom.get_fullname().strip()

        # Pad atom name if:
        #     - smaller than 4 characters
        # AND - is not C, N, O, S, H, F, P, ..., one letter elements
        # AND - first character is NOT numeric (funky hydrogen naming rules)
        if len(name) < 4 and name[:1].isalpha() and len(element.strip()) < 2:
            name = " " + name

        altloc = atom.get_altloc()
        x, y, z = atom.get_coord()
        bfactor = atom.get_bfactor()
        occupancy = atom.get_occupancy()
        try:
            occupancy_str = "%6.2f" % occupancy
        except TypeError:
            if occupancy is None:
                occupancy_str = " " * 6
                import warnings
                from Bio import BiopythonWarning
                warnings.warn("Missing occupancy in atom %s written as blank" %
                              repr(atom.get_full_id()), BiopythonWarning)
            else:
                raise TypeError("Invalid occupancy %r in atom %r"
                                % (occupancy, atom.get_full_id()))

        args = (record_type, atom_number, name, altloc, resname, chain_id,
                resseq, icode, x, y, z, occupancy_str, bfactor, segid,
                element, charge)
        return _ATOM_FORMAT_STRING % args

    # Public methods

    def set_structure(self, pdb_object):
        # MODBY: Alexandre Fassio
        # Get the Structure object.
        parent_struct = pdb_object.get_parent_by_level(level='S')
        # Save a copy of the CONECT records from structure.
        conects = parent_struct.conects

        # Check what the user is providing and build a structure appropriately
        if pdb_object.level == "S":
            structure = pdb_object
        else:
            sb = StructureBuilder()
            # MODBY: Alexandre
            # Keeps the id from the original structure.
            # MODBY: Alexandre
            # Now it pass the PDB file as a parameter
            sb.init_structure(parent_struct.id, parent_struct.pdb_file)
            sb.init_seg(' ')
            # Build parts as necessary
            if pdb_object.level == "M":
                sb.structure.add(pdb_object)
                self.structure = sb.structure
            else:
                sb.init_model(0)
                if pdb_object.level == "C":
                    sb.structure[0].add(pdb_object)
                else:
                    sb.init_chain('A')
                    if pdb_object.level == "R":
                        parent_id = pdb_object.parent.id
                        try:
                            sb.structure[0]['A'].id = parent_id
                        except Exception:
                            pass

                        # MODBY: Alexandre Fassio
                        # The old code (sb.structure[0]['A'].add(pdb_object)) used to try to access a chain whose id is always 'A'.
                        # However, when the informed residue has a different chain id, this access will raise a KeyError exception.
                        # Thus, we should always access the new chain id by using the parent_id variable.
                        sb.structure[0][parent_id].add(pdb_object)
                    else:
                        # Atom
                        sb.init_residue('DUM', ' ', 1, ' ')
                        parent_id = pdb_object.parent.parent.id
                        try:
                            sb.structure[0]['A'].id = parent_id
                        except Exception:
                            pass

                        # MODBY: Alexandre Fassio
                        # The old code (sb.structure[0]['A'].child_list[0].add(pdb_object)) used to try to access a chain whose id
                        # is always 'A'. However, when the informed atom has a different chain id, this access will raise a
                        # KeyError exception. Thus, we should always access the new chain id by using the parent_id variable.
                        sb.structure[0][parent_id].child_list[0].add(pdb_object)

            # Return structure
            structure = sb.structure

        # MODBY: Alexandre Fassio
        # Set the copied CONECTs to the new Structure object.
        structure.conects = conects
        self.structure = structure

    def _get_list(self, entity, sort=False):
        # S -> M
        # M -> C
        if entity.get_level() in ["S", "M"]:
            the_list = entity.get_list()

            # Sort Models/Chains by id.
            if sort:
                return sorted(the_list, key=lambda x: x.id)

        # C -> R
        elif entity.get_level() == "C":
            the_list = entity.get_unpacked_list()

            # Sort residues by index, i.e., by their order in a PDB chain.
            if sort:
                return sorted(the_list, key=lambda x: x.idx)

        # R -> A
        #
        #   Always keep atoms order as in PDB.
        #
        elif entity.get_level() == "R":
            the_list = entity.get_unpacked_list()

        return the_list

    def save(self, file, select=Select(), write_conects=False,
             write_end=True, preserve_atom_numbering=False, sort=False):
        """
        @param select: selects which entities will be written.
        @type select: object

        Typically select is a subclass of L{Select}, it should
        have the following methods:

         - accept_model(model)
         - accept_chain(chain)
         - accept_residue(residue)
         - accept_atom(atom)

        These methods should return 1 if the entity is to be
        written out, 0 otherwise.

        Typically select is a subclass of L{Select}.
        """
        get_atom_line = self._get_atom_line
        if isinstance(file, basestring):
            fp = open(file, "w")
            close_file = 1
        else:
            # filehandle, I hope :-)
            fp = file
            close_file = 0

        # MODBY: Alexandre Fassio
        # Atom serial numbers to save on the PDB, i.e., not filtered atoms.
        valid_serial_numbers = set()
        serial_number_mapping = {}

        # multiple models?
        if len(self.structure) > 1 or self.use_model_flag:
            model_flag = 1
        else:
            model_flag = 0

        for model in self._get_list(self.structure, sort):
            if not select.accept_model(model):
                continue
            # necessary for ENDMDL
            # do not write ENDMDL if no residues were written
            # for this model
            model_residues_written = 0
            if not preserve_atom_numbering:
                atom_number = 1
            if model_flag:
                fp.write("MODEL      %s\n" % model.serial_num)
            for chain in self._get_list(model, sort):
                if not select.accept_chain(chain):
                    continue
                chain_id = chain.get_id()
                # necessary for TER
                # do not write TER if no residues were written
                # for this chain
                chain_residues_written = 0
                for residue in self._get_list(chain, sort):
                    if not select.accept_residue(residue):
                        continue
                    hetfield, resseq, icode = residue.get_id()
                    resname = residue.get_resname()
                    segid = residue.get_segid()
                    for atom in self._get_list(residue, sort):
                        if select.accept_atom(atom):
                            chain_residues_written = 1
                            model_residues_written = 1

                            if preserve_atom_numbering:
                                atom_number = atom.get_serial_number()

                            # MODBY: Alexandre Fassio
                            # Create a mapping between the current serial number and the new serial number.
                            # If preserve_atom_numbering is set to True, this mapping is irrelevant as both
                            # serial numbers will be the same.
                            serial_number_mapping[atom.get_serial_number()] = atom_number

                            s = get_atom_line(atom, hetfield, segid,
                                              atom_number, resname,
                                              resseq, icode, chain_id)

                            # MODBY: Alexandre Fassio
                            # This atom will not be filtered.
                            valid_serial_numbers.add(atom_number)
                            fp.write(s)
                            if not preserve_atom_numbering:
                                atom_number += 1

                if chain_residues_written:
                    # The serial number of the TER record is always one number greater than the serial number of the ATOM/HETATM
                    # preceding the TER.
                    last_serial_num = atom_number
                    if not preserve_atom_numbering:
                        # Now, increment the current serial number in case there is a new chain.
                        atom_number += 1

                    fp.write("TER   %5i      %3s %c%4i%c                                                      \n"
                             % (last_serial_num, resname, chain_id, resseq, icode))

            if model_flag and model_residues_written:
                fp.write("ENDMDL\n")

        # MODBY: Alexandre Fassio
        # Print CONECT records
        if write_conects:
            conects = self.structure.conects
            for serial_number in sorted(conects):
                # Substitutes the old serial number by the new serial number
                new_serial_number = serial_number_mapping.get(serial_number, None)

                if new_serial_number not in valid_serial_numbers:
                    continue

                # It may return a list of lists as some atoms may have more than one CONECT line.
                bonded_atoms = conects[serial_number]
                max_num_fields = 4
                bonded_atoms_sets = [bonded_atoms[i:i + max_num_fields]
                                     for i in range(0, len(bonded_atoms),
                                                    max_num_fields)]

                for bonded_atoms in bonded_atoms_sets:
                    valid_bonded_atoms = []
                    for tmp_serial_number in bonded_atoms:
                        # Substitutes the old serial number by the new serial number
                        tmp_serial_number = serial_number_mapping.get(tmp_serial_number, None)

                        if tmp_serial_number in valid_serial_numbers:
                            valid_bonded_atoms.append(tmp_serial_number)

                    if len(valid_bonded_atoms) == 0:
                        continue

                    valid_bonded_atoms = [str(x) for x in valid_bonded_atoms]
                    missing_values = max_num_fields - len(valid_bonded_atoms)
                    valid_bonded_atoms += [''] * missing_values

                    record = _CONECT_FORMAT_STRING % (str(new_serial_number),
                                                      *valid_bonded_atoms)

                    fp.write(record)
        if write_end:
            fp.write('END\n')
        if close_file:
            fp.close()


if __name__ == "__main__":

    from luna.MyBio.PDB.PDBParser import PDBParser

    import sys

    p = PDBParser(PERMISSIVE=True)

    s = p.get_structure("test", sys.argv[1])

    io = PDBIO()
    io.set_structure(s)
    io.save("out1.pdb")

    with open("out2.pdb", "w") as fp:
        s1 = p.get_structure("test1", sys.argv[1])
        s2 = p.get_structure("test2", sys.argv[2])
        io = PDBIO(1)
        io.set_structure(s1)
        io.save(fp)
        io.set_structure(s2)
        io.save(fp, write_end=1)
