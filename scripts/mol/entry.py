import re
import logging
from operator import xor
from os.path import exists


from rdkit.Chem import MolToPDBBlock
from pybel import readfile
from pybel import informats as OB_FORMATS


from mol.wrappers.rdkit import RDKIT_FORMATS, read_multimol_file
from mol.interaction.filter import InteractionFilter
from util.default_values import ACCEPTED_MOL_OBJ_TYPES, ENTRY_SEPARATOR
from util.file import get_file_format
from util.exceptions import InvalidEntry, IllegalArgumentError, MoleculeObjectError, MoleculeNotFoundError
from MyBio.PDB.PDBParser import PDBParser
from MyBio.PDB.Entity import Entity


logger = logging.getLogger()

# Source: https://richjenks.com/filename-regex/
FILENAME_REGEX = r"(?!.{256,})(?!(aux|clock\$|con|nul|prn|com[1-9]|lpt[1-9])(?:$|\.))[^ ][ \.\w\-$()+=[\];#@~,&']+[^\. ]"

PDB_PARSER = PDBParser(PERMISSIVE=True, QUIET=True, FIX_ATOM_NAME_CONFLICT=False, FIX_OBABEL_FLAGS=False)

REGEX_RESNUM_ICODE = re.compile(r'^(\d+)([a-zA-z]?)$')

PCI_ENTRY_REGEX = re.compile(r'^%s:\w:\w{1,3}:\-?\d{1,4}[a-zA-z]?$' % FILENAME_REGEX)
PPI_ENTRY_REGEX = re.compile(r'^%s:\w$' % FILENAME_REGEX)


class Entry:

    def __init__(self, pdb_id, chain_id, comp_name=None, comp_num=None, comp_icode=None, is_hetatm=True, sep=ENTRY_SEPARATOR,
                 inter_filter=None):

        if xor(comp_name is None, comp_num is None):
            raise IllegalArgumentError("You tried to define a compound, so you must inform its name and number.")

        if comp_num is not None:
            try:
                assert float(comp_num).is_integer()
                comp_num = int(comp_num)
            except (ValueError, AssertionError):
                raise IllegalArgumentError("The informed residue number '%s' is invalid. It must be an integer." % str(comp_num))

        if comp_icode is not None:
            comp_icode = str(comp_icode)
            if comp_icode.isdigit() or len(comp_icode) > 1:
                raise IllegalArgumentError("The informed residue icode '%s' is invalid. It must be a character." % str(comp_icode))

        if inter_filter is not None and isinstance(inter_filter, InteractionFilter) is False:
            raise IllegalArgumentError("The informed interaction filter must be an instance of '%s'." % InteractionFilter)


        self.pdb_id = pdb_id
        self.chain_id = chain_id
        self.comp_name = comp_name
        self.comp_num = comp_num
        self._comp_icode = comp_icode
        self.is_hetatm = is_hetatm
        self.sep = sep
        self.inter_filter = inter_filter

        if not self.is_valid():
            raise InvalidEntry("Entry '%s' does not match the PDB format." % self.to_string())

    @classmethod
    def from_string(cls, entry_str, is_hetatm=True, sep=ENTRY_SEPARATOR, inter_filter=None):
        entries = entry_str.split(sep)
        if len(entries) >= 2 and len(entries) <= 4:
            if len(entries) == 4:
                # Separate ligand number from insertion code.
                matched = REGEX_RESNUM_ICODE.match(entries[3])
                if matched:
                    comp_num = matched.group(1)
                    icode = None if matched.group(2) == "" else matched.group(2)
                    entries = entries[0:3] + [comp_num, icode]
                else:
                    raise IllegalArgumentError("The field residue/ligand number and its insertion code (if applicable) '%s' is invalid. "
                                               "It must be an integer followed by one insertion code character when applicable."
                                               % entries[3])
            return cls(*entries, is_hetatm=is_hetatm, sep=sep, inter_filter=inter_filter)
        else:
            raise IllegalArgumentError("The number of fields in the informed string '%s' is incorrect. A valid string must contain "
                                       "two obligatory fields (PDB and chain id) and may contain two optional fields (residue name "
                                       "and residue number followed by its insertion code when applicable)." % entry_str)

    @property
    def comp_icode(self):
        if isinstance(self._comp_icode, str):
            return self._comp_icode
        else:
            return ' '

    @property
    def full_id(self):
        return (self.pdb_id, self.chain_id, self.comp_name, self.comp_num, self.comp_icode)

    def to_string(self, sep=None):
        full_id = self.full_id

        # An entry object will always have a PDB and chain id.
        entry = list(full_id[0:2])
        # If it contains additional information about the residue it will also include them.
        if full_id[2] and full_id[3]:
            comp_name = str(full_id[2]) if full_id[2] else ""

            comp_num_and_icode = str(full_id[3]) if full_id[3] else ""
            comp_num_and_icode += str(full_id[4]) if full_id[4].strip() else ""
            entry += [comp_name, comp_num_and_icode]

        sep = sep or self.sep

        return sep.join(entry)

    def is_valid(self):
        full_id = self.full_id

        # If it contains additional information about the residue it will also include them.
        if full_id[2] and full_id[3]:
            regex = PCI_ENTRY_REGEX
        else:
            regex = PPI_ENTRY_REGEX

        return regex.match(self.to_string(":")) is not None

    def get_biopython_key(self):
        if self.comp_name and self.comp_num:
            if self.comp_name == 'HOH' or self.comp_name == 'WAT':
                return ('W', self.comp_num, self.comp_icode)
            elif self.is_hetatm:
                return ('H_%s' % self.comp_name, self.comp_num, self.comp_icode)
            else:
                return (' ', self.comp_num, self.comp_icode)

        return self.chain_id

    def __repr__(self):
        return '<%s: %s>' % (self.__class__.__name__, self.to_string(self.sep))


class DBEntry(Entry):

    def __init__(self, ligand_entry_id, pdb_id, chain_id, comp_name=None, comp_num=None, comp_icode=None):

        self.id = ligand_entry_id
        super().__init__(pdb_id, chain_id, comp_name, comp_num, comp_icode)


class CompoundEntry(Entry):

    def __init__(self, pdb_id, chain_id, comp_name, comp_num, comp_icode=None, sep=ENTRY_SEPARATOR, inter_filter=None):

        super().__init__(pdb_id, chain_id, comp_name, comp_num, comp_icode, is_hetatm=True, sep=sep, inter_filter=inter_filter)

    @classmethod
    def from_string(cls, entry_str, sep=ENTRY_SEPARATOR, inter_filter=None):
        entries = entry_str.split(sep)
        if len(entries) == 4:
            # Separate ligand number from insertion code.
            matched = REGEX_RESNUM_ICODE.match(entries[3])
            if matched:
                comp_num = matched.group(1)
                icode = None if matched.group(2) == "" else matched.group(2)
                entries = entries[0:3] + [comp_num, icode]
            else:
                raise IllegalArgumentError("The field residue/ligand number and its insertion code (if applicable) '%s' is invalid. "
                                           "It must be an integer followed by one insertion code character when applicable." % entries[3])
            return cls(*entries, sep=sep, inter_filter=inter_filter)
        else:
            raise IllegalArgumentError("The number of fields in the informed string '%s' is incorrect. A valid ligand entry must contain "
                                       "four obligatory fields: PDB, chain id, residue name, and residue number followed by its insertion "
                                       "code when applicable)." % entry_str)


class ChainEntry(Entry):

    def __init__(self, pdb_id, chain_id, sep=ENTRY_SEPARATOR, inter_filter=None):
        super().__init__(pdb_id, chain_id, is_hetatm=False, sep=sep, inter_filter=inter_filter)

    @classmethod
    def from_string(cls, entry_str, sep=ENTRY_SEPARATOR, inter_filter=None):
        entries = entry_str.split(sep)
        if len(entries) == 2:
            return cls(*entries, sep=sep, inter_filter=inter_filter)
        else:
            raise IllegalArgumentError("The number of fields in the informed string '%s' is incorrect. A valid string must contain "
                                       "two obligatory fields: PDB and chain id." % entry_str)


class MolEntry(Entry):

    def __init__(self, pdb_id, mol_id, mol_file=None, mol_file_ext=None, mol_obj_type='rdkit', sep=ENTRY_SEPARATOR, inter_filter=None):

        if mol_obj_type not in ACCEPTED_MOL_OBJ_TYPES:
            raise IllegalArgumentError("Objects of type '%s' are not currently accepted. "
                                       "The available options are: %s." % (mol_obj_type, ", ".join(ACCEPTED_MOL_OBJ_TYPES)))

        self.mol_id = mol_id
        self.mol_file = mol_file
        self.mol_file_ext = mol_file_ext
        self.mol_obj_type = mol_obj_type
        self._mol_obj = None

        super().__init__(pdb_id, "z", "LIG", 9999, sep=sep, inter_filter=inter_filter)

    @property
    def full_id(self):
        return (self.pdb_id, self.mol_id)

    @property
    def mol_obj(self):
        if self._mol_obj is None:
            self.load_mol_from_file()

        return self._mol_obj

    @mol_obj.setter
    def mol_obj(self, obj):
        self._mol_obj = obj

    def is_mol_obj_loaded(self):
        return self._mol_obj is not None

    def load_mol_from_file(self):
        logger.info("It will try to load the molecule '%s'." % self.mol_id)

        if self.mol_file is None:
            raise IllegalArgumentError("It cannot load the molecule as no molecular file was provided.")

        ext = self.mol_file_ext or get_file_format(self.mol_file)

        available_formats = OB_FORMATS if self.mol_obj_type == "openbabel" else RDKIT_FORMATS
        tool = "Open Babel" if self.mol_obj_type == "openbabel" else "RDKit"
        if ext not in available_formats:
            raise IllegalArgumentError("Extension '%s' informed or assumed from the filename is not a format "
                                       "recognized by %s." % (ext, tool))

        if not exists(self.mol_file):
            raise FileNotFoundError("The file '%s' was not found." % self.mol_file)

        try:
            if self.mol_obj_type == "openbabel":
                for ob_mol in readfile(ext, self.mol_file):
                    if self.mol_id == ob_mol.OBMol.GetTitle():
                        self._mol_obj = ob_mol
                        break
            else:
                for rdk_mol in read_multimol_file(self.mol_file, mol_file_format=ext, targets=[self.mol_id], removeHs=False):
                    # It returns None if the molecule parsing generated errors.
                    self._mol_obj = rdk_mol
                    break
        except Exception as e:
            logger.exception(e)
            raise MoleculeObjectError("An error occurred while parsing the molecular file with %s and the molecule "
                                      "object could not be created. Check the logs for more information." % tool)

        if self._mol_obj is None:
            raise MoleculeNotFoundError("Ligand '%s' not found in the input file or generated errors while parsing it with %s."
                                        % (self.mol_id, tool))

        logger.info("Molecule '%s' was successfully loaded." % self.mol_id)

    def get_biopython_structure(self, entity=None, parser=PDB_PARSER):
        if self.mol_obj_type == "openbabel":
            pdb_block = self.mol_obj.write('pdb')
            chain_id = "A"
            comp_name = " "
        else:
            pdb_block = MolToPDBBlock(self.mol_obj)
            chain_id = " "
            comp_name = "H_UNL"

        comp_structure = PDB_PARSER.get_structure_from_pdb_block(self.pdb_id, pdb_block)

        chain = comp_structure[0][chain_id]
        if self.chain_id != chain.id:
            chain.id = self.chain_id

        lig = chain[(comp_name, 1, " ")]
        lig.id = ("H_%s" % self.comp_name, self.comp_num, " ")
        lig.resname = self.comp_name

        if entity is not None:
            if isinstance(entity, Entity):
                structure = entity.get_parent_by_level('S')
                if self.chain_id not in structure[0].child_dict:
                    chain = chain.copy()
                    structure[0].add(chain)
                else:
                    lig = lig.copy()
                    structure[0][self.chain_id].add(lig)
            else:
                raise IllegalArgumentError("The informed entity is not a valid Biopython object.")
        else:
            entity = comp_structure

        return entity

    def __repr__(self):
        return '<MolEntry: %s%s%s>' % (self.pdb_id, self.sep, self.mol_id)


def recover_entries_from_entity(entity, get_ligands=True, get_chains=True, sep=ENTRY_SEPARATOR):

    if entity.level == "S":
        if get_ligands:
            residues = entity[0].get_residues()
        if get_chains:
            chains = entity[0].get_chains()

    elif entity.level == "M":
        if get_ligands:
            residues = entity.get_residues()
        if get_chains:
            chains = entity.get_chains()
    else:
        if get_ligands:
            # If the entity is already a Chain, get_parent_by_level() returns the same object.
            # But, if the entity is a Residue or Atom, it will return its corresponding chain parent.
            residues = entity.get_parent_by_level("C").get_residues()
        if get_chains:
            chains = entity.get_parent_by_level("M").get_chains()

    if get_ligands:
        pdb_id = entity.get_parent_by_level("S").id
        for res in residues:
            if res.is_hetatm():
                entry = sep.join([pdb_id, res.parent.id, res.resname, "%d%s" % res.id[1:]])
                yield entry

    if get_chains:
        pdb_id = entity.get_parent_by_level("S").id
        for chain in chains:
            entry = sep.join([pdb_id, chain.id])
            yield entry
