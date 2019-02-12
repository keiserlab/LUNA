from pybel import (readfile, informats)
from util.file import get_file_format
from util.exceptions import (IllegalArgumentError, MoleculeObjectError, MoleculeNotFoundError)
from os.path import exists
from MyBio.PDB.PDBParser import PDBParser
from MyBio.PDB.Entity import Entity

import re
import logging


logger = logging.getLogger(__name__)

# Source: https://richjenks.com/filename-regex/
FILENAME_REGEX = r"(?!.{256,})(?!(aux|clock\$|con|nul|prn|com[1-9]|lpt[1-9])(?:$|\.))[^ ][ \.\w\-$()+=[\];#@~,&']+[^\. ]"

PDB_PARSER = PDBParser(PERMISSIVE=True, QUIET=True, FIX_ATOM_NAME_CONFLICT=False, FIX_OBABEL_FLAGS=False)


class Entry:

    def __init__(self, pdb_id, chain_id, lig_name=None, lig_num=None, lig_icode=None, separator=':'):
        self.pdb_id = pdb_id
        self.chain_id = chain_id
        self.lig_name = lig_name
        self.lig_num = lig_num
        self._lig_icode = lig_icode
        self.separator = separator

    @property
    def lig_icode(self):
        if isinstance(self._lig_icode, str):
            return self._lig_icode
        else:
            return ' '

    @property
    def full_id(self):
        return (self.pdb_id, self.chain_id, self.lig_name, self.lig_num, self.lig_icode)

    def to_string(self, sep=':'):
        entry = [self.pdb_id, self.chain_id]
        if self.lig_name:
            entry.append(str(self.lig_name))
        if self.lig_num:
            lig_num = str(self.lig_num)
            if self._lig_icode:
                lig_num += self._lig_icode
            entry.append(lig_num)
        return sep.join(entry)

    def is_hetatm(self):
        return self.lig_name and self.lig_num

    def get_biopython_key(self):
        if self.is_hetatm():
            if self.lig_name == 'HOH' or self.lig_name == 'WAT':
                key = ('W', self.lig_num, self.lig_icode)
            else:
                key = ('H_%s' % self.lig_name, self.lig_num, self.lig_icode)
        else:
            key = self.chain_id
        return key

    def __repr__(self):
        return '<Entry: %s>' % self.to_string(self.separator)


class DBEntry(Entry):

    def __init__(self, ligand_entry_id, pdb_id, chain_id, lig_name=None, lig_num=None, lig_icode=None):
        self.id = ligand_entry_id
        super().__init__(pdb_id, chain_id, lig_name, lig_num, lig_icode)


class MolEntry(Entry):

    def __init__(self, pdb_id, mol_id, mol_file, mol_file_type=None):
        self.mol_id = mol_id
        self.mol_file = mol_file
        self.mol_file_type = mol_file_type
        self._mol_obj = None

        super().__init__(pdb_id, "z", "LIG", 9999)

    @property
    def full_id(self):
        return (self.pdb_id, self.mol_id)

    @property
    def mol_obj(self):
        if self._mol_obj is None:
            ext = self.mol_file_type or get_file_format(self.mol_file)

            if ext not in informats:
                raise IllegalArgumentError("Extension '%s' informed or assumed by the filename is not a "
                                           "recognised Open Babel format." % ext)

            if not exists(self.mol_file):
                raise FileNotFoundError("The file '%s' was not found." % self.mol_file)

            try:
                for ob_mol in readfile(ext, self.mol_file):
                    if self.mol_id == ob_mol.OBMol.GetTitle():
                        self._mol_obj = ob_mol
                        break
            except Exception as e:
                logger.exception(e)
                raise MoleculeObjectError("An error occurred while parsing the file with Open Babel and the molecule "
                                          "object could not be created. Check the logs for more information.")

            if self._mol_obj is None:
                raise MoleculeNotFoundError("Ligand '%s' not found in the input file." % self.mol_id)
        return self._mol_obj

    @mol_obj.setter
    def mol_obj(self, obj):
        self._mol_obj = obj

    def get_biopython_structure(self, entity=None, parser=PDB_PARSER):
        pdb_block = self.mol_obj.write('pdb')
        lig_structure = PDB_PARSER.get_structure_from_pdb_block(self.pdb_id, pdb_block)

        chain = lig_structure[0]['A']
        if self.chain_id != chain.id:
            chain.id = self.chain_id

        lig = chain[(" ", 1, " ")]
        lig.id = ("H_%s" % self.lig_name, self.lig_num, " ")
        lig.resname = self.lig_name

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
            entity = lig_structure

        return entity

    def __repr__(self):
        return '<MolEntry: %s%s%s>' % (self.pdb_id, self.separator, self.mol_id)


class EntryValidator:

    def __init__(self, pattern):
        self.pattern = pattern
        self.regex = re.compile(pattern, flags=0)

    def is_valid(self, entry):
        return self.regex.match(entry) is not None


class PLIEntryValidator(EntryValidator):

    def __init__(self, free_filename_format=False):
        if free_filename_format:
            pattern = r'^%s:\w:\w{1,3}:\-?\d{1,4}[a-zA-z]?$' % FILENAME_REGEX
        else:
            pattern = r'^\w{4}:\w:\w{1,3}:\-?\d{1,4}[a-zA-z]?$'

        super().__init__(pattern)


class PPIEntryValidator(EntryValidator):

    def __init__(self, free_filename_format=False):
        if free_filename_format:
            pattern = '^\\w+(\\w|\\-)*:\\w$'
        else:
            pattern = '^\\w{4}:\\w$'
        super().__init__(pattern)


