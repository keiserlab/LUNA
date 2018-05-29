import re


class LigandEntry():

    def __init__(self, pdb_id, chain_id,
                 lig_name=None, lig_num=None,
                 lig_icode=None, separator=":"):

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
            return " "

    def to_string(self, sep=":"):
        entry = [self.pdb_id, self.chain_id]

        if self.lig_name:
            entry.append(str(self.lig_name))
        if self.lig_num:
            lig_num = str(self.lig_num)
            if self._lig_icode:
                lig_num += self._lig_icode

            entry.append(lig_num)

        return sep.join(entry)

    def __repr__(self):
        return '<LigandEntry: %s>' % self.to_string(self.separator)


class DBLigandEntry(LigandEntry):

    def __init__(self, ligand_entry_id, pdb, chain,
                 lig_name=None, lig_num=None, lig_icode=None):

        self.id = ligand_entry_id

        super().__init__(pdb, chain, lig_name, lig_num, lig_icode)


class EntryValidator:

    def __init__(self, pattern):
        self.regex = re.compile(pattern, flags=0)

    def is_valid(self, entry):
        return self.regex.match(entry) is not None


class PLIEntryValidator(EntryValidator):

    def __init__(self, uploadFormat=False):
        pattern = '^\w{4}:\w:\w{1,3}:\-?\d+$'

        if (uploadFormat):
            pattern = '^\w+(\w|\-)*:\w:\w{1,3}:\-?\d+$'

        super(PLIEntryValidator, self).__init__(pattern)


class PPIEntryValidator(EntryValidator):

    def __init__(self, uploadFormat=False):
        pattern = '^\w{4}:\w$'

        if (uploadFormat):
            pattern = '^\w+(\w|\-)*:\w$'

        super(PPIEntryValidator, self).__init__(pattern)
