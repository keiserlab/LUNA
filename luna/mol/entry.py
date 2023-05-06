import re
import logging
from operator import xor
from os.path import exists
from collections import defaultdict
import ast

from rdkit.Chem import Mol as RDMol
from openbabel import OBMol
from openbabel.pybel import readfile
from openbabel.pybel import Molecule as PybelMol
from openbabel.pybel import informats as OB_FORMATS

from luna.wrappers.rdkit import (RDKIT_FORMATS, read_multimol_file,
                                 read_mol_from_file)
from luna.wrappers.base import MolWrapper
from luna.util.default_values import (ACCEPTED_MOL_OBJ_TYPES, ENTRY_SEPARATOR,
                                      ARTIFACTS_LIST)
from luna.util.file import get_file_format, get_filename
from luna.util.exceptions import (InvalidEntry,
                                  IllegalArgumentError,
                                  MoleculeObjectError,
                                  MoleculeObjectTypeError,
                                  MoleculeNotFoundError)
from luna.MyBio.PDB.PDBParser import PDBParser, WATER_NAMES, DEFAULT_CHAIN_ID
from luna.MyBio.PDB.Entity import Entity


logger = logging.getLogger()

PCI_ENTRY_REGEX = re.compile(r'^.{1,255}:\w:\w[\w+\-]{0,2}:'
                             r'\-?\d{1,4}[a-zA-z]?$')
PPI_ENTRY_REGEX = re.compile(r'^.{1,255}:\w$')

REGEX_RESNUM_ICODE = re.compile(r'^([\-\+]?\d+)([a-zA-z]?)$')


class Entry:
    """Entries determine the target molecule to which interactions and other
    properties will be calculated. They can be ligands, chains, etc, and can
    be defined in a number of ways. Each entry has an associated PDB file that
    may contain macromolecules (protein, RNA, DNA) and other small molecules,
    water, and ions. The PDB file provides the context to where the
    interactions with the target molecule will be calculated.

    Parameters
    ----------
    pdb_id : str
        A 4-symbols structure id from PDB or a local PDB filename.
        Example: '3QL8' or 'file1'.
    chain_id : str
        A 1-symbol chain id. Example: 'A'.
    comp_name : str, optional
        A 1 to 3-symbols molecule name (residue name in the PDB format).
        Obligatory if ``is_hetatm`` is True. Example: 'X01'.
    comp_num : int, optional
        A valid 4-digits integer (residue sequence number in the PDB format).
        Obligatory if ``is_hetatm`` is True. Example: 300 or -1.
    comp_icode : str, optional
        A 1-character molecule insertion code (residue insertion code in the
        PDB format). Example: 'A'.
    is_hetatm : bool
        If the molecule is a ligand or not. The default value is True.
    sep : str
        A separator character to format the entry string.
        The default value is ':'.
    parser : :class:`~luna.MyBio.PDB.PDBParser.PDBParser` or \
                :class:`~luna.MyBio.PDB.FTMapParser.FTMapParser`, optional
            Define a PDB parser object. If not provided, the default parser
            will be used.

    Raises
    ------
    IllegalArgumentError
        If ``is_hetatm`` is True, but the molecule name and number are
        not provided. If the molecule number is provided but it is not
        an integer. If ``comp_icode`` is provided but it is not a
        valid character.
    InvalidEntry
        If the provided information does not match the PDB format.

    Examples
    --------

    Chain entry: can be used to calculate interactions with a given chain.

    >>> from luna.mol.entry import Entry
    >>> e = Entry(pdb_id="3QL8", chain_id="A")
    >>> print(e)
    <Entry: 3QL8:A>

    Molecule entry: can be used to calculate interactions with a given molecule
    (residue or nucleotide).

    >>> from luna.mol.entry import Entry
    >>> e = Entry(pdb_id="3QL8", chain_id="A", comp_name="HIS", comp_num=125, \
is_hetatm=False)
    >>> print(e)
    <Entry: 3QL8:A:HIS:125>

    Ligand entry: can be used to calculate interactions with a given ligand.

    >>> from luna.mol.entry import Entry
    >>> e = Entry(pdb_id="3QL8", chain_id="A", comp_name="X01", comp_num=300, \
is_hetatm=True)
    >>> print(e)
    <Entry: 3QL8:A:X01:300>

    You can use a different character separator for the entries. For example:

    >>> from luna.mol.entry import Entry
    >>> e = Entry(pdb_id="3QL8", chain_id="A", comp_name="X01", comp_num=300, \
is_hetatm=True, sep="/")
    >>> print(e)
    <Entry: 3QL8/A/X01/300>
    """

    def __init__(self, pdb_id, chain_id, comp_name=None,
                 comp_num=None, comp_icode=None,
                 is_hetatm=True, sep=ENTRY_SEPARATOR, parser=None):

        if xor(comp_name is None, comp_num is None):
            raise IllegalArgumentError("You tried to define a molecule, so "
                                       "you must inform its name and number.")

        if comp_num is not None:
            try:
                assert float(comp_num).is_integer()
                comp_num = int(comp_num)
            except (ValueError, AssertionError):
                error_msg = ("The informed molecule number '%s' is invalid. "
                             "It must be an integer." % str(comp_num))
                raise IllegalArgumentError(error_msg)

        if comp_icode is not None:
            comp_icode = str(comp_icode)
            if comp_icode.isdigit() or len(comp_icode) > 1:
                error_msg = ("The informed molecule icode '%s' is invalid. "
                             "It must be a character." % str(comp_icode))
                raise IllegalArgumentError(error_msg)

        self._pdb_id = pdb_id
        self._chain_id = chain_id
        self._comp_name = comp_name
        self._comp_num = comp_num
        self._comp_icode = comp_icode
        self.is_hetatm = is_hetatm
        self.sep = sep
        self.parser = parser

        if not self.is_valid():
            raise InvalidEntry("Entry '%s' does not match the PDB format."
                               % self.to_string())

    @classmethod
    def from_string(cls, entry_str, is_hetatm=True, sep=ENTRY_SEPARATOR):
        """Initialize from a string.

        Parameters
        ----------
        entry_str : str
            A string representing the entry. Example: '3QL8:A:X01:300'.
        is_hetatm : bool
            Defines if the molecule is a ligand or not.
            The default value is True.
        sep : str
            The separator character used in ``entry_str``.
            The default value is ':'. For example: if ``entry_str`` is set to
            '3QL8|A|X01|300', then ``sep`` should be defined as '|'.

        Returns
        -------
         : `Entry`

        Raises
        ------
        IllegalArgumentError
            If the fields in ``entry_str`` do not match the format expected to
            define a chain (`ChainEntry`) or a molecule (`MolEntry`).

        Examples
        --------

        Chain entry: can be used to calculate interactions with a given chain.

        >>> from luna.mol.entry import Entry
        >>> e = Entry.from_string("3QL8:A", sep=":")
        >>> print(e)
        <Entry: 3QL8:A>

        Molecule entry: can be used to calculate interactions with a given
        molecule (residue or nucleotide).

        >>> from luna.mol.entry import Entry
        >>> e = Entry.from_string("3QL8:A:HIS:125", sep=":")
        >>> print(e)
        <Entry: 3QL8:A:HIS:125>

        Ligand entry: can be used to calculate interactions with a given
        ligand.

        >>> from luna.mol.entry import Entry
        >>> e = Entry.from_string("3QL8:A:X01:300", sep=":")
        >>> print(e)
        <Entry: 3QL8:A:X01:300>
        """

        entries = entry_str.split(sep)

        # Try to initialize a new ChainEntry.
        if len(entries) == 2:
            if any([str(i).strip() == "" for i in entries]):
                error_msg = ("The number of fields in the informed string "
                             "'%s' is incorrect. A valid ChainEntry must "
                             "contain two obligatory fields: PDB and chain "
                             "id." % entry_str)
                raise IllegalArgumentError(error_msg)

            return cls(*entries, is_hetatm=False, sep=sep)

        # Try to initialize a new MolEntry.
        elif len(entries) == 4:
            if any([str(i).strip() == "" for i in entries]):
                error_msg = ("The number of fields in the informed string "
                             "'%s' is incorrect. A valid MolEntry must "
                             "contain four obligatory fields: PDB, chain id, "
                             "molecule name, and molecule number followed by "
                             "its insertion code when applicable." % entry_str)
                raise IllegalArgumentError(error_msg)

            # Separate ligand number from insertion code.
            matched = REGEX_RESNUM_ICODE.match(entries[3])
            if matched:
                comp_num = matched.group(1)
                try:
                    assert float(comp_num).is_integer()
                    comp_num = int(comp_num)
                except (ValueError, AssertionError):
                    error_msg = ("The informed molecule number '%s' is "
                                 "invalid. It must be an integer."
                                 % str(comp_num))
                    raise IllegalArgumentError(error_msg)

                icode = None if matched.group(2) == "" else matched.group(2)
                entries = entries[0:3] + [comp_num, icode]
            else:
                error_msg = ("The molecule number and its insertion code "
                             "(if applicable) '%s' is invalid. It must be an "
                             "integer followed by one insertion code "
                             "character when applicable." % entries[3])
                raise IllegalArgumentError(error_msg)

            return cls(*entries, is_hetatm=is_hetatm, sep=sep)
        else:
            error_msg = ("The number of fields in the informed string '%s' is "
                         "incorrect. A valid string must contain two "
                         "obligatory fields (PDB and chain id) and it may "
                         "contain two optional fields (molecule name and "
                         "molecule number followed by its insertion code when "
                         "applicable)." % entry_str)
            raise IllegalArgumentError(error_msg)

    @classmethod
    def from_file(cls, input_file, pdb_id=None, mol_file=None,
                  entries_sep=":", fields_sep=",", **kwargs):
        """Initialize from a file containing a list of strings
        representing entries.

        Parameters
        ----------
        input_file : str
            The file from where the list of strings (one per line) will
            be read from.

            .. note::
                The list of strings can be a mix of complexes whose
                molecules lie all in the same PDB file, or whose protein
                structure is in a PDB file and whose ligand is on a
                molecular file. See ``entries_sep`` and ``fields_sep``.
        pdb_id : str
            A 4-symbols structure id from PDB or a local PDB filename.
            Example: '3QL8' or 'file1'.

            .. note::
                Mandatory if there is any line containing only the ligand name.
                This only applies to ligands that will be read from
                multimolecular files.
        mol_file : str
            Pathname of the molecular file.

            .. note::
                Mandatory if there is any line containing only the ligand name.
                This only applies to ligands that will be read from
                multimolecular files.
        entries_sep : str
            a separator for entries whose molecules lie all in the same PDB
            file, not requiring, therefore, an additional molecular file.
            Such entries cause the function to yield `ChainEntry` or
            `MolEntry` objects. Default: ':'. E.g.: 3QQK:A:X02:497.
        fields_sep : str
            a character used to separate fields in ``input_file`` for
            complexes whose protein structure is in a PDB file and
            whose ligand is on a molecular file. Such entries cause the
            function to yield `MolFileEntry` objects. Default: ','

            Each line should contain either the ligand name only or the
            following fields:

                * ``pdb_id`` (str): a 4-symbols structure id from \
PDB or a local PDB filename. Example: '3QL8' or 'file1'.
                * ``mol_id`` (str): the ligand id in the molecular file.

                * ``mol_file`` (str): pathname of the molecular file.

                * ``is_multimol_file`` (bool): if ``mol_file`` contains \
multiple molecules or not. If True, ``mol_id`` should match some ligand name \
in ``mol_file``.

            .. note::
                If a line contains only the ligand name, providing a value to
                the arguments ``pdb_id`` and ``mol_file`` is mandatory.
        kwargs : dict, optional
            Extra arguments to lines that yield `MolFileEntry` objects.
            Refer to the documentation for a list of all possible arguments.

        Yields
        ------
        `ChainEntry`, `MolEntry`, `MolFileEntry`
            Define a chain, a molecule from a PDB file, or a molecule from
            a molecular file, respectively, according to a string
            from ``input_file``.

        Raises
        ------
        IllegalArgumentError
            If a given line contain a string that do not match the format
            expected to define a chain (`ChainEntry`) or a molecule from a PDB
            file (`MolEntry`), or a molecule from a molecular file
            (`MolFileEntry`).

        Examples
        --------

        Example of a file containing a mix of entry types:

        >>> 3QQK:A:X02:497
        >>> ZINC000065293174
        >>> protein,ZINC000012442563,inputs/ligands.mol2,True
        >>> protein,ZINC000007786517,inputs/ZINC000007786517.mol,False
        """
        if not exists(input_file):
            raise OSError("File '%s' does not exist." % input_file)

        with open(input_file, "r") as IN:
            c = 1
            for row in IN:
                row = row.strip()
                if row == "":
                    continue

                args = row.split(fields_sep)

                has_error = False
                if len(args) == 1:
                    args = row.split(entries_sep)

                    if len(args) == 1:
                        if pdb_id is None or mol_file is None:
                            msg = ("It seems a molecule name was provided "
                                   f"in line #{c}. In this case, "
                                   "'pdb_id' and 'mol_file' are mandatory.")
                            raise IllegalArgumentError(msg)

                        yield MolFileEntry.from_mol_file(pdb_id, args[0],
                                                         mol_file,
                                                         is_multimol_file=True,
                                                         **kwargs)

                    elif len(args) == 2:
                        yield ChainEntry.from_string(row, sep=entries_sep)

                    elif len(args) == 4:
                        yield MolEntry.from_string(row, sep=entries_sep)

                    else:
                        has_error = True

                elif len(args) == 4:
                    (curr_pdb_id, curr_mol_id,
                        curr_mol_file, is_multimol) = args
                    is_multimol = ast.literal_eval(is_multimol)

                    cls = MolFileEntry
                    yield cls.from_mol_file(curr_pdb_id, curr_mol_id,
                                            curr_mol_file,
                                            is_multimol_file=is_multimol,
                                            **kwargs)
                else:
                    has_error = True

                if has_error:
                    msg = (f"Invalid number of arguments found in line #{c}. "
                           "You should provide either one or four arguments. "
                           "If only one argument is provided, it can be an "
                           "entry string (e.g., 3QQK:A:X02:497) or a "
                           "molecule name.")
                    raise IllegalArgumentError(msg)

                c += 1

    @property
    def pdb_id(self):
        """str, read-only: the pdb id."""
        return self._pdb_id

    @property
    def chain_id(self):
        """str, read-only: the chain id."""
        return self._chain_id

    @property
    def comp_name(self):
        """str, read-only: the molecule name."""
        return self._comp_name

    @property
    def comp_num(self):
        """int, read-only: the molecule number."""
        return self._comp_num

    @property
    def comp_icode(self):
        """str, read-only: the molecule insertion code."""
        if isinstance(self._comp_icode, str):
            return self._comp_icode
        else:
            return ' '

    @property
    def full_id(self):
        """tuple, read-only: The full id of the entry is the tuple
        (PDB id or filename, chain id) for entries representing chains and
        (PDB id or filename, chain id, molecule name, molecule number,
        insertion code) for entries representing molecules."""

        entry = [self.pdb_id, self.chain_id]
        if self.comp_name is not None and self.comp_num is not None:
            entry.append(self.comp_name)
            entry.append(self.comp_num)
            entry.append(self.comp_icode)
        return tuple(entry)

    def to_string(self, sep=None):
        """ Convert the entry to a string using ``sep`` as a
        separator character.

        Parameters
        ----------
        sep : str or None
            If None (the default), use the separator character defined during
            the entry object creation. Otherwise, uses ``sep`` as the
            separator character.

        Examples
        --------

        >>> from luna.mol.entry import Entry
        >>> e = Entry(pdb_id="3QL8", chain_id="A", comp_name="X01", \
comp_num=300, is_hetatm=True, sep=":")
        >>> print(e.to_string("/"))
        3QL8/A/X01/300
        """

        full_id = self.full_id

        # An entry object will always have a PDB and chain id.
        entry = list(full_id[0:2])

        # If it contains additional information about the molecule it
        # will also include them.
        if len(full_id) > 2:
            if full_id[2] is not None and full_id[3] is not None:
                comp_name = str(full_id[2]).strip()
                comp_num_and_icode = \
                    str(full_id[3]).strip() + str(full_id[4]).strip()
                entry += [comp_name, comp_num_and_icode]

        sep = sep or self.sep

        return sep.join(entry)

    def is_valid(self):
        """Check if the entry matches the expected format for protein-protein
        or protein-molecule complexes.

        Returns
        -------
         : bool
        """
        full_id = self.full_id

        # Regex for ChainEntry (pdb_id, chain_id).
        if len(full_id) == 2:
            return PPI_ENTRY_REGEX.match(self.to_string(":")) is not None

        # Regex for MolEntry (pdb_id, chain_id, comp_name, comp_num, icode).
        elif len(full_id) == 5:
            return PCI_ENTRY_REGEX.match(self.to_string(":")) is not None

        # Return False for anything else
        return False

    def get_biopython_key(self, full_id=False):
        """Represent the entry as a key to select chains or molecules from
        Biopython Entity objects.

        Parameters
        full_id : bool
            If True, return the full id of a chain or ligand.
            For chains, it consists of a tuple containing the PDB and the
            chain id. For ligands, it consists of a tuple containing the PDB,
            the chain, and the ligand id. The default value is False.

        Returns
        -------
         : str or tuple
            Return str if the entry represents a chain and if ``full_id``
            is False. Otherwise, return a tuple.

        Examples
        --------

        >>> from luna.mol.entry import Entry
        >>> e = Entry(pdb_id="3QL8", chain_id="A", comp_name="X01", \
comp_num=300, is_hetatm=True, sep=":")
        >>> print(e.get_biopython_key())
        ('H_X01', 300, ' ')

        """
        key = []
        if full_id:
            key = [self.pdb_id, 0, self.chain_id]

        if self.comp_name is not None and self.comp_num is not None:
            if self.comp_name == 'HOH' or self.comp_name == 'WAT':
                comp_id = ('W', self.comp_num, self.comp_icode)
            elif self.is_hetatm:
                comp_id = ('H_%s' % self.comp_name, self.comp_num,
                           self.comp_icode)
            else:
                comp_id = (' ', self.comp_num, self.comp_icode)

            if full_id:
                key.append(comp_id)
                return tuple(key)
            return comp_id

        if full_id:
            return tuple(key)

        return self.chain_id

    def __repr__(self):
        return '<%s: %s>' % (self.__class__.__name__, self.to_string(self.sep))


class ChainEntry(Entry):
    """Define a chain.

    Parameters
    ----------
    pdb_id : str
        A 4-symbols structure id from PDB or a local PDB filename.
        Example: '3QL8' or 'file1'.
    chain_id : str
        A 1-symbol chain id. Example: 'A'.
    sep : str
        A separator character to format the entry string.
        The default value is ':'.

    Raises
    ------
    InvalidEntry
        If the provided information does not match the PDB format.

    Examples
    --------
    >>> from luna.mol.entry import ChainEntry
    >>> e = ChainEntry(pdb_id="3QL8", chain_id="A")
    >>> print(e)
    <ChainEntry: 3QL8:A>

    """

    def __init__(self, pdb_id, chain_id, sep=ENTRY_SEPARATOR, parser=None):
        super().__init__(pdb_id, chain_id, is_hetatm=False, sep=sep,
                         parser=parser)

    @classmethod
    def from_string(cls, entry_str, sep=ENTRY_SEPARATOR):
        """Initialize from a string.

        Parameters
        ----------
        entry_str : str
            A string representing the entry. Example: '3QL8:A'.
        sep : str
            The separator character used in ``entry_str``.
            The default value is ':'. For example: if ``entry_str`` is set to
            '3QL8|A', then ``sep`` should be defined as '|'.

        Returns
        -------
         : `Entry`

        Raises
        ------
        IllegalArgumentError
            If the fields in ``entry_str`` do not match the format expected
            to define a chain.

        Examples
        --------
        >>> from luna.mol.entry import ChainEntry
        >>> e = ChainEntry.from_string("3QL8:A", sep=":")
        >>> print(e)
        <ChainEntry: 3QL8:A>

        """
        entries = entry_str.split(sep)
        if len(entries) == 2:
            return cls(*entries, sep=sep)
        else:
            error_msg = ("The number of fields in the informed string '%s' is "
                         "incorrect. A valid string must contain two "
                         "obligatory fields: PDB and chain id." % entry_str)
            raise IllegalArgumentError(error_msg)

    @property
    def full_id(self):
        """tuple, read-only: The full id of the entry is the tuple
        (PDB id or filename, chain id)."""
        return (self.pdb_id, self.chain_id)


class MolEntry(Entry):
    """Define a molecule from a PDB file, which can be a residue,
    nucleotide, ligand, water, or ion.

    Parameters
    ----------
    pdb_id : str
        A 4-symbols structure id from PDB or a local PDB filename.
        Example: '3QL8' or 'file1'.
    chain_id : str
        A 1-symbol chain id. Example: 'A'.
    comp_name : str
        A 1 to 3-symbols molecule name (residue name in the PDB format).
        Example: 'X01'.
    comp_num : int
        A valid 4-digits integer (residue sequence number in the PDB format).
        Example: 300 or -1.
    comp_icode : str, optional
        A 1-character molecule insertion code (residue insertion code in the
        PDB format). Example: 'A'.
    sep : str
        A separator character to format the entry string.
        The default value is ':'.

    Raises
    ------
    InvalidEntry
        If the provided information does not match the PDB format.

    Examples
    --------

    Molecule entry: can be used to calculate interactions with a given
    molecule (residue or nucleotide).

    >>> from luna.mol.entry import MolEntry
    >>> e = MolEntry(pdb_id="3QL8", chain_id="A", comp_name="HIS", \
comp_num=125, is_hetatm=False)
    >>> print(e)
    <MolEntry: 3QL8:A:HIS:125>

    Ligand entry: can be used to calculate interactions with a given ligand.

    >>> from luna.mol.entry import MolEntry
    >>> e = MolEntry(pdb_id="3QL8", chain_id="A", comp_name="X01", \
comp_num=300, is_hetatm=True)
    >>> print(e)
    <MolEntry: 3QL8:A:X01:300>

    """

    def __init__(self, pdb_id, chain_id, comp_name, comp_num,
                 comp_icode=None, sep=ENTRY_SEPARATOR, parser=None):

        super().__init__(pdb_id, chain_id, comp_name, comp_num, comp_icode,
                         is_hetatm=True, sep=sep, parser=parser)

    @classmethod
    def from_string(cls, entry_str, sep=ENTRY_SEPARATOR):
        entries = entry_str.split(sep)
        if len(entries) == 4:
            # Separate ligand number from insertion code.
            matched = REGEX_RESNUM_ICODE.match(entries[3])
            if matched:
                comp_num = matched.group(1)
                icode = None if matched.group(2) == "" else matched.group(2)
                entries = entries[0:3] + [comp_num, icode]
            else:
                error_msg = ("The molecule number and its insertion code "
                             "(if applicable) '%s' is invalid. It must be an "
                             "integer followed by one insertion code "
                             "character when applicable." % entries[3])
                raise IllegalArgumentError(error_msg)
            return cls(*entries, sep=sep)
        else:
            error_msg = ("The number of fields in the informed string '%s' is "
                         "incorrect. A valid molecule entry must contain four "
                         "obligatory fields: PDB, chain id, molecule name, "
                         "and molecule number followed by its insertion code "
                         "when applicable." % entry_str)
            raise IllegalArgumentError(error_msg)


class MolFileEntry(Entry):

    """Define a ligand from a molecular file.
    This class should be used for docking and molecular dynamics campaigns
    where usually one has the protein structure in the PDB format and the
    ligand structure in a separate molecular file.

    Parameters
    ----------
    pdb_id : str
        A 4-symbols structure id from PDB or a local PDB filename.
        Example: '3QL8' or 'file1'.
    mol_id : str
        The ligand id in the molecular file.
    sep : str
        A separator character to format the entry string.
        The default value is ':'.

    Attributes
    ----------
    mol_id : str
        The ligand id.
    mol_file : str
        Pathname of the molecular file.
    mol_file_ext : str
        The molecular file format.
        If not provided, try to recover the molecular file extension
        directly from ``mol_file``.
    mol_obj_type : {'rdkit', 'openbabel'}
        Define which library (RDKit or Open Babel) to use to parse
        the molecular file.
    overwrite_mol_name : bool
        If True, substitute the ligand name in the parsed molecular object
        with ``mol_id``. Only works for single-molecule files
        (``is_multimol_file`` = False) as in these cases ``mol_id`` does not
        need to match the ligand name in the molecular file.
    is_multimol_file : bool
        If ``mol_file`` contains multiple molecules or not.
        If True, ``mol_id`` should match some ligand name in ``mol_file``.
    """

    def __init__(self, pdb_id, mol_id, sep=ENTRY_SEPARATOR):

        self.mol_id = mol_id

        #
        # Initialize empty properties.
        #
        self._mol_obj = None
        self.mol_file = None
        # TODO: Find a way to assume the Mol file type when not provided
        self.mol_file_ext = None
        self.mol_obj_type = None
        self.overwrite_mol_name = None
        self.is_multimol_file = None

        super().__init__(pdb_id, DEFAULT_CHAIN_ID, "LIG", 9999,
                         is_hetatm=True, sep=sep)

    @classmethod
    def from_mol_obj(cls, pdb_id, mol_id, mol_obj, sep=ENTRY_SEPARATOR):
        """Initialize from an already loaded molecular object.

        This function is useful in cases where a molecular object is parsed and
        pre-processed using a different protocol defined by the user.

        Parameters
        ----------
        pdb_id : str
            A 4-symbols structure id from PDB or a local PDB filename.
            Example: '3QL8' or 'file1'.
        mol_id : str
            The ligand id.
            As the molecular object is already provided, the ligand id does
            not need to match the ligand name in the molecular object.
        mol_obj : :class:`~luna.wrappers.base.MolWrapper`, \
                    :class:`rdkit.Chem.rdchem.Mol`, or \
                    :class:`openbabel.pybel.Molecule`
            The molecular object.
        sep : str
            A separator character to format the entry string.
            The default value is ':'.

        Returns
        -------
         : `MolFileEntry`

        Raises
        ------
        MoleculeObjectTypeError
            If the molecular object is not an instance
            of :class:`~luna.wrappers.base.MolWrapper`,
            :class:`rdkit.Chem.rdchem.Mol`,
            or :class:`openbabel.pybel.Molecule`.
        IllegalArgumentError
            If ``entity`` is not a valid Biopython object.

        Examples
        --------

        In this example, we will initialize a new MolFileEntry with the ligand
        'ZINC000007786517' and the structure located in a PDB file of name
        'D4', which is the structure used for docking the molecule.

        First, let's parse the molecular file.

        >>> from luna.wrappers.rdkit import read_mol_from_file
        >>> mol_file = "tutorial/inputs/ZINC000007786517.mol"
        >>> mol_obj = read_mol_from_file(mol_file, mol_format="mol")

        Now, we create the new MolFileEntry object as follows:

        >>> e = MolFileEntry.from_mol_obj("D4", "ZINC000007786517",
        ...                               mol_obj, sep=ENTRY_SEPARATOR)
        >>> print(e)
        <MolFileEntry: D4:ZINC000007786517>
        >>> print(e.mol_obj.to_smiles())
        Cc1cccc(NC(=O)C[N@@H+](C)C2CCCCC2)c1C
        """

        if mol_obj is not None:
            if isinstance(mol_obj, MolWrapper):
                mol_obj = mol_obj.unwrap()
            elif isinstance(mol_obj, PybelMol):
                mol_obj = mol_obj.OBMol

            if isinstance(mol_obj, RDMol):
                mol_obj_type = "rdkit"
            elif isinstance(mol_obj, OBMol):
                mol_obj_type = "openbabel"
            else:
                logger.exception("Objects of type '%s' are not currently "
                                 "accepted." % mol_obj.__class__)
                raise MoleculeObjectTypeError("Objects of type '%s' are not "
                                              "currently accepted."
                                              % mol_obj.__class__)
        else:
            if mol_obj_type not in ACCEPTED_MOL_OBJ_TYPES:
                error_msg = ("Objects of type '%s' are not currently "
                             "accepted. The available options are: %s."
                             % (mol_obj_type,
                                ", ".join(ACCEPTED_MOL_OBJ_TYPES)))
                raise IllegalArgumentError(error_msg)

        entry = cls(pdb_id, mol_id, sep)
        entry.mol_obj = mol_obj
        entry.mol_obj_type = mol_obj_type

        # TODO: Find a way to assume the Mol file type when not provided
        # mol_file_ext

        return entry

    @classmethod
    def from_mol_file(cls, pdb_id, mol_id, mol_file, is_multimol_file,
                      mol_file_ext=None, mol_obj_type='rdkit',
                      autoload=False, overwrite_mol_name=False,
                      sep=ENTRY_SEPARATOR):
        """Initialize from a molecular file.

        Parameters
        ----------
        pdb_id : str
            A 4-symbols structure id from PDB or a local PDB filename.
            Example: '3QL8' or 'file1'.
        mol_id : str
            The ligand id in the molecular file.
        mol_file : str
            Pathname of the molecular file.
        is_multimol_file : bool
            If ``mol_file`` contains multiple molecules or not.
            If True, ``mol_id`` should match some ligand name in ``mol_file``.
        mol_file_ext : str, optional
            The molecular file format.
            If not provided, try to recover the molecular file extension
            directly from ``mol_file``.
        mol_obj_type : {'rdkit', 'openbabel'}
            If "rdkit", parse the converted molecule with RDKit and return an
            instance of :class:`rdkit.Chem.rdchem.Mol`. If "openbabel", parse
            the converted molecule with Open Babel and return an instance of
            :class:`openbabel.pybel.Molecule`. The default value is 'rdkit'.
        autoload : bool
            If True, parse the ligand from the molecular file during the entry
            initialization. Otherwise, only load the ligand when first used.
        overwrite_mol_name : bool
            If True, substitute the ligand name in the parsed molecular object
            with ``mol_id``. Only works for single-molecule files
            (``is_multimol_file`` = False) as in these cases ``mol_id`` does
            not need to match the ligand name in the molecular file.
        sep : str
            A separator character to format the entry string.
            The default value is ':'.

        Returns
        -------
         : `MolFileEntry`

        Raises
        ------
        FileNotFoundError
            If ``mol_file`` does not exist.
        IllegalArgumentError
            If ``mol_obj_type`` is not either 'rdkit' nor 'openbabel'.
        MoleculeObjectError
            If any errors occur while parsing the molecular file. Detailed
            information about the errors can be found in the logging outputs.
        MoleculeNotFoundError
            If the ligand ``mol_id`` was not found in the input file and
            ``is_multimol_file`` is True.

        Examples
        --------

        In this first example, we will read the ligand 'ZINC000007786517' from
        a single-molecule file. As we are working with a single-molecule file,
        ``mol_id`` can be any value you prefer.

        >>> from luna.mol.entry import MolFileEntry
        >>> mol_file = "tutorial/inputs/ZINC000007786517.mol"
        >>> e = MolFileEntry.from_mol_file(pdb_id="D4", mol_id="Ligand",
        ...                                mol_file=mol_file, \
mol_obj_type='rdkit',
        ...                                is_multimol_file=False)
        >>> print(e)
        D4:Ligand
        >>> print(e.mol_obj.to_smiles())
        Cc1cccc(NC(=O)C[N@@H+](C)C2CCCCC2)c1C

        Now, let's say we need to read the ligand 'ZINC000096459890' from a
        multi-molecular file and that we want to use Open Babel to parse the
        molecule. To do so, remember that it should exist a ligand with the
        name ``mol_id`` in ``mol_file``. Otherwise, it will raise the
        exception MoleculeNotFoundError.

        >>> from luna.mol.entry import MolFileEntry
        >>> mol_file = "tutorial/inputs/ligands.mol2"
        >>> e = MolFileEntry.from_mol_file(pdb_id="D4", \
mol_id="ZINC000096459890",
        ...                                mol_file=mol_file, \
mol_obj_type='openbabel',
        ...                                is_multimol_file=True)
        >>> print(e)
        <MolFileEntry: D4:ZINC000096459890>
        >>> print(e.mol_obj.to_smiles())
        O=C(OCCCN1C=CC=CC1=O)c1ccc2ccc(Cl)cc2n1

        Below, we show what happens if ``mol_id`` does not exist in
        ``mol_file``. Observe we set ``autoload`` to True to parse the
        molecule right away.

        >>> from luna.mol.entry import MolFileEntry
        >>> mol_file = "tutorial/inputs/ligands.mol2"
        >>> e = MolFileEntry.from_mol_file(pdb_id="D4", mol_id="Ligand",
        ...                                mol_file=mol_file, \
mol_obj_type='openbabel',
        ...                                 is_multimol_file=True, \
autoload=True)
        luna.util.exceptions.MoleculeNotFoundError: "The ligand 'Ligand' was \
not found in the
        input file or generated errors while parsing it with Open Babel."
        """

        entry = cls(pdb_id, mol_id, sep)

        entry.mol_file = mol_file
        entry.is_multimol_file = is_multimol_file
        entry.mol_file_ext = mol_file_ext or get_file_format(mol_file)
        entry.mol_obj_type = mol_obj_type
        entry.overwrite_mol_name = overwrite_mol_name

        if autoload:
            entry._load_mol_from_file()

        return entry

    @classmethod
    def from_file(cls, input_file, pdb_id=None, mol_file=None,
                  fields_sep=",", **kwargs):

        if not exists(input_file):
            raise OSError("File '%s' does not exist." % input_file)

        with open(input_file, "r") as IN:
            c = 1
            for row in IN:
                row = row.strip()
                if row == "":
                    continue

                args = row.split(fields_sep)

                if len(args) == 1:
                    if pdb_id is None or mol_file is None:
                        msg = (f"Only one field was found in line #{c}. "
                               "In these cases, 'pdb_id' and 'mol_file' "
                               "are mandatory.")
                        raise IllegalArgumentError(msg)

                    curr_pdb_id = pdb_id
                    curr_mol_id = args[0]
                    curr_mol_file = mol_file
                    is_multimol_file = True

                elif len(args) == 4:
                    (curr_pdb_id, curr_mol_id,
                        curr_mol_file, is_multimol_file) = args
                    is_multimol_file = ast.literal_eval(is_multimol_file)

                else:
                    msg = (f"Invalid number of arguments found in line #{c}. "
                           "You should provide either one or four arguments.")
                    raise IllegalArgumentError(msg)

                yield cls.from_mol_file(curr_pdb_id, curr_mol_id,
                                        curr_mol_file,
                                        is_multimol_file=is_multimol_file,
                                        **kwargs)
                c += 1

    @property
    def full_id(self):
        """tuple, read-only: The full id of the entry is the tuple
        (PDB id or filename, ligand id)."""
        return (self.pdb_id, self.mol_id)

    @property
    def mol_obj(self):
        """:class:`~luna.wrappers.base.MolWrapper`, \
        :class:`rdkit.Chem.rdchem.Mol`, \
        or :class:`openbabel.pybel.Molecule`: The molecule."""
        if self._mol_obj is None and self.mol_file is not None:
            self._load_mol_from_file()
        return self._mol_obj

    @mol_obj.setter
    def mol_obj(self, mol_obj):
        self._mol_obj = MolWrapper(mol_obj)

    def is_valid(self):
        """Check if the entry represents a valid protein-ligand complex.

        Returns
        -------
         : bool
        """
        return True

    def is_mol_obj_loaded(self):
        """Check if the molecular object has already been loaded.

        Returns
        -------
         : bool
        """
        return self._mol_obj is not None

    def _load_mol_from_file(self):
        logger.debug("It will try to load the molecule '%s'." % self.mol_id)

        if self.mol_file is None:
            raise IllegalArgumentError("It cannot load the molecule as no "
                                       "molecular file was provided.")

        available_formats = (OB_FORMATS if self.mol_obj_type == "openbabel"
                             else RDKIT_FORMATS)
        tool = "Open Babel" if self.mol_obj_type == "openbabel" else "RDKit"
        if self.mol_file_ext not in available_formats:
            raise IllegalArgumentError("Extension '%s' informed or assumed "
                                       "from the filename is not a format "
                                       "recognized by %s."
                                       % (self.mol_file_ext, tool))

        if not exists(self.mol_file):
            raise FileNotFoundError("The file '%s' was not found."
                                    % self.mol_file)

        try:
            if self.mol_obj_type == "openbabel":
                mols = readfile(self.mol_file_ext, self.mol_file)
                # If it is a multimol file, then we need to loop over the
                # molecules to find the target one. Note that in this case,
                # the ids must match.
                if self.is_multimol_file:
                    for ob_mol in mols:
                        filename = get_filename(ob_mol.OBMol.GetTitle())
                        if self.mol_id == filename:
                            self.mol_obj = ob_mol
                            break
                else:
                    self.mol_obj = mols.__next__()
            else:
                if self.mol_file_ext == "pdb":
                    self.mol_obj = \
                        read_mol_from_file(self.mol_file,
                                           mol_format=self.mol_file_ext,
                                           removeHs=False)
                else:
                    # If 'targets' is None, then the entire Mol file
                    # will be read.
                    targets = None
                    # If it is a multimol file than loop through it until the
                    # informed molecule (by its mol_id) is found.
                    if self.is_multimol_file:
                        targets = [self.mol_id]

                    for rdk_mol, mol_id \
                        in read_multimol_file(self.mol_file,
                                              mol_format=self.mol_file_ext,
                                              targets=targets,
                                              removeHs=False):
                        # It returns None if the molecule parsing
                        # generated errors.
                        self.mol_obj = rdk_mol
                        break
        except Exception as e:
            logger.exception(e)

            error_msg = ("An error occurred while parsing the molecular file "
                         "with %s and the molecule object for the entry '%s' "
                         "could not be created. Check the logs for more "
                         "information." % (tool, self.to_string()))
            raise MoleculeObjectError(error_msg)

        if self._mol_obj is None:
            error_msg = ("The ligand '%s' was not found in the input file or "
                         "generated errors while parsing it with %s."
                         % (self.mol_id, tool))
            raise MoleculeNotFoundError(error_msg)
        else:
            if not self.mol_obj.has_name() or self.overwrite_mol_name:
                self.mol_obj.set_name(self.mol_id)

        logger.debug("Molecule '%s' was successfully loaded." % self.mol_id)

    def get_biopython_structure(self, entity=None, parser=None):
        """Transform the molecular object into a Biopython Entity object.

        If ``entity`` is provided, the molecular object is appended to it,
        i.e., this function can be used to join a ligand and the structure
        used during docking or molecular dynamics.

        By default, the ligand is added to a chain of id `z`.

        Parameters
        ----------
        entity : :class:`~luna.MyBio.PDB.Entity.Entity`, optional
            Append the molecular object to ``entity``.
            If not provided, a new :class:`~luna.MyBio.PDB.Entity.Entity`
            is created.
        parser : :class:`~luna.MyBio.PDB.PDBParser.PDBParser`, optional
            Define a PDB parser object. If not provided, the default parser
            will be used.

        Returns
        -------
         : :class:`~luna.MyBio.PDB.Entity.Entity`

        Raises
        ------
        IllegalArgumentError
            If ``entity`` is not a valid Biopython object.

        Examples
        --------

        In this example, we will demonstrate how to join a protein structure
        and a ligand docked against it.

        First, let's parse the PDB file.

        >>> from luna.MyBio.PDB.PDBParser import PDBParser
        >>> pdb_parser = PDBParser(PERMISSIVE=True, QUIET=True)
        >>> pdb_file = "tutorial/inputs/D4.pdb"
        >>> structure = pdb_parser.get_structure("Protein", pdb_file)

        Observe that the list of chains in the parsed structure contains
        only one element.

        >>> print(structure[0].child_list)
        [<Chain id=A>]

        Now, we will read the ligand and append it to the existing protein
        structure.

        >>> from luna.mol.entry import MolFileEntry
        >>> mol_file = "tutorial/inputs/ZINC000007786517.mol"
        >>> e = MolFileEntry.from_mol_file("D4", "ZINC000007786517", mol_file,
        ...                                mol_obj_type='rdkit', \
is_multimol_file=False)
        >>> joined_structure = e.get_biopython_structure(structure)

        Observe that now the list of chains contains chains 'A' and 'z',
        which is the default chain where ligands are added.

        >>> print(joined_structure[0].child_list)
        [<Chain id=A>, <Chain id=z>]

        If we loop over the residues in chain 'z', we will find our ligand.

        >>> for r in structure[0]["z"]:
        >>>     print(r)
        <Residue LIG het=H_LIG resseq=9999 icode= >
        """

        if parser is None:
            parser = PDBParser(PERMISSIVE=True, QUIET=True,
                               FIX_EMPTY_CHAINS=True,
                               FIX_ATOM_NAME_CONFLICT=True,
                               FIX_OBABEL_FLAGS=False)

        mol_file_ext = self.mol_file_ext
        if mol_file_ext is None and self.mol_file is not None:
            mol_file_ext = get_file_format(self.mol_file)

        if self.mol_obj_type == "openbabel":
            pdb_block = self.mol_obj.to_pdb_block()

            atm = self.mol_obj.unwrap().GetFirstAtom()
            residue_info = atm.GetResidue()

            # When the PDBParser finds an empty chain, it automatically
            # replace it by 'z'.
            chain_id = (residue_info.GetChain()
                        if residue_info.GetChain().strip() != ""
                        else self.chain_id)
            comp_num = residue_info.GetNum()

            if residue_info.GetName() in WATER_NAMES:
                comp_name = "W"
            elif residue_info.IsHetAtom(atm):
                comp_name = "H_%s" % residue_info.GetName()
            else:
                comp_name = " "

            if mol_file_ext == "pdb":
                self.chain_id = chain_id
                self.comp_name = residue_info.GetName()
                self.comp_num = comp_num
                self.is_hetatm = residue_info.IsHetAtom(atm)
        else:
            pdb_block = self.mol_obj.to_pdb_block()

            if mol_file_ext == "pdb":
                residue_info = (self.mol_obj.unwrap()
                                .GetAtoms()[0].GetPDBResidueInfo())
                # When the PDBParser finds an empty chain, it automatically
                # replace it by 'z'.
                chain_id = (residue_info.GetChainId()
                            if residue_info.GetChainId().strip() != ""
                            else self.chain_id)
                comp_num = residue_info.GetResidueNumber()

                if residue_info.GetResidueName() in WATER_NAMES:
                    comp_name = "W"
                elif residue_info.GetIsHeteroAtom():
                    comp_name = "H_%s" % residue_info.GetResidueName()
                else:
                    comp_name = " "

                self.chain_id = chain_id
                self.comp_name = residue_info.GetResidueName()
                self.comp_num = comp_num
                self.is_hetatm = residue_info.GetIsHeteroAtom()
            else:
                # When the PDBParser finds an empty chain,
                # it automatically replace it by 'z'.
                chain_id = self.chain_id
                comp_name = "H_UNL"
                comp_num = 1

        comp_structure = parser.get_structure_from_pdb_block(self.pdb_id,
                                                             pdb_block)

        chain = comp_structure[0][chain_id]
        if self.chain_id != chain.id:
            chain.id = self.chain_id

        lig = chain[(comp_name, comp_num, " ")]

        # It only substitutes the ligand id if it is different from the id
        # defined by the MolFileEntry object property. This update will never
        # happen when the ligand file is a PDB file as the ids are guaranteed
        # to be equal.
        if lig.id != ("H_%s" % self.comp_name, self.comp_num, " "):
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
                    # Update the ligand index according to the number of
                    # molecules already present in the chain.
                    lig.idx = len(structure[0][self.chain_id].child_list)
                    structure[0][self.chain_id].add(lig)
            else:
                raise IllegalArgumentError("The informed entity is not a "
                                           "valid Biopython object.")
        else:
            entity = comp_structure

        return entity

    def __repr__(self):
        return '<MolFileEntry: %s%s%s>' % (self.pdb_id, self.sep, self.mol_id)

    def __getstate__(self):
        if self._mol_obj is not None:
            self.mol_obj = MolWrapper(self.mol_obj)
        return self.__dict__

    def __setstate__(self, state):
        self.__dict__.update(state)


def recover_entries_from_entity(entity, get_small_molecules=True,
                                get_chains=False, ignore_artifacts=True,
                                by_cluster=False, sep=ENTRY_SEPARATOR):
    """ Search for chains and small molecules in ``entity`` and return
    them as strings.

    Parameters
    ----------
    entity : :class:`~luna.MyBio.PDB.Entity.Entity`
        An entity from where chains or small molecules will be recovered.
    get_small_molecules : bool
        If True, identify small molecules and return them as `MolEntry`
        objects. The default value is True.
    get_chains : bool
        If True, identify chains and return them as `ChainEntry` objects.
        The default value is False.
    ignore_artifacts : bool
        If True, ignore the following crystallography artifacts: ACE, ACT,
        BME, CSD, CSW, EDO, FMT, GOL, MSE, NAG, NO3, PO4, SGM, SO4, or
        TPO. The default value is True.
    by_cluster : bool
        If True, aggregate entries by cluster. Cluster ids are exclusive
        to :class:`~luna.MyBio.PDB.Residue.Residue` instances and are
        automatically set by :class:`~luna.MyBio.PDB.FTMapParser.FTMapParser`,
        a parser for FTMap results. By default, the cluster id of
        :class:`~luna.MyBio.PDB.Residue.Residue` instances are set to None,
        therefore, if the cluster id is not explicitly defined, all entries
        will be aggregated to the same key ``None``.
    sep : str
        A separator character to format the entry string.
        The default value is ':'.

    Returns
    ------
     : list or dict
        If ``by_cluster`` is set to False, a list of `ChainEntry` or
        `MolEntry` objects is returned. Otherwise, a dict is returned,
        in which keys are clusters and values are lists of `ChainEntry`
        or `MolEntry` objects. When no cluster information is available, all
        entries are aggregated in a key of value ``None``. Cluster ids are
        exclusive to :class:`~luna.MyBio.PDB.Residue.Residue` instances,
        therefore, `ChainEntry` objects are always placed in a key of
        value ``None``.

    Examples
    --------

    First, let's parse a PDB file.

    >>> from luna.MyBio.PDB.PDBParser import PDBParser
    >>> pdb_parser = PDBParser(PERMISSIVE=True, QUIET=True)
    >>> pdb_file = "tutorial/inputs/3QQK.pdb"
    >>> structure = pdb_parser.get_structure("Protein", pdb_file)

    Now, we can recover entries from the parsed PDB file:

    >>> from luna.mol.entry import recover_entries_from_entity
    >>> entries = recover_entries_from_entity(structure, get_chains=True)
    >>> for e in entries:
    >>>     print(e)
    <MolEntry: Protein:A:X02:497>
    <ChainEntry: Protein:A>

    """

    clusters = defaultdict(list)
    entries = []

    if entity.level == "S":
        if get_small_molecules:
            residues = entity[0].get_residues()
        if get_chains:
            chains = entity[0].get_chains()

    elif entity.level == "M":
        if get_small_molecules:
            residues = entity.get_residues()
        if get_chains:
            chains = entity.get_chains()
    else:
        if get_small_molecules:
            # If the entity is already a Chain, get_parent_by_level() returns
            # the same object. But, if the entity is a Residue or Atom, it will
            # return its corresponding chain parent.
            residues = entity.get_parent_by_level("C").get_residues()
        if get_chains:
            chains = entity.get_parent_by_level("M").get_chains()

    if get_small_molecules:
        pdb_id = entity.get_parent_by_level("S").id
        for res in residues:
            if res.is_hetatm():
                if ignore_artifacts and res.resname in ARTIFACTS_LIST:
                    continue

                comp_num_and_icode = ""
                if isinstance(res.id[1], int):
                    comp_num_and_icode = str(res.id[1])
                comp_num_and_icode += \
                    str(res.id[2]) if res.id[2].strip() else ""

                entry = MolEntry.from_string(sep.join([pdb_id, res.parent.id,
                                                       res.resname,
                                                       comp_num_and_icode]))
                if by_cluster:
                    clusters[res.cluster_id].append(entry)
                else:
                    entries.append(entry)

    if get_chains:
        pdb_id = entity.get_parent_by_level("S").id
        for chain in chains:
            entry = ChainEntry.from_string(sep.join([pdb_id, chain.id]))
            if by_cluster:
                clusters[None].append(entry)
            else:
                entries.append(entry)

    if by_cluster:
        return clusters
    return entries
