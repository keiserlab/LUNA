from itertools import product

from luna.util.default_values import (COV_SEARCH_RADIUS, BOUNDARY_CONFIG)
from luna.util.exceptions import EntityLevelError, IllegalArgumentError
from luna.interaction.cov import is_covalently_bound
from luna.MyBio.util import ENTITY_LEVEL_NAME
from luna.MyBio.PDB.NeighborSearch import NeighborSearch
from luna.MyBio.PDB import Selection
from luna.MyBio.PDB.Residue import Residue

import logging

logger = logging.getLogger()


def get_all_contacts(entity,
                     radius=BOUNDARY_CONFIG["bsite_cutoff"],
                     level='A'):
    """Recover all residue-residue or atom-atom contacts in ``entity``.

    Parameters
    ----------
    entity : :class:`~luna.MyBio.PDB.Entity.Entity`
        The PDB object from where atoms and residues will be recovered.
    radius : float
        The cutoff distance (in Å) for defining contacts.
        The default value is 6.2.
    level : {'R', 'A'}
        Return residues ('R') or atoms ('A') in contact with ``source``.

    Returns
    -------
     : list of tuple of (:class:`~luna.MyBio.PDB.Residue.Residue` or \
                :class:`~luna.MyBio.PDB.Atom.Atom`, \
                :class:`~luna.MyBio.PDB.Residue.Residue` or \
                :class:`~luna.MyBio.PDB.Atom.Atom`)
        Each tuple contains either a pair of residues or atoms in contact.

    Raises
    ------
    EntityLevelError
        If ``level`` is neither 'R' nor 'A'.

    Examples
    --------

    In this example, we will identify all residue-residue contacts
    within 2.5 Å in the PDB 3QQK. So, let's first parse the PDB file.

    >>> from luna.util.default_values import LUNA_PATH
    >>> from luna.MyBio.PDB.PDBParser import PDBParser
    >>> pdb_parser = PDBParser(PERMISSIVE=True, QUIET=True)
    >>> structure = pdb_parser.get_structure("Protein", \
f"{LUNA_PATH}/tutorial/inputs/3QQK.pdb")

    Then, to recover all residue-residue contacts within 2.5 Å use
    :meth:`get_all_contacts` with ``level`` set to 'R' and radius set to 2.5.

    >>> from luna.interaction.contact import get_all_contacts
    >>> contacts = get_all_contacts(structure, radius=2.5, level="R")
    >>> print(len(contacts))
    314

    """
    try:
        if level == 'S' or level == 'M' or level == 'C':
            raise EntityLevelError("Maximum entity level to be chosen "
                                   "is: R (residues) or A (atoms)")

        if level not in ENTITY_LEVEL_NAME:
            raise EntityLevelError("The defined level '%s' does not exist"
                                   % level)

        logger.debug("Trying to select all contacts in the PDB file %s."
                     % entity.get_parent_by_level('S').id)

        all_atoms = list(entity.get_atoms())
        ns = NeighborSearch(all_atoms)
        pairs = ns.search_all(radius, level)

        logger.debug("Number of nearby %s(s) found: %d."
                     % (ENTITY_LEVEL_NAME[level].lower(), len(pairs)))

        return pairs
    except Exception as e:
        logger.exception(e)
        raise


def get_contacts_with(source,
                      target=None,
                      entity=None,
                      radius=BOUNDARY_CONFIG["bsite_cutoff"],
                      level='A'):
    """Recover atoms or residues in contact with ``source``.

    Parameters
    ----------
    source : :class:`~luna.MyBio.PDB.Entity.Entity`
        The reference, which can be any :class:`~luna.MyBio.PDB.Entity.Entity`
        instance (structure, model, chain, residue, or atom).
    target : :class:`~luna.MyBio.PDB.Entity.Entity`, optional
        If provided, only contacts with the ``target`` will be considered.
    entity : :class:`~luna.MyBio.PDB.Entity.Entity`
        The PDB object from where atoms will be recovered.
        If not provided (the default), the model object that contains
        ``source`` will be used instead.
    radius : float
        The cutoff distance (in Å) for defining contacts.
        The default value is 6.2.
    level : {'R', 'A'}
        Return residues ('R') or atoms ('A') in contact with ``source``.

    Returns
    -------
     : set of tuple of (:class:`~luna.MyBio.PDB.Residue.Residue` or \
            :class:`~luna.MyBio.PDB.Atom.Atom`, \
            :class:`~luna.MyBio.PDB.Residue.Residue` or \
            :class:`~luna.MyBio.PDB.Atom.Atom`)
        Each tuple contains two items:\
            the first corresponds to a residue/atom from the ``source``,\
            and the second corresponds to a residue/atom in contact
            with ``source``.

    Raises
    ------
    EntityLevelError
        If ``level`` is neither 'R' nor 'A'.

    Examples
    --------

    **Example 1)** In this example, we will identify residue-residue contacts
    between a ligand and nearby residues.

    First, let's parse a PDB file to work with.

    >>> from luna.util.default_values import LUNA_PATH
    >>> from luna.MyBio.PDB.PDBParser import PDBParser
    >>> pdb_parser = PDBParser(PERMISSIVE=True, QUIET=True)
    >>> structure = pdb_parser.get_structure("Protein", \
f"{LUNA_PATH}/tutorial/inputs/3QQK.pdb")

    Now, we define the target ligand.

    >>> ligand = structure[0]["A"][('H_X02', 497, ' ')]

    Then, to recover contacts between this ligand and nearby residues we use
    :meth:`get_contacts_with` with ``level`` set to 'R'.

    >>> from luna.interaction.contact import get_contacts_with
    >>> contacts = sorted(get_contacts_with(ligand,
    ...                                     entity=structure,
    ...                                     radius=3,
    ...                                     level="R"))
    >>> for pair in contacts:
    ...     print(pair)
    (<Residue X02 het=H_X02 resseq=497 icode= >, \
<Residue GLU het=  resseq=81 icode= >)
    (<Residue X02 het=H_X02 resseq=497 icode= >, \
<Residue LEU het=  resseq=83 icode= >)
    (<Residue X02 het=H_X02 resseq=497 icode= >, \
<Residue X02 het=H_X02 resseq=497 icode= >)
    (<Residue X02 het=H_X02 resseq=497 icode= >, \
<Residue HOH het=W resseq=321 icode= >)


    **Example 2)** In this example, we will identify atom-atom contacts
    between a ligand and a given residue.

    First, let's parse a PDB file to work with.

    >>> from luna.util.default_values import LUNA_PATH
    >>> from luna.MyBio.PDB.PDBParser import PDBParser
    >>> pdb_parser = PDBParser(PERMISSIVE=True, QUIET=True)
    >>> structure = pdb_parser.get_structure("Protein", \
f"{LUNA_PATH}/tutorial/inputs/3QQK.pdb")

    Now, we define the target ligand and residue.

    >>> ligand = structure[0]["A"][('H_X02', 497, ' ')]
    >>> residue = structure[0]["A"][(' ', 81, ' ')]

    Then, to recover contacts between these compounds we use
    :meth:`get_contacts_with`. As we need atom-wise contacts, set ``level``
    to 'A'.

    >>> from luna.interaction.contact import get_contacts_with
    >>> contacts = sorted(get_contacts_with(ligand,
    ...                   target=residue,
    ...                   entity=structure,
    ...                   radius=4,
    ...                   level="A"))
    >>> for pair in contacts:
    ...     print(pair)
    (<Atom C7>, <Atom O>)
    (<Atom N10>, <Atom C>)
    (<Atom N10>, <Atom O>)
    """
    try:
        if level == 'S' or level == 'M' or level == 'C':
            raise EntityLevelError("Maximum entity level to be chosen is: "
                                   "R (residues) or A (atoms)")

        if level not in ENTITY_LEVEL_NAME:
            raise EntityLevelError("The defined level '%s' does not exist"
                                   % level)

        entity = entity or source.get_parent_by_level("M")

        source_atoms = Selection.unfold_entities([source], 'A')

        target_atoms = []
        if target is None:
            target_atoms = list(entity.get_atoms())
        else:
            if target.level == "A":
                target_atoms = [target]
            else:
                target_residues = Selection.unfold_entities(target, 'R')
                target_atoms = [a for r in target_residues
                                for a in r.get_unpacked_list()]

        ns = NeighborSearch(target_atoms)
        entities = set()
        for atom in source_atoms:
            entity = atom.get_parent_by_level(level)
            nb_entities = ns.search(atom.coord, radius, level)
            pairs = set(product([entity], nb_entities))
            entities.update(pairs)

        logger.debug("Number of nearby %s(s) found: %d."
                     % (ENTITY_LEVEL_NAME[level].lower(), len(entities)))
        return entities
    except Exception as e:
        logger.exception(e)
        raise


def get_proximal_compounds(source, radius=COV_SEARCH_RADIUS):
    """Recover proximal compounds to ``source``.

    Parameters
    ----------
    source : :class:`~luna.MyBio.PDB.Residue.Residue`
        The reference compound.
    radius : float
        The cutoff distance (in Å) for defining proximity.
        The default value is 2.2, which may recover potential residues bound to
        ``source`` through covalent bonds.

    Returns
    -------
     : list of :class:`~luna.MyBio.PDB.Residue.Residue`
        The list of proximal compounds always include ``source``.

    Raises
    ------
    IllegalArgumentError
        If ``source`` is not a :class:`~luna.MyBio.PDB.Residue.Residue`

    Examples
    --------

    In this example, we will identify all proximal compounds to a given residue
    in the PDB 3QQK.
    So, let's first parse the PDB file.

    >>> from luna.util.default_values import LUNA_PATH
    >>> from luna.MyBio.PDB.PDBParser import PDBParser
    >>> pdb_parser = PDBParser(PERMISSIVE=True, QUIET=True)
    >>> structure = pdb_parser.get_structure("Protein", \
f"{LUNA_PATH}/tutorial/inputs/3QQK.pdb")

    Now, we select the residue of our interest.

    >>> residue = structure[0]["A"][(' ', 81, ' ')]

    Finally, we call :meth:`get_proximal_compounds`.

    >>> from luna.interaction.contact import get_proximal_compounds
    >>> compounds = get_proximal_compounds(residue)
    >>> for c in compounds:
    ...    print(c)
    <Residue PHE het=  resseq=80 icode= >
    <Residue GLU het=  resseq=81 icode= >
    <Residue PHE het=  resseq=82 icode= >
    """
    if not isinstance(source, Residue):
        raise IllegalArgumentError("Invalid object type. "
                                   "A Residue object is expected instead.")

    model = source.get_parent_by_level('M')
    proximal = get_contacts_with(source, entity=model,
                                 radius=radius, level='R')

    # Sorted by the compound order as in the PDB.
    return sorted(list(set([p[1] for p in proximal])),
                  key=lambda r: (r.parent.parent.id, r.parent.id, r.idx))


def get_cov_contacts_with(source, target=None, entity=None):
    """Recover potential covalent bonds with ``source``.

    Covalent bonds between two nearby atoms A and B are determined
    as in Open Babel:

    .. math:: 0.4 <= \overrightarrow{\|AB\|} <= A_{cov} + B_{cov} + 0.45

    Where :math:`\overrightarrow{\|AB\|}` is the distance
    between atoms A and B, and :math:`A_{cov}` and :math:`B_{cov}`
    are the covalent radius of atoms A and B, respectively.

    Parameters
    ----------
    source : :class:`~luna.MyBio.PDB.Entity.Entity`
        The reference, which can be any :class:`~luna.MyBio.PDB.Entity.Entity`
        instance (structure, model, chain, residue, or atom).
    target : :class:`~luna.MyBio.PDB.Entity.Entity`, optional
        If provided, only covalent bonds with the ``target``
        will be considered.
    entity : :class:`~luna.MyBio.PDB.Entity.Entity`
        The PDB object from where atoms will be recovered.
        If not provided (the default), the model object that contains
        ``source`` will be used instead.

    Returns
    -------
     : set of tuple of (:class:`~luna.MyBio.PDB.Atom.Atom`, \
            :class:`~luna.MyBio.PDB.Atom.Atom`)
        Pairs of atoms with potential covalent bonds.


    Examples
    --------

    In this example, we will identify all covalent bonds involving a given
    residue in the PDB 3QQK. So, let's first parse the PDB file.

    >>> from luna.util.default_values import LUNA_PATH
    >>> from luna.MyBio.PDB.PDBParser import PDBParser
    >>> pdb_parser = PDBParser(PERMISSIVE=True, QUIET=True)
    >>> structure = pdb_parser.get_structure("Protein", \
f"{LUNA_PATH}/tutorial/inputs/3QQK.pdb")

    Now, we select the residue of our interest.

    >>> residue = structure[0]["A"][(' ', 81, ' ')]

    Finally, let`s recover potential covalent bonds involving this residue
    with :meth:`get_cov_contacts_with`. In the below snippet, pairs of atoms
    are sorted by atoms` serial number, and the residue information is printed
    together with the atom name.

    >>> from luna.interaction.contact import get_cov_contacts_with
    >>> cov_bonds = sorted(get_cov_contacts_with(residue, entity=structure),
    ...                    key=lambda x: (x[0].serial_number,
    ...                                   x[1].serial_number))
    >>> for atm1, atm2 in cov_bonds:
    >>>     pair = ("%s%d/%s" % (atm1.parent.resname,
    ...                          atm1.parent.id[1],
    ...                          atm1.name),
    ...             "%s%d/%s" % (atm2.parent.resname,
    ...                          atm2.parent.id[1],
    ...                          atm2.name))
    >>>     print(pair)
    ('PHE80/C', 'GLU81/N')
    ('GLU81/N', 'GLU81/CA')
    ('GLU81/CA', 'GLU81/C')
    ('GLU81/CA', 'GLU81/CB')
    ('GLU81/C', 'GLU81/O')
    ('GLU81/C', 'PHE82/N')
    ('GLU81/CB', 'GLU81/CG')
    ('GLU81/CG', 'GLU81/CD')
    ('GLU81/CD', 'GLU81/OE1')
    ('GLU81/CD', 'GLU81/OE2')
    """

    entity = entity or source.get_parent_by_level("M")
    entities = get_contacts_with(source,
                                 target=target,
                                 entity=entity,
                                 radius=COV_SEARCH_RADIUS,
                                 level='A')

    cov_bonds = set()
    for atm1, atm2 in entities:
        if is_covalently_bound(atm1, atm2):
            cov_bonds.add(tuple(sorted([atm1, atm2],
                                       key=lambda x: x.serial_number)))

    return cov_bonds
