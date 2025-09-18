import logging
from io import StringIO

from Bio.PDB.PDBIO import Select

from luna.pdb.io.base import PDBIO
from luna.util.exceptions import FileNotCreated


logger = logging.getLogger()


def save_to_file(entity, 
                 output_file, 
                 select=Select(), 
                 write_conects=True,
                 write_end=True, 
                 preserve_atom_numbering=True, 
                 sort=False):
    """
    Write a structure, model, chain, residue, or atom object to a PDB file.

    Parameters
    ----------
    entity : :class:`~luna.pdb.core.entity.Entity`
         The structure, model, chain, residue, or atom to be written.
    output_file : str
        Path to the output PDB file.
    select : :class:`Bio.PDB.PDBIO.Select`
        Selector for filtering atoms/residues/chains. 
        Defaults to selecting all.
    write_conects : bool
        If True, write CONECT records (if available).
    write_end : bool
        If True, write the END record.
    preserve_atom_numbering : bool
        If True, keep the original atom serial numbers.
        If False, atoms are renumbered sequentially.
    sort : bool, default=False
        If True, sort chains and residues before writing.
    """
    try:
        io = PDBIO()
        io.set_structure(entity)
        io.save(output_file, 
                select=select,
                write_conects=write_conects,
                write_end=write_end,
                preserve_atom_numbering=preserve_atom_numbering,
                sort=sort)

    except Exception as e:
        logger.exception(e)
        raise FileNotCreated(f"PDB file '{output_file}' could not be created.")


def entity_to_string(entity, 
                     select=Select(),
                     write_conects=True,
                     write_end=True,
                     preserve_atom_numbering=True):
    """
    Convert a structure, model, chain, residue, or atom to string.

    This function works on a structural level. That means if ``entity`` is not
    a :class:`~luna.pdb.core.structure.Structure` object, the structure will
    be recovered directly from ``entity``. Therefore, use ``select`` to select
    specific chains, residues, and atoms from the structure object.

    Parameters
    ----------
    entity : :class:`~luna.pdb.core.entity.Entity`
        The object to convert.
    select : :class:`Bio.PDB.PDBIO.Select`
        Selector for filtering atoms/residues/chains. 
        Defaults to selecting all.
    write_conects : bool
        If True, write CONECT records.
    write_end : bool
        If True, write the END record.
    preserve_atom_numbering : bool
        If True, keep the original atom serial numbers.
        If False, atoms are renumbered sequentially.

    Returns
    -------
    str
        The serialized PDB content as a string.
    """
    fh = StringIO()
    io = PDBIO()
    io.set_structure(entity)
    io.save(fh, 
            select=select,
            write_conects=write_conects,
            write_end=write_end,
            preserve_atom_numbering=preserve_atom_numbering)
    return fh.getvalue()
