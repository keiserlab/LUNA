import tempfile
import logging

from Bio.PDB import Selection
from Bio.PDB.PDBIO import Select

from openbabel.pybel import readfile
from openbabel.pybel import Molecule as PybelWrapper

from rdkit.Chem import (MolFromMolBlock, SanitizeFlags,
                        SanitizeMol, MolToMolFile)


from luna.pdb.io.helpers import save_to_file
from luna.pdb.io.selector import AtomSelector
from luna.mol.validator import MolValidator
from luna.mol.standardiser import Standardizer
from luna.wrappers.base import MolWrapper
from luna.wrappers.rdkit import read_mol_from_file
from luna.wrappers.obabel import convert_molecule
from luna.util.file import new_unique_filename, remove_files
from luna.util.exceptions import (MoleculeSizeError,
                                  MoleculeObjectError)
from luna.util.default_values import OPENBABEL


logger = logging.getLogger(__name__)


def biopython_entity_to_mol(entity,
                            select=Select(),
                            amend_mol=True,
                            template=None,
                            add_h=False,
                            ph=None,
                            metals_coord=None,
                            wrapped=True,
                            openbabel=OPENBABEL,
                            tmp_path=None,
                            keep_tmp_files=False):
    """Convert an object :class:`~luna.pdb.core.entity.Entity` to a
    molecular object (:class:`~luna.wrappers.base.MolWrapper`,
    :class:`rdkit.Chem.rdchem.Mol`, or :class:`openbabel.pybel.Molecule`).

    Parameters
    ----------
    entity : :class:`~luna.pdb.core.entity.Entity`
        The entity to be converted.
    select : :class:`Bio.PDB.PDBIO.Select`
        Decide which atoms will be considered.
        By default, all atoms are accepted.
    amend_mol : bool
        If True, validate and standardize molecules.
        The validation is employed with
        :class:`~luna.mol.validator.MolValidator`.
        Currently, only residues are standardized.
    template : :class:`luna.mol.template.Template`, optional
        Fix the converted molecule's bond order based on the bond orders
        in a template molecule. The ``template`` should implement
        assign_bond_order().
    add_h : bool
        If True, add hydrogen to the converted molecule.
    ph : float, optional
        Add hydrogens considering pH ``ph``.
    wrapped : bool
        If True, wrap the molecular object with
        :class:`~luna.wrappers.base.MolWrapper`.
    openbabel : str
        Pathname to Open Babel.
    tmp_path : str
        A temporary directory to where temporary files will be saved.
        If not provided, the system's default temporary directory will
        be used instead.
    keep_tmp_files : bool
        If True, keep all temporary files. Otherwise, removes them in the end.

    Returns
    -------
    mol_obj : :class:`~luna.wrappers.base.MolWrapper`, \
                :class:`rdkit.Chem.rdchem.Mol`, \
                or :class:`openbabel.pybel.Molecule`
        The converted molecule.
    ignored_atoms : list
        List of ignored atoms. Currently, ignored atoms may
        contain only metals.
    """

    tmp_path = tmp_path or tempfile.gettempdir()

    logger.debug("It will try to create a new MOL object from the provided "
                 "entity.")
    logger.debug("Temporary files will be saved at '%s'." % tmp_path)

    # First it saves the selection into a PDB file and then it converts the
    # file to .mol. I had to do it because the OpenBabel 2.4.1 had a
    # problem with some molecules containing aromatic rings. In such cases,
    # the aromatic ring was being wrongly perceived and some atoms received
    # more double bonds than it was expected. The version 2.3.2 works better.
    # Therefore, I recomend using Open Babel 2.3.3 instead.
    filename = new_unique_filename(tmp_path)
    pdb_file = '%s_pdb-file.pdb' % filename

    logger.debug("First: it will try to create a new PDB file (%s) "
                 "from the provided entity." % pdb_file)

    # Apparently, Open Babel creates a bug when it tries to parse a file with
    # CONECTS containing serial numbers with more than 4 digits.
    # E.g.: 1OZH:A:HE3:1406, line CONECT162811627916282.
    # By setting preserve_atom_numbering to False, it solves the problem.
    save_to_file(entity,
                 pdb_file,
                 select,
                 preserve_atom_numbering=False,
                 sort=True)

    ini_input_file = pdb_file
    if template is not None:
        if entity.level == "R" and entity.is_hetatm():
            # Note that the template molecule should have no explicit hydrogens
            # else the algorithm will fail.
            rdmol = read_mol_from_file(pdb_file, mol_format="pdb",
                                       removeHs=True)
            new_rdmol = template.assign_bond_order(rdmol, entity.resname)

            ini_input_file = '%s_tmp-mol-file.mol' % filename
            MolToMolFile(new_rdmol, ini_input_file)

            if not keep_tmp_files:
                remove_files([pdb_file])
        else:
            logger.warning("It cannot apply a template on the provided entity "
                           "because it should be a single compound "
                           "(Residue class).")

    # Convert the PDB file to Mol file with the proper protonation
    # and hydrogen addition if required.
    mol_file = '%s_mol-file.mol' % filename
    ob_opt = {"error-level": 5}
    logger.debug("Next: it will try to convert the PDB file to "
                 " .mol using Open Babel.")
    if add_h:
        logger.debug("Hydrogens will be added to the molecule.")
        if ph is not None:
            ob_opt["p"] = ph
        else:
            ob_opt["h"] = ""
    convert_molecule(ini_input_file, output_file=mol_file,
                     opts=ob_opt, openbabel=openbabel)

    # Currently, ignored atoms are only metals.
    ignored_atoms = []

    mol_obj = None
    if amend_mol:
        logger.debug("Now: a validation will be performed and "
                     "it will try to fix some errors.")

        try:
            mol_obj = next(readfile("mol", mol_file))
        except Exception:
            error_msg = ("An error occurred while parsing the file '%s' with "
                         "Open Babel and the molecule object could not be "
                         "created. Check the logs for more information."
                         % mol_file)
            raise MoleculeObjectError(error_msg)

        # Standardize a specific set of atoms if provided,
        # otherwise use all atoms from entity.
        if isinstance(select, AtomSelector):
            target_atoms = select.entries
        else:
            target_atoms = Selection.unfold_entities(entity, 'A')

        # If the add_h property is set to False, the code will not remove
        # any existing hydrogens from the PDB structure.
        # In these situations, the list of atoms may contain hydrogens.
        # But, we do not need to attribute properties to hydrogens.
        # We just need them to correctly set properties to heavy atoms.
        # So let's just ignore them.
        target_atoms = [atm for atm in target_atoms
                        if select.accept_atom(atm) and atm.element != "H"]

        mol_obj = MolWrapper(mol_obj)
        if mol_obj.get_num_heavy_atoms() != len(target_atoms):
            error_msg = ("The number of heavy atoms in the PDB selection "
                         "and in the MOL file are different.")
            raise MoleculeSizeError(error_msg)

        # Ignore hydrogen atoms.
        atm_obj_list = [atm for atm in mol_obj.get_atoms()
                        if atm.get_atomic_num() != 1]

        # Pairs of AtomWrapper and Atom (Biopython) objects
        # to standardize.
        atom_pairs = []
        pdb_mapping = {}
        for i, atm_obj in enumerate(atm_obj_list):
            atom_pairs.append((atm_obj, target_atoms[i]))

            pdb_mapping[atm_obj.get_idx()] = target_atoms[i]

        updated_metals_coord = {}
        # If there is something to be standardised.
        if atom_pairs:
            rs = Standardizer()
            mol_obj.unwrap().DeleteHydrogens()
            rs.standardize(atom_pairs, metals_coord=metals_coord)

            updated_metals_coord = rs.metals_coord

            # After standardizing residues, we need to recreate the Mol
            # file, otherwise implicit hydrogens will not be included
            # in the MOL object and, therefore, their coordinates could
            # not be accessed. If you try to generate coordinates directly
            # from the object, hydrogens will be incorrectly placed.
            mol_obj = PybelWrapper(mol_obj.unwrap())
            new_mol_file = '%s_tmp-mol-file.mol' % filename
            mol_obj.write("mol", new_mol_file, overwrite=True)

            # Overwrite mol_file by converting the new molecular file using
            # the user specified parameters. Note that right now it will add
            # explicit hydrogens to the molecules according to the provided pH.
            convert_molecule(new_mol_file, output_file=mol_file, 
                             opts=ob_opt, openbabel=OPENBABEL)

            # Let's finally read the correct and standardized molecular file.
            try:
                mol_obj = next(readfile("mol", mol_file))
            except Exception:
                error_msg = ("An error occurred while parsing the file "
                             "'%s' with Open Babel and the molecule object"
                             "could not be created. Check the logs for "
                             "more information." % mol_file)
                raise MoleculeObjectError(error_msg)

            # Remove temporary files.
            if not keep_tmp_files:
                remove_files([new_mol_file])

        mv = MolValidator(metals_coord=updated_metals_coord)
        is_valid = mv.validate_mol(mol_obj, pdb_mapping)
        logger.debug('Validation finished!!!')

        if not is_valid:
            logger.warning("The molecular file '%s' contain invalid atoms. "
                           "Check the logs for more information." % mol_file)

        # Validate molecule using RDKit sanitization methods.
        try:
            aux_mol = (PybelWrapper(mol_obj.unwrap())
                       if isinstance(mol_obj, MolWrapper)
                       else mol_obj)

            # The sanitization is set off. We will apply it in the next step.
            aux_mol = MolFromMolBlock(aux_mol.write('mol'),
                                      sanitize=False, removeHs=False)
            # Sanitize molecule is applied now, so we will be able to catch
            # the exceptions raised by RDKit, otherwise it would
            # not be possible.
            SanitizeMol(aux_mol, SanitizeFlags.SANITIZE_ALL)
        except Exception:
            error_msg = ("An error occurred while parsing the molecular block "
                         "with RDKit. The block was generated by Open Babel "
                         "from the file '%s'. Check the logs for more "
                         "information." % mol_file)
            raise MoleculeObjectError(error_msg)
    else:
        try:
            # Create a new Mol object.
            mol_obj = next(readfile("mol", mol_file))
        except Exception:
            error_msg = ("An error occurred while parsing the file '%s' and "
                         "the molecule object could not be created. "
                         "Check the logs for more information."
                         % mol_file)
            raise MoleculeObjectError(error_msg)

    # Remove temporary files.
    if not keep_tmp_files:
        remove_files([ini_input_file, mol_file])

    if wrapped:
        mol_obj = MolWrapper(mol_obj)
    elif isinstance(mol_obj, MolWrapper):
        mol_obj = mol_obj.unwrap()

    return mol_obj, ignored_atoms