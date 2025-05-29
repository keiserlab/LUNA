"""Consumer class that builds a Structure object.

This is used by the PDBParser and MMCIFparser classes.
"""

import warnings

from Bio.PDB.PDBExceptions import PDBConstructionException
from Bio.PDB.PDBExceptions import PDBConstructionWarning

# SMCRA hierarchy
from luna.pdb.core.structure import Structure
from luna.pdb.core.model import Model
from luna.pdb.core.chain import Chain
from luna.pdb.core.residue import Residue, DisorderedResidue
from luna.pdb.core.atom import Atom, DisorderedAtom


class StructureBuilder(object):
    """
    A Biopython-compatible structure builder for LUNA.

    Responsible for constructing the SMCRA hierarchy (Structure, Model, Chain,
    Residue, Atom) from parsed data. Includes additional features like:

    - CONECT record storage
    - Residue merging logic (for OpenBabel cleanup)
    - Disordered residue and atom handling
    - Residue indexing (idx) and line tracking
    """

    def __init__(self):
        self.line_counter = 0
        self.header = {}

        # StructureBuilder now keeps CONECT records.
        self.conects = {}

    def _is_completely_disordered(self, residue):
        """Return 1 if all atoms in the residue have a non blank altloc."""
        atom_list = residue.get_unpacked_list()
        for atom in atom_list:
            altloc = atom.get_altloc()
            if altloc == " ":
                return 0
        return 1

    # Public methods called by the Parser classes

    def set_header(self, header):
        """
        Set the PDB header block.

        Parameters
        ----------
        header : dict
            Header dictionary parsed from the PDB file.
        """
        self.header = header

    def set_conects(self, conects):
        """
        Store CONECT records parsed from the PDB file.

        Parameters
        ----------
        conects : dict
            Mapping of atom serial numbers to bonded atom serial numbers.
        """
        self.conects = conects

    def set_line_counter(self, line_counter):
        """
        Track the current line being parsed (for error reporting).

        Parameters
        ----------
        line_counter : int
            Current PDB file line number.
        """
        self.line_counter = line_counter

    def init_structure(self, structure_id, pdb_file=None):
        """
        Initialize a new Structure object.

        Parameters
        ----------
        structure_id : str
            Unique identifier for the structure.
        """
        self.structure = Structure(structure_id)

    def init_model(self, model_id, serial_num=None):
        """
        Initialize a new Model object and add it to the structure.

        Parameters
        ----------
        model_id : int
            Model index.
        serial_num : int, optional
            PDB serial number.
        """
        self.model = Model(model_id, serial_num)
        self.structure.add(self.model)

    def init_chain(self, chain_id):
        """
        Initialize a new Chain object and add it to the current model.

        Parameters
        ----------
        chain_id : str
            Chain identifier.
        """
        if self.model.has_id(chain_id):
            self.chain = self.model[chain_id]
            warnings.warn("WARNING: Chain %s is discontinuous at line %i."
                          % (chain_id, self.line_counter),
                          PDBConstructionWarning)
        else:
            self.chain = Chain(chain_id)
            self.model.add(self.chain)

    def init_seg(self, segid):
        """
        Set the current segment ID.

        Parameters
        ----------
        segid : str
            Segment ID (used by PDB files).
        """
        self.segid = segid

    def init_residue(self, resname, field, resseq, icode):
        """
        Initialize a new Residue (or DisorderedResidue) object.

        Merges with existing residue if applicable, handles alternate conformations
        and OpenBabel artifacts.

        Parameters
        ----------
        resname : str
            Residue name (e.g., "GLY").
        field : str
            Hetero flag (" ", "W", or "H" / "H_"...).
        resseq : int
            Residue sequence number.
        icode : str
            Insertion code.
        """
        if field != " ":
            if field == "H":
                # The hetero field consists of H_ + the residue name (e.g. H_FUC)
                field = "H_" + resname
        res_id = (field, resseq, icode)

        if field == " ":
            if self.chain.has_id(res_id):
                # There already is a residue with the id (field, resseq, icode).
                # This only makes sense in the case of a point mutation.
                warnings.warn("WARNING: Residue ('%s', %i, '%s') "
                              "redefined at line %i."
                              % (field, resseq, icode, self.line_counter),
                              PDBConstructionWarning)
                duplicate_residue = self.chain[res_id]
                if duplicate_residue.is_disordered() == 2:
                    # The residue in the chain is a DisorderedResidue object.
                    # So just add the last Residue object.
                    if duplicate_residue.disordered_has_id(resname):
                        # The residue was already made
                        self.residue = duplicate_residue
                        duplicate_residue.disordered_select(resname)
                    else:
                        # Make a new residue and add it to the already
                        # present DisorderedResidue
                        #
                        # Create a new residue with the new property (idx) added to the Residue class.
                        new_residue = Residue(res_id, resname, self.segid, len(self.chain.child_list), self.line_counter)
                        duplicate_residue.disordered_add(new_residue)
                        self.residue = duplicate_residue
                        return
                else:
                    if resname == duplicate_residue.resname:
                        warnings.warn("WARNING: Residue ('%s', %i, '%s','%s')"
                                      " already defined with the same name at line  %i."
                              % (field, resseq, icode, resname, self.line_counter),
                              PDBConstructionWarning)
                        self.residue = duplicate_residue
                        return
                    # Make a new DisorderedResidue object and put all
                    # the Residue objects with the id (field, resseq, icode) in it.
                    # These residues each should have non-blank altlocs for all their atoms.
                    # If not, the PDB file probably contains an error.
                    if not self._is_completely_disordered(duplicate_residue):
                        # if this exception is ignored, a residue will be missing
                        self.residue = None
                        raise PDBConstructionException(
                            "Blank altlocs in duplicate residue %s ('%s', %i, '%s')"
                            % (resname, field, resseq, icode))
                    self.chain.detach_child(res_id)

                    # Create a new residue with the new property (idx) added to the Residue class.
                    new_residue = Residue(res_id, resname, self.segid, len(self.chain.child_list), self.line_counter)
                    disordered_residue = DisorderedResidue(res_id)
                    self.chain.add(disordered_residue)
                    disordered_residue.disordered_add(duplicate_residue)
                    disordered_residue.disordered_add(new_residue)
                    self.residue = disordered_residue
                    return
        else:
            # Identify hetero residues already added to a chain.
            if self.chain.has_id(res_id):
                # There already is a hetero residue with the id (field, resseq, icode).
                # It can occur when hydrogen atoms are added in the final of the PDB file.
                warnings.warn("WARNING: Residue ('%s', %i, '%s') "
                              "redefined at line %i."
                              % (field, resseq, icode, self.line_counter),
                              PDBConstructionWarning)
                duplicate_residue = self.chain[res_id]
                if resname == duplicate_residue.resname:
                    warnings.warn("WARNING: Residue ('%s', %i, '%s','%s')"
                                  " already defined with the same name at line  %i."
                          % (field, resseq, icode, resname, self.line_counter),
                          PDBConstructionWarning)
                    self.residue = duplicate_residue
                    return

        # Create a new residue with the new property (idx) added to the Residue class.
        self.residue = Residue(res_id, resname, self.segid, len(self.chain.child_list), self.line_counter)
        self.chain.add(self.residue)

    def init_atom(self, name, coord, b_factor, occupancy, altloc, fullname,
                  serial_number=None, element=None):
        """
        Initialize a new Atom or DisorderedAtom object and add it to the current residue.

        Parameters
        ----------
        name : str
            Atom name (e.g., "CA").
        coord : np.ndarray
            Cartesian coordinates.
        b_factor : float
            B-factor.
        occupancy : float
            Occupancy value.
        altloc : str
            Alternate location indicator.
        fullname : str
            Atom name with spaces (e.g., " CA ").
        serial_number : int, optional
            Atom serial number.
        element : str, optional
            Atomic element.
        """
        residue = self.residue
        
        # if residue is None, an exception was generated during
        # the construction of the residue
        if residue is None:
            return
            
        # First check if this atom is already present in the residue.
        # If it is, it might be due to the fact that the two atoms have atom
        # names that differ only in spaces (e.g. "CA.." and ".CA.",
        # where the dots are spaces). If that is so, use all spaces
        # in the atom name of the current atom.
        if residue.has_id(name):
            duplicate_atom = residue[name]
            # atom name with spaces of duplicate atom
            duplicate_fullname = duplicate_atom.get_fullname()
            if duplicate_fullname != fullname:
                # name of current atom now includes spaces
                name = fullname
                warnings.warn("Atom names %r and %r differ "
                              "only in spaces at line %i."
                              % (duplicate_fullname, fullname,
                                 self.line_counter),
                              PDBConstructionWarning)
                
        self.atom = Atom(name, coord, b_factor, occupancy, altloc,
                         fullname, serial_number, element)
        
        if altloc != " ":
            # The atom is disordered
            if residue.has_id(name):
                # Residue already contains this atom
                duplicate_atom = residue[name]
                if duplicate_atom.is_disordered() == 2:
                    duplicate_atom.disordered_add(self.atom)
                else:
                    # This is an error in the PDB file:
                    # a disordered atom is found with a blank altloc
                    # Detach the duplicate atom, and put it in a
                    # DisorderedAtom object together with the current
                    # atom.
                    residue.detach_child(name)
                    disordered_atom = DisorderedAtom(name)
                    residue.add(disordered_atom)
                    disordered_atom.disordered_add(self.atom)
                    disordered_atom.disordered_add(duplicate_atom)
                    residue.flag_disordered()
                    warnings.warn("WARNING: disordered atom found "
                                  "with blank altloc before line %i.\n"
                                  % self.line_counter,
                                  PDBConstructionWarning)
            else:
                # The residue does not contain this disordered atom
                # so we create a new one.
                disordered_atom = DisorderedAtom(name)
                residue.add(disordered_atom)
                # Add the real atom to the disordered atom, and the
                # disordered atom to the residue
                disordered_atom.disordered_add(self.atom)
                residue.flag_disordered()
        else:
            # The atom is not disordered
            residue.add(self.atom)

    def set_anisou(self, anisou_array):
        """
        Set anisotropic B-factor for the current atom.

        Parameters
        ----------
        anisou_array : np.ndarray
            Anisotropic B-factor.
        """
        self.atom.set_anisou(anisou_array)

    def set_siguij(self, siguij_array):
        """
        Set standard deviation of anisotropic B-factor for current atom.

        Parameters
        ----------
        siguij_array : np.ndarray
        """
        self.atom.set_siguij(siguij_array)

    def set_sigatm(self, sigatm_array):
        """
        Set standard deviation of atomic positions for current atom.

        Parameters
        ----------
        sigatm_array : np.ndarray
        """
        self.atom.set_sigatm(sigatm_array)

    def get_structure(self):
        """
        Finalize and return the constructed Structure.

        Returns
        -------
        Bio.PDB.Structure
            Structure object populated from parsed data.
        """
        # Add the header dict
        self.structure.header = self.header

        # Add the conects dict.
        self.structure.conects = self.conects

        return self.structure

    def set_symmetry(self, spacegroup, cell):
        pass
