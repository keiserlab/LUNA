import warnings
import logging

from Bio.File import as_handle
from Bio.PDB.PDBExceptions import PDBConstructionException
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from luna.pdb.parser.base import PDBParser
from luna.pdb.builder import StructureBuilder


logger = logging.getLogger()

PROBES_SET = ("ACD", "ACN", "ACT", "ADY", "AMN", "BDY", 
              "BEN", "BUT", "CHX", "DFO", "DME", "EOL", 
              "ETH", "PHN", "THS", "URE")


class FTMapParser(PDBParser):
    """
    Parser for FTMap PDB output.

    FTMap generates multiple `HEADER`-delimited chunks representing protein,
    clusters, and probe molecules. This parser merges them into a single Biopython
    Structure and performs:
      - optional filtering of clusters (`only_clusters`)
      - optional filtering of probes/residues (`only_compounds`)
      - serial number remapping
      - CONECT remapping
      - probe-to-ligand transformation
    """

    def __init__(self, 
                 PERMISSIVE=True, 
                 QUIET=False,
                 FIX_ATOM_NAME_CONFLICT=False,
                 FIX_OBABEL_FLAGS=False,
                 FIX_ALL_OBABEL_ERRORS=False):
        """
        Initialize an FTMapParser instance.

        Parameters
        ----------
        PERMISSIVE : bool
            Allow recovery from malformed PDB files.
            If False, exceptions in constructing the SMCRA data structure are fatal. 
            If true (DEFAULT), the exceptions are caught, but some residues or atoms 
            will be missing.
        QUIET : bool
            Suppress warnings during parsing.
        FIX_ATOM_NAME_CONFLICT : bool
            Fix conflicting atom names (e.g., hydrogens added by OpenBabel).
        FIX_OBABEL_FLAGS : bool
            Fix OpenBabel-added atom flags.
        FIX_ALL_OBABEL_ERRORS : bool
            Enable both atom name and flag corrections.
        """
        self.PERMISSIVE = bool(PERMISSIVE)
        self.QUIET = bool(QUIET)

        if FIX_ALL_OBABEL_ERRORS:
            FIX_ATOM_NAME_CONFLICT = True
            FIX_OBABEL_FLAGS = True
        
        # Correct conflicts in atom names.
        self.FIX_ATOM_NAME_CONFLICT = FIX_ATOM_NAME_CONFLICT
        # Do not fix empty chains using the default method. 
        # Instead, FTMapParser implements its own method to fix empty chains.
        self.FIX_EMPTY_CHAINS = False
        # Correct flags that were incorrectly set by OpenBabel.
        self.FIX_OBABEL_FLAGS = FIX_OBABEL_FLAGS
        # Correct all OpenBabel errors.
        self.FIX_ALL_OBABEL_ERRORS = FIX_ALL_OBABEL_ERRORS

    def get_structure(self, id, file, only_clusters=None, only_compounds=None):
        """
        Return a Biopython Structure object parsed from FTMap output.

        Parameters
        ----------
        id : str
            Structure ID.
        file : str or handle
            PDB file path or open handle.
        only_clusters : list of int, optional
            Limit parsing to specified FTMap cluster IDs.
        only_compounds : list of residue IDs, optional
            Limit output to specific residues.

        Returns
        -------
        Structure
        """
        with warnings.catch_warnings():
            if self.QUIET:
                warnings.filterwarnings("ignore", category=PDBConstructionWarning)

            with as_handle(file, mode='rU') as handle:
                self._parse_chunks(handle.readlines(), 
                                   id, 
                                   file,
                                   only_clusters=only_clusters,
                                   only_compounds=only_compounds)

            # Return the Structure instance
            return self.structure_builder.get_structure()

    def _parse_chunks(self, all_lines, id, file, only_clusters=None, only_compounds=None):
        """
        Parse all HEADER-separated chunks from an FTMap PDB file.

        This builds and merges partial structures using a temporary StructureBuilder.

        Parameters
        ----------
        lines : list of str
            Full PDB file lines.
        id : str
            Structure ID.
        file : str
            PDB file path.
        only_clusters : list of int, optional
        only_compounds : list of residue IDs, optional
        """
        chunks = []
        new_chunk = []

        # Track selected clusters/compounds
        self.only_clusters = set(only_clusters or [])
        self.only_compounds = self._get_valid_residue_ids(only_compounds or [])

        # Split file into chunks (on each HEADER)
        for line in lines:
            if line.startswith("HEADER"):
                if not line.startswith("HEADER    ") and line.startswith("HEADER "):
                    line = line.replace("HEADER ", "HEADER    ")
                if new_chunk:
                    chunks.append(new_chunk)
                    new_chunk = []
            new_chunk.append(line)

        if new_chunk:
            chunks.append(new_chunk)

        final_structure_builder = StructureBuilder()
        final_structure_builder.init_structure(id, file)
        final_structure_builder.init_model(0)

        curr_chain_id = "a"
        curr_serial_number = 1
        global_line_counter = 0

        for cid, chunk_lines in enumerate(chunks):
            self.structure_builder = StructureBuilder()
            self.structure_builder.init_structure(f"Chunk{cid}")

            self.header = None
            self.conects = None
            self.trailer = None
            self.line_counter = 0

            self._parse(chunk_lines)

            self.structure_builder.set_header(self.header)
            self.structure_builder.set_conects(self.conects)

            structure = self.structure_builder.get_structure()

            if structure.header["head"] == "protein":
                final_structure_builder.set_header(structure.header)
                if self.structure_builder.atom:
                    curr_serial_number = self.structure_builder.atom.serial_number + 1

            elif structure.header["head"].startswith("crosscluster.") and self.only_clusters:
                cluster_id = int(structure.header["head"].split(".")[1])
                if self.only_clusters and cluster_id not in self.only_clusters:
                    continue

            remapped_conects = {}
            for model_id in structure.child_dict:
                if not final_structure_builder.structure.has_id(model_id):
                    final_structure_builder.init_model(model_id)

                # Set the current model to be the one with id `model_id`, which can be an existing
                # model or a model just initialized.
                final_structure_builder.model = final_structure_builder.structure[model_id]

                empty_chains = []
                for chain_id in list(structure[model_id].child_dict):

                    # If some chain id is empty, then replace it with the current default chain id.
                    if chain_id.strip() == "":
                        structure[model_id][chain_id].id = curr_chain_id
                        chain_id = curr_chain_id
                        curr_chain_id = chr(ord(curr_chain_id) + 1)

                    if not final_structure_builder.model.has_id(chain_id):
                        final_structure_builder.model.add(structure[model_id][chain_id])

                        # Set the current chain to be the one given by the id `chain_id`.
                        final_structure_builder.chain = final_structure_builder.model[chain_id]

                        # Only modify residues if the structure contains a cluster of probes.
                        if structure.header["head"].startswith("crosscluster."):
                            remove_residues = []
                            # Loop over each residue and set probes to be considered ligands.
                            for res in final_structure_builder.chain:
                                # Modify probes to be hetatm.
                                if res.resname in PROBES_SET:
                                    res.id = ("H_%s" % res.resname, res.id[1], res.id[2])
                                    res.cluster_id = int(structure.header["head"].split(".")[1])

                                # Update the line number where the residue starts.
                                res.at_line = global_line_counter + res.at_line

                                # Ignore compounds that are not in the set ``only_compounds``.
                                if not self._accept_residue(res):
                                    remove_residues.append(res.id)
                                else:
                                    for atm in res.get_atoms():
                                        remapped_conects[atm.serial_number] = curr_serial_number
                                        atm.serial_number = curr_serial_number
                                        curr_serial_number += 1

                            # Remove residues not defined in ``only_compounds``.
                            for res_id in remove_residues:
                                del final_structure_builder.chain[res_id]
                    else:
                        # Set the current chain to be the one given by the id `chain_id`.
                        final_structure_builder.chain = final_structure_builder.model[chain_id]

                        # Loop over each residue and set probes to be considered ligands.
                        for res in structure[model_id][chain_id]:
                            # Only modify residues if the structure contains a cluster of probes.
                            if structure.header["head"].startswith("crosscluster."):
                                # Modify probes to be hetatm.
                                if res.resname in PROBES_SET:
                                    res.id = ("H_%s" % res.resname, res.id[1], res.id[2])
                                    res.cluster_id = int(structure.header["head"].split(".")[1])

                            # Update the line number where the residue starts.
                            res.at_line = global_line_counter + res.at_line

                            try:
                                # Ignore compounds that are not in the set ``only_compounds``.
                                if not self._accept_residue(res):
                                    continue
                                else:
                                    if curr_serial_number is not None:
                                        for atm in res.get_atoms():
                                            remapped_conects[atm.serial_number] = curr_serial_number
                                            atm.serial_number = curr_serial_number
                                            curr_serial_number += 1

                                final_structure_builder.chain.add(res)
                            except PDBConstructionException as message:
                                self._handle_PDB_exception(message, res.at_line)

                    if len(final_structure_builder.chain) == 0:
                        empty_chains.append(final_structure_builder.chain.id)

                # Remove empty chains.
                for chain_id in empty_chains:
                    del final_structure_builder.model[chain_id]

            # Remap the conects based on the atoms' new numbering.
            if remapped_conects:
                new_conects = {}
                for atm1, atoms in self.structure_builder.conects.items():
                    # Ignore atoms that were removed due to ``only_clusters`` and ``only_compounds``.
                    if atm1 not in remapped_conects:
                        continue

                    new_atoms_list = []
                    for atm2 in atoms:
                        # Ignore atoms that were removed due to ``only_clusters`` and ``only_compounds``.
                        if atm2 not in remapped_conects:
                            continue
                        new_atoms_list.append(remapped_conects[atm2])
                    new_conects[remapped_conects[atm1]] = new_atoms_list

                self.structure_builder.set_conects(new_conects)

                # Update the conects in final_structure_builder with the
                # remapped conects from the current structure
                final_structure_builder.conects.update(self.structure_builder.conects)
            elif structure.header["head"] == "protein":
                # Update the conects in final_structure_builder with the
                # conects from the original protein structure
                final_structure_builder.conects.update(self.structure_builder.conects)

            global_line_counter += self.structure_builder.line_counter

        self.structure_builder = final_structure_builder

    def _get_valid_residue_ids(self, res_ids):
        """
        Validate compound residue ID filters.

        Parameters
        ----------
        res_ids : iterable
            Residue ID tuples.

        Returns
        -------
        set
            Set of valid 3-tuple or 4-tuple residue IDs.
        """
        valid_residue_ids = set()

        for res_id in res_ids:
            if not isinstance(res_id, tuple) or len(res_id) not in [3, 4]:
                logger.warning("The residue id '%s' is not valid. "
                               "A valid residue id consists of a tuple containing three (hetflag, resseq, icode) or "
                               "four (structure id, model id, chain id, (hetflag, resseq, icode)) elements." % str(res_id))
                continue
            if len(res_id) == 4:
                if not isinstance(res_id[3], tuple) or len(res_id[3]) != 3:
                    logger.warning("The residue id '%s' is not valid. "
                                   "A valid residue id consists of a tuple containing three (hetflag, resseq, icode) or "
                                   "four (structure id, model id, chain id, (hetflag, resseq, icode)) elements." % str(res_id))
                    continue
            valid_residue_ids.add(res_id)
        return valid_residue_ids

    def _accept_residue(self, res):
        """
        Check whether a residue passes the compound filter.

        Parameters
        ----------
        res : Residue
            A Biopython Residue object.

        Returns
        -------
        bool
            True if accepted.
        """
        if self.only_compounds:
            for res_id in self.only_compounds:
                if len(res_id) == 4:
                    if res_id[1:] == res.get_full_id()[1:]:
                        return True
                elif len(res_id) == 3:
                    if res_id == res.id:
                        return True
            return False
        return True
