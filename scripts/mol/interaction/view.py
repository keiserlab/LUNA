from mol.entry import MolEntry
from mol.wrappers.pymol import PymolWrapper, PymolSessionManager
from util.exceptions import PymolSessionNotInitialized
from util.file import get_filename


class InteractionViewer(PymolSessionManager):

    def set_view(self, inter_tuples):

        if not isinstance(self.wrapper, PymolWrapper):
            raise PymolSessionNotInitialized("No session was initialized.")

        for uid, (target_entry, interactions_mngr) in enumerate(inter_tuples):
            pdb_file = target_entry.pdb_file
            main_grp = "OBJ%d_%s" % (uid, get_filename(pdb_file))

            mol_obj = target_entry.mol_obj if isinstance(target_entry, MolEntry) else None

            # Load PDB and extract hetatm.
            self.set_pdb_view(pdb_file, main_grp, mol_obj)

            # Add interactions and styles.
            self.set_interactions_view(interactions_mngr.interactions, main_grp, uid)

        self.set_last_details_to_view()
