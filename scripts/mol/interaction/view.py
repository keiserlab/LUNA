from mol.wrappers.pymol import PymolWrapper, PymolSessionManager
from util.exceptions import PymolSessionNotInitialized
from util.file import get_filename


class InteractionViewer(PymolSessionManager):

    def set_view(self, inter_tuples):

        if not isinstance(self.wrapper, PymolWrapper):
            raise PymolSessionNotInitialized("No session was initialized.")

        pdb_files_read = set()

        for uid, (target_entry, interactions_mngr) in enumerate(inter_tuples):
            pdb_file = target_entry.pdb_file
            main_grp = "OBJ%d_%s" % (uid, get_filename(pdb_file))

            if pdb_file not in pdb_files_read:
                # Load PDB and extract hetatm.
                self.set_pdb_view(pdb_file, main_grp)
                pdb_files_read.add(pdb_file)

            # Add interactions and styles.
            self.set_interactions_view(interactions_mngr.interactions, main_grp, uid)

        self.set_last_details_to_view()
