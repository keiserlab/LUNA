from mol.entry import MolEntry
from mol.wrappers.pymol import PymolWrapper, PymolSessionManager, mybio_to_pymol_selection
from util.exceptions import PymolSessionNotInitialized


class InteractionViewer(PymolSessionManager):

    def set_view(self, inter_tuples):

        if not isinstance(self.wrapper, PymolWrapper):
            raise PymolSessionNotInitialized("No session was initialized.")

        for (target_entry, interactions_mngr) in inter_tuples:
            pdb_file = target_entry.pdb_file
            main_grp = target_entry.to_string(sep="-")

            mol_obj = target_entry.mol_obj if isinstance(target_entry, MolEntry) else None

            # Load PDB and extract hetatm.
            self.set_pdb_view(pdb_file, main_grp, mol_obj)

            # Add interactions and styles.
            self.set_interactions_view(interactions_mngr.interactions, main_grp)

        self.set_last_details_to_view()

        # for (target_entry, interactions_mngr) in inter_tuples:
        #     pdb_file = target_entry.pdb_file
        #     main_grp = target_entry.to_string(sep="-")

        #     self.wrapper.color([("white", "%s and visible and !name PS*" % main_grp)])

        #     interacting_atms = set()
        #     for inter in interactions_mngr.interactions:
        #         if inter.type != "Hydrophobic":
        #             continue

        #         for atm in inter.src_grp.atoms + inter.trgt_grp.atoms:
        #             if atm not in interacting_atms:
        #                 self.wrapper.color([("pink", "%s and %s" % (main_grp, mybio_to_pymol_selection(atm)))])

        #         for atm in inter.src_interacting_atms + inter.trgt_interacting_atms:
        #             self.wrapper.color([("hotpink", "%s and %s" % (main_grp, mybio_to_pymol_selection(atm)))])
        #             interacting_atms.add(atm)
