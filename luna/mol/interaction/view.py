from luna.mol.entry import MolEntry
from luna.mol.interaction.calc import InteractionsManager
from luna.mol.wrappers.pymol import PymolWrapper, PymolSessionManager, mybio_to_pymol_selection
from luna.util.exceptions import PymolSessionNotInitialized
from luna.MyBio.util import entity_to_string


class InteractionViewer(PymolSessionManager):

    def set_view(self, inter_tuples):

        if not isinstance(self.wrapper, PymolWrapper):
            raise PymolSessionNotInitialized("No session was initialized.")

        for target_entry, inter_data in inter_tuples:
            pdb_file = target_entry.pdb_file
            main_grp = target_entry.to_string(sep="-")

            mol_block = (entity_to_string(target_entry.get_biopython_structure())
                         if isinstance(target_entry, MolEntry) else None)

            # Load PDB and extract hetatm.
            self.set_pdb_view(pdb_file, main_grp, mol_block)

            if isinstance(inter_data, InteractionsManager):
                interactions = inter_data.interactions
            else:
                interactions = inter_data

            # Add interactions and styles.
            self.set_interactions_view(interactions, main_grp)

            if self.show_hydrop_surface:
                self.wrapper.color([("white", "%s and !name PS*" % main_grp)])

                interacting_atms = set()
                for inter in interactions:
                    if inter.type != "Hydrophobic":
                        continue

                    for atm in inter.src_grp.atoms + inter.trgt_grp.atoms:
                        if atm not in interacting_atms:
                            self.wrapper.color([("pink", "%s and %s" % (main_grp, mybio_to_pymol_selection(atm)))])

                    for atm in inter.src_interacting_atms + inter.trgt_interacting_atms:
                        self.wrapper.color([("hotpink", "%s and %s" % (main_grp, mybio_to_pymol_selection(atm)))])
                        interacting_atms.add(atm)

        self.set_last_details_to_view()
