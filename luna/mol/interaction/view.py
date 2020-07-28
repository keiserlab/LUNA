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
            main_grp = target_entry.to_string(sep="-").replace("'", "-")

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

                # It will display all hydrophobic groups if a AtomGroupsManager is available.
                atm_grp_mngr = interactions[0].src_grp.manager
                if atm_grp_mngr is not None:
                    interacting_atms = set()
                    for atm_grp in atm_grp_mngr.filter_by_types(["Hydrophobe", "Hydrophobic"], must_contain_all=False):
                        inters = [i for i in atm_grp.interactions if i.type == "Hydrophobic"]

                        for comp in atm_grp.compounds:
                            comp_sel = mybio_to_pymol_selection(comp)
                            self.wrapper.show([("sticks", comp_sel)])

                        if len(inters) > 0:
                            for inter in inters:
                                for atm in inter.src_grp.atoms + inter.trgt_grp.atoms:
                                    if atm not in interacting_atms:
                                        self.wrapper.color([("pink", "%s and %s" % (main_grp, mybio_to_pymol_selection(atm)))])

                                for atm in inter.src_interacting_atms + inter.trgt_interacting_atms:
                                    self.wrapper.color([("hotpink", "%s and %s" % (main_grp, mybio_to_pymol_selection(atm)))])
                                    interacting_atms.add(atm)

                        else:
                            for atm in atm_grp.atoms:
                                self.wrapper.color([("wheat", "%s and %s" % (main_grp, mybio_to_pymol_selection(atm)))])

                # Otherwise, it will display only hydrophobic groups comprising the interacting groups.
                else:
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
