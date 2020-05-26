from luna.mol.entry import MolEntry
from luna.mol.wrappers.pymol import PymolWrapper, PymolSessionManager, mybio_to_pymol_selection
from luna.util.exceptions import PymolSessionNotInitialized


import logging

logger = logging.getLogger()


class ShellViewer(PymolSessionManager):

    def set_view(self, shell_tuples):

        if not isinstance(self.wrapper, PymolWrapper):
            raise PymolSessionNotInitialized("No session was initialized.")

        for target_entry, shells in shell_tuples:
            pdb_file = target_entry.pdb_file
            main_grp = target_entry.to_string(sep="-")

            mol_obj = target_entry.mol_obj if isinstance(target_entry, MolEntry) else None

            # Load PDB and extract hetatm.
            self.set_pdb_view(pdb_file, main_grp, mol_obj)

            for index, shell in enumerate(shells):
                sphere_obj = "%s.spheres.s%d_lvl%d_center%d" % (main_grp, index, shell.level, hash(shell.central_atm_grp))
                centroid_obj = "%s.grps.centroid" % (sphere_obj)

                self.wrapper.add_pseudoatom(centroid_obj, {"color": "white", "pos": list(shell.central_atm_grp.centroid)})
                self.wrapper.hide([("nonbonded", centroid_obj)])
                self.wrapper.show([("spheres", centroid_obj), ("nb_spheres", centroid_obj), ("dots", centroid_obj)])
                self.wrapper.set("dot_color", "red")
                self.wrapper.set("sphere_scale", shell.radius, {"selection": centroid_obj})
                self.wrapper.set("sphere_transparency", 0.85, {"selection": centroid_obj})
                self.wrapper.run_cmds([("center", {"selection": centroid_obj})])

                # TMP
                self.wrapper.label([(centroid_obj, '"%s"' % "+".join([a.name for a in sorted(shell.central_atm_grp.atoms)]))])

                # Add interactions and styles.
                self.set_interactions_view(shell.interactions, sphere_obj)

        self.set_last_details_to_view()

        self.wrapper.label([("visible and !name PS*", "name")])

        for target_entry, shells in shell_tuples:
            for index, shell in enumerate(shells):
                for c in shell.central_atm_grp.compounds:
                    self.wrapper.hide([("sticks", mybio_to_pymol_selection(c))])
                    self.wrapper.show([("lines", mybio_to_pymol_selection(c))])

                    # Show residue label if required.
                    if self.show_comp_labels:
                        if c.is_residue():
                            self.wrapper.label([("%s AND name CA" % mybio_to_pymol_selection(c), '"%s-%s" % (resn, resi)')])
                        else:
                            any_atm = next(c.get_atoms())
                            atm_sel = mybio_to_pymol_selection(any_atm)
                            self.wrapper.label([(atm_sel, '"%s-%s" % (resn, resi)')])
