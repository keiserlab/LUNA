from mol.wrappers.pymol import PymolWrapper, PymolSessionManager, mybio_to_pymol_selection
from util.exceptions import PymolSessionNotInitialized
from util.file import get_filename


import logging

logger = logging.getLogger()


class ShellViewer(PymolSessionManager):

    def set_view(self, shell_tuples):

        if not isinstance(self.wrapper, PymolWrapper):
            raise PymolSessionNotInitialized("No session was initialized.")

        pdb_files_read = set()

        for uid, (pdb_file, shell) in enumerate(shell_tuples):
            main_grp = "OBJ_%s" % get_filename(pdb_file)

            if pdb_file not in pdb_files_read:
                # Load PDB and extract hetatm.
                self.set_pdb_view(pdb_file, main_grp)
                pdb_files_read.add(pdb_file)

            sphere_obj = "%s.spheres.s%d_lvl%d_center%d" % (main_grp, uid, shell.level, hash(shell.central_atm_grp))

            centroid_obj = "%s.grps.s%d_center" % (sphere_obj, uid)
            self.wrapper.add_pseudoatom(centroid_obj, {"color": "white", "pos": list(shell.central_atm_grp.centroid)})
            self.wrapper.hide([("nonbonded", centroid_obj)])
            self.wrapper.show([("sphere", centroid_obj), ("nb_spheres", centroid_obj), ("dots", centroid_obj)])
            self.wrapper.set("dot_color", "red")
            self.wrapper.set("sphere_scale", shell.radius, {"selection": centroid_obj})
            self.wrapper.set("sphere_transparency", 0.85, {"selection": centroid_obj})
            self.wrapper.run_cmds([("center", {"selection": centroid_obj})])

            # Compound view (residue, ligand, etc)
            for compound in shell.central_atm_grp.compounds:
                comp_repr = "sphere" if compound.is_water() else "sticks"
                carb_color = "green" if compound.is_target() else "gray"
                self.wrapper.show([(comp_repr, mybio_to_pymol_selection(compound))])
                self.wrapper.color([(carb_color, mybio_to_pymol_selection(compound) + " AND elem C")])

            # Add interactions and styles.
            self.set_interactions_view(shell.interactions, main_grp, uid)

        self.set_last_details_to_view()
