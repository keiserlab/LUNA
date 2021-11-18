from luna.mol.entry import MolEntry
from luna.mol.wrappers.pymol import PymolWrapper, PymolSessionManager, mybio_to_pymol_selection
from luna.util.exceptions import PymolSessionNotInitialized
from luna.MyBio.util import entity_to_string


import logging

logger = logging.getLogger()


class ShellViewer(PymolSessionManager):

    def set_view(self, shell_tuples):

        if not isinstance(self.wrapper, PymolWrapper):
            raise PymolSessionNotInitialized("No session was initialized.")

        for target_entry, shells, pdb_path in shell_tuples:
            pdb_file = "%s/%s.pdb" % (pdb_path, target_entry.pdb_id)
            main_grp = target_entry.to_string(sep="-")

            mol_block = (entity_to_string(target_entry.get_biopython_structure())
                         if isinstance(target_entry, MolEntry) else None)

            # Load PDB and extract hetatm.
            self.set_pdb_view(pdb_file, main_grp, mol_block)

            residue_selections = set()
            for index, shell in enumerate(shells):
                sphere_obj = "%s.spheres.s%d_lvl%d_center%d" % (main_grp, index, shell.level, hash(shell.central_atm_grp))
                centroid_obj = "%s.centroid" % (sphere_obj)

                self.wrapper.add_pseudoatom(centroid_obj, {"color": "white", "pos": list(shell.central_atm_grp.centroid)})
                self.wrapper.hide([("nonbonded", centroid_obj)])
                self.wrapper.show([("spheres", centroid_obj), ("nb_spheres", centroid_obj), ("dots", centroid_obj)])
                self.wrapper.set("dot_color", "red")
                self.wrapper.set("sphere_scale", shell.radius, {"selection": centroid_obj})
                self.wrapper.set("sphere_transparency", 0.85, {"selection": centroid_obj})
                self.wrapper.run_cmds([("center", {"selection": centroid_obj})])

                self.wrapper.label([(centroid_obj, '"%s"' % "+".join([a.name for a in sorted(shell.central_atm_grp.atoms)]))])

                # Add interactions and styles.
                interacting_residue_sels = self.set_interactions_view(shell.interactions, main_grp, sphere_obj)

                residue_selections.update(interacting_residue_sels)

            self.wrapper.select(name="%s.inter_residues" % main_grp, selection=" or ".join(residue_selections))

        self.set_last_details_to_view()
