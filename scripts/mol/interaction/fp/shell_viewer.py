from mol.wrappers.pymol import PymolWrapper, mybio_to_pymol_selection
from util.exceptions import PymolSessionNotInitialized
from util.default_values import PYMOL_INTERACTION_COLOR, INTERACTION_SHORT_NAMES
from util.file import get_filename


import logging

logger = logging.getLogger()


class PymolShellViewer:

    def __init__(self, show_cartoon=False, bg_color="white", pharm_color=None,
                 inter_color=PYMOL_INTERACTION_COLOR, pse_export_version="1.8"):

        self.show_cartoon = show_cartoon
        self.bg_color = bg_color
        self.pse_export_version = pse_export_version
        self.inter_color = inter_color
        self.wrapper = None

    def new_session(self, shell_tuples, output_file):
        self.start_session()
        self.set_view(shell_tuples)
        self.save_session(output_file)
        self.finish_session()

    def start_session(self):
        self.wrapper = PymolWrapper()
        self.wrapper.set("pse_export_version", self.pse_export_version)
        self.wrapper.set("transparency_mode", 3)
        self.wrapper.set("group_auto_mode", 2)
        self.wrapper.run_cmds([("bg_color", {"color": self.bg_color})])
        self.wrapper.set("internal_gui_width", 350)

    def set_view(self, shell_tuples):

        if not isinstance(self.wrapper, PymolWrapper):
            raise PymolSessionNotInitialized("No session was initialized.")

        directional_inters = ["Hydrogen bond", "Water-bridged hydrogen bond", "Weak hydrogen bond", "Halogen bond", "Halogen-pi",
                              "Chalcogen bond", "Chalcogen-pi"]

        pdb_files_read = set()

        for t, (pdb_file, shell) in enumerate(shell_tuples):
            main_grp = "OBJ_%s" % get_filename(pdb_file)
            pdb_obj = "%s.prot" % main_grp

            if pdb_file not in pdb_files_read:
                self.wrapper.load(pdb_file, pdb_obj)
                self.wrapper.extract([("%s.hets" % main_grp, "hetatm and %s" % pdb_obj)])

                self.wrapper.color_by_element([main_grp])
                self.wrapper.hide([("everything", main_grp)])
                if self.show_cartoon:
                    self.wrapper.show([("cartoon", main_grp)])

                pdb_files_read.add(pdb_file)

            sphere_obj = "%s.spheres.s%d_lvl%d_center%d" % (main_grp, t, shell.level, hash(shell.central_atm_grp))

            centroid_obj = "%s.grps.s%d_center" % (sphere_obj, t)
            self.wrapper.add_pseudoatom(centroid_obj, {"color": "white", "pos": list(shell.central_atm_grp.centroid)})
            self.wrapper.hide([("nonbonded", centroid_obj)])
            self.wrapper.show([("sphere", centroid_obj), ("nb_spheres", centroid_obj), ("dots", centroid_obj)])
            self.wrapper.set("dot_color", "red")
            self.wrapper.set("sphere_scale", shell.radius, {"selection": centroid_obj})
            self.wrapper.set("sphere_transparency", 0.7, {"selection": centroid_obj})
            self.wrapper.run_cmds([("center", {"selection": centroid_obj})])

            # Compound view (residue, ligand, etc)
            for compound in shell.central_atm_grp.compounds:
                comp_repr = "sphere" if compound.is_water() else "sticks"
                carb_color = "green" if compound.is_target() else "gray"
                self.wrapper.show([(comp_repr, mybio_to_pymol_selection(compound))])
                self.wrapper.color([(carb_color, mybio_to_pymol_selection(compound) + " AND elem C")])

            for i, inter in enumerate(shell.interactions):
                if inter.type == "Proximal":
                    continue

                obj1_name = "obj%d_grp%s" % (t, hash(inter.src_grp))
                obj2_name = "obj%d_grp%s" % (t, hash(inter.trgt_grp))

                # Pseudoatoms
                if not self.wrapper.obj_exists(obj1_name):
                    self.wrapper.add_pseudoatom(obj1_name, {"vdw": 1, "pos": list(inter.src_grp.centroid)})
                if not self.wrapper.obj_exists(obj2_name):
                    self.wrapper.add_pseudoatom(obj2_name, {"vdw": 1, "pos": list(inter.trgt_grp.centroid)})

                # Set the representation for each compound in the groups involved in the interaction.
                for compound in inter.src_grp.compounds.union(inter.trgt_grp.compounds):
                    comp_repr = "sphere" if compound.is_water() else "sticks"
                    comp_sel = mybio_to_pymol_selection(compound)
                    self.wrapper.show([(comp_repr, comp_sel)])

                    carb_color = "green" if compound.is_target() else "gray"
                    self.wrapper.color([(carb_color, comp_sel + " AND elem C")])

                # Check if the interaction involves the same compound: intramolecular interactions.
                inter_grp = "intra" if inter.is_intramol_interaction() else "inter"
                inter_name = "%s.all_inters.%s.%s.obj%d_inter%d.line%d" % (sphere_obj, inter_grp,
                                                                           INTERACTION_SHORT_NAMES[inter.type], t, hash(inter), hash(inter))
                self.wrapper.distance(inter_name, obj1_name, obj2_name)
                self.wrapper.hide([("label", inter_name)])

                # Set styles to the interactions.
                self.wrapper.color([(self.inter_color.get_color(inter.type), inter_name)])

                # Add arrows over the interaction lines to represent directional interactions
                if inter.type in directional_inters:
                    arrow_name = "%s.all_inters.%s.%s.obj%d_inter%d.arrow%d" % (sphere_obj, inter_grp, INTERACTION_SHORT_NAMES[inter.type],
                                                                                t, hash(inter), hash(inter))
                    arrow_opts = {"radius": 0.07, "gap": 0.9, "hlength": 0.4, "hlength": 0.4,
                                  "hlength": 0.4, "color": self.inter_color.get_color(inter.type)}
                    self.wrapper.arrow(arrow_name, obj1_name, obj2_name, arrow_opts)

                # If a group object contains more than one atom, the centroid object will be displayed.
                if len(inter.src_grp.atoms) > 1:
                    # Add the centroids to the group "grps" and append them to the main group
                    self.wrapper.group("%s.grps" % sphere_obj, [obj1_name])
                    self._set_centroid_style(obj1_name)
                # Otherwise, just remove the centroid as it will not add any new information (the atom represented
                # by the centroid is already been displayed).
                else:
                    self.wrapper.delete([obj1_name])

                # If a group object contains more than one atom, the centroid object will be displayed.
                if len(inter.trgt_grp.atoms) > 1:
                    # Add the centroids to the group "grps" and append them to the main group
                    self.wrapper.group("%s.grps" % sphere_obj, [obj2_name])
                    self._set_centroid_style(obj2_name)
                # Otherwise, just remove the centroid as it will not add any new information (the atom represented
                # by the centroid is already been displayed).
                else:
                    self.wrapper.delete([obj2_name])

        self.wrapper.set("sphere_scale", "0.3", {"selection": "visible and resn hoh"})
        self.wrapper.hide([("everything", "elem H")])

    def _set_centroid_style(self, centroid):
        # Set styles to the centroid.
        self.wrapper.hide([("nonbonded", centroid)])
        self.wrapper.show([("sphere", centroid)])
        self.wrapper.set("sphere_scale", 0.2, {"selection": centroid})

    def save_session(self, output_file):
        if not isinstance(self.wrapper, PymolWrapper):
            raise PymolSessionNotInitialized("No session was initialized.")

        self.wrapper.save_session(output_file)

    def finish_session(self):
        if not isinstance(self.wrapper, PymolWrapper):
            raise PymolSessionNotInitialized("No session was initialized.")

        self.wrapper.reinitialize()
        self.wrapper = None
