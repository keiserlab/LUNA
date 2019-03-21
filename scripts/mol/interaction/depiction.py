from mol.wrappers.pymol import PymolWrapper, mybio_to_pymol_selection
from util.default_values import PYMOL_INTERACTION_COLOR, INTERACTION_SHORT_NAMES
from util.file import get_filename
from util.exceptions import PymolSessionNotInitialized


class PymolInteractionViewer:

    def __init__(self, show_cartoon=False, bg_color="white", pharm_color=None,
                 inter_color=PYMOL_INTERACTION_COLOR, pse_export_version="1.8"):

        self.show_cartoon = show_cartoon
        self.bg_color = bg_color
        self.pse_export_version = pse_export_version
        self.inter_color = inter_color
        self.wrapper = None

    def new_session(self, inter_tuples, output_file):
        self.start_session()
        self.set_view(inter_tuples)
        self.save_session(output_file)
        self.finish_session()

    def start_session(self):
        self.wrapper = PymolWrapper()
        self.wrapper.set("pse_export_version", self.pse_export_version)
        self.wrapper.set("transparency_mode", 3)
        self.wrapper.set("group_auto_mode", 2)
        self.wrapper.run_cmds([("bg_color", {"color": self.bg_color})])

    def set_view(self, inter_tuples):

        for t, (pdb_file, interactions) in enumerate(inter_tuples):
            main_grp = "OBJ_%s" % get_filename(pdb_file)

            pdb_obj = "%s.prot" % main_grp
            self.wrapper.load(pdb_file, pdb_obj)
            self.wrapper.extract([("%s.hets" % main_grp, "hetatm and %s" % pdb_obj)])

            self.wrapper.color_by_element([main_grp])
            self.wrapper.hide([("everything", main_grp)])

            if self.show_cartoon:
                self.wrapper.show([("cartoon", main_grp)])

            for i, inter in enumerate(interactions):
                if inter.type == "Proximal":
                    continue

                obj1_name = "obj%d_grp%s" % (t, hash(inter.atm_grp1))
                obj2_name = "obj%d_grp%s" % (t, hash(inter.atm_grp2))

                # Set the representation for each compound in the groups involved in the interaction.
                for compound in inter.atm_grp1.compounds.union(inter.atm_grp2.compounds):
                    comp_repr = "sphere" if compound.is_water() else "sticks"
                    comp_sel = mybio_to_pymol_selection(compound)
                    self.wrapper.show([(comp_repr, comp_sel)])
                    carb_color = "green" if compound.is_target() else "gray"
                    self.wrapper.color([(carb_color, comp_sel + " AND elem C")])

                self.wrapper.add_pseudoatom(obj1_name, {"vdw": 1, "pos": list(inter.atm_grp1.centroid)})
                self.wrapper.add_pseudoatom(obj2_name, {"vdw": 1, "pos": list(inter.atm_grp2.centroid)})

                # Interactions
                inter_name = "%s.inters.%s.obj%d_inter%d" % (main_grp, INTERACTION_SHORT_NAMES[inter.type], t, i)
                self.wrapper.distance(inter_name, obj1_name, obj2_name)
                self.wrapper.hide([("label", inter_name)])
                # Set styles to the interactions.
                self.wrapper.color([(self.inter_color.get_color(inter.type), inter_name)])

                # If a group object contains more than one atom, the centroid object will be displayed.
                if len(inter.atm_grp1.atoms) > 1:
                    # Add the centroids to the group "grps" and append them to the main group
                    self.wrapper.group("%s.grps" % main_grp, [obj1_name])
                    self._set_centroid_style(obj1_name)
                # Otherwise, just remove the centroid as it will not add any new information (the atom represented
                # by the centroid is already been displayed).
                else:
                    self.wrapper.delete([obj1_name])

                # If a group object contains more than one atom, the centroid object will be displayed.
                if len(inter.atm_grp2.atoms) > 1:
                    # Add the centroids to the group "grps" and append them to the main group
                    self.wrapper.group("%s.grps" % main_grp, [obj2_name])
                    self._set_centroid_style(obj2_name)
                # Otherwise, just remove the centroid as it will not add any new information (the atom represented
                # by the centroid is already been displayed).
                else:
                    self.wrapper.delete([obj2_name])

        self.wrapper.set("sphere_scale", "0.3", {"selection": "visible and resn hoh"})

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
