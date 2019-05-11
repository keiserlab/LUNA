from mol.wrappers.pymol import PymolWrapper, mybio_to_pymol_selection
from util.default_values import PYMOL_INTERACTION_COLOR, INTERACTION_SHORT_NAMES
from util.file import get_filename
from util.exceptions import PymolSessionNotInitialized


class PymolInteractionViewer:

    def __init__(self, show_cartoon=False, bg_color="white", pharm_color=None,
                 add_directional_arrows=True, inter_color=PYMOL_INTERACTION_COLOR, pse_export_version="1.8"):

        self.show_cartoon = show_cartoon
        self.bg_color = bg_color
        self.pse_export_version = pse_export_version
        self.inter_color = inter_color
        self.wrapper = None
        self.add_directional_arrows = add_directional_arrows

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
        self.wrapper.set("internal_gui_width", 350)

    def set_view(self, inter_tuples):

        if not isinstance(self.wrapper, PymolWrapper):
            raise PymolSessionNotInitialized("No session was initialized.")

        directional_inters = ["Hydrogen bond", "Water-bridged hydrogen bond", "Weak hydrogen bond", "Halogen bond", "Halogen-pi",
                              "Chalcogen bond", "Chalcogen-pi", "Parallel multipolar", "Antiparallel multipolar", "Tilted multipolar",
                              "Orthogonal multipolar", "Multipolar"]

        nucleophile_inters = ["Orthogonal multipolar", "Parallel multipolar", "Antiparallel multipolar", "Tilted multipolar", "Multipolar",
                              "Cation-nucleophile", "Unfavorable anion-nucleophile", "Unfavorable nucleophile-nucleophile"]
        electrophile_inters = ["Orthogonal multipolar", "Parallel multipolar", "Antiparallel multipolar", "Tilted multipolar", "Multipolar",
                               "Anion-electrophile", "Unfavorable cation-electrophile", "Unfavorable electrophile-electrophile"]

        unfavorable_inters = ["Repulsive", "Unfavorable anion-nucleophile", "Unfavorable cation-electrophile",
                              "Unfavorable nucleophile-nucleophile", "Unfavorable electrophile-electrophile"]

        pdb_files_read = set()

        for t, (pdb_file, interactions) in enumerate(inter_tuples):
            main_grp = "OBJ%d_%s" % (t, get_filename(pdb_file))
            pdb_obj = "%s.prot" % main_grp

            if pdb_file not in pdb_files_read:
                self.wrapper.load(pdb_file, pdb_obj)
                self.wrapper.extract([("%s.hets" % main_grp, "hetatm and %s" % pdb_obj)])

                self.wrapper.color_by_element([main_grp])
                self.wrapper.hide([("everything", main_grp)])

                if self.show_cartoon:
                    self.wrapper.show([("cartoon", main_grp)])
            else:
                pdb_files_read.add(pdb_file)

            for i, inter in enumerate(interactions):
                if inter.type == "Proximal":
                    continue

                obj1_name = "obj%d_grp%s" % (t, hash(inter.atm_grp1))
                centroid_obj1 = inter.atm_grp1.centroid
                # Define the centroid in a nucleophile with two atoms as the position of its more electronegative atom.
                # Remember that the position in the interaction object matters. We have defined that the first group is always
                # the nucleophile for both dipole-dipole and ion-dipole interactions.
                if inter.type in nucleophile_inters and len(inter.atm_grp1.atoms) == 2:
                    dipole_atm = inter.atm_grp1.atoms[0] if (inter.atm_grp1.atoms[0].electronegativity >
                                                             inter.atm_grp1.atoms[1].electronegativity) else inter.atm_grp1.atoms[1]
                    obj1_name += "_%s" % hash(dipole_atm.name)
                    centroid_obj1 = dipole_atm.coord
                # For unfavorable multipolar interactions, it may happen that the first atom group is an electrophile as well.
                elif inter.type == "Unfavorable electrophile-electrophile" and len(inter.atm_grp1.atoms) == 2:
                    dipole_atm = inter.atm_grp1.atoms[0] if (inter.atm_grp1.atoms[0].electronegativity <
                                                             inter.atm_grp1.atoms[1].electronegativity) else inter.atm_grp1.atoms[1]
                    obj1_name += "_%s" % hash(dipole_atm.name)
                    centroid_obj1 = dipole_atm.coord

                obj2_name = "obj%d_grp%s" % (t, hash(inter.atm_grp2))
                centroid_obj2 = inter.atm_grp2.centroid
                # Define the centroid in an electrophile with two atoms as the position of its less electronegative atom.
                # Remember that the position in the interaction object matters. We have defined that the second group is always
                # the electrophile for both dipole-dipole and ion-dipole interactions.
                if inter.type in electrophile_inters and len(inter.atm_grp2.atoms) == 2:
                    dipole_atm = inter.atm_grp2.atoms[0] if (inter.atm_grp2.atoms[0].electronegativity <
                                                             inter.atm_grp2.atoms[1].electronegativity) else inter.atm_grp2.atoms[1]
                    obj2_name += "_%s" % hash(dipole_atm.name)
                    centroid_obj2 = dipole_atm.coord
                # For unfavorable multipolar interactions, it may happen that the second atom group is an nucleophile as well.
                elif inter.type == "Unfavorable nucleophile-nucleophile" and len(inter.atm_grp2.atoms) == 2:
                    dipole_atm = inter.atm_grp2.atoms[0] if (inter.atm_grp2.atoms[0].electronegativity >
                                                             inter.atm_grp2.atoms[1].electronegativity) else inter.atm_grp2.atoms[1]
                    obj2_name += "_%s" % hash(dipole_atm.name)
                    centroid_obj2 = dipole_atm.coord

                # Add pseudoatoms
                if not self.wrapper.obj_exists(obj1_name):
                    self.wrapper.add_pseudoatom(obj1_name, {"vdw": 1, "pos": list(centroid_obj1)})

                if not self.wrapper.obj_exists(obj2_name):
                    self.wrapper.add_pseudoatom(obj2_name, {"vdw": 1, "pos": list(centroid_obj2)})

                # Set the representation for each compound in the groups involved in the interaction.
                for compound in inter.atm_grp1.compounds.union(inter.atm_grp2.compounds):
                    comp_repr = "sphere" if compound.is_water() else "sticks"
                    comp_sel = mybio_to_pymol_selection(compound)
                    self.wrapper.show([(comp_repr, comp_sel)])
                    carb_color = "green" if compound.is_target() else "gray"
                    self.wrapper.color([(carb_color, comp_sel + " AND elem C")])

                # Check if the interaction involves the same compound: intramolecular interactions.
                inter_grp = "intra" if inter.is_intramol_interaction() else "inter"

                inter_name = "%s.all_inters.%s.%s.obj%d_inter%d.line" % (main_grp, inter_grp, INTERACTION_SHORT_NAMES[inter.type], t, i)
                self.wrapper.distance(inter_name, obj1_name, obj2_name)
                self.wrapper.hide([("label", inter_name)])

                # Set styles to the interactions.
                self.wrapper.color([(self.inter_color.get_color(inter.type), inter_name)])

                if self.add_directional_arrows:
                    # Add arrows over the interaction lines to represent directional interactions
                    if inter.type in directional_inters:
                        arrow_name = "%s.all_inters.%s.%s.obj%d_inter%d.arrow" % (main_grp, inter_grp,
                                                                                  INTERACTION_SHORT_NAMES[inter.type], t, i)
                        arrow_opts = {"radius": 0.03, "gap": 0.9, "hlength": 0.5, "hradius": 0.2,
                                      "color": self.inter_color.get_color(inter.type)}
                        self.wrapper.arrow(arrow_name, obj1_name, obj2_name, arrow_opts)

                    if inter.type in unfavorable_inters:
                        arrow_name1 = "%s.all_inters.%s.%s.obj%d_inter%d.arrow1" % (main_grp, inter_grp,
                                                                                    INTERACTION_SHORT_NAMES[inter.type], t, i)
                        arrow_name2 = "%s.all_inters.%s.%s.obj%d_inter%d.arrow2" % (main_grp, inter_grp,
                                                                                    INTERACTION_SHORT_NAMES[inter.type], t, i)
                        square_name = "%s.all_inters.%s.%s.obj%d_inter%d.block" % (main_grp, inter_grp,
                                                                                   INTERACTION_SHORT_NAMES[inter.type], t, i)

                        arrow_opts = {"radius": 0.03, "gap": 0.9, "hlength": 0.5, "hradius": 0.2,
                                      "color": self.inter_color.get_color(inter.type)}
                        square_opts = {"radius": 0.3, "gap": 1.5, "hlength": 0, "hradius": 0,
                                       "color": self.inter_color.get_color(inter.type)}

                        # Two arrows in different directions
                        self.wrapper.arrow(arrow_name1, obj1_name, obj2_name, arrow_opts)
                        self.wrapper.arrow(arrow_name2, obj2_name, obj1_name, arrow_opts)
                        # Add a square-like object
                        self.wrapper.arrow(square_name, obj1_name, obj2_name, square_opts)

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
