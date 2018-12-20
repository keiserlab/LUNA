from mol.wrappers.pymol import (PymolWrapper, PymolColorMap, mybio_to_pymol_selection)
from util.exceptions import PymolSessionNotInitialized

import logging
logger = logging.getLogger(__name__)


DEFAULT_INTERACTIONS_COLOR = {
    "Proximal": "gray60",
    "Hydrogen bond": "tv_blue",
    "Attractive": "forest",
    "Salt bridge": "palegreen",
    "Cation-pi": "salmon",
    "Edge-to-face pi-stacking": "tv_red",
    "Face-to-face pi-stacking": "tv_red",
    "Parallel-displaced pi-stacking": "tv_red",
    "Hydrophobic": "orange",
    "Halogen bond": "aquamarine",
    "Repulsive": "violetpurple",
    "Water-bridged hydrogen bond": "lightblue",
    "Pi-stacking": "tv_red"
}


class PymolShellViewer:

    def __init__(self, input_file, show_cartoon=False, bg_color="white", pharm_color=None,
                 inter_color=None, pse_export_version="1.8"):

        self.input_file = input_file
        self.show_cartoon = show_cartoon
        self.bg_color = bg_color
        self.pse_export_version = pse_export_version
        self.wrapper = None

        self.inter_color = inter_color or PymolColorMap(DEFAULT_INTERACTIONS_COLOR, "white")

    def new_session(self, shells, output_file):
        self.start_session()
        self.set_view(shells)
        self.save_session(output_file)
        self.finish_session()

    def start_session(self):
        self.wrapper = PymolWrapper()

        self.wrapper.set("pse_export_version", self.pse_export_version)
        self.wrapper.set("transparency_mode", 3)
        self.wrapper.run_cmds([("bg_color", {"color": self.bg_color})])

        self.wrapper.load(self.input_file)
        self.wrapper.color_by_element(["all"])
        self.wrapper.hide_all()

        if self.show_cartoon:
            self.wrapper.show([("cartoon", "all")])

    def set_view(self, shells):
        if not isinstance(self.wrapper, PymolWrapper):
            raise PymolSessionNotInitialized("No session was initialized.")

        for (s, shell) in enumerate(shells):
            centroid_name = "sphere_%d" % s

            self.wrapper.add_pseudoatom(centroid_name, {"color": "white", "pos": list(shell.central_atm_grp.centroid)})

            self.wrapper.hide([("nonbonded", centroid_name)])
            self.wrapper.show([("sphere", centroid_name), ("nb_spheres", centroid_name), ("dots", centroid_name)])
            self.wrapper.set("dot_color", "red")

            # Compound view (residue, ligand, etc)
            comp_repr = "sphere" if shell.central_atm_grp.compound.is_water() else "sticks"
            carb_color = "green" if shell.central_atm_grp.compound.is_target() else "gray"

            self.wrapper.show([(comp_repr, mybio_to_pymol_selection(shell.central_atm_grp.compound))])
            self.wrapper.color([(carb_color, mybio_to_pymol_selection(shell.central_atm_grp.compound) + " AND elem C")])

            # if shell.central_atm_grp.compound.is_water():
                # self.wrapper.set("sphere_scale", "0.3", {"selection": centroid_name})

            self.wrapper.set("sphere_scale", shell.radius, {"selection": centroid_name})
            self.wrapper.set("sphere_transparency", 0.7, {"selection": centroid_name})
            self.wrapper.run_cmds([("center", {"selection": centroid_name})])

            for (i, inter) in enumerate(shell.interactions):
                if inter.type == "Proximal":
                    continue

                obj1_name = "Group_%s" % hash(inter.atm_grp1)
                obj2_name = "Group_%s" % hash(inter.atm_grp2)

                comp1_sel = mybio_to_pymol_selection(inter.atm_grp1.compound)
                comp2_sel = mybio_to_pymol_selection(inter.atm_grp2.compound)

                comp1_repr = "sphere" if inter.atm_grp1.compound.is_water() else "sticks"
                comp2_repr = "sphere" if inter.atm_grp2.compound.is_water() else "sticks"
                self.wrapper.show([(comp1_repr, comp1_sel), (comp2_repr, comp2_sel)])

                # if inter.atm_grp1.compound.is_water():
                #     self.wrapper.set("sphere_scale", "0.3", {"selection": comp1_sel})
                # if inter.atm_grp2.compound.is_water():
                #     self.wrapper.set("sphere_scale", "0.3", {"selection": comp2_sel})

                carb1_color = "green" if inter.atm_grp1.compound.is_target() else "gray"
                carb2_color = "green" if inter.atm_grp2.compound.is_target() else "gray"
                self.wrapper.color([(carb1_color, comp1_sel + " AND elem C")])
                self.wrapper.color([(carb2_color, comp2_sel + " AND elem C")])

                self.wrapper.add_pseudoatom(obj1_name, {"vdw": 1, "pos": list(inter.atm_grp1.centroid)})
                self.wrapper.add_pseudoatom(obj2_name, {"vdw": 1, "pos": list(inter.atm_grp2.centroid)})

                self.wrapper.hide([("nonbonded", obj1_name), ("nonbonded", obj2_name)])
                self.wrapper.show([("sphere", obj1_name), ("sphere", obj2_name)])

                self.wrapper.set("sphere_scale", 0.4, {"selection": obj1_name})
                self.wrapper.set("sphere_scale", 0.4, {"selection": obj2_name})

                self.wrapper.distance("inter_%s_%d" % (s, i), obj1_name, obj2_name)
                self.wrapper.color([(self.inter_color.get_color(inter.type), "inter_%s_%d" % (s, i))])

        self.wrapper.set("sphere_scale", "0.3", {"selection": "visible and resn hoh"})

    def save_session(self, output_file):
        if not isinstance(self.wrapper, PymolWrapper):
            raise PymolSessionNotInitialized("No session was initialized.")

        self.wrapper.save_session(output_file)

    def finish_session(self):
        if not isinstance(self.wrapper, PymolWrapper):
            raise PymolSessionNotInitialized("No session was initialized.")

        self.wrapper.reinitialize()
        self.wrapper = None
