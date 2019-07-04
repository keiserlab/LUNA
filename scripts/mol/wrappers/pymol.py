from pymol import cmd
from pymol import util

from mol.wrappers.cgo_arrow import cgo_arrow
from util.exceptions import PymolSessionNotInitialized
from util.default_values import PYMOL_INTERACTION_COLOR, INTERACTION_SHORT_NAMES
from util.file import get_filename, get_file_format

import logging

logger = logging.getLogger()

NUCLEOPHILE_INTERS = ["Orthogonal multipolar", "Parallel multipolar", "Antiparallel multipolar", "Tilted multipolar", "Multipolar",
                      "Cation-nucleophile", "Unfavorable anion-nucleophile", "Unfavorable nucleophile-nucleophile"]

ELECTROPHILE_INTERS = ["Orthogonal multipolar", "Parallel multipolar", "Antiparallel multipolar", "Tilted multipolar", "Multipolar",
                       "Anion-electrophile", "Unfavorable cation-electrophile", "Unfavorable electrophile-electrophile"]

UNFAVORABLE_INTERS = ["Repulsive", "Unfavorable anion-nucleophile", "Unfavorable cation-electrophile",
                      "Unfavorable nucleophile-nucleophile", "Unfavorable electrophile-electrophile"]


class PymolWrapper:

    def load(self, input_file, obj_name=None):
        self.input_file = input_file
        if not obj_name:
            obj_name = get_filename(input_file)
        cmd.load(input_file, obj_name)

    def show(self, tuples):
        for representation, selection in tuples:
            cmd.show(representation, selection)

    def hide(self, tuples):
        for representation, selection in tuples:
            cmd.hide(representation, selection)

    def hide_all(self):
        self.hide([('everything', '')])

    def label(self, tuples):
        for selection, expression in tuples:
            cmd.label(selection, expression)

    def add_pseudoatom(self, name, opts=None):
        opts = opts or {}
        cmd.pseudoatom(name, **opts)

    def sel_exists(self, name):
        return name in cmd.get_names("selections")

    def obj_exists(self, name):
        return name in cmd.get_names("objects")

    def group_exists(self, name):
        return name in cmd.get_names("group_objects")

    def get_coords(self, selection):
        return cmd.get_coords(selection)

    def distance(self, name, sel1, sel2):
        cmd.distance(name, sel1, sel2)

    def arrow(self, name, atm_sel1, atm_sel2, opts=None):
        opts = opts or {}
        cgo_arrow(atm_sel1, atm_sel2, name=name, **opts)

    def save_png(self, output_file, width=1200, height=1200, dpi=100, ray=1):
        cmd.png(output_file, width, height, dpi, ray)

    def save_session(self, output_file):
        if get_file_format(output_file) != 'pse':
            output_file += '.pse'
        cmd.save(output_file)

    def color(self, tuples):
        for color, selection in tuples:
            cmd.color(color, selection)

    def color_by_element(self, selections):
        for selection in selections:
            util.cbag(selection)

    def group(self, name, members, action=None):
        for member in members:
            if action:
                cmd.group(name, member, action)
            else:
                cmd.group(name, member)

    def set(self, name, value, opts=None):
        opts = opts or {}
        cmd.set(name, value, **opts)

    def align(self, mobile, target, opts=None):
        opts = opts or {}
        cmd.align(mobile, target, **opts)

    def delete(self, selections):
        for selection in selections:
            cmd.delete(selection)

    def remove(self, selections):
        for selection in selections:
            cmd.remove(selection)

    def extract(self, tuples):
        for name, selection in tuples:
            cmd.extract(name, selection)

    def reinitialize(self):
        cmd.reinitialize('everything')
        self.input_file = None

    def quit(self):
        cmd.quit()

    def run(self, func_name, opts):
        return getattr(cmd, func_name)(**opts)

    def run_cmds(self, commands):
        for func_name, opts in commands:
            getattr(cmd, func_name)(**opts)


class PymolSessionManager:

    def __init__(self, show_cartoon=False, bg_color="white", pharm_color=None,
                 add_directional_arrows=True, inter_color=PYMOL_INTERACTION_COLOR, pse_export_version="1.8"):

        self.show_cartoon = show_cartoon
        self.bg_color = bg_color
        self.pse_export_version = pse_export_version
        self.inter_color = inter_color
        self.wrapper = None
        self.add_directional_arrows = add_directional_arrows

    def new_session(self, data, output_file):
        self.start_session()
        self.set_view(data)
        self.save_session(output_file)
        self.finish_session()

    def start_session(self):
        self.wrapper = PymolWrapper()
        self.wrapper.set("pse_export_version", self.pse_export_version)
        self.wrapper.set("transparency_mode", 3)
        self.wrapper.set("group_auto_mode", 2)
        self.wrapper.run_cmds([("bg_color", {"color": self.bg_color})])
        self.wrapper.set("internal_gui_width", 350)

    def set_view(self, data):
        raise NotImplementedError("Use a class that implements this method.")

    def set_pdb_view(self, pdb_file, pdb_obj):
        prot_obj = "%s.prot" % pdb_obj
        self.wrapper.load(pdb_file, prot_obj)

        self.wrapper.extract([("%s.hets" % pdb_obj, "hetatm and %s" % prot_obj)])

        self.wrapper.color_by_element([pdb_obj])
        self.wrapper.hide([("everything", pdb_obj)])

        if self.show_cartoon:
            self.wrapper.show([("cartoon", pdb_obj)])

    def set_interactions_view(self, interactions, add_to_grp, uid):

        for i, inter in enumerate(interactions):
            if inter.type == "Proximal":
                continue

            obj1_name = "obj%d_grp%s" % (uid, hash(inter.src_grp))
            centroid_obj1 = inter.src_grp.centroid
            # Define the centroid in a nucleophile with two atoms as the position of its more electronegative atom.
            # Remember that the position in the interaction object matters. We have defined that the first group is always
            # the nucleophile for both dipole-dipole and ion-dipole interactions.
            if inter.type in NUCLEOPHILE_INTERS and len(inter.src_grp.atoms) == 2:
                dipole_atm = inter.src_grp.atoms[0] if (inter.src_grp.atoms[0].electronegativity >
                                                        inter.src_grp.atoms[1].electronegativity) else inter.src_grp.atoms[1]
                obj1_name += "_%s" % hash(dipole_atm.name)
                centroid_obj1 = dipole_atm.coord
            # For unfavorable multipolar interactions, it may happen that the first atom group is an electrophile as well.
            elif inter.type == "Unfavorable electrophile-electrophile" and len(inter.src_grp.atoms) == 2:
                dipole_atm = inter.src_grp.atoms[0] if (inter.src_grp.atoms[0].electronegativity <
                                                        inter.src_grp.atoms[1].electronegativity) else inter.src_grp.atoms[1]
                obj1_name += "_%s" % hash(dipole_atm.name)
                centroid_obj1 = dipole_atm.coord

            obj2_name = "obj%d_grp%s" % (uid, hash(inter.trgt_grp))
            centroid_obj2 = inter.trgt_grp.centroid
            # Define the centroid in an electrophile with two atoms as the position of its less electronegative atom.
            # Remember that the position in the interaction object matters. We have defined that the second group is always
            # the electrophile for both dipole-dipole and ion-dipole interactions.
            if inter.type in ELECTROPHILE_INTERS and len(inter.trgt_grp.atoms) == 2:
                dipole_atm = inter.trgt_grp.atoms[0] if (inter.trgt_grp.atoms[0].electronegativity <
                                                         inter.trgt_grp.atoms[1].electronegativity) else inter.trgt_grp.atoms[1]
                obj2_name += "_%s" % hash(dipole_atm.name)
                centroid_obj2 = dipole_atm.coord
            # For unfavorable multipolar interactions, it may happen that the second atom group is an nucleophile as well.
            elif inter.type == "Unfavorable nucleophile-nucleophile" and len(inter.trgt_grp.atoms) == 2:
                dipole_atm = inter.trgt_grp.atoms[0] if (inter.trgt_grp.atoms[0].electronegativity >
                                                         inter.trgt_grp.atoms[1].electronegativity) else inter.trgt_grp.atoms[1]
                obj2_name += "_%s" % hash(dipole_atm.name)
                centroid_obj2 = dipole_atm.coord

            # Add pseudoatoms
            if not self.wrapper.obj_exists(obj1_name):
                self.wrapper.add_pseudoatom(obj1_name, {"vdw": 1, "pos": list(centroid_obj1)})

            if not self.wrapper.obj_exists(obj2_name):
                self.wrapper.add_pseudoatom(obj2_name, {"vdw": 1, "pos": list(centroid_obj2)})

            # Set the representation for each compound in the groups involved in the interaction.
            for compound in inter.src_grp.compounds.union(inter.trgt_grp.compounds):
                comp_repr = "sphere" if compound.is_water() else "sticks"
                comp_sel = mybio_to_pymol_selection(compound)
                self.wrapper.show([(comp_repr, comp_sel)])
                carb_color = "green" if compound.is_target() else "gray"
                self.wrapper.color([(carb_color, comp_sel + " AND elem C")])

            # Check if the interaction involves the same compound: intramolecular interactions.
            inter_grp = "intra" if inter.is_intramol_interaction() else "inter"

            inter_name = "%s.all_inters.%s.%s.obj%d_inter%d.line" % (add_to_grp, inter_grp, INTERACTION_SHORT_NAMES[inter.type], uid, i)
            self.wrapper.distance(inter_name, obj1_name, obj2_name)
            self.wrapper.hide([("label", inter_name)])

            # Set styles to the interactions.
            self.wrapper.color([(self.inter_color.get_color(inter.type), inter_name)])

            if self.add_directional_arrows:
                if inter.type in UNFAVORABLE_INTERS:
                    arrow_name1 = "%s.all_inters.%s.%s.obj%d_inter%d.arrow1" % (add_to_grp, inter_grp,
                                                                                INTERACTION_SHORT_NAMES[inter.type], uid, i)
                    arrow_name2 = "%s.all_inters.%s.%s.obj%d_inter%d.arrow2" % (add_to_grp, inter_grp,
                                                                                INTERACTION_SHORT_NAMES[inter.type], uid, i)
                    square_name = "%s.all_inters.%s.%s.obj%d_inter%d.block" % (add_to_grp, inter_grp,
                                                                               INTERACTION_SHORT_NAMES[inter.type], uid, i)

                    arrow_opts = {"radius": 0.03, "gap": 0.9, "hlength": 0.5, "hradius": 0.2,
                                  "color": self.inter_color.get_color(inter.type)}
                    square_opts = {"radius": 0.3, "gap": 1.5, "hlength": 0, "hradius": 0,
                                   "color": self.inter_color.get_color(inter.type)}

                    # Two arrows in different directions
                    self.wrapper.arrow(arrow_name1, obj1_name, obj2_name, arrow_opts)

                    if not inter.is_directional():
                        self.wrapper.arrow(arrow_name2, obj2_name, obj1_name, arrow_opts)

                    # Add a square-like object
                    self.wrapper.arrow(square_name, obj1_name, obj2_name, square_opts)
                # Add arrows over the interaction lines to represent directional interactions
                elif inter.is_directional():
                    arrow_name = "%s.all_inters.%s.%s.obj%d_inter%d.arrow" % (add_to_grp, inter_grp,
                                                                              INTERACTION_SHORT_NAMES[inter.type], uid, i)
                    arrow_opts = {"radius": 0.03, "gap": 0.9, "hlength": 0.5, "hradius": 0.2,
                                  "color": self.inter_color.get_color(inter.type)}
                    self.wrapper.arrow(arrow_name, obj1_name, obj2_name, arrow_opts)

            # If a group object contains more than one atom, the centroid object will be displayed.
            if len(inter.src_grp.atoms) > 1:
                # Add the centroids to the group "grps" and append them to the main group
                self.wrapper.group("%s.grps" % add_to_grp, [obj1_name])
                self._set_centroid_style(obj1_name)
            # Otherwise, just remove the centroid as it will not add any new information (the atom represented
            # by the centroid is already been displayed).
            else:
                self.wrapper.delete([obj1_name])

            # If a group object contains more than one atom, the centroid object will be displayed.
            if len(inter.trgt_grp.atoms) > 1:
                # Add the centroids to the group "grps" and append them to the main group
                self.wrapper.group("%s.grps" % add_to_grp, [obj2_name])
                self._set_centroid_style(obj2_name)
            # Otherwise, just remove the centroid as it will not add any new information (the atom represented
            # by the centroid is already been displayed).
            else:
                self.wrapper.delete([obj2_name])

    def _set_centroid_style(self, centroid):
        # Set styles to the centroid.
        self.wrapper.hide([("nonbonded", centroid)])
        self.wrapper.show([("sphere", centroid)])
        self.wrapper.set("sphere_scale", 0.2, {"selection": centroid})

    def set_last_details_to_view(self):
        self.wrapper.set("sphere_scale", "0.3", {"selection": "visible and resn hoh"})
        self.wrapper.hide([("everything", "elem H")])

    def save_session(self, output_file):
        if not isinstance(self.wrapper, PymolWrapper):
            raise PymolSessionNotInitialized("No session was initialized.")

        self.wrapper.save_session(output_file)

    def finish_session(self):
        if not isinstance(self.wrapper, PymolWrapper):
            raise PymolSessionNotInitialized("No session was initialized.")

        self.wrapper.reinitialize()
        self.wrapper = None


def mybio_to_pymol_selection(entity):
    params = {}
    if entity.level == 'S' or entity.level == 'M':
        params[''] = 'all'
    else:
        if entity.level == 'C':
            params['chain'] = entity.id
        else:
            if entity.level == 'R':
                params['resn'] = entity.resname
                params['res'] = str(entity.id[1]) + entity.id[2].strip()
                params['chain'] = entity.get_parent().id
            else:
                if entity.level == 'A':
                    residue = entity.get_parent()
                    params['id'] = entity.serial_number
                    params['name'] = entity.name
                    params['resn'] = residue.resname
                    params['res'] = str(residue.id[1]) + residue.id[2].strip()
                    params['chain'] = residue.get_parent().id

            return (' AND ').join(['%s %s' % (k, str(v)) for k, v in params.items()])
