from pymol import cmd
from pymol import util


from luna.wrappers.cgo_arrow import cgo_arrow
from luna.util.exceptions import (PymolSessionNotInitialized,
                                  IllegalArgumentError)
from luna.util.default_values import (PYMOL_INTERACTION_COLOR,
                                      INTERACTION_SHORT_NAMES)
from luna.util.file import get_filename, get_file_format

from Bio.Data.SCOPData import protein_letters_3to1

import logging

logger = logging.getLogger()

NUCLEOPHILE_INTERS = ["Orthogonal multipolar", "Parallel multipolar",
                      "Antiparallel multipolar", "Tilted multipolar",
                      "Multipolar", "Cation-nucleophile",
                      "Unfavorable anion-nucleophile",
                      "Unfavorable nucleophile-nucleophile"]

ELECTROPHILE_INTERS = ["Orthogonal multipolar", "Parallel multipolar",
                       "Antiparallel multipolar", "Tilted multipolar",
                       "Multipolar", "Anion-electrophile",
                       "Unfavorable cation-electrophile",
                       "Unfavorable electrophile-electrophile"]

UNFAVORABLE_INTERS = ["Repulsive", "Unfavorable anion-nucleophile",
                      "Unfavorable cation-electrophile",
                      "Unfavorable nucleophile-nucleophile",
                      "Unfavorable electrophile-electrophile"]


class PymolWrapper:
    """This class provides functions to provide easy access
    to common functions from Pymol."""

    def get_cmd(self):
        """Expose the Pymol ``cmd`` object, so that one can
        call Pymol functions directly."""
        return cmd

    def load(self, input_file, obj_name=None):
        """Load a molecular file (e.g., PDB files).

        Parameters
        ----------
        input_file : str
            The pathname of the molecular file to be loaded.
        obj_name : str, optional
            Pymol object to store the loaded structure.
            If not provided, the filename will be used instead.
        """
        self.input_file = input_file
        if not obj_name:
            obj_name = get_filename(input_file)
        cmd.load(input_file, obj_name)

    def show(self, tuples):
        """Display atom and bond representations for certain selections.

        Parameters
        ----------
        tuples : iterable of tuple
            Each tuple should contain a Pymol representation (e.g., 'sticks')
            and a selection (e.g., 'hetatm').
        """
        for representation, selection in tuples:
            cmd.show(representation, selection)

    def hide(self, tuples):
        """Hide atom and bond representations for certain selections.

        Parameters
        ----------
        tuples : iterable of tuple
            Each tuple should contain a Pymol representation (e.g., 'sticks')
             and a selection (e.g., 'hetatm').
        """
        for representation, selection in tuples:
            cmd.hide(representation, selection)

    def hide_all(self):
        """Hide all representations."""
        self.hide([('everything', '')])

    def center(self, selection):
        """Translate the window, the clipping slab, and the origin to a point
        centered within the selection."""
        cmd.center(selection)

    def label(self, tuples):
        """Draw text labels for PyMOL objects.

        Parameters
        ----------
        tuples : iterable of tuple
            Each tuple should contain a selection (e.g., 'hetatm') and some
            string to label the given selection.
        """
        for selection, expression in tuples:
            cmd.label(selection, expression)

    def add_pseudoatom(self, name, opts=None):
        """Create a molecular object with a pseudoatom or add a pseudoatom to
        a molecular object if the specified object already exists.

        Parameters
        ----------
        name : str
            The object name to create or modify.
        opts : dict
            A set of options to create the pseudoatom.
            Check `Pymol <https://pymolwiki.org/index.php/Pseudoatom>`_ to
            discover which options are available.

        """
        opts = opts or {}
        cmd.pseudoatom(name, **opts)

    def select(self, selection, name="sele", enable=0):
        """Create a named selection from an atom selection.

        Parameters
        ----------
        selection : str
            The expression to select atoms.
        name : str
            The selection name, which by default is 'sele'.
        enable : {0, 1}
            If ``1``, activate the selection, i.e., show selection indicators.
            The default value is 0, which implies the selection
            indicators won't be shown.
        """
        cmd.select(name, selection, enable=enable)

    def get_names(self, obj_type):
        """Get names of objects, grouped objects, or selections.

        Parameters
        ----------
        obj_type : {'objects', 'selections', 'all', 'public_objects', \
            'public_selections', 'public_nongroup_objects', \
            'public_group_objects', 'nongroup_objects', 'group_objects'}
                The target object type.

        Returns
        -------
         : list of str
        """
        return cmd.get_names(obj_type)

    def sel_exists(self, name):
        """Check if a selection exists given by its name ``name``."""
        return name in cmd.get_names("selections")

    def obj_exists(self, name):
        """Check if an object exists given by its name ``name``."""
        return name in cmd.get_names("objects")

    def group_exists(self, name):
        """Check if a group of objects exists given by its name ``name``."""
        return name in cmd.get_names("group_objects")

    def get_coords(self, selection):
        """Get atomic coordinates for a given atom selection.

        Parameters
        ----------
        selection : str
            The expression to select atoms.

        Returns
        -------
         : array_like of float (size 3)
            Atomic coordinates (x, y, z) of each atom selected.
        """
        return cmd.get_coords(selection)

    def distance(self, name, sel1, sel2):
        """Create a new distance object between two atoms given by their
        selection-expressions.

        Parameters
        ----------
        name : str
            Name of the distance object to create.
        sel1 : str
            The expression to select the first atom.
        sel2 : str
            The expression to select the second atom.
        """
        cmd.distance(name, sel1, sel2)

    def arrow(self, name, atm_sel1, atm_sel2, opts=None):
        """Draw an arrow object between two atoms given by their
        selection-expressions.

        Parameters
        ----------
        name : str
            Name of the arrow object to create.
        sel1 : str
            The expression to select the first atom.
        sel2 : str
            The expression to select the second atom.
        opts : dict
            A set of options to create the arrow.
            Check `Pymol <https://pymolwiki.org/index.php/Cgo_arrow>`_ to
            discover which options are available.
        """
        opts = opts or {}
        cgo_arrow(atm_sel1, atm_sel2, name=name, **opts)

    def save_png(self, output_file, width=1200, height=1200, dpi=100, ray=1):
        """Save the current Pymol session as a PNG format image file.

        Parameters
        ----------
        output_file : str
            The output image pathname.
        width : int
            The width in pixels. The default value is 1,200.
        height : int or str
            The height in pixels. The default value is 1,200.
        dpi : float
            Dots-per-inch. The default value is 100.
        ray : {0, 1}
            If ``1`` (the default), run ray first to make high-resolution
            photos.
        """
        cmd.png(output_file, width, height, dpi, ray)

    def save_session(self, output_file):
        """Save the current PyMOL state to a PSE format file to later use.

        Parameters
        ----------
        output_file : str
            The output pathname.
        """
        if get_file_format(output_file) != 'pse':
            output_file += '.pse'
        cmd.save(output_file)

    def color(self, tuples):
        """Color objects and atoms.

        Parameters
        ----------
        tuples : iterable of tuple
            Each tuple should contain a color (e.g., 'red') and a selection
            (e.g., 'hetatm').
        """
        for color, selection in tuples:
            cmd.color(color, selection)

    def color_by_element(self, selections, c_color="green"):
        """Color atoms by their default element color (e.g., oxygen in red).

        Parameters
        ----------
        selections : iterable of str
            A sequence of selections to define which atoms will be colored by
            element.
        c_color : {'green', 'cyan', 'light magenta', 'yellow', 'salmon', \
                    'white', 'slate', 'bright orange', 'purple', 'pink'}
            The carbon color. The default value is 'green'.
        """

        valid_colors = ['green', 'cyan', 'light magenta', 'yellow', 'salmon',
                        'white', 'slate', 'bright orange', 'purple', 'pink']
        if c_color.lower() not in valid_colors:
            raise IllegalArgumentError("Invalid color '%s'. The accepted "
                                       "colors are: %s."
                                       % (c_color, ", ".join(valid_colors)))

        c_color = c_color.lower()

        for selection in selections:
            if c_color == 'green':
                util.cbag(selection)
            elif c_color == 'cyan':
                util.cbac(selection)
            elif c_color == 'light magenta':
                util.cbam(selection)
            elif c_color == 'yellow':
                util.cbay(selection)
            elif c_color == 'salmon':
                util.cbas(selection)
            elif c_color == 'white':
                util.cbaw(selection)
            elif c_color == 'slate':
                util.cbab(selection)
            elif c_color == 'bright orange':
                util.cbao(selection)
            elif c_color == 'purple':
                util.cbap(selection)
            elif c_color == 'pink':
                util.cbak(selection)

    def group(self, name, members, action=None):
        """Create or update a group object.

        Parameters
        ----------
        name : str
            The group name to create or update.
        members : iterable of str
            The objects to include in the group.
        action : {'add', 'remove', 'open', 'close', 'toggle', 'auto', \
                    'empty', 'purge', 'excise'}, optional
            An action to take. If not provided, the default value
            'auto' will be used instead. The description of the actions
            are described below (source: Pymol documentation):
                * add:     add members to group.
                * remove:  remove members from group (members will be \
                    ungrouped).
                * empty:   remove all members from group.
                * purge:   remove all members from group and delete them.
                * excise:  remove all members from group and delete group.
                * open:    expand group display in object menu panel.
                * close:   collapse group display in object menu panel.
                * toggle:  toggle group display in object menu panel.
                * auto:    add or toggle.
        """
        for member in members:
            if action:
                cmd.group(name, member, action)
            else:
                cmd.group(name, member)

    def create(self, name, selection):
        """Create a new molecular object from a selection.

        Note that the selected atoms won't be extracted from the original
        object. Instead, a copy of them will be created in the new object.

        Parameters
        ----------
        name : str
            The object name to be created.
        selection : str
            The expression to select atoms.
        """
        cmd.create(name, selection)

    def alter(self, selection, expression):
        """Modify atomic properties.

        Parameters
        ----------
        selection : str
            The expression to select atoms.
        expression : str
            Expression in Python language to define which properties should be
            modified. This can be used, for instance, to rename an atom or
            chain.

        Examples
        --------

        Alter the name of a ligand carbon from 'C1' to 'CA'.

        >>> pw_obj.alter("hetatm and name C1", "name='CA'")

        Alter the chain A name to 'Z'.

        >>> pw_obj.alter("chain A", "chain='Z'")
        """
        cmd.alter(selection, expression)

    def set(self, name, value, opts=None):
        """Modify global, object, object-state, or per-atom settings.

        Parameters
        ----------
        name : str
            The setting name to modify.
        value : str
            The new setting value.
        opts : dict
            A set of options.
            Check `Pymol <https://pymolwiki.org/index.php/Pseudoatom>`_
            to discover which options are available.
        """
        opts = opts or {}
        cmd.set(name, value, **opts)

    def align(self, mobile, target, opts=None):
        """Align the structure ``mobile`` to the ``target`` structure.

        Parameters
        ----------
        mobile : str
            The structure to be aligned given by an atomic selection.
        target : str
            The target structure given by an atomic selection.
        opts : dict
            Alignment options.
            Check `Pymol <https://pymolwiki.org/index.php/Align>`_
            to discover which options are available.
        """
        opts = opts or {}
        cmd.align(mobile, target, **opts)

    def delete(self, selections):
        """Delete the provided selections.

        Parameters
        ----------
        selections : iterable of str
            A sequence of selections to be deleted.
            Wildcards can be used to define object or selection names.
        """
        for selection in selections:
            cmd.delete(selection)

    def remove(self, selections):
        """Remove a selection of atoms from models.

        Parameters
        ----------
        selections : iterable of str
            A sequence of selections to define which atoms will be removed.
        """
        for selection in selections:
            cmd.remove(selection)

    def extract(self, tuples):
        """Perform multiple extractions, i.e., extract atoms from an object to
        another object.

        Parameters
        ----------
        tuples : iterable of tuple
            Each tuple should contain the object name to where atoms will be
            added and the selection itself that defines which atoms will be
            extracted (e.g., 'hetatm').
        """
        for name, selection in tuples:
            cmd.extract(name, selection)

    def load_mol_from_pdb_block(self, pdb_block, obj_name):
        """Load a molecular file from a PDB block string.

        Parameters
        ----------
        pdb_block :str
            The PDB block string.
        obj_name : str
            Pymol object to store the loaded structure.
        """
        cmd.read_pdbstr(pdb_block, obj_name, state=1)

    def reinitialize(self):
        """Clear all objects and resets all parameters to default."""
        cmd.reinitialize('everything')
        self.input_file = None

    def quit(self):
        """Terminate Pymol."""
        cmd.quit()

    def run(self, func_name, opts):
        """Run a Pymol command ``func_name`` with parameters ``opts``.

        Parameters
        ----------
        func_name : str
            The Pymol command name.
        opts : dict
            Parameters to pass to the command.
        """
        return getattr(cmd, func_name)(**opts)

    def run_cmds(self, commands):
        """Run a set of Pymol commands.

        Parameters
        ----------
        commands : iterable of tuple
            Each tuple should contain a Pymol command and its parameters.
            See :meth:`run` to more details.
        """

        for func_name, opts in commands:
            getattr(cmd, func_name)(**opts)


class PymolSessionManager:

    """Class to start, manage, and save Pymol sessions.
    This class provides useful built-in functions to load PDB/Mol files and
    show interactions.

    .. note::
        This class is not intended to be used directly because :meth:`set_view`
        is not implemented by default. Instead, you should use a class that
        inherits from `PymolSessionManager` and implements :meth:`set_view`.
        An example is the class
        :class:`~luna.interaction.view.InteractionViewer` that implements a
        custom :meth:`set_view` to show interactions. Therefore, you should
        define your own logic beyond :meth:`set_view` to save a Pymol session
        that meets your goals.

    Parameters
    ----------
    show_cartoon : bool
        If True, show the protein structure as cartoons.
    bg_color : str
        The background color. The default value is "white".
        Check `Pymol <https://pymolwiki.org/index.php/Color_Values>`_
        to discover which colors are available.
    add_directional_arrows : bool
        If True, show arrows for directional interactions (e.g., hydrogen bonds
        and multipolar interactions).
    show_res_labels : bool
        If True (the default), show residue labels.
    inter_color : :class:`~luna.util.ColorPallete`
        A Pymol-compatible color scheme for interactions.
        The default value is
        :const:`~luna.util.default_values.PYMOL_INTERACTION_COLOR`.
    pse_export_version : str
        Define a legacy format for saving Pymol sessions (PSE files).
        The default value os '1.8'.
    """

    def __init__(self,
                 show_cartoon=False,
                 bg_color="white",
                 add_directional_arrows=True,
                 show_res_labels=True,
                 inter_color=PYMOL_INTERACTION_COLOR,
                 pse_export_version="1.8"):

        self.show_cartoon = show_cartoon
        self.bg_color = bg_color
        self.pse_export_version = pse_export_version
        self.inter_color = inter_color
        self.add_directional_arrows = add_directional_arrows
        self.show_res_labels = show_res_labels
        self.wrapper = None

    def new_session(self, data, output_file):
        """Start a new session, which includes the following steps:

            * Start a new Pymol session (:meth:`start_session`);
            * Set the view (:meth:`set_view`);
            * Save the Pymol session to ``output_file``;
            * Finish the Pymol session.

        Parameters
        ----------
        data : iterable
            Data to be processed by :meth:`set_view`.
        output_file : str
            The pathname to where the Pymol session will be saved.
        """
        self.start_session()
        self.set_view(data)
        self.save_session(output_file)
        self.finish_session()

    def start_session(self):
        """Start a new session and set Pymol settings, including the
        background color and the PSE export version.
        """
        self.wrapper = PymolWrapper()
        self.wrapper.set("pse_export_version", self.pse_export_version)
        self.wrapper.set("transparency_mode", 3)
        self.wrapper.set("group_auto_mode", 2)
        self.wrapper.run_cmds([("bg_color", {"color": self.bg_color})])
        self.wrapper.set("internal_gui_width", 370)

    def set_view(self, data):
        """Set the session view. However, this method is not implemented by
        default. Instead, you should use a class that inherits from
        `PymolSessionManager` and implements :meth:`set_view`. An example is
        the class :class:`~luna.interaction.view.InteractionViewer` that
        implements a custom :meth:`set_view` to show interactions.
        Therefore, you should define your own logic beyond :meth:`set_view` to
        save a Pymol session that meets your goals.

        Parameters
        ----------
        data : iterable
            The data that will be used to set the Pymol view.
        """
        raise NotImplementedError("Use a class that implements this method.")

    def load_pdb(self,
                 pdb_file,
                 pdb_obj,
                 mol_block=None,
                 is_ftmap_output=False):
        """Load molecules from PDB files to the current Pymol session.

        Optionally, ligands can also be loaded from a separate molecular string
        block. This is especially useful when working with docked molecules in
        which the protein structure is in a PDB file and ligands are in a
        separate molecular file.

        Parameters
        ----------
        pdb_file : str
            The pathname of the PDB file to be loaded.
        pdb_obj : str
            Pymol object to store the loaded structure.
        mol_block : str, optional
            A molecular string block to load together with the PDB file.
        is_ftmap_output : bool
            If the PDB file is an FTMap output.
            If so, an additional processing step is performed to
            standardize the loaded Pymol objects.
        """
        prot_obj = "%s.prot" % pdb_obj

        self.wrapper.load(pdb_file, prot_obj)

        if is_ftmap_output:
            objs = self.wrapper.get_names("objects")
            to_merge = []
            for obj in objs:
                if obj == "protein":
                    to_merge.append(obj)

                    sel = \
                        ("protein and not resn %s"
                         % " and not resn ".join(protein_letters_3to1.keys()))
                    self.wrapper.alter(sel, "type='HETATM'")

                elif obj.endswith(".pdb"):
                    to_merge.append(obj)
                    self.wrapper.alter(obj, "type='HETATM'")

            self.wrapper.create(prot_obj, " | ".join(to_merge))
            self.wrapper.delete(objs)

        self.wrapper.extract([("%s.hets" % pdb_obj,
                               "hetatm and %s" % prot_obj)])

        if mol_block is not None:
            self.wrapper.load_mol_from_pdb_block(mol_block,
                                                 "%s.hets" % pdb_obj)

        self.wrapper.color_by_element([pdb_obj])
        self.wrapper.hide([("everything", pdb_obj)])

        if self.show_cartoon:
            self.wrapper.show([("cartoon", pdb_obj)])

    def set_interactions_view(self,
                              interactions,
                              main_grp,
                              secondary_grp=None):
        """Display molecular interactions.

        Parameters
        ----------
        interactions : iterable of \
                :class:`~luna.interaction.type.InteractionType`
            A sequence of interactions to show.
        main_grp : str
            Main Pymol object to store atom groups.
        secondary_grp : str, optional
            Secondary Pymol object to store interactions.
            If not provided, ``main_grp`` will be used instead.
        """

        residue_selections = set()

        secondary_grp = secondary_grp or main_grp

        for i, inter in enumerate(interactions):

            #
            # Centroid 1
            #
            centroid_hash1 = hash(tuple(sorted(inter.src_interacting_atms)))
            obj1_name = "%s.centroids.%s" % (main_grp, centroid_hash1)
            centroid_obj1 = inter.src_centroid
            centroid_obj1_visible = True
            # Define the centroid in a nucleophile with two atoms as the
            # position of its more electronegative atom. Remember that
            # the position in the interaction object matters. We have
            # defined that the first group is always the nucleophile for
            # both dipole-dipole and ion-dipole interactions.
            if (inter.type in NUCLEOPHILE_INTERS
                    and len(inter.src_grp.atoms) == 2):
                dipole_atm = \
                    (inter.src_grp.atoms[0]
                     if (inter.src_grp.atoms[0].electronegativity
                         > inter.src_grp.atoms[1].electronegativity)
                     else inter.src_grp.atoms[1])
                obj1_name += "_%s" % hash(dipole_atm.name)
                centroid_obj1 = dipole_atm.coord
                centroid_obj1_visible = False
            # For unfavorable multipolar interactions, it may happen that the
            # first atom group is an electrophile as well.
            elif (inter.type == "Unfavorable electrophile-electrophile"
                    and len(inter.src_grp.atoms) == 2):
                dipole_atm = \
                    (inter.src_grp.atoms[0]
                     if (inter.src_grp.atoms[0].electronegativity
                         < inter.src_grp.atoms[1].electronegativity)
                     else inter.src_grp.atoms[1])
                obj1_name += "_%s" % hash(dipole_atm.name)
                centroid_obj1 = dipole_atm.coord
                centroid_obj1_visible = False

            #
            # Centroid 2
            #
            centroid_hash2 = hash(tuple(sorted(inter.trgt_interacting_atms)))
            obj2_name = "%s.centroids.%s" % (main_grp, centroid_hash2)
            centroid_obj2 = inter.trgt_centroid
            centroid_obj2_visible = True
            # Define the centroid in an electrophile with two atoms as
            # the position of its less electronegative atom. Remember that
            # the position in the interaction object matters. We have defined
            # that the second group is always the electrophile for both
            # dipole-dipole and ion-dipole interactions.
            if (inter.type in ELECTROPHILE_INTERS
                    and len(inter.trgt_grp.atoms) == 2):
                dipole_atm = \
                    (inter.trgt_grp.atoms[0]
                     if (inter.trgt_grp.atoms[0].electronegativity
                         < inter.trgt_grp.atoms[1].electronegativity)
                     else inter.trgt_grp.atoms[1])
                obj2_name += "_%s" % hash(dipole_atm.name)
                centroid_obj2 = dipole_atm.coord
                centroid_obj2_visible = False
            # For unfavorable multipolar interactions, it may happen that the
            # second atom group is a nucleophile as well.
            elif (inter.type == "Unfavorable nucleophile-nucleophile"
                    and len(inter.trgt_grp.atoms) == 2):
                dipole_atm = \
                    (inter.trgt_grp.atoms[0]
                     if (inter.trgt_grp.atoms[0].electronegativity
                         > inter.trgt_grp.atoms[1].electronegativity)
                     else inter.trgt_grp.atoms[1])
                obj2_name += "_%s" % hash(dipole_atm.name)
                centroid_obj2 = dipole_atm.coord
                centroid_obj2_visible = False

            # Add pseudoatoms
            if not self.wrapper.obj_exists(obj1_name):
                self.wrapper.add_pseudoatom(obj1_name,
                                            {"vdw": 1,
                                             "pos": list(centroid_obj1)})
            if not self.wrapper.obj_exists(obj2_name):
                self.wrapper.add_pseudoatom(obj2_name,
                                            {"vdw": 1,
                                             "pos": list(centroid_obj2)})

            # Set the representation for each compound in the groups involved
            # in the interaction.
            compounds = inter.src_grp.compounds.union(inter.trgt_grp.compounds)
            for compound in compounds:
                if compound.is_water() or compound.is_metal():
                    comp_repr = "sphere"
                elif compound.is_hetatm():
                    if (len(compound.child_list) == 1
                            or len([atm for atm in compound.child_list
                                    if atm.element != "H"]) == 1):
                        comp_repr = "sphere"
                    else:
                        comp_repr = "sticks"
                else:
                    comp_repr = "sticks"

                comp_sel = ("%s and %s"
                            % (main_grp, mybio_to_pymol_selection(compound)))
                self.wrapper.show([(comp_repr, comp_sel)])
                carb_color = "green" if compound.is_target() else "gray"
                self.wrapper.color([(carb_color, comp_sel + " AND elem C")])

                if compound.is_residue():
                    residue_selections.add(comp_sel)

                    # Show residue label if required.
                    if self.show_res_labels:
                        self.wrapper.label([("%s AND name CA" % comp_sel,
                                             '"%s-%s" % (resn, resi)')])

            # Check if the interaction involves the same
            # compound: intramolecular interactions.
            inter_grp = "intra" if inter.is_intramol_interaction() else "inter"

            src_grp_name = "+".join(["%s-%s-%d%s"
                                     % (c.parent.id, c.resname,
                                        c.id[1], c.id[2].strip())
                                     for c in sorted(inter.src_grp.compounds)])

            trgt_grp_name = \
                "+".join(["%s-%s-%d%s" % (c.parent.id, c.resname,
                                          c.id[1], c.id[2].strip())
                          for c in sorted(inter.trgt_grp.compounds)])

            inter_grp = ("%s.all_inters.%s.%s.i%d_%s_and_%s"
                         % (secondary_grp, inter_grp,
                            INTERACTION_SHORT_NAMES[inter.type],
                            i, src_grp_name, trgt_grp_name))

            inter_name = "%s.line" % inter_grp

            self.wrapper.distance(inter_name, obj1_name, obj2_name)
            self.wrapper.hide([("label", inter_name)])

            inter_color = self.inter_color.get_color(inter.type)

            # Set styles to the interactions.
            self.wrapper.color([(inter_color,
                                 inter_name)])

            if self.add_directional_arrows:

                if inter.type in UNFAVORABLE_INTERS:
                    arrow_name1 = "%s.arrow1" % inter_grp
                    arrow_name2 = "%s.arrow2" % inter_grp

                    square_name = "%s.block" % inter_grp

                    arrow_opts = {"radius": 0.03, "gap": 0.9,
                                  "hlength": 0.5, "hradius": 0.2,
                                  "color": inter_color}
                    square_opts = {"radius": 0.3,
                                   "gap": 1.5, "hlength": 0,
                                   "hradius": 0, "color": inter_color}

                    # Two arrows in different directions
                    self.wrapper.arrow(arrow_name1,
                                       obj1_name,
                                       obj2_name,
                                       arrow_opts)

                    if not inter.is_directional():
                        self.wrapper.arrow(arrow_name2,
                                           obj2_name,
                                           obj1_name,
                                           arrow_opts)
                    self.wrapper.arrow(arrow_name2,
                                       obj2_name,
                                       obj1_name,
                                       arrow_opts)

                    # Add a square-like object
                    self.wrapper.arrow(square_name, obj1_name,
                                       obj2_name, square_opts)

                # Add arrows over the interaction lines to
                # represent directional interactions
                elif inter.is_directional():
                    arrow_name = "%s.arrow" % inter_grp
                    arrow_opts = {"radius": 0.03, "gap": 0.9,
                                  "hlength": 0.5, "hradius": 0.2,
                                  "color": inter_color}
                    self.wrapper.arrow(arrow_name,
                                       obj1_name,
                                       obj2_name,
                                       arrow_opts)

            # If a group object contains more than one atom.
            if inter.src_grp.size > 1 and centroid_obj1_visible:
                # Add the centroids to the group "centroids" and append them
                # to the main group
                self._set_centroid_style(obj1_name)
            # Otherwise, just remove the centroid as it will not add any new
            # information (the atom represented by the centroid is already the
            # atom itself).
            else:
                self.wrapper.delete([obj1_name])

            # If a group object contains more than one atom.
            if inter.trgt_grp.size > 1 and centroid_obj2_visible:
                # Add the centroids to the group "centroids" and append them
                # to the main group
                self._set_centroid_style(obj2_name)
            # Otherwise, just remove the centroid as it will not add any new
            # information (the atom represented by the centroid is already
            # been displayed).
            else:
                self.wrapper.delete([obj2_name])

        return residue_selections

    def _set_centroid_style(self, centroid):
        # Set styles to the centroid.
        self.wrapper.hide([("nonbonded", centroid)])
        self.wrapper.show([("sphere", centroid)])
        self.wrapper.set("sphere_scale", 0.2, {"selection": centroid})

    def set_last_details_to_view(self):
        """This method can be called to apply final modifications to the
        Pymol session. In its default version, the following modifications
        are applied:
            * Dash radius for interactions is set to 0.08;
            * Labels are set to bold and their size is set to 20;
            * Atomic spheres' scale is set to 0.3;
            * Hydrogen atoms are hidden;
            * The view is centered within the visible objects.
        """
        self.wrapper.set("dash_radius", 0.08)
        self.wrapper.set("label_font_id", "13")
        self.wrapper.set("label_size", "20")
        self.wrapper.set("sphere_scale", "0.3",
                         {"selection": "visible and not name PS*"})
        self.wrapper.hide([("everything", "elem H+D")])
        self.wrapper.center("visible")

    def save_session(self, output_file):
        """Save the Pymol session as a PSE file of name ``output_file``."""
        if not isinstance(self.wrapper, PymolWrapper):
            raise PymolSessionNotInitialized("No session was initialized.")

        self.wrapper.save_session(output_file)

    def finish_session(self):
        """Clear all objects and resets all parameters to default."""
        if not isinstance(self.wrapper, PymolWrapper):
            raise PymolSessionNotInitialized("No session was initialized.")

        self.wrapper.reinitialize()
        self.wrapper = None


def mybio_to_pymol_selection(entity):
    """Transform an :class:`~luna.MyBio.PDB.Entity.Entity` instance into a
    Pymol selection-expression, which can then be used to select atoms in a
    Pymol session.

    Parameters
    ----------
    entity : :class:`~luna.MyBio.PDB.Entity.Entity`
        An entity to be transformed into a Pymol selection-expression.

    Returns
    -------
     : str
        The Pymol selection-expression.

    Examples
    --------

    First, let's parse a PDB file to work with.

    >>> from luna.util.default_values import LUNA_PATH
    >>> from luna.MyBio.PDB.PDBParser import PDBParser
    >>> pdb_parser = PDBParser(PERMISSIVE=True, QUIET=True)
    >>> structure = pdb_parser.get_structure("Protein",
    ...                                      f"{LUNA_PATH}/tutorial/\
inputs/3QQK.pdb")

    Now, let's get the Pymol selection-expression for the chain A.

    >>> from luna.wrappers.pymol import mybio_to_pymol_selection
    >>> print(mybio_to_pymol_selection(structure[0]['A']))
    chain A

    Finally, we can get the Pymol selection-expression for the ligand X02.

    >>> from luna.wrappers.pymol import mybio_to_pymol_selection
    >>> print(mybio_to_pymol_selection(structure[0]["A"][('H_X02', 497, ' ')]))
    resn X02 AND res 497 AND chain A
    """
    params = {}
    if entity.level == 'S' or entity.level == 'M':
        params[''] = 'all'
    elif entity.level == 'C':
        params['chain'] = entity.id
    elif entity.level == 'R':
        params['resn'] = entity.resname
        params['res'] = str(entity.id[1]) + entity.id[2].strip()
        params['chain'] = entity.get_parent().id
    elif entity.level == 'A':
        residue = entity.get_parent()
        params['id'] = entity.serial_number
        params['name'] = entity.name
        params['resn'] = residue.resname
        params['res'] = str(residue.id[1]) + residue.id[2].strip()
        params['chain'] = residue.get_parent().id
    else:
        return {}

    if "resn" in params:
        # Escape characters that can generate problems with Pymol
        params['resn'] = params['resn'].replace("+", "\\+")
        params['resn'] = params['resn'].replace("-", "\\-")

    return (' AND ').join(['%s %s' % (k, str(v)) for k, v in params.items()])
