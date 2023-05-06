from luna.mol.entry import MolFileEntry
from luna.interaction.calc import InteractionsManager
from luna.wrappers.pymol import (PymolWrapper, PymolSessionManager,
                                 mybio_to_pymol_selection)
from luna.util.file import is_file_valid
from luna.util.exceptions import PymolSessionNotInitialized
from luna.MyBio.util import entity_to_string


class InteractionViewer(PymolSessionManager):

    """Class that inherits from
    :class:`~luna.wrappers.pymol.PymolSessionManager` and implements
    :meth:`set_view` to depict interactions in Pymol and save the view
    as a Pymol session.

    This class can be used to visualize multiple complexes into the same Pymol
    session, for instance, to compare binding modes. To do so, it is
    recommended that the protein structures are in the same coordinate system,
    i.e., they should be aligned first. This can be achieved with
    :func:`luna.align.tmalign.align_2struct` or any other tool of your
    preference.

    Parameters
    ----------
    show_hydrop_surface : bool
        If True, highlight hydrophobic surfaces. The default value is False.
    **kwargs : dict, optional
        Extra arguments to `InteractionViewer`.
        Refer to :class:`~luna.wrappers.pymol.PymolSessionManager`
        documentation for a list of all possible arguments.

    Examples
    --------

    In the below examples, we will assume a LUNA project object named
    ``proj_obj`` already exists.

    **Example 1)** In the first example, we will visualize all interactions
    identified in the first protein-ligand complex. To do so, we will provide
    a tuple containing an :class:`~luna.interaction.calc.InteractionsManager`
    from where the interactions will be recovered.

    First access the property ``interactions_mngrs`` to get an iterable of
    :class:`~luna.interaction.calc.InteractionsManager` objects and get the
    first one.

    >>> interactions_mngr = list(proj_obj.interactions_mngrs)[0]

    As `InteractionViewer` expects a list of tuples, we will define it now.
    The first item of the tuple is the :class:`~luna.mol.entry.Entry` instance,
    which represents a ligand. This can be obtained directly from the
    :class:`~luna.interaction.calc.InteractionsManager` object.
    The :class:`~luna.mol.entry.Entry` instance is necessary because the second
    item in the tuple (interactions) may be an iterable of
    :class:`~luna.interaction.type.InteractionType` from where such
    information cannot be recovered. Finally, the third item can be either a
    PDB file or a directory. In this example, we will use the PDB directory
    defined during the LUNA project initialization.

    >>> inter_tuples = [(interactions_mngr.entry,
    ...                  interactions_mngr,
    ...                  proj_obj.pdb_path)]

    Now, we create a new `InteractionViewer` object and call
    :meth:`~luna.wrappers.pymol.PymolSessionManager.new_session`, which will
    initialize the session, depict the interactions, and save the session to
    an output PSE file.

    >>> inter_view = InteractionViewer()
    >>> inter_view.new_session(inter_tuples, "output.pse")


    **Example 2)** In this example, we will create a new Pymol session where
    a given set of interactions will be shown. Let's say, for instance, we
    only want to visualize hydrogen bonds.

    To do so, instead of defining a tuple with an
    :class:`~luna.interaction.calc.InteractionsManager` instance, we will
    define a list of :class:`~luna.interaction.type.InteractionType` objects,
    which will contain only hydrogen bonds.

    Let's start with the selection of hydrogen bonds. Here, we will use the
    built-in method
    :meth:`luna.interaction.calc.InteractionsManager.filter_by_types`, which
    permits to filter interactions by type.

    >>> interactions_mngr = list(proj_obj.interactions_mngrs)[0]
    >>> hydrogen_bonds = interactions_mngr.filter_by_types(['Hydrogen bond'])

    Now, we just create the tuple as we did before and define the list of
    interactions we want to depict in the Pymol session.

    >>> inter_tuples = [(interactions_mngr.entry, \
    ...                  hydrogen_bonds,
    ...                  proj_obj.pdb_path)]

    Finally, we create a new `InteractionViewer` object and call
    :meth:`InteractionViewer.new_session`, which will initialize the session,
    depict the interactions, and save the session to an output PSE file.

    >>> inter_view = InteractionViewer()
    >>> inter_view.new_session(inter_tuples, "output.pse")

    **Example 3)** In this final example, we will create a new Pymol session to
    visualize all complexes in a LUNA project. To do so, we create a list with
    one tuple for each complex:

    >>> inter_tuples = [(im.entry, im, proj_obj.pdb_path) \
for im in proj_obj.interactions_mngrs]

    Then, as we did before, just create a new `InteractionViewer` object and
    call :meth:`InteractionViewer.new_session`.

    >>> inter_view = InteractionViewer()
    >>> inter_view.new_session(inter_tuples, "output.pse")

    **Note:** it is recommended that all complexes have the same atomic
    coordinates to make the analysis and comparisons easier. If that's
    not the case, you may want to align the structures first. To do so,
    you can use :func:`luna.align.tmalign.align_2struct` or any other tool
    of your preference.
    """

    def __init__(self, show_hydrop_surface=False, **kwargs):
        self.show_hydrop_surface = show_hydrop_surface
        super().__init__(**kwargs)

    def set_view(self, inter_tuples):
        """Depict interactions into the current Pymol session.

        Parameters
        ----------
        inter_tuples : iterable of tuple
            Each tuple must contain three items: an
            :class:`~luna.mol.entry.Entry` instance, an iterable of
            :class:`~luna.interaction.type.InteractionType` or an
            :class:`~luna.interaction.calc.InteractionsManager`,
            and a PDB file or the directory where the PDB file is located.
        """

        if not isinstance(self.wrapper, PymolWrapper):
            raise PymolSessionNotInitialized("No session was initialized.")

        for target_entry, inter_data, pathname in inter_tuples:

            if is_file_valid(pathname):
                pdb_file = pathname
            else:
                pdb_file = "%s/%s.pdb" % (pathname, target_entry.pdb_id)

            main_grp = target_entry.to_string(sep="-").replace("'", "-")

            mol_block = \
                (entity_to_string(target_entry.get_biopython_structure())
                 if isinstance(target_entry, MolFileEntry) else None)

            # Load PDB and extract hetatm.
            self.load_pdb(pdb_file, main_grp, mol_block)

            if isinstance(inter_data, InteractionsManager):
                interactions = inter_data.interactions
            else:
                interactions = inter_data

            # Add interactions and styles.
            interacting_residue_sels = \
                self.set_interactions_view(interactions, main_grp)
            if interacting_residue_sels:
                sel = " or ".join(interacting_residue_sels)
                self.wrapper.select(name="%s.inter_residues" % main_grp,
                                    selection=sel)

            if self.show_hydrop_surface:
                self.wrapper.color([("white", "%s and !name PS*" % main_grp)])

                # It will display all hydrophobic groups if an
                # AtomGroupsManager is available.
                atm_grp_mngr = interactions[0].src_grp.manager
                if atm_grp_mngr is not None:
                    interacting_atms = set()
                    for atm_grp in \
                        atm_grp_mngr.filter_by_types(["Hydrophobe",
                                                      "Hydrophobic"],
                                                     must_contain_all=False):

                        inters = [i for i in atm_grp.interactions
                                  if i.type == "Hydrophobic"]

                        for comp in atm_grp.compounds:
                            comp_sel = ("%s and %s"
                                        % (main_grp,
                                           mybio_to_pymol_selection(comp)))
                            self.wrapper.show([("sticks", comp_sel)])

                        if len(inters) > 0:
                            for inter in inters:
                                for atm in (inter.src_grp.atoms
                                            + inter.trgt_grp.atoms):
                                    sel = mybio_to_pymol_selection(atm)
                                    if atm not in interacting_atms:
                                        self.wrapper.color([("pink",
                                                             "%s and %s"
                                                             % (main_grp,
                                                                sel))])

                                for atm in (inter.src_interacting_atms
                                            + inter.trgt_interacting_atms):
                                    sel = mybio_to_pymol_selection(atm)
                                    self.wrapper.color([("hotpink",
                                                         "%s and %s"
                                                         % (main_grp, sel))])
                                    interacting_atms.add(atm)

                        else:
                            for atm in atm_grp.atoms:
                                sel = mybio_to_pymol_selection(atm)
                                self.wrapper.color([("wheat",
                                                     "%s and %s"
                                                     % (main_grp, sel))])

                    self.wrapper.hide([("sticks",
                                        "%s and not hetatm" % main_grp)])
                    self.wrapper.show([("sticks",
                                        "%s.inter_residues" % main_grp)])
                    self.wrapper.color([("white",
                                         "%s and not hetatm and not "
                                         "%s.inter_residues and !name PS*"
                                         % (main_grp, main_grp))])

                # Otherwise, it will display only hydrophobic groups comprising
                # the interacting groups.
                else:
                    interacting_atms = set()
                    for inter in interactions:
                        if inter.type != "Hydrophobic":
                            continue

                        for atm in inter.src_grp.atoms + inter.trgt_grp.atoms:
                            sel = mybio_to_pymol_selection(atm)
                            if atm not in interacting_atms:
                                self.wrapper.color([("pink", 
                                                     "%s and %s"
                                                     % (main_grp, sel))])

                        for atm in inter.src_interacting_atms + inter.trgt_interacting_atms:
                            sel = mybio_to_pymol_selection(atm)
                            self.wrapper.color([("hotpink",
                                                 "%s and %s"
                                                 % (main_grp, sel))])
                            interacting_atms.add(atm)

        self.set_last_details_to_view()
