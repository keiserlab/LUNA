from luna.mol.entry import MolFileEntry
from luna.wrappers.pymol import PymolWrapper, PymolSessionManager
from luna.MyBio.util import entity_to_string
from luna.util.exceptions import PymolSessionNotInitialized
from luna.util.file import is_file_valid


import logging

logger = logging.getLogger()


class ShellViewer(PymolSessionManager):
    """Class that inherits from
    :class:`~luna.wrappers.pymol.PymolSessionManager` and implements
    :meth:`set_view` to depict shells (IFP features) and interactions
    in Pymol and save the view as a Pymol session.

    This class can be used to visualize multiple complexes into the same Pymol
    session, for instance, to compare binding modes and analyze similar shells.
    To do so, it is recommended that the protein structures are in the same
    coordinate system, i.e., they should be aligned first. This can be achieved
    with :func:`luna.align.tmalign.align_2struct` or any other tool of
    your preference.

    Examples
    --------

    In the below example, we will assume a LUNA project object named
    ``proj_obj`` already exists.

    To visualize shells, we first need to generate them. So, let's define a
    :class:`~luna.interaction.fp.shell.ShellGenerator` object that will create
    shells over 2 iterations (levels). At each iteration, the shell radius will
    be increased by 3 and substructural information will be encoded following
    EIFP definition. Here, as an example, we will generate shells for the first
    :class:`~luna.mol.groups.AtomGroupsManager` object at ``proj_obj``.

    >>> from luna.interaction.fp.shell import ShellGenerator
    >>> from luna.interaction.fp.type import IFPType
    >>> num_levels, radius_step = 2, 3
    >>> sg = ShellGenerator(num_levels, radius_step, ifp_type=IFPType.EIFP)
    >>> atm_grps_mngr = list(proj_obj.atm_grps_mngrs)[0]
    >>> sm = sg.create_shells(atm_grps_mngr)

    The function
    :meth:`~luna.interaction.fp.shell.ShellGenerator.create_shells` returns a
    :class:`~luna.interaction.fp.shell.ShellManager` object, which provides
    built-in methods to access the created shells. Therefore, you can interact
    with it to select the shells you want to visualize in Pymol. As an example,
    let's select all unique shells at the last level:

    >>> shells = sm.get_shells_by_level(num_levels - 1, unique_shells=True)

    .. note::
        Levels are 0-indexed. So, the first level is 0, second is 1, etc.
        That means if ``num_levels`` is 5, the last level will be 4.

    As `ShellViewer` expects a list of tuples, we will define it now. The first
    item of the tuple is the :class:`~luna.mol.entry.Entry` instance, which
    represents a ligand. This can be obtained directly from the
    :class:`~luna.mol.groups.AtomGroupsManager` object. The second item is the
    list of shells you want to visualize. Finally, the third item can be either
    a PDB file or a directory. In this example, we will use the PDB directory
    defined during the LUNA project initialization.

    >>> shell_tuples = [(atm_grps_mngr.entry, shells, proj_obj.pdb_path)]

    To finish, we now create a new `ShellViewer` object and call
    :meth:`~luna.wrappers.pymol.PymolSessionManager.new_session`, which will
    initialize the session, depict shells and interactions, and save the
    session to an output PSE file.

    >>> from luna.interaction.fp.view import ShellViewer
    >>> sv = ShellViewer()
    >>> sv.new_session(shell_tuples, "example.pse")
    """

    def set_view(self, shell_tuples):
        """Depict shells (IFP features) and interactions into the current
        Pymol session.

        Parameters
        ----------
        shell_tuples : iterable of tuple
            Each tuple must contain three items: an
            :class:`~luna.mol.entry.Entry` instance,
            an iterable of :class:`~luna.interaction.fp.shell.Shell`,
            and a PDB file or the directory where the PDB file is located.
        """

        if not isinstance(self.wrapper, PymolWrapper):
            raise PymolSessionNotInitialized("No session was initialized.")

        for target_entry, shells, pathname in shell_tuples:

            if is_file_valid(pathname):
                pdb_file = pathname
            else:
                pdb_file = "%s/%s.pdb" % (pathname, target_entry.pdb_id)

            main_grp = target_entry.to_string(sep="-")

            mol_block = \
                (entity_to_string(target_entry.get_biopython_structure())
                 if isinstance(target_entry, MolFileEntry) else None)

            # Load PDB and extract hetatm.
            self.load_pdb(pdb_file, main_grp, mol_block)

            residue_selections = set()
            for index, shell in enumerate(shells):
                sphere_obj = ("%s.spheres.s%d_lvl%d_center%d"
                              % (main_grp, index, shell.level,
                                 hash(shell.central_atm_grp)))
                centroid_obj = "%s.centroid" % (sphere_obj)

                centroids = list(shell.central_atm_grp.centroid)
                self.wrapper.add_pseudoatom(centroid_obj,
                                            {"color": "white",
                                             "pos": centroids})
                self.wrapper.hide([("nonbonded", centroid_obj)])
                self.wrapper.show([("spheres", centroid_obj),
                                   ("nb_spheres", centroid_obj),
                                   ("dots", centroid_obj)])
                self.wrapper.set("dot_color", "red")
                self.wrapper.set("sphere_scale",
                                 shell.radius,
                                 {"selection": centroid_obj})
                self.wrapper.set("sphere_transparency", 0.85,
                                 {"selection": centroid_obj})
                self.wrapper.run_cmds([("center",
                                        {"selection": centroid_obj})])

                atm_names = [a.name
                             for a in sorted(shell.central_atm_grp.atoms)]
                self.wrapper.label([(centroid_obj,
                                     '"%s"' % "+".join(atm_names))])

                # Add interactions and styles.
                interacting_residue_sels = \
                    self.set_interactions_view(shell.interactions,
                                               main_grp, sphere_obj)

                residue_selections.update(interacting_residue_sels)

            if residue_selections:
                self.wrapper.select(name="%s.inter_residues" % main_grp,
                                    selection=" or ".join(residue_selections))

        self.set_last_details_to_view()
