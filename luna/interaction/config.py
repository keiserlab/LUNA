import logging
logger = logging.getLogger()


class InteractionConfig(dict):
    """Generic class to define parameters for interactions.

    Parameters
    ----------
    config : dict, optional
        A dict containing parameters for calculating interactions.

    Examples
    --------

    >>> config = {"max_ha_dist_hb_inter": 2.5}
    >>> inter_config = InteractionConfig(config)
    >>> print(inter_config["max_ha_dist_hb_inter"])
    2.5
    inter_config["max_ha_dist_hb_inter"] = 3
    >>> print(inter_config["max_ha_dist_hb_inter"])
    3
    """

    def __init__(self, config=None):
        config = config or {}
        super().__init__(config)

    @property
    def params(self):
        """list: The list of parameters."""
        return sorted(self.keys())


class DefaultInteractionConfig(InteractionConfig):
    """Default parameters for calculating interactions in LUNA."""

    def __init__(self):
        config = {}

        # Hydrogen bond
        config["max_da_dist_hb_inter"] = 3.9
        config["max_ha_dist_hb_inter"] = 2.8
        config["min_dha_ang_hb_inter"] = 90
        config["min_har_ang_hb_inter"] = 90
        config["min_dar_ang_hb_inter"] = 90

        # Weak hydrogen bond
        # Ref: Panigrahi, S. K. & Desiraju, G. R. (2007).
        # Ref: Desiraju, G. R. & Steiner, T. (2001).
        config["max_da_dist_whb_inter"] = 4
        config["max_ha_dist_whb_inter"] = 3
        config["min_dha_ang_whb_inter"] = 110
        config["min_har_ang_whb_inter"] = 90
        config["min_dar_ang_whb_inter"] = 90

        # Weak hydrogen bonds: hydrogen bonds involving aromatic rings
        # Ref: Hydrogen bonds with π-acceptors in proteins: frequencies and role in stabilizing local 3D structures [Steiner, 2001]
        # Ref: Strong and Weak Hydrogen Bonds [Panigrahi, 2007]
        config["max_dc_dist_whb_inter"] = 4.5
        config["max_hc_dist_whb_inter"] = 3.5
        config["min_dhc_ang_whb_inter"] = 120
        config["max_disp_ang_whb_inter"] = 40

        # Ionic interactions
        # Ref: Barlow, D. J.; Thornton, J. M. Ion-Pairs in Proteins. Journal of Molecular Biology 1983.
        config["max_dist_repuls_inter"] = 6
        config["max_dist_attract_inter"] = 6

        # Aromatic stacking:
        # Ref: Bhattacharyya et al. 2003. Geometry of Interaction of the Histidine Ring with Other Planar and Basic Residues.
        # Ref: Chakrabarti, P.; Bhattacharyya, R. 2007. Geometry of Nonbonded Interactions Involving Planar Groups in Proteins.
        config["max_cc_dist_pi_pi_inter"] = 6
        # Define if a ring is tilted, same plane or in a T-format
        config["min_dihed_ang_slope_pi_pi_inter"] = 30
        config["max_dihed_ang_slope_pi_pi_inter"] = 60
        # Define if the ring is offset or aligned with another ring center.
        config["min_disp_ang_offset_pi_pi_inter"] = 30
        config["max_disp_ang_offset_pi_pi_inter"] = 60

        # Amide-aromatic stacking
        # [1] A systematic analysis of atomic protein–ligand interactions in the PDB [Freitas, 2007].
        # [2] Efficient Stacking on Protein Amide Fragments [Harder, 2013].
        # [3] The environment of amide groups in protein–ligand complexes: H-bonds and beyond [Cotesta, 2006].
        # [4] Hydrogen bonds with π-acceptors in proteins: frequencies and role in stabilizing local 3D structures [Steiner, 2001]
        # [5] Example of amide-pi: Systematic Investigation of Halogen Bonding in Protein–Ligand Interactions [Hardegger, 2011]
        #
        # Ref: [1], [2], [4], and [5].
        config["max_cc_dist_amide_pi_inter"] = 4.5
        # Ref: [1] and [3].
        config["max_dihed_ang_amide_pi_inter"] = 30
        # Ref: [3]: I use the centroid of an amide, while they use the nitrogen.
        config["max_disp_ang_pi_pi_inter"] = 30

        # Hydrophobic interaction
        config["max_dist_hydrop_inter"] = 4.5
        config["min_surf_size"] = 1
        config["min_inter_atom_in_surf"] = 1

        # Cation-pi interaction
        # Ref: Gallivan, J. P.; Dougherty, D. A. Cation-π Interactions in Structural Biology. Proc Natl Acad Sci U S A 1999, 96 (17), 9459–9464.
        config["max_dist_cation_pi_inter"] = 6

        # Halogen bond.
        # Interaction model: C-X ---- A-R, where C is a carbon, X a halogen, A an acceptor and R is an atom bonded to A.
        # Ref: Halogen bonds in biological molecules [Auffinger, 2004]
        # Ref: The Important Role of Halogen Bond in Substrate Selectivity of Enzymatic Catalysis [Jiang, 2016]
        #
        # Distance X-A when A is a single atom.
        config["max_xa_dist_xbond_inter"] = 4
        # Distance X-A when A is an aromatic ring, so C stands for Centroid.
        config["max_xc_dist_xbond_inter"] = 4.5
        config["min_cxa_ang_xbond_inter"] = 120
        config["min_xar_ang_xbond_inter"] = 80
        config["max_disp_ang_xbond_inter"] = 60

        # Chalcogen bond.
        # Interaction model: R-Y ---- A-N,
        # Where R is a carbon/sulfur, Y a chalcogen (S, Se, Te), A an acceptor, and N a covalently bonded atom to A.
        # Ref: Mining and Structural Characterization of S···X Chalcogen Bonds in Protein Database [Iwaoka, Michio, and Natsuki Babe, 2015].
        # Ref: Chalcogen Bonding ‘2S–2N Squares’ versus Competing Interactions: Exploring the Recognition Properties of Sulfur [Ams et al, 2018].
        # Ref: S···O and S···N Sulfur Bonding Interactions in Protein–Ligand Complexes:
        #       Empirical Considerations and Scoring Function [Koebel, 2016].
        #
        # Distance Y-A when A is a single atom.
        config["max_ya_dist_ybond_inter"] = 4
        # Distance Y-A, when A (acceptor) is an aromatic ring, so C stands for Centroid.
        config["max_yc_dist_ybond_inter"] = 4.5
        # Angles
        config["min_rya_ang_ybond_inter"] = 120
        config["min_yan_ang_ybond_inter"] = 80
        # When A (acceptor) is an aromatic ring.
        config["max_disp_ang_ybond_inter"] = 60

        # Orthogonal multipolar interaction (dipole-dipole)
        # All rules were obtained departing from Paulini R, Müller K, Diederich F. 2005.
        #       Orthogonal multipolar interactions in structural chemistry and biology.
        # Model: A-N ... E-Y, where N and E are the nucleophilic and electrophilic atom, respectively.
        # Distance between the nucleophilic and electrophilic atom.
        config["max_ne_dist_multipolar_inter"] = 4
        # The angle NEY have a min and a max angle.
        config["min_ney_ang_multipolar_inter"] = 70
        config["max_ney_ang_multipolar_inter"] = 110
        config["max_disp_ang_multipolar_inter"] = 40
        # Orthogonal multipolar
        config["min_an_ey_ang_ortho_multipolar_inter"] = 70
        config["max_an_ey_ang_ortho_multipolar_inter"] = 110
        # Parallel and antiparallel multipolar
        config["max_an_ey_ang_para_multipolar_inter"] = 25
        config["min_an_ey_ang_antipara_multipolar_inter"] = 155

        # Ion-multipole interaction (ion-dipole)
        # Model: I ... D-Y, where I is the ion, D the dipole atom of interest (the electrophile or nucleophile) and Y is its counterpart.
        config["max_id_dist_ion_multipole_inter"] = 4.5
        config["min_idy_ang_ion_multipole_inter"] = 60
        config["max_disp_ang_ion_multipole_inter"] = 40

        # Proximal interactions
        config["max_dist_proximal"] = 6
        config["min_dist_proximal"] = 2

        # Van der Waals interactions. From: Arpeggio.
        config["vdw_tolerance"] = 0.1

        # From Chimera: https://www.cgl.ucsf.edu/chimerax/docs/user/commands/clashes.html
        config["vdw_clash_tolerance"] = 0.6

        # From Chimera: https://www.cgl.ucsf.edu/chimerax/docs/user/commands/clashes.html
        # Controls all interactions, except clashes that has its own parameter.
        config["min_bond_separation"] = 3
        # Controls only clashes.
        config["min_bond_separation_for_clash"] = 4

        config["boundary_cutoff"] = 6.2

        super().__init__(config)
