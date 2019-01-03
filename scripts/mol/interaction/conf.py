import logging
logger = logging.getLogger(__name__)


class InteractionConf():

    def __init__(self, conf):
        self._conf = conf
        self._expand_dict()

    @property
    def conf(self):
        return self._conf

    @property
    def keys(self):
        return [k for k in self._conf]

    def add(self, key, val):
        if key not in self._conf:
            self._conf[key] = val
            self.__dict__[key] = val
        else:
            logger.info("Key '%s' already exists." % key)

    def alter(self, key, val):
        if key in self._conf:
            self.conf[key] = val
        else:
            logger.info("Key '%s' does not exist." % key)

    def get_value(self, key):
        if key in self._conf:
            return self.conf[key]
        else:
            logger.info("Key '%s' does not exist." % key)

    def _expand_dict(self):
        for key in self._conf:
            self.__dict__[key] = self._conf[key]

    def __getattr__(self, key):
        if key in self._conf:
            return self._conf[key]
        else:
            logger.info("Key '%s' does not exist." % key)
            return None


class DefaultInteractionConf(InteractionConf):

    def __init__(self):

        conf = {}

        # Hydrogen bond
        conf["max_da_dist_hb_inter"] = 3.9
        conf["max_ha_dist_hb_inter"] = 2.8
        conf["min_dha_ang_hb_inter"] = 90

        # Weak hydrogen bond
        # Ref: Panigrahi, S. K. & Desiraju, G. R. (2007).
        # Ref: Desiraju, G. R. & Steiner, T. (2001).
        conf["max_da_dist_whb_inter"] = 4
        conf["max_ha_dist_whb_inter"] = 3
        conf["min_dha_ang_whb_inter"] = 110
        conf["min_har_ang_whb_inter"] = 90
        conf["min_dar_ang_whb_inter"] = 90

        # Ionic interactions
        conf["max_dist_repuls_inter"] = 6
        conf["max_dist_attract_inter"] = 6

        # Aromatic stacking
        conf["max_cc_dist_pi_pi_inter"] = 6
        conf["min_dihed_ang_pi_pi_inter"] = 30
        conf["max_disp_ang_pi_pi_inter"] = 20

        # Hydrophobic interaction
        conf["max_dist_hydrop_inter"] = 4.5

        # Cation-pi interaction
        conf["max_dist_cation_pi_inter"] = 6

        # Halogen bond.
        # Interaction model: C-X ---- A-R,
        # Where C is a carbon, X a halogen, A an acceptor and
        # R is an atom bonded to A.
        # Distance X-A when A is an single atom.
        conf["max_xa_dist_xbond_inter"] = 4
        # Distance X-A when A is an aromatic ring, so C comes from Centroid.
        conf["max_xc_dist_xbond_inter"] = 4.5
        # Ref: Halogen bonds in biological molecules [Auffinger, 2004]
        conf["min_cxa_ang_xbond_inter"] = 120
        conf["min_xar_ang_xbond_inter"] = 80
        conf["max_disp_ang_xbond_inter"] = 60

        # Proximal interactions
        conf["max_dist_proximal"] = 6

        # Covalent interactions
        conf["cov_dist_tolerance"] = 0.45

        conf["boundary_cutoff"] = 7

        super().__init__(conf)
