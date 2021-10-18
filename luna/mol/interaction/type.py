import numpy as np

from luna.mol.interaction.math import atom_coordinates, centroid
from luna.util.default_values import PYMOL_INTERACTION_COLOR_AS_RGB
from luna.util import rgb2hex


NUCLEOPHILE_INTERS = ["Orthogonal multipolar", "Parallel multipolar", "Antiparallel multipolar", "Tilted multipolar", "Multipolar",
                      "Cation-nucleophile", "Unfavorable anion-nucleophile", "Unfavorable nucleophile-nucleophile"]

ELECTROPHILE_INTERS = ["Orthogonal multipolar", "Parallel multipolar", "Antiparallel multipolar", "Tilted multipolar", "Multipolar",
                       "Anion-electrophile", "Unfavorable cation-electrophile", "Unfavorable electrophile-electrophile"]

UNFAVORABLE_INTERS = ["Repulsive", "Unfavorable anion-nucleophile", "Unfavorable cation-electrophile",
                      "Unfavorable nucleophile-nucleophile", "Unfavorable electrophile-electrophile"]


class InteractionType():

    def __init__(self, src_grp, trgt_grp, inter_type, src_interacting_atms=None, trgt_interacting_atms=None,
                 src_centroid=None, trgt_centroid=None, directional=False, params=None):

        self._src_grp = src_grp
        self._trgt_grp = trgt_grp

        src_interacting_atms = src_interacting_atms or []
        self._src_interacting_atms = list(src_interacting_atms)

        trgt_interacting_atms = trgt_interacting_atms or []
        self._trgt_interacting_atms = list(trgt_interacting_atms)

        self._src_centroid = np.array(src_centroid) if src_centroid is not None else None
        self._trgt_centroid = np.array(trgt_centroid) if trgt_centroid is not None else None

        self._type = inter_type
        self.directional = directional
        self._params = params or {}
        self._hash_cache = None

        self.apply_refs()
        self._expand_dict()

    @property
    def src_grp(self):
        return self._src_grp

    @src_grp.setter
    def src_grp(self, atm_grp):
        self._src_grp = atm_grp

        # Reset hash.
        self._hash_cache = None
        self.apply_refs()

    @property
    def trgt_grp(self):
        return self._trgt_grp

    @trgt_grp.setter
    def trgt_grp(self, atm_grp):
        self._trgt_grp = atm_grp

        # Reset hash.
        self._hash_cache = None
        self.apply_refs()

    @property
    def src_interacting_atms(self):
        return self._src_interacting_atms or self.src_grp.atoms

    @property
    def trgt_interacting_atms(self):
        return self._trgt_interacting_atms or self.trgt_grp.atoms

    @property
    def src_centroid(self):
        if self._src_centroid is None:
            if self._src_interacting_atms:
                self._src_centroid = centroid(atom_coordinates(self._src_interacting_atms))
            else:
                src_centroid = self._src_grp.centroid
                # Define the centroid in a nucleophile with two atoms as the position of its more electronegative atom.
                # Remember that the position in the interaction object matters. We have defined that the first group is always
                # the nucleophile for both dipole-dipole and ion-dipole interactions.
                if self.type in NUCLEOPHILE_INTERS and len(self.src_grp.atoms) == 2:
                    dipole_atm = self.src_grp.atoms[0] if (self.src_grp.atoms[0].electronegativity
                                                           > self.src_grp.atoms[1].electronegativity) else self.src_grp.atoms[1]
                    src_centroid = dipole_atm.coord
                # For unfavorable multipolar interactions, it may happen that the first atom group is an electrophile as well.
                elif self.type == "Unfavorable electrophile-electrophile" and len(self.src_grp.atoms) == 2:
                    dipole_atm = self.src_grp.atoms[0] if (self.src_grp.atoms[0].electronegativity
                                                           < self.src_grp.atoms[1].electronegativity) else self.src_grp.atoms[1]
                    src_centroid = dipole_atm.coord

                self._src_centroid = src_centroid
        return self._src_centroid

    @src_centroid.setter
    def src_centroid(self, centroid):
        if centroid is None:
            self._src_centroid = None
        else:
            self._src_centroid = np.array(centroid)

    @property
    def trgt_centroid(self):
        if self._trgt_centroid is None:
            if self._trgt_interacting_atms:
                self._trgt_centroid = centroid(atom_coordinates(self._trgt_interacting_atms))
            else:
                trgt_centroid = self._trgt_grp.centroid

                # Define the centroid in an electrophile with two atoms as the position of its less electronegative atom.
                # Remember that the position in the interaction object matters. We have defined that the second group is always
                # the electrophile for both dipole-dipole and ion-dipole interactions.
                if self.type in ELECTROPHILE_INTERS and len(self.trgt_grp.atoms) == 2:
                    dipole_atm = self.trgt_grp.atoms[0] if (self.trgt_grp.atoms[0].electronegativity
                                                            < self.trgt_grp.atoms[1].electronegativity) else self.trgt_grp.atoms[1]
                    trgt_centroid = dipole_atm.coord
                # For unfavorable multipolar interactions, it may happen that the second atom group is a nucleophile as well.
                elif self.type == "Unfavorable nucleophile-nucleophile" and len(self.trgt_grp.atoms) == 2:
                    dipole_atm = self.trgt_grp.atoms[0] if (self.trgt_grp.atoms[0].electronegativity
                                                            > self.trgt_grp.atoms[1].electronegativity) else self.trgt_grp.atoms[1]
                    trgt_centroid = dipole_atm.coord

                self._trgt_centroid = trgt_centroid
        return self._trgt_centroid

    @trgt_centroid.setter
    def trgt_centroid(self, centroid):
        if centroid is None:
            self._trgt_centroid = None
        else:
            self._trgt_centroid = np.array(centroid)

    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, new_type):
        self._type = new_type

        # Reset hash.
        self._hash_cache = None

    @property
    def params(self):
        return self._params

    @property
    def required_interactions(self):
        interactions = []
        if "depends_on" in self._params:
            interactions = self._params["depends_on"]

        return interactions

    def get_partner(self, comp):
        if comp == self.src_grp:
            return self.trgt_grp
        elif comp == self.trgt_grp:
            return self.src_grp
        return None

    def is_directional(self):
        return self.directional

    def is_intramol_interaction(self):
        comps1 = self.src_grp.compounds
        comps2 = self.trgt_grp.compounds
        return len(comps1) == 1 and len(comps2) == 1 and comps1 == comps2

    def is_intermol_interaction(self):
        return not self.is_intramol_interaction()

    def show_src_centroid(self):
        show_centroid = True

        # Define the centroid in a nucleophile with two atoms as the position of its more electronegative atom.
        # Remember that the position in the interaction object matters. We have defined that the first group is always
        # the nucleophile for both dipole-dipole and ion-dipole interactions.
        if self.type in NUCLEOPHILE_INTERS and len(self.src_grp.atoms) == 2:
            show_centroid = False
        # For unfavorable multipolar interactions, it may happen that the first atom group is an electrophile as well.
        elif self.type == "Unfavorable electrophile-electrophile" and len(self.src_grp.atoms) == 2:
            show_centroid = False

        return show_centroid

    def show_trgt_centroid(self):
        show_centroid = True

        # Define the centroid in an electrophile with two atoms as the position of its less electronegative atom.
        # Remember that the position in the interaction object matters. We have defined that the second group is always
        # the electrophile for both dipole-dipole and ion-dipole interactions.
        if self.type in ELECTROPHILE_INTERS and len(self.trgt_grp.atoms) == 2:
            show_centroid = False
        # For unfavorable multipolar interactions, it may happen that the second atom group is a nucleophile as well.
        elif self.type == "Unfavorable nucleophile-nucleophile" and len(self.trgt_grp.atoms) == 2:
            show_centroid = False

        return show_centroid

    def apply_refs(self):
        self.src_grp.add_interactions([self])
        self.trgt_grp.add_interactions([self])

    def clear_refs(self):
        self.src_grp.remove_interactions([self])
        self.trgt_grp.remove_interactions([self])

    def _expand_dict(self):
        for key in self._params:
            self.__dict__[key] = self._params[key]

    def as_json(self):
        inter_obj = {}
        inter_obj["type"] = self.type

        inter_obj["color"] = rgb2hex(*PYMOL_INTERACTION_COLOR_AS_RGB.get_unnormalized_color(self.type))
        inter_obj["is_directional"] = self.is_directional()
        inter_obj["is_intramol_interaction"] = self.is_intramol_interaction()
        inter_obj["is_intermol_interaction"] = self.is_intermol_interaction()

        src_grp_obj = self.src_grp.as_json()
        src_grp_obj["add_pseudo_group"] = len(self.src_grp.atoms) > 1
        src_grp_obj["centroid"] = self.src_centroid.tolist()
        src_grp_obj["show_centroid"] = self.show_src_centroid()

        if src_grp_obj["add_pseudo_group"]:
            src_grp_obj["pseudo_group_name"] = "+".join([a.name for a in sorted(self.src_grp.atoms)])

        inter_obj["src_grp"] = src_grp_obj

        #
        # Target atom group
        #
        trgt_grp_obj = self.trgt_grp.as_json()
        trgt_grp_obj["add_pseudo_group"] = len(self.trgt_grp.atoms) > 1
        trgt_grp_obj["centroid"] = self.trgt_centroid.tolist()
        trgt_grp_obj["show_centroid"] = self.show_trgt_centroid()

        if trgt_grp_obj["add_pseudo_group"]:
            trgt_grp_obj["pseudo_group_name"] = "+".join([a.name for a in sorted(self.trgt_grp.atoms)])

        inter_obj["trgt_grp"] = trgt_grp_obj

        return inter_obj

    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(self, other.__class__):
            is_equal_compounds = ((self.src_grp == other.src_grp and self.trgt_grp == other.trgt_grp)
                                  or (self.src_grp == other.trgt_grp and self.trgt_grp == other.src_grp))

            is_equal_interactions = self.type == other.type
            has_equal_params = self.params == other.params

            return is_equal_compounds and is_equal_interactions and has_equal_params
        return False

    def __ne__(self, other):
        """Overrides the default implementation"""
        return not self.__eq__(other)

    def __hash__(self):
        """Overrides the default implementation"""

        if self._hash_cache is None:
            # First, it flats the dictionary by transforming it to a list.
            # Then, it transforms the list into an immutable data structure (tuple).
            params_values = []
            for key in sorted(self.params):
                if type(self.params[key]) is list:
                    val = tuple(self.params[key])
                else:
                    val = self.params[key]
                params_values.append(val)
            params_as_tuple = tuple(params_values)

            # The properties src_grp and trgt_grp properties makes an InteractionType object order dependent.
            # For example, Class(X,Y) would be considered different from Class(Y,X).
            # However, in both cases the interactions should be considered the same.
            # Then, the next line turns the order dependent arguments into an independent order data.
            comp_values_as_tuple = tuple(sorted([self.src_grp, self.trgt_grp]))
            self._hash_cache = hash(tuple([comp_values_as_tuple, self.type, params_as_tuple]))

        return self._hash_cache

    def __repr__(self):
        return ('<InteractionType: compounds=(%s, %s) type=%s>' % (self.src_grp, self.trgt_grp, self.type))
