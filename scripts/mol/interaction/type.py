import numpy as np

from mol.interaction.math import atom_coordinates, centroid


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

        self._expand_dict()

    @property
    def src_grp(self):
        return self._src_grp

    @property
    def trgt_grp(self):
        return self._trgt_grp

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
                self._src_centroid = self._src_grp.centroid
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
                self._trgt_centroid = self._trgt_grp.centroid
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
        partner = None
        if comp == self.src_grp:
            partner = self.trgt_grp
        elif comp == self.trgt_grp:
            partner = self.src_grp

        return partner

    def is_directional(self):
        return self.directional

    def is_intramol_interaction(self):
        comps1 = self.src_grp.compounds
        comps2 = self.trgt_grp.compounds
        return len(comps1) == 1 and len(comps2) == 1 and comps1 == comps2

    def is_intermol_interaction(self):
        return not self.is_intramol_interaction()

    def clear_refs(self):
        self.src_grp.remove_interactions([self])
        self.trgt_grp.remove_interactions([self])

    def _expand_dict(self):
        for key in self._params:
            self.__dict__[key] = self._params[key]

    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(self, other.__class__):
            is_equal_compounds = ((self.src_grp == other.src_grp and self.trgt_grp == other.trgt_grp) or
                                  (self.src_grp == other.trgt_grp and self.trgt_grp == other.src_grp))

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
            comp_values_as_tuple = tuple(sorted([self.src_grp, self.trgt_grp], key=hash))
            self._hash_cache = hash(tuple([comp_values_as_tuple, self.type, params_as_tuple]))

        return self._hash_cache

    def __repr__(self):
        return ('<InteractionType: compounds=(%s, %s) type=%s>' % (self.src_grp, self.trgt_grp, self.type))
