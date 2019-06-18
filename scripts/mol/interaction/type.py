class InteractionType():

    def __init__(self, src_grp, trgt_grp, inter_type, directional=False, params=None, recursive=True):
        self._src_grp = src_grp
        self._trgt_grp = trgt_grp
        self._type = inter_type
        self.directional = directional
        self._params = params or {}
        self._recursive = recursive
        self._hash_cache = None

        if recursive:
            self.src_grp.add_interactions([self])
            self.trgt_grp.add_interactions([self])

        self._expand_dict()

    @property
    def src_grp(self):
        return self._src_grp

    @property
    def trgt_grp(self):
        return self._trgt_grp

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
        if self._recursive:
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
