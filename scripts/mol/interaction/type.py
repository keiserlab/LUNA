class InteractionType():

    def __init__(self, atm_grp1, atm_grp2, inter_type, params=None, recursive=True):
        self._atm_grp1 = atm_grp1
        self._atm_grp2 = atm_grp2
        self._type = inter_type
        self._params = params or {}

        self._hash_chache = None

        if recursive:
            self.atm_grp1.add_interaction(self)
            self.atm_grp2.add_interaction(self)
        self._recursive = recursive

        self._expand_dict()

    @property
    def atm_grp1(self):
        return self._atm_grp1

    @property
    def atm_grp2(self):
        return self._atm_grp2

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
        if comp == self.atm_grp1:
            partner = self.atm_grp2
        elif comp == self.atm_grp2:
            partner = self.atm_grp1

        return partner

    def clear_refs(self):
        if self._recursive:
            self.atm_grp1.remove_interaction(self)
            self.atm_grp2.remove_interaction(self)

    def _expand_dict(self):
        for key in self._params:
            self.__dict__[key] = self._params[key]

    def __eq__(self, other):
        """Overrides the default implementation"""
        is_equal = False
        if isinstance(self, other.__class__):
            is_equal_compounds = ((self.atm_grp1 == other.atm_grp1 and self.atm_grp2 == other.atm_grp2) or
                                  (self.atm_grp1 == other.atm_grp2 and self.atm_grp2 == other.atm_grp1))

            is_equal_interactions = self.type == other.type
            has_equal_params = self.params == other.params

            is_equal = is_equal_compounds and is_equal_interactions and has_equal_params

        return is_equal

    def __ne__(self, other):
        """Overrides the default implementation (unnecessary in Python 3)"""
        return not self.__eq__(other)

    def __hash__(self):
        """Overrides the default implementation"""

        if self._hash_chache is None:
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

            # The properties atm_grp1 and atm_grp2 properties makes an InteractionType object order dependent.
            # For example, Class(X,Y) would be considered different from Class(Y,X).
            # However, in both cases the interactions should be considered the same.
            # Then, the next line turns the order dependent arguments into an independent order data.
            comp_values_as_tuple = tuple(sorted([self.atm_grp1, self.atm_grp2], key=id))
            self._hash_chache = hash(tuple([comp_values_as_tuple, self.type, params_as_tuple]))

        return self._hash_chache

    def __repr__(self):
        return ('<InteractionType: compounds=(%s, %s) type=%s>' % (self.atm_grp1, self.atm_grp2, self.type))
