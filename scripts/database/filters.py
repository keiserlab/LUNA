class FilterRules:

    def __init__(self, label=None, rules=None):
        self.label = label

        if (rules is None):
            rules = []

        self.rules = rules

    def add_rule(self, rule):
        self.rules.append(rule)
