class EntryValidator:

    def __init__(self, pattern):
        import re
        self.regex = re.compile(pattern, flags=0)

    def is_valid(self, entry):
        return self.regex.match(entry) is not None


class PliEntryValidator(EntryValidator):

    def __init__(self, uploadFormat=False):
        pattern = '^\w{4}:\w:\w{1,3}:\-?\d+$'

        if (uploadFormat):
            pattern = '^\w+(\w|\-)*:\w:\w{1,3}:\-?\d+$'

        super(PliEntryValidator, self).__init__(pattern)


class PpiEntryValidator(EntryValidator):

    def __init__(self, uploadFormat=False):
        pattern = '^\w{4}:\w$'

        if (uploadFormat):
            pattern = '^\w+(\w|\-)*:\w$'

        super(PpiEntryValidator, self).__init__(pattern)
