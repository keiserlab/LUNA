class Simple:

    status = 0
    description = None
    errors = []
    data = None

    def __init__(self, status=0, description=None, errors=[], data=None):
        self.status = status
        self.description = description
        self.errors = errors
        self.data = data
