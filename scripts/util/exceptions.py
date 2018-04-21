class EntityLevelError(KeyError):
    pass


class ChainNotFoundError(KeyError):
    pass


class MoleculeNotFoundError(KeyError):
    pass


class IllegalArgumentError(ValueError):
    pass


class PDBNotReadError(IOError):
    pass


class InvalidSuperpositionFileError(ValueError):
    pass


class FileNotCreated(IOError):
    pass


class FingerprintNotCreated(RuntimeError):
    pass


class InvalidFileFormat(ValueError):
    pass


class InvalidNapoliEntry(ValueError):
    pass


class DirectoryAlreadyExists(OSError):
    pass
