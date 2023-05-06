class LUNAError(Exception):
    """Base class for LUNA-specific errors.

       This class is provided for future LUNA-specific functionality.
    """
    pass


class CompatibilityError(RuntimeError):
    pass


class EntityLevelError(KeyError):
    pass


class MissingAtomsError(ValueError):
    pass


class AtomsListError(ValueError):
    pass


class ChainNotFoundError(KeyError):
    pass


class MoleculeNotFoundError(KeyError):
    pass


class IllegalArgumentError(ValueError):
    pass


class PDBNotReadError(IOError):
    pass


class PKLNotReadError(IOError):
    pass


class InvalidSuperpositionFileError(ValueError):
    pass


class FileNotCreated(IOError):
    pass


class FingerprintNotCreated(RuntimeError):
    pass


class InvalidFileFormat(ValueError):
    pass


class InvalidEntry(ValueError):
    pass


class DirectoryAlreadyExists(OSError):
    pass


class DuplicateEntry(ValueError):
    pass


class ProcessingFailed(RuntimeError):
    pass


class ShellCenterNotFound(KeyError):
    pass


class PymolSessionNotInitialized(NameError):
    pass


class InvalidFingerprintType(TypeError):
    pass


class BitsValueError(ValueError):
    pass


class FingerprintCountsError(ValueError):
    pass


class MoleculeSizeError(ValueError):
    pass


class SanitizationError(ValueError):
    pass


class AtomObjectTypeError(TypeError):
    pass


class BondObjectTypeError(TypeError):
    pass


class MoleculeObjectTypeError(TypeError):
    pass


class MoleculeObjectError(RuntimeError):
    pass


class LicenseNotFoundError(RuntimeError):
    pass
