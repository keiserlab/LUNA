from enum import Enum, auto


class IFPType(Enum):
    """An enumeration of Interaction FingerPrints (IFPs) available in LUNA."""

    # auto() creates automatic identifiers for each type as this value is not important.
    # Therefore, always remember to compare IFP types using their names and not their values.
    EIFP = auto()
    HIFP = auto()  # Hybrid version: atomic invariants for atoms and pharmacophore for groups
    FIFP = auto()
