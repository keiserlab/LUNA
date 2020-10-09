from enum import Enum, auto


class IFPType(Enum):
    # auto() creates automatic identifiers for each type as this value is not important.
    # Therefore, always remember to compare IFP types using their names and not their values.
    EIFP = auto()
    EIFP_WITH_PHARM_FOR_GROUPS = auto()
    FIFP = auto()
