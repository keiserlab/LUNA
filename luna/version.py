"""
    Pattern: MAJOR.MINOR.PATCH as in https://semver.org/, where
       - MAJOR stands for big scientific updates;

       - MINOR stands for updates that "break" the compatibility
            between pharmacophoric perception, interactions,
            and interaction fingerprints;

       - PATCH stands for backward compatible bug fixes.
"""
version_info = (0, 13, 0)
version = '.'.join(str(c) for c in version_info)
__version__ = version


def has_version_compatibility(v):
    version_to_check = v
    if isinstance(v, str):
        version_to_check = [int(i) for i in v.split(".")]
    try:
        if tuple(version_to_check[0:2]) == version_info[0:2]:
            return True
        else:
            return False
    except TypeError:
        raise
