from os.path import basename, exists, isdir, isfile, splitext
from os import makedirs, remove, listdir
from shutil import rmtree
import string
import logging
import pickle
import gzip

from luna.util.exceptions import FileNotCreated, PKLNotReadError
from luna.util import new_random_string


logger = logging.getLogger()


def detect_compression_format(filename):
    """
    Attempts to detect compression format from the filename extension.
    Returns None if no format could be detected.
    """
    if filename.endswith('.bz2'):
        return "bz2"
    elif filename.endswith('.xz'):
        return "xz"
    elif filename.endswith('.gz'):
        return "gz"
    else:
        return None


def get_file_format(path, max_split=None, ignore_compression=False):
    """Detect the file format from pathname ``path``.

    Parameters
    ----------
    path : str
        The pathname.
    max_split : int or None
        Specifies how many splits to do. The default value is None,
        which is "all occurrences".
    ignore_compression: bool, optional
        Ignore compression format. The default value is False,
        which does not ignore compression format.

    Returns
    -------
    filename : str
        The file format.
    """
    if ignore_compression is True:
        comp_format = detect_compression_format(path)
        if comp_format:
            path = path[:-(len(comp_format) + 1)]
    return generic_splitext(path, max_split)[1][1:]


def get_filename(path, max_split=None):
    """Detect the filename from pathname ``path``.

    Parameters
    ----------
    path : str
        The pathname.
    max_split : int or None
        Specifies how many splits to do. The default value is None,
        which is "all occurrences".

    Returns
    -------
    filename : str
        The filename.
    """
    return generic_splitext(path, max_split)[0]


def generic_splitext(path, max_split=None):
    """Split the pathname ``path`` into a pair (``filename``, ``ext``).

    Parameters
    ----------
    path : str
        The pathname.
    max_split : int or None
        Specifies how many splits to do. The default value is None,
        which is "all occurrences".

    Returns
    -------
    filename : str
        The filename.
    ext: str
        The file extension.
    """
    filename = basename(path)

    num_ext = filename.count('.') + 1
    if max_split is None or max_split < 1:
        max_split = num_ext
    else:
        max_split = min(max_split, num_ext)

    filename, ext = splitext(filename)
    max_split -= 1
    while max_split > 0:
        filename, tmp_ext = splitext(filename)
        ext = tmp_ext + ext
        max_split -= 1

    return (filename, ext)


def generate_json_file(json_data, output_file, indent=4, sort_keys=True):
    """Serialize ``json_data`` to a JSON formatted string and save it
    at ``output_file``.

    Parameters
    ----------
    json_data : object
        The data to be serialized.
    output_file : str
        The output file where the serialized data will be saved.
    indent : int
        If indent is a non-negative integer or string, then JSON array
        elements and object members will be pretty-printed with that indent
        level. An indent level of 0, negative, or "" will only insert newlines.
        None selects the most compact representation. Using a positive integer
        indent indents that many spaces per level. If indent is a string
        (such as "\t"), that string is used to indent each level.
        The default value is 4.

    sort_keys : bool
        If ``sort_keys`` is True, the output of dictionaries will be
        sorted by key.
    """
    try:
        import simplejson as json
        logger.warning("Module 'simplejson' imported.")
    except ImportError:
        logger.warning("Module 'simplejson' not available. "
                       "Built-in module 'json' will be imported.")
        import json

    try:
        with open(output_file, "w") as IN:
            json.dump(json_data, IN, indent=indent, sort_keys=sort_keys)
    except Exception as e:
        logger.exception(e)
        raise


def parse_json_file(json_file):
    """Deserialize the JSON file ``json_file``.

    Parameters
    ----------
    json_file : str
        The input JSON file.

    Returns
    -------
    json_data : object
        The deserialized data.
    """
    try:
        import simplejson as json
        logger.warning("Module 'simplejson' imported.")
    except ImportError:
        logger.warning("Module 'simplejson' not available. "
                       "Built-in module 'json' will be imported.")
        import json

    try:
        with open(json_file, "r") as IN:
            return json.load(IN)
    except Exception as e:
        logger.exception(e)
        raise


def create_directory(path, clear=False):
    """Create the directory pathaname ``path``.

    Parameters
    ----------
    path : str
        The directory pathname to be created.
    clear : bool
        If True, clear the directory if already exists.
        The default value is False.
    """
    try:
        if not exists(path):
            makedirs(path)
        elif clear:
            logger.info("The directory '%s' already exists, and it will be "
                        "cleared before the program continues." % path)
            remove_directory(path)
        else:
            logger.info("The directory '%s' already exists, but it will not "
                        "be cleared." % path)
    except OSError as e:
        logger.exception(e)
        raise


def remove_files(files):
    """Remove the provided files.

    Parameters
    ----------
    files : iterable
        An iterable object containing a sequence of files to be removed.
    """
    for f in files:
        if exists(f):
            remove(f)
        else:
            logger.info("File '%s' does not exist." % f)


def remove_directory(path, only_empty_paths=False):
    """Remove the directory given by the pathname ``path``.

    Parameters
    ----------
    path : str
        The pathname.
    only_empty_paths: bool, optional
        If True, do not remove non-empty directories.
        The default value is False, which removes all directories.
    """
    if isdir(path):
        # Do nothing if the directory is not empty and if only
        # empty paths must be removed.
        if only_empty_paths and len(listdir(path)) != 0:
            return

        try:
            rmtree(path, ignore_errors=True)
        except OSError:
            raise


def new_unique_filename(path,
                        size=32,
                        chars=string.ascii_uppercase + string.digits,
                        retries=5):
    """Generate a new unique random pathname.

    Parameters
    ----------
    path : str
        The target pathname.
    size : int, optional
        The size of the new filename. The default value is 32.
    chars : iterable, optional
        A sequence of characters to choose from.
        The default value is 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'.
    retries : int, optional
        The function will keep trying to generate a unique filename
        (not exist in the pathname ``path``) until the maximum number
        of retries ``retries`` is reached.

    Returns
    -------
    unique_pathname : str
        A new random unique pathname (path + random filename).
    """
    for r in range(retries):
        filename = '%s/%s' % (path, new_random_string(size, chars))
        if not exists(filename):
            return filename


def validate_filesystem(path, type):
    """Validate if a file or directory exists and if it has the expected type.

    Parameters
    ----------
    path : str
        The pathname of the file or directory to be validated.
    type : {'file', 'directory'}
        The filesystem type.

    Raises
    -------
    FileNotFoundError
        If the file or directory given by the pathname ``path`` does not exist.
    IsADirectoryError
        If a file is provided but a directory is found instead at the
        pathname ``path``.
    NotADirectoryError
        If the filesystem given by the pathname ``path`` is not a directory.
    """
    if exists(path) is False:
        raise FileNotFoundError("File or directory '%s' does not exist" % path)

    if type == "file" and isfile(path) is False:
        raise IsADirectoryError("'%s' is not a file." % path)
    elif type == "directory" and isdir(path) is False:
        raise NotADirectoryError("'%s' is not a directory." % path)


def is_directory_valid(path):
    """Check if ``path`` exists and if it is in fact a directory."""
    try:
        validate_filesystem(path, "directory")
        return True
    except Exception:
        return False


def is_file_valid(path):
    """Check if ``path`` exists and if it is in fact a file."""
    try:
        validate_filesystem(path, "file")
        return True
    except Exception:
        return False


def validate_directory(path):
    """Validate ``path`` as a directory."""
    validate_filesystem(path, "directory")


def validate_file(path):
    """Validate ``path`` as a file."""
    validate_filesystem(path, "file")


def pickle_data(data, output_file, compressed=True):
    """Write the pickled representation of the object ``data`` to
    the file ``output_file``.

    Parameters
    ----------
    data : object
        The object to be pickled.
    output_file : str
        The output file where the pickled representation will be saved.
    compressed : bool, optional
        If True (the default), compress the pickled representation as
        a gzip file (.gz).

    Raises
    -------
    FileNotCreated
        If the file could not be created.
    """
    open_func = open
    if compressed:
        open_func = gzip.open
        if output_file.endswith(".gz") is False:
            output_file += ".gz"

    try:
        with open_func(output_file, "wb") as OUT:
            pickle.dump(data, OUT, pickle.HIGHEST_PROTOCOL)
    except OSError as e:
        raise FileNotCreated("File '%s' could not be created."
                             % output_file) from e
    except Exception:
        raise


def unpickle_data(input_file):
    """Read the pickled representation of an object from the file
    ``input_file`` and return the reconstituted object hierarchy
    specified therein. ``input_file`` can be a gzip-compressed file.

    Raises
    -------
    PKLNotReadError
        If the file could not be loaded.
    """
    try:
        # Try to decompress and unpickle the data first.
        with gzip.open(input_file, "rb") as IN:
            return pickle.load(IN)
    except Exception:
        pass

    # If the decompression failed, let's try to unpickle the file directly.
    # Maybe it is not a compressed file.
    try:
        with open(input_file, "rb") as IN:
            return pickle.load(IN)
    except OSError as e:
        raise PKLNotReadError("File '%s' could not be loaded."
                              % input_file) from e
