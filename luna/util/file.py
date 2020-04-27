from os.path import basename, exists, isdir, isfile, splitext
from os import makedirs, remove
from shutil import rmtree
import string
import random
import logging
import pickle
import gzip

from luna.util.exceptions import FileNotCreated, PKLNotReadError


logger = logging.getLogger()


def detect_compression_format(filename):
    """
    Attempts to detect file format from the filename extension.
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
    if ignore_compression is True:
        comp_format = detect_compression_format(path)
        if comp_format:
            path = path[:-(len(comp_format) + 1)]
    return generic_splitext(path, max_split)[1][1:]


def get_filename(path, max_split=None):
    return generic_splitext(path, max_split)[0]


def generic_splitext(path, max_split=None):
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


def generate_json_file(json_data, output_file):
    try:
        import simplejson as json
        logger.info("Module 'simplejson' imported.")
    except ImportError:
        logger.info("Module 'simplejson' not available. Built-in module 'json' will be imported.")
        import json

    try:
        with open(output_file, "w") as IN:
            json.dump(json_data, IN, indent=4, sort_keys=True)
    except Exception as e:
        logger.exception(e)
        raise


def parse_json_file(json_file):
    try:
        import simplejson as json
        logger.info("Module 'simplejson' imported.")
    except ImportError:
        logger.info("Module 'simplejson' not available. Built-in module 'json' will be imported.")
        import json

    try:
        with open(json_file, "r") as IN:
            return json.load(IN)
    except Exception as e:
        logger.exception(e)
        raise


def create_directory(path, clear=False):
    try:
        if not exists(path):
            makedirs(path)
        elif clear:
            logger.info("The directory '%s' already exists, and it will be cleared before the program continues." % path)
            clear_directory(path)
        else:
            logger.info("The directory '%s' already exists, but it will not be cleared." % path)
    except OSError as e:
        logger.exception(e)
        raise


def remove_files(files):
    for f in files:
        if exists(f):
            remove(f)
        else:
            logger.info("File %s does not exist." % f)


def clear_directory(path):
    if isdir(path):
        try:
            rmtree(path, ignore_errors=True)
        except OSError as e:
            logger.exception(e)
            raise


def new_filename(size=32, chars=string.ascii_uppercase + string.digits):
    return ('').join((random.choice(chars) for i in range(size)))


def get_unique_filename(path, size=32, chars=string.ascii_uppercase + string.digits, retries=5):
    for r in range(retries):
        filename = '%s/%s' % (path, new_filename(size, chars))
        if not exists(filename):
            return filename
    return None


def validate_filesystem(path, type):
    if exists(path) is False:
        raise FileNotFoundError("File or directory '%s' does not exist" % path)

    if type == "file" and isfile(path) is False:
        raise IsADirectoryError("'%s' is not a file" % path)
    elif type == "directory" and isdir(path) is False:
        raise NotADirectoryError("'%s' is not a directory" % path)


def is_directory_valid(path):
    try:
        validate_filesystem(path, "directory")
        return True
    except Exception:
        return False


def is_file_valid(path):
    try:
        validate_filesystem(path, "file")
        return True
    except Exception:
        return False


def try_validate_directory(path):
    validate_filesystem(path, "directory")


def try_validate_file(path):
    validate_filesystem(path, "file")


def pickle_data(data, output_file, compressed=True):
    open_func = open
    if compressed:
        open_func = gzip.open
        if output_file.endswith(".gz") is False:
            output_file += ".gz"

    try:
        with open_func(output_file, "wb") as OUT:
            pickle.dump(data, OUT, pickle.HIGHEST_PROTOCOL)
    except OSError as e:
        logger.exception(e)
        raise FileNotCreated("File '%s' could not be created." % output_file)
    except Exception as e:
        logger.exception(e)
        raise


def unpickle_data(input_file):

    # Try to decompress and unpickle the data first.
    try:
        with gzip.open(input_file, "rb") as IN:
            return pickle.load(IN)
    except Exception:
        pass

    # If the decompression failed, let's try to unpickle the file directly. Maybe it is not a compressed file.
    try:
        with open(input_file, "rb") as IN:
            return pickle.load(IN)
    except OSError as e:
        logger.exception(e)
        raise PKLNotReadError("File '%s' could not be loaded." % input_file)
    return None
