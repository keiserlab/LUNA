from os.path import basename, exists, isdir, isfile, splitext
from os import makedirs, remove
from shutil import rmtree

import string
import random
import logging


logger = logging.getLogger()


def get_file_format(path, max_split=None):
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
    except ImportError as e:
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
    except ImportError as e:
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
            logger.warning("The directory '%s' already exists, and it will be cleaned before the program continues." % path)
            clear_directory(path)
        else:
            logger.warning("The directory '%s' already exists, but it will not be cleaned before the program continues." % path)
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
