from os.path import (basename, exists, isdir, splitext)
from os import makedirs
from shutil import rmtree

from util.exceptions import DirectoryAlreadyExists

import logging
logger = logging.getLogger(__name__)


def get_file_format(path, maxsplit=None):
    return generic_splitext(path, maxsplit)[1][1:]


def get_filename(path, maxsplit=None):
    return generic_splitext(path)[0]


def generic_splitext(path, maxsplit=None):
    filename = basename(path)

    numExt = filename.count('.') + 1
    if (maxsplit is None or maxsplit < 1):
        maxsplit = numExt
    else:
        maxsplit = min(maxsplit, numExt)

    filename, fileExt = splitext(filename)
    maxsplit -= 1
    while maxsplit > 0:
        filename, tmpFileExt = splitext(filename)
        fileExt = tmpFileExt + fileExt
        maxsplit -= 1

    return (filename, fileExt)


def generate_json_file(json_data, outputFile):
    try:
        import simplejson as json
        logger.info("Module 'simplejson' imported.")
    except ImportError as e:
        logger.exception(e)
        logger.warning("Module 'simplejson not available.")
        logger.warning("Built-in module 'json' will be imported.")
        import json

    try:
        with open(outputFile, "w") as f:
            json.dump(json_data, f, indent=4, sort_keys=True)
    except Exception as e:
        logger.exception(e)
        raise


def create_directory(path, clear=False):
    try:
        if not exists(path):
            makedirs(path)
        elif clear:
            clear_directory(path)
        else:
            raise DirectoryAlreadyExists("The directory '%s' already exists."
                                         % path)
    except OSError as e:
        logger.exception(e)
        raise


def clear_directory(path):
    if isdir(path):
        try:
            rmtree(path)
        except OSError as e:
            logger.exception(e)
            raise


def validate_filesystem(path, type):
    import os.path

    if os.path.exists(path) is False:
        raise FileNotFoundError("File or directory '%s' does not exist" % path)

    if type == "file" and os.path.isfile(path) is False:
        raise IsADirectoryError("'%s' is not a file" % path)
    elif type == "directory" and os.path.isdir(path) is False:
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
