from os.path import (basename, splitext)


def get_file_format(path):
    return generic_splitext(path)[1][1:]


def get_filename(path):
    return generic_splitext(path)[0]


def generic_splitext(path):
    filename, fileExt = splitext(basename(path))
    while '.' in filename:
        filename, tḿpFileExt = splitext(filename)
        fileExt = tḿpFileExt + fileExt

    return (filename, fileExt)


def generate_json_file(json_data, outputFile):
    import logging
    logger = logging.getLogger(__name__)

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


def clear_directory(path):
    import os, shutil
    
    if (os.path.isdir(path)):
        shutil.rmtree(path)
