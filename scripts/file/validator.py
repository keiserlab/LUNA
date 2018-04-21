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
    except Exception as e:
        return False


def is_file_valid(path):
    try:
        validate_filesystem(path, "file")
        return True
    except Exception as e:
        return False


def try_validate_directory(path):
    validate_filesystem(path, "directory")


def try_validate_file(path):
    validate_filesystem(path, "file")
