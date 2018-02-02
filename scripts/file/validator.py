def validate_filesystem(path, type):
    import os.path

    if (os.path.exists(path) is False):
        raise FileNotFoundError("File or directory '%s' does not exist" % path)

    if (type == "file" and os.path.isfile(path) is False):
        raise IsADirectoryError("'%s' is not a file" % path)        
    elif (type == "directory" and os.path.isdir(path) is False):        
        raise NotADirectoryError("'%s' is not a directory" % path)    

    return True


def is_directory_valid(path):
    return validate_filesystem(path, "directory")


def is_file_valid(path):
    return validate_filesystem(path, "file")
