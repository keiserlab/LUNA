from os.path import (basename, splitext)


def get_file_format(filename):
    return splitext(basename(filename))[1][1:]
