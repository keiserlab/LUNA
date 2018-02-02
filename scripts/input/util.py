def entry_format(entries, sep=":"):
    return sep.join([str(x) for x in entries if x is not None])


def format_selector(selTuple):
    resname, resnum, icode = selTuple

    if (resname == " "):
        resname = None
    elif (resname == "W"):
        resname = "HOH"
    else:
        resname = resname[2:]

    return entry_format([resname, str(resnum) + icode])
