"""
Functions to convert strings.
"""

import re


def camelcase(string):
    """ Convert string into camel case.

    Parameters
    ----------
        string : str
            String to convert.

    Returns
    -------
        string : str
            Camel case string.
    """

    string = re.sub(r"\w[\s\W]+\w", '', str(string))
    if not string:
        return string
    return lowercase(string[0]) + re.sub(r"[\-_\.\s]([a-z])", lambda matched: uppercase(matched.group(1)), string[1:])


def capitalcase(string):
    """Convert string into capital case.
    First letters will be uppercase.

    Parameters
    ----------
        string : str
            String to convert.

    Returns
    -------
        string : str
            Capital case string.

    Examples
    --------
    >>> import luna.util.stringcase as sc
    >>> new_str = sc.camelcase("This is an interesting example.")
    >>> print(new_str)
    'this is an interesting example.'

    """
    string = str(string)
    if not string:
        return string
    return uppercase(string[0]) + string[1:]


def constcase(string):
    """Convert string into upper snake case.
    Join punctuation with underscore and convert letters into uppercase.

    Parameters
    ----------
        string : str
            String to convert.

    Returns
    -------
        string : str
            Const cased string.
    """

    return uppercase(snakecase(string))


def lowercase(string):
    """Convert string into lower case.

    Parameters
    ----------
        string : str
            String to convert.

    Returns
    -------
        string : str
            Lowercase case string.

    Examples
    --------
    >>> import luna.util.stringcase as sc
    >>> new_str = sc.camelcase("This is an Interesting EXAMPLE.")
    >>> print(new_str)
    'this is an interesting example.'
    """

    return str(string).lower()


def pascalcase(string):
    """Convert string into pascal case.

    Parameters
    ----------
        string : str
            String to convert.

    Returns
    -------
        string : str
            Pascal case string.
    """

    return capitalcase(camelcase(string))


def pathcase(string):
    """Convert string into path case.
    Join punctuation with slash.

    Parameters
    ----------
        string : str
            String to convert.

    Returns
    -------
        string : str
            Path cased string.
    """
    string = snakecase(string)
    if not string:
        return string
    return re.sub(r"_", "/", string)


def backslashcase(string):
    """Convert string into spinal case.
    Join punctuation with backslash.

    Parameters
    ----------
        string : str
            String to convert.

    Returns
    -------
        string : str
            Spinal cased string.
    """
    str1 = re.sub(r"_", r"\\", snakecase(string))

    return str1
    # return re.sub(r"\\n", "", str1))  # TODO: make regex fot \t ...


def sentencecase(string):
    """Convert string into sentence case.
    First letter capped and each punctuations are joined with space.

    Parameters
    ----------
        string : str
            String to convert.

    Returns
    -------
        string : str
            Sentence cased string.
    """
    joiner = ' '
    string = re.sub(r"[\-_\.\s]", joiner, str(string))
    if not string:
        return string
    return capitalcase(trimcase(
        re.sub(r"[A-Z]", lambda matched: joiner +
               lowercase(matched.group(0)), string)
    ))


def snakecase(string):
    """Convert string into snake case.
    Join punctuation with underscore.

    Parameters
    ----------
        string : str
            String to convert.

    Returns
    -------
        string : str
            Snake cased string.
    """

    string = re.sub(r"[\-\.\s]", '_', str(string))
    if not string:
        return string
    return lowercase(string[0]) + re.sub(r"[A-Z]", lambda matched: '_' + lowercase(matched.group(0)), string[1:])


def spinalcase(string):
    """Convert string into spinal case.
    Join punctuation with hyphen.

    Parameters
    ----------
        string : str
            String to convert.

    Returns
    -------
        string : str
            Spinal cased string.
    """

    return re.sub(r"_", "-", snakecase(string))


def dotcase(string):
    """Convert string into dot case.
    Join punctuation with dot.

    Parameters
    ----------
        string : str
            String to convert.

    Returns
    -------
        string : str
            Dot cased string.
    """

    return re.sub(r"_", ".", snakecase(string))


def titlecase(string):
    """Convert string into sentence case.
    First letter capped while each punctuations is capitalsed
    and joined with space.

    Parameters
    ----------
        string : str
            String to convert.

    Returns
    -------
        string : str
            Title cased string.
    """

    return ' '.join(
        [capitalcase(word) for word in snakecase(string).split("_")]
    )


def trimcase(string):
    """Convert string into trimmed string.

    Parameters
    ----------
        string : str
            String to convert.

    Returns
    -------
        string : str
            Trimmed case string.

    Examples
    --------
    >>> import luna.util.stringcase as sc
    >>> new_str = sc.camelcase("      This is an Interesting EXAMPLE.     ")
    >>> print(new_str)
    'This is an Interesting EXAMPLE.'
    """
    return str(string).strip()


def uppercase(string):
    """Convert string into upper case.

    Parameters
    ----------
        string : str
            String to convert.

    Returns
    -------
        string : str
            Uppercase case string.

    Examples
    --------
    >>> import luna.util.stringcase as sc
    >>> new_str = sc.camelcase("This is an Interesting EXAMPLE.")
    >>> print(new_str)
    'THIS IS AN INTERESTING EXAMPLE.'
    """
    return str(string).upper()


def alphanumcase(string):
    """Cuts all non-alphanumeric symbols,
    i.e. cuts all expect except 0-9, a-z and A-Z.

    Parameters
    ----------
        string : str
            String to convert.

    Returns
    -------
        string : str
            String with cutted non-alphanumeric symbols.
    """
    # return filter(str.isalnum, str(string))
    return re.sub("\W+", "", string)
