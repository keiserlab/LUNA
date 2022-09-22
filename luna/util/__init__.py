import string
import random
import inspect
import warnings

import logging
logger = logging.getLogger()


def new_random_string(size=32, chars=string.ascii_uppercase + string.digits):
    """Generate a new random string of size ``size`` containing only the
    characters provided in ``chars``.

    Parameters
    ----------
    size : int, optional
            The size of the new string. The default value is 32.
    chars : iterable, optional
            A sequence of characters to choose from. The default value is
            'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'.

    Returns
    -------
    random_string : str
            The new random string.
    """
    return ('').join((random.choice(chars) for i in range(size)))


def rgb2hex(r, g, b):
    return "#{:02x}{:02x}{:02x}".format(r, g, b)


def hex2rgb(hexcode):
    return tuple(map(ord, hexcode[1:].decode('hex')))


# TODO: Transform hex to rgb and vice versa
class ColorPallete:
    def __init__(self, color_map=None, default_color=(255, 255, 255)):
        self.color_map = color_map or {}
        self.default_color = default_color

    def add_color(self, key, color):
        self.color_map[key] = color

    def get_normalized_color(self, key):
        color = self.get_color(key)

        if isinstance(color, str):
            return color
        elif all([x < 1 for x in color]):
            return color
        return tuple(c / 255 for c in color)

    def get_unnormalized_color(self, key):
        color = self.get_color(key)

        if isinstance(color, str):
            return color
        elif all([x > 1 for x in color]):
            return color
        return tuple(int(c * 255) for c in color)

    def get_color(self, key):
        if key in self.color_map:
            return self.color_map[key]
        elif self.default_color:
            logger.warning("Key '%s' does not exist. Then, the default color "
                           "will be used: %s." % (key, self.default_color))
            return self.default_color
        else:
            raise KeyError("Key '%s' does not exist and no default "
                           "color was defined." % key)

    def __contains__(self, key):
        return key in self.color_map


def func_call_to_str(func, *args, **kwargs):
    arg_names = func.__code__.co_varnames[:func.__code__.co_argcount]
    args = args[:len(arg_names)]
    defaults = func.__defaults__ or ()
    args = (args + defaults[len(defaults) - (func.__code__.co_argcount - len(args)):])
    params = zip(arg_names, args)
    args = args[len(arg_names):]
    if args:
        params.append(('args', args))
    if kwargs:
        params.append(('kwargs', kwargs))

    args_as_str = ', '.join('%s=%r' % p for p in params) + ' )'
    func_as_str = "%s(%s)" % (func.__name__, args_as_str)

    return func_as_str


def iter_to_chunks(l, n):
    return [l[i:i + n] for i in range(0, len(l), n)]


class LUNAWarning(Warning):
    """Base LUNA warning class.

    Unlike normal warnings, these are by default always set to on.
    """


class LUNADeprecationWarning(LUNAWarning, DeprecationWarning):
    """A warning class for a deprecated method or class."""


# Always show custom warnings for this package
warnings.filterwarnings("always", category=LUNAWarning)


class deprecated(object):
    """Decorator to mark a function as deprecated.

    Issue a deprecation warning when a function is called, and update the
    documentation. A deprecation version must be provided.

    Examples
    --------
    >>> from luna.util import deprecated
    >>> @deprecated("1.1", remove_version="1.3",
    ...             msg="Function no longer needed")
    ... def my_function():
    ...     pass

    Notes
    -----
    Adapted from https://wiki.python.org/moin/PythonDecoratorLibrary
    """

    def __init__(self, deprecated_version, remove_version=None, msg=None):
        """Constructor.

        Parameters
        ----------
        deprecated_version : str
            Version in which object was deprecated (e.g. '1.1')
        remove_version : str, optional
            Version in which object will be removed (e.g. '1.2'). If not
            specified, it is assumed the object will be removed in the next
            release (e.g. '1.2' if `deprecated_version` is '1.1')
        msg : str, optional
            Message to include with deprecation warning, to explain deprecation
            or point to newer version.
        """
        self.deprecated_version = deprecated_version
        if remove_version is None:
            version_info = deprecated_version.split(".")
            version_info[1] = str(int(version_info[1]) + 1)
            for i in range(2, len(version_info)):
                version_info[i] = "0"
            remove_version = ".".join(version_info)
        self.remove_version = remove_version
        if msg is None:
            self.extra = ""
        else:
            self.extra = " {0}".format(msg)

    def __call__(self, obj):
        if inspect.isfunction(obj):
            return self.deprecate_function(obj)
        else:
            raise ValueError("Deprecated object is not a function.")

    def deprecate_function(self, f):
        """Return the decorated function."""
        msg = (
            "Function `{0}` was deprecated in {1} and will be removed "
            "in {2}.{3}"
        ).format(
            f.__name__,
            self.deprecated_version,
            self.remove_version,
            self.extra,
        )

        def new_func(*args, **kwargs):
            warnings.warn(msg, category=LUNADeprecationWarning, stacklevel=2)
            return f(*args, **kwargs)

        new_func.__name__ = f.__name__
        new_func.__dict__ = f.__dict__
        new_func.__doc__ = f.__doc__
        self.update_docstring(new_func)
        return new_func

    def update_docstring(self, obj):
        """Add deprecation note to docstring."""
        msg = (
            ".. note:: Deprecated in LUNA {0}.\n"
            "   `{1}` will be removed in LUNA {2}.{3}"
        ).format(
            self.deprecated_version,
            obj.__name__,
            self.remove_version,
            self.extra,
        )
        obj.__doc__ = "{0}\n\n{1}".format(msg, obj.__doc__)
        return obj
