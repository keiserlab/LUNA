import logging

logger = logging.getLogger()


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

    def get_color(self, key):
        if key in self.color_map:
            return self.color_map[key]
        elif self.default_color:
            logger.warning("Key '%s' does not exist. Then, the default color will be used: %s." % (key, self.default_color))
            return self.default_color
        else:
            raise KeyError("Key '%s' does not exist and no default color was defined." % key)

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
