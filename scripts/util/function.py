def func_call_2str(func, *args, **kwargs):
    arg_names = func.__code__.co_varnames[:func.__code__.co_argcount]
    args = args[:len(arg_names)]
    defaults = func.__defaults__ or ()
    args = (args + defaults[len(defaults) -
            (func.__code__.co_argcount - len(args)):])
    params = zip(arg_names, args)
    args = args[len(arg_names):]
    if args:
        params.append(('args', args))
    if kwargs:
        params.append(('kwargs', kwargs))

    args_as_str = ', '.join('%s=%r' % p for p in params) + ' )'
    func_as_str = "%s(%s)" % (func.__name__, args_as_str)

    return func_as_str
