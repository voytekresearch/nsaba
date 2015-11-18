"""
nsabatools.py
Contains tools, classes and data structures for nsaba submodules.
"""


class ReadOnlyDict(dict):
    def __init__(self, *args, **kwargs):
        super(ReadOnlyDict, self).__init__(*args, **kwargs)

    def __setitem__(self, key, value):
        if key in self.keys():
            raise UserWarning("Read Only Dictionary")
        else:
            super(ReadOnlyDict, self).__setitem__(key, value)


def not_operational(func):
    """Decorater: Used for methods that may be imported and accidentally used."""
    def func_wrap(*args, **kwargs):
        raise ImportError("'%s': is still in development and not operational." % func.func_name)
    return func_wrap


class preprint(object):
    """Decorater: Prints string prior to executing function."""
    def __init__(self, string):
        self.string = string

    def __call__(self, f):
        def f_wrapper(*args, **kwargs):
            print self.string
            return f(*args, **kwargs)
        return f_wrapper
