"""
Contains tools for testing and debugging
nsaba submodules.
"""


# Uses for methods that may be imported and accidentally used.
def not_operational(func):
    def func_wrap(*args, **kwargs):
        raise ImportError("'%s': is still in development and not operational." % func.func_name)
    return func_wrap


class preprint(object):
    def __init__(self, string):
        self.string = string

    def __call__(self, f):
        def f_wrapper(*args, **kwargs):
            print self.string
            return f(*args, **kwargs)
        return f_wrapper
