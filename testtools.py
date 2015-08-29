"""
Contains tools for testing and debugging
nsaba submodules.
"""


# Uses for methods that may be imported and accidentally used.
def not_operational(func):
    def func_wrap(*args, **kwargs):
        raise ImportError("'%s': is still in development and not operational." % func.func_name)
    return func_wrap
