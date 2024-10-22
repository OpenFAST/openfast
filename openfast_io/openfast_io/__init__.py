try:
    from ._version import __version__, __version_tuple__

except ImportError:
    __version__ = "undefined"
    __version_tuple__ = None  # type: ignore