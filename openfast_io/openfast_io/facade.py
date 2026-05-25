"""
Facade aliases -- re-exports from the canonical locations.

``InputReader_Facade`` and ``InputWriter_Facade`` are convenience aliases
that point to the same classes now living in ``FAST_reader.py`` and
``FAST_writer.py`` (which themselves are thin facades over the new IO layer).
"""
from openfast_io.FAST_reader import InputReader_OpenFAST as InputReader_Facade   # noqa: F401
from openfast_io.FAST_writer import InputWriter_OpenFAST as InputWriter_Facade   # noqa: F401

__all__ = ["InputReader_Facade", "InputWriter_Facade"]
