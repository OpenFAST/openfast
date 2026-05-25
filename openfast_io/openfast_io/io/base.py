from abc import ABC, abstractmethod
from pathlib import Path


class ModuleIO(ABC):
    """Base class for all module-level IO classes.

    Each subclass handles reading and writing one OpenFAST module's input files.
    read() returns a plain dict (the module's fst_vt slice).
    write() accepts that same dict and produces the file(s).

    base_dir is always passed explicitly — module IOs do not assume a working
    directory. This makes them testable in isolation.
    """

    @abstractmethod
    def read(self, file_path: Path, base_dir: Path) -> dict:
        """Read module input file(s). Returns dict matching the fst_vt module slice."""
        ...

    @abstractmethod
    def write(self, data: dict, file_path: Path, base_dir: Path) -> None:
        """Write module input file(s) from data dict."""
        ...
