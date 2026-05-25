from pathlib import Path
import pytest
from openfast_io.io.base import ModuleIO


def test_module_io_is_abstract():
    """ModuleIO cannot be instantiated directly."""
    with pytest.raises(TypeError):
        ModuleIO()


def test_module_io_subclass_must_implement_read():
    class NoRead(ModuleIO):
        def write(self, data, file_path, base_dir):
            pass
    with pytest.raises(TypeError):
        NoRead()


def test_module_io_subclass_must_implement_write():
    class NoWrite(ModuleIO):
        def read(self, file_path, base_dir):
            return {}
    with pytest.raises(TypeError):
        NoWrite()


def test_concrete_subclass_is_instantiable():
    class ConcreteIO(ModuleIO):
        def read(self, file_path: Path, base_dir: Path) -> dict:
            return {}
        def write(self, data: dict, file_path: Path, base_dir: Path) -> None:
            pass
    io = ConcreteIO()
    assert io.read(Path("."), Path(".")) == {}
