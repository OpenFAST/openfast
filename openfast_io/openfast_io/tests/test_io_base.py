from pathlib import Path
import pytest
from openfast_io.io.base import ModuleIO

from openfast_io.io.aerodyn import AeroDynIO
from openfast_io.io.elastodyn import ElastoDynIO
from openfast_io.io.simple_elastodyn import SimpleElastoDynIO
from openfast_io.io.beamdyn import BeamDynIO
from openfast_io.io.inflowwind import InflowWindIO
from openfast_io.io.aerodisk import AeroDiskIO
from openfast_io.io.servodyn import ServoDynIO
from openfast_io.io.hydrodyn import HydroDynIO
from openfast_io.io.seastate import SeaStateIO
from openfast_io.io.subdyn import SubDynIO
from openfast_io.io.moordyn import MoorDynIO
from openfast_io.io.map_io import MAPIO
from openfast_io.io.extptfm import ExtPtfmIO


# Every concrete module IO class — replaces the per-module
# *_is_module_io / *_implements_module_io boilerplate.
ALL_IO_CLASSES = [
    AeroDynIO, ElastoDynIO, SimpleElastoDynIO, BeamDynIO, InflowWindIO,
    AeroDiskIO, ServoDynIO, HydroDynIO, SeaStateIO, SubDynIO, MoorDynIO,
    MAPIO, ExtPtfmIO,
]


@pytest.mark.parametrize("cls", ALL_IO_CLASSES, ids=lambda c: c.__name__)
def test_io_class_is_module_io(cls):
    """Each concrete IO class instantiates and is a ModuleIO subclass."""
    assert isinstance(cls(), ModuleIO)


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
