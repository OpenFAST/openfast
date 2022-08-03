
import sys
import argparse
import numpy as np
from pathlib import Path
interface_path = Path(__file__).parent.parent.parent.parent / "glue-codes" / "python"
sys.path.insert(0, str(interface_path))
import openfast_library

def test_hub_position(library_path, input_file):
    openfastlib = openfast_library.FastLibAPI(library_path, input_file)
    openfastlib.fast_init()
    absolute_position, rotational_velocity, orientation_dcm = openfastlib.get_hub_position()

    # Initial hub position is at -5, 0, 90.55
    np.testing.assert_allclose(
        absolute_position,
        np.array([-5.0, 0.0, 90.55]),
        rtol=1e-5,
        atol=1e-8,
        verbose=True
    )

    # This case is initially still
    # Velocities should be 0 and the DCM should be identity
    np.testing.assert_array_equal( rotational_velocity, np.zeros(3) )
    np.testing.assert_array_equal( np.reshape( orientation_dcm[:], (3,3) ), np.eye(3) )


if __name__=="__main__":

    parser = argparse.ArgumentParser(description="Executes Python-based tests for OpenFAST Library.")
    parser.add_argument("input_file", metavar="Input-File", type=str, nargs=1, help="Path to an input file.")

    args = parser.parse_args()
    input_file = args.input_file[0]

    library_path = Path(__file__).parent.parent.parent.parent / "build" / "modules" / "openfast-library" / "libopenfastlib"
    if sys.platform == "linux" or sys.platform == "linux2":
        library_path = library_path.with_suffix(".so")
    elif sys.platform == "darwin":
        library_path = library_path.with_suffix(".dylib")
    elif sys.platform == "win32":
        # TODO
        pass

    test_hub_position(library_path, input_file)
