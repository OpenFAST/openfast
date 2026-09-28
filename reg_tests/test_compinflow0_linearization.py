from pathlib import Path
import re
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
import rtestlib as rtl
from fast_linearization_file import FASTLinearizationFile


def header_count(header, label):
    for line in header:
        if label in line:
            return int(line.split(label, 1)[1])
    raise AssertionError(f"Missing {label!r} from linearization header")


def validate_linearization(path):
    linearization = FASTLinearizationFile(str(path))
    nx = header_count(linearization["header"], "Number of continuous states:")
    nu = header_count(linearization["header"], "Number of inputs:")
    ny = header_count(linearization["header"], "Number of outputs:")
    expected_shapes = {
        "A": (nx, nx),
        "B": (nx, nu),
        "C": (ny, nx),
        "D": (ny, nu),
    }
    for key, shape in expected_shapes.items():
        if shape[0] and shape[1]:
            assert key in linearization, f"{path}: missing {key} matrix"
            assert linearization[key].shape == shape, f"{path}: unexpected {key} matrix shape"
            assert np.isfinite(linearization[key]).all(), f"{path}: non-finite {key} values"
    for key in ("x", "xdot", "u", "y"):
        if key in linearization:
            assert np.isfinite(linearization[key]).all(), f"{path}: non-finite {key} operating point"
    for key, count in (("x", nx), ("u", nu), ("y", ny)):
        if count:
            assert key in linearization, f"{path}: missing {key} operating point"
            assert linearization[key].size == count, f"{path}: unexpected {key} operating point size"

    descriptions = linearization.get("u_info", {}).get("Description", [])
    assert len(descriptions) == nu, f"{path}: input channel count does not match header"
    return descriptions


def main(executable, fixture_root):
    executable = Path(executable).resolve()
    fixture_root = Path(fixture_root).resolve()
    case_name = "5MW_Land_Linear_Aero"
    with tempfile.TemporaryDirectory(prefix="openfast-compinflow0-") as temporary:
        run_root = Path(temporary)
        case_dir = run_root / case_name
        rtl.copyTree(str(fixture_root / case_name), str(case_dir))
        rtl.copyTree(str(fixture_root / "5MW_Baseline"), str(run_root / "5MW_Baseline"))
        input_path = case_dir / f"{case_name}.fst"
        input_text = input_path.read_text()
        input_text, replacements = re.subn(r"(?m)^(\s*)1(\s+CompInflow\b)", r"\g<1>0\2", input_text)
        assert replacements == 1, f"Expected one CompInflow switch in {input_path}"
        input_path.write_text(input_text)
        for path in case_dir.glob("*.lin"):
            path.unlink()

        result = subprocess.run(
            [str(executable), input_path.name],
            cwd=case_dir,
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 0, f"OpenFAST exited {result.returncode}\n{result.stdout}\n{result.stderr}"

        expected = {
            f"{case_name}.1.AD.lin",
            f"{case_name}.1.ED.lin",
            f"{case_name}.1.SrvD.lin",
            f"{case_name}.1.lin",
        }
        actual = {path.name for path in case_dir.glob("*.lin")}
        assert actual == expected, f"Expected {sorted(expected)}, found {sorted(actual)}"

        ad_file = case_dir / f"{case_name}.1.AD.lin"
        descriptions = validate_linearization(ad_file)
        for channel in (
            "Extended input: horizontal wind speed",
            "Extended input: vertical power-law shear",
            "Extended input: propagation direction",
        ):
            assert not any(channel in description for description in descriptions), f"Unexpected AD input channel: {channel}"
        for filename in sorted(expected - {f"{case_name}.1.AD.lin"}):
            validate_linearization(case_dir / filename)
        print("CompInflow=0 linearization passed: four expected files, no IfW file, no extended AD wind channels, finite consistent module and system outputs")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
