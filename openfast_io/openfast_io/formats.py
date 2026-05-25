"""JSON and YAML serialization helpers for fst_vt dicts."""
import json
import yaml
from .FileTools import remove_numpy


def fst_vt_to_json(fst_vt: dict, indent: int = 2) -> str:
    """Serialize fst_vt to JSON. Numpy types are converted to Python natives."""
    return json.dumps(remove_numpy(fst_vt), indent=indent, default=str)


def fst_vt_from_json(json_str: str) -> dict:
    """Deserialize fst_vt from JSON.

    Note: numpy arrays become plain Python lists. InputWriter_OpenFAST
    handles both — it tests for list/array equivalence, not exact type.
    """
    return json.loads(json_str)


def fst_vt_to_yaml(fst_vt: dict) -> str:
    """Serialize fst_vt to YAML string."""
    return yaml.dump(remove_numpy(fst_vt), default_flow_style=False, sort_keys=False)


def fst_vt_from_yaml(yaml_str_or_path: str) -> dict:
    """Deserialize fst_vt from YAML string or file path."""
    try:
        with open(yaml_str_or_path) as f:
            return yaml.safe_load(f)
    except (OSError, TypeError):
        return yaml.safe_load(yaml_str_or_path)
