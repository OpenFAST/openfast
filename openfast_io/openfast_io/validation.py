"""Validation layer for fst_vt dicts.

Performs schema-based checks, cross-module consistency checks, and
file-reference existence checks.
"""
import os
from dataclasses import dataclass, field
from pathlib import Path
from .schema import get_schema, FILE_REF_PARAMS


@dataclass
class ValidationIssue:
    severity: str        # 'ERROR' | 'WARNING' | 'INFO'
    modules: list[str] = field(default_factory=list)
    parameter: str | None = None
    message: str = ''
    source_file: str | None = None
    source_line: int | None = None


def validate_fst_vt(fst_vt: dict, version: str = '5.0.0',
                    check_files: bool = True,
                    base_dir: 'Path | str | None' = None) -> list[ValidationIssue]:
    """Validate a loaded fst_vt. Returns list of issues ordered by severity.

    Args:
        fst_vt: the loaded variable tree to validate.
        version: target OpenFAST version, used for removed-parameter checks.
        check_files: if True, verify that file-reference parameters point to
            existing files.
        base_dir: directory that relative file-reference paths in ``fst_vt``
            are resolved against (fst_vt stores paths relative to the case
            directory, not the current working directory). If ``None``,
            falls back to resolving against the current working directory.
    """
    issues = []
    issues.extend(_check_removed_params(fst_vt, version))
    issues.extend(_check_cross_module(fst_vt))
    if check_files:
        issues.extend(_check_file_refs(fst_vt, base_dir=base_dir))
    severity_order = {'ERROR': 0, 'WARNING': 1, 'INFO': 2}
    return sorted(issues, key=lambda i: severity_order.get(i.severity, 3))


def _check_removed_params(fst_vt: dict, version: str) -> list[ValidationIssue]:
    """Warn if fst_vt contains parameters removed in the target version."""
    issues = []
    from .schema import _SCHEMA
    if version == '5.0.0':
        removed_in_v5 = _SCHEMA.get('4.0.0', {})
        for module, params in removed_in_v5.items():
            present = fst_vt.get(module, {})
            if isinstance(present, dict):
                for param in params:
                    if param in present:
                        issues.append(ValidationIssue(
                            severity='WARNING', modules=[module], parameter=param,
                            message=f"{module}.{param} was removed in v5.0.0 and will be ignored. "
                                    f"See schema for migration guidance."
                        ))
    return issues


def _check_cross_module(fst_vt: dict) -> list[ValidationIssue]:
    """Physics-level cross-module consistency checks."""
    issues = []
    fst = fst_vt.get('Fst', {})
    ed = fst_vt.get('ElastoDyn', {})

    # Blade count consistency
    n_blades_ed = ed.get('NumBl', 3)
    ad_blades = fst_vt.get('AeroDynBlade', [])
    if isinstance(ad_blades, list) and ad_blades and len(ad_blades) != n_blades_ed:
        # Only flag if AeroDynBlade is non-empty and list (i.e., AeroDyn was read)
        has_data = any(bool(b) for b in ad_blades)
        if has_data:
            issues.append(ValidationIssue(
                severity='ERROR', modules=['ElastoDyn', 'AeroDyn'],
                parameter='NumBl',
                message=f"Blade count mismatch: ElastoDyn.NumBl={n_blades_ed} "
                        f"but AeroDynBlade has {len(ad_blades)} entries."
            ))

    # ServoDyn DLL not set when CompServo=1
    if fst.get('CompServo', 0) == 1:
        dll = fst_vt.get('ServoDyn', {}).get('DLL_FileName', '')
        if not dll:
            issues.append(ValidationIssue(
                severity='WARNING', modules=['Fst', 'ServoDyn'],
                parameter='DLL_FileName',
                message="CompServo=1 but ServoDyn.DLL_FileName is not set."
            ))

    # HydroDyn data present but disabled
    if fst.get('CompHydro', 0) == 0 and fst_vt.get('HydroDyn'):
        issues.append(ValidationIssue(
            severity='INFO', modules=['Fst', 'HydroDyn'], parameter=None,
            message="HydroDyn data loaded but CompHydro=0 — it will be ignored at runtime."
        ))

    return issues


def _check_file_refs(fst_vt: dict, base_dir: 'Path | str | None' = None) -> list[ValidationIssue]:
    """Check that file-reference parameters point to existing files.

    Relative paths are resolved against ``base_dir`` when given (fst_vt
    stores paths relative to the case directory). When ``base_dir`` is
    ``None``, relative paths are resolved against the current working
    directory instead, matching the historical behavior.
    """
    issues = []
    for module, param_names in FILE_REF_PARAMS.items():
        data = fst_vt.get(module, {})
        if not isinstance(data, dict):
            continue
        for param in param_names:
            val = data.get(param, '')
            if not val:
                continue
            paths = val if isinstance(val, list) else [val]
            for p in paths:
                p_str = str(p).strip('"').strip("'")
                if p_str and p_str.lower() not in ('unused', 'default', ''):
                    candidate = Path(p_str)
                    if base_dir is not None and not candidate.is_absolute():
                        candidate = Path(os.path.normpath(Path(base_dir) / candidate))
                    if not candidate.exists():
                        issues.append(ValidationIssue(
                            severity='ERROR', modules=[module], parameter=param,
                            message=f"{module}.{param} references a file that does not exist: {p_str}"
                        ))
    return issues
