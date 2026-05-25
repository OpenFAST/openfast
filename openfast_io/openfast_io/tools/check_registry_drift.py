"""
Developer tool: compare Fortran Registry files against openfast_io reader coverage.

Run after modifying any *_Registry.txt file to find parameters added, removed, or renamed.

Usage:
    python -m openfast_io.tools.check_registry_drift --openfast-root /path/to/openfast [--version 5.0.0]
"""
import re
import argparse
from dataclasses import dataclass
from pathlib import Path


REGISTRY_MAP = {
    'ElastoDyn':  ('modules/elastodyn/src/ElastoDyn_Registry.txt', 'ED_InputFile'),
    'AeroDyn':    ('modules/aerodyn/src/AeroDyn_Registry.txt', 'AD_InputFile'),
    'InflowWind': ('modules/inflowwind/src/InflowWind_Registry.txt', 'InflowWind_InputFile'),
    'ServoDyn':   ('modules/servodyn/src/ServoDyn_Registry.txt', 'SrvD_InputFile'),
    'HydroDyn':   ('modules/hydrodyn/src/HydroDyn_Registry.txt', 'HD_InputFile'),
    'SubDyn':     ('modules/subdyn/src/SubDyn_Registry.txt', 'SD_InputFile'),
    'BeamDyn':    ('modules/beamdyn/src/BeamDyn_Registry.txt', 'BD_InputFile'),
    'MoorDyn':    ('modules/moordyn/src/MoorDyn_Registry.txt', 'MD_InputFile'),
    'SeaState':   ('modules/seastate/src/SeaState_Registry.txt', 'SeaSt_InputFile'),
}


@dataclass
class RegistryParam:
    name: str
    fortran_type: str
    dims: str
    description: str
    units: str | None


def parse_registry(registry_path: Path, type_name: str) -> list[RegistryParam]:
    """Extract all fields from a specific typedef block in a registry file."""
    params = []
    current_module = None
    with open(registry_path, errors='replace') as f:
        for line in f:
            line = line.rstrip()
            if not line or line.lstrip().startswith(('#', '!')):
                continue
            parts = re.split(r'\t+', line)
            if len(parts) < 6 or parts[0].strip() != 'typedef':
                continue
            module = parts[1].strip()
            if module == '^':
                module = current_module
            else:
                current_module = module
            if parts[2].strip() != type_name:
                continue
            desc = parts[8].strip().strip('"') if len(parts) > 8 else ''
            units = parts[9].strip() if len(parts) > 9 else None
            params.append(RegistryParam(
                name=parts[4].strip(),
                fortran_type=parts[3].strip(),
                dims=parts[5].strip(),
                description=desc,
                units=units if units and units != '-' else None,
            ))
    return params


def scan_reader_params(reader_path: Path, module_key: str) -> set[str]:
    """Scan FAST_reader.py for fst_vt['ModuleKey']['ParamName'] assignments."""
    source = reader_path.read_text(errors='replace')
    pattern = r"""fst_vt\[['"]""" + re.escape(module_key) + r"""['"]\]\s*\[['"](\w+)['"]\]"""
    return set(re.findall(pattern, source))


def scan_schema_params(module_key: str) -> set[str]:
    """Get parameter names currently defined in schema.py for this module."""
    try:
        from openfast_io.schema import get_schema
        return set(get_schema(module_key).keys())
    except ImportError:
        return set()


def check_drift(openfast_root: Path, version: str = '5.0.0'):
    reader_path = openfast_root / 'openfast_io/openfast_io/FAST_reader.py'
    if not reader_path.exists():
        print(f"ERROR: FAST_reader.py not found at {reader_path}")
        return

    total_missing_reader = 0
    total_missing_schema = 0
    total_stale = 0

    for module, (reg_rel, type_name) in REGISTRY_MAP.items():
        reg_path = openfast_root / reg_rel
        if not reg_path.exists():
            print(f"\n{module}: ✗ registry not found ({reg_path})")
            continue

        reg_params = parse_registry(reg_path, type_name)
        reg_names = {p.name for p in reg_params}
        reader_names = scan_reader_params(reader_path, module)
        schema_names = scan_schema_params(module)

        missing_reader = reg_names - reader_names
        missing_schema = reg_names - schema_names
        stale = reader_names - reg_names

        total_missing_reader += len(missing_reader)
        total_missing_schema += len(missing_schema)
        total_stale += len(stale)

        status = '✓' if not missing_reader and not stale else '⚠'
        print(f"\n{module} ({type_name}: {len(reg_names)} params in registry) {status}")

        if not missing_reader and not missing_schema and not stale:
            print("  ✓ Fully covered")
            continue
        if not reader_names:
            print("  ✗ NO READER — entire module unimplemented in openfast_io")
        elif missing_reader:
            print(f"  ⚠ READER MISSING ({len(missing_reader)}):")
            for n in sorted(missing_reader)[:20]:
                p = next((x for x in reg_params if x.name == n), None)
                desc = f'  [{p.fortran_type}] "{p.description}"' if p else ''
                print(f"       - {n}{desc}")
            if len(missing_reader) > 20:
                print(f"       ... and {len(missing_reader) - 20} more")
        if missing_schema:
            print(f"  ⚠ SCHEMA MISSING ({len(missing_schema)})")
        if stale:
            print(f"  ✗ STALE IN READER ({len(stale)}) — may have been removed/renamed:")
            for n in sorted(stale)[:10]:
                print(f"       - {n}")
            if len(stale) > 10:
                print(f"       ... and {len(stale) - 10} more")

    print(f"\n{'─'*60}")
    print(f"Summary: {total_missing_reader} reader gaps, "
          f"{total_missing_schema} schema gaps, "
          f"{total_stale} stale entries")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Check openfast_io reader vs. registry drift')
    parser.add_argument('--openfast-root', required=True, type=Path)
    parser.add_argument('--version', default='5.0.0')
    args = parser.parse_args()
    check_drift(args.openfast_root, args.version)
