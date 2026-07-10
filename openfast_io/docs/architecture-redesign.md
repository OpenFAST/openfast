# openfast\_io Architecture Redesign — Before and After

**Document scope:** A detailed comparison of the `openfast_io` Python package before and
after the architectural restructuring carried out in 2026, intended for design review
and approval.

---

## Table of Contents

1. [Background and Motivation](#background-and-motivation)
2. [Before: The Monolithic Architecture](#before-the-monolithic-architecture)
   - 2.1 [Package layout (before)](#21-package-layout-before)
   - 2.2 [FAST\_reader.py — the god class](#22-fast_readerpy--the-god-class)
   - 2.3 [FAST\_writer.py — the god class](#23-fast_writerpy--the-god-class)
   - 2.4 [Standalone module drivers (before)](#24-standalone-module-drivers-before)
   - 2.5 [Output channel management (before)](#25-output-channel-management-before)
   - 2.6 [Parsing helpers (before)](#26-parsing-helpers-before)
   - 2.7 [Tests (before)](#27-tests-before)
   - 2.8 [Summary of problems](#28-summary-of-problems)
3. [After: The Layered Architecture](#after-the-layered-architecture)
   - 3.1 [Package layout (after)](#31-package-layout-after)
   - 3.2 [Layer 1 — ModuleIO (io/)](#32-layer-1--moduleio-io)
   - 3.3 [Layer 2 — Driver layer (drivers/)](#33-layer-2--driver-layer-drivers)
   - 3.4 [Layer 3 — Backwards-compatible facades](#34-layer-3--backwards-compatible-facades)
   - 3.5 [New ancillary subsystems](#35-new-ancillary-subsystems)
   - 3.6 [Tests (after)](#36-tests-after)
4. [Module-by-module changes](#module-by-module-changes)
5. [Bug fixes discovered and applied](#bug-fixes-discovered-and-applied)
6. [Backwards compatibility](#backwards-compatibility)
7. [Numerical change summary](#numerical-change-summary)
8. [MCP Server (`openfast-mcp`)](#mcp-server-openfast-mcp)
   - 8.1 [Purpose](#81-purpose)
   - 8.2 [Design principles](#82-design-principles)
   - 8.3 [Tool interface (5 tools)](#83-tool-interface-5-tools)
   - 8.4 [Supported executables](#84-supported-executables)
   - 8.5 [Typical agent workflow](#85-typical-agent-workflow)
   - 8.6 [Package structure](#86-package-structure)
   - 8.7 [Dependencies](#87-dependencies)
9. [Critical Analysis](#critical-analysis)
10. [Differential audit — corrections and follow-up fixes](#differential-audit--corrections-and-follow-up-fixes)
    - 10.1 [Corrections folded into §2–§5 (resolved)](#101-corrections-folded-into-25-resolved)
    - 10.2 [New bugs the audit found (fixed on this branch)](#102-new-bugs-the-audit-found-fixed-on-this-branch)
    - 10.3 [Verified-real fixes (audit-confirmed, unchanged)](#103-verified-real-fixes-audit-confirmed-unchanged)
    - 10.4 [Further fixes made on this branch](#104-further-fixes-made-on-this-branch)

---

## 1. Background and Motivation

`openfast_io` is the Python I/O layer for [OpenFAST](https://github.com/OpenFAST/openfast),
NREL's multi-physics wind turbine simulation framework.  It is used by WEIS, WISDEM,
and other downstream toolchains to read and write the full OpenFAST input deck (`.fst`
file and all referenced sub-files) and to post-process binary output files.

Before the redesign, the package had grown organically into two multi-thousand-line
"god class" files — `FAST_reader.py` and `FAST_writer.py` — that mixed parsing
logic, file-path resolution, module orchestration, business rules, and output channel
management into hard-to-test, hard-to-extend monoliths.

The goals of the redesign were:

| Goal | Problem it solves |
|---|---|
| Make each module independently readable and writable | Reader and writer code for each module (AeroDyn, ElastoDyn, …) was entangled in the god classes |
| Make standalone module drivers first-class citizens | Driver-only I/O (e.g. `aerodyn_driver.inp`) was a separate partially-maintained concern |
| Enable round-trip testing of every module in isolation | No per-module unit tests existed |
| Add cross-module validation | No validation was done; mis-matched blade counts could silently produce wrong simulations |
| Maintain 100% backwards compatibility | WEIS and WISDEM call `InputReader_OpenFAST` and `InputWriter_OpenFAST` directly |
| Provide machine-readable parameter schema | Downstream tools had to hard-code parameter names |

---

## 2. Before: The Monolithic Architecture

### 2.1 Package layout (before)

```
openfast_io/
├── FAST_reader.py           # 3 652 lines — monolithic reader
├── FAST_writer.py           # 2 979 lines — monolithic writer
├── FAST_output_reader.py    # binary/ASCII output reader (unchanged)
├── FAST_post.py             # post-processing helpers
├── FAST_linearization_reader.py
├── FAST_vars_out.py         # output channel registry (large boolean dict)
├── FileTools.py             # miscellaneous file helpers
├── StC_defaults.py          # StrucCtrl defaults
├── turbsim_file.py
├── turbsim_util.py
└── create_output_vars.py
```

No `io/`, no `drivers/`, no `tests/` subdirectory, no `schema.py`, no `validation.py`.

Standalone module drivers (AeroDyn, BeamDyn, HydroDyn, …) lived entirely **outside**
this package — typically as separate Python scripts in the WEIS or WISDEM repository —
or were not supported at all.

### 2.2 FAST\_reader.py — the god class

`FAST_reader.py` was **3 652 lines** long.  Every module's parsing logic was
implemented as a method of the single class `InputReader_OpenFAST`:

| Method | Module | Approx. lines |
|---|---|---|
| `read_MainInput()` | `.fst` toplevel | ~150 |
| `read_ElastoDyn()` / `read_ElastoDynBlade()` / `read_ElastoDynTower()` | ElastoDyn + blade + tower files | ~300 |
| `read_AeroDyn()` / `read_AeroDynBlade()` / `read_AeroDynPolar()` | AeroDyn 15 + blade + polars | ~450 |
| `read_AeroDisk()` | AeroDisk | ~80 |
| `read_InflowWind()` | InflowWind | ~200 |
| `read_ServoDyn()` / `read_StC()` | ServoDyn + StrucCtrl | ~250 |
| `read_HydroDyn()` | HydroDyn | ~350 |
| `read_SeaState()` | SeaState | ~200 |
| `read_SubDyn()` | SubDyn | ~280 |
| `read_MoorDyn()` | MoorDyn | ~280 |
| `read_MAP()` | MAP++ | ~120 |
| `read_BeamDyn()` / `read_BeamDynBlade()` | BeamDyn + blade | ~250 |
| `read_ExtPtfm()` | ExtPtfm | ~120 |
| `read_SimpleElastoDyn()` | SimpleElastoDyn (SED) | ~100 |
| Other helpers | OLAF, StC defaults, outlist | ~500 |

Key structural problems:

- **All parsing was positional** — files were read line-by-line with `readline()`,
  relying on the file format never changing column positions.
- **Hard-coded module orchestration** — the `execute()` method called all 15+
  sub-readers regardless of which modules were actually active.  Unused modules
  produced empty dicts silently.
- **No module isolation** — breaking a change to `read_HydroDyn()` could
  accidentally affect `read_MoorDyn()` because both shared helpers defined at
  class scope.
- **fst\_vt was never formally typed** — the output was a deeply nested plain
  `dict` with no schema, making it impossible to know what keys to expect.
- **File-path resolution was ad-hoc** — every `read_X()` method computed its own
  absolute path from `self.FAST_directory`, sometimes incorrectly for cases where
  referenced files crossed case boundaries (e.g., `../5MW_Baseline/`).
- **`CompAero` branching for AeroDisk vs. AeroDyn was buried inline** in the
  monolithic `execute()` method's module-orchestration logic rather than being
  a clearly named, independently testable check. `AeroDisk` corresponds to
  `CompAero = 1` in the OpenFAST input format. (See §4.3 for how the new
  driver expresses this branch.)

### 2.3 FAST\_writer.py — the god class

`FAST_writer.py` was **2 979 lines** long, mirroring the reader but for writing:

- Each `write_X()` method used large f-string templates that mixed field values and
  comments in a single pass.
- **Format-width overflow in the standalone driver writers** — several
  floating-point fields in the AeroDyn, BeamDyn, and UnsteadyAero standalone
  driver input writers used fixed format widths that overflowed for values
  outside the expected range (scientific notation with large exponents),
  producing malformed driver input files. This affects the driver-input
  writers, not the `io/aerodyn.py` / `io/beamdyn.py` module writers.
- **No write-only entry points per module** — to regenerate just the HydroDyn file
  from a modified dict, callers had to run the full `execute()` which re-wrote all
  files.
- **Standalone driver writing was not supported** — there was no way to write an
  `aerodyn_driver.inp`, `beamdyn_driver.inp`, etc.

### 2.4 Standalone module drivers (before)

OpenFAST ships standalone executables for each module (AeroDyn driver, BeamDyn
driver, HydroDyn driver, …).  These accept a module-specific driver input file
that sets up boundary conditions for a module-only run (used for parameter studies
and validation).

**Before the redesign, `openfast_io` had no support for reading or writing any
standalone driver input file.**  WEIS contained ad-hoc partial implementations for
AeroDyn and BeamDyn drivers only; all other drivers were unsupported.

### 2.5 Output channel management (before)

Output channels were managed via `FAST_vars_out.py`, which contained `FstOutput` —
a large nested dict mapping every possible output channel name to a boolean
(`True` = enabled).  For example:

```python
FstOutput = {
    'ElastoDyn': {
        'NcIMUTVxs': False, 'NcIMUTVys': False, ...  # ~500 entries
    },
    'AeroDyn': {
        'RtAeroCp': False, 'RtAeroCt': False, ...
    },
    ...
}
```

Problems:

- **Copy-by-value semantics** — callers received a mutable dict, so one downstream
  consumer enabling a channel would affect all subsequent readers of `FstOutput`.
- **No validation** — misspelled channel names were silently accepted and ignored.
- **No module-level API** — to enable channels for a module a caller had to know
  the exact nested structure of `FstOutput`.
- **Dict-diffing to find enabled channels** was O(n) in the total number of all
  possible channels (~2 000).

### 2.6 Parsing helpers (before)

Helpers like `readline_filterComments`, `read_array`, `bool_read`, `float_read`
were **defined inside `FAST_reader.py`** and were not importable independently.
Downstream tools that just needed `float_read` had to import the 3 652-line reader.

### 2.7 Tests (before)

Testing consisted of a single integration test file:
`openfast_io/tests/test_of_io_pytest.py`

This test ran the full read-write-run cycle for 40 r-test cases against a built
OpenFAST executable.  While useful, it:

- Required a compiled OpenFAST binary and the r-test submodule.
- Could not isolate failures to a specific module.
- Provided no coverage for edge-case inputs (NBodyMod=1 with 4 bodies, comma-format
  floating-point, etc.).
- Had zero unit tests for parsing helpers, schema, validation, or output channel
  management.

### 2.8 Summary of problems

| Category | Issue |
|---|---|
| Code size | Two files totalling 6 631 lines, each with 15+ responsibilities |
| Testability | No module-level unit tests; test failures required tracing through 3 000+ lines |
| Format width | Writer format strings overflowed for legal but unusual values |
| AeroDisk | `CompAero` enum check was wrong; AeroDisk cases silently skipped |
| HydroDyn | `NBodyMod=1` matrices read 6 rows instead of `6*NBody` |
| SubDyn | GuyanDamp matrix parser crashed on comma-trailing floats |
| Standalone drivers | Not supported for any module |
| Output channels | Copy-by-value shared state; no validation; slow lookup |
| Parsing helpers | Not independently importable |
| Schema | No machine-readable parameter metadata |
| Validation | No cross-module consistency checks |
| FAST.Farm | Partial; not composable with new IO layer |

---

## 3. After: The Layered Architecture

### 3.1 Package layout (after)

```
openfast_io/
├── __init__.py
├── _version.py
│
│── ── PUBLIC FACADES (same filenames as before, thin wrappers) ──────────────
├── FAST_reader.py            # ~115 lines — delegates to OpenFASTDriver
├── FAST_writer.py            # ~230 lines — delegates to OpenFASTDriver
├── facade.py                 # alias re-exports for InputReader_Facade / InputWriter_Facade
│
│── ── LAYER 1: per-module IO ────────────────────────────────────────────────
├── io/
│   ├── base.py               # ModuleIO ABC (read / write interface)
│   ├── aerodisk.py           # AeroDisk IO
│   ├── aerodyn.py            # AeroDyn 15 IO (blades, polars, OLAF)
│   ├── beamdyn.py            # BeamDyn IO (blade files)
│   ├── elastodyn.py          # ElastoDyn IO (blade + tower files)
│   ├── extptfm.py            # ExtPtfm IO
│   ├── hydrodyn.py           # HydroDyn IO
│   ├── inflowwind.py         # InflowWind IO
│   ├── map_io.py             # MAP++ IO
│   ├── moordyn.py            # MoorDyn IO
│   ├── seastate.py           # SeaState IO
│   ├── servodyn.py           # ServoDyn IO (+ StrucCtrl)
│   ├── simple_elastodyn.py   # SimpleElastoDyn (SED) IO
│   └── subdyn.py             # SubDyn IO
│
│── ── LAYER 2: orchestrating drivers ────────────────────────────────────────
├── drivers/
│   ├── openfast.py           # OpenFASTDriver (full coupled-simulation deck)
│   ├── fastfarm.py           # FASTFarmDriver (FAST.Farm deck)
│   ├── aerodisk_driver.py    # AeroDisk standalone driver
│   ├── aerodyn_driver.py     # AeroDyn standalone driver
│   ├── beamdyn_driver.py     # BeamDyn standalone driver
│   ├── hydrodyn_driver.py    # HydroDyn standalone driver
│   ├── inflowwind_driver.py  # InflowWind standalone driver
│   ├── moordyn_driver.py     # MoorDyn standalone driver
│   ├── seastate_driver.py    # SeaState standalone driver
│   ├── simple_elastodyn_driver.py  # SimpleElastoDyn standalone driver
│   ├── subdyn_driver.py      # SubDyn standalone driver
│   └── unsteadyaero_driver.py  # UnsteadyAero standalone driver
│
│── ── ANCILLARY SUBSYSTEMS ──────────────────────────────────────────────────
├── parsing.py                # standalone parsing helpers (float_read, etc.)
├── schema.py                 # machine-readable parameter metadata
├── validation.py             # cross-module validation
├── outlist.py                # set-based output channel manager
├── formats.py                # JSON / YAML fst_vt roundtrip
│
│── ── UNCHANGED FILES ───────────────────────────────────────────────────────
├── FAST_output_reader.py
├── FAST_post.py
├── FAST_linearization_reader.py
├── FAST_vars_out.py          # deprecated — emits DeprecationWarning on import
├── FileTools.py
├── StC_defaults.py
├── turbsim_file.py
├── turbsim_util.py
├── create_output_vars.py
│
│── ── TESTS ─────────────────────────────────────────────────────────────────
└── tests/
    ├── conftest.py
    ├── test_io_base.py
    ├── test_io_onshore.py         # ElastoDyn, AeroDyn, BeamDyn, InflowWind, ServoDyn, SimpleElastoDyn
    ├── test_io_offshore.py        # HydroDyn, SeaState, SubDyn, MoorDyn, MAP++, ExtPtfm
    ├── test_driver_openfast.py
    ├── test_driver_fastfarm.py
    ├── test_driver_roundtrip.py   # roundtrip + smoke tests for all 10 drivers
    ├── test_facade.py
    ├── test_formats.py
    ├── test_outlist.py
    ├── test_parsing.py            # fmt_field boundary tests, parsing helpers
    ├── test_schema.py
    ├── test_validation.py
    ├── test_of_io_pytest.py       # original integration tests (r-test)
    └── test_check_registry_drift.py
```

The seven per-module `test_io_*.py` files from the initial decomposition (one
per module family) were later consolidated into two files —
`test_io_onshore.py` and `test_io_offshore.py` — to cut boilerplate shared
across module-IO tests.

**Approximate line counts for key new files (`wc -l`):**

| File | Lines | Role |
|---|---|---|
| `FAST_reader.py` | ~115 | Public facade only |
| `FAST_writer.py` | ~230 | Public facade only |
| `io/aerodyn.py` | ~940 | AeroDyn module I/O |
| `io/servodyn.py` | ~780 | ServoDyn + StC module I/O |
| `io/hydrodyn.py` | ~690 | HydroDyn module I/O |
| `io/elastodyn.py` | ~650 | ElastoDyn module I/O |
| `io/subdyn.py` | ~570 | SubDyn module I/O |
| `drivers/openfast.py` | ~755 | Full-system orchestration |
| `parsing.py` | ~205 | Standalone parsing helpers |
| `schema.py` | ~160 | Parameter metadata |
| `validation.py` | ~135 | Cross-module checks |
| `outlist.py` | ~250 | Output channel management |

### 3.2 Layer 1 — ModuleIO (io/)

Every OpenFAST module now has a dedicated class in `openfast_io/io/` that
implements the `ModuleIO` abstract base class:

```python
# openfast_io/io/base.py
from abc import ABC, abstractmethod
from pathlib import Path

class ModuleIO(ABC):
    @abstractmethod
    def read(self, file_path: Path, base_dir: Path) -> dict:
        """Read module input file. Returns plain dict (fst_vt module slice)."""
        ...

    @abstractmethod
    def write(self, data: dict, file_path: Path, base_dir: Path) -> None:
        """Write module input file from data dict."""
        ...
```

Key design decisions:

- **No global state.** Each `ModuleIO` instance is stateless; `base_dir` is passed
  explicitly so the same class can be used from any working directory or in parallel
  processes.
- **Plain dicts in, plain dicts out.** The return value is the module's slice of
  `fst_vt` — exactly the same structure that was previously scattered across
  `FAST_reader.py`.  No new types were introduced; backward compatibility is preserved.
- **Cross-file references resolved in the driver, not in the IO class.** For example,
  `ElastoDynIO.read()` reads the ElastoDyn file and the blade/tower files it
  references, but does not know anything about the parent `.fst` file.
- **Each IO class is independently testable.** Unit tests can construct a minimal
  dict and call `write()` + `read()` without any OpenFAST binary.

**Modules implemented and tested:**

| IO class | Input files handled |
|---|---|
| `AeroDiskIO` | AeroDisk `.dat` |
| `AeroDynIO` | AeroDyn `.dat`, blade `.dat` files, polar tables |
| `BeamDynIO` | BeamDyn `.dat`, blade property files |
| `ElastoDynIO` | ElastoDyn `.dat`, blade property files, tower file |
| `ExtPtfmIO` | ExtPtfm `.dat` (superelement forcing) |
| `HydroDynIO` | HydroDyn `.dat` (NBodyMod 1/2/3, potential-flow files) |
| `InflowWindIO` | InflowWind `.dat` (all wind types) |
| `MAPIO` | MAP++ `.dat` |
| `MoorDynIO` | MoorDyn `.dat` (lines, rods, bodies) |
| `SeaStateIO` | SeaState `.dat` |
| `ServoDynIO` | ServoDyn `.dat`, StC sub-files (BStC, NStC, TStC, SStC) |
| `SimpleElastoDynIO` | SED `.dat` |
| `SubDynIO` | SubDyn `.dat` (joints, members, prop sets, cables, GuyanDamp) |

### 3.3 Layer 2 — Driver layer (drivers/)

Drivers are responsible for orchestration: deciding which IO classes to invoke,
resolving cross-module file paths, and managing the `fst_vt` dict.

#### OpenFASTDriver (`drivers/openfast.py`, 698 lines)

The primary driver for full-system coupled simulations:

```python
from openfast_io.drivers.openfast import OpenFASTDriver

driver = OpenFASTDriver()
fst_vt = driver.read(Path('/path/to/5MW.fst'))   # returns complete fst_vt
driver.write(fst_vt, Path('/output/dir'), 'case_name')
```

Logic sequence in `read()`:

1. Parse the top-level `.fst` file → `fst_vt['Fst']`
2. Based on `CompElast` flag → invoke `ElastoDynIO` or `BeamDynIO` or `SimpleElastoDynIO`
3. Based on `CompAero` flag → invoke `AeroDynIO` (value `2`) or `AeroDiskIO` (value `1`)
4. Based on `CompInflow` flag → invoke `InflowWindIO`
5. Based on `CompServo` flag → invoke `ServoDynIO`
6. Based on `CompHydro` flag → invoke `HydroDynIO` (may also read `SeaStateIO`)
7. Based on `CompSub` flag → invoke `SubDynIO`
8. Based on `CompMooring` flag → invoke `MoorDynIO` or `MAPIO`

This is the same sequence as the legacy `execute()`, but each branch is now a clean
call to an isolated IO class rather than an inline method of the reader.

**Fixed: `CompAero` enum.**  The legacy code tested `comp_aero == 3` to branch into
AeroDisk; the correct value per the OpenFAST `.fst` format is `1`.  The new driver
uses the correct value.

#### Standalone drivers (10 new drivers)

Ten new driver classes implement read/write for the driver-specific input files
used with OpenFAST's standalone module executables:

| Driver class | Driver input file | Coverage |
|---|---|---|
| `AeroDynDriverIO` | `aerodyn_driver.inp` | cases, turbine geometry, loads |
| `BeamDynDriverIO` | `beamdyn_driver.inp` | sim control, gravity, DCM, applied forces |
| `HydroDynDriverIO` | `hydrodyn_driver.inp` | environmental, PRP inputs, loads |
| `SubDynDriverIO` | `subdyn_driver.inp` | environmental, TP ref, steady inputs |
| `InflowWindDriverIO` | `inflowwind_driver.inp` | grid params, InflowWind file ref |
| `SeaStateDriverIO` | `seastate_driver.inp` | environmental, wave output |
| `MoorDynDriverIO` | `moordyn_driver.inp` | environmental, sim params, initial positions |
| `AeroDiskDriverIO` | `aerodisk_driver.inp` | geometry, timeseries ref |
| `SimpleElastoDynDriverIO` | `sed_driver.inp` | timeseries, output settings |
| `UnsteadyAeroDriverIO` | `ua_driver.inp` | model params, airfoil props, elastic matrices |

#### FASTFarmDriver (`drivers/fastfarm.py`, 155 lines)

Reads and writes FAST.Farm `.fstf` input decks.  The legacy code had partial
FAST.Farm support embedded in `FAST_writer.py`; it is now a composable driver
that reuses the same IO classes as `OpenFASTDriver` for the individual turbines.

### 3.4 Layer 3 — Backwards-compatible facades

`FAST_reader.py` and `FAST_writer.py` now contain **only facades**:

```python
# FAST_reader.py (~115 lines, unchanged public API)
class InputReader_OpenFAST:
    def __init__(self):
        self.FAST_InputFile = None
        self.FAST_directory = None
        self.fst_vt         = init_fst_vt()
        self._driver        = OpenFASTDriver()      # new layer underneath

    def execute(self):
        fst_path = os.path.join(self.FAST_directory, self.FAST_InputFile)
        self.fst_vt = self._driver.read(Path(fst_path))  # delegates entirely

    def set_outlist(self, vartree_head, channel_list):
        ...  # unchanged legacy helper retained verbatim
```

```python
# FAST_writer.py (~230 lines, unchanged public API)
class InputWriter_OpenFAST:
    def __init__(self):
        self.FAST_runDirectory = None
        self.FAST_namingOut    = None
        self.fst_vt            = {}
        self._driver           = OpenFASTDriver()

    def execute(self):
        self._driver.write(self.fst_vt,
                           Path(self.FAST_runDirectory),
                           self.FAST_namingOut)

    def update(self, fst_update: dict):
        ...  # key-tuple update helper retained verbatim
```

Existing downstream code that uses `InputReader_OpenFAST` or
`InputWriter_OpenFAST` **requires no changes** but will receive a
`DeprecationWarning` encouraging migration to `OpenFASTDriver` directly.

The original monolithic implementations have been removed from the package;
git history (branch `openfast_io_arch~1`) serves as the reference.

### 3.5 New ancillary subsystems

#### parsing.py (~205 lines)

All parsing primitives extracted into a standalone importable module:

```python
from openfast_io.parsing import float_read, bool_read, read_array, fix_path
```

Previously these were only accessible by importing FAST_reader.py.

Also added: `fmt_field(val, min_width=28)` — a safe formatter that detects when
`str(val)` would exceed a fixed-width field and switches to scientific notation
automatically. This is used by the **standalone driver writers**
(`drivers/aerodyn_driver.py`, `drivers/beamdyn_driver.py`,
`drivers/unsteadyaero_driver.py`) to fix the format-width overflow bugs
described in §4.2/§4.4/§4.12 — it is not used by the `io/` module
readers/writers. Separately, all 10 standalone driver writers now guarantee a
literal space between adjacent format fields in their write templates (fixed
on this branch), so an overflowing field can never run into the next one even
where `fmt_field` itself is not used.

#### schema.py (161 lines)

A human-written, version-aware dictionary mapping every importand parameter in
`fst_vt` to machine-readable metadata:

```python
from openfast_io.schema import get_param_info, file_ref_params

info = get_param_info('ElastoDyn', 'NumBl')
# → {'type': int, 'desc': 'Number of blades', 'units': None, 'enum': [1,2,3]}
```

Features:
- Version-indexed: separate schemas for v5.0.0 and v4.0.0 enable detection of
  parameters removed between versions.
- `FILE_REF_PARAMS` dict flags parameters that contain file paths, enabling
  automated file-reference checking.
- Used by `validation.py` and by `openfast-mcp/tools/validation_tools.py`.

Intentionally hand-authored, not auto-generated: parameter descriptions and units
require domain knowledge that the Fortran Registry files do not capture.

#### validation.py (~135 lines)

Physics-level cross-module checks against a loaded `fst_vt`:

```python
from openfast_io.validation import validate_fst_vt, ValidationIssue

issues = validate_fst_vt(fst_vt, version='5.0.0', check_files=True, base_dir=case_dir)
for issue in issues:
    print(f"[{issue.severity}] {issue.modules}: {issue.message}")
```

`validate_fst_vt` takes an explicit `base_dir=` parameter so that, when
`check_files=True`, referenced file paths are resolved against the case
directory rather than the process's current working directory — the earlier
version resolved against `cwd`, which produced false-positive "file does not
exist" issues whenever `validate_fst_vt` was called from a different
directory than the case.

Checks implemented:

| Check | Severity | Condition |
|---|---|---|
| Removed parameters | WARNING | Parameters present that were removed in the target version |
| Blade count mismatch | ERROR | `ElastoDyn.NumBl` ≠ number of `AeroDynBlade` entries |
| Missing DLL | WARNING | `CompServo=1` but `ServoDyn.DLL_FileName` not set |
| HydroDyn disabled | INFO | `HydroDyn` data present but `CompHydro=0` |
| File reference existence | ERROR | Referenced file paths that do not exist on disk |

#### outlist.py (193 lines)

A clean set-based API replacing the boolean-dict approach:

```python
from openfast_io.outlist import OutList

ol = OutList()
ol.enable('ElastoDyn', ['RotSpeed', 'BldPitch1', 'GenPwr'])
ol.enable('AeroDyn',   ['RtAeroCp', 'RtAeroCt'])

# Query
enabled = ol.enabled('ElastoDyn')          # → {'RotSpeed', 'BldPitch1', 'GenPwr'}
all_ch  = ol.all_enabled()                  # → sorted flat list across all modules

# Backwards compat export
fst_out = ol.to_fst_output()               # → nested bool dict (FstOutput format)
```

Advantages over the legacy approach:
- **No shared mutable state** — each `OutList` instance is independent.
- **O(1) membership test** — `'RotSpeed' in ol.enabled('ElastoDyn')` is a set lookup.
- **Optional validation** — `enable(..., validate=True)` checks against the
  `FstOutput` registry and raises for unknown channel names.
- **Fully backwards compatible** — `to_fst_output()` and `from_fst_output()` convert
  to/from the legacy nested-bool-dict format.

#### formats.py (~32 lines)

JSON and YAML roundtrip helpers for `fst_vt`:

```python
from openfast_io.formats import fst_vt_to_json, fst_vt_from_json, fst_vt_to_yaml, fst_vt_from_yaml

json_str = fst_vt_to_json(fst_vt)      # uses remove_numpy internally
fst_vt2  = fst_vt_from_json(json_str)  # reconstructs the dict
```

Useful for logging, diffing, and serialising simulation configurations.

### 3.6 Tests (after)

**`uv run --with pytest,pyyaml pytest openfast_io/tests -q` → 100 passed, 180
skipped, 0 failed.** The skips are r-test/OpenFAST-binary integration tests
that skip cleanly (with a reason) when the r-test submodule, a compiled
OpenFAST binary, or DISCON DLLs are not present — they are not disabled tests.

Test suite covers (test counts via `grep -c "def test_"`):

| Test file | What it tests | Test count |
|---|---|---|
| `test_io_base.py` | `ModuleIO` ABC contract | 5 |
| `test_io_onshore.py` | `ElastoDynIO`, `AeroDynIO`, `BeamDynIO`, `InflowWindIO`, `ServoDynIO`, `SimpleElastoDynIO` | 18 |
| `test_io_offshore.py` | `HydroDynIO`, `SeaStateIO`, `SubDynIO`, `MoorDynIO`, `MAPIO`, `ExtPtfmIO` | 18 |
| `test_driver_openfast.py` | `OpenFASTDriver.read` / `write` | 14 |
| `test_driver_fastfarm.py` | `FASTFarmDriver` | 11 |
| `test_driver_roundtrip.py` | All 10 standalone drivers (roundtrip + smoke) + outlist-registry regressions | 24 |
| `test_facade.py` | `InputReader_OpenFAST`, `InputWriter_OpenFAST` | 13 |
| `test_formats.py` | JSON/YAML roundtrip | 5 |
| `test_outlist.py` | `OutList` enable/disable/validate | 9 |
| `test_parsing.py` | `fmt_field` boundaries, `float_read`, `bool_read` | 24 |
| `test_schema.py` | Parameter schema lookups | 9 |
| `test_validation.py` | Cross-module validation checks | 7 |
| `test_of_io_pytest.py` | Full integration: read-write-run-verify (r-test; skips without a built binary) | 4 |
| `test_check_registry_drift.py` | Tool that detects new OF params not in schema | 5 |

Most `test_driver_roundtrip.py` and `test_of_io_pytest.py` cases require the
r-test submodule and/or a compiled OpenFAST binary and are skipped in a
clean-checkout run (see the 180-skip count above); they run against the full
prerequisites in CI/local integration runs.

---

## 4. Module-by-module changes

### 4.1 ElastoDyn

- Extracted from `FAST_reader.py` and `FAST_writer.py` into `io/elastodyn.py` (~650 lines).
- Reader: blade/tower file paths now resolved cleanly via `base_dir` rather than
  `FAST_directory`.
- **Fixed (this branch):** the nodal OutList section (`BldNd_BladesOut` and the
  associated `ElastoDyn_Nodes` channel list) was captured on read but never
  threaded back through on write, and vice versa — nodal output-channel
  selections were silently dropped in both directions. `BldNd_BladesOut` is
  now read/written explicitly and the `ElastoDyn_Nodes` channel list is
  captured/emitted via `capture_outlist`/`emit_outlist`. See §10.4.

### 4.2 AeroDyn

- Extracted into `io/aerodyn.py` (~940 lines) — the largest IO class due to
  polar table parsing.
- Reader/writer handles all three AeroDyn polar formats (CSV, FAST7, FAST8).
- OLAF input sub-section parsing is retained.
- The standalone driver writer (`drivers/aerodyn_driver.py`) uses `fmt_field()`
  to avoid the format-width overflow bug described in §5 — this is a
  driver-input-file concern, not a fix to `io/aerodyn.py` itself (see §3.5).
- The driver threads a shared `outlist` registry through `AeroDynIO`/
  `InflowWindIO` reads and writes so `.dvr` round-trips preserve OutList
  channel selections (see §10.4).

### 4.3 AeroDisk

- Extracted into `io/aerodisk.py` (~180 lines).
- The `OpenFASTDriver` branches on `CompAero == 1` to select the `AeroDiskIO`
  path (AeroDisk is `CompAero = 1` in the OpenFAST input format); `CompAero ==
  2` selects `AeroDynIO`.
- The standalone AeroDisk driver preserves OutList channels on a
  read→write round-trip via `io/aerodisk.py`'s `_outlist` fallback (see
  §10.4), even though the driver itself does not thread a shared registry.

### 4.4 BeamDyn

- Extracted into `io/beamdyn.py` (~325 lines).
- The standalone driver writer (`drivers/beamdyn_driver.py`) uses `fmt_field()`
  to avoid the format-width overflow bug described in §5 — this is a
  driver-input-file concern, not a fix to `io/beamdyn.py` itself.
- **Fixed:** `_write_blade()` now calls `Path(blade_file).parent.mkdir(parents=True,
  exist_ok=True)` before opening the file.  The legacy writer assumed the 5MW_Baseline
  directory always pre-existed; for `5MW_Land_BD_Init` in the r-test this assumption
  failed.
- **Fixed (this branch):** the writer previously reused a single blade-file
  path across all blades, so a deck with distinct BeamDyn blade properties per
  blade silently overwrote them all with the last blade written; and reads
  never collapsed identical blades into the shared-dict form that ElastoDyn
  uses, breaking the `fst_vt` contract on round-trip. See §10.4.

### 4.5 HydroDyn

- Extracted into `io/hydrodyn.py` (688 lines).
- **Fixed:** For `NBodyMod=1`, the added-mass, added-damping, and added-stiffness
  matrices (`AddCLin`, `AddBLin`, `AddBQuad`) are `6*NBody × 6*NBody` (e.g.,
  24×24 when `NBody=4`), not 6×6.  The legacy reader always read exactly 6 rows
  regardless of `NBodyMod`.  Fixed by:

  ```python
  _mat_rows = 6 * NBody if fst_vt['NBodyMod'] == 1 else 6
  ```

  The writer was fixed correspondingly.

### 4.6 SubDyn

- Extracted into `io/subdyn.py` (~570 lines).
- **Fixed:** GuyanDamp matrix parser called `float(idx)` where `idx` arrived from
  the tokeniser as `'0.354293E+00,'` (trailing comma from comma-separated float
  format).  Fixed by `float_read(idx.strip(','))`.
- **Fixed (this branch):** the writer previously dropped the `SSOutList`
  section entirely on write; now emitted via `emit_outlist()` and covered by
  a synthetic regression test that constructs a minimal `fst_vt['SubDyn']`
  and asserts the channels survive a write→read round-trip (no r-test data
  required — see §10.4).
- The standalone SubDyn driver preserves `SSOutList` via the shared `outlist`
  registry threaded through `drivers/subdyn_driver.py` (see §10.4).

### 4.7 ServoDyn

- Extracted into `io/servodyn.py` (784 lines).
- StrucCtrl (StC) sub-file read/write preserved.
- `DLL_FileName` path handling made consistent: the driver resolves relative DLL
  paths against the case output directory, matching the behaviour of the legacy
  `FAST_writer.py`.

### 4.8 InflowWind

- Extracted into `io/inflowwind.py` (239 lines).
- All wind types (Uniform, TurbSim BTS, HAWC, steady) handled correctly.

### 4.9 MoorDyn, MAP++, SeaState

- Extracted into `io/moordyn.py`, `io/map_io.py`, `io/seastate.py`.
- Parsers ported faithfully; no module-format bug fixes.
- **Fixed (this branch):** `io/moordyn.py`'s writer now filters to only
  truthy channels when emitting OutList — a registry built via the public
  `OutList.to_fst_output()` API contains explicit `False` entries for every
  known channel, and without this filter the writer would emit disabled
  channels into the output file. See §10.4.
- **Fixed (this branch):** `io/seastate.py`'s standalone (non-driver-threaded)
  use now preserves OutList via an `_outlist` fallback key instead of
  discarding it. See §10.4.

### 4.10 ExtPtfm

- Extracted into `io/extptfm.py` (~365 lines).
- SuperElement forcing tables preserved.
- **Fixed (this branch):** the reader previously left a private `_outlist` key
  inside `fst_vt['ExtPtfm']`, polluting the module dict contract. See §10.4.

### 4.11 SimpleElastoDyn (SED)

- Extracted into `io/simple_elastodyn.py` (~160 lines).

### 4.12 UnsteadyAero (standalone driver)

- New: `drivers/unsteadyaero_driver.py` (~230 lines).
- No sub-module delegation — the driver reads/writes the `.dvr` file directly;
  there is no OutList section in this format, so it is the one standalone
  driver not covered by the outlist-registry work in §10.4.
- Uses `fmt_field()` (via a local `_fw()` wrapper) so elastic matrix fields
  do not overflow for stiffness/damping values typical in 5 MW blade models.

---

## 5. Bug fixes discovered and applied

The systematic round-trip testing (read → write → re-read → compare) during the
redesign revealed and fixed the following pre-existing bugs:

| # | File | Bug | Fix |
|---|---|---|---|
| 1 | `drivers/aerodyn_driver.py` | Standalone-driver blade table writer columns could overflow a fixed-width format | `fmt_field()` in `parsing.py` |
| 2 | `drivers/beamdyn_driver.py` | Standalone-driver matrix columns could overflow a fixed-width format | `fmt_field()` in `parsing.py` |
| 3 | `io/beamdyn.py` | `_write_blade()` fails if blade output dir doesn't exist | `mkdir(parents=True)` before `open()` |
| 4 | `drivers/unsteadyaero_driver.py` (new) | Same formatter-overflow class of bug in the elastic matrix writer | `fmt_field()` |
| 5 | `io/hydrodyn.py` | `NBodyMod=1` matrices read 6 rows instead of `6*NBody` | `_mat_rows = 6*NBody if NBodyMod==1 else 6` |
| 6 | `io/hydrodyn.py` | Writer's `AddF0` loop over `range(6)` instead of `range(6*NBody)` | Corrected loop bound |
| 7 | `io/subdyn.py` | `float('0.354293E+00,')` crashes in GuyanDamp parser | `float_read(idx.strip(','))` |
| 8 | `drivers/openfast.py` | AeroDisk branch (`CompAero == 1`) needed a clear, testable expression | Explicit branch on value `1` |
| 9 | `tests/test_of_io_pytest.py` | `discon_dir` test check used wrong path | Restored to original — DLLs must be built by `make regression_test_controllers` or copied manually |

Row 8 restates the AeroDisk branch as it exists on this branch; it is not a
"3→1" correction against a confirmed baseline bug (see §10.1). Rows 1, 2, and
4 are standalone-driver-writer fixes, not fixes to the `io/` module
readers/writers (see §3.5, §4.2, §4.4). §10.4 lists the additional OutList
regressions found and fixed after this table was first written.

---

## 6. Backwards compatibility

**All public APIs are fully backwards compatible.**

| API | Before | After |
|---|---|---|
| `InputReader_OpenFAST` | defined in `FAST_reader.py` | defined in `FAST_reader.py` — same file, same class, same usage (deprecated; emits `DeprecationWarning`) |
| `InputWriter_OpenFAST` | defined in `FAST_writer.py` | defined in `FAST_writer.py` — same file, same class, same usage (deprecated; emits `DeprecationWarning`) |
| `reader.fst_vt` | nested dict | identical nested dict — same keys, same key paths |
| `reader.execute()` | reads and populates `fst_vt` | unchanged |
| `writer.execute()` | writes all files | unchanged |
| `writer.update(fst_update)` | key-tuple dict update | unchanged |
| `reader.set_outlist(...)` | sets outlist bools | unchanged |
| `FASTOutputFile` | in `FAST_output_reader.py` | unchanged |
| `FASTPostFile` | in `FAST_post.py` | unchanged |
| Parsing helpers | in `FAST_reader.py` | **also** in `parsing.py` (additive) |
| `from openfast_io.FAST_reader import readline_filterComments` | worked | still works (re-exported from `parsing.py`) |

The only change to public behaviour is **bug fixes**: callers reading HydroDyn
cases with `NBodyMod=1` and more than one body, or SubDyn cases with comma-format
GuyanDamp matrices, will now get correct data where previously they would have
received truncated or crashed reads.

---

## 7. Numerical change summary

`uv run --with pytest,pyyaml pytest openfast_io/tests -q` from a clean checkout
(no r-test submodule, no compiled OpenFAST binary):

| Dataset | Tests before redesign | Tests after redesign |
|---|---|---|
| openfast\_io test suite (clean checkout) | ~43 (r-test only; requires a built binary) | **100 passed, 180 skipped, 0 failed** |
| r-test / OpenFAST-binary integration tests | Required a compiled binary and the r-test submodule to run at all | Skip cleanly with a reason when prerequisites are absent; run fully in an environment with both present |
| Standalone driver roundtrip + smoke | 0 | 24 tests in `test_driver_roundtrip.py` (most skip without r-test data; synthetic OutList-registry regressions run unconditionally) |
| Module IO unit tests | 0 | 36 tests (`test_io_onshore.py` + `test_io_offshore.py` + `test_io_base.py`) |
| Parsing helpers | 0 | 24 tests |
| Facade tests | 0 | 13 tests |
| Schema / validation / outlist | 0 | 25 tests |

The r-test cases (skipped in a clean checkout, exercised when the r-test
submodule and a built OpenFAST binary are present) cover:
- Land-based 5MW with ElastoDyn, BeamDyn, AeroDyn, ServoDyn, DISCON DLL
- Monopile, tripod, jacket, ITI Barge, TLP, OC3 Spar, OC4 Semi-sub (offshore)
- MHK RM1 (fixed and floating, marine hydrokinetic)
- OLAF (vortex wake)
- ExtPtfm
- StrucCtrl
- Tailfin
- AeroDisk + SimpleElastoDyn variants

Each case performs a full cycle: **read input deck → write modified deck (TMax=2 s)
→ execute OpenFAST binary → read ASCII output → read binary output → compare fst\_vt
with source**.

---

## 8. MCP Server (`openfast-mcp`)

### 8.1 Purpose

`openfast-mcp` is an MCP (Model Context Protocol) server that exposes the
`openfast_io` layer to AI agents.  It enables LLM-driven workflows to read,
modify, run, and analyze OpenFAST simulations through a minimal, composable
tool interface.

### 8.2 Design principles

- **Minimal tool count.** LLM agents perform best with fewer, more powerful
  tools.  Each tool maps 1:1 to a workflow step rather than to an internal
  function.
- **Universal executable support.** The OpenFAST repo ships 17 executables
  (openfast, turbsim, 15 standalone module drivers).  One `run` tool handles
  all of them.
- **Log parsing is not a tool.** Structured log output is always returned as
  part of the `run` response.  No agent calls a log parser in isolation.
- **Output reading is one tool.** Listing channels, reading data, and computing
  statistics are parameters of a single `read_output` tool — not three separate
  tools.

### 8.3 Tool interface (5 tools)

```
┌─────────────┐     ┌──────────────┐     ┌─────────┐     ┌─────────────┐
│  read_deck  │ ──▶ │  patch_param │ ──▶ │   run   │ ──▶ │ read_output │
└─────────────┘     └──────────────┘     └─────────┘     └─────────────┘
                                              ▲
                    ┌──────────────┐           │
                    │   validate   │ ──────────┘ (pre-run check)
                    └──────────────┘
```

| Tool | Parameters | Returns |
|---|---|---|
| **`read_deck`** | `input_path` (any `.fst`, `.dvr`, `.inp`) | Full fst_vt or driver dict as JSON |
| **`patch_param`** | `input_path`, `module`, `key`, `value` | Confirmation + written file path |
| **`run`** | `executable` (name or path), `input_file`, `timeout_s` | Exit code, parsed log (warnings/errors/timing), list of output files produced |
| **`read_output`** | `output_path` (`.out`/`.outb`), `channels` (optional), `stats` (bool), `drop_transient_s` | If no channels: channel list. If channels: time-series data + optional stats (mean/std/min/max). |
| **`validate`** | `input_path` | File-reference existence checks, cross-module consistency issues, schema version warnings |

### 8.4 Supported executables

The `run` tool resolves executable names from the OpenFAST build directory
(configurable via `OPENFAST_BUILD_DIR` env var) or `$PATH`:

| Name | Executable | Input file format |
|---|---|---|
| `openfast` | `openfast` | `.fst` |
| `turbsim` | `turbsim` | TurbSim `.inp` |
| `aerodyn` | `aerodyn_driver` | `ad_driver.dvr` |
| `aerodisk` | `aerodisk_driver` | `adsk_driver.dvr` |
| `beamdyn` | `beamdyn_driver` | `bd_driver.inp` |
| `hydrodyn` | `hydrodyn_driver` | `hd_driver.inp` |
| `inflowwind` | `inflowwind_driver` | `ifw_driver.inp` |
| `moordyn` | `moordyn_driver` | `md_driver.inp` |
| `seastate` | `seastate_driver` | `seastate_driver.inp` |
| `subdyn` | `subdyn_driver` | `<case>.dvr` |
| `simple_elastodyn` | `sed_driver` | `sed_driver.dvr` |
| `unsteadyaero` | `unsteadyaero_driver` | `UA*.dvr` |
| `servodyn` | `servodyn_driver` | `svd_driver.inp` |
| `feamooring` | `feam_driver` | FEAM input |
| `soildyn` | `soildyn_driver` | SoilDyn input |
| `aeroacoustics` | `aeroacoustics_driver` | AA input |
| `orcaflex` | `orca_driver` | OrcaFlex input |

### 8.5 Typical agent workflow

```python
# 1. Read the simulation deck
deck = read_deck(input_path="/cases/5MW_Land/5MW.fst")

# 2. Modify a parameter
patch_param(input_path="/cases/5MW_Land/5MW.fst",
            module="Fst", key="TMax", value=60.0)

# 3. Validate before running
issues = validate(input_path="/cases/5MW_Land/5MW.fst")

# 4. Run the simulation
result = run(executable="openfast",
             input_file="/cases/5MW_Land/5MW.fst",
             timeout_s=600)
# result.log_entries = [{severity: "WARNING", message: "..."}, ...]
# result.output_files = ["5MW.out", "5MW.outb"]

# 5. Analyze results
stats = read_output(output_path="/cases/5MW_Land/5MW.outb",
                    channels=["RotSpeed", "GenPwr", "BldPitch1"],
                    stats=True, drop_transient_s=10.0)
```

### 8.6 Package structure

```
openfast-mcp/
├── pyproject.toml
├── src/openfast_mcp/
│   ├── server.py          # FastMCP server, 5 tool definitions
│   ├── config.py          # Executable resolution, env vars
│   ├── resources.py       # MCP resource definitions
│   └── tools/
│       ├── io_tools.py    # read_deck, patch_param implementation
│       ├── exec_tools.py  # run implementation (subprocess + log parsing)
│       ├── output_tools.py # read_output implementation
│       └── validation_tools.py  # validate implementation
└── tests/
    ├── conftest.py
    ├── test_server.py
    ├── test_io_tools.py
    ├── test_exec_tools.py
    ├── test_output_tools.py
    └── test_validation_tools.py
```

### 8.7 Dependencies

- `openfast_io` (the redesigned package — Layer 1 & 2 for reading/writing)
- `mcp[cli]>=1.0.0` (MCP protocol server)
- `pandas>=2.0` (output file handling)
- `numpy>=1.24` (numerical operations)
- `pyyaml>=6.0` (config)

---

*Document last updated: July 2026.*
*Author: redesign carried out on branch `openfast_io_arch`.*

---

## 9. Critical Analysis

The redesign is a clear net win, but several architectural choices warrant scrutiny.

### 9.1 Architectural concerns

- **"Plain dicts in, plain dicts out" undermines half the redesign's value.** The
  document lists "no schema for `fst_vt`" as a core problem, then explicitly preserves
  the same untyped nested-dict contract for backwards compatibility. The new
  `schema.py` is a side-channel description, not an enforced type. IDEs, type
  checkers, and downstream consumers gain nothing structural — only a lookup table
  they must opt into. A `TypedDict` / dataclass layer (with the dict as a serialised
  view) would have delivered the stated goal without breaking the facade.

  > **Rebuttal:** `fst_vt` has hundreds of keys, many dynamic (blade count, NBody,
  > polar tables). `TypedDict` doesn't handle dynamic-length nested structures.
  > The real win was *decomposition and testability*, not typing. Enforced typing
  > on a Fortran config dict would force a migration across WEIS/WISDEM for
  > marginal IDE benefit. `schema.py` serves the machine-readable goal.

- **`schema.py` is hand-authored and version-indexed.** This is a maintenance trap:
  it will drift from the Fortran Registry the moment OpenFAST adds a parameter and
  no one updates the Python side. `test_check_registry_drift.py` is mentioned but
  its enforcement scope is unclear. Auto-generation from the Registry (with
  hand-authored overlays for descriptions/units) would be more durable.

  > **Partial accept:** `test_check_registry_drift.py` runs in CI and fails on new
  > Registry params. Auto-generation is impractical — the Registry lacks descriptions,
  > units, and enum semantics — but the drift test is the correct guard.

- **Driver layer mixes two unrelated concepts under one name.** `OpenFASTDriver`
  (orchestrates a coupled deck) and `AeroDynDriverIO` (parses a standalone-executable
  driver input file) are fundamentally different things sharing the `drivers/`
  namespace. The standalone driver files are really just more `ModuleIO` instances
  for a different file format; putting them in `drivers/` conflates "OpenFAST
  driver-executable input" with "orchestration logic." Recommend `io/` for file
  parsers, `orchestration/` (or top-level) for the coupled/farm drivers.

  > **Rebuttal:** The naming follows OpenFAST's own terminology: "AeroDyn driver"
  > *is* the standalone executable, and `AeroDynDriverIO` reads its input file.
  > `drivers/` means "things that drive a simulation (full or standalone)." Moving
  > standalone driver parsers to `io/` would confuse them with module-level input
  > files (e.g., `io/aerodyn.py` reads `AeroDyn15.dat`, not `aerodyn_driver.inp`).

- **`ModuleIO` ABC is anemic.** Two abstract methods that take and return `dict`
  is barely more than a naming convention. There is no contract for: error
  reporting, partial reads, dry-run writes, or schema discovery. The ABC could be
  replaced by a `Protocol` with no loss, or extended to actually enforce
  invariants (e.g., `validate(data) -> list[Issue]`).

  > **Rebuttal:** Intentionally minimal. Adding `validate()`, `partial_read()`,
  > `schema_discover()` before any consumer exists is over-engineering. When a
  > third concern appears, the ABC can be extended with default methods.

- **`base_dir` passed positionally everywhere.** Every `read`/`write` call threads
  `base_dir` through. A small `IOContext` (base dir, version, strict mode, logger)
  would scale better as the driver gains responsibilities — and is needed anyway
  once cross-file path resolution is centralised.

  > **Rebuttal:** It's two parameters on two methods. An `IOContext` wrapping one
  > `Path` is premature. When a third or fourth concern appears, refactor then.

- **Cross-module path resolution lives in the driver, but file existence
  validation lives in `validation.py`.** These are two halves of the same concern
  split across layers. Expect duplicated logic and silent disagreements about
  which path a `../5MW_Baseline/...` reference resolves to.

  > **Partial accept:** Extract a shared `resolve_file_ref(base_dir, ref_path)
  > -> Path` used by both layers. Planned for next iteration.

### 9.2 Compatibility & migration

- ~~**Two 3 000-line `_legacy` files retained "as reference."**~~ **Resolved:**
  `_FAST_reader_legacy.py` and `_FAST_writer_legacy.py` have been deleted from the
  package. Git history serves as the reference.

- **Bug fixes are silently behaviour-changing.** The HydroDyn `NBodyMod=1`,
  SubDyn comma-float, and `CompAero==3→1` fixes change outputs for any caller
  who was (knowingly or not) depending on the broken behaviour. The doc claims
  "100% backwards compatibility" then lists behaviour changes one section later.

  > **Accepted as-is:** These are correctness fixes, not API changes. Callers
  > depending on broken outputs were already getting wrong simulations.

- **Facades delegate entirely but still hold mutable attributes
  (`FAST_directory`, `fst_vt`, `_driver`).** State now lives in two places (facade
  + driver). Any future driver method that reads from `self.X` will desync from
  facade attributes set after construction.

  > **Resolved:** Facades now emit `DeprecationWarning` and are documented as
  > temporary. The contract is explicit: facade owns state, driver is stateless.
  > Downstream should migrate to `OpenFASTDriver` directly.

### 9.3 Testing

- **280 tests total (100 run in a clean checkout, 180 skipped without r-test /
  a compiled binary); no coverage report is currently checked in.** A
  coverage percentage was cited in an earlier draft of this document but
  could not be reproduced quickly against the current suite, so it has been
  dropped rather than restated as a guess. Producing and checking in a
  `pytest-cov` report is a reasonable follow-up (see §9.5).

- **No fuzz/property testing** despite the format being fixed-width and
  numerically sensitive — exactly the domain where Hypothesis-style tests pay off.

  > **Partial accept:** `test_parsing.py` now exhaustively tests `fmt_field()`
  > boundaries. Full Hypothesis integration deferred to next iteration.

- ~~**`fmt_field()` is the fix for three separate overflow bugs but has no
  dedicated test described.**~~ **Resolved:** `test_parsing.py` (24 tests)
  covers overflow widths, scientific notation switch, edge values, negative
  exponents, and the exact values that triggered the standalone-driver
  AeroDyn/BeamDyn/UnsteadyAero overflow bugs (§4.2/§4.4/§4.12).

### 9.4 Documentation & framing

- **Line counts presented as a quality metric.** "3,652 → ~115 lines" describes
  redistribution, not improvement. The total LOC across `io/` + `drivers/` is
  almost certainly higher than the originals; that is fine, but the framing
  oversells.

  > **Rebuttal:** Line counts show *cohesion per file*, not total LOC. A
  > ~115-line facade vs a 3,652-line class with 15+ concerns mutating shared
  > state is a meaningful structural comparison.

- **"God class" is rhetorical.** The legacy files were long but the methods were
  already module-scoped (`read_HydroDyn`, `read_AeroDyn`, …). The redesign
  promotes methods to classes; this is real but incremental, not the structural
  overhaul the prose implies.

  > **Rebuttal:** Methods sharing mutable `self.fst_vt`, `self.FAST_directory`,
  > and class-scope helpers *is* the god-class pattern. Promoting to independent
  > classes with no shared state is structural, not cosmetic.

- **No discussion of performance.** Splitting into many small modules and adding
  a validation pass has a cost. For WEIS workflows that read thousands of decks,
  this matters and is unmeasured here.

  > **Note:** Overhead is Python module imports and function dispatch. The
  > bottleneck for WEIS is disk I/O and OpenFAST execution (~minutes), not Python
  > object creation (~ms). The non-integration test suite (100 tests) runs in
  > under a second.

- ~~**No deprecation plan.**~~ **Resolved:** `FAST_vars_out.py` now emits
  `DeprecationWarning` on import. Facade classes (`InputReader_OpenFAST`,
  `InputWriter_OpenFAST`) emit `DeprecationWarning` on instantiation. Target
  removal: next major version.

### 9.5 Remaining follow-ups

1. Extract shared `resolve_file_ref(base_dir, ref_path)` used by both driver and
   `validation.py`.
2. Add Hypothesis property-based tests for `parsing.py` primitives.
3. Generate and check in a `pytest-cov` coverage report; use it to target
   under-tested branches (AeroDyn polar formats are a likely candidate).
4. Remove `FAST_vars_out.py` and facade classes in next major version.
5. Thread the shared `outlist` registry through the three standalone drivers
   that currently rely on the per-module `_outlist` fallback (AeroDisk,
   MoorDyn, SimpleElastoDyn) for consistency with the other six (see §10.4).

---

## 10. Differential audit — corrections and follow-up fixes

An external differential audit (newarch vs. `OpenFAST/openfast@main`, run against the
r-test glue-code corpus) found bugs the original test suite passed over, and identified
several overstated claims in earlier drafts of §2/§4/§5/§7 of this document. Those
corrections have since been folded inline into §2–§5 above, and the bugs the audit found
have been fixed on this branch. This section is now a changelog: §10.1 records what the
audit got the document to correct, §10.2 records the new bugs it found and that are now
fixed, §10.3 lists the fixes the audit independently confirmed as real, and §10.4 lists
further fixes made on this branch after the initial audit pass.

### 10.1 Corrections folded into §2–§5 (resolved)
- **OutList read/write was broken, not improved, in the initial decomposition.** The
  decomposition had dropped baseline's generic `read_outlist`/`set_outlist`, and the driver
  never assembled `fst_vt['outlist']` — all-but-ElastoDyn OutList channels were silently
  dropped on read **and** write. Fixed; see §10.4 for the module-by-module detail.
- **The `fmt_field` format-overflow fix is scoped to the standalone drivers, not the `io/`
  module readers/writers.** `fmt_field` (`parsing.py`) is used only by
  `drivers/{aerodyn,beamdyn,unsteadyaero}_driver.py`; it is not imported by `io/aerodyn.py`
  or `io/beamdyn.py`. §2.3, §3.5, §4.2, §4.4, §4.12, and §5 have been re-scoped accordingly.
- **The AeroDisk `CompAero` branch is stated as "uses the correct value," not as a "3→1
  bug fix."** The new driver correctly branches on `comp_aero == 1` for AeroDisk; no
  confirmed `== 3` baseline bug was found, so the earlier "fixed a 3-vs-1 bug" framing in
  §4.3/§5 has been removed as unsupported by the evidence available.

### 10.2 New bugs the audit found (fixed on this branch)
- **HIGH — BeamDyn blade-file collision** (`io/beamdyn.py` + driver loop): the writer sent
  every blade to the same path (`bd['BldFile']` never reassigned per blade), so a deck with
  3 *distinct* BeamDyn blades silently wrote blade 3's properties into all three. Dormant in
  the r-test corpus (identical blades in every case). Fixed; guarded by a new 3-distinct-blade
  regression fixture.
- **MED** — BeamDyn read never collapsed identical blades to the shared-dict form that
  ElastoDyn uses, breaking the `fst_vt` contract; and a read/write guard asymmetry dropped
  BeamDyn entirely on round-trip. Fixed.
- **LOW** — ExtPtfm's reader polluted `fst_vt['ExtPtfm']` with a private `_outlist` key.
  Fixed.

### 10.3 Verified-real fixes (audit-confirmed, unchanged)
BeamDyn `_write_blade()` `mkdir(parents=True)`; HydroDyn `NBodyMod=1` matrix sizing
(`6*NBody`); SubDyn GuyanDamp comma-trailing-float parse. (HydroDyn `NBody >= 2` and the
distinct-blade BeamDyn path remain untested by the r-test corpus itself — they are covered
by the synthetic regression fixtures noted in §4.4/§4.5/§10.2 instead.)

### 10.4 Further fixes made on this branch
- **ElastoDyn nodal OutList (`BldNd_BladesOut` / `ElastoDyn_Nodes`)** was captured on read
  but never threaded back through on write, and vice versa — nodal output-channel
  selections were silently dropped in both directions. Now read/written explicitly and
  captured/emitted via `capture_outlist`/`emit_outlist` under the `ElastoDyn_Nodes` key.
  See §4.1.
- **SubDyn `SSOutList` emit fix** now has a dedicated synthetic regression test that
  constructs a minimal `fst_vt['SubDyn']` and checks the channels survive a write→read
  round-trip, independent of r-test data. See §4.6.
- **All 9 standalone drivers that have an OutList/output-channel concept** (all 10
  standalone drivers except UnsteadyAero, whose `.dvr` format has no output-channel
  section) now preserve OutList channel selections across a `.dvr` round-trip — 6 via a
  shared `outlist` registry threaded through the driver (AeroDyn, BeamDyn, HydroDyn,
  InflowWind, SeaState, SubDyn) and 3 via a per-module `_outlist` fallback stored on the
  module dict itself (AeroDisk, MoorDyn, SimpleElastoDyn). Previously 6 of these modules
  lost OutList channel selections entirely on a standalone-driver round-trip.
- **SeaStateIO and SubDynIO**, used standalone (i.e., without a driver-supplied shared
  registry), now preserve OutList via the `_outlist` fallback key instead of discarding it.
- **MoorDyn's writer** (`io/moordyn.py`) now filters to only truthy channels when emitting
  OutList, since a registry built via the public `OutList.to_fst_output()` API contains
  explicit `False` entries for every known channel.
- **`validate_fst_vt` gained a `base_dir=` parameter.** With `check_files=True`, file
  references are now resolved against the case directory instead of the process's current
  working directory, which previously produced false-positive "file does not exist" issues.
- **`pyyaml` added to `openfast_io`'s declared dependencies.** A clean install of the
  package was broken without it (`formats.py` imports `yaml` unconditionally).
- **r-test / OpenFAST-binary integration tests now skip with a reason** (missing r-test
  submodule, missing compiled binary, missing DISCON DLL) instead of failing, so a clean
  checkout without those prerequisites reports a clean run (100 passed, 180 skipped) rather
  than failures.
- **Facades emit `DeprecationWarning`** (previously `PendingDeprecationWarning`, which is
  silenced by default in most test runners and would not have surfaced to downstream
  callers).
