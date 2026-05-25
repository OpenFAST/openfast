"""OutList: clean interface for OpenFAST output channel management.

Replaces the per-channel bool-dict pattern in FAST_vars_out.FstOutput with
a set-based approach. The legacy bool-dict is still available via to_fst_output()
for backwards compatibility.
"""
from __future__ import annotations
import copy
from typing import Iterable

try:
    from .FAST_vars_out import FstOutput
except ImportError:
    FstOutput = {}


class OutList:
    """Manages output channel selections for an OpenFAST simulation.

    Internally stores enabled channels per module as sets of channel name strings.
    The FstOutput registry is used for validation only, not as the primary data store.

    Usage::

        ol = OutList()
        ol.enable('ElastoDyn', ['RotSpeed', 'BldPitch1', 'GenPwr'])
        ol.enable('AeroDyn', ['RtAeroCp', 'RtAeroCt'])
        print(ol.enabled('ElastoDyn'))  # {'RotSpeed', 'BldPitch1', 'GenPwr'}
        print(ol.all_enabled())          # flat list of all enabled channel names
    """

    def __init__(self):
        # module_name → set of enabled channel names
        self._channels: dict[str, set[str]] = {}

    def enable(self, module: str, channels: Iterable[str]) -> None:
        """Enable output channels for a module."""
        if module not in self._channels:
            self._channels[module] = set()
        self._channels[module].update(channels)

    def disable(self, module: str, channels: Iterable[str]) -> None:
        """Disable specific channels for a module."""
        if module in self._channels:
            self._channels[module] -= set(channels)

    def enabled(self, module: str) -> set[str]:
        """Get set of enabled channel names for a module."""
        return set(self._channels.get(module, set()))

    def all_enabled(self) -> list[str]:
        """Get flat sorted list of ALL enabled channel names across all modules."""
        result = []
        for module in sorted(self._channels):
            result.extend(sorted(self._channels[module]))
        return result

    def enabled_by_module(self) -> dict[str, list[str]]:
        """Get enabled channels grouped by module, sorted."""
        return {m: sorted(chs) for m, chs in sorted(self._channels.items()) if chs}

    def is_enabled(self, channel: str) -> bool:
        """Check if a channel is enabled in any module."""
        return any(channel in chs for chs in self._channels.values())

    def clear(self, module: str | None = None) -> None:
        """Clear all channels for a module, or all modules if None."""
        if module is None:
            self._channels.clear()
        elif module in self._channels:
            self._channels[module].clear()

    def validate(self) -> list[str]:
        """Return list of enabled channels not found in the FstOutput registry.

        These channels may cause OpenFAST to error at runtime.
        The FstOutput registry is auto-generated from
        openfast/docs/OtherSupporting/OutListParameters.xlsx.
        """
        registry = FstOutput  # read-only check — no deepcopy
        unknown = []
        for module, channels in self._channels.items():
            reg_module = registry.get(module, {})
            if not reg_module:
                unknown.extend(f"{module}.{ch}" for ch in sorted(channels))
            else:
                known = _flatten_keys(reg_module)
                for ch in sorted(channels):
                    if ch not in known:
                        unknown.append(f"{module}.{ch}")
        return unknown

    # ── Conversion to/from legacy bool-dict format ──

    def to_fst_output(self) -> dict:
        """Convert to the legacy FstOutput bool-dict format for backwards compatibility.

        Returns a fresh deepcopy of FstOutput with all channels set to False,
        then only the enabled channels set to True.
        """
        out = copy.deepcopy(FstOutput) if FstOutput else {}
        # Reset everything to False first
        for module in out:
            if isinstance(out[module], dict):
                _set_all_false(out[module])
        # Now enable only the channels in our sets
        for module, channels in self._channels.items():
            if module in out:
                _set_channels_in_dict(out[module], channels)
        return out

    @classmethod
    def from_fst_output(cls, fst_output: dict) -> 'OutList':
        """Create OutList from a legacy FstOutput bool-dict (e.g., fst_vt['outlist'])."""
        ol = cls()
        for module, vartree in fst_output.items():
            if isinstance(vartree, dict):
                enabled = _extract_true_channels(vartree)
                if enabled:
                    ol._channels[module] = enabled
        return ol

    # ── File I/O (used by drivers) ──

    @staticmethod
    def read_from_file(f, module: str) -> set[str]:
        """Read an OutList section from an open file handle.

        Works for both structured and free-form OutList formats.
        Returns the set of channel names found.
        """
        channels = set()
        data = f.readline()
        while data.strip() == '':
            data = f.readline()
        while data and not data.strip().startswith('END'):
            line = data.split('!')[0]  # strip comment
            line = line.split('-')[0] if '-' in line and '"' not in line.split('-')[0] else line.split('!')[0]
            for delim in ['"', "'", ',', ';', '\t']:
                line = line.replace(delim, ' ')
            tokens = [w.strip() for w in line.split() if w.strip()]
            channels.update(tokens)
            data = f.readline()
            while data.strip() == '':
                data = f.readline()
        return channels

    def write_to_file(self, f, module: str) -> None:
        """Write an OutList section to an open file handle."""
        channels = sorted(self._channels.get(module, set()))
        for ch in channels:
            f.write(f'"{ch}"\n')
        f.write('END of OutList section (the word "END" must appear in the first 3 columns of the last OutList line)\n')


def _flatten_keys(d: dict) -> set[str]:
    """Recursively collect all leaf keys from a nested dict."""
    keys = set()
    for k, v in d.items():
        if isinstance(v, dict):
            keys.update(_flatten_keys(v))
        else:
            keys.add(k)
    return keys


def _set_all_false(vartree: dict) -> None:
    """Recursively set all leaf bool values to False in a nested dict."""
    for k, v in vartree.items():
        if isinstance(v, dict):
            _set_all_false(v)
        elif isinstance(v, bool):
            vartree[k] = False


def _set_channels_in_dict(vartree: dict, channels: set[str]) -> None:
    """Recursively set matching channel keys to True in a nested dict."""
    for k, v in vartree.items():
        if isinstance(v, dict):
            _set_channels_in_dict(v, channels)
        elif k in channels:
            vartree[k] = True


def _extract_true_channels(vartree: dict) -> set[str]:
    """Recursively extract channel names set to True from a nested dict."""
    result = set()
    for k, v in vartree.items():
        if isinstance(v, dict):
            result.update(_extract_true_channels(v))
        elif v is True:
            result.add(k)
    return result
