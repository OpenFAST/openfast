"""
FAST_reader.py -- backwards-compatible facade over the new composable IO layer.

This file preserves the public API that WEIS and other downstream code depends
on (``InputReader_OpenFAST``, parsing helpers, ``set_outlist``, etc.) while
delegating all heavy lifting to ``OpenFASTDriver`` and the per-module IO
classes under ``openfast_io.io``.

The original monolithic implementation is available in git history
(branch ``openfast_io_arch~1``) for reference.
"""
from __future__ import annotations

import copy
import os
from functools import reduce
from pathlib import Path

import operator

# -- Re-exports from parsing.py (backwards-compat for external importers) -----
from openfast_io.parsing import (        # noqa: F401
    readline_filterComments,
    readline_ignoreComments,
    read_array,
    fix_path,
    bool_read,
    float_read,
    int_read,
    quoted_read,
)

# -- New driver layer ----------------------------------------------------------
from openfast_io.drivers.openfast import OpenFASTDriver, init_fst_vt


class InputReader_OpenFAST:
    """OpenFAST input file reader -- backwards-compatible facade.

    .. deprecated::
        This facade class is deprecated and will be removed in a future version.
        Use ``openfast_io.drivers.openfast.OpenFASTDriver`` directly instead.

    Usage is identical to the legacy reader::

        reader = InputReader_OpenFAST()
        reader.FAST_InputFile = '5MW.fst'
        reader.FAST_directory = '/path/to/case'
        reader.execute()
        print(reader.fst_vt['Fst']['TMax'])
    """

    def __init__(self):
        import warnings
        warnings.warn(
            "InputReader_OpenFAST is deprecated. "
            "Use openfast_io.drivers.openfast.OpenFASTDriver directly.",
            PendingDeprecationWarning,
            stacklevel=2,
        )
        self.FAST_InputFile = None   # FAST input file (ext=.fst)
        self.FAST_directory = None   # Path to fst directory files
        self.path2dll       = None   # Path to controller DLL
        self.fst_vt         = init_fst_vt()
        self._driver        = OpenFASTDriver()

    # -- Core API --------------------------------------------------------------

    def execute(self):
        """Read the full simulation deck.  Populates ``self.fst_vt``."""
        fastdir = '' if self.FAST_directory is None else self.FAST_directory
        fst_path = os.path.join(fastdir, self.FAST_InputFile)
        self.fst_vt = self._driver.read(Path(fst_path))

    # -- Outlist helpers (unchanged from legacy -- pure dict manipulation) ------

    def set_outlist(self, vartree_head, channel_list):
        """Recursively set output channel names to True in the nested outlist dict."""

        def get_dict(vartree, branch):
            return reduce(operator.getitem, branch, vartree_head)

        def set_dict(vartree, branch, val):
            get_dict(vartree, branch[:-1])[branch[-1]] = val

        def loop_dict(vartree, search_var, branch):
            for var in vartree.keys():
                branch_i = copy.copy(branch)
                branch_i.append(var)
                if isinstance(vartree[var], dict):
                    loop_dict(vartree[var], search_var, branch_i)
                else:
                    if var == search_var:
                        set_dict(vartree_head, branch_i, True)

        for var in channel_list:
            var = var.replace(' ', '')
            loop_dict(vartree_head, var, [])


if __name__ == "__main__":
    from openfast_io.FileTools import check_rtest_cloned

    parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))) + os.sep

    fast = InputReader_OpenFAST()
    fast.FAST_InputFile = '5MW_Land_BD_DLL_WTurb.fst'
    fast.FAST_directory = os.path.join(
        parent_dir, 'reg_tests', 'r-test',
        'glue-codes', 'openfast',
        '5MW_Land_BD_DLL_WTurb',
    )
    check_rtest_cloned(fast.FAST_directory)
    fast.execute()
