"""
FAST_writer.py -- backwards-compatible facade over the new composable IO layer.

This file preserves the public API that WEIS and other downstream code depends
on (``InputWriter_OpenFAST``, ``update``, ``get_outlist``, ``update_outlist``,
helper functions, etc.) while delegating all heavy lifting to
``OpenFASTDriver.write()`` and the per-module IO classes under ``openfast_io.io``.

The original monolithic implementation is available in git history
(branch ``openfast_io_arch~1``) for reference.
"""
from __future__ import annotations

import copy
import os
import operator
from functools import reduce
from pathlib import Path

import numpy as np


# -- Re-exported helper functions (external code may import these) -------------

def auto_format(f, var):
    """Error handling for variables with 'Default' options."""
    if isinstance(var, str):
        f.write('{:}\n'.format(var))
    elif isinstance(var, int):
        f.write('{:3}\n'.format(var))
    elif isinstance(var, float):
        f.write('{: 2.15e}\n'.format(var))


def float_default_out(val, trim=False):
    """Formatted float output when 'default' is an option."""
    if type(val) is float:
        return '{:.4f}'.format(val) if trim else '{: 22f}'.format(val)
    else:
        return '{:}'.format(val) if trim else '{:<22}'.format(val)


def int_default_out(val, trim=False):
    """Formatted int output when 'default' is an option."""
    if type(val) is float:
        return '{:d}'.format(int(val)) if trim else '{:<22d}'.format(int(val))
    else:
        return '{:}'.format(val) if trim else '{:<22}'.format(val)


def get_dict(vartree, branch):
    """Given a list of nested dictionary keys, return the dict at that point."""
    return reduce(operator.getitem, branch, vartree)


# -- New driver layer ----------------------------------------------------------
from openfast_io.drivers.openfast import OpenFASTDriver


class InputWriter_OpenFAST:
    """OpenFAST input file writer -- backwards-compatible facade.

    .. deprecated::
        This facade class is deprecated and will be removed in a future version.
        Use ``openfast_io.drivers.openfast.OpenFASTDriver`` directly instead.

    Usage is identical to the legacy writer::

        writer = InputWriter_OpenFAST()
        writer.fst_vt = reader.fst_vt
        writer.FAST_runDirectory = '/output/path'
        writer.FAST_namingOut = 'my_case'
        writer.update(fst_update={'Fst': {'TMax': 20.0}})
        writer.execute()
    """

    def __init__(self):
        import warnings
        warnings.warn(
            "InputWriter_OpenFAST is deprecated. "
            "Use openfast_io.drivers.openfast.OpenFASTDriver directly.",
            PendingDeprecationWarning,
            stacklevel=2,
        )
        self.FAST_namingOut = None     # Base name for output files
        self.FAST_runDirectory = None  # Output directory
        self.fst_vt = {}
        self.fst_update = {}
        self._driver = OpenFASTDriver()

    # -- Core API --------------------------------------------------------------

    def execute(self):
        """Write all enabled module input files.  Delegates to ``OpenFASTDriver.write()``."""
        if self.FAST_runDirectory is None:
            raise ValueError('FAST_runDirectory must be set before calling execute()')
        if self.FAST_namingOut is None:
            raise ValueError('FAST_namingOut must be set before calling execute()')

        if not os.path.exists(self.FAST_runDirectory):
            os.makedirs(self.FAST_runDirectory)

        self._driver.write(
            self.fst_vt,
            Path(self.FAST_runDirectory),
            self.FAST_namingOut,
        )

    # -- Update helpers (unchanged from legacy -- pure dict manipulation) ------

    def update(self, fst_update={}):
        """Apply user-supplied overrides to ``self.fst_vt``."""
        if fst_update:
            self.fst_update = fst_update

        def loop_dict(vartree, branch):
            for var in vartree.keys():
                branch_i = copy.copy(branch)
                branch_i.append(var)
                if type(vartree[var]) is dict:
                    loop_dict(vartree[var], branch_i)
                else:
                    try:
                        get_dict(self.fst_vt, branch_i[:-1])[branch_i[-1]] = (
                            get_dict(self.fst_update, branch_i[:-1])[branch_i[-1]]
                        )
                    except (KeyError, TypeError):
                        pass

        if self.fst_update:
            # Check if keys are tuples (WEIS-style: {('Fst','TMax'): 20.0})
            first_key = next(iter(self.fst_update))
            if isinstance(first_key, tuple):
                fst_update_orig = copy.copy(self.fst_update)
                self.fst_update = {}
                for var_list in fst_update_orig.keys():
                    branch = []
                    for i, var in enumerate(var_list[0:-1]):
                        if var not in get_dict(self.fst_update, branch).keys():
                            get_dict(self.fst_update, branch)[var] = {}
                        branch.append(var)
                    get_dict(self.fst_update, branch)[var_list[-1]] = fst_update_orig[var_list]

            loop_dict(self.fst_update, [])

    # -- Outlist helpers (unchanged from legacy) -------------------------------

    def get_outlist(self, vartree_head, channel_list=[]):
        """Recursively find values set to True in the nested outlist dict."""

        def loop_dict(vartree, outlist_i):
            for var in vartree.keys():
                if type(vartree[var]) is dict:
                    loop_dict(vartree[var], outlist_i)
                else:
                    if vartree[var]:
                        outlist_i.append(var)
            return outlist_i

        if not channel_list:
            channel_list = vartree_head.keys()

        outlist = []
        for var in channel_list:
            var = var.replace(' ', '')
            outlist_i = loop_dict(vartree_head[var], [])
            if outlist_i:
                outlist.append(sorted(outlist_i))

        return outlist

    def update_outlist(self, channels):
        """Set output channels to specified boolean values."""

        def get_outlist_dict(vartree, branch):
            return reduce(operator.getitem, branch, self.fst_vt['outlist'])

        def set_outlist_dict(vartree, branch, val):
            get_outlist_dict(vartree, branch[:-1])[branch[-1]] = val

        def loop_dict(vartree, search_var, val, branch):
            for var in vartree.keys():
                branch_i = copy.copy(branch)
                branch_i.append(var)
                if type(vartree[var]) is dict:
                    loop_dict(vartree[var], search_var, val, branch_i)
                else:
                    if var == search_var:
                        set_outlist_dict(self.fst_vt['outlist'], branch_i, val)

        channel_list = channels.keys()
        for var in channel_list:
            val = channels[var]
            var = var.replace(' ', '')
            loop_dict(self.fst_vt['outlist'], var, val, [])


if __name__ == "__main__":
    from openfast_io.FAST_reader import InputReader_OpenFAST
    from openfast_io.FileTools import check_rtest_cloned

    fst_update = {}
    fst_update['Fst', 'TMax'] = 20.
    fst_update['AeroDyn', 'TwrAero'] = False

    parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))) + os.sep
    build_of_io_dir = os.path.join(parent_dir, 'build_ofio')
    Path(build_of_io_dir).mkdir(parents=True, exist_ok=True)

    # Read the model
    fast = InputReader_OpenFAST()
    fast.FAST_InputFile = '5MW_Land_BD_DLL_WTurb.fst'
    fast.FAST_directory = os.path.join(
        parent_dir, 'reg_tests', 'r-test',
        'glue-codes', 'openfast',
        '5MW_Land_BD_DLL_WTurb',
    )
    check_rtest_cloned(fast.FAST_directory)
    fast.execute()

    # Write out the model
    fastout = InputWriter_OpenFAST()
    fastout.fst_vt = fast.fst_vt
    fastout.FAST_runDirectory = os.path.join(build_of_io_dir, 'fast_write_main_test')
    fastout.FAST_namingOut = '5MW_Land_BD_DLL_WTurb_write'
    fastout.update(fst_update=fst_update)
    fastout.execute()
