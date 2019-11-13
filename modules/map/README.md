# Map++ Module
This is an externally developed module with further information
available on the developer's documentation site:
[map++ readthedocs](https://map-plus-plus.readthedocs.io/en/latest).

The legacy version of FAST's information regarding this module
are available at the [NWTC Software Portal](https://nwtc.nrel.gov/MAP/).

## Overview
The Mooring Analysis Program (MAP++) is a library designed to be used in
parallel with other CAE tools to model the steady-state forces on a
Multi-Segmented, Quasi-Static (MSQS) mooring line. The MSQS model is developed
based on an extension of conventional single line static solutions.
Conceptually, MAP++'s MSQS module solves the algebraic equations for all
elements simultaneously with the condition that the total force at connection
points sum to zero. Seabed contact, seabed friction, and externally applied
forces can be modeled with this tool. This allows multi-element mooring lines
with arbitrary connection configurations to be analyzed.

Because MAP++ is compiled as a library, it can be linked with other programs at
run-time. Alternatively, one may use MAP++'s native Python binding routines to
access the program through Python. These features give users the option to
execute MAP++ as a stand-alone design or simulation tool. The entry points into
MAP++ follow the function calling conventions outlined by the NWTC FAST
modularization framework.
