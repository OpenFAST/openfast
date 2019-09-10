# OrcaFlex Interface Module
The legacy version of this module and additional documentation are available
the [NWTC Software Portal](https://nwtc.nrel.gov/OrcaFlexInterface/).

## Overview
OrcaFlex is a commercial software package developed by Orcina for the design
and analysis of marine systems. When the OrcaFlexInterface module is used in
OpenFAST, all hydrodynamic and mooring loads will be computed using OrcaFlex,
while the turbine, tower, and floating platform structural dynamics;
aerodynamics; and control and electrical-drive dynamics will be computed by
OpenFAST.

To use this module with OpenFAST, you will need the following:
- OpenFAST for Windows®
- A valid OrcaFlex license
- FASTlinkDLL.dll; This DLL is compiled by Orcina and is called by OpenFAST
  during the simulation to compute the loads on the platform by OrcaFlex.
  Both 32- and 64-bit versions of FASTlinkDLL.dll, which are compatible with
  the 32- and 64-bit Windows executable versions of OpenFAST, are available
  at https://orcina.com/Support/FASTlink.zip.

## Sample Models
Sample models for OpenFAST and OrcaFlexInterface can be downloaded
[here](https://nwtc.nrel.gov/enduser).

This self-extracting archive for Windows contains documentation on using the
OpenFAST-OrcaFlex interface as well as several sample models set to call
FASTlinkDLL.
