# IceFloe Module
The legacy version of this module and additional documentation are available
the [NWTC Software Portal](https://nwtc.nrel.gov/IceFloe/).

## Overview
The U.S. Department of Energy (DOE) awarded DNV GL a project to create a model
for interaction of bottom-fixed offshore wind turbines with surface ice for use
with common simulation tools. The IceFloe module, which conforms to the
standards of the FAST Modularization Framework, was created from this project.

The IceFloe module includes the option to apply loads from several models
(listed in table below) to a monopole structure or a multi-leg structure of
either three or four legs of the same diameter. For multi-leg support
structures, the ice loads are calculated independently for each leg; however
factors based on sheltering of one leg by another are also used, which can be
automatically applied or user specified. The project's final technical report
provides more details.

| Ice Load Model in IceFloe     | Source                 |
| ----------------------------- | ---------------------- |
| Continuous random crushing    | ISO 19906, Karna       |
| Intermittent crushing per ISO | ISO 19906              |
| Lock-in crushing per ISO      | ISO 19906              |
| Lock-in crushing per IEC      | IEC 61400-3, Korzhavin |
| Coupled crushing              | Määtänen               |
| Flexural failure per ISO      | ISO 19906, Croasdale   |
| Flexural failure per IEC      | IEC 61400-3, Ralston   |

## Manual
IceFloe documentation is available
[here](https://nwtc.nrel.gov/system/files/DDRP0133-IceLoadFinalReport2014_10_30.pdf).
