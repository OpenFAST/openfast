.. _AA-conclusions:

Conclusions
-----------

This document describes a set of frequency-based aeroacoustics models coupled to
the open-source aeroservoelastic solver OpenFAST. The goal of these models is to
predict the aeroacoustics emissions of wind turbine rotors. The document shows a
code-to-code comparison between the models coupled to OpenFAST and the models
implemented at the Technical University of Munich and coupled to the
aeroservoelastic solver Cp-Lambda. The comparison is performed simulating the
aeroacoustics emissions of the IEA Wind Task 37 land-based reference wind
turbine. The results show a good agreement between the two implementations. The
same turbine model is later used to exercise the aeroacoustics model showcasing
its capabilities. Finally, the appendices describe the entries of the input
files of OpenFAST to run the aeroacoustics analysis.

Future work will focus on the validation of the aeroacoustics models. 
In parallel, propagation models will be investigated and implemented. 
Finally, attention will be dedicated to infrasound noise and to the
time-domain models that can simulate it.

