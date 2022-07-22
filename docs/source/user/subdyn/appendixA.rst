.. _sd_appendix_A:

Appendix A. OC4 Jacket Input File
=================================

SubDyn's primary input file 
:download:`(OC4 Jacket SubDyn's Input File) <./examples/OC4_Jacket_SD_Input.dat>`: 

This file includes information on the integration method (e.g., Adams-Bashforth 4 :sup:`th`} order), 
numerical-solution parameters (e.g., integration time interval, static solver flag, numer of modes to retain within the Craig-Bampton reduction), 
finite element analysis information (beam element model, number of elements per member),
and the geometric definition of the beam members via joints, member connectivity, and member cross-sectional properties.  
This file also specifies any SSI input files (soil/pile stiffness and mass matrices).
