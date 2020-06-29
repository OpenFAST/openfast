.. _bd_appendix:

Appendix
========

.. _bd_input_files:

BeamDyn Input Files
-------------------

In this appendix we describe the BeamDyn input-file structure and provide examples for the NREL 5MW Reference Wind Turbine.

OpenFAST+BeamDyn and stand-alone BeamDyn (static and dynamic) simulations all require two files:

1) BeamDyn primary input file 
:download:`(NREL 5MW static example) <examples/bd_primary_nrel_5mw.inp>`: This file includes information on the numerical-solution parameters (e.g., numerical damping, quadrature rules), and the geometric definition of the beam reference line via "members" and "key points".  This file also specifies the "blade input file."

2) BeamDyn blade input file :download:`(NREL 5MW example) <examples/nrel_5mw_blade.inp>`: 

Stand-alone BeamDyn simulation also require a driver input file; we list here examples for static and dynamic simulations:

3a) BeamDyn driver for dynamic simulations :download:`(NREL 5MW example) <examples/bd_driver_dynamic_nrel_5mw.inp>`: This file specifies the inputs for a single blade (e.g., forces, orientations, root velocity) and specifies the BeamDyn primary input file.

3b) BeamDyn driver for static simulations :download:`(NREL 5MW example) <examples/bd_driver_static_nrel_5mw.inp>`: Same as above but for static analysis.


.. _app-output-channel:

BeamDyn List of Output Channels
-------------------------------

This is a list of all possible output parameters for the BeamDyn module.
The names are grouped by meaning, but can be ordered in the OUTPUTS
section of the BeamDyn primary input file as the user sees fit.
N\ :math:`\beta`, refers to output node :math:`\beta`, where
:math:`\beta` is a number in the range [1,9], corresponding to entry
:math:`\beta` in the ``OutNd`` list. When coupled to FAST,
“:math:`B\alpha`” is prefixed to each output name, where :math:`\alpha`
is a number in the range [1,3], corresponding to the blade number. The
outputs are expressed in one of the following three coordinate systems:

-  **r**: a floating reference coordinate system fixed to the root of the
   moving beam; when coupled to FAST for blades, this is equivalent to
   the IEC blade (b) coordinate system.

-  **l**: a floating coordinate system local to the deflected beam.

-  **g**: the global inertial frame coordinate system; when coupled to FAST,
   this is equivalent to FAST’s global inertial frame (i) coordinate
   system.

.. _bd-output-channel:

.. figure:: figs/bd_output_channel.pdf
   :width: 500px
   :align: center

   BeamDyn Output Channel List
