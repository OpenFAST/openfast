.. _FF:FutureWork:

Future Work
===========

This list contains features that could be implemented in future
releases:

-  Develop more efficient methods of generating/processing ambient wind
   from a high-fidelity precursor simulation, including:

   -  Propagate 2D planes of ambient wind data using Taylor’s frozen
      turbulence hypothesis as an alternative to 3D volumes

   -  Allow for nonuniform grids in Turbsim

   -  Use Dynamic Mode Decomposition to compress the file size of the
      low-resolution domains

   -  Implement Gabor mode enrichment to replace the high-resolution
      domains

   -  Develop a more efficient ABLSolver based on a simple rectangular
      (rather than a generally unstructured) grid.

-  Improve the eddy-viscosity formulation with additional physics.

-  Pursue additional wake-modeling approaches, including:

   -  Introduce simpler wake-deficit models, e.g., the Gaussian wake
      model by Bastankhah and Porté-Agel and the super-Gaussian model by
      Blondel and Cathelain

   -  Introduce simpler wake-deflection models, e.g., the model by
      Jiménez or the model by Qian and Ishihara

   -  Apply a free-vortex method for the near wake

   -  Incorporate a kidney-shaped wake under skewed-flow conditions,
      e.g., by incorporating opposing vortices from the skew-induced
      horseshoe vortex

   -  Deform the base-wake deficit (introduce asymmetry) as a result of
      background turbulence (in addition to wake meandering)

   -  Incorporate wake-added turbulence

   -  Improve the treatment of complex terrain (beyond specifying
      ambient wind data as NaN in VTK format)

   -  Include wakes from the nacelle and support structure

   -  Reflect wakes off of the ground.

-  Address deep-array effects for large wind farms and account for flow
   speedup around the edges of the wind farm – i.e., account for the
   wind-farm blockage effect – e.g., by mimicking the wind farm-induced
   boundary layer with surface roughness in the LES ambient wind
   precursor.

-  Implement a model to mimic the measurements taken from a LIDAR and
   other remote sensing technologies.

-  Incorporate MPI to support the modeling of large wind farms by taking
   advantage of memory parallelization and parallelization between nodes
   of an HPC.

-  Allow for a more general module form, e.g.:

   -  Support continuous states

   -  Support direct feedthrough of input to output

   -  Support full-system linearization.

-  Support an interface to Simulink for super and individual wind
   turbine controllers.

-  Implement checkpoint-restart capability.

-  Enable binary wind data input and output formats and binary
   time-series results output format.

-  Add ability to output disturbed wind in VTK format on 2D slices that
   need not be parallel to the *X-Y, Y-Z* and/or *X-Z* planes of the
   global inertial-frame coordinate system.

-  Rename the ambient wind data input files in VTK format following the
   naming convention used for the FAST.Farm-generated visualization
   output files in VTK format (with leading zeros and without the *t*).

-  Support super controller-, inflow-, and wake-related output channels
   for more than the first 9 wind turbines in the wind farm.

-  Interface FAST.Farm to the Wind-Plant Integrated System Design &
   Engineering Model
   (`WISDEM <https://github.com/NREL/WISDEM>`__\ :math:`^\text{TM}`) for
   systems-engineering applications (multidisciplinary design, analysis,
   and optimization; uncertainty quantification; and so on).

-  Develop a wrapper for stand-alone
   AeroDyn – the aerodynamics module
   of OpenFAST (or an equivalent BEM tool) – as an alternative to
   OpenFAST to support advanced performance-only wind-farm analysis that
   is much more computationally efficient than FAST.Farm analysis using
   OpenFAST.

-  Address unique offshore wind energy challenges, e.g.:

   -  Ensure consistent waves across an offshore wind farm

   -  Support the air-water interface

   -  Consider shared mooring and anchoring arrangements (for floating
      offshore wind farms).

-  Adopt the capability to support undersea marine turbine arrays (which
   may require supporting direct feedthrough of input to output to
   handle the added-mass effects).
