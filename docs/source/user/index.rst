.. _user_guide:

User Documentation
==================

We are in the process of transitioning legacy FAST v8 documentation, which can be found at https://www.nrel.gov/wind/nwtc.html.

.. note::

    Much of the documentation here is legacy documentation from FAST v8. While most of it is still
    directly applicable to OpenFAST, portions may be out of date.


.. _general-reference-docs:

General
~~~~~~~
.. toctree::
   :maxdepth: 1

   fast_to_openfast.rst
   api_change.rst
   input_file_overview.rst

Workshop material, legacy documentation, and other resources are listed below.

- `Workshop Presentations <https://drive.google.com/drive/folders/1BDDfcnIyvmZCwf7eFo0ISI7aF_FMAOvt?usp=sharing>`_
- :download:`Old FAST v6 User’s Guide <../../OtherSupporting/Old_FAST6_UsersGuide.pdf>`
- :download:`FAST v8 README <../../OtherSupporting/FAST8_README.pdf>`
- `Implementation of Substructure Flexibility and Member-Level Load Capabilities for Floating Offshore Wind Turbines in OpenFAST <https://www.nrel.gov/docs/fy20osti/76822.pdf>`_
- `FAST modularization framework for wind turbine simulation: full-system linearization <https://www.nrel.gov/docs/fy17osti/67015.pdf>`_
- `Full-System Linearization for Floating Offshore Wind Turbines in OpenFAST <https://www.nrel.gov/docs/fy19osti/71865.pdf>`_ 
- :download:`FAST with Labview <../../OtherSupporting/UsingFAST4Labview.pdf>`
- :download:`OutListParameters.xlsx <../../OtherSupporting/OutListParameters.xlsx>` - Contains the full list of outputs for each module.


Module Documentation
~~~~~~~~~~~~~~~~~~~~
This section contains documentation for the OpenFAST module-coupling environment and its underlying modules.
Documentation covers usage of models, underlying theory, and in some cases module verification.

.. toctree::
   :maxdepth: 1

   AeroDyn <aerodyn/index.rst>
   OLAF <aerodyn-olaf/index.rst>
   Aeroacoustics <aerodyn-aeroacoustics/index.rst>
   BeamDyn <beamdyn/index.rst>
   SubDyn <subdyn/index.rst>
   ExtPtfm <extptfm//index.rst>
   ElastoDyn <elastodyn/index.rst>
   HydroDyn <hydrodyn/index.rst>
   InflowWind <inflowwind/index.rst>
   ServoDyn <servodyn/index.rst>
   Structural Control <servodyn-stc/StC_index.rst>
   C++ API <cppapi/index.rst>
   FAST.Farm <fast.farm/index.rst>
   

Modularization Framework
~~~~~~~~~~~~~~~~~~~~~~~~

Information specific to the modularization framework of OpenFAST is provided here. These are a collection
of publications, presentations, and past studies on the subject.

- `The New Modularization Framework for the FAST Wind Turbine CAE Tool <https://www.nrel.gov/docs/fy13osti/57228.pdf>`_
- :download:`Example Module Implementation Plans <../../OtherSupporting/ModulePlan_GasmiPaperExamples.doc>`
- :download:`Module and Mesh-Mapping Linearization Implementation Plan <../../OtherSupporting/LinearizationOfMeshMapping_Rev18_Rev2.doc>`
- :download:`Interpolation of DCMs <../../OtherSupporting/DCM_Interpolation/DCM_Interpolation.pdf>` - A summary of the mathematics used in the interpolation of DCM (direction cosine matrices) using logarithmic mapping and matrix exponentials.
- :download:`Set-point Linearization Development Plan <../../OtherSupporting/DevelopmentPlan-SetPoint-Linearization.pdf>`

.. - :download:`OpenFAST Steady State Solution <../../OtherSupporting/OpenFASTSteadyStateSolution_Rev7.doc>`


Glue Code and Mesh Mapping
~~~~~~~~~~~~~~~~~~~~~~~~~~

- `FAST Modular Wind Turbine CAE Tool: Nonmatching Spatial and Temporal Meshes <https://www.nrel.gov/docs/fy14osti/60742.pdf>`_
- `FAST Modular Framework for Wind Turbine Simulation: New Algorithms and Numerical Examples <https://dx.doi.org/10.2514/6.2015-1461>`_
- :download:`OpenFAST Algorithms <../../OtherSupporting/OpenFAST_Algorithms/OpenFAST_Algorithms.pdf>` - A summary of the solve method used in the glue code.
- :download:`Predictor-Corrector Approach <../../OtherSupporting/ProposedPCApproach_Rev4.docx>`
