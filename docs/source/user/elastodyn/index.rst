ElastoDyn Users Guide and Theory Manual
=======================================

This document offers a quick reference guide for the ElastoDyn software
program. It is intended to be used by the general user in combination
with other OpenFAST manuals. The manual will be updated as new releases are
issued and as needed to provide further information on advancements or
modifications to the software. For general information on OpenFAST,
users should refer to section :numref:`general-reference-docs`.

.. note::

   We are in the process of migrating the documentation from FAST 8 to
   OpenFAST. For reference, various portions of old documentation are
   provided here. While most of it is still directly applicable
   to OpenFAST, portions may be out of date.


The following documents are a detailed derivation of the equations of
motion of ElastoDyn. These documents have not been compiled into a report,
so they contain mostly equations and little explanatory text. A reader
with a background in kinematics and dynamics may be able to comprehend
the equations. The documents make the most sense if studied in the following order:

1. :download:`FASTDOFs.xls <../../../OtherSupporting/ElastoDyn/FASTDOFs.xls>`:
   Contains a listing of the DOF indices used by the equations of motion in FAST.
2. :download:`FASTCoordinateSystems.doc <../../../OtherSupporting/ElastoDyn/FASTCoordinateSystems.doc>`:
   Documents the transformation matrices relating each coordinate system in FAST. Unfortunately, there are no pictures in this document that diagram these coordinate systems. They can hopefully be visualized by means of the transformation matrices.
3. :download:`FASTKinematics.doc <../../../OtherSupporting/ElastoDyn/FASTKinematics.doc>`:
   Documents the linear position, velocity, and acceleration vectors of each "important" point in the system and documents the angular velocity and acceleration vectors of each "important" reference frame in the system.  Also included is documentation of the partial velocity vectors needed by Kane's dynamics.
4. :download:`FASTKinetics.doc <../../../OtherSupporting/ElastoDyn/FASTKinetics.doc>`:
   Documents the derivation of the equations of motion using Kane's dynamics.
5. :download:`FASTLoads.doc <../../../OtherSupporting/ElastoDyn/FASTLoads.doc>`:
   Documents how the output loads are computed using terms from the equations of motion.
6. :download:`FASTMotions.doc <../../../OtherSupporting/ElastoDyn/FASTMotions.doc>`:
   Documents how the output motions are computed using variables from the equations of motion.
7. :download:`FASTLogicFlow.doc <../../../OtherSupporting/ElastoDyn/FASTLogicFlow.doc>`:
   Contains a listing of the subroutine names used by FAST. The names are listed in the order they are called within the program.

There a few minor errors in the equations documented in these papers
that may be clear after understanding the equations. The implemented code does
not have these errors. The papers do not describe the Fortran source code and
variable naming conventions, but a source code comparison is possible
with careful study.

Note that the "unofficial FAST Theory Manual" applies to the structural
equations of FAST v7 and the ElastoDyn module of FAST v8 and OpenFAST.
 
.. toctree::
   
   input.rst
