.. _other_docs:

Other documentation
~~~~~~~~~~~~~~~~~~~
Additional documentation exists that may be useful for developers seeking deeper
understanding of the solver and mathematics.  This documentation is not generally
necessary for most development efforts.

- :download:`DCM_Interpolation.pdf    <../../OtherSupporting/DCM_Interpolation/DCM_Interpolation.pdf>`
   This is a summary of the mathematics used in the interpolation of 
   DCM (direction cosine matrices) using logarithmic mapping and matrix exponentials.
- :download:`OpenFAST_Algorithms.pdf  <../../OtherSupporting/OpenFAST_Algorithms/OpenFAST_Algorithms.pdf>`
   This is a summary of the solve method used in the glue code.
- :download:`OutListParameters.xlsx   <../../OtherSupporting/OutListParameters.xlsx>`
   This Excel file contains the full list of outputs for each module.  It is used to generate the
   Fortran code for the output channel list handling for each module (this code is generally in
   the _IO.f90 files).  The MATLAB script available in the
   `matlab-toolbox <https://github.com/OpenFAST/matlab-toolbox>`__ repository at *Utilities/GetOutListParameters.m*.
