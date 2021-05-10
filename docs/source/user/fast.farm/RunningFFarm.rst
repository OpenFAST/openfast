.. _FF:Running:

Running FAST.Farm
=================

As FAST.Farm is a module of OpenFAST, the process of downloading, compiling,
and running FAST.Farm is the same as that for OpenFAST. Such instructions are
available in the :ref:`installation` documentation.

.. note::
   To improve the speed of OpenFAST compiled with FAST.Farm enabled, the user
   may wish to compile in single precision with `OpenMP`.  To do so, add the
   `-DDOUBLE_PRECISION:BOOL=OFF -DOPENMP=ON` options with CMake.

.. note::
   Checkpoint-restart capability has not yet been implemented within FAST.Farm.
