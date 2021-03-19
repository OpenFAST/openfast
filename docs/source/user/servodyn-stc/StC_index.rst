.. _StC:

Structural Control (SrvD)
=========================

.. only:: html

   The Structural Control module within ServoDyn is used to simulate tuned mass
   dampers and tuned liquid column dampers. This module also provides an option
   for applying a force time series load at the location of the Structural
   Control node.

   The location of the StC node is given in the main ServoDyn input file (see
   :numref:`SrvD-Stc-inputs` and :numref:`StC-Locations`).  These may be mounted
   at the nacelle, tower, blade, or platform.  Multiple StC's may be placed at
   each location, each with it's own input file.  Output channels from the StC
   module are handled by ServoDyn (see :numref:`SrvD-Outputs` for details).



.. toctree::
   :maxdepth: 2

   StC_input.rst
   StC_Theory.rst
   StC_TLCD_Theory.rst
   zrefs.rst

