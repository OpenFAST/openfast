.. _glue-code:

Glue Code
=========

The OpenFAST *glue code* is the layer of software that initializes each physics
module, manages the flow of data between them, orchestrates the time-stepping
loop, and—optionally—linearizes the assembled system.  This section documents
the glue code from a user and module-developer perspective.

.. toctree::
   :maxdepth: 2

   overview
   modvar
   modglue
   solver
   linearization
