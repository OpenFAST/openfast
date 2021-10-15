
.. _ep-usage:



Usage
-----

The *ExtPtfm* module uses superelement properties of the support structure
provided by the user (e.g., mass, stiffness, damping, and time series of excitation forces)
to compute the reaction of the support-structure at the interface.
The module uses a linear formulation and internally keeps track
of Guyan and Craig-Bampton degrees of freedom.


Typical sequentially coupled workflow with *ExtPtfm*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The overall workflow includes the following steps:

-  The substructure designer performs a time-domain simulation of the
   isolated substructure under a given sea state, using an external 
   tool such as a finite-element tool that has the capability to generate
   a superelement model. The underlying model behind the external tool 
   and the time series of loads are reduced using the CB technique described in :ref:`ep-theory`,
   where the leader DOF are selected as the ones at the substructure
   interface node. Results from the reduction are written to a file
   containing the reduced system matrices, :math:`\boldsymbol{M}_r`,
   :math:`\boldsymbol{C}_r`, :math:`\boldsymbol{K}_r`, and the time
   series of reduced loads, :math:`\boldsymbol{f}_r`.

-  The file is imported in *OpenFAST* by the *ExtPtfm* module, and a
   time-domain simulation of the full wind turbine is run with the
   reduced representation of the substructure. At every time step, the
   *ExtPtfm* module receives as inputs the displacement, velocity, and acceleration
   of the interface point and returns the loads at this point.

-  *OpenFAST* exports times series of loads and displacements at the
   interface, which are then returned to the substructure designer.
   These inputs are used as boundary conditions to the external 
   tool and then another time-domain simulation of the substructure is
   run. Stress concentrations are computed, and code checks are
   performed.


Using the module
~~~~~~~~~~~~~~~~

The ``ExtPtfm`` module is activated by setting the flag ``CompSub`` to
:math:`2` in the main OpenFAST input file. The variable ``SubFile`` in
this same file needs to be set to a valid input file for the module (see
:ref:`ep_input-files` for the input file specifications).

.. code::

   ---  [...]
     2  CompSub - Compute sub-structural dynamics (switch) {0=None; 1=SubDyn; 2=ExtPtfm}
   ---  [...]
   "Turbine_ExtPtfm.dat" SubFile - Sub-structure input file (SubDyn or ExtPtfm)



Recommendations
~~~~~~~~~~~~~~~

**Time step** superelements may contain high frequencies. 
As a rule of thumb, it is recommended to use a time step 
satisfying the following criteria: :math:`\Delta t=\frac{1}{10 f_\text{max}}`, where :math:`f_\text{max}` is the maximum frequency present in the superelement or full OpenFAST model.


**Number of Craig-Bampton modes** It is recommended to perform a sensitivity analyses of the results
for an increasing number of Craig-Bampton modes to get an idea of how many modes are needed to reach convergence 
and see which modes contributes the most to the sytem response.





