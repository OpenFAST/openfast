
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

OpenFAST 5.0
~~~~~~~~~~~~

With the release of OpenFAST 5.0, significant changes were made to the input files of ExtPtfm and its internal formulation. 
Couplings between ExtPtfm and HydroDyn and between ExtPtfm and the mooring modules have been enabled to more conveniently 
simulate a flexible floating structure. Finally, linearization with ExtPtfm is now possible.

When simulating a floating structure in ExtPtfm, HydroDyn can be enabled in the glue-code to compute and provide 
potential-flow wave excitation and radiation loads and hydrostatic loads for both rigid-body modes and additional 
generalized degrees of freedom through the `NAddDOF` option (see HydroDyn documentation). With this setup, we need 
the number of active modes in ExtPtfm to be equal to :math:`6+NAddDOF` in HydroDyn, with the first 6 modes being 
rigid-body modes, followed by a number of elastic modes. 

Note that if HydroDyn also contains strip-theory members, the resulting strip-theory loads will all be mapped to 
the ExtPtfm rigid-body modes only, not the elastic modes. Therefore, users are discouraged from using strip-theory-only 
members in HydroDyn. The same applies to marine growth and ballast/flooding in HydroDyn. However, hybrid members 
with drag only can still be used in many cases to help obtain more relastic floater motions. The drag force tends 
to be small and does not have a strong direct impact on the structural loading. In this case, mapping the drag force 
to the rigid-body modes only can be an acceptable simplification. The same consideration also applies to the 
second-order potential-flow options in HydroDyn that do no support generalized degrees of freedom, yet. 

When coupling to mooring, a set of connection points must be defined in the ExtPtfm module. These connections points 
are coupled to the fairleads in the mooring models. Generally, the connection points in ExtPtfm should be the same 
set of points used as fairleads in the mooring model; however, they need not be exactly the same as OpenFAST automatically 
introduces nearest-neighbor mapping between the two assuming rigid connections.

Note that the above workflow with ExtPtfm coupling to HydroDyn and mooring is still a work-in-progress. 


