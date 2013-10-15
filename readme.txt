HydroDyn v2.00.01a-gjh  02-Oct-2013
  

Version 2.00 has been updated to work within the FAST Modularization Framework. 
 
The major changes in version 2.00 are:

*    Version 2.00 integrates with FAST v8.XX.
*    Version 2.00 does not integrate with FAST v7.xx.
*    Mooring lines have been removed from HydroDyn. The mooring line modeling is now being handled in a different FAST module (MAP). You can 
       however use the additional preload, stiffness and damping matrices to model a mooring line system within HydroDyn.
*    Multi-member structures can now be modeled within HydroDyn (applicable e.g. to fixed-bottom tripod and jacket structures, as well as 
       thin members (e.g., braces/spokes) of floating platforms). This modeling implement's Morison's equation for multiple slender cylinders, 
       inclined and tapered members, inertia/drag/added mass/buoyancy/dynamic pressure loads, filled/flooded effects, and marine growth effects.
*    Hybrid models with both linear radiation/diffraction calculations and Morison Equation calculations are possible.
*    Added a linear state-space-based radiation formulation to be used together with SS_Fitting.
*    The new implementation has well-defined data exchange interfaces (following the FAST modularization framework) which should make integration 
       of HydroDyn into other multi-physics software packages much simpler. 

------------- Known Issues -------------
* Heave Coefs are zero at marine growth boundary unless user inputs a joint at that boundary.  
* Check whether FSLoc can be set to DEFAULT in which case it is set to the value of MSL2SWL.
* WaveStMod must be equal to 0, no wave stretching
* MSL2SWL must be set to 0
* OutAll must be set to False
* JointOvlp switch must be set to 0, no super members will be generated
* Floating Platform DOF flags must be TRUE
* Input data error checking is incomplete
* If a member output falls onto a joint which was created by the software, due to element splitting rules, 
    only the distributed load from one side of the joint will be reported, not a sum of all distributed loads at that joint.
* Each joint ID must appear at least once in the Members Table
* Driver File: Morison section is currently ignored

------------- Future Work  -------------

* Further verify through data from the IEA Wind Task 23/30 OC3/OC4 projects
* Further validate with test data
* Write the HydroDyn user and theory manual 
* Add threshold to A, K, C to minimize computations
* Add support for wave stretching
* Add support for offset between mean sea level and still water level (MSL2SWL)
* Add implementation of OutAll which will generate a complete set of outputs
* Add full support for Floating platform DOF flags
* Add ability to report the added mass forces.  Currently, HydroDyn doesn't have access to structural accelerations.
* Add super member element support
                * Add check to make sure requested member output locations are not within a super member element
* Add time-varying Directional Cosine matrices for the meshes
* Move any time-step array allocations to OtherStates so there are no allocates on a per time step basis
* See if using vector and matrix library routines instead of DO loops speeds-up code
* Add final LumpedMesh and DistribMesh information to the Summary file outputs
* Need to implement additional input data checking
* Manually DEALLOCATE all arrays which were allocated via ALLOCATE
* Source code needs additional comments and formatting for clarity
* Need to consider using double precision for any length/distance and interpolation factor calculations
* Add second-order hydrodynamic loads
* Add wave directional spreading
* Add ability to print multiple input file errors/mistakes instead of returning/failing upon the first encountered error.