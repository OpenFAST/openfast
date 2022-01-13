
.. _AD_theory:

AeroDyn Theory
==============

This theory manual is work in progress, please refer to the AeroDyn 14 manual for more details :cite:`ad-AeroDyn:manual`. Many changes have occured since AeroDyn 14 (e.g. BEM formulation, coordinate system used in the BEM equations, dynamic stall, dynamic BEM), but these changes are not yet documented here.



Steady BEM
~~~~~~~~~~

The steady blade element momentum (BEM) equations are solved as a constrained equation, and the formulation follows the description from Ning :cite:`ad-Ning:2014`.


.. _AD_twr_shadow:

Tower shadow models
~~~~~~~~~~~~~~~~~~~

Powles tower shadow model (**TwrShadow=1**) is given by:

.. math::
   u_{TwrShadow} = - \frac{C_d}{  \sqrt{\overline{r}}  }
               \cos\left( \frac{\pi/2 \overline{y}}{\sqrt{\overline{r}}}\right)^2

where :math:`\overline{r} = \sqrt{ \overline{x}^2 + \overline{y}^2 }`.


Eames tower shadow model (**TwrShadow=2**) is given by:

.. math::
   u_{TwrShadow} = -\frac{C_d}{ TI \: \overline{x} \, \sqrt{2 \pi }  }
               \exp{\left(  -\frac{1}{2}  \left(\frac{ \overline{y}}{ TI \: \overline{x} } \right)^2 \right) }

where :math:`TI` is the turbulence intensity at the tower node. 


.. _AD_buoyancy:

Buoyancy
~~~~~~~~

When a solid object is submerged in a fluid, it experiences a net force, buoyancy, 
from the hydrostatic fluid pressure acting on its surface. This force can often 
be neglected in less dense fluids, such as air, but can be significant in denser
fluids, such as water. To capture the effects of this force on MHK turbines, 
buoyant loads are calculated for the turbine blades, tower, hub, and nacelle. 
Marine growth is neglected for all components. :numref:`AD_buoy_coords` - 
:numref:`AD_buoy_hubnacelle` detail the coordinate systems and blade, tower, hub,
and nacelle buoyancy calculations.

.. _AD_buoy_coords:

Coordinate Systems
------------------
The buoyant force acting on an element depends on its instantaneous orientation 
and depth. The orientation is defined by heading and inclination angles, which
are calculated for each element at every time step. Total water depth is defined
by the user, relative to the mean sea level. The instantaneous depth of each
element is based on the total water depth and element position and is calculated
at each time step.

.. _AD_buoy_bladestower:

Blades and Tower
----------------
To allow for an efficient analytical solution, the blades and tower are modeled
as tapered cylinders. The cross-sectional area of the tapered cylinders is set
equal to the blade or tower cross-sectional area. Loads are estimated by breaking
the blade or tower into elements of a given length and integrating the hydrostatic
pressure over the wetted area of each element. For the blades, loads are applied
at a user-specified center of buoyancy. For the tower, loads are applied at the
centerline. When applicable, end effects are accounted for by calculating the
fluid pressure on the exposed axial face of the element. The tower is assumed to
be embedded into the seabed, such that no end effects at the tower base are needed.

The buoyancy calculation for the blades and tower is completed according to the following steps:

1.	Calculate parameters related to element geometry that do not change with time
2.	Check that no elements cross the free surface or go beneath the seabed
3.	Calculate the instantaneous orientation and depth of each element
4.	Integrate hydrostatic fluid pressure over the wetted surface of each element and express as a force acting at the center of buoyancy
5.	For blades, calculate the buoyant force on the axial face of the blade root and tip; add the tip force to the adjacent element and store the root force
6.	For the tower, calculate and store the buoyant force on the axial face of the tower top
7.	Move buoyant loads from the center of buoyancy to the aerodynamic center
8.	Express buoyant loads in the form expected by OpenFAST
9.	Add buoyant loads to aerodynamic loads

.. _AD_buoy_hubnacelle:

Hub and Nacelle
---------------
The hub and nacelle are treated as separate components. The buoyant force is 
determined by the volume of either the hub or nacelle and applied at its 
user-specified center of buoyancy. Corrections are made to account for the joints
between the hub and blades and the nacelle and tower, as the joint locations are
not exposed to fluid pressure. No correction is made for the joint between the 
hub and nacelle.

The buoyancy calculation for the hub and nacelle is completed according to the following steps:

1.	Check that the component does not cross the free surface or go beneath the seabed
2.	Calculate the instantaneous depth of the component
3.	Calculate the buoyant force from the volume of the component
4.	Move buoyant loads from the center of buoyancy to the aerodynamic center
5.	For the hub, correct loads to account for the joints with each blade
6.	For the nacelle, correct loads to account for the joint with the tower
