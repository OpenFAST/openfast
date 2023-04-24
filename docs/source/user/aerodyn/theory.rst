
.. _AD_theory:

AeroDyn Theory
==============

This theory manual is work in progress, please refer to the AeroDyn 14 manual for more details :cite:`ad-AeroDyn:manual`. Many changes have occured since AeroDyn 14 (e.g. BEM formulation, coordinate system used in the BEM equations, dynamic stall, dynamic BEM), but these changes are not yet documented here.



Steady BEM
~~~~~~~~~~

The steady blade element momentum (BEM) equations are solved as a constrained equation, and the formulation follows the description from Ning :cite:`ad-Ning:2014`.



.. _AD_DBEMT:

Dynamic BEM Theory (DBEMT)
~~~~~~~~~~~~~~~~~~~~~~~~~~



Two equivalent versions of Oye's dynamic inflow model are implemented in AeroDyn.
The first one uses discrete time, it can be used with the constant-tau1 model 
(``DBEMT_Mod=1``) or the varying-tau1 model (``DBEMT_Mod=2``), but it cannot be used for linearization.
The second version uses a continuous-time state-space formulation  (``DBEMT_Mod=1``), it assumes a constant-tau1, and can be used for linearization.  
For a same value of :math:`\tau_1`, the discrete-time and continuous-time formulations returns exactly the same results.





Oye's dynamic inflow model consists of two first-order differential equations (see :cite:`ad-Branlard:book`):

.. math::
   \begin{aligned}
       \boldsymbol{W}_\text{int}+\tau_1    \boldsymbol{\dot{W}}_\text{int} &= \boldsymbol{W}_\text{qs} + k \tau_1 \boldsymbol{\dot{W}}_\text{qs} \\
       \boldsymbol{W}+\tau_2 \boldsymbol{\dot{W}} &= \boldsymbol{W}_\text{int}
   \end{aligned}

where 
:math:`\boldsymbol{W}` is the dynamic induction vector at the rotor (at a given blade position and radial position), 
:math:`\boldsymbol{W}_\text{qs}` is the quasi-steady induction, 
:math:`\boldsymbol{W}_\text{int}` is an intermediate value coupling the quasi-steady and the actual inductions (may be discontinuous if the quasi-steady indution is discontinuous).
and
:math:`(\dot{\ })` represents the time derivative.
The coupling constant :math:`k`, with values between 0 and 1, is usually chosen as :math:`k=0.6`.
Oye's dynamic inflow model relies on two time constants, :math:`\tau_1` and :math:`\tau_2` :

.. math::
        \tau_1=\frac{1.1}{1-1.3 \min(\overline{a},0.5)} \frac{R}{\overline{U}_0}
        , \qquad
        \tau_2 =\left[ 0.39-0.26\left(\frac{r}{R}\right)^2\right] \tau_1

where :math:`R` is the rotor radius, :math:`\overline{U}_0` is the average wind speed over the rotor, :math:`\overline{a}` is the average axial induction over the rotor, and :math:`r` is the radial position along the blade.
For ``DBEMT_Mod=1`` or ``DBEMT_Mod=3``, the user needs to provide the value of :math:`\tau_1`.




The continuous-time state-space formulation of the dynamic inflow model (``DBEMT_Mod=3``) was derived in :cite:`ad-Branlard:2022`.

.. math::
   \begin{align}
      \begin{bmatrix}
      \boldsymbol{\dot{W}}_\text{red}\\
      \boldsymbol{\dot{W}}\\
      \end{bmatrix}
      =
      \begin{bmatrix}
      -\frac{1}{\tau_1}\boldsymbol{I}_2 & \boldsymbol{0} \\
       \frac{1}{\tau_2}\boldsymbol{I}_2 &
      -\frac{1}{\tau_2}\boldsymbol{I}_2 \\
      \end{bmatrix}
      \begin{bmatrix}
      \boldsymbol{W}_\text{red}\\
      \boldsymbol{W}\\
      \end{bmatrix}
      +
      \begin{bmatrix}
       \frac{1-k}{\tau_1} \\
       \frac{k}{\tau_2}\\
      \end{bmatrix}
     \boldsymbol{W}_\text{qs}
   \end{align}

where 
:math:`\boldsymbol{I}_2` is the 2x2 identity matrix,
:math:`\boldsymbol{W}_\text{red}` is the reduced induction which is a continuous, scaled, and lagged version of the quasi-steady induction, defined as:

.. math::
    \boldsymbol{W}_\text{int} = \boldsymbol{W}_\text{red} + k \boldsymbol{W}_\text{qs} 


The discrete-time version of the model is documented in the unpublished manual of DBEMT.
The current discrete-time formulation is complex and in the future it can be simplified by using :math:`\boldsymbol{W}_\text{red}`.






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
by the user, relative to the still water level (or relative to the mean sea 
level when running AeroDyn in standalone mode with the AeroDyn driver). The
instantaneous depth of each element is based on its position in global coordinates
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
be either embedded into the seabed or attached to another support structure member,
such that no end effects at the tower base are needed. For MHK turbines with a support
structure (i.e., any structure other than a simple tower embedded in the seabed), it
is currently recommended to model the entire support structure, including the tower, in HydroDyn.
Future releases will include the ability to neglect fluid loads at the interface between a tower
modeled in AeroDyn and a platform modeled in HydroDyn.


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

Although the blade and tower buoyant loads are not based on volume, the volumes of these components are
written to the AeroDyn summary file for reference. The blade and tower volumes are calculated by summing
the volume of each element, assumed to be a tapered cylinder. The volume of a single element :math:`(V_{elem})`
is given by:

.. math::
   V_{elem} = \frac{\pi}{3} (r_i^2 + r_i r_{i+1} + r_{i+1}^2) dl

where :math:`r_i` is the element radius at node :math:`i`, :math:`r_{i+1}` is the element radius at node 
:math:`i+1`, and :math:`dl` is the element length.


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
