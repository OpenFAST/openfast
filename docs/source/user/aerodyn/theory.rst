
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
   \begin{align}
       \boldsymbol{W}_\text{int}+\tau_1    \boldsymbol{\dot{W}}_\text{int}
           &=
       \boldsymbol{W}_\text{qs} + k \tau_1 \boldsymbol{\dot{W}}_\text{qs} \\
       \boldsymbol{W}+\tau_2 \boldsymbol{\dot{W}}
           &=
       \boldsymbol{W}_\text{int}
   \end{align}

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




