
.. _AD_theory:

AeroDynTheory
=============

This theory manual is work in progress, please refer to the AeroDyn manual for more details. 


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




.. _AD_UA:

Unsteady aerodynamics
---------------------

Beddoes-Leishman type models (UAMod=2,3)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Beddoes-Leishman type dynamic stall models are currently described in the document: 
The Unsteady Aerodynamics Module for FAST 8, from  Rick Damiani and Greg Hayman (2017)


Beddoes-Leishman state space models (UAMod=4,5)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TODO


Oye model (UAMod=6)
~~~~~~~~~~~~~~~~~~~


See Hansen book

Boeing-Vertol model (UAMod=7)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Boeing-Vertol is used mentioned in the following paper: The Development of CACTUS, a Wind and Marine Turbine Performance Simulation Code from Jonathan C. Murray  and Matthew Barone (2011).
The documentation presented here was inspired from the implementation done in the vortex code CACTUS.



.. math::

   \alpha_{dyn} = \alpha_qs - k_1 \gamma \sqrt{\left| \frac{c\dot{\alpha}}{2U}\right|}







