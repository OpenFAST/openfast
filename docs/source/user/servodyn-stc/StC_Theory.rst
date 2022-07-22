.. _StC-Theory:

==========================================================
Theory Manual for the Tuned Mass Damper Module in OpenFAST
==========================================================

:Author: William La Cava & Matthew A. Lackner
   Department of Mechanical and Industrial Engineering
   University of Massachusetts Amherst
   Amherst, MA 01003
   ``wlacava@umass.edu``, ``lackner@ecs.umass.edu``

This document was edited by Jason M. Jonkman of NREL
to include an independent vertically oriented TMD in OpenFAST.
``jason.jonkman@nrel.gov``


This manual describes updated functionality in OpenFAST that simulates
the addition of tuned mass dampers (TMDs) for structural control. The
dampers can be added to the blades, nacelle, tower, or substructure. For
application studies of these systems, refer to
:cite:`stc-lackner_passive_2011,stc-lackner_structural_2011,stc-namik_active_2013,stc-stewart_effect_2011,stc-stewart_impact_2014,stc-stewart_optimization_2013`.
The TMDs are three independent, 1 DOF, linear mass spring damping
elements that act in the local :math:`x`, :math:`y`, and :math:`z`
coordinate systems of each component. The other functionality of the
structural control (StC) module, including an omnidirectional TMD and
TLCD are not documented herein. We first present the theoretical
background and then describe the code changes.

Theoretical Background
======================

Definitions
-----------

.. container::
   :name: tab:defs

   .. table:: Definitions

      +-----------------+--------------------+
      | Variable        | Description        |
      +=================+====================+
      | |O_eq|          | |O_desc|           |
      +-----------------+--------------------+
      | |P_eq|          | |P_desc|           |
      +-----------------+--------------------+
      | |TMD_eq|        | |TMD_desc|         |
      +-----------------+--------------------+
      | |G_eq|          | |G_desc|           |
      +-----------------+--------------------+
      | |N_eq|          | |N_desc|           |
      +-----------------+--------------------+
      | |TMD_OG_eq|     | |TMD_OG_desc|      |
      +-----------------+--------------------+
      | |TMD_PN_eq|     | |TMD_PN_desc|      |
      +-----------------+--------------------+
      | |TMD_X_eq|      | |TMD_X_desc|       |
      +-----------------+--------------------+
      | |TMD_Y_eq|      | |TMD_Y_desc|       |
      +-----------------+--------------------+
      | |TMD_Z_eq|      | |TMD_Z_desc|       |
      +-----------------+--------------------+
      | |P_OG_eq|       | |P_OG_desc|        |
      +-----------------+--------------------+
      | |R_OG_eq|       | |R_OG_desc|        |
      +-----------------+--------------------+
      | |R_GN_eq|       | |R_GN_desc|        |
      +-----------------+--------------------+
      | |Omega_NON_eq|  | |Omega_NON_desc|   |
      +-----------------+--------------------+
      | |OmegaD_NON_eq| | |OmegaD_NON_desc|  |
      +-----------------+--------------------+
      | |a_GOG_eq|      | |a_GOG_desc|       |
      +-----------------+--------------------+
      | |a_GON_eq|      | |a_GON_desc|       |
      +-----------------+--------------------+

.. |O_eq|            replace:: :math:`O`
.. |O_desc|          replace:: origin point of global inertial reference frame
.. |P_eq|            replace:: :math:`P`
.. |P_desc|          replace:: origin point of non-inertial reference frame fixed to  component (blade, nacelle, tower, substructure) where TMDs are at rest
.. |TMD_eq|          replace:: :math:`TMD`
.. |TMD_desc|        replace:: location point of a TMD
.. |G_eq|            replace:: :math:`G`
.. |G_desc|          replace:: axis orientation of global reference frame
.. |N_eq|            replace:: :math:`N`
.. |N_desc|          replace:: axis orientation of local component reference frame with  unit vectors :math:`\hat{\imath}, \hat{\jmath}, \hat{k}`
.. |TMD_OG_eq|       replace:: :math:`\vec{r}_{_{_{TMD/O_G}}} = \left[ \begin{array}{c} x \\ y\\ z \end{array} \right]_{_{TMD/O_G}}`
.. |TMD_OG_desc|     replace::  position of a TMD with respect to (w.r.t.) :math:`O` with orientation :math:`G`
.. |TMD_PN_eq|       replace:: :math:`\vec{r}_{_{_{TMD/P_N}}} = \left[ \begin{array}{c} x \\ y\\ z \end{array} \right]_{_{TMD/P_N}}`
.. |TMD_PN_desc|     replace:: position of a TMD w.r.t. :math:`P_N`
.. |TMD_X_eq|        replace:: :math:`\vec{r}_{_{_{TMD_X}}}`
.. |TMD_X_desc|      replace:: position vector for :math:`TMD_X`
.. |TMD_Y_eq|        replace:: :math:`\vec{r}_{_{_{TMD_Y}}}`
.. |TMD_Y_desc|      replace:: position vector for :math:`TMD_Y`
.. |TMD_Z_eq|        replace:: :math:`\vec{r}_{_{_{TMD_Z}}}`
.. |TMD_Z_desc|      replace:: position vector for :math:`TMD_Z`
.. |P_OG_eq|         replace:: :math:`\vec{r}_{_{P/O_G}} =\left[ \begin{array}{c} x \\ y\\ z \end{array} \right]_{_{P/O_G}}`
.. |P_OG_desc|       replace:: position vector of component w.r.t. :math:`O_G`
.. |R_OG_eq|         replace:: :math:`R_{_{N/G}}`
.. |R_OG_desc|       replace:: 3 x 3 rotation matrix transforming orientation :math:`G` to :math:`N`
.. |R_GN_eq|         replace:: :math:`R_{_{G/N}} = R_{_{N/G}}^T`
.. |R_GN_desc|       replace:: transformation from :math:`N` to :math:`G` 
.. |Omega_NON_eq|    replace:: :math:`\vec{\omega}_{_{N/O_N}} = \dot{\left[ \begin{array}{c} \theta \\ \phi \\ \psi \end{array} \right]}_{_{N/O_N}}`
.. |Omega_NON_desc|  replace:: angular velocity of component in orientation :math:`N`; defined likewise for  :math:`G`
.. |OmegaD_NON_eq|   replace:: :math:`\dot{\vec{\omega}}_{_{N/O_N}} = \vec{\alpha}_{_{N/O_N}}`
.. |OmegaD_NON_desc| replace:: angular acceleration of component
.. |a_GOG_eq|        replace:: :math:`\vec{a}_{G/O_G} = \left[ \begin{array}{c}0 \\ 0\\ -g \end{array} \right]_{/O_G}`
.. |a_GOG_desc|      replace:: gravitational acceleration in global coordinates
.. |a_GON_eq|        replace:: :math:`\vec{a}_{G/O_N} = R_{_{N/G}} \vec{a}_{G/O_G} = \left[ \begin{array}{c}a_{_{G_X}} \\ a_{_{G_Y}}\\ a_{_{G_Z}} \end{array} \right]_{/O_N}`
.. |a_GON_desc|      replace:: gravity w.r.t. :math:`O_N`


Equations of motion
-------------------

The position vectors of the TMDs in the two reference frames :math:`O`
and :math:`P` are related by

.. math:: \vec{r}_{_{TMD/O_G}} =  \vec{r}_{_{P/O_G}} +  \vec{r}_{_{TMD/P_G}}

Expressed in orientation :math:`N`,

.. math:: \vec{r}_{_{TMD/O_N}} =  \vec{r}_{_{P/O_N}} +  \vec{r}_{_{TMD/P_N}}

.. math:: \Rightarrow \vec{r}_{_{TMD/P_N}} =  \vec{r}_{_{TMD/O_N}} -  \vec{r}_{_{P/O_N}}

Differentiating, [1]_

.. math::
   \dot{\vec{r}}_{_{TMD/P_N}}= \dot{\vec{r}}_{_{TMD/O_N}} 
      - \dot{\vec{r}}_{_{P/O_N}} 
      - \vec{\omega}_{_{N/O_N}} \times \vec{r}_{_{TMD/P_N}}

differentiating again gives the acceleration of the TMD w.r.t. :math:`P`
(the nacelle position), oriented with :math:`N`:

.. math::
   \begin{array}{cc}
      \ddot{\vec{r}}_{_{TMD/P_N}} =
               &  \ddot{\vec{r}}_{_{TMD/O_N}}
                  - \ddot{\vec{r}}_{_{P/O_N}} - \vec{\omega}_{_{N/O_N}} 
                  \times (\vec{\omega}_{_{N/O_N}} \times \vec{r}_{_{TMD/P_N}}) \\[1.1em] 
               &- \vec{\alpha}_{_{N/O_N}} \times \vec{r}_{_{TMD/P_N}} 
                  - 2 \vec{\omega}_{_{N/O_N}} \times \dot{\vec{r}}_{_{TMD/P_N}}
   \end{array}
   :label: accel

The right-hand side contains the following terms:

.. container::
   :name: tab:

   .. table:: RHS terms 

      +--------------------+-----------------------+
      | |Rddot_TMD_ON_eq|  | |Rddot_TMD_ON_desc|   |
      +--------------------+-----------------------+
      | |Rddot_P_ON_eq|    | |Rddot_P_ON_desc|     |
      +--------------------+-----------------------+
      | |Omega_N_ON_eq|    | |Omega_N_ON_desc|     |
      +--------------------+-----------------------+
      | |CentripAcc_eq|    | |CentripAcc_desc|     |
      +--------------------+-----------------------+
      | |TangentAcc_eq|    | |TangentAcc_desc|     |
      +--------------------+-----------------------+
      | |Coriolus_eq|      | |Coriolus_desc|       |
      +--------------------+-----------------------+

.. |Rddot_TMD_ON_eq|   replace:: :math:`\ddot{\vec{r}}_{_{TMD/O_N}}`
.. |Rddot_TMD_ON_desc| replace:: acceleration of the TMD in the *inertial* frame :math:`O_N`
.. |Rddot_P_ON_eq|   replace:: :math:`\ddot{\vec{r}}_{_{P/O_N}} = R_{_{N/G}} \ddot{\vec{r}}_{_{P/O_G}}`
.. |Rddot_P_ON_desc| replace:: acceleration of the Nacelle origin :math:`P` w.r.t. :math:`O_N`
.. |Omega_N_ON_eq|   replace:: :math:`\vec{\omega}_{_{N/O_N}} = R_{_{N/G}} \vec{\omega}_{_{N/O_G}}`
.. |Omega_N_ON_desc| replace:: angular velocity of nacelle w.r.t. :math:`O_N`
.. |CentripAcc_eq|   replace:: :math:`\vec{\omega}_{_{N/O_N}} \times (\vec{\omega}_{_{N/O_N}} \times \vec{r}_{_{TMD/P_N}})`
.. |CentripAcc_desc| replace:: Centripetal acceleration
.. |TangentAcc_eq|   replace:: :math:`\vec{\alpha}_{_{N/O_N}} \times \vec{r}_{_{TMD/P_N}}`
.. |TangentAcc_desc| replace:: Tangential acceleration
.. |Coriolus_eq|   replace:: :math:`2\vec{\omega}_{_{N/O_N}} \times \dot{\vec{r}}_{_{TMD/P_N}}`
.. |Coriolus_desc| replace:: Coriolis acceleration


The acceleration in the inertial frame
:math:`\ddot{\vec{r}}_{_{TMD/O_N}}` can be replaced with a force balance

.. math::
   \begin{aligned}
      \ddot{\vec{r}}_{_{TMD/O_N}} = \left[ 
         \begin{array}{c} \ddot{x} \\
            \ddot{y} \\
            \ddot{z}
         \end{array}
      \right]_{_{TMD/O_N}} = \frac{1}{m} \left[ 
         \begin{array}{c} 
            \sum{F_X} \\
            \sum{F_Y} \\
            \sum{F_Z} 
         \end{array}
      \right]_{_{TMD/O_N}} = \frac{1}{m} \vec{F}_{_{TMD/O_N}}
    \end{aligned}

Substituting the force balance into Equation :eq:`accel` gives
the general equation of motion for a TMD:

.. math::
   \begin{array}{cc}
      \ddot{\vec{r}}_{_{TMD/P_N}} = & \frac{1}{m} \vec{F}_{_{TMD/O_N}}
         - \ddot{\vec{r}}_{_{P/O_N}}
         - \vec{\omega}_{_{N/O_N}} \times (\vec{\omega}_{_{N/O_N}}
               \times \vec{r}_{_{TMD/P_N}}) \\[1.1em]
      & - \vec{\alpha}_{_{N/O_N}} \times \vec{r}_{_{TMD/P_N}}
         - 2 \vec{\omega}_{_{N/O_N}} \times \dot{\vec{r}}_{_{TMD/P_N}}
   \end{array}
   :label: EOM

We will now solve the equations of motion for :math:`TMD_X`,
:math:`TMD_Y`, and :math:`TMD_Z`.

TMD_X :
~~~~~~~

The external forces :math:`\vec{F}_{_{TMD_X/O_N}}` are given by

.. math::
   \vec{F}_{_{TMD_X/O_N}} = \left[
      \begin{array}{c}
         - c_x \dot{x}_{_{TMD_X/P_N}}
         - k_x x_{_{TMD_X/P_N}}
         + m_x a_{_{G_X/O_N}}
         + F_{ext_x}
         + F_{StopFrc_{X}} \\
         F_{Y_{_{TMD_X/O_N}}}
         + m_x a_{_{G_Y/O_N}}  \\
         F_{Z_{_{TMD_X/O_N}}}
         + m_x a_{_{G_Z/O_N}}
      \end{array}
   \right]

:math:`TMD_X` is fixed to frame :math:`N` in the :math:`y` and :math:`z`
directions so that

.. math::
   {r}_{_{TMD_X/P_N}} = \left[
      \begin{array}{c}
         x_{_{TMD_X/P_N}} \\
         0 \\
         0 
      \end{array}
   \right]

The other components of Eqn. :eq:`EOM` are:

.. math::
   \vec{\omega}_{_{N/O_N}} \times (\vec{\omega}_{_{N/O_N}} \times \vec{r}_{_{TMD_X/P_N}})
         = x_{_{TMD_X/P_N}} \left[
      \begin{array}{c}
         - (\dot{\phi}_{_{N/O_N}}^2 + \dot{\psi}_{_{N/O_N}}^2) \\
         \dot{\theta}_{_{N/O_N}}\dot{\phi}_{_{N/O_N}} \\
         \dot{\theta}_{_{N/O_N}}\dot{\psi}_{_{N/O_N}}
       \end{array}
   \right]

.. math::
   2\vec{\omega}_{_{N/O_N}} \times \dot{\vec{r}}_{_{TMD_X/P_N}}
         = \dot{x}_{_{TMD_X/P_N}} \left[
      \begin{array}{c} 0 \\
         2\dot{\psi}_{_{N/O_N}} \\
         -2\dot{\phi}_{_{N/O_N}}
      \end{array}
   \right]

.. math:: \vec{\alpha}_{_{N/O_N}} \times \vec{r}_{_{TMD_X/P_N}} = x_{_{TMD_X/P_N}} \left[ \begin{array}{c} 0 \\ \ddot{\psi}_{_{N/O_N}} \\ -\ddot{\phi}_{_{N/O_N}}\end{array} \right]

Therefore :math:`\ddot{x}_{_{TMD_X/P_N}}` is governed by the equations

.. math::
   \begin{aligned}
      \ddot{x}_{_{TMD_X/P_N}} =& (\dot{\phi}_{_{N/O_N}}^2 
         + \dot{\psi}_{_{N/O_N}}^2-\frac{k_x}{m_x}) x_{_{TMD_X/P_N}}
         - (\frac{c_x}{m_x}) \dot{x}_{_{TMD_X/P_N}} 
         -\ddot{x}_{_{P/O_N}}+a_{_{G_X/O_N}} \\ 
      &+ \frac{1}{m_x} ( F_{ext_X} + F_{StopFrc_{X}})
   \end{aligned}
   :label: EOM_Xx

The forces :math:`F_{Y_{_{TMD_X/O_N}}}` and :math:`F_{Z_{_{TMD_X/O_N}}}`
are solved noting
:math:`\ddot{y}_{_{TMD_X/P_N}} = \ddot{z}_{_{TMD_X/P_N}} = 0`:

.. math::
   F_{Y_{_{TMD_X/O_N}}} = m_x \left( - a_{_{G_Y/O_N}} +\ddot{y}_{_{P/O_N}} 
      + (\ddot{\psi}_{_{N/O_N}}
      + \dot{\theta}_{_{N/O_N}}\dot{\phi}_{_{N/O_N}} ) x_{_{TMD_X/P_N}}
      + 2\dot{\psi}_{_{N/O_N}} \dot{x}_{_{TMD_X/P_N}} \right)
   :label: EOM_Xy

.. math::
   F_{Z_{_{TMD_X/O_N}}} = m_x \left( - a_{_{G_Z/O_N}} +\ddot{z}_{_{P/O_N}}
      - (\ddot{\phi}_{_{N/O_N}}
      - \dot{\theta}_{_{N/O_N}}\dot{\psi}_{_{N/O_N}} ) x_{_{TMD_X/P_N}}
      - 2\dot{\phi}_{_{N/O_N}} \dot{x}_{_{TMD_X/P_N}} \right)
   :label: EOM_Xz
    
TMD_Y:
~~~~~~

The external forces :math:`\vec{F}_{_{TMD_Y/P_N}}` on :math:`TMD_Y` are
given by

.. math::
   \vec{F}_{_{TMD_Y/P_N}} =  \left[
      \begin{array}{c}
         F_{X_{_{TMD_Y/O_N}}} + m_y a_{_{G_X/O_N}}\\
         - c_y \dot{y}_{_{TMD_Y/P_N}} - k_y y_{_{TMD_Y/P_N}}
         + m_y a_{_{G_Y/O_N}} + F_{ext_y} + F_{StopFrc_{Y}} \\
         F_{Z_{_{TMD_Y/O_N}}}+ m_y a_{_{G_Z/O_N}}
      \end{array}
   \right]

:math:`TMD_Y` is fixed to frame :math:`N` in the :math:`x` and :math:`z`
directions so that

.. math::
   {r}_{_{TMDYX/P_N}} = \left[
      \begin{array}{c}
         0 \\
         y_{_{TMD_Y/P_N}} \\
         0
      \end{array}
   \right]

The other components of Eqn. :eq:`EOM` are:

.. math::
   \vec{\omega}_{_{N/O_N}} \times (\vec{\omega}_{_{N/O_N}}
         \times \vec{r}_{_{TMD_Y/P_N}})
      = y_{_{TMD_Y/P_N}}
      \left[
         \begin{array}{c}
            \dot{\theta}_{_{N/O_N}}\dot{\phi}_{_{N/O_N}} \\
            -(\dot{\theta}_{_{N/O_N}}^2 + \dot{\psi}_{_{N/O_N}}^2)  \\
            \dot{\phi}_{_{N/O_N}}\dot{\psi}_{_{N/O_N}}
         \end{array}
      \right]

.. math::
   2\vec{\omega}_{_{N/O_N}} \times \dot{\vec{r}}_{_{TMD_Y/P_N}}
      = \dot{y}_{_{TMD_Y/P_N}} \left[
         \begin{array}{c}
            - 2 \dot{\psi}_{_{N/O_N}} \\
            0 \\
            2 \dot{\theta}_{_{N/O_N}}
         \end{array}
      \right]

.. math::
   \vec{\alpha}_{_{N/O_N}} \times \vec{r}_{_{TMD_Y/P_N}}
      = y_{_{TMD_Y/P_N}} \left[
         \begin{array}{c}
            - \ddot{\psi}_{_{N/O_N}} \\
            0 \\
            \ddot{\theta}_{_{N/O_N}}
         \end{array}
      \right]

Therefore :math:`\ddot{y}_{_{TMD_Y/P_N}}` is governed by the equations

.. math::
   \begin{aligned}
      \ddot{y}_{_{TMD_Y/P_N}}
         = & (\dot{\theta}_{_{N/O_N}}^2
            + \dot{\psi}_{_{N/O_N}}^2-\frac{k_y}{m_y}) y_{_{TMD_Y/P_N}}
            - (\frac{c_y}{m_y}) \dot{y}_{_{TMD_Y/P_N}} 
            -\ddot{y}_{_{P/O_N}} + a_{_{G_Y/O_N}}\\ 
         &+ \frac{1}{m_y} (F_{ext_Y} + F_{StopFrc_{Y}})
   \end{aligned}
   :label: EOM_Yy

The forces :math:`F_{X_{_{TMD_Y/O_N}}}` and :math:`F_{Z_{_{TMD_Y/O_N}}}`
are solved noting
:math:`\ddot{x}_{_{TMD_Y/P_N}} = \ddot{z}_{_{TMD_Y/P_N}} = 0`:

.. math::
   F_{X_{_{TMD_Y/O_N}}} = m_y \left( - a_{_{G_X/O_N}} + \ddot{x}_{_{P/O_N}}
      - (\ddot{\psi}_{_{N/O_N}}
      - \dot{\theta}_{_{N/O_N}}\dot{\phi}_{_{N/O_N}}) y_{_{TMD_Y/P_N}}
      - 2\dot{\psi}_{_{N/O_N}} \dot{y}_{_{TMD_Y/P_N}} \right)
   :label: EOM_Yx

.. math::
   F_{Z_{_{TMD_Y/O_N}}} = m_y \left( - a_{_{G_Z/O_N}} + \ddot{z}_{_{P/O_N}}
      + (\ddot{\theta}_{_{N/O_N}}
      + \dot{\phi}_{_{N/O_N}}\dot{\psi}_{_{N/O_N}}) y_{_{TMD_Y/P_N}}
      + 2\dot{\theta}_{_{N/O_N}} \dot{y}_{_{TMD_Y/P_N}} \right)
   :label: EOM_Yz


TMD_Z :
~~~~~~~

The external forces :math:`\vec{F}_{_{TMD_Z/O_N}}` are given by

.. math::
   \vec{F}_{_{TMD_Z/O_N}} = \left[
      \begin{array}{c}
         F_{X_{_{TMD_Z/O_N}}} + m_z a_{_{G_X/O_N}} \\
         F_{Y_{_{TMD_Z/O_N}}} + m_z a_{_{G_Y/O_N}} \\
         - c_z \dot{z}_{_{TMD_Z/P_N}} - k_z z_{_{TMD_Z/P_N}}
         + m_z a_{_{G_Z/O_N}} + F_{ext_z} + F_{StopFrc_{Z}} + F_{Z_{PreLoad}}
      \end{array}
   \right]

where :math:`F_{Z_{PreLoad}}` is a spring pre-load to shift the neutral position
when gravity acts upon the mass for the :math:`TMD_Z`.
:math:`TMD_Z` is fixed to frame :math:`N` in the :math:`x` and :math:`y`
directions so that

.. math::
   {r}_{_{TMD_Z/P_N}} = \left[
      \begin{array}{c}
         0 \\
         0 \\
         z_{_{TMD_Z/P_N}}
      \end{array}
   \right]

The other components of Eqn. :eq:`EOM` are:

.. math::
   \vec{\omega}_{_{N/O_N}} \times (\vec{\omega}_{_{N/O_N}} \times \vec{r}_{_{TMD_Z/P_N}})
      = z_{_{TMD_Z/P_N}} \left[
         \begin{array}{c}
            \dot{\theta}_{_{N/O_N}}\dot{\psi}_{_{N/O_N}} \\
            \dot{\phi}_{_{N/O_N}}\dot{\psi}_{_{N/O_N}} \\
            -(\dot{\theta}_{_{N/O_N}}^2 + \dot{\phi}_{_{N/O_N}}^2)
         \end{array}
      \right]

.. math::
   2\vec{\omega}_{_{N/O_N}} \times \dot{\vec{r}}_{_{TMD_Z/P_N}}
      = \dot{z}_{_{TMD_Z/P_N}} \left[
         \begin{array}{c}
            2\dot{\phi}_{_{N/O_N}} \\
            -2\dot{\theta}_{_{N/O_N}} \\
            0
         \end{array}
      \right]

.. math::
   \vec{\alpha}_{_{N/O_N}} \times \vec{r}_{_{TMD_Z/P_N}}
      = z_{_{TMD_Z/P_N}} \left[
         \begin{array}{c}
            \ddot{\phi}_{_{N/O_N}} \\
            -\ddot{\theta}_{_{N/O_N}} \\
            0
         \end{array}
      \right]

Therefore :math:`\ddot{z}_{_{TMD_Z/P_N}}` is governed by the equations

.. math::
   \begin{aligned}
      \ddot{z}_{_{TMD_Z/P_N}}
         = & (\dot{\theta}_{_{N/O_N}}^2
            + \dot{\phi}_{_{N/O_N}}^2-\frac{k_z}{m_z}) z_{_{TMD_Z/P_N}}
            - (\frac{c_z}{m_z}) \dot{z}_{_{TMD_Z/P_N}} 
            -\ddot{z}_{_{P/O_N}} + a_{_{G_Z/O_N}}\\
         &+ \frac{1}{m_z} (F_{ext_Z} + F_{StopFrc_{Z}} + F_{Z_{PreLoad}})
   \end{aligned}
   :label: EOM_Zz



The forces :math:`F_{X_{_{TMD_Z/O_N}}}` and :math:`F_{Z_{_{TMD_Z/O_N}}}`
are solved noting
:math:`\ddot{x}_{_{TMD_Z/P_N}} = \ddot{y}_{_{TMD_Z/P_N}} = 0`:

.. math::
   F_{X_{_{TMD_Z/O_N}}} = m_z \left( - a_{_{G_X/O_N}} + \ddot{x}_{_{P/O_N}}
      + (\ddot{\phi}_{_{N/O_N}}
      + \dot{\theta}_{_{N/O_N}}\dot{\psi}_{_{N/O_N}}) z_{_{TMD_Z/P_N}}
      + 2\dot{\phi}_{_{N/O_N}} \dot{z}_{_{TMD_Z/P_N}} \right)
   :label: EOM_Zx

.. math::
   F_{Y_{_{TMD_Z/O_N}}} = m_z \left( - a_{_{G_Y/O_N}} + \ddot{y}_{_{P/O_N}}
      - (\ddot{\theta}_{_{N/O_N}}
      - \dot{\phi}_{_{N/O_N}}\dot{\psi}_{_{N/O_N}}) z_{_{TMD_Z/P_N}}
      - 2\dot{\theta}_{_{N/O_N}} \dot{z}_{_{TMD_Z/P_N}} \right)
   :label: EOM_Zy
    

State Equations
---------------

Inputs:
~~~~~~~

The inputs are the component linear acceleration and angular position,
velocity and acceleration:

.. math::
   \vec{u} = \left[
      \begin{array}{c}
         \ddot{\vec{r}}_{_{P/O_G}} \\
         \vec{R}_{_{N/G}} \\
         \vec{\omega}_{_{N/O_G}} \\
         \vec{\alpha}_{_{N/O_G}}
      \end{array}
   \right]
   \Rightarrow \left[
      \begin{array}{c}
         \ddot{\vec{r}}_{_{P/O_N}} \\
         \vec{\omega}_{_{N/O_N}} \\
         \vec{\alpha}_{_{N/O_N}}
      \end{array}
    \right]
    = \left[
      \begin{array}{c}
         \vec{R}_{_{N/G}} \ddot{\vec{r}}_{_{P/O_G}} \\
         \vec{R}_{_{N/G}} \vec{\omega}_{_{N/O_G}} \\
         \vec{R}_{_{N/G}} \vec{\alpha}_{_{N/O_G}
      }\end{array}
   \right]

States:
~~~~~~~

The states are the position and velocity of the TMDs along their
respective DOFs in the component reference frame:

.. math::
   \vec{R}_{_{TMD/P_N}} = \left[
      \begin{array}{c}
         x \\
         \dot{x} \\
         y \\
         \dot{y} \\
         z \\
         \dot{z}
      \end{array}
   \right]_{_{TMD/P_N}} 
   = \left[
      \begin{array}{c}
         {x}_{_{TMD_X/P_N}} \\
         \dot{x}_{_{TMD_X/P_N}} \\
         {y}_{_{TMD_Y/P_N}} \\
         \dot{y}_{_{TMD_Y/P_N}} \\
         {z}_{_{TMD_Z/P_N}} \\
         \dot{z}_{_{TMD_Z/P_N}}
      \end{array}
   \right]

The equations of motion can be re-written as a system of non-linear
first-order equations of the form

.. math::
   \dot{\vec{R}}_{_{TMD}} = A \vec{R}_{_{TMD}} + B

\ where

.. math::
   A(\vec{u}) = \left[
      \begin{array}{cccccc}
      0& 1 &0&0&0&0 \\
      (\dot{\phi}_{_{P/O_N}}^2 + \dot{\psi}_{_{P/O_N}}^2-\frac{k_x}{m_x}) & - (\frac{c_x}{m_x}) &0&0&0&0 \\
      0&0&0& 1 &0&0 \\
      0&0& (\dot{\theta}_{_{P/O_N}}^2 + \dot{\psi}_{_{P/O_N}}^2-\frac{k_y}{m_y}) & - (\frac{c_y}{m_y}) &0&0 \\
      0&0&0&0&0& 1 \\
      0&0&0&0& (\dot{\theta}_{_{P/O_N}}^2 + \dot{\phi}_{_{P/O_N}}^2-\frac{k_z}{m_z}) & - (\frac{c_z}{m_z}) \\ 
   \end{array} \right]

and

.. math::
   B(\vec{u}) = \left[
      \begin{array}{l}
         0 \\
         -\ddot{x}_{_{P/O_N}}+a_{_{G_X/O_N}} + \frac{1}{m_x} ( F_{ext_X} + F_{StopFrc_{X}}) \\
         0 \\
         -\ddot{y}_{_{P/O_N}}+a_{_{G_Y/O_N}} + \frac{1}{m_y} (F_{ext_Y}+ F_{StopFrc_{Y}}) \\
         0 \\
         -\ddot{z}_{_{P/O_N}}+a_{_{G_Z/O_N}} + \frac{1}{m_z} (F_{ext_Z}+ F_{StopFrc_{Z}} + F_{Z_{PreLoad}})
      \end{array}
   \right]
   :label: Bu

The inputs are coupled to the state variables, resulting in A and B as
:math:`f(\vec{u})`.

Outputs
-------

The output vector :math:`\vec{Y}` is

.. math::
   \vec{Y} = \left[
      \begin{array}{c}
         \vec{F}_{_{P_G}} \\
         \vec{M}_{_{P_G}}
      \end{array}
   \right]

The output includes reaction forces corresponding to
:math:`F_{Y_{_{TMD_X/O_N}}}`, :math:`F_{Z_{_{TMD_X/O_N}}}`,
:math:`F_{X_{_{TMD_Y/O_N}}}`, :math:`F_{Z_{_{TMD_Y/O_N}}}`,
:math:`F_{X_{_{TMD_Z/O_N}}}`, and :math:`F_{Y_{_{TMD_Z/O_N}}}` from Eqns.
:eq:`EOM_Xy`, :eq:`EOM_Xz`, :eq:`EOM_Yx`, :eq:`EOM_Yz`, :eq:`EOM_Zx`, and
:eq:`EOM_Zy`. The resulting forces :math:`\vec{F}_{_{P_G}}` and moments
:math:`\vec{M}_{_{P_G}}` acting on the component are

.. math::
   \begin{aligned}
      \vec{F}_{_{P_G}} = R^T_{_{N/G}} & \left[
         \begin{array}{l}
            k_x {x}_{_{TMD/P_N}} + c_x \dot{x}_{_{TMD/P_N}} - F_{StopFrc_{X}} - F_{ext_x} - F_{X_{_{TMD_Y/O_N}}} - F_{X_{_{TMD_Z/O_N}}} \\ 
            k_y {y}_{_{TMD/P_N}} + c_y \dot{y}_{_{TMD/P_N}} - F_{StopFrc_{Y}} - F_{ext_y} - F_{Y_{_{TMD_X/O_N}}} - F_{Y_{_{TMD_Z/O_N}}} \\ 
            k_z {z}_{_{TMD/P_N}} + c_z \dot{z}_{_{TMD/P_N}} - F_{StopFrc_{Z}} - F_{ext_z} - F_{Z_{_{TMD_X/O_N}}} - F_{Z_{_{TMD_Y/O_N}}} - F_{Z_{PreLoad}}
         \end{array}
      \right]
   \end{aligned}
   :label: OutputForces

and

.. math::
   \vec{M}_{_{P_G}} = R^T_{_{N/G}} \left[
      \begin{array}{c}
         M_{_X} \\
         M_{_Y} \\
         M_{_Z}
      \end{array}
   \right]_{_{N/N}} = R^T_{_{N/G}} \left[
      \begin{array}{c}
         -(F_{Z_{_{TMD_Y/O_N}}}) y_{_{TMD/P_N}} + (F_{Y_{_{TMD_Z/O_N}}} ) z_{_{TMD/P_N}} \\
          (F_{Z_{_{TMD_X/O_N}}}) x_{_{TMD/P_N}} - (F_{X_{_{TMD_Z/O_N}}} ) z_{_{TMD/P_N}} \\
         -(F_{Y_{_{TMD_X/O_N}}}) x_{_{TMD/P_N}} + ( F_{X_{_{TMD_Y/O_N}}}) y_{_{TMD/P_N}}
      \end{array}
   \right]

Stop Forces
~~~~~~~~~~~

The extra forces :math:`F_{StopFrc_{X}}`, :math:`F_{StopFrc_{Y}}`, and
:math:`F_{StopFrc_{Z}}` are added to output forces in the case that the
movement of TMD_X, TMD_Y, or TMD_Z exceeds the maximum track length for
the mass. Otherwise, they equal zero. The track length has limits on the
positive and negative ends in the TMD direction (X_PSP and X_NSP, Y_PSP
and Y_NSP, and Z_PSP and Z_NSP). If we define a general maximum and
minimum displacements as :math:`x_{max}` and :math:`x_{min}`,
respectively, the stop forces have the form

.. math::
   F_{StopFrc} = -\left\{
      \begin{array}{lr}
         \begin{aligned}
            k_S \Delta x  & \quad : ( x > x_{max} \wedge \dot{x}<=0) \vee ( x < x_{min} \wedge \dot{x}>=0)\\
            k_S \Delta x + c_S \dot{x} & \quad : ( x > x_{max} \wedge \dot{x}>0) \vee ( x < x_{min} \wedge \dot{x}<0)\\
            0 & \quad : \text{otherwise}
         \end{aligned}
      \end{array}
   \right.

where :math:`\Delta x` is the distance the mass has traveled beyond the
stop position and :math:`k_S` and :math:`c_S` are large stiffness and
damping constants.


.. _SrvD-StCz-PreLoad:

Pre-Load Forces
~~~~~~~~~~~~~~~

The extra force :math:`F_{Z_{PreLoad}}` is added to the output forces as a
method to shift the at rest position of the TMD_Z when gravity is acting on it.
This is particularly useful for substructure mounted StCs when very large masses
with soft spring constants are used. This appears in the term
:math:`\vec{F}_{_{TMD_Z/O_N}}` and in eq equations of motion given by :eq:`Bu`
and resulting forces in :eq:`OutputForces`.


Code Modifications
==================

The Structural Control (StC) function is a submodule linked into ServoDyn. In
addition to references in ServoDyn.f90 and ServoDyn.txt, new files that contain
the StC module are listed below.

New Files
---------

-  StrucCtrl.f90 : the structural control module

-  StrucCtrl.txt : registry file include files, inputs, states, parameters,
   and outputs shown in Tables `1 <#tbl2>`__ and `2 <#tbl1>`__

-  StrucCtrl_Types.f90: automatically generated

Variables
---------

.. container::
   :name: tbl2

   .. table:: Summary of field definitions in the StC registry. Note that state vector :math:`\vec{tmd_x}` corresponds to :math:`\vec{R}_{_{TMD/P_N}}`, and that the outputs :math:`\vec{F}_{_{P_G}}` and :math:`\vec{M}_{_{P_G}}` are contained    in the MeshType object (y.Mesh). :math:`X_{DSP}`, :math:`Y_{DSP}`, and :math:`Z_{DSP}` are initial displacements of the TMDs.

      +----------------------+------------------------------------------------------------------------------+
      + DataType             + Variable name                                                                +
      +======================+==============================================================================+
      | **InitInput**        |                                                                              |
      +----------------------+------------------------------------------------------------------------------+
      |                      | InputFile                                                                    |
      +----------------------+------------------------------------------------------------------------------+
      |                      | Gravity                                                                      |
      +----------------------+------------------------------------------------------------------------------+
      |                      | :math:`\vec{r}_{_{N/O_G}}`                                                   |
      +----------------------+------------------------------------------------------------------------------+
      | **Input u**          |                                                                              |
      +----------------------+------------------------------------------------------------------------------+
      |                      | :math:`\ddot{\vec{r}}_{_{P/O_G}}`                                            |
      +----------------------+------------------------------------------------------------------------------+
      |                      | :math:`\vec{R}_{_{N/O_G}}`                                                   |
      +----------------------+------------------------------------------------------------------------------+
      |                      | :math:`\vec{\omega}_{_{N/O_G}}`                                              |
      +----------------------+------------------------------------------------------------------------------+
      |                      | :math:`\vec{\alpha}_{_{N/O_G}}`                                              |
      +----------------------+------------------------------------------------------------------------------+
      | **Parameter p**      |                                                                              |
      +----------------------+------------------------------------------------------------------------------+
      |                      | :math:`m_x`                                                                  |
      +----------------------+------------------------------------------------------------------------------+
      |                      | :math:`c_x`                                                                  |
      +----------------------+------------------------------------------------------------------------------+
      |                      | :math:`k_x`                                                                  |
      +----------------------+------------------------------------------------------------------------------+
      |                      | :math:`m_y`                                                                  |
      +----------------------+------------------------------------------------------------------------------+
      |                      | :math:`c_y`                                                                  |
      +----------------------+------------------------------------------------------------------------------+
      |                      | :math:`k_y`                                                                  |
      +----------------------+------------------------------------------------------------------------------+
      |                      | :math:`m_z`                                                                  |
      +----------------------+------------------------------------------------------------------------------+
      |                      | :math:`c_z`                                                                  |
      +----------------------+------------------------------------------------------------------------------+
      |                      | :math:`k_z`                                                                  |
      +----------------------+------------------------------------------------------------------------------+
      |                      | :math:`K_S = \left[ k_{SX}\hspace{1em}k_{SY}\hspace{1em}k_{SZ}\right]`       |
      +----------------------+------------------------------------------------------------------------------+
      |                      | :math:`C_S = \left[c_{SX}\hspace{1em}c_{SY}\hspace{1em}c_{SZ}\right]`        |
      +----------------------+------------------------------------------------------------------------------+
      |                      | :math:`P_{SP}=\left[X_{PSP}\hspace{1em}Y_{PSP}\hspace{1em}Z_{PSP}\right]`    |
      +----------------------+------------------------------------------------------------------------------+
      |                      | :math:`P_{SP}=\left[X_{NSP}\hspace{1em}Y_{NSP}\hspace{1em}Z_{NSP}\right]`    |
      +----------------------+------------------------------------------------------------------------------+
      |                      | :math:`F{ext}`                                                               |
      +----------------------+------------------------------------------------------------------------------+
      |                      | :math:`Gravity`                                                              |
      +----------------------+------------------------------------------------------------------------------+
      |                      | TMDX_DOF                                                                     |
      +----------------------+------------------------------------------------------------------------------+
      |                      | TMDY_DOF                                                                     |
      +----------------------+------------------------------------------------------------------------------+
      |                      | TMDZ_DOF                                                                     |
      +----------------------+------------------------------------------------------------------------------+
      |                      | :math:`X_{DSP}`                                                              |
      +----------------------+------------------------------------------------------------------------------+
      |                      | :math:`Y_{DSP}`                                                              |
      +----------------------+------------------------------------------------------------------------------+
      |                      | :math:`Z_{DSP}`                                                              |
      +----------------------+------------------------------------------------------------------------------+
      | **State x**          |                                                                              |
      +----------------------+------------------------------------------------------------------------------+
      |                      | :math:`\vec{tmd_x}`                                                          |
      +----------------------+------------------------------------------------------------------------------+
      | **Output y**         |                                                                              |
      +----------------------+------------------------------------------------------------------------------+
      |                      | Mesh                                                                         |
      +----------------------+------------------------------------------------------------------------------+


The input, parameter, state and output definitions are summarized in
Table `1 <#tbl2>`__. The inputs from file are listed in Table
`2 <#tbl1>`__.

.. container::
   :name: tbl1

   .. table:: Data read in from TMDInputFile.

      +------------+------------+------------------------------------------------------+
      | Field Name | Field Type | Description                                          |
      +============+============+======================================================+
      | TMD_CMODE  | int        | Control Mode (1:passive, 2:active)                   |
      +------------+------------+------------------------------------------------------+
      | TMD_X_DOF  | logical    | DOF on or off                                        |
      +------------+------------+------------------------------------------------------+
      | TMD_Y_DOF  | logical    | DOF on or off                                        |
      +------------+------------+------------------------------------------------------+
      | TMD_Z_DOF  | logical    | DOF on or off                                        |
      +------------+------------+------------------------------------------------------+
      | TMD_X_DSP  | real       | TMD_X initial displacement                           |
      +------------+------------+------------------------------------------------------+
      | TMD_Y_DSP  | real       | TMD_Y initial displacement                           |
      +------------+------------+------------------------------------------------------+
      | TMD_Z_DSP  | real       | TMD_Z initial displacement                           |
      +------------+------------+------------------------------------------------------+
      | TMD_X_M    | real       | TMD mass                                             |
      +------------+------------+------------------------------------------------------+
      | TMD_X_K    | real       | TMD stiffness                                        |
      +------------+------------+------------------------------------------------------+
      | TMD_X_C    | real       | TMD damping                                          |
      +------------+------------+------------------------------------------------------+
      | TMD_Y_M    | real       | TMD mass                                             |
      +------------+------------+------------------------------------------------------+
      | TMD_Y_K    | real       | TMD stiffness                                        |
      +------------+------------+------------------------------------------------------+
      | TMD_Y_C    | real       | TMD damping                                          |
      +------------+------------+------------------------------------------------------+
      | TMD_Z_M    | real       | TMD mass                                             |
      +------------+------------+------------------------------------------------------+
      | TMD_Z_K    | real       | TMD stiffness                                        |
      +------------+------------+------------------------------------------------------+
      | TMD_Z_C    | real       | TMD damping                                          |
      +------------+------------+------------------------------------------------------+
      | TMD_X_PSP  | real       | positive stop position (maximum X mass displacement) |
      +------------+------------+------------------------------------------------------+
      | TMD_X_NSP  | real       | negative stop position (minimum X mass displacement) |
      +------------+------------+------------------------------------------------------+
      | TMD_X_K_SX | real       | stop spring stiffness                                |
      +------------+------------+------------------------------------------------------+
      | TMD_X_C_SX | real       | stop spring damping                                  |
      +------------+------------+------------------------------------------------------+
      | TMD_Y_PSP  | real       | positive stop position (maximum Y mass displacement) |
      +------------+------------+------------------------------------------------------+
      | TMD_Y_NSP  | real       | negative stop position (minimum Y mass displacement) |
      +------------+------------+------------------------------------------------------+
      | TMD_Y_K_S  | real       | stop spring stiffness                                |
      +------------+------------+------------------------------------------------------+
      | TMD_Y_C_S  | real       | stop spring damping                                  |
      +------------+------------+------------------------------------------------------+
      | TMD_Z_PSP  | real       | positive stop position (maximum Z mass displacement) |
      +------------+------------+------------------------------------------------------+
      | TMD_Z_NSP  | real       | negative stop position (minimum Z mass displacement) |
      +------------+------------+------------------------------------------------------+
      | TMD_Z_K_S  | real       | stop spring stiffness                                |
      +------------+------------+------------------------------------------------------+
      | TMD_Z_C_S  | real       | stop spring damping                                  |
      +------------+------------+------------------------------------------------------+
      | TMD_P_X    | real       | x origin of P in nacelle coordinate system           |
      +------------+------------+------------------------------------------------------+
      | TMD_P_Y    | real       | y origin of P in nacelle coordinate system           |
      +------------+------------+------------------------------------------------------+
      | TMD_P_Z    | real       | z origin of P in nacelle coordinate system           |
      +------------+------------+------------------------------------------------------+

Acknowledgements
================

The authors would like to thank Dr. Jason Jonkman for reviewing this
manual.

.. [1]
   Note that :math:`( R a ) \times ( Rb ) = R( a \times b )`.
