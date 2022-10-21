.. _ed_theory:


ElastoDyn Theory
================

Note this document is work in progress and is greatly incomplete. 
This documentation was started to document some code changes to the the tail furl and rotor furl part of ElastoDyn. 
Please refer to the different ressources provided in :numref:`ed_intro` for additional documents.



Notations
---------

**Points**

The following (partial) list of points are defined by ElastoDyn:

- ``Z``: the platform reference point
- ``O``: the tower-top/base plate point
- ``W``: the specified point on the tail-furl axis
- ``I``: the tail boom center of mass
- ``J``: the tail fin center of mass

**Bodies**

The following (partial) list of bodies are defined by ElastoDyn:

- ``E``: the earth/inertial frame
- ``X``: the platform body
- ``N``: the nacelle body
- ``A``: the tail-furl body





Kinematics
----------


ElastoDyn computes the position, velocity and accelerations of key points of the structure, starting from the platform reference point `Z` and going up in the structure.

The different position vectors are available in the data stucture ``RtHSdat``.
For instance, the global position of point J is given by:

.. math::  :label: TFPointJPos

   \boldsymbol{r}_J =  \boldsymbol{r}_Z +  \boldsymbol{r}_{ZO} +  \boldsymbol{r}_{OW} +  \boldsymbol{r}_{WJ}

The translational displacement vector (how much a point has moved compared to its reference position)  is calculated as follows: :math:`\boldsymbol{r}_J-\boldsymbol{r}_{J,\text{ref}}`.


The coordinate systems of ElastoDyn are stored in the variable ``CoordSys``.
The orientation matrix of a given coordinate system can be formed using the unit vectors (assumed to be column vectors) of a given coordinate system expressed in the inertial frame.
For instance for the tailfin coordinate system:

.. math::  :label: TFPointJOrientation

   \boldsymbol{R}_{Ai} = \begin{bmatrix}
      \left.\boldsymbol{\hat{x}_\text{tf}^t}\right|_i \\
      \left.\boldsymbol{\hat{y}_\text{tf}^t}\right|_i \\
      \left.\boldsymbol{\hat{z}_\text{tf}^t}\right|_i \\
   \end{bmatrix}


Angular velocities are stored in variables ``RtHSdat%AngVelE*`` with respect to the initial frame ("Earth", `E`). 
For instance, the angular velocity of the tail-furl body (body `A`) is:

.. math::  :label: TFPointJAngVel

   \boldsymbol{\omega}_{A/E} =  \boldsymbol{\omega}_{X/E} + \boldsymbol{\omega}_{N/X} + \boldsymbol{\omega}_{A/N}

where :math:`\boldsymbol{\omega}_{N/X}=\boldsymbol{\omega}_{B/X}+\boldsymbol{\omega}_{N/B}`



Linear (translational) velocities of the different points are found in the variables ``RtHSdat%LinVelE*``, and are computed based on Kane's partial velocities (which are Jacobians of the velocity with respect to the time derivatives of the degrees of freedom).
For instance, the linear velocity of point J is computed as:

.. math:: :label: TFPointJVel

    \boldsymbol{v}_J = \sum_{j} \frac{\partial v_J}{\partial \dot{q}_j} \dot{q}_j

where the Jacobians :math:`\frac{\partial v_J}{\partial \dot{q}_j}` are stored in ``RtHSdat%PLinVelEJ(:,0)`` 




Translational accelerations are computed as the sum of contribution from the first and second time derivatives of the degrees of freedom.
For instance, the acceleration of point `J` is computed as:

.. math:: :label: TFPointJAng

    \boldsymbol{\tilde{a}}_J &= \sum_{j\in PA} \frac{\partial a_J}{\partial \dot{q}_j} \dot{q}_j

    \boldsymbol{a}_J &= \boldsymbol{\tilde{a}}_J + \sum_{j\in PA} \frac{\partial v_J}{\partial \dot{q}_j} \ddot{q}_j


where  :math:`\frac{\partial a_J}{\partial \dot{q}_j}` are stored in ``RtHSdat%PLinVelEJ(:,1)``

Angular accelerations requires similar computations currently not documented.



.. _ed_rtfrl_theory:

Rotor and tail furl
-------------------

The user can select linear spring and damper models, together with 
up- and down-stop springs, and up- and down-stop dampers. 

The torque applied from the linear spring and damper is:

.. math::  :label: TFLinTorque

   Q_\text{lin} = - k \theta  - d \dot{\theta}

where :math:`\theta` is the degree of freedom (rotor or tail furl), 
:math:`k` is the linear spring constant (``RFrlSpr`` or ``TFrlSpr``)
:math:`d` is the linear damping constant (``RFrlDmp`` or ``TFrlDmp``).

The up-/down- stop spring torque is defined as:

.. math::  :label: TFStopTorqueSpring

   Q_\text{stop, spr} = \begin{cases} 
      - k_{US} (\theta-\theta_{k_{US}}),&\text{if } \theta>\theta_{k_{US}}  \\ 
      - k_{DS} (\theta-\theta_{k_{DS}}),&\text{if } \theta<\theta_{k_{DS}}  \\ 
        0 ,&\text{otherwise}
        \end{cases}

where 
:math:`k_{US}` is the up-stop spring constant (``RFrlUSSpr`` or ``TFrlUSSpr``),
:math:`\theta_{k_{US}}` is the up-stop spring angle (``RFrlUSSP`` or ``TFrlUSSP``),
and similar notations are used for the down-stop spring.

The up-/down- stop damping torque is defined as:

.. math::  :label: TFStopTorqueDamp

   Q_\text{stop, dmp} = \begin{cases} 
      - d_{US} \dot{\theta},&\text{if } \theta>\theta_{d_{US}}  \\ 
      - d_{DS} \dot{\theta},&\text{if } \theta<\theta_{d_{DS}}  \\ 
        0 ,&\text{otherwise}
        \end{cases}

where similar nnotations are used.
The total moment on the given degree of freedom is:

.. math::  :label: TFTotTorque

   Q = Q_\text{lin} + Q_\text{stop,spr} + Q_\text{stop,dmp}
   
     


   




.. _ed_dev_notes:


Developer notes
===============


**Internal coordinate systems**

The different coordinate systems of ElastoDyn are stored in the variable ``CoordSys``.
The coordinate systems used internally by ElastoDyn are using a different convention than the OpenFAST input/output coordinate system. 

For instance, for the coordinate system of the nacelle, with unit axes noted :math:`x_n,y_n,z_n` in OpenFAST, and :math:`d_1,d_2,d_3` in ElastoDyn, the following conversions apply:
:math:`d_1 = x_n`,  
:math:`d_2 =z_n` and 
:math:`d_3 =-y_n`.

The following (partial) list of coordinate systems are defined internally by ElastoDyn:

-  `z` : inertial coordinate system 
-  `a` : tower base coordinate system 
-  `t` : tower-node coordinate system (one per node)
-  `d` : nacelle coordinate system 
-  `c` : shaft-tilted coordinate system 
-  `rf` : rotor furl coordinate system 
-  `tf` : tail furl coordinate system 
-  `g` : hub coordinate system 
