.. _beamdyn-theory:

BeamDyn Theory
==============

This section focuses on the theory behind the BeamDyn module. The
theoretical foundation, numerical tools, and some special handling in
the implementation will be introduced. References will be provided in
each section detailing the theories and numerical tools.

In this chapter, matrix notation is used to denote vectorial or
vectorial-like quantities. For example, an underline denotes a vector
:math:`\underline{u}`, an over bar denotes unit vector :math:`\bar{n}`,
and a double underline denotes a tensor
:math:`\underline{\underline{\Delta}}`. Note that sometimes the
underlines only denote the dimension of the corresponding matrix.

Coordinate Systems
------------------

:numref:`blade-geometry` (in :numref:`bd-input-files`) and
:numref:`bd-frame` show the coordinate system used in BeamDyn.

.. _bd-frame:

.. figure:: figs/bd_frame.pdf
   :width: 100%
   :align: center

   Global, blade reference, and internal coordinate systems in BeamDyn. Illustration by Al Hicks, NREL.


Global Coordinate System
~~~~~~~~~~~~~~~~~~~~~~~~

The global coordinate system is denoted as ``X``, ``Y``, and ``Z``
in :numref:`bd-frame`. This is an inertial frame and in FAST its
origin is usually placed at the bottom of the tower as shown.

BD Coordinate System
~~~~~~~~~~~~~~~~~~~~

The BD coordinate system is denoted as :math:`x_1`, :math:`x_2`, and
:math:`x_3` respectively in :numref:`bd-frame`. This is an inertial
frame used internally in BeamDyn (i.e., doesn’t rotate with the rotor)
and its origin is placed at the initial position of the blade root
point.

Blade Reference Coordinate System
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The blade reference coordinate system is denoted as
:math:`X_{rt}`, :math:`Y_{rt}`, and
:math:`Z_{rt}` in :numref:`bd-frame` at initialization
(:math:`t = 0`). The blade reference coordinate system is a floating
frame that attaches at the blade root and is rotating with the blade.
Its origin is at the blade root and the directions of axes following the
IEC standard, i.e., :math:`Z_r` is pointing along the blade axis
from root to tip; :math:`Y_r` pointing nominally towards the
trailing edge of the blade and parallel with the chord line at the
zero-twist blade station; and :math:`X_r` is orthogonal with the
:math:`Y_r` and :math:`Z_r` axes, such that they form a
right-handed coordinate system (pointing nominally downwind). We note
that the initial blade reference coordinate system, denoted by subscript
:math:`r0`, coincides with the BD coordinate system, which is used
internally in BeamDyn and introduced in the previous section. The axis
convention relations between the initial blade reference coordinate
system and the BD coordinate system can be found in :numref:`IECBD`.

.. _IECBD:

.. table:: Transformation between blade coordinate system and BD coordinate system.

   +---------------+------------------+------------------+------------------+
   | Blade Frame   | :math:`X_{r0}`   | :math:`Y_{r0}`   | :math:`Z_{r0}`   |
   +---------------+------------------+------------------+------------------+
   | BD Frame      | :math:`x_2`      | :math:`x_3`      | :math:`x_1`      |
   +---------------+------------------+------------------+------------------+

Local blade coordinate system
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The local blade coordinate system is used for some input and output
quantities, for example, the cross-sectional mass and stiffness matrices
and the the sectional force and moment resultants. This coordinate
system is different from the blade reference coordinate system in that
its :math:`Z_l` axis is always tangent to the blade axis as the blade
deflects. Note that a subscript :math:`l` denotes the local blade
coordinate system.

Geometrically Exact Beam Theory
-------------------------------

The theoretical foundation of BeamDyn is the geometrically exact beam
theory. This theory features the capability of beams that are initially
curved and twisted and subjected to large displacement and rotations.
Along with a proper two-dimensional (2D) cross-sectional analysis, the
coupling effects between all six DOFs, including extension, bending,
shear, and torsion, can be captured by GEBT as well . The term,
“geometrically exact” refer to the fact that there is no approximation
made on the geometries, including both initial and deformed geometries,
in formulating the equations :cite:`HodgesBeamBook`.

The governing equations of motion for geometrically exact beam theory
can be written as :cite:`Bauchau:2010`

.. math::
   	:label: GovernGEBT-1-2

   	\dot{\underline{h}} - \underline{F}^\prime &= \underline{f} \\
   	\dot{\underline{g}} + \dot{\tilde{u}} \underline{h} - \underline{M}^\prime + (\tilde{x}_0^\prime + \tilde{u}^\prime)^T \underline{F} &= \underline{m}

where :math:`{\underline{h}}` and :math:`{\underline{g}}` are the
linear and angular momenta resolved in the inertial coordinate system,
respectively; :math:`{\underline{F}}` and :math:`{\underline{M}}` are
the beam’s sectional force and moment resultants, respectively;
:math:`{\underline{u}}` is the one-dimensional (1D) displacement of a
point on the reference line; :math:`{\underline{x}}_0` is the position
vector of a point along the beam’s reference line; and
:math:`{\underline{f}}` and :math:`{\underline{m}}` are the distributed
force and moment applied to the beam structure. The notation
:math:`(\bullet)^\prime` indicates a derivative with respect to beam
axis :math:`x_1` and :math:`\dot{(\bullet)}` indicates a derivative with
respect to time. The tilde operator :math:`({\widetilde{\bullet}})`
defines a skew-symmetric tensor corresponding to the given vector. In
the literature, it is also termed as “cross-product matrix”. For
example,

.. math::

   {\widetilde{n}} =
   		\begin{bmatrix}
		0 & -n_3 & n_2 \\
		n_3 & 0 & -n_1 \\
		-n_2 & n_1 & 0\\
		\end{bmatrix}

The constitutive equations relate the velocities to the momenta and the
1D strain measures to the sectional resultants as

.. math::
   	:label: ConstitutiveMass-Stiff

   	\begin{Bmatrix}
   	\underline{h} \\
   	\underline{g}
   	\end{Bmatrix}
   	= \underline{\underline{\mathcal{M}}} \begin{Bmatrix}
   	\dot{\underline{u}} \\
   	\underline{\omega}
   	\end{Bmatrix} \\

   	\begin{Bmatrix}
   	\underline{F} \\
   	\underline{M}
   	\end{Bmatrix}
   	= \underline{\underline{\mathcal{C}}} \begin{Bmatrix}
   	\underline{\epsilon} \\
   	\underline{\kappa}
   	\end{Bmatrix}

where :math:`\underline{\underline{\mathcal{M}}}` and
:math:`\underline{\underline{\mathcal{C}}}` are the :math:`6 \times 6`
sectional mass and stiffness matrices, respectively (note that they are
not really tensors); :math:`\underline{\epsilon}` and
:math:`\underline{\kappa}` are the 1D strains and curvatures,
respectively; and, :math:`\underline{\omega}` is the angular velocity
vector that is defined by the rotation tensor
:math:`\underline{\underline{R}}` as :math:`\underline{\omega} =
axial(\dot{\underline{\underline{R}}}~\underline{\underline{R}}^T)`. The
axial vector :math:`{\underline{a}}` associated with a second-order
tensor :math:`{\underline{\underline{A}}}` is denoted
:math:`{\underline{a}}=axial({\underline{\underline{A}}})` and its
components are defined as

.. math::
   :label: axial

   {\underline{a}} = axial({\underline{\underline{A}}})=\begin{Bmatrix}
   a_1 \\
   a_2 \\
   a_3
   \end{Bmatrix}
   =\frac{1}{2}
   \begin{Bmatrix}
   A_{32}-A_{23} \\
   A_{13}-A_{31} \\
   A_{21}-A_{12}
   \end{Bmatrix}

The 1D strain measures are defined as

.. math::
   :label: 1DStrain

   \begin{Bmatrix}
      {\underline{\epsilon}} \\
      {\underline{\kappa}}
   \end{Bmatrix}
   =
   \begin{Bmatrix}
           {\underline{x}}^\prime_0 + {\underline{u}}^\prime - ({\underline{\underline{R}}} ~{\underline{\underline{R}}}_0) \bar{\imath}_1 \\
           {\underline{k}}
   \end{Bmatrix}

where :math:`{\underline{k}} = axial [({\underline{\underline{R R_0}}})^\prime ({\underline{\underline{R R_0}}})^T]` is the sectional curvature vector resolved in the inertial basis; :math:`{\underline{\underline{R}}}_0` is the initial rotation tensor; and :math:`\bar{\imath}_1` is the unit vector along :math:`x_1` direction in the inertial basis. These three sets of equations, including equations of motion Eq. :eq:`GovernGEBT-1-2`, constitutive equations Eq. :eq:`ConstitutiveMass-Stiff`, and kinematical equations Eq. :eq:`1DStrain`, provide a full mathematical description of the beam elasticity problems.

.. _num-imp:

Numerical Implementation with Legendre Spectral Finite Elements
---------------------------------------------------------------

For a displacement-based finite element implementation, there are six
degree-of-freedoms at each node: three displacement components and three
rotation components. Here we use :math:`{\underline{q}}` to denote the
elemental displacement array as :math:`\underline{q}=\left[
\underline{u}^T~~\underline{c}^T\right]` where :math:`{\underline{u}}`
is the displacement and :math:`{\underline{c}}` is the
rotation-parameter vector. The acceleration array can thus be defined as
:math:`\underline{a}=\left[ \ddot{\underline{u}}^T~~ \dot{\underline{\omega}}^T \right]`.
For nonlinear finite-element analysis, the discretized and incremental
forms of displacement, velocity, and acceleration are written as

.. math::
     :label: Discretized

     \underline{q} (x_1) &= \underline{\underline{N}} ~\hat{\underline{q}}~~~~\Delta \underline{q}^T = \left[ \Delta \underline{u}^T~~\Delta \underline{c}^T \right] \\
     \underline{v}(x_1) &= \underline{\underline{N}}~\hat{\underline{v}}~~~~\Delta \underline{v}^T = \left[\Delta \underline{\dot{u}}^T~~\Delta \underline{\omega}^T \right] \\
     \underline{a}(x_1) &= \underline{\underline{N}}~ \hat{\underline{a}}~~~~\Delta \underline{a}^T = \left[ \Delta \ddot{\underline{u}}^T~~\Delta \dot{\underline{\omega}}^T \right]

where :math:`{\underline{\underline{N}}}` is the shape function matrix
and :math:`(\hat{\cdot})` denotes a column matrix of nodal values.

The displacement fields in an element are approximated as

.. math::
       :label: InterpolateDisp

       {\underline{u}}(\xi) &=  h^k(\xi) {\underline{\hat{u}}}^k \\
       {\underline{u}}^\prime(\xi) &=  h^{k\prime}(\xi) {\underline{\hat{u}}}^k

where :math:`h^k(\xi)`, the component of shape function matrix
:math:`{\underline{\underline{N}}}`, is the :math:`p^{th}`-order
polynomial Lagrangian-interpolant shape function of node :math:`k`,
:math:`k=\{1,2,...,p+1\}`, :math:`{\underline{\hat{u}}}^k` is the
:math:`k^{th}` nodal value, and :math:`\xi \in \left[-1,1\right]` is the
element natural coordinate. However, as discussed in
:cite:`Bauchau-etal:2008`, the 3D rotation field cannot
simply be interpolated as the displacement field in the form of

.. math::
       :label: InterpolateRot

       {\underline{c}}(\xi) &= h^k(\xi) {\underline{\hat{c}}}^k \\
       {\underline{c}}^\prime(\xi) &= h^{k \prime}(\xi) {\underline{\hat{c}}}^k

where :math:`{\underline{c}}` is the rotation field in an element and
:math:`{\underline{\hat{c}}}^k` is the nodal value at the :math:`k^{th}`
node, for three reasons:

1) rotations do not form a linear space so that they must be “composed” rather than added;
2) a rescaling operation is needed to eliminate the singularity existing in the vectorial rotation parameters;
3) the rotation field lacks objectivity, which, as defined by :cite:`Crisfield1999`, refers to the invariance of strain measures computed through interpolation to the addition of a rigid-bodymotion.

Therefore, we adopt the more robust interpolation approach
proposed by :cite:`Crisfield1999` to deal with the finite
rotations. Our approach is described as follows

Step 1:
    Compute the nodal relative rotations,
    :math:`{\underline{\hat{r}}}^k`, by removing the reference rotation,
    :math:`{\underline{\hat{c}}}^1`, from the finite rotation at each
    node,
    :math:`{\underline{\hat{r}}}^k = ({\underline{\hat{c}}}^{1-}) \oplus
    {\underline{\hat{c}}}^k`. It is noted that the minus sign on
    :math:`{\underline{\hat{c}}}^1` denotes that the relative rotation
    is calculated by removing the reference rotation from each node. The
    composition in that equation is an equivalent of
    :math:`{\underline{\underline{R}}}({\underline{\hat{r}}}^k) = {\underline{\underline{R}}}^T({\underline{\hat{c}}}^1)~{\underline{\underline{R}}}({\underline{{\underline{c}}}}^k).`

Step 2:
    Interpolate the relative-rotation field:
    :math:`{\underline{r}}(\xi) = h^k(\xi) {\underline{\hat{r}}}^k` and
    :math:`{\underline{r}}^\prime(\xi) = h^{k \prime}(\xi) {\underline{\hat{r}}}^k`.
    Find the curvature field
    :math:`{\underline{\kappa}}(\xi) = {\underline{\underline{R}}}({\underline{\hat{c}}}^1) {\underline{\underline{H}}}({\underline{r}}) {\underline{r}}^\prime`,
    where :math:`{\underline{\underline{H}}}` is the tangent tensor that
    relates the curvature vector :math:`{\underline{k}}` and rotation
    vector :math:`{\underline{c}}` as

    .. math::
       :label: Tensor

           {\underline{k}} = {\underline{\underline{H}}}~ {\underline{c}}^\prime

Step 3:
    Restore the rigid-body rotation removed in Step 1:
    :math:`{\underline{c}}(\xi) = {\underline{\hat{c}}}^1 \oplus {\underline{r}}(\xi)`.

Note that the relative-rotation field can be computed with respect to
any of the nodes of the element; we choose node 1 as the reference node
for convenience. In the LSFE approach, shape functions (i.e., those
composing :math:`{\underline{\underline{N}}}`) are :math:`p^{th}`-order
Lagrangian interpolants, where nodes are located at the :math:`p+1`
Gauss-Lobatto-Legendre (GLL) points in the :math:`[-1,1]` element
natural-coordinate domain. :numref:`N4_lsfe` shows representative
LSFE basis functions for fourth- and eighth-order elements. Note that
nodes are clustered near element endpoints. More details on the LSFE and
its applications can be found in
References :cite:`Patera:1984,Ronquist:1987,Sprague:2003,Sprague:2004`.


.. _N4_lsfe:

.. figure:: figs/n4.pdf
   :width: 47%
   :align: center

   Representative :math:`p+1` Lagrangian-interpolant shape functions in the element natural coordinates for a fourth-order LSFEs, where nodes are located at the Gauss-Lobatto-Legendre points.

.. _N8_lsfe:

.. figure:: figs/n8.pdf
   :width: 47%
   :align: center

   Representative :math:`p+1` Lagrangian-interpolant shape functions in the element natural coordinates for a eighth-order LSFEs, where nodes are located at the Gauss-Lobatto-Legendre points.



Wiener-Milenković Rotation Parameter
------------------------------------

In BeamDyn, the 3D rotations are represented as Wiener-Milenković
parameters defined in the following equation:

.. math::
   :label: WMParameter

   {\underline{c}} = 4 \tan\left(\frac{\phi}{4} \right) \bar{n}

where :math:`\phi` is the rotation angle and :math:`\bar{n}` is the
unit vector of the rotation axis. It can be observed that the valid
range for this parameter is :math:`|\phi| < 2 \pi`. The singularities
existing at integer multiples of :math:`\pm 2 \pi` can be removed by a
rescaling operation at :math:`\pi` as:

.. math::
   :label: RescaledWM

   {\underline{r}} = \begin{cases}
   4(q_0{\underline{p}} + p_0 {\underline{q}} + \tilde{p} {\underline{q}} ) / (\Delta_1 + \Delta_2), & \text{if } \Delta_2 \geq 0 \\
   -4(q_0{\underline{p}} + p_0 {\underline{q}} + \tilde{p} {\underline{q}} ) / (\Delta_1 - \Delta_2), & \text{if } \Delta_2 < 0
   \end{cases}

where :math:`{\underline{p}}`, :math:`{\underline{q}}`, and
:math:`{\underline{r}}` are the vectorial parameterization of three
finite rotations such that
:math:`{\underline{\underline{R}}}({\underline{r}}) = {\underline{\underline{R}}}({\underline{p}}) {\underline{\underline{R}}}({\underline{q}})`,
:math:`p_0 = 2 - {\underline{p}}^T {\underline{p}}/8`,
:math:`q_0 = 2 - {\underline{q}}^T {\underline{q}}/8`,
:math:`\Delta_1 = (4-p_0)(4-q_0)`, and
:math:`\Delta_2 = p_0 q_0 - {\underline{p}}^T {\underline{q}}`. It is
noted that the rescaling operation could cause a discontinuity of the
interpolated rotation field; therefore a more robust interpolation
algorithm has been introduced in Section :ref:`num-imp` where the
rescaling-independent relative-rotation field is interpolated.

The rotation tensor expressed in terms of Wiener-Milenković parameters is

.. math::
      :label: eqn:RotTensorWM

      {\underline{\underline{R}}} ({\underline{c}}) = \frac{1}{(4-c_0)^2}
      \begin{bmatrix}
      c_0^2 + c_1^2 - c_2^2 - c_3^2 & 2(c_1 c_2 - c_0 c_3) & 2(c_1 c_3 + c_0 c_2) \\
      2(c_1 c_2 + c_0 c_3) & c_0^2 - c_1^2 + c_2^2 - c_3^2 & 2(c_2 c_3 - c_0 c_1) \\
      2(c_1 c_3 - c_0 c_2)  & 2(c_2 c_3 + c_0 c_1) & c_0^2 - c_1^2 - c_2^2 + c_3^2 \\
      \end{bmatrix}

where :math:`{\underline{c}} = \left[ c_1~~c_2~~c_3\right]^T` is the
Wiener-Milenković parameter and
:math:`c_0 = 2 - \frac{1}{8}{\underline{c}}^T {\underline{c}}`. The
relation between rotation tensor and direction cosine matrix (DCM) is

.. math::
   :label: RT2DCM

   {\underline{\underline{R}}} = ({\underline{\underline{DCM}}})^T

Interested users are referred to :cite:`Bauchau-etal:2008`
and :cite:`Wang:GEBT2013` for more details on the rotation
parameter and its implementation with GEBT.

Linearization Process
---------------------

The nonlinear governing equations introduced in the previous section are
solved by Newton-Raphson method, where a linearization process is
needed. The linearization of each term in the governing equations are
presented in this section.

According to :cite:`Bauchau:2010`, the linearized governing
equations in Eq.  :eq:`GovernGEBT-1-2` are in the form of

.. math::
   :label: LinearizedEqn

   \hat{\underline{\underline{M}}} \Delta \hat{\underline{a}} +\hat{\underline{\underline{G}}} \Delta \hat{\underline{v}}+ \hat{\underline{\underline{K}}} \Delta \hat{\underline{q}} = \hat{\underline{F}}^{ext} - \hat{\underline{F}}

where the :math:`\hat{{\underline{\underline{M}}}}`,
:math:`\hat{{\underline{\underline{G}}}}`, and
:math:`\hat{{\underline{\underline{K}}}}` are the elemental mass,
gyroscopic, and stiffness matrices, respectively;
:math:`\hat{{\underline{F}}}` and :math:`\hat{{\underline{F}}}^{ext}`
are the elemental forces and externally applied loads, respectively.
They are defined for an element of length :math:`l` along :math:`x_1` as
follows

.. math::
   	:label: hatMGKFFext

   	\hat{{\underline{\underline{M}}}}&= \int_0^l \underline{\underline{N}}^T \mathcal{\underline{\underline{M}}} ~\underline{\underline{N}} dx_1 \\
   	\hat{{\underline{\underline{G}}}} &= \int_0^l {\underline{\underline{N}}}^T {\underline{\underline{\mathcal{G}}}}^I~{\underline{\underline{N}}} dx_1\\
   	\hat{{\underline{\underline{K}}}}&=\int_0^l \left[ {\underline{\underline{N}}}^T ({\underline{\underline{\mathcal{K}}}}^I + \mathcal{{\underline{\underline{Q}}}})~ {\underline{\underline{N}}} + {\underline{\underline{N}}}^T \mathcal{{\underline{\underline{P}}}}~ {\underline{\underline{N}}}^\prime + {\underline{\underline{N}}}^{\prime T} \mathcal{{\underline{\underline{C}}}}~ {\underline{\underline{N}}}^\prime + {\underline{\underline{N}}}^{\prime T} \mathcal{{\underline{\underline{O}}}}~ {\underline{\underline{N}}} \right] d x_1 \\
   	\hat{{\underline{F}}} &= \int_0^l ({\underline{\underline{N}}}^T {\underline{\mathcal{F}}}^I + {\underline{\underline{N}}}^T \mathcal{{\underline{F}}}^D + {\underline{\underline{N}}}^{\prime T} \mathcal{{\underline{F}}}^C)dx_1 \\
   	\hat{{\underline{F}}}^{ext}& = \int_0^l {\underline{\underline{N}}}^T \mathcal{{\underline{F}}}^{ext} dx_1

where :math:`\mathcal{{\underline{F}}}^{ext}` is the applied load
vector. The new matrix notations in Eqs. :eq:`hatMGKFFext` to are briefly introduced
here. :math:`\mathcal{{\underline{F}}}^C` and
:math:`\mathcal{{\underline{F}}}^D` are elastic forces obtained from
Eq. :eq:`GovernGEBT-1-2` as

.. math::
   	:label: FCD

   	\mathcal{{\underline{F}}}^C &= \begin{Bmatrix}
            {\underline{F}} \\
   	{\underline{M}}
   	\end{Bmatrix} = {\underline{\underline{\mathcal{C}}}} \begin{Bmatrix}
   	{\underline{\epsilon}} \\
   	{\underline{\kappa}}
   	\end{Bmatrix} \\
   	\mathcal{{\underline{F}}}^D & = \begin{bmatrix}
   	\underline{\underline{0}} & \underline{\underline{0}}\\
   	(\tilde{x}_0^\prime+\tilde{u}^\prime)^T & \underline{\underline{0}}
   	\end{bmatrix}
   	\mathcal{{\underline{F}}}^C \equiv {\underline{\underline{\Upsilon}}}~ \mathcal{{\underline{F}}}^C

where :math:`\underline{\underline{0}}` denotes a :math:`3 \times 3`
null matrix. The :math:`{\underline{\underline{\mathcal{G}}}}^I`,
:math:`{\underline{\underline{\mathcal{K}}}}^I`,
:math:`\mathcal{{\underline{\underline{O}}}}`,
:math:`\mathcal{{\underline{\underline{P}}}}`,
:math:`\mathcal{{\underline{\underline{Q}}}}`, and
:math:`{\underline{\mathcal{F}}}^I` in Eqs. :eq:`hatMGKFFext`  are defined as

.. math::
      :label: mathcalGKOPFI

      {\underline{\underline{\mathcal{G}}}}^I &= \begin{bmatrix}
      {\underline{\underline{0}}} & (\widetilde{\tilde{\omega} m {\underline{\eta}}})^T+\tilde{\omega} m \tilde{\eta}^T  \\
      {\underline{\underline{0}}} & \tilde{\omega} {\underline{\underline{\varrho}}}-\widetilde{{\underline{\underline{\varrho}}} {\underline{\omega}}}
      \end{bmatrix} \\
      {\underline{\underline{\mathcal{K}}}}^I &= \begin{bmatrix}
      {\underline{\underline{0}}} & \dot{\tilde{\omega}}m\tilde{\eta}^T + \tilde{\omega} \tilde{\omega}m\tilde{\eta}^T  \\
      {\underline{\underline{0}}} & \ddot{\tilde{u}}m\tilde{\eta} + {\underline{\underline{\varrho}}} \dot{\tilde{\omega}}-\widetilde{{\underline{\underline{\varrho}}} {\underline{\dot{\omega}}}}+\tilde{\omega} {\underline{\underline{\varrho}}} \tilde{\omega} - \tilde{\omega}  \widetilde{{\underline{\underline{\varrho}}} {\underline{\omega}}}
      \end{bmatrix}\\
      \mathcal{{\underline{\underline{O}}}} &= \begin{bmatrix}
      {\underline{\underline{0}}} & {\underline{\underline{C}}}_{11} \tilde{E_1} - \tilde{F} \\
      {\underline{\underline{0}}}& {\underline{\underline{C}}}_{21} \tilde{E_1} - \tilde{M}
      \end{bmatrix} \\
      \mathcal{{\underline{\underline{P}}}} &= \begin{bmatrix}
      {\underline{\underline{0}}} & {\underline{\underline{0}}} \\
      \tilde{F} +  ({\underline{\underline{C}}}_{11} \tilde{E_1})^T & ({\underline{\underline{C}}}_{21} \tilde{E_1})^T
      \end{bmatrix}  \\
      \mathcal{{\underline{\underline{Q}}}} &= {\underline{\underline{\Upsilon}}}~ \mathcal{{\underline{\underline{O}}}} \\
      {\underline{\mathcal{F}}}^I &= \begin{Bmatrix}
      m \ddot{{\underline{u}}} + (\dot{\tilde{\omega}} + \tilde{\omega} \tilde{\omega})m {\underline{\eta}} \\
      m \tilde{\eta} \ddot{{\underline{u}}} +{\underline{\underline{\varrho}}}\dot{{\underline{\omega}}}+\tilde{\omega}{\underline{\underline{\varrho}}}{\underline{\omega}}
      \end{Bmatrix}

where :math:`m` is the mass density per unit length,
:math:`{\underline{\eta}}` is the location of the sectional center of
mass, :math:`{\underline{\underline{\varrho}}}` is the moment of inertia
tensor, and the following notations were introduced to simplify the
above expressions

.. math::
       :label: E1-PartC

       {\underline{E}}_1 &= {\underline{x}}_0^\prime + {\underline{u}}^\prime \\
       {\underline{\underline{\mathcal{C}}}} &= \begin{bmatrix}
       {\underline{\underline{C}}}_{11} & {\underline{\underline{C}}}_{12} \\
       {\underline{\underline{C}}}_{21} & {\underline{\underline{C}}}_{22}
       \end{bmatrix}

Damping Forces and Linearization
--------------------------------

A viscous damping model has been implemented into BeamDyn to account for
the structural damping effect. The damping force is defined as

.. math::
      :label: Damping

      {\underline{f}}_d = {\underline{\underline{\mu}}}~ {\underline{\underline{\mathcal{C}}}} \begin{Bmatrix}
      \dot{\epsilon} \\
      \dot{\kappa}
      \end{Bmatrix}

where :math:`{\underline{\underline{\mu}}}` is a user-defined
damping-coefficient diagonal matrix. The damping force can be recast in
two separate parts, like :math:`{\underline{\mathcal{F}}}^C` and
:math:`{\underline{\mathcal{F}}}^D` in the elastic force, as

.. math::
      :label: DampingForce-1-2

      {\underline{\mathcal{F}}}^C_d &= \begin{Bmatrix}
      {\underline{F}}_d \\
      {\underline{M}}_d
      \end{Bmatrix} \\
      {\underline{\mathcal{F}}}^D_d &= \begin{Bmatrix}
       {\underline{0}} \\
       (\tilde{x}^\prime_0 + \tilde{u}^\prime)^T \underline{F}_d
       \end{Bmatrix}

The linearization of the structural damping forces are as follows:

.. math::
       :label: DampingForceLinear-1-2

       \Delta {\underline{\mathcal{F}}}^C_d &= {\underline{\underline{\mathcal{S}}}}_d \begin{Bmatrix}
       \Delta {\underline{u}}^\prime \\
       \Delta {\underline{c}}^\prime
       \end{Bmatrix} + {\underline{\underline{\mathcal{O}}}}_d \begin{Bmatrix}
       \Delta {\underline{u}} \\
       \Delta {\underline{c}}
       \end{Bmatrix} + {\underline{\underline{\mathcal{G}}}}_d \begin{Bmatrix}
       \Delta {\underline{\dot{u}}} \\
       \Delta {\underline{\omega}}
       \end{Bmatrix}     + {\underline{\underline{\mu}}} ~{\underline{\underline{C}}} \begin{Bmatrix}
       \Delta {\underline{\dot{u}}}^\prime \\
       \Delta {\underline{\omega}}^\prime
       \end{Bmatrix} \\
       \Delta {\underline{\mathcal{F}}}^D_d &= {\underline{\underline{\mathcal{P}}}}_d \begin{Bmatrix}
       \Delta {\underline{u}}^\prime \\
       \Delta {\underline{c}}^\prime
       \end{Bmatrix} + {\underline{\underline{\mathcal{Q}}}}_d \begin{Bmatrix}
       \Delta {\underline{u}} \\
       \Delta {\underline{c}}
       \end{Bmatrix} + {\underline{\underline{\mathcal{X}}}}_d \begin{Bmatrix}
       \Delta {\underline{\dot{u}}} \\
       \Delta {\underline{\omega}}
       \end{Bmatrix}     + {\underline{\underline{\mathcal{Y}}}}_d \begin{Bmatrix}
       \Delta {\underline{\dot{u}}}^\prime \\
       \Delta {\underline{\omega}}^\prime
       \end{Bmatrix}

where the newly introduced matrices are defined as

.. math::
       :label: DampingSd-Od-Gd-Pd-Qd-Xd-Yd

       {\underline{\underline{\mathcal{S}}}}_d &=
       {\underline{\underline{\mu}}} {\underline{\underline{\mathcal{C}}}} \begin{bmatrix}
       \tilde{\omega}^T & {\underline{\underline{0}}} \\
       {\underline{\underline{0}}} & \tilde{\omega}^T
       \end{bmatrix} \\
       {\underline{\underline{\mathcal{O}}}}_d &=
       \begin{bmatrix}
       {\underline{\underline{0}}} & {\underline{\underline{\mu}}} {\underline{\underline{C}}}_{11} (\dot{\tilde{u}}^\prime - \tilde{\omega} \tilde{E}_1) - \tilde{F}_d \\
       {\underline{\underline{0}}} &{\underline{\underline{\mu}}} {\underline{\underline{C}}}_{21} (\dot{\tilde{u}}^\prime - \tilde{\omega} \tilde{E}_1) - \tilde{M}_d
       \end{bmatrix} \\
       {\underline{\underline{\mathcal{G}}}}_d &=
       \begin{bmatrix}
       {\underline{\underline{0}}} & {\underline{\underline{C}}}_{11}^T {\underline{\underline{\mu}}}^T \tilde{E}_1 \\
       {\underline{\underline{0}}} & {\underline{\underline{C}}}_{12}^T {\underline{\underline{\mu}}}^T \tilde{E}_1
       \end{bmatrix} \\
       {\underline{\underline{\mathcal{P}}}}_d &=
       \begin{bmatrix}
       {\underline{\underline{0}}} & {\underline{\underline{0}}}  \\
       \tilde{F}_d + \tilde{E}_1^T {\underline{\underline{\mu}}} {\underline{\underline{C}}}_{11} \tilde{\omega}^T &  \tilde{E}_1^T {\underline{\underline{\mu}}} {\underline{\underline{C}}}_{12} \tilde{\omega}^T
       \end{bmatrix} \\
       {\underline{\underline{\mathcal{Q}}}}_d &=
       \begin{bmatrix}
       {\underline{\underline{0}}} & {\underline{\underline{0}}}  \\
       {\underline{\underline{0}}} &  \tilde{E}_1^T {\underline{\underline{O}}}_{12}
       \end{bmatrix} \\
       {\underline{\underline{\mathcal{X}}}}_d &=
       \begin{bmatrix}
       {\underline{\underline{0}}} & {\underline{\underline{0}}}  \\
        {\underline{\underline{0}}} &  \tilde{E}_1^T {\underline{\underline{G}}}_{12}
       \end{bmatrix} \\
       {\underline{\underline{\mathcal{Y}}}}_d &=
       \begin{bmatrix}
       {\underline{\underline{0}}} & {\underline{\underline{0}}}  \\
         \tilde{E}_1^T {\underline{\underline{\mu}}} {\underline{\underline{C}}}_{11} &   \tilde{E}_1^T {\underline{\underline{\mu}}} {\underline{\underline{C}}}_{12}
       \end{bmatrix} \\

where :math:`{\underline{\underline{O}}}_{12}` and
:math:`{\underline{\underline{G}}}_{12}` are the :math:`3 \times 3` sub
matrices of :math:`\mathcal{{\underline{\underline{O}}}}` and
:math:`\mathcal{{\underline{\underline{G}}}}` as
:math:`{\underline{\underline{C}}}_{12}` in Eq. :eq:`E1-PartC`.

.. _convergence-criterion:

Convergence Criterion and Generalized-\ :math:`\alpha` Time Integrator
----------------------------------------------------------------------

The system of nonlinear equations in Eqs. :eq:`GovernGEBT-1-2` are solved using the
Newton-Raphson method with the linearized form in Eq. :eq:`LinearizedEqn`. In the present
implementation, an energy-like stopping criterion has been chosen, which
is calculated as

.. math::
       :label: StoppingCriterion

       | \Delta \mathbf{U}^{(i)T} \left( {^{t+\Delta t}} \mathbf{R} -  {^{t+\Delta t}} \mathbf{F}^{(i-1)}  \right) | \leq | \epsilon_E \left( \Delta \mathbf{U}^{(1)T} \left( {^{t+\Delta t}} \mathbf{R} - {^t}\mathbf{F} \right) \right) |

where :math:`|\cdot|` denotes the absolute value,
:math:`\Delta \mathbf{U}` is the incremental displacement vector,
:math:`\mathbf{R}` is the vector of externally applied nodal point
loads, :math:`\mathbf{F}` is the vector of nodal point forces
corresponding to the internal element stresses, and :math:`\epsilon_E`
is the user-defined energy tolerance. The superscript on the left side
of a variable denotes the time-step number (in a dynamic analysis),
while the one on the right side denotes the Newton-Raphson iteration
number. As pointed out by :cite:`Bathe-Cimento:1980`, this
criterion provides a measure of when both the displacements and the
forces are near their equilibrium values.

Time integration is performed using the generalized-\ :math:`\alpha`
scheme in BeamDyn, which is an unconditionally stable (for linear
systems), second-order accurate algorithm. The scheme allows for users
to choose integration parameters that introduce high-frequency numerical
dissipation. More details regarding the generalized-\ :math:`\alpha`
method can be found in :cite:`Chung-Hulbert:1993,Bauchau:2010`.

Calculation of Reaction Loads
-----------------------------

Since the root motion of the wind turbine blade, including displacements
and rotations, translational and angular velocities, and translational
and angular accelerates, are prescribed as inputs to BeamDyn either by
the driver (in stand-alone mode) or by FAST glue code (in FAST-coupled
mode), the reaction loads at the root are needed to satisfy equality of
the governing equations. The reaction loads at the root are also the
loads passing from blade to hub in a full turbine analysis.

The governing equations in Eq. :eq:`GovernGEBT-1-2` can be recast in a compact form

.. math::
   :label: CompactGovern

   {\underline{\mathcal{F}}}^I - {\underline{\mathcal{F}}}^{C\prime} + {\underline{\mathcal{F}}}^D = {\underline{\mathcal{F}}}^{ext}

with all the vectors defined in Section [sec:LinearProcess]. At the
blade root, the governing equation is revised as

.. math::
   :label: CompactGovernRoot

   {\underline{\mathcal{F}}}^I - {\underline{\mathcal{F}}}^{C\prime} + {\underline{\mathcal{F}}}^D = {\underline{\mathcal{F}}}^{ext}+{\underline{\mathcal{F}}}^R

where :math:`{\underline{\mathcal{F}}}^R = \left[ {\underline{F}}^R~~~{\underline{M}}^R\right]^T`
is the reaction force vector and it can be solved from
Eq. :eq:`CompactGovernRoot` given that the motion fields are known at this
point.

Calculation of Blade Loads
--------------------------

BeamDyn can also calculate the blade loads at each finite element node
along the blade axis. The governing equation in Eq. :eq:`CompactGovern` are
recast as

.. math::
   :label: GovernBF

   {\underline{\mathcal{F}}}^A + {\underline{\mathcal{F}}}^V - {\underline{\mathcal{F}}}^{C\prime} + {\underline{\mathcal{F}}}^D = {\underline{\mathcal{F}}}^{ext}

where the inertial force vector :math:`{\underline{\mathcal{F}}}^I` is
split into :math:`{\underline{\mathcal{F}}}^A` and :math:`{\underline{\mathcal{F}}}^V`:

.. math::
       :label: mathcalFA-FV

       {\underline{\mathcal{F}}}^A &= \begin{Bmatrix}
       m \ddot{{\underline{u}}} + \dot{\tilde{\omega}}m {\underline{\eta}}\\
       m \tilde{\eta} \ddot{{\underline{u}}} + {\underline{\underline{\rho}}} \dot{{\underline{\omega}}}
       \end{Bmatrix} \\
       {\underline{\mathcal{F}}}^V &= \begin{Bmatrix}
       \tilde{\omega} \tilde{\omega} m {\underline{\eta}}\\
        \tilde{\omega} {\underline{\underline{\rho}}} {\underline{\omega}}
       \end{Bmatrix} \\

The blade loads are thus defined as

.. math::
   :label: BladeForce

   {\underline{\mathcal{F}}}^{BF} \equiv {\underline{\mathcal{F}}}^V - {\underline{\mathcal{F}}}^{C\prime} + {\underline{\mathcal{F}}}^D

We note that if structural damping is considered in the analysis, the
:math:`{\underline{\mathcal{F}}}^{C}_d` and
:math:`{\underline{\mathcal{F}}}^D_d` are incorporated into the internal
elastic forces, :math:`{\underline{\mathcal{F}}}^C` and
:math:`{\underline{\mathcal{F}}}^D`, for calculation.
