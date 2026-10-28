.. _MD-SlipSupport-Theory:

Slip Support Theory
===================

This page describes the formulation of the slip supports in MoorDyn-F. Input and usage are described in
:ref:`MD-SlipSupport`.


Notation
--------

A line has nodes :math:`j = 0,\dots,N` with positions :math:`\mathbf{r}_j` and velocities
:math:`\dot{\mathbf{r}}_j`. Segment :math:`k = 1,\dots,N` connects nodes :math:`k-1` and :math:`k` and has
unstretched length :math:`l_k`. The unstretched arc-length coordinate of node :math:`j` is

.. math::

   s_0 = 0, \qquad s_j = s_{j-1} + l_j, \qquad L = s_N .

The hoop is located at the position :math:`\mathbf{x}_h` of its Point and moves with velocity
:math:`\dot{\mathbf{x}}_h`. The contact location on the line is the unstretched arc length :math:`s_c`
measured from end A. All directions along the line are oriented from end A toward end B. :math:`\Delta t`
denotes the requested mooring time step ``dtM``.


Contact Kernel
--------------

The hoop acts on a virtual contact point that is a weighted average of nearby line nodes. The weights come
from a cubic B-spline kernel of fixed width in unstretched arc length,

.. math::

   B(t) =
   \begin{cases}
   \left(4 - 6t^2 + 3|t|^3\right)/6, & |t| < 1 \\
   \left(2 - |t|\right)^3/6,         & 1 \le |t| < 2 \\
   0,                                & |t| \ge 2
   \end{cases}

with kernel scale :math:`W = \mathrm{Width}\cdot h`, where :math:`h` is the longest segment length among
all lines through the same hoop Point. The unnormalized weight of node :math:`j` is

.. math::

   a_j(s_c) = B\!\left(\frac{s_j - s_c}{W}\right)\Delta s_j ,

where :math:`\Delta s_j = (l_j + l_{j+1})/2` is the line length represented by node :math:`j`, with
:math:`\Delta s_0 = l_1` and :math:`\Delta s_N = l_N`.

Near the line ends the kernel extends beyond the line. There, virtual nodes are created by point
reflection about the end node: a virtual node at arc length :math:`-s_j` with position
:math:`2\mathbf{r}_0 - \mathbf{r}_j` near end A, and a virtual node at :math:`2L - s_{N-j}` with position
:math:`2\mathbf{r}_N - \mathbf{r}_{N-j}` near end B. The weight :math:`a` of a virtual node is folded back
onto the real nodes, as :math:`+2a` on the end node and :math:`-a` on the reflected node. The normalized
weights are

.. math::

   w_j(s_c) = \frac{\tilde a_j(s_c)}{A(s_c)}, \qquad A(s_c) = \sum_j \tilde a_j(s_c) ,

where :math:`\tilde a_j` includes the folded contributions. The contact point, its material velocity, and
the unit tangent are

.. math::

   \mathbf{C}(s_c) = \sum_j w_j \mathbf{r}_j, \qquad
   \mathbf{V}(s_c) = \sum_j w_j \dot{\mathbf{r}}_j, \qquad
   \hat{\mathbf{t}} = \frac{\mathbf{C}'}{\lVert\mathbf{C}'\rVert}, \quad
   \mathbf{C}' = \sum_j \frac{\mathrm{d}w_j}{\mathrm{d}s_c}\mathbf{r}_j ,

with :math:`\mathrm{d}w_j/\mathrm{d}s_c = (\tilde a_j' - w_j A')/A`.

For fixed node positions, :math:`\mathbf{C}(s)` is a smooth curve through the line nodes parameterized by
the unstretched arc length :math:`s`, and :math:`\mathbf{C}'` is its derivative with respect to :math:`s`
at the current node positions. Only the weights depend on :math:`s`; the node positions are held fixed.
:math:`\hat{\mathbf{t}}` is therefore the tangent of the smoothed curve at the contact point rather than
the direction of any single segment, and it varies continuously as the contact moves past nodes. The
weight derivatives follow from the kernel,
:math:`a_j' = -B'\!\left((s_j - s_c)/W\right)\Delta s_j/W`, with the same folding at the line ends. If
:math:`\lVert\mathbf{C}'\rVert` vanishes, the direction of the contact segment
:math:`\mathbf{r}_k - \mathbf{r}_{k-1}` is used instead.

Properties of the kernel:

- The weights sum to one, so a rigid translation of the line moves the contact point by the same amount,
  and the reaction on the hoop balances the forces applied to the line.
- With the reflected virtual nodes, :math:`\mathbf{C}(0) = \mathbf{r}_0` and :math:`\mathbf{C}(L) = \mathbf{r}_N`.
- The kernel has an RMS spread of :math:`W/\sqrt{3}` along the line.
- :math:`\mathbf{C}` is a smoothed version of the line. At a bend with curvature :math:`\kappa`, it lies
  inside the line by approximately :math:`W^2\kappa/6`, so the line passes slightly outside the hoop.
- On a line with uniform segments away from the ends, the unnormalized weight sum :math:`A` and the mean and
  variance of the weight distribution do not depend on the contact position when :math:`W/h` is an integer.
  The node-passing ripple of the contact force is governed mainly by the variation of :math:`\sum_j w_j^2`
  as the kernel shifts along the evenly spaced nodes. This variation repeats every segment, and its relative
  range over one segment decreases rapidly with :math:`W/h`:

  .. list-table::
     :header-rows: 1

     * - :math:`W/h`
       - 1
       - 1.5
       - 2
       - 3
       - 6
     * - Variation of :math:`\sum_j w_j^2` (max - min)/mean
       - 8.3%
       - 0.41%
       - 0.27%
       - 0.043%
       - 0.0024%


Contact Location
----------------

The contact location minimizes the distance between the hoop and the contact curve,

.. math::

   s_c = \arg\min_{s}\ \lVert \mathbf{C}(s) - \mathbf{x}_h \rVert^2 .

It is parameterized by segment :math:`k` and a local coordinate :math:`u \in [0,1]` as
:math:`s_c = s_{k-1} + u\,l_k`. For each candidate segment, :math:`u` is initialized by projecting
:math:`\mathbf{x}_h` onto the chord :math:`\mathbf{r}_{k-1}\mathbf{r}_k` and refined by Gauss-Newton
iterations

.. math::

   u \leftarrow \operatorname{clamp}_{[0,1]}\!\left(u - \frac{(\mathbf{C} - \mathbf{x}_h)\cdot\mathbf{C}_u}
   {\mathbf{C}_u\cdot\mathbf{C}_u}\right), \qquad \mathbf{C}_u = l_k\,\mathbf{C}' .

The candidate segments are the segment found in the previous evaluation and two segments on either side
of it (all segments on the first evaluation). The candidate with the smallest distance is used. The contact
location is not a state; it is recomputed from the line and hoop positions in every evaluation.

If the contact moves by more than two segments between evaluations, the search stops at the edge of the
window, so the contact lags behind the hoop and the contact force is applied at the wrong location. The
window follows the contact by up to two segments per evaluation, so the lag clears after a few evaluations.
This is unlikely in practice: with an explicit time step limited by the axial wave speed :math:`c` of the
line (:math:`\Delta t \lesssim h/c`), it would require slip speeds comparable to :math:`c`.


Normal Contact Force
--------------------

The contact offset and relative velocity perpendicular to the line are

.. math::

   \boldsymbol{\delta} = (\mathbf{C} - \mathbf{x}_h) - \left[(\mathbf{C} - \mathbf{x}_h)\cdot\hat{\mathbf{t}}\right]\hat{\mathbf{t}},
   \qquad
   \mathbf{v}_{rel} = \mathbf{V} - \dot{\mathbf{x}}_h,
   \qquad
   \mathbf{v}_\perp = \mathbf{v}_{rel} - (\mathbf{v}_{rel}\cdot\hat{\mathbf{t}})\,\hat{\mathbf{t}} .

The transverse penalty force acting on the contact point is

.. math::

   \mathbf{F}_\perp = -K\boldsymbol{\delta} - c\,\mathbf{v}_\perp, \qquad c = 2\zeta\sqrt{K m_C},

where :math:`m_C` is the inertia felt at the contact point when a force is shared among the nodes,

.. math::

   m_C = \frac{m_{node}}{\sum_j w_j^2}, \qquad
   m_{node} = \frac{L}{N}\left(m' + \rho_w C_{an} A_d\right),

with line mass per unit length :math:`m'`, water density :math:`\rho_w`, transverse added-mass
coefficient :math:`C_{an}`, and :math:`A_d = \pi d^2/4`.

The expression for :math:`m_C` follows from treating the nodes as independent masses on the short time scale
of the contact. A force :math:`\mathbf{F}` at the contact point is applied to node :math:`j` as
:math:`w_j\mathbf{F}`, which accelerates the node by :math:`w_j\mathbf{F}/m_{node}`. The contact point moves
with the weighted node motion, so its acceleration is
:math:`\sum_j w_j\,(w_j\mathbf{F}/m_{node}) = \mathbf{F}\sum_j w_j^2/m_{node}`, which gives the effective mass
:math:`m_C`. Because the weights are nonnegative and sum to one, :math:`\sum_j w_j^2 \le 1`, and therefore

.. math::

   m_C \ge m_{node} .

Equality holds when the whole weight sits on a single node, such as at a line end. Spreading the force over
more nodes makes the contact point heavier; for a uniform line,
:math:`\sum_j w_j^2 \approx 0.48\,h/W`, so :math:`m_C \approx 4\,m_{node}` for :math:`W/h = 2`. Near the line
ends, the negative folded weights can make :math:`\sum_j w_j^2` slightly exceed one, so :math:`m_C` can be
slightly smaller than :math:`m_{node}` there.

:math:`m_C` is recomputed from the current weights in every evaluation, so the damping coefficient :math:`c`
and the friction regularization and bounds that depend on :math:`m_C` (see below) follow the contact as it
slides. :math:`K`, in contrast, is fixed at initialization. When ``K`` = 0 is input,

.. math::

   K = m_{node}\left(\frac{0.3}{\Delta t}\right)^2 ,

which keeps the contact frequency :math:`\sqrt{K/m_C}` at or below about :math:`0.3/\Delta t`, including at
the line ends where the weight concentrates on the end node.


Force Distribution to Nodes
---------------------------

Each node receives its share of the contact force split along its own tangent. The node tangent
:math:`\hat{\mathbf{q}}_j` is the unit bisector of the adjacent segment directions,

.. math::

   \hat{\mathbf{q}}_j = \frac{\hat{\mathbf{e}}_j + \hat{\mathbf{e}}_{j+1}}{\lVert\hat{\mathbf{e}}_j + \hat{\mathbf{e}}_{j+1}\rVert},
   \qquad \hat{\mathbf{e}}_k = \frac{\mathbf{r}_k - \mathbf{r}_{k-1}}{\lVert\mathbf{r}_k - \mathbf{r}_{k-1}\rVert},

using the single adjacent segment at the end nodes (and :math:`\hat{\mathbf{t}}` if the bisector vanishes).
The normal force applied to node :math:`j` is the part of its share perpendicular to :math:`\hat{\mathbf{q}}_j`,

.. math::

   \mathbf{f}^N_j = w_j\left[\mathbf{F}_\perp - (\mathbf{F}_\perp\cdot\hat{\mathbf{q}}_j)\,\hat{\mathbf{q}}_j\right],
   \qquad
   F_N = \left\lVert \sum_j \mathbf{f}^N_j \right\rVert .

:math:`F_N` is the normal force actually delivered to the line and is used for the friction capacity. With
the friction force :math:`F_f` defined below, the total force on node :math:`j` and the reaction on the hoop
are

.. math::

   \mathbf{f}_j = \mathbf{f}^N_j + w_j F_f\,\hat{\mathbf{q}}_j, \qquad
   \mathbf{F}_h = -\sum_j \mathbf{f}_j .

The node forces :math:`\mathbf{f}_j` are added to the net forces of the line nodes; forces on end nodes are
passed to the attached objects as usual. :math:`\mathbf{F}_h` is added to the net force of the hoop Point.

Splitting the force along each node's own tangent ensures that a frictionless hoop applies no force along the
line at any node. Without it, the parts of the kernel acting on nodes away from a sharp bend would push the
line along its length, and the hoop could hold tension without friction. The split makes the contact force
non-conservative (it is not the gradient of a contact potential), but the contact damping likely far
outweighs the resulting energy change.


Slip Velocity and Friction Coefficient
--------------------------------------

The slip velocity of the line relative to the hoop is

.. math::

   v = \mathbf{v}_{rel}\cdot\hat{\mathbf{t}} ,

positive when the line material moves toward end B. Because :math:`s_c` is a material coordinate,
:math:`\dot s_c \approx -v/(1+\varepsilon)` with :math:`\varepsilon` the local line strain. The friction
coefficient transitions from static to kinetic friction with the Stribeck velocity :math:`v_S`,

.. math::

   \mu(v) = \mu_k + (\mu_s - \mu_k)\,e^{-(v/v_S)^2} .


.. _MD-SlipSupport-Coulomb:

Regularized Coulomb Friction
----------------------------

When ``StickDisp`` = 0,

.. math::

   F_f = -\mu(v)\,F_N \tanh\!\left(\frac{v}{\tilde v_{tol}}\right),
   \qquad
   \tilde v_{tol} = \max\!\left(v_{tol},\ \frac{\mu_s F_N \Delta t}{m_C}\right).

Near zero slip, the friction acts as a viscous damper with coefficient :math:`\mu_s F_N/\tilde v_{tol}`. The
lower bound on :math:`\tilde v_{tol}`, evaluated with the current :math:`F_N` and :math:`m_C`, keeps this
damping within the stability limit of the explicit time
integration, :math:`(\mu_s F_N/\tilde v_{tol})\,\Delta t/m_C \le 1`.


LuGre Friction
--------------

When ``StickDisp`` :math:`= d > 0`, friction follows the LuGre model (Canudas de Wit et al., 1995) with the
friction normalized by the normal force. The bristle deflection :math:`z` is a continuous state of the
slip support. The bristle stiffness per unit normal force is

.. math::

   \sigma_0 = \frac{\mu_s}{d}, \qquad
   \tilde\sigma_0 = \min\!\left(\sigma_0,\ \frac{m_C}{F_N}\left(\frac{0.5}{\Delta t}\right)^2\right),

where the cap, evaluated with the current :math:`F_N` and :math:`m_C`, keeps the bristle frequency
:math:`\sqrt{\tilde\sigma_0 F_N/m_C}` at or below
:math:`0.5/\Delta t`. The bristle dynamics and friction force are

.. math::

   \dot z = v - \kappa z, \qquad
   \kappa = \min\!\left(\frac{|v|\,\tilde\sigma_0}{\mu(v)},\ \frac{0.5}{\Delta t}\right),

.. math::

   F_f^{*} = -\left(F_N\tilde\sigma_0 z + \sigma_1 F_N\,\dot z\right), \qquad
   \sigma_1 F_N = 2\sqrt{\tilde\sigma_0 F_N\, m_C}\; e^{-(v/v_S)^2} ,

with no viscous term (:math:`\sigma_2 = 0`). The bristle damping :math:`\sigma_1` is critical damping of the
bristle stiffness acting on :math:`m_C`. It is faded out as the slip speed exceeds :math:`v_S`, which keeps the
model dissipative during sliding transients. The friction force is bounded,

.. math::

   F_f = \operatorname{sign}(F_f^{*})\,\min\!\left(|F_f^{*}|,\ F_{max}\right), \qquad
   F_{max} =
   \begin{cases}
   \mu_s F_N,    & \kappa < 0.5/\Delta t \\
   \mu(v)\,F_N,  & \kappa = 0.5/\Delta t
   \end{cases}

The :math:`\mu_s F_N` bound is a safeguard that preserves the frictional lag of the LuGre model at breakaway.
When the relaxation rate :math:`\kappa` is capped (fast sliding with a small ``StickDisp``), :math:`z` can
exceed its steady value and the :math:`\mu(v) F_N` bound restores the kinetic friction level.

Limiting behavior:

- Steady sliding: :math:`z \to \mu(v)\operatorname{sign}(v)/\tilde\sigma_0` and
  :math:`F_f \to -\mu(v) F_N \operatorname{sign}(v)`.
- Sticking (:math:`|z| \ll \mu/\tilde\sigma_0`): :math:`\dot z \approx v`, so :math:`z` is the elastic slip
  and the contact acts as a spring along the line with stiffness :math:`\tilde\sigma_0 F_N = \mu_s F_N/d`.
  Breakaway occurs as :math:`|z|` approaches :math:`d`.
- Like the original LuGre model, small oscillations of :math:`v` relax :math:`z` through the :math:`|v|`
  term, which can produce slow drift below breakaway when the line vibrates (Dupont et al., 2002).


.. _MD-SlipSupport-IC:

Initial Conditions
------------------

**Line layout.** For a line with :math:`n` slip supports, the hoop positions :math:`\mathbf{H}_1,\dots,\mathbf{H}_n`
(in the order the supports are listed) and the line end positions define the polyline
:math:`\mathbf{r}_0, \mathbf{H}_1, \dots, \mathbf{H}_n, \mathbf{r}_N` with cumulative chord length
:math:`c_i`. Internal node :math:`j` is placed on this polyline at chord distance
:math:`(s_j/L)\,c_{n+1}`, and the initial node velocities are zero.

**Dynamic relaxation.** During initial condition generation, the friction force is replaced by a viscous
slip damper and the LuGre state is frozen,

.. math::

   F_f = -c_{IC}\,v, \qquad \dot z = 0, \qquad
   c_{IC} = \min\!\left(2\sqrt{2\,w_{wet}\,M_L},\ \frac{0.5\,m_{node}}{\Delta t}\right),

where :math:`w_{wet} = |m' - \rho_w A_d|\,g` is the wet weight per unit length and
:math:`M_L = L\,(m' + \rho_w C_{at} A_d)` is the line mass including axial added mass.

The first term is an ad-hoc damping that slows sliding at the rate :math:`c_{IC}/M_L`, which is
proportional to :math:`\sqrt{g/L}` for lines much heavier than water. The second term keeps the damper
within the explicit stability limit, :math:`c_{IC}\,\Delta t/m_{node} \le 0.5`. The value of
:math:`c_{IC}` affects only the convergence rate of dynamic relaxation, not the equilibrium it converges to.
:math:`z` is set to zero at the start and at the end of dynamic relaxation.


Release
-------

When the contact point sits at a line end (:math:`k = 1, u = 0` or :math:`k = N, u = 1`), the overshoot of
the hoop beyond the end is

.. math::

   o = -(\mathbf{x}_h - \mathbf{C})\cdot\hat{\mathbf{t}} \ \ \text{(end A)}, \qquad
   o = (\mathbf{x}_h - \mathbf{C})\cdot\hat{\mathbf{t}} \ \ \text{(end B)} .

The slip support releases the line permanently if :math:`o > 0.1\,l_k` in every force evaluation within an
accepted time step. Requiring this in every evaluation prevents releases triggered by intermediate
Runge-Kutta stages or linearization perturbations.


References
----------

- C. Canudas de Wit, H. Olsson, K. J. Åström, and P. Lischinsky, "A new model for control of systems with
  friction," *IEEE Transactions on Automatic Control*, 40(3):419-425, 1995.
- P. Dupont, V. Hayward, B. Armstrong, and F. Altpeter, "Single state elastoplastic friction models,"
  *IEEE Transactions on Automatic Control*, 47(5):787-792, 2002.
