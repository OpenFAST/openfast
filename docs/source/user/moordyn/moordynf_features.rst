.. _MD-F-Features:

MoorDyn-F Specific Features
===========================

This page describes input options that are available in MoorDyn-F (the Fortran implementation of MoorDyn
used by OpenFAST) in addition to those documented in the
`MoorDyn input file documentation <https://moordyn.readthedocs.io/en/latest/inputs.html>`_.
All other sections of the MoorDyn input file (line types, points, lines, failures, options, and outputs) are
unchanged and follow the external documentation.


.. _MD-SlipSupport:

Slip Supports
-------------

A slip support models a line threaded through a hoop, eyelet, or guide that is carried by a MoorDyn Point.
The hoop constrains the line transversely but lets it slide along its length, resisted by friction.
The hoop has no physical size: there is no clearance between the line and the hoop, and no bend radius.
Typical applications include lines running through fairlead guides, mid-line guides on a structure, or a
free weight or buoy riding on a line.

The line keeps its original discretization while it slides: no nodes are added, removed, or moved along
the line. Instead, the hoop acts on a virtual contact point that moves continuously along the line, and
the contact force is spread to the nearby line nodes with a smooth weighting kernel. As a result, the
contact force varies smoothly as the line nodes slide past. The formulation is given in
:ref:`MD-SlipSupport-Theory`.

Main capabilities:

- The hoop can be any Point type: fixed, coupled, free, or attached to a Body. The reaction force from the
  line is added to the Point, so it is included in coupled loads passed to the glue code and in the
  dynamics of free Points and Bodies.
- Several slip supports may act on the same line, and several lines may pass through the same Point.
- Friction is either a regularized Coulomb model or a LuGre model with true sticking. Both use separate
  static and kinetic friction coefficients with a smooth (Stribeck) transition between them.
- If the line slides completely through a hoop, the support releases the line permanently and the
  simulation continues.
- Linearization is supported (see :ref:`MD-SlipSupport-Linearization`).


.. _MD-SlipSupport-Input:

Input
~~~~~

Slip supports are defined in an optional ``SLIP SUPPORTS`` section of the MoorDyn input file. The section
header is a line containing ``---`` and the text ``SLIP SUPPORTS``, followed by a column-name line, a units
line, and one row per slip support. **The section must come after the** ``LINES`` **section**, because it
refers to Line properties.

Example (one hoop on Point 2 threading Line 1, with LuGre friction):

.. code-block:: none

   ---------------------- SLIP SUPPORTS ------------------------------------------------
   ID   Point   Line   Width   K       Zeta   MuS    MuK    VTol    VS      StickDisp
   (#)  (#)     (#)    (h)     (N/m)   (-)    (-)    (-)    (m/s)   (m/s)   (m)
   1    2       1      2       0       0.8    0.6    0.3    0.01    0.1     0.005

.. list-table:: Slip support input columns
   :header-rows: 1
   :widths: 12 10 78

   * - Column
     - Units
     - Description
   * - ``ID``
     - (-)
     - Slip support number. Must be sequential starting from 1.
   * - ``Point``
     - (-)
     - ID of the Point carrying the hoop. Any Point type may be used. A free Point used as a hoop
       should have a nonzero mass (see :ref:`MD-SlipSupport-Usage`).
   * - ``Line``
     - (-)
     - ID of the Line threaded through the hoop.
   * - ``Width``
     - (-)
     - Contact kernel width, given as a multiple of the segment length :math:`h`. Here :math:`h` is the
       longest unstretched segment length among all lines that pass through the same Point, so all lines
       through one hoop share the same physical contact width :math:`W = \mathrm{Width}\cdot h`.
       Must be at least 1. Values below 2 produce a warning because the contact force may ripple as the
       line slides past nodes. **2 is recommended.** :math:`4W` must also be shorter than the
       unstretched length of the line.
   * - ``K``
     - (N/m)
     - Transverse contact stiffness (penalty stiffness) that keeps the line centered in the hoop. Enter 0
       for an automatic value tied to the line node mass and the time step ``dtM`` (recommended).
   * - ``Zeta``
     - (-)
     - Transverse contact damping ratio, relative to critical damping of the contact stiffness. A value of
       about 0.8 is recommended.
   * - ``MuS``
     - (-)
     - Static friction coefficient. Must satisfy ``MuS`` :math:`\ge` ``MuK``. Enter 0 together with
       ``MuK`` = 0 for a frictionless hoop.
   * - ``MuK``
     - (-)
     - Kinetic (sliding) friction coefficient, :math:`\ge 0`.
   * - ``VTol``
     - (m/s)
     - Regularization velocity of the Coulomb friction model. Must be positive when ``MuS`` > 0 and
       ``StickDisp`` = 0. Ignored when LuGre friction is used. The line creeps through the hoop at a speed of
       order ``VTol`` when it should be stuck, so it should be small compared with the sliding speeds of
       interest. MoorDyn raises it internally when needed for numerical stability
       (see :ref:`MD-SlipSupport-Coulomb`).
   * - ``VS``
     - (m/s)
     - Stribeck velocity: the slip speed over which friction drops from the static to the kinetic level.
       Must be positive when ``MuS`` > ``MuK`` or when LuGre friction is used.
   * - ``StickDisp``
     - (m)
     - Pre-sliding displacement of the LuGre friction model: the elastic slip the contact can take before
       it breaks loose, typically a few millimeters. Enter 0 to use the regularized Coulomb model instead.

Friction model selection:

- ``StickDisp`` = 0: regularized Coulomb friction. No additional states. The line can creep slowly through
  the hoop while it should be stuck.
- ``StickDisp`` > 0: LuGre friction. One additional continuous state per slip support. The line truly
  sticks until the static friction limit is reached. If ``StickDisp`` is too small for the time step,
  MoorDyn stiffens the contact less than requested to stay stable and prints a one-time warning. Either
  increase ``StickDisp`` or reduce the time step to remove the warning.


.. _MD-SlipSupport-Outputs:

Output Channels
~~~~~~~~~~~~~~~

Slip support output channels are requested in the ``OUTPUTS`` section with the prefix ``SLIP`` followed by
the slip support ID and a quantity code, for example ``SLIP1FN``.

.. list-table:: Slip support output channels
   :header-rows: 1
   :widths: 18 10 72

   * - Channel
     - Units
     - Description
   * - ``SLIP<n>S``
     - (m)
     - Contact location: unstretched arc length from line end A to the contact point.
   * - ``SLIP<n>VS``
     - (m/s)
     - Slip velocity of the line relative to the hoop, along the line. Positive when the line moves toward
       end B (``SLIP<n>S`` decreasing). This is a physical velocity of the stretched line, not the rate of
       change of ``SLIP<n>S``, which is measured in unstretched arc length:
       d(``SLIP<n>S``)/dt :math:`\approx -` ``SLIP<n>VS`` :math:`/(1+\varepsilon)`, with :math:`\varepsilon`
       the line strain at the contact point.
   * - ``SLIP<n>FN``
     - (N)
     - Magnitude of the normal (transverse) contact force applied to the line. Always :math:`\ge 0`.
   * - ``SLIP<n>FF``
     - (N)
     - Friction force applied to the line, along the line. Positive pushes the line toward end B. It
       normally has the opposite sign of ``SLIP<n>VS``.
   * - ``SLIP<n>FX``, ``SLIP<n>FY``, ``SLIP<n>FZ``
     - (N)
     - Total force exerted by the line on the hoop in global coordinates. This is the reaction added to the
       Point.
   * - ``SLIP<n>Z``
     - (m)
     - LuGre bristle deflection. Zero when Coulomb friction is used.

After a slip support is released, ``SLIP<n>VS``, ``SLIP<n>FN``, ``SLIP<n>FF``, and ``SLIP<n>FX/FY/FZ`` are
zero, while ``SLIP<n>S`` and ``SLIP<n>Z`` keep their last values.


.. _MD-SlipSupport-Usage:

Usage Notes
~~~~~~~~~~~

**Initial line shape.** A line with slip supports is not initialized with the catenary solution. Its
internal nodes are placed at equal distances along a polyline that runs from end A through each slip
support Point, in the order the slip supports of that line are listed, to end B. The initial contact
location therefore follows the ratio of the straight-line distances between the line ends and the hoops.
A hoop does not need to coincide with a node. List the slip supports of a line in order from end A to
end B.

**Initial condition generation.** During dynamic relaxation (``TmaxIC`` > 0), friction is replaced by a
viscous slip damper so that the line can slide to its frictionless equilibrium. The damping coefficient is
set automatically from the line properties (see :ref:`MD-SlipSupport-IC`). LuGre states are held at
zero and reset after dynamic relaxation. Note that the convergence check of the dynamic relaxation only
monitors fairlead tensions, so a line may still be sliding slowly when it ends; friction then holds the
remaining imbalance. Set ``TmaxIC`` = 0 to start a free hoop (for example, a bead on a line) from its
input position.

**Release.** A slip support releases its line permanently once the hoop is more than 0.1 segment
length beyond either line end throughout a whole time step. A message is written to the screen and log
file.

**Choosing Width.** Larger values give smoother contact forces when the line nodes slide past the support,
but spread the load over a longer length and place the line slightly outside the hoop at bends.
``Width`` = 2 removes node-passing ripple for practical purposes. Because ``Width`` is relative to the
segment length, mesh refinement also narrows the contact and shrinks the offset between the line and the
hoop; to keep a fixed physical contact width while refining, scale ``Width`` with the number of
segments.

**Bend angle.** The model is intended for moderate bend angles at the hoop. It can work up to a hairpin
turn, but likely with degraded accuracy. To model a line wrapping around a large-radius sheave or bollard,
use several slip supports distributed around its circumference.

**Free hoop points.** The automatic stiffness ``K`` = 0 is sized from the line node mass. A free hoop Point
that is much lighter than a line node may need an explicit, smaller ``K`` for numerical stability.

**Active tensioning.** Changes of the last segment length by line control channels are not considered when
checking ``Width``.


.. _MD-SlipSupport-Linearization:

Linearization
~~~~~~~~~~~~~

Slip supports may be used with linearization. Each LuGre slip support adds one first-order state (the
bristle deflection) named ``Slip <n> bristle z, m``. The linearized model reflects the numerical model of
the contact:

- With Coulomb friction, near-zero slip the friction acts as a viscous damper with coefficient
  :math:`\mu_s F_N / \tilde v_{tol}` along the line, where :math:`\tilde v_{tol} \ge` ``VTol`` is the
  adjusted regularization velocity (see :ref:`MD-SlipSupport-Coulomb`).
- The transverse contact stiffness adds a high-frequency transverse mode. With the automatic stiffness, its
  frequency is about :math:`0.3/dt_M` rad/s.
- With LuGre friction, the stuck contact appears as a stiff spring along the line (the bristle stiffness)
  with an associated fast first-order mode.
