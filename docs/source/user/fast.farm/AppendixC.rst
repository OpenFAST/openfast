.. _FF:App:Output:

List of Output Channels
=======================

This is a list of all possible output parameters available within
FAST.Farm (except those that are available from OpenFAST, which are
specified within the OpenFAST input file(s) and output separately for
each turbine). The names are grouped by meaning, but can be ordered in
the OUTPUTS section of the FAST.Farm primary input file as you see fit.

T\ :math:`\alpha` refers to turbine :math:`\alpha`, where :math:`\alpha`
is a one-digit number in the range [1,9], corresponding to row
:math:`\alpha` in the wind turbine input table. If **NumTurbines** > 9,
only values for the first 9 turbines can be output. Setting
:math:`\alpha` > **NumTurbines** yields invalid output.

In\ :math:`\zeta` and Ot\ :math:`\zeta` refer to super-controller input
and output :math:`\zeta`, respectively, where :math:`\zeta` is a
one-digit number in the range [1,9], corresponding to element
:math:`\zeta` in the input and output arguments of the super-controller
source code. If there are more than 9 elements, only values for the
first 9 inputs and outputs can be output. Setting :math:`\zeta` greater
than the number of elements yields invalid output.

N\ :math:`\beta` refers to radial output node :math:`\beta`, where
:math:`\beta` is a two-digit number in the range [01,20], corresponding
to entry :math:`\beta` in the **OutRadii** list, where node
:math:`\beta` is at radius **dr** :math:`\times`
**OutRadii**\ [:math:`\beta`]. Setting :math:`\beta` > **NOutRadii**
yields invalid output.

W\ :math:`\eta` refers to wind point :math:`\eta`, where :math:`\eta` is
a one-digit number in the range [1,9], corresponding to entry
:math:`\eta`\ in the **WindVelX**, **WindVelY**, and **WindVelZ** lists.
Setting :math:`\eta` > **NWindVel** yields invalid output. Setting
**WindVelX**, **WindVelY**, and **WindVelZ** outside the low-resolution
wind domain also yields invalid output.

:math:`\delta` refers to the X, Y, or Z coordinate axis.

D\ :math:`\gamma` refers to downstream distance :math:`\gamma`, where
:math:`\gamma` is a one-digit number in the range [1,9], corresponding
to entry :math:`\gamma` in the **OutDist** list. Setting :math:`\gamma`
> **NOutDist** yields invalid output. The output is also invalid if
**OutDist** is a distance further downstream than the wake has been
calculated or for any distance where the wake from the turbine has
overlapped itself.

.. container::
   :name: Tab:FF:Outputs

   .. table:: List of Available FAST.Farm Output Channels

      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
      | Channel Name                                                 | Units             | Description                                     |
      +==============================================================+===================+=================================================+
      | *Super Controller*                                                                                                                 |
      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
      | SCGblIn\ :math:`\zeta`                                       | (user)            | Global (turbine independent) super              |
      |                                                              |                   | controller input :math:`\zeta`                  |
      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
      | SCT\ :math:`\alpha`\ In\ :math:`\zeta`                       | (user)            | Turbine-dependent super controller input        |
      |                                                              |                   | :math:`\zeta` for turbine :math:`\alpha`        |
      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
      | SCGblOt\ :math:`\zeta`                                       | (user)            | Global (turbine independent) super              |
      |                                                              |                   | controller output :math:`\zeta`                 |
      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
      | SCT\ :math:`\alpha`\ Ot\ :math:`\zeta`                       | (user)            | Turbine-dependent super controller input        |
      |                                                              |                   | :math:`\zeta` for turbine :math:`\alpha`        |
      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
      | *Wind Turbine and Inflow*                                                                                                          |
      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
      | RtAxs\ :math:`\delta`\ T\ :math:`\alpha`                     | (-)               | Orientation of the rotor centerline for turbine |
      |                                                              |                   | :math:`\alpha` in the global coordinate system  |
      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
      | RtPos\ :math:`\delta`\ T\ :math:`\alpha`                     | \(m\)             | Position of the rotor (hub) center for turbine  |
      |                                                              |                   | :math:`\alpha` in the global coordinate system  |
      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
      | RtDiamT\ :math:`\alpha`                                      | \(m\)             | Rotor diameter for turbine :math:`\alpha`       |
      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
      | YawErrT\ :math:`\alpha`                                      | (deg)             | Nacelle-yaw error for turbine :math:`\alpha`    |
      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
      | TIAmbT\ :math:`\alpha`                                       | (%)               | Ambient turbulence intensity of the wind at the |
      |                                                              |                   | the rotor disk for  turbine :math:`\alpha`. The |
      |                                                              |                   | ambient turbulence  intensity is based on a     |
      |                                                              |                   | spatial-average of the three vector components, |
      |                                                              |                   | instead of just the axial component.            |
      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
      | RtVAmbT\ :math:`\alpha`                                      | (m/s)             | Rotor-disk-averaged ambient wind speed (normal  |
      |                                                              |                   | to disk, not including structural motion, local |
      |                                                              |                   | induction or wakes from upstream turbines) for  |
      |                                                              |                   | turbine :math:`\alpha`                          |
      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
      | RtVRelT\ :math:`\alpha`                                      | (m/s)             | Rotor-disk-averaged relative wind speed (normal |
      |                                                              |                   | to disk, including structural motion and wakes  |
      |                                                              |                   | from upstream turbines, but not including local |
      |                                                              |                   | induction) for turbine :math:`\alpha`           |
      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
      | CtT\ :math:`\alpha`\ N\ :math:`\beta`                        | (-)               | Azimuthally averaged thrust force coefficient   |
      |                                                              |                   | (normal to disk) for radial output node         |
      |                                                              |                   | :math:`\beta` of turbine :math:`\alpha`         |
      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
      | *Wake (for an Individual Rotor)*                                                                                                   |
      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
      | WkAxs\ :math:`\delta`\ T\ :math:`\alpha`\ D\ :math:`\gamma`  | (-)               | Orientation of the wake centerline for          |
      |                                                              |                   | downstream distance :math:`\gamma`  of turbine  |
      |                                                              |                   | :math:`\alpha` in the global coordinate system  |
      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
      | WkPos\ :math:`\delta`\ T\ :math:`\alpha`\ D\ :math:`\gamma`  | \(m\)             | Center position of the wake centerline for      |
      |                                                              |                   | downstream distance :math:`\gamma` of turbine   |
      |                                                              |                   | :math:`\alpha` in the global coordinate system  |
      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
      | WkVel\ :math:`\delta`\ T\ :math:`\alpha`\ D\ :math:`\gamma`  | (m/s)             | Advection, deflection, and meandering velocity  |
      |                                                              |                   | (not including the horizontal wake-deflection   |
      |                                                              |                   | correction or low-pass time-filtering) of the   |
      |                                                              |                   | wake for downstream distance :math:`\gamma` of  |
      |                                                              |                   | turbine :math:`\alpha` in the global coordinate |
      |                                                              |                   | system                                          |
      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
      | WkDiamT\ :math:`\alpha`\ D\ :math:`\gamma`                   | \(m\)             | Wake diameter for downstream distance           |
      |                                                              |                   | :math:`\gamma` of turbine :math:`\alpha`        |
      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
      | WkDfVxT\ :math:`\alpha`\ N\ :math:`\beta`\ D\ :math:`\gamma` | (m/s)             | Axial wake velocity deficits for radial output  |
      |                                                              |                   | node :math:`\beta` and downstream distance      |
      |                                                              |                   | :math:`\gamma` of turbine :math:`\alpha`        |
      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
      | WkDfVrT\ :math:`\alpha`\ N\ :math:`\beta`\ D\ :math:`\gamma` | (m/s)             | Radial wake velocity deficits for radial output |
      |                                                              |                   | node :math:`\beta` and downstream distance      |
      |                                                              |                   | :math:`\gamma` of turbine :math:`\alpha`        |
      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
      | EddVisT\ :math:`\alpha`\ N\ :math:`\beta`\ D\ :math:`\gamma` | (m\ :math:`^2`/s) | Total eddy viscosity for radial output node     |
      |                                                              |                   | :math:`\beta` and downstream distance           |
      |                                                              |                   | :math:`\gamma` of turbine :math:`\alpha`        |
      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
      | EddAmbT\ :math:`\alpha`\ N\ :math:`\beta`\ D\ :math:`\gamma` | (m\ :math:`^2`/s) | Individual contribution to the eddy viscosity   |
      |                                                              |                   | from ambient turbulence for radial output node  |
      |                                                              |                   | :math:`\beta` and downstream distance           |
      |                                                              |                   | :math:`\gamma` of turbine :math:`\alpha`        |
      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
      | EddShrT\ :math:`\alpha`\ N\ :math:`\beta`\ D\ :math:`\gamma` | (m\ :math:`^2`/s) | Individual contributions to the eddy viscosity  |
      |                                                              |                   | from the shear layer for radial output node     |
      |                                                              |                   | :math:`\beta` and downstream distance           |
      |                                                              |                   | :math:`\gamma` of turbine :math:`\alpha`        |
      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
      | *Ambient Wind and Array Effects*                                                                                                   |
      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
      | W\ :math:`\eta`\ VAmb\ :math:`\delta`                        | (m/s)             | Ambient wind velocity (not including wakes) for |
      |                                                              |                   | point :math:`\eta` in the global coordinate     |
      |                                                              |                   | system (from the low-resolution domain)         |
      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
      | W\ :math:`\eta`\ VDis\ :math:`\delta`                        | (m/s)             | Disturbed wind velocity (including wakes) for   |
      |                                                              |                   | point :math:`\eta` in the global coordinate     |
      |                                                              |                   | system (from the low-resolution domain)         |
      +--------------------------------------------------------------+-------------------+-------------------------------------------------+
