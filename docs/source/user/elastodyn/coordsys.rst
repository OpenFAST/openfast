.. _ed_coordsys:


Coordinate systems
==================

For the coordinates system not detailed in subsections below, please refer to the following references:

-  `FAST 7 Manual <https://www.nrel.gov/docs/fy06osti/38230.pdf>`_

- :download:`FASTCoordinateSystems.doc <../../../OtherSupporting/ElastoDyn/FASTCoordinateSystems.doc>`:
  Documents the transformation matrices relating each coordinate system in FAST. Unfortunately, there are no pictures in this document that diagram these coordinate systems. They can hopefully be visualized by means of the transformation matrices.


Tail-Furl coordinate system
---------------------------
The tail-furl DOF allows the user to model the unusual
configuration of a bearing that permits the tail to rotate
about the yawing-portion of the structure atop the
tower. In order to include tail-furling in a model,
the user must designate the turbine as a furling machine by
setting the input ``Furling`` from the ElastoDyn input file to
True. Then you must assemble the furling input file,
``FurlFile``, and use the tail-furl flag, ``TFrlDOF``, to enable
this feature.
The angular tail-furl motion takes place about the
tail-furl axis defined by inputs ``TFrlPntxn``, ``TFrlPntyn``,
``TFrlPntzn``, ``TFrlSkew``, and ``TFrlTilt`` available in
``FurlFile``.
Inputs ``TFrlPntxn``, ``TFrlPntyn``, and ``TFrlPntzn`` locate an arbitrary point on the tail-furl axis
relative to the tower-top. 
Three angles ``TFinAngles`` define the angular orientation of the tail- furl axis passing through this point.  
See :numref:`figTFAxes` for a schematic.

The tail-furl bearing can be an ideal bearing with
no friction by setting ``TFrlMod`` to 0; by setting
``TFrlMod`` to 1, it also has a standard model that
includes a linear spring, linear damper and Coulomb
damper, as well as up- and down-stop springs, and up-
and down-stop dampers. OpenFAST models the stop
springs with a linear function of tail-furl deflection.
The tail-furl stops start at a specified angle and work as
a linear spring based on the deflection past the stop
angles. The tail-furl dampers are linear functions of
the furl rate and start at the specified up-stop and
down-stop angles. These dampers are bidirectional,
resisting motion equally in both directions once past
the stop angle.

A user-defined tail-furl spring and damper model
is also available. To use it, set ``TFrlMod`` to 2 and
create a subroutine entitled ``UserTFrl()`` with the
arguments ``TFrlDef``, ``TFrlRate``, ``ZTime``, ``DirRoot``, and
``TFrlMom``:

- ``TFrlDef``: Current tail-furl angular deflection in radians (input)
- ``TFrlRate``: Current tail-furl angular rate in rad/sec (input)
- ``ZTime``: Current simulation time in sec (input)
- ``DirRoot``: Simulation root name including the full path to the current working directory (input)
- ``TFrlMom``: Tail-furl moment in N.m (output)

The source file UserSubs.f90 contains a dummy
``UserTFrl()`` routine; replace it with your own and
rebuild OpenFAST. 
Argument ``DirRoot`` may be used to write a record of
what is done in ``UserTFrl()`` to be stored with the
simulation results.

The geometries of the tail boom mass center, tail
fin mass center, and tail fin aerodynamic surface,
which are all components of the furling-tail assembly,
are defined relative to the tower-top as shown in :numref:`figTFGeom`.
This definition was chosen in order to avoid
having to define a coordinate system in the furling-tail
assembly since such a coordinate system would most
likely have an obscure orientation, making it difficult
for users to input configuration information relative to
it. This definition also avoids the complications
involved in having to define geometries differently,
depending on whether or not a tail-furl assembly exists
separately from the nacelle, which depends on whether
tail-furl is present or absent in the turbine.

Since the component geometry of the furling-tail
assembly is defined relative to the tower-top, this
geometry naturally changes with the tail-furl angle. In
order to avoid having to define different geometries for
different tail-furl positions (for example, variations in
the initial tail-furl angle), ElastoDyn expects the component
geometry of the furling-tail assembly to be
defined/input at a tail-furl angle of zero. As such, the
initial tail-furl angle does not affect the specification of
any other tail-furl geometry. Stated another way, the
input geometries for the tail-furl assembly components
define the tail configuration when the tail-furl angle is
zero regardless of initial tail-furl position. Users
should be clear of this convention when assembling
their furling input file. Further clarification on this
furling geometry convention is provided in the Rotor-
Furl section above.



.. _figTFAxes:
.. figure:: figs/TailFinAxes.png
   :width: 70%
           
   Layout of a three-bladed, upwind, furling turbine: furl axes


.. _figTFFurl:
.. figure:: figs/TailFinFurl.png
   :width: 70%
           
   Layout of a three-bladed, upwind, furling turbine: rotor-furl structure


.. _figTFGeom:
.. figure:: figs/TailFinGeom.png
   :width: 70%
           
   Layout of a three-bladed, upwind, furling turbine: tail-furl structure.  
   NOTE: The tail fin "CP" (center of pressure) parameters are now replaced by the location of the reference point.



