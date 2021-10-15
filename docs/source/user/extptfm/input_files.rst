.. _ep_input-files:

Input Files
===========

The different input files used by ExtPtfm are described in this section.  ExtPtfm uses the SI system (kg, m, s, N)
No lines should be added or removed from the input files.



.. _ep_main-input-file:

ExtPftm Module Input File
-------------------------

Prior to OpenFAST 2.3, the ExtPtfm module had no "module" input file, and the Guyan ASCII input file was given directly. This is no longer supported, and a module input file is required.
An example of `OpenFAST setup with ExtPtfm is available here <https://github.com/OpenFAST/r-test/blob/main/glue-codes/openfast/5MW_OC4Jckt_ExtPtfm/>`_.
An example of `ExtPtfm module input file is available here <https://github.com/OpenFAST/r-test/blob/main/glue-codes/openfast/5MW_OC4Jckt_ExtPtfm/ExtPtfm.dat>`_.
The format is similar to other *OpenFAST* module.
The input parameters are:

-  ``DT``: Time step for numerical integration (s). The user may specify
   a time step here or use “default” to rely on the glue-code time-step.

-  ``IntMethod``: numerical method for the time integration. The
   Runge-Kutta, Adams–Bashforth and Adams–Bashforth-Moutlon methods are
   available.

-  ``FileFormat``: file format used for the reduction inputs. Available
   formats are ``GuyanASCII`` (see :ref:`epGuyanInputFile`) and ``FlexASCII`` (see :ref:`epSuperelementInputFile`)

-  ``Red_Filename``: path of file containing the Guyan/Craig-Bampton
   inputs.

-  ``RedCst_Filename``: path of file containing the Guyan/Craig-Bampton
   inputs that are constant. This input is not used yet but may be
   introduced in the future to accommodate for reduction file formats
   that use two files: one that contains the constants that are
   structure dependent (static loads from e.g. gravity and matrices),
   and one that contain the time-varying signals that are simulation
   dependent (e.g. loads (on top of the constants loads) and wave
   elevation)

-  ``ActiveDOFList``: list of size ``NActiveDOFList`` containing the CB
   modes indices that are active. This list is not read if
   ``NActiveDOFList<=0``. When specified, all system matrices are reshaped as
   :math:`\boldsymbol{M}_\text{new}=\boldsymbol{M}(I,I)` where :math:`I`
   is the list of indices, potentially unsorted and noncontiguous.
   Setting ``NActiveDOFList=0`` is equivalent to a Guyan reduction.
   Setting :math:`\texttt{NActiveCBDOF}=-1` uses all the CB DOF provided
   in the input file.

-  ``InitPosList``: list of size ``NInitPosList`` containing the initial
   positions for the CB modes. This list is not read if
   ``NInitPosList<=0``, in which case all the CB DOF positions are
   initialized to 0.

-  ``InitVelList``: list of size ``NInitVelList`` containing the initial
   velocities for the CB modes. This list is not read if
   ``NInitVelList<=0``, in which case all the CB DOF velocities are
   initialized to 0.

 - ``SumPrint``:  Print summary data to <RootName>.sum 

 - ``OutFile`` , ``TabDelim``, ``OutFmt``, ``TStart``: Output flags, currently unused

-  ``OutList``: Specifies the list of outputs that the user requests.
   These outputs are described in :ref:`epOutputChannels`.



.. _epOutputChannels:

Output channels
---------------



Outputs are written to disk via the ’.out’ or ’.outb’ files exported 
by *OpenFAST*.
The time series of loads and displacements at the
interface can be selected in ElastoDyn (e.g. ``PtfmPitch``)
Additional “write outputs” are 
implemented in *ExtPtfm*, according to the
list given below. The symbols used in the theory section (:ref:`ep-theory`) are also given in the table.

.. table:: Output channels for the *ExtPtfm* module

   ================ ======================================================================== ==================================== =========
   **Channel name** **Description**                                                          **Symbol**                           **Units**
   ================ ======================================================================== ==================================== =========
   ``IntrfFx``      Platform interface force - Directed along the x-direction                :math:`f_{C}[1]`                     (N)
   ``IntrfFy``      Platform interface force - Directed along the y-direction                :math:`f_{C}[2]`                     (N)
   ``IntrfFz``      Platform interface force - Directed along the z-direction                :math:`f_{C}[3]`                     (N)
   ``IntrfMx``      Platform interface moment - Directed along the x-direction               :math:`f_{C}[4]`                     (Nm)
   ``IntrfMy``      Platform interface moment - Directed along the y-direction               :math:`f_{C}[5]`                     (Nm)
   ``IntrfMz``      Platform interface moment - Directed along the z-direction               :math:`f_{C}[6]`                     (Nm)
   ``InpF_Fx``      Reduced input force at interface point - Directed along the x-direction  :math:`f_{r1}[1]`                    (N)
   ``InpF_Fy``      Reduced input force at interface point - Directed along the y-direction  :math:`f_{r1}[2]`                    (N)
   ``InpF_Fz``      Reduced input force at interface point - Directed along the z-direction  :math:`f_{r1}[3]`                    (N)
   ``InpF_Mx``      Reduced input moment at interface point - Directed along the x-direction :math:`f_{r1}[4]`                    (Nm)
   ``InpF_My``      Reduced input moment at interface point - Directed along the y-direction :math:`f_{r1}[5]`                    (Nm)
   ``InpF_Mz``      Reduced input moment at interface point - Directed along the z-direction :math:`f_{r1}[6]`                    (Nm)
   ``CBQ_XXX``      Displacement of CB DOF number XXX (e.g. ``CBQ_001``)                     :math:`\boldsymbol{x}_2[XXX]`        (-)
   ``CBQD_XXX``     Velocity of CB DOF number XXX     (e.g. ``CBQD_001``)                    :math:`\boldsymbol{\dot{x}}_2[XXX]`  (-)
   ``CBQD2_XXX``    Acceleration of CB DOF number XXX                                        :math:`\boldsymbol{\ddot{x}}_2[XXX]` (-)
   ``CBF_XXX``      Reduced input modal force in CB DOF number XXX                           :math:`\boldsymbol{f}_{r2}[XXX]`     (-)
   ``WaveElevExt``  Wave elevation provided in the external file                             :math:`\eta`                         (m)
   ================ ======================================================================== ==================================== =========



.. _epGuyanInputFile:

Guyan input file (``GuyanASCII``)
---------------------------------

The Guyan input files format is a legacy file format used for superelements
that only contains 6 interface degrees of freedom.

Example
^^^^^^^

An example of ASCII file that was used for Guyan-Reduced sub-structure
is given below, where numerical values are implied instead of ``M11``,
``t1``, ``Fx1`` etc.

.. code::

   Comment
   #Mass
     M11 ... M16
        [...]
     M61 ... M66
   #Damping
     [6 x 6 matrix]
   #Stiffness
     [6 x 6 matrix]
   # time-varying force
   # Time  Fx  Fy  Fz   Mx   My    Mz 
   # s    (N) (N) (N) (N-m) (N-m) (N-m)
     t1   Fx1 Fy1 Fz1  Mx1   My1   Mz1
                   [...]
     tN   FxN FyN FzN  MxN   MyN   MzN

Specifications
^^^^^^^^^^^^^^

The format is a fixed form format
where the line number are assumed. The format specifications are defined
below:

-  ASCII file

-  Line 1 is an arbitrary comment

-  Line 2 must contain ‘\ ``#mass``\ ‘ (case insensitive). If not, the
   file format is invalidated. This requirement is important to
   differentiate between this format and other ASCII formats.

-  Lines: 9 and 16 are comment lines that are ignored

-  Lines 3-8, 10-15, 17-22 contain six float values, forming the
   elements of the mass, damping and stiffness matrices respectively.
   These values corresponds to :math:`\boldsymbol{M}_r`,
   :math:`\boldsymbol{K}_r` and :math:`\boldsymbol{D}_r`.

-  Lines 23-25 are comment lines and are ignored

-  The remaining lines of the files contain 7 float values,
   corresponding to the values: :math:`t`, :math:`[F_x(t)`,
   :math:`F_y(t)`, :math:`F_z(t)`, :math:`M_x(t)`, :math:`M_y(t)`,
   :math:`M_z(t)] = \boldsymbol{f}_{r1}(t)`. The number of time steps is
   here noted :math:`N` but it is not specified in the file.

In particular the same file format
should be used for Guyan and Craig-Bampton reduced substructures. The
following sections define formats that can serve for both purposes.

.. _epSuperelementInputFile:

Superelement input file (``FlexASCII``)
---------------------------------------

This superelement input file is used to provide the system matrices and time series of loads 
for superelements with an arbitrary number of Craig-Bampton modes.

Example
^^^^^^^


An example of superelement file is available `here <https://github.com/OpenFAST/r-test/blob/main/glue-codes/openfast/5MW_OC4Jckt_ExtPtfm/ExtPtfm_SE.dat>`_.
A "dummy" example is given below, where numerical values are
implied for: ``n``, ``dt``, ``t``, ``M11``, ``F1`` etc.

.. code::

   !Comment
   !Comment Flex 5 Format
   !Dimension:                         n
   !Time increment in simulation:      dt
   !Total simulation time in file:     T
   !Mass Matrix (Units (kg,m))
   !Dimension:                         n
     M11 ... M1n
        [...]
     Mn1 ... Mnn
   !Stiffness Matrix (Units (N,m))
   !Dimension:                         n
     [n x n matrix]
   !Damping Matrix (Units (N,m,kg))
   !Dimension:                         n
     [n x n matrix]
   !Loading and Wave Elevation (Units (N,m))
   !Dimension: 1 time column -  n force columns - 1 wave elevation column
     t1   F11 ... F1n eta1
            [...]
     tN   F1N ... FnN etaN


Specifications
^^^^^^^^^^^^^^

The file follows the following specifications:

-  ASCII file

-  Line 1: arbitrary comment that needs to start with an exclamation
   mark ‘\ ``!``\ ‘

-  Line 2: comment which must contain the string ‘\ ``Flex 5 format``\ ‘
   (case insensitive)

-  The following lines are header lines that should start with an
   exclamation mark.

-  The header lines are either comments or lines containing
   keyword/value pairs

-  The following (case insensitive) keywords are currently supported for
   the header:

   -  ‘\ ``!dimension``:‘ followed by the integer
      ``n``\ :math:`=6+n_{CB}`

   -  ‘\ ``!time increment in simulation:``\ ‘: followed by the time
      step :math:`dt`

   -  ‘\ ``!total simulation time in file:``\ ‘: following by the
      simulation length :math:`T`

-  The remaining lines consists of the following special (case
   insensitive) keywords:

   -  ‘\ ``!mass matrix``\ ‘: followed by some text. The next line
      provide a dimension, but it is ignored. The dimension line is then
      followed by :math:`n` lines each containing :math:`n` float values
      These values corresponds to :math:`\boldsymbol{M}_r`.

   -  ‘\ ``!stiffness matrix``\ ‘: similar to the mass matrix, the
      values corresponds to :math:`\boldsymbol{K}_r`.

   -  ‘\ ``!damping matrix``\ ‘: similar to the mass matrix, the values
      corresponds to :math:`\boldsymbol{D}_r`.

   -  ‘\ ``!loading``\ ‘: followed by some text. The next line contains
      the dimensions but is ignored. The remaining lines of the file
      after this keyword should each contain ``n``\ +2 values,
      corresponding to the time :math:`t`, the loads
      :math:`\boldsymbol{f}_r(t)` and the wave elevation
      :math:`\eta(t)` NOTE: the wave elevation is intended for outputs
      only, but it is not outputed yet.
      The number of lines :math:`N` should be
      consistent with the definition of :math:`dt` and :math:`T` from
      the header. The inputs are linearly  interpolated if the time step is
      different from the time step of *ExtPtfm*.

-  For now, the units information and the dimension information after
   the keywords are ignored. Only the dimension provided in the header
   is read and should be respected throughout the file. The reason for
   discarding these information is that at the time of writing, there is
   no guarantee that this information is always provided, and the format
   specifications of the units and dimension were not specified.






