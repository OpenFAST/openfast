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

   ``RBMod``: method for handling rigid-body motion of floating structures (switch). 
   0: No special handling for rigid-body motion (fixed-bottom structures); 
   1: Transform to rigid-body frame of reference; 
   2: Transform to rigid-body frame of reference and add fictitious forces and exact self-weight.

-  ``Red_Filename``: path of file containing the Guyan/Craig-Bampton
   inputs.

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

   ``Connections``: flag for including connection points on the structure. 
   If true, a file specifiying the connections must be provided through 
   ``Conn_FileName``. Currently, this feature is only used to couple with 
   one of the mooring modules.

   ``UserForcing``: flag for user-defined modal forcing. If true, a file 
   containing the force time series to be applied must be provided through 
   ``Force_FileName``.

   ``ConnForcing``: flag for user-defined external force at connection points.
   If true, a file containing the force time series to be applied at the connection 
   points must be provided through ``FConn_FileName``. This option requires 
   ``Connections`` to be true. Application of moments is not supported.

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
   ``ExtrnFx``      Reduced input force at interface point - Directed along the x-direction  :math:`f_{r1}[1]`                    (N)
   ``ExtrnFy``      Reduced input force at interface point - Directed along the y-direction  :math:`f_{r1}[2]`                    (N)
   ``ExtrnFz``      Reduced input force at interface point - Directed along the z-direction  :math:`f_{r1}[3]`                    (N)
   ``ExtrnMx``      Reduced input moment at interface point - Directed along the x-direction :math:`f_{r1}[4]`                    (Nm)
   ``ExtrnMy``      Reduced input moment at interface point - Directed along the y-direction :math:`f_{r1}[5]`                    (Nm)
   ``ExtrnMz``      Reduced input moment at interface point - Directed along the z-direction :math:`f_{r1}[6]`                    (Nm)
   ``CBDXXX``       Displacement of CB DOF number XXX (e.g. ``CBD001``)                      :math:`\boldsymbol{x}_2[XXX]`        (-)
   ``CBVXXX``       Velocity of CB DOF number XXX     (e.g. ``CBV001``)                      :math:`\boldsymbol{\dot{x}}_2[XXX]`  (-)
   ``CBAXXX``       Acceleration of CB DOF number XXX (e.g. ``CBA001``)                      :math:`\boldsymbol{\ddot{x}}_2[XXX]` (-)
   ``CBFXXX``       Reduced input modal force in CB DOF number XXX (e.g. ``CBF001``)         :math:`\boldsymbol{f}_{r2}[XXX]`     (-)
   ================ ======================================================================== ==================================== =========

.. _epSuperelementInputFile:

Guyan/Craig-Bampton superelement input file (provided through ``Red_Filename``)
-------------------------------------------------------------------------------

This superelement input file is used to provide the system matrices.

Example
^^^^^^^

An example of superelement file is available `here <https://github.com/OpenFAST/r-test/blob/main/glue-codes/openfast/5MW_OC4Jckt_ExtPtfm/ExtPtfm_SE.dat>`_.
A "dummy" example is given below, where numerical values are
implied for: ``n``, ``M(i,j)``, ``K(i,j)``, ``C(i,j)``, etc.

.. code::

   !Comment
   !Comment
   !Dimension: n
   
   !Mass Matrix (Units (kg,m))
     M(1,1) ... M(1,n)
           [...]
     M(n,1) ... M(n,n)
   
   !Stiffness Matrix (Units (N,m))
     K(1,1) ... K(1,n)
           [...]
     K(n,1) ... K(n,n)
   
   !Damping Matrix (Units (N,m,kg))
     C(1,1) ... C(1,n)
           [...]
     C(n,1) ... C(n,n)

   !Weight constant (Units (N,m))
     W_0(1) W_0(2) ... W_0(n)

   !Weight stiffness (Units (N,m))
     K_W(1,1) ... K_W(1,n)
             [...]
     K_W(n,1) ... K_W(n,n)

Specifications
^^^^^^^^^^^^^^

The file follows the following specifications:

-  ASCII file

-  The file can start with an arbitrary number of comment lines that 
   start with an exclamation mark ‘\ ``!``\ ‘

-  The following (case insensitive) keyword must be provided next:

   -  ‘\ ``!dimension``:‘ followed by the integer
      ``n``\ :math:`=6+n_{CB}`

-  The remaining lines consist of the following special (case
   insensitive) keywords:

   -  ‘\ ``!mass matrix``\ ‘: followed by some text. The next :math:`n` 
      lines each contain :math:`n` float values. These values correspond 
      to :math:`\boldsymbol{M}_r`. Note that when modeling a floating 
      structure with ``RBMod`` > 0, the first 6 modes, the interface modes,
      also serve as rigid-body modes. Therefore, the first 6-by-6 entries 
      of :math:`\boldsymbol{M}_r` must be a self-consistent rigid-body mass 
      matrix. ExtPtfm will give a warning if this is not the case. 
      Internally, ExtPtfm uses this information to determine the mass, 
      rigid-body moments of inertia, and CoG location.

   -  ‘\ ``!stiffness matrix``\ ‘: similar to the mass matrix, the
      values correspond to :math:`\boldsymbol{K}_r`. For a floating 
      structure, there should not be any structural stiffness associated 
      with rigid-body motion; therefore, the first 6 rows and 6 columns of 
      :math:`\boldsymbol{K}_r` should all be zeros.

   -  ‘\ ``!weight constant``\ ‘: followed by some text. The next line 
      should contain :math:`n` float values. The values correspond to the 
      constant self-weight :math:`\boldsymbol{W}_0`. For a floating 
      structure, the 1st, 2nd, and 6th entries of :math:`\boldsymbol{W}_0`
      must be zeros. The constant roll and pitch moments due to self-weight 
      must be consistent with the CoG location derived from the mass matrix.

   -  ‘\ ``!weight stiffness``\ ‘: similar to the mass matrix, the values 
      correspond to :math:`\boldsymbol{K}_W`. For a floating structure, 
      the adopted convention requires 
      :math:`\boldsymbol{W}_0` - :math:`\boldsymbol{K}_W` * (modal displacment) 
      to return the self-weight in a frame of reference following the 
      rigid-body/interface motion. This implies that entry (1,5) of 
      :math:`\boldsymbol{K}_W` must be equal to :math:`\boldsymbol{W}_0(3)`. 
      Entry (2,4) must be equal to :math:`-\boldsymbol{W}_0(3)`. Entries (4,4) 
      and (5,5) must be equal to :math:`\boldsymbol{W}_0(3) * z_{CG}`. 
      Entry (6,4) must be equal to :math:`-\boldsymbol{W}_0(3) * x_{CG}`, and 
      entry (6,5) must be equal to :math:`-\boldsymbol{W}_0(3) * y_{CG}`, where 
      :math:`(x_{CG},y_{CG},z_{CG})` are the location of the rigid-body 
      center of mass measured from the platform reference point defined in 
      ElastoDyn, i.e., (*PtfmRefxt*, *PtfmRefyt*, *PtfmRefzt*). All entries in the 
      1st, 2nd, and 6th columns of :math:`\boldsymbol{K}_W` should be 
      zeros. Again, if the input matrix does not follow this structure, a 
      warning will be given.

.. _epSuperelementInputFile:

Connections input file (provided through ``Conn_Filename``)
-----------------------------------------------------------

The connection input file is used to provide the number and locations of connection 
points on the structure. Currently, the connection points are only used to couple  
with mooring fairleads when any of the available mooring modules is enabled. This file 
also needs to provide the structural motion/deflection at each connection point for 
each Guyan mode and Craig-Bampton mode.

Example
^^^^^^^

A "dummy" example is given below, where numerical values are
implied for ``m``, ``X1``, ``Y1``, ``Z1``, etc.

.. code::

   !Comment
   !Comment
   !nConn: m
   
   !Connections (m)
     X1  Y1  Z1
        [...]
     Xm  Ym  Zm
   
   !Displacement (m)
     [3m x n matrix]

Specifications
^^^^^^^^^^^^^^

The file follows the following specifications:

-  ASCII file

-  The file can start with an arbitrary number of comment lines that 
   start with an exclamation mark ‘\ ``!``\ ‘

-  The following (case insensitive) keywords must be provided next 
   in the order given:

   -  ‘\ ``!nConn``:‘ followed by the integer
      ``m`` giving the number of connection points on the structure

   -  ‘\ ``!Connections``\ ‘: followed by some text. The next :math:`m` 
      lines each containing 3 float values give the :math:`x`, :math:`y`, 
      and :math:`z` coordinates of the :math:`m` connection points. Note 
      that the coordinates are in reference to the same global origin used 
      in, e.g., HydroDyn and MoorDyn input files, and are not relative to 
      the platform reference point defined in ElastoDyn.

   -  ‘\ ``!Displacement``\ ‘: followed by some text. The next :math:`3m` 
      lines each contain :math:`n` float values. The columns of this matrix 
      correspond to the :math:`n` modes defined in the superelement input 
      file and give the :math:`x`,  :math:`y`, and :math:`z` displacements 
      of the first connection point in the first 3 rows, followed by those 
      of the second connection point in the next 3 rows, and so on, for unit 
      motion of each mode. When modeling a floating body with ``RBMod`` > 0, 
      the first 6 modes are the Guyan interface modes that also serve as 
      rigid-body modes. Therefore, the first 6 columns of this matrix should 
      simply describe the linearized rigid-body motion of each connection 
      point following the platform reference point defined in ElastoDyn, 
      i.e., (*PtfmRefxt*, *PtfmRefyt*, *PtfmRefzt*).

   
User modal and connection force file (provided through ``Force_Filename`` and ``FConn_Filename`` )
--------------------------------------------------------------------------------------------------

These two files allow the user to prescribe external forcing to the system modes 
and the connection points. The two files have similar formatting. The only difference 
is in the number of columns.

Example
^^^^^^^

A "dummy" example is given below, where numerical values are
implied for ``k``, ``t1``, ``F1``, etc.

.. code::

   !Comment
   !Comment
   !nSteps: k
   
   !Forcing (N/-)
     t(1)  F_1(1)  F_2(1) ...
           [...]
     t(k)  F_1(k)  F_2(k) ...
   

Specifications
^^^^^^^^^^^^^^

The files follow the following specifications:

-  ASCII file

-  The file can start with an arbitrary number of comment lines that 
   start with an exclamation mark ‘\ ``!``\ ‘

-  The following (case insensitive) keywords must be provided next 
   in the order they are given:

   -  ‘\ ``!NSteps``:‘ followed by the integer
      ``k`` giving the number of time steps in the force time series

   -  ‘\ ``!Forcing``\ ‘: followed by some text. The next :math:`k` 
      lines provide the user-defined modal or connection force time 
      series. The number of columns depends on whether this is for 
      modal forcing or connection forcing. In both cases, the first 
      column is time. For modal forcing, there should be :math:`n` 
      columns following the time column, with each column corresponding 
      to one of the modes defined in the superelement input file. 
      For connection forcing, there should be :math:`3m` columns 
      following the time column, giving the :math:`x`, :math:`y`, and 
      :math:`z` force components at the first connection point, 
      followed by that at the second connection point, and so on. Note 
      that the times in the first column need not match the simulation 
      time step or be evenly spaced. ExtPtfm will perform linear 
      interpolation between the provided times as needed. 
