!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2021 National Renewable Energy Laboratory
! Author: Nicole Mendoza
!
! This file is part of MoorDyn.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
!**********************************************************************************************************************************
MODULE MoorDyn_C

   USE ISO_C_BINDING
   USE MoorDyn
   USE MoorDyn_Types
   USE NWTC_Library
   USE VersionInfo

IMPLICIT NONE

PUBLIC :: MD_C_Init
PUBLIC :: MD_C_UpdateStates
PUBLIC :: MD_C_CalcOutput
PUBLIC :: MD_C_End

!------------------------------------------------------------------------------------
!  Error handling
!     This must exactly match the value in the python-lib. If ErrMsgLen changes at
!     some point in the nwtc-library, this should be updated, but the logic exists
!     to correctly handle different lengths of the strings
integer(IntKi),   parameter            :: ErrMsgLen_C = 1025
integer(IntKi),   parameter            :: IntfStrLen  = 1025       ! length of other strings through the C interface


!------------------------------------------------------------------------------------
!  Version info for display
TYPE(ProgDesc), PARAMETER              :: version   = ProgDesc( 'MoorDyn library', '', '' )


!--------------------------------------------------------------------------------------------------------------------------------------------------------
!  Data storage
!     All MoorDyn data is stored within the following data structures inside this
!     module.  No data is stored within MoorDyn itself, but is instead passed in
!     from this module.  This data is not available to the calling code unless
!     explicitly passed through the interface (derived types such as these are
!     non-trivial to pass through the c-bindings).
TYPE(MD_InitInputType)                  :: InitInp             !< Input data for initialization routine
TYPE(MD_InputType), ALLOCATABLE         :: u(:)                !< An initial guess for the input; input mesh must be defined
TYPE(MD_ParameterType)                  :: p                   !< Parameters
TYPE(MD_ContinuousStateType)            :: x(0:2)              !< Initial continuous states
TYPE(MD_DiscreteStateType)              :: xd(0:2)             !< Initial discrete states
TYPE(MD_ConstraintStateType)            :: z(0:2)              !< Initial constraint states
TYPE(MD_OtherStateType)                 :: other(0:2)          !< Initial other states
TYPE(MD_OutputType)                     :: y                   !< Initial system outputs (outputs are not calculated; only the output mesh is initialized)
TYPE(MD_MiscVarType)                    :: m                   !< Initial misc/optimization variables
TYPE(MD_InitOutputType)                 :: InitOutData         !< Output for initialization routine

!--------------------------------------------------------------------------------------------------------------------------------------------------------
! Time tracking
!     For the solver in MD, previous timesteps input must be stored for extrapolation
!     to the t+dt timestep.  This can be either linear (1) quadratic (2).  The
!     InterpOrder variable tracks what this is and sets the size of the inputs `u`
!     passed into MD. Inputs `u` will be sized as follows:
!        linear    interp     u(2)  with inputs at T,T-dt
!        quadratic interp     u(3)  with inputs at T,T-dt,T-2*dt
!  Time tracking
!     When we are performing a correction step, time information of previous
!     calls is needed to decide how to apply correction logic or cycle the inputs
!     and resave the previous timestep states.
!  Correction steps
!     OpenFAST has the ability to perform correction steps.  During a correction
!     step, new input values are passed in but the timestep remains the same.
!     When this occurs the new input data at time t is used with the state
!     information from the previous timestep (t) to calculate new state values
!     time t+dt in the UpdateStates routine.  In OpenFAST this is all handled by
!     the glue code.  However, here we do not pass state information through the
!     interface and therefore must store it here analogously to how it is handled
!     in the OpenFAST glue code.
INTEGER(IntKi)                          :: InterpOrder         !< Interpolation order: must be 1 (linear) or 2 (quadratic)
REAL(DbKi), DIMENSION(:), ALLOCATABLE   :: InputTimes(:)       !< InputTimes array
REAL(DbKi)                              :: InputTimePrev       !< input time of last UpdateStates call
real(DbKi)                              :: dT_Global           !< dT of the code calling this module
integer(IntKi)                          :: N_Global            !< global timestep
real(DbKi)                              :: T_Initial           !< initial Time of simulation

! We are including the previous state info here (not done in OpenFAST this way)
INTEGER(IntKi),   PARAMETER             :: STATE_LAST = 0      !< Index for previous state (not needed in OF, but necessary here)
INTEGER(IntKi),   PARAMETER             :: STATE_CURR = 1      !< Index for current state
INTEGER(IntKi),   PARAMETER             :: STATE_PRED = 2      !< Index for predicted state

! Note the indexing is different on inputs (no clue why, but thats how OF handles it)
INTEGER(IntKi),   PARAMETER             :: INPUT_LAST = 3      !< Index for previous  input at t-dt
INTEGER(IntKi),   PARAMETER             :: INPUT_CURR = 2      !< Index for current   input at t
INTEGER(IntKi),   PARAMETER             :: INPUT_PRED = 1      !< Index for predicted input at t+dt

!--------------------------------------------------------------------------------------------------------------------------------------------------------
! Meshes for motions and loads
!     Meshes are used within MD to handle all motions and loads. Rather than directly
!     map to those nodes, we will create a mapping to go between the array of node
!     positions passed into this module and what is used inside MD.  This is done
!     through a pair of meshes for the motion and loads corresponding to the node
!     positions passed in.
!------------------------------
!  Meshes for external nodes
!     These point meshes are merely used to simplify the mapping of motions/loads
!     to/from MD using the library mesh mapping routines.  These meshes may contain
!     one or multiple points.
!        - 1 point   -- rigid floating body assumption
!        - N points  -- flexible structure (either floating or fixed bottom)
TYPE(MeshType)                          :: MD_MotionMesh       !< mesh for motions of external nodes
TYPE(MeshType)                          :: MD_LoadMesh         !< mesh for loads  for external nodes
TYPE(MeshMapType)                       :: Map_Motion_2_MD     !< Mesh mapping between input motion mesh and MD
TYPE(MeshMapType)                       :: Map_MD_2_Load       !< Mesh mapping between MD output loads mesh and external nodes mesh

! Motions input (so that we don't have to reallocate all the time)
REAL(ReKi)                              :: tmpPositions(6,1)   !< temp array.  Probably don't need this, but makes conversion from C clearer.
REAL(ReKi)                              :: tmpVelocities(6,1)  !< temp array.  Probably don't need this, but makes conversion from C clearer.
REAL(ReKi)                              :: tmpAccelerations(6,1) !< temp array.  Probably don't need this, but makes conversion from C clearer.
REAL(ReKi)                              :: tmpForces(6,1)      !< temp array.  Probably don't need this, but makes conversion to   C clearer.

CONTAINS

!> This routine sets the error status in C_CHAR for export to calling code.
!! Make absolutely certain that we do not overrun the end of ErrMsg_C.  That is hard coded to 1025,
!! but ErrMsgLen is set in the nwtc_library, and could change without updates here.  We don't want an
!! inadvertant buffer overrun -- that can lead to bad things.
subroutine SetErr(ErrStat, ErrMsg, ErrStat_C, ErrMsg_C)
   integer,                intent(in   )  :: ErrStat                 !< aggregated error message (fortran type)
   character(ErrMsgLen),   intent(in   )  :: ErrMsg                  !< aggregated error message (fortran type)
   integer(c_int),         intent(  out)  :: ErrStat_C
   character(kind=c_char), intent(  out)  :: ErrMsg_C(ErrMsgLen_C)
   integer                                :: i
   ErrStat_C = ErrStat     ! We will send back the same error status that is used in OpenFAST
   if (ErrMsgLen > ErrMsgLen_C-1) then   ! If ErrMsgLen is > the space in ErrMsg_C, do not copy everything over
      ErrMsg_C = TRANSFER( trim(ErrMsg(1:ErrMsgLen_C-1))//C_NULL_CHAR, ErrMsg_C )
   else
      ErrMsg_C = TRANSFER( trim(ErrMsg)//C_NULL_CHAR, ErrMsg_C )
   endif
end subroutine SetErr

!===============================================================================================================
!---------------------------------------------- MD INIT --------------------------------------------------------
!===============================================================================================================
SUBROUTINE MD_C_Init(InputFileString_C, InputFileStringLength_C, DT_C, G_C, RHO_C, DEPTH_C, PtfmInit_C, InterpOrder_C, NumChannels_C, OutputChannelNames_C, OutputChannelUnits_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='MD_C_Init')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: MD_C_Init
!GCC$ ATTRIBUTES DLLEXPORT :: MD_C_Init
#endif
   TYPE(C_PTR)                                    , INTENT(IN   )   :: InputFileString_C        !< Input file as a single string with lines deliniated by C_NULL_CHAR
   INTEGER(C_INT)                                 , INTENT(IN   )   :: InputFileStringLength_C  !< length of the input file string
   REAL(C_DOUBLE)                                 , INTENT(IN   )   :: DT_C
   REAL(C_FLOAT)                                  , INTENT(IN   )   :: G_C
   REAL(C_FLOAT)                                  , INTENT(IN   )   :: RHO_C
   REAL(C_FLOAT)                                  , INTENT(IN   )   :: DEPTH_C
   REAL(C_FLOAT)                                  , INTENT(IN   )   :: PtfmInit_C(6)
   INTEGER(C_INT)                                 , INTENT(IN   )   :: InterpOrder_C
   INTEGER(C_INT)                                 , INTENT(  OUT)   :: NumChannels_C
   CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: OutputChannelNames_C(100000)
   CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: OutputChannelUnits_C(100000)
   INTEGER(C_INT)                                 , INTENT(  OUT)   :: ErrStat_C
   CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: ErrMsg_C(ErrMsgLen_C)

   ! Local Variables
   CHARACTER(KIND=C_char, LEN=InputFileStringLength_C), POINTER     :: InputFileString          !< Input file as a single string with NULL chracter separating lines
   INTEGER(IntKi)                                                   :: ErrStat, ErrStat2
   CHARACTER(ErrMsgLen)                                             :: ErrMsg,  ErrMsg2
   INTEGER                                                          :: I, J, K
   character(*), parameter                                          :: RoutineName = 'MD_C_Init'


   ! Initialize library and display info on this compile
   ErrStat = ErrID_None
   ErrMsg = ''

   CALL NWTC_Init( ProgNameIn=version%Name )
   CALL DispCopyrightLicense( version%Name )
   CALL DispCompileRuntimeInfo( version%Name )



   ! Convert the MD input file to FileInfoType
   !----------------------------------------------------------------------------------------------------------------------------------------------

   ! Get fortran pointer to C_NULL_CHAR deliniated input file as a string 
   CALL C_F_pointer(InputFileString_C, InputFileString)

   ! Convert string inputs to FileInfoType
   CALL InitFileInfo(InputFileString, InitInp%PassedPrimaryInputData, ErrStat2, ErrMsg2); if (Failed()) return

   ! Set other inputs for calling MD_Init
   !----------------------------------------------------------------------------------------------------------------------------------------------
   
   ! Check the interpolation order
   IF (InterpOrder_C .EQ. 1 .OR. InterpOrder_C .EQ. 2) THEN
      InterpOrder = INT(InterpOrder_C, IntKi)
      call AllocAry( InputTimes, InterpOrder+1, 'InputTimes', ErrStat2, ErrMsg2); if (Failed()) return
   ELSE
      ErrStat2 = ErrID_Fatal
      ErrMsg2  = 'InterpOrder must be 1 (linear) or 2 (quadratic)'
      if (Failed()) return
   END IF

   dT_Global                = REAL(DT_C, DbKi)
   N_Global                 = 0_IntKi                     ! Assume we are on timestep 0 at start
   InitInp%FileName         = 'notUsed'
   InitInp%RootName         = 'MDroot'
   InitInp%UsePrimaryInputFile = .FALSE.

   ! Environment variables -- These should be passed in from C.
   InitInp%g                = REAL(G_C, ReKi)
   InitInp%rhoW             = REAL(RHO_C, ReKi)
   InitInp%WtrDepth         = REAL(DEPTH_C, ReKi)

   ! Platform position (x,y,z,Rx,Ry,Rz) -- where rotations are small angle assumption in radians.
   ! This data is used to set the CoupledKinematics mesh that will be used at each timestep call
   CALL AllocAry (InitInp%PtfmInit, 6, 1, 'InitInp%PtfmInit', ErrStat2, ErrMsg2 ); if (Failed()) return
   DO I = 1,6
       InitInp%PtfmInit(I,1)  = REAL(PtfmInit_C(I),ReKi)
   END DO

   ALLOCATE(u(InterpOrder+1), STAT=ErrStat2)
   if (ErrStat2 /= 0) then
      ErrStat2 = ErrID_Fatal
      ErrMsg2  = 'Failed to allocate Inputs type for MD'
      if (Failed()) return
   endif
   
   !-------------------------------------------------
   ! Call the main subroutine MD_Init
   !-------------------------------------------------
   CALL MD_Init(InitInp, u(1), p, x(STATE_CURR), xd(STATE_CURR), z(STATE_CURR), other(STATE_CURR), y, m, dT_Global, InitOutData, ErrStat2, ErrMsg2); if (Failed()) return

   !-------------------------------------------------
   !  Set output channel information for driver code
   !-------------------------------------------------

   ! Number of channels
   NumChannels_C = size(InitOutData%WriteOutputHdr)

   ! Transfer the output channel names and units to c_char arrays for returning
   k=1
   DO i=1,NumChannels_C
      DO j=1,ChanLen    ! max length of channel name.  Same for units
         OutputChannelNames_C(k)=InitOutData%WriteOutputHdr(i)(j:j)
         OutputChannelUnits_C(k)=InitOutData%WriteOutputUnt(i)(j:j)
         k=k+1
      END DO
   END DO

   ! Null terminate the string
   OutputChannelNames_C(k) = C_NULL_CHAR
   OutputChannelUnits_C(k) = C_NULL_CHAR

   !-------------------------------------------------------------
   ! Set the interface  meshes for motion inputs and loads output
   !-------------------------------------------------------------
   DO i = 1,6
       tmpPositions(i,1)     = REAL(PtfmInit_C(i),ReKi)
   END DO
   tmpVelocities     = 0_ReKi
   tmpAccelerations  = 0_ReKi
   CALL SetMotionLoadsInterfaceMeshes(ErrStat2,ErrMsg2); if (Failed()) return

   DO i=2,InterpOrder+1
      CALL MD_CopyInput (u(1),  u(i),  MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
   END DO
   InputTimePrev = -dT_Global    ! Initialize for MD_C_UpdateStates

   !-------------------------------------------------------------
   ! Initial setup of other pieces of x,xd,z,other
   !-------------------------------------------------------------
   CALL MD_CopyContState  ( x(    STATE_CURR), x(    STATE_PRED), MESH_NEWCOPY, ErrStat2, ErrMsg2);   if (Failed())  return
   CALL MD_CopyDiscState  ( xd(   STATE_CURR), xd(   STATE_PRED), MESH_NEWCOPY, ErrStat2, ErrMsg2);   if (Failed())  return
   CALL MD_CopyConstrState( z(    STATE_CURR), z(    STATE_PRED), MESH_NEWCOPY, ErrStat2, ErrMsg2);   if (Failed())  return
   CALL MD_CopyOtherState ( other(STATE_CURR), other(STATE_PRED), MESH_NEWCOPY, ErrStat2, ErrMsg2);   if (Failed())  return

   !-------------------------------------------------------------
   ! Setup the previous timestep copies of states
   !-------------------------------------------------------------
   CALL MD_CopyContState  ( x(    STATE_CURR), x(    STATE_LAST), MESH_NEWCOPY, ErrStat2, ErrMsg2);   if (Failed())  return
   CALL MD_CopyDiscState  ( xd(   STATE_CURR), xd(   STATE_LAST), MESH_NEWCOPY, ErrStat2, ErrMsg2);   if (Failed())  return
   CALL MD_CopyConstrState( z(    STATE_CURR), z(    STATE_LAST), MESH_NEWCOPY, ErrStat2, ErrMsg2);   if (Failed())  return
   CALL MD_CopyOtherState ( other(STATE_CURR), other(STATE_LAST), MESH_NEWCOPY, ErrStat2, ErrMsg2);   if (Failed())  return

   !-------------------------------------------------
   ! Clean up variables and set up for MD_C_CalcOutput
   !------------------------------------------------- 
   CALL MD_DestroyInitInput( InitInp, ErrStat2, ErrMsg2 );        if (Failed())  return
   CALL MD_DestroyInitOutput( InitOutData, ErrStat2, ErrMsg2 );   if (Failed())  return

   call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)


CONTAINS
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed) then
         call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
      endif
   end function Failed
END SUBROUTINE MD_C_Init

!===============================================================================================================
!---------------------------------------------- MD UPDATE STATES -----------------------------------------------
!===============================================================================================================
SUBROUTINE MD_C_UpdateStates(Time_C, TimeNext_C, POSITIONS_C, VELOCITIES_C, ACCELERATIONS_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='MD_C_UpdateStates')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: MD_C_UpdateStates
!GCC$ ATTRIBUTES DLLEXPORT :: MD_C_UpdateStates
#endif
   real(c_double),                                  intent(in   )  :: Time_C
   real(c_double),                                  intent(in   )  :: TimeNext_C
   REAL(C_FLOAT)                                  , INTENT(IN   )   :: POSITIONS_C(1,6)
   REAL(C_FLOAT)                                  , INTENT(IN   )   :: VELOCITIES_C(1,6)
   REAL(C_FLOAT)                                  , INTENT(IN   )   :: ACCELERATIONS_C(1,6)
   INTEGER(C_INT)                                 , INTENT(  OUT)   :: ErrStat_C
   CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: ErrMsg_C(ErrMsgLen_C)

   ! Local Variables
   INTEGER(IntKi)                                                   :: ErrStat, ErrStat2, J
   CHARACTER(ErrMsgLen)                                             :: ErrMsg,  ErrMsg2
   LOGICAL                                                          :: CorrectionStep
   character(*), parameter                                          :: RoutineName = 'MD_C_UpdateStates'

   ! Set up error handling for MD_C_CalcOutput
   ErrStat = ErrID_None
   ErrMsg = ''
   CorrectionStep = .FALSE.

   !-------------------------------------------------------
   ! Check the time for current timestep and next timestep
   !-------------------------------------------------------
   !     These inputs are used in the time stepping algorithm within MD_UpdateStates
   !     For quadratic interpolation (InterpOrder==2), 3 timesteps are used.  For
   !     linear (InterOrder==1), 2 timesteps (the MD code can handle either).
   !        u(1)  inputs at t + dt        ! Next timestep
   !        u(2)  inputs at t             ! This timestep
   !        u(3)  inputs at t - dt        ! previous timestep (quadratic only)
   !
   !  NOTE: Within MD, the Radiation calculations can be done at an integer multiple of the
   !        timestep.  This is checked at each UpdateStates call.  However, if we compile
   !        in double precision, the values of Time_C and TimeNext_C are in double precison,
   !        but InputTimes is in DbKi (which is promoted quad precision when compiling in
   !        double precision) and the check may fail.  So we are going to set the times we
   !        we pass over to UpdateStates using the global timestep and the stored DbKi value
   !        for the timestep rather than the lower precision (when compiled double) time
   !        values passed in.  It is a bit of a clumsy workaround for this precision loss,
   !        but should not affect any results.

   !  Check if we are repeating an UpdateStates call (for example in a predictor/corrector loop)
   if ( EqualRealNos( real(Time_C,DbKi), InputTimePrev ) ) then
      CorrectionStep = .true.
   else ! Setup time input times array
      InputTimePrev          = real(Time_C,DbKi)            ! Store for check next time
      if (InterpOrder>1) then ! quadratic, so keep the old time
         InputTimes(INPUT_LAST) = ( N_Global - 1 ) * dT_Global    ! u(3) at T-dT
      endif
      InputTimes(INPUT_CURR) =   N_Global       * dT_Global       ! u(2) at T
      InputTimes(INPUT_PRED) = ( N_Global + 1 ) * dT_Global       ! u(1) at T+dT
      N_Global = N_Global + 1_IntKi                               ! increment counter to T+dT
   endif


   IF (CorrectionStep) THEN
       ! Step back to previous state because we are doing a correction step
       !     -- repeating the T -> T+dt update with new inputs at T+dt
       !     -- the STATE_CURR contains states at T+dt from the previous call, so revert those
       CALL MD_CopyContState   (x(    STATE_LAST), x(    STATE_CURR), MESH_UPDATECOPY, ErrStat2, ErrMsg2);  IF (Failed())  RETURN
       CALL MD_CopyDiscState   (xd(   STATE_LAST), xd(   STATE_CURR), MESH_UPDATECOPY, ErrStat2, ErrMsg2);  IF (Failed())  RETURN
       CALL MD_CopyConstrState (z(    STATE_LAST), z(    STATE_CURR), MESH_UPDATECOPY, ErrStat2, ErrMsg2);  IF (Failed())  RETURN
       CALL MD_CopyOtherState  (other(STATE_LAST), other(STATE_CURR), MESH_UPDATECOPY, ErrStat2, ErrMsg2);  IF (Failed())  RETURN
   ELSE
       ! Cycle inputs back one timestep since we are moving forward in time.
       IF (InterpOrder>1) THEN ! quadratic, so keep the old time
           CALL MD_CopyInput( u(INPUT_CURR), u(INPUT_LAST), MESH_UPDATECOPY, ErrStat2, ErrMsg2);   IF (Failed())  RETURN
       END IF
       ! Move inputs from previous t+dt (now t) to t
       CALL MD_CopyInput( u(INPUT_PRED), u(INPUT_CURR), MESH_UPDATECOPY, ErrStat2, ErrMsg2);       IF (Failed())  RETURN
   END IF

   ! Reshape position and velocity (transposing from a row vector to a column vector)
   DO J = 1,6
       tmpPositions(J,1)     = REAL(POSITIONS_C(1,J),ReKi)
       tmpVelocities(J,1)    = REAL(VELOCITIES_C(1,J),ReKi)
       tmpAccelerations(J,1) = REAL(ACCELERATIONS_C(1,J),ReKi)
   END DO

   ! Transfer motions to input meshes
   CALL Set_MotionMesh( ErrStat2, ErrMsg2 );          IF (Failed())  RETURN
   CALL MD_SetInputMotion( u(INPUT_PRED), ErrStat2, ErrMsg2 ); IF (Failed())  RETURN

   ! Set copy the current state over to the predicted state for sending to UpdateStates
   !     -- The STATE_PREDicted will get updated in the call.
   !     -- The UpdateStates routine expects this to contain states at T at the start of the call (history not passed in)
   CALL MD_CopyContState   (x(    STATE_CURR), x(    STATE_PRED), MESH_UPDATECOPY, ErrStat2, ErrMsg2);  IF (Failed())  RETURN
   CALL MD_CopyDiscState   (xd(   STATE_CURR), xd(   STATE_PRED), MESH_UPDATECOPY, ErrStat2, ErrMsg2);  IF (Failed())  RETURN
   CALL MD_CopyConstrState (z(    STATE_CURR), z(    STATE_PRED), MESH_UPDATECOPY, ErrStat2, ErrMsg2);  IF (Failed())  RETURN
   CALL MD_CopyOtherState  (other(STATE_CURR), other(STATE_PRED), MESH_UPDATECOPY, ErrStat2, ErrMsg2);  IF (Failed())  RETURN

   !-------------------------------------------------
   ! Call the main subroutine MD_UpdateStates
   !-------------------------------------------------
   CALL MD_UpdateStates( InputTimes(INPUT_CURR), N_Global, u, InputTimes, p, x(STATE_PRED), xd(STATE_PRED), z(STATE_PRED), other(STATE_PRED), m, ErrStat2, ErrMsg2);  IF (Failed())  RETURN

   !-------------------------------------------------------
   ! Cycle the states
   !-------------------------------------------------------
   ! Move current state at T to previous state at T-dt
   !     -- STATE_LAST now contains info at time T
   !     -- this allows repeating the T --> T+dt update
   CALL MD_CopyContState   (x(    STATE_CURR), x(    STATE_LAST), MESH_UPDATECOPY, ErrStat2, ErrMsg2);  IF (Failed())  RETURN
   CALL MD_CopyDiscState   (xd(   STATE_CURR), xd(   STATE_LAST), MESH_UPDATECOPY, ErrStat2, ErrMsg2);  IF (Failed())  RETURN
   CALL MD_CopyConstrState (z(    STATE_CURR), z(    STATE_LAST), MESH_UPDATECOPY, ErrStat2, ErrMsg2);  IF (Failed())  RETURN
   CALL MD_CopyOtherState  (other(STATE_CURR), other(STATE_LAST), MESH_UPDATECOPY, ErrStat2, ErrMsg2);  IF (Failed())  RETURN
   ! Update the predicted state as the new current state
   !     -- we have now advanced from T to T+dt.  This allows calling with CalcOuput to get the outputs at T+dt
   CALL MD_CopyContState   (x(    STATE_PRED), x(    STATE_CURR), MESH_UPDATECOPY, ErrStat2, ErrMsg2);  IF (Failed())  RETURN
   CALL MD_CopyDiscState   (xd(   STATE_PRED), xd(   STATE_CURR), MESH_UPDATECOPY, ErrStat2, ErrMsg2);  IF (Failed())  RETURN
   CALL MD_CopyConstrState (z(    STATE_PRED), z(    STATE_CURR), MESH_UPDATECOPY, ErrStat2, ErrMsg2);  IF (Failed())  RETURN
   CALL MD_CopyOtherState  (other(STATE_PRED), other(STATE_CURR), MESH_UPDATECOPY, ErrStat2, ErrMsg2);  IF (Failed())  RETURN



   call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)

CONTAINS
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)  call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
   end function Failed
END SUBROUTINE MD_C_UpdateStates

!===============================================================================================================
!---------------------------------------------- MD CALC OUTPUT -------------------------------------------------
!===============================================================================================================
!> Calculate the moordyn results iven the current set of states and inputs
SUBROUTINE MD_C_CalcOutput(Time_C, POSITIONS_C, VELOCITIES_C, ACCELERATIONS_C, FORCES_C, OUTPUTS_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='MD_C_CalcOutput')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: MD_C_CalcOutput
!GCC$ ATTRIBUTES DLLEXPORT :: MD_C_CalcOutput
#endif
   REAL(C_DOUBLE)                                 , INTENT(IN   )   :: Time_C
   REAL(C_FLOAT)                                  , INTENT(IN   )   :: POSITIONS_C(1,6)
   REAL(C_FLOAT)                                  , INTENT(IN   )   :: VELOCITIES_C(1,6)
   REAL(C_FLOAT)                                  , INTENT(IN   )   :: ACCELERATIONS_C(1,6)
   REAL(C_FLOAT)                                  , INTENT(  OUT)   :: FORCES_C(1,6)
   REAL(C_FLOAT)                                  , INTENT(  OUT)   :: OUTPUTS_C(p%NumOuts)
   INTEGER(C_INT)                                 , INTENT(  OUT)   :: ErrStat_C
   CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: ErrMsg_C(ErrMsgLen_C)

   ! Local Variables
   REAL(DbKi)                                                       :: t
   INTEGER(IntKi)                                                   :: ErrStat, ErrStat2, J
   CHARACTER(ErrMsgLen)                                             :: ErrMsg,  ErrMsg2
   character(*), parameter                                          :: RoutineName = 'MD_C_CalcOutput'

   ! Set up error handling for MD_C_CalcOutput
   ErrStat = ErrID_None
   ErrMsg = ''

   ! Set up inputs to MD_CalcOutput
   !-----------------------------------------------------------------------------------------------------------

   ! Time
   t = REAL(Time_C, DbKi)

   ! Reshape position and velocity (from row vector to a column vector)
   DO J = 1,6
       tmpPositions(J,1)       = REAL(POSITIONS_C(1,J),ReKi)
       tmpVelocities(J,1)      = REAL(VELOCITIES_C(1,J),ReKi)
       tmpAccelerations(J,1)   = REAL(ACCELERATIONS_C(1,J),ReKi)
   END DO

   ! Transfer motions to input meshes
   CALL Set_MotionMesh(ErrStat2, ErrMsg2 );  if (Failed()) return;

   ! transfer input motion mesh to u(1) meshes
   CALL MD_SetInputMotion( u(1), ErrStat2, ErrMsg2 );  if (Failed()) return;

   !-------------------------------------------------
   ! Call the main subroutine MD_CalcOutput
   !-------------------------------------------------
   CALL MD_CalcOutput( t, u(1), p, x(STATE_CURR), xd(STATE_CURR), z(STATE_CURR), other(STATE_CURR), y, m, ErrStat2, ErrMsg2 );  if (Failed()) return;

   !-------------------------------------------------
   ! Convert the outputs of MD_calcOutput back to C
   !-------------------------------------------------
   ! Transfer resulting load meshes to intermediate mesh
   CALL MD_TransferLoads( u(1), y, ErrStat2, ErrMsg2 );  if (Failed()) return;

   ! Set output force/moment array
   CALL Set_OutputLoadArray( )

   ! Reshape for return
   DO J = 1,6
       FORCES_C(1,J) = REAL(tmpForces(J,1), c_float)
   END DO

   OUTPUTS_C = REAL(y%WriteOutput, C_FLOAT)

CONTAINS
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed) then
         call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
      endif
   end function Failed
END SUBROUTINE MD_C_CalcOutput

!===============================================================================================================
!----------------------------------------------- MD END --------------------------------------------------------
!===============================================================================================================
!> Cleanup memory
!! NOTE: the error handling here is slightly different than in other routines
SUBROUTINE MD_C_End(ErrStat_C,ErrMsg_C) BIND (C, NAME='MD_C_End')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: MD_C_End
!GCC$ ATTRIBUTES DLLEXPORT :: MD_C_End
#endif
   INTEGER(C_INT)                , INTENT(  OUT)      :: ErrStat_C
   CHARACTER(KIND=C_CHAR)        , INTENT(  OUT)      :: ErrMsg_C(ErrMsgLen_C)

   ! Local variables
   INTEGER(IntKi)                                     :: ErrStat, ErrStat2, i
   CHARACTER(ErrMsgLen)                               :: ErrMsg,  ErrMsg2
   character(*), parameter                            :: RoutineName = 'MD_C_End'

   ! Set up error handling for MD_C_End
   ErrStat = ErrID_None
   ErrMsg = ''

   ! Call the main subroutine MD_End
   CALL MD_End(u(1), p, x(1), xd(1), z(1), other(1), y, m, ErrStat2, ErrMsg2)
   call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   !  NOTE: MoorDyn_End only takes 1 instance of u, not the array.  So extra
   !        logic is required here (this isn't necessary in the fortran driver
   !        or in openfast, but may be when this code is called from C, Python,
   !        or some other code using the c-bindings)
   IF (allocated(u)) THEN
      DO i=2,size(u)
         CALL MD_DestroyInput( u(i), ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
      IF (allocated(u))             deallocate(u)
   END IF

   ! Destroy any other copies of states (rerun on (STATE_CURR) is ok)
   call MD_DestroyContState(   x(    STATE_LAST), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call MD_DestroyContState(   x(    STATE_CURR), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call MD_DestroyContState(   x(    STATE_PRED), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call MD_DestroyDiscState(   xd(   STATE_LAST), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call MD_DestroyDiscState(   xd(   STATE_CURR), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call MD_DestroyDiscState(   xd(   STATE_PRED), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call MD_DestroyConstrState( z(    STATE_LAST), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call MD_DestroyConstrState( z(    STATE_CURR), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call MD_DestroyConstrState( z(    STATE_PRED), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call MD_DestroyOtherState(  other(STATE_LAST), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call MD_DestroyOtherState(  other(STATE_CURR), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call MD_DestroyOtherState(  other(STATE_PRED), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   ! if deallocate other items now
   if (allocated(InputTimes))    deallocate(InputTimes)

   ! Clear out mesh related data storage
   call ClearMesh()

   call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
CONTAINS
   !> Don't leave junk in memory.  So destroy meshes and mappings.
   subroutine ClearMesh()
      call MeshDestroy( MD_MotionMesh, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      call MeshDestroy( MD_LoadMesh, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      ! Destroy mesh mappings
      call NWTC_Library_Destroymeshmaptype( Map_Motion_2_MD, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      call NWTC_Library_Destroymeshmaptype( Map_MD_2_Load, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   end subroutine ClearMesh
END SUBROUTINE MD_C_End

!===============================================================================================================
!----------------------------------------- ADDITIONAL SUBROUTINES ----------------------------------------------
!===============================================================================================================
!! This subroutine sets the interface meshes to map to the input motions to the MD
!! meshes.  This subroutine also sets the meshes for the output loads.
SUBROUTINE SetMotionLoadsInterfaceMeshes(ErrStat,ErrMsg)
   INTEGER(IntKi),         INTENT(  OUT)  :: ErrStat    !< temporary error status
   CHARACTER(ErrMsgLen),   INTENT(  OUT)  :: ErrMsg     !< temporary error message
   REAL(ReKi)                             :: InitPos(3)
   REAL(R8Ki)                             :: theta(3)
   REAL(R8Ki)                             :: Orient(3,3)
   !-------------------------------------------------------------
   ! Set the interface  meshes for motion inputs and loads output
   !-------------------------------------------------------------
   ! Motion mesh
   !     This point mesh may contain more than one point. Mapping will be used to map
   !     this to the input mesh.
   CALL MeshCreate(  MD_MotionMesh                       ,  &
                     IOS              = COMPONENT_INPUT  ,  &
                     Nnodes           = 1                ,  &
                     ErrStat          = ErrStat          ,  &
                     ErrMess          = ErrMsg           ,  &
                     TranslationDisp  = .TRUE.,    Orientation = .TRUE., &
                     TranslationVel   = .TRUE.,    RotationVel = .TRUE., &
                     TranslationAcc   = .TRUE.,    RotationAcc = .TRUE.  )
   IF (ErrStat >= AbortErrLev) RETURN

   ! initial position and orientation of node
   InitPos  = tmpPositions(1:3,1)
   theta    = REAL(tmpPositions(4:6,1),DbKi)    ! convert ReKi to DbKi to avoid roundoff
   CALL SmllRotTrans( 'InputRotation', theta(1), theta(2), theta(3), Orient, 'Orient', ErrStat, ErrMsg )
   CALL MeshPositionNode(  MD_MotionMesh            , &
                           1                        , &
                           InitPos                  , &  ! position
                           ErrStat, ErrMsg          , &
                           Orient                     )  ! orientation
   IF (ErrStat >= AbortErrLev) RETURN
    
   CALL MeshConstructElement ( MD_MotionMesh, ELEMENT_POINT, ErrStat, ErrMsg, 1 )
   IF (ErrStat >= AbortErrLev) RETURN

   CALL MeshCommit ( MD_MotionMesh, ErrStat, ErrMsg )
   IF (ErrStat >= AbortErrLev) RETURN

   MD_MotionMesh%RemapFlag  = .TRUE.

   ! For checking the mesh, uncomment this.
   !     note: CU is is output unit (platform dependent).
   !call MeshPrintInfo( CU, MD_MotionMesh )
 
   !-------------------------------------------------------------
   ! Loads mesh
   !     This point mesh may contain more than one point. Mapping will be used to map
   !     the loads from output mesh
   CALL MeshCopy( SrcMesh  = MD_MotionMesh      ,&
                  DestMesh = MD_LoadMesh        ,&
                  CtrlCode = MESH_SIBLING       ,&
                  IOS      = COMPONENT_OUTPUT   ,&
                  ErrStat  = ErrStat            ,&
                  ErrMess  = ErrMsg             ,&
                  Force    = .TRUE.             ,&
                  Moment   = .TRUE.             )
       IF (ErrStat >= AbortErrLev) RETURN
   
   MD_LoadMesh%RemapFlag  = .TRUE.

   ! For checking the mesh, uncomment this.
   !     note: CU is is output unit (platform dependent).
   !call MeshPrintInfo( CU, MD_LoadMesh )

   !-------------------------------------------------------------
   ! Set the mapping meshes
   IF ( allocated(u(1)%CoupledKinematics) ) THEN      ! input motions
       CALL MeshMapCreate( MD_MotionMesh, u(1)%CoupledKinematics(1), Map_Motion_2_MD, ErrStat, ErrMsg )
       IF (ErrStat >= AbortErrLev) RETURN
   END IF
   IF ( allocated(y%CoupledLoads) ) THEN              ! output loads
       CALL MeshMapCreate( y%CoupledLoads(1), MD_LoadMesh, Map_MD_2_Load, ErrStat, ErrMsg )
       IF (ErrStat >= AbortErrLev) RETURN
   END IF
END SUBROUTINE SetMotionLoadsInterfaceMeshes

!---------------------------------------------------------------------------------------------------------------
!> This routine is operating on module level data, hence few inputs
SUBROUTINE Set_MotionMesh(ErrStat, ErrMsg)
   REAL(R8Ki)                                :: theta(3)
   REAL(R8Ki)                                :: Orient(3,3)
   INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat
   CHARACTER(ErrMsgLen),      INTENT(  OUT)  :: ErrMsg
   ! Set mesh corresponding to input motions
   theta = REAL(tmpPositions(4:6,1),DbKi)    ! convert ReKi to DbKi to avoid roundoff
   CALL SmllRotTrans( 'InputRotation', theta(1), theta(2), theta(3), Orient, 'Orient', ErrStat, ErrMsg )
   MD_MotionMesh%TranslationDisp(1:3,1) = tmpPositions(1:3,1) - MD_MotionMesh%Position(1:3,1)  ! relative displacement only
   MD_MotionMesh%Orientation(1:3,1:3,1) = Orient
   MD_MotionMesh%TranslationVel( 1:3,1) = tmpVelocities(1:3,1)
   MD_MotionMesh%RotationVel(    1:3,1) = tmpVelocities(4:6,1)
   MD_MotionMesh%TranslationAcc( 1:3,1) = tmpAccelerations(1:3,1)
   MD_MotionMesh%RotationAcc(    1:3,1) = tmpAccelerations(4:6,1)
END SUBROUTINE Set_MotionMesh

!---------------------------------------------------------------------------------------------------------------
!> Map the motion of the intermediate input mesh over to the input meshes
!! This routine is operating on module level data, hence few inputs
SUBROUTINE MD_SetInputMotion( u_local, ErrStat, ErrMsg )
   TYPE(MD_InputType),        INTENT(INOUT)  :: u_local
   INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat
   CHARACTER(ErrMsgLen),      INTENT(  OUT)  :: ErrMsg

   ErrStat = ErrID_None
   ErrMsg = ''

   IF ( allocated(u_local%CoupledKinematics) ) THEN
      CALL Transfer_Point_to_Point( MD_MotionMesh, u_local%CoupledKinematics(1), Map_Motion_2_MD, ErrStat, ErrMsg )
   END IF
END SUBROUTINE MD_SetInputMotion

!---------------------------------------------------------------------------------------------------------------
!> Map the loads of the output meshes to the intermediate output mesh.  Since
!! we are mapping two meshes over to a single one, we use an intermediate
!! temporary mesh -- This step is not currently necessary in MD since only one set
!! of load points is used with MD.
SUBROUTINE MD_TransferLoads( u_local, y_local, ErrStat, ErrMsg )
   TYPE(MD_InputType),        INTENT(IN   )  :: u_local     ! Only one input (probably at T)
   TYPE(MD_OutputType),       INTENT(IN   )  :: y_local     ! Only one input (probably at T)
   INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat
   CHARACTER(ErrMsgLen),      INTENT(  OUT)  :: ErrMsg

   ErrStat = ErrID_None
   ErrMsg = ''
   MD_LoadMesh%Force            = 0.0_ReKi
   MD_LoadMesh%Moment           = 0.0_ReKi

   !  mesh
   IF ( allocated(y_local%CoupledLoads) ) THEN
       CALL Transfer_Point_to_Point( y_local%CoupledLoads(1), MD_LoadMesh, Map_MD_2_Load, ErrStat, ErrMsg, u_local%CoupledKinematics(1), MD_MotionMesh )
       IF (ErrStat >= AbortErrLev)  RETURN
   END IF
END SUBROUTINE MD_TransferLoads

!---------------------------------------------------------------------------------------------------------------
!> Transfer the loads from the load mesh to the temporary array for output
!! This routine is operating on module level data, hence few inputs
SUBROUTINE Set_OutputLoadArray()
   ! Set mesh corresponding to input motions
   tmpForces(1:3,1)          = MD_LoadMesh%Force (1:3,1)
   tmpForces(4:6,1)          = MD_LoadMesh%Moment(1:3,1)
END SUBROUTINE Set_OutputLoadArray

END MODULE
