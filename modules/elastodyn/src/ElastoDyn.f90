!**********************************************************************************************************************************
! The ElastoDyn.f90, ElastoDyn_IO.f90, and ElastoDyn_Types.f90 make up the ElastoDyn module of the
! FAST Modularization Framework. ElastoDyn_Types is auto-generated based on FAST_Registry.txt.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2012-2016  National Renewable Energy Laboratory
! Copyright (C) 2017-2018  Envision Energy USA, LTD.
!
!    This file is part of ElastoDyn.
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
MODULE ElastoDyn

   USE ElastoDyn_IO
   USE NWTC_LAPACK

   USE ED_UserSubs         ! <- module not in the FAST Framework!

   USE ElastoDyn_AllBldNdOuts_IO

   IMPLICIT NONE

   PRIVATE
   
   REAL(R8Ki)  :: MinPerturb = SQRT( EPSILON( 1.0_R8Ki ) ) ! minimum value for perturbation in ED jacobians


      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: ED_Init                           ! Initialization routine
   PUBLIC :: ED_End                            ! Ending routine (includes clean up)

   PUBLIC :: ED_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                               !   continuous states, and updating discrete states
   PUBLIC :: ED_CalcOutput                     ! Routine for computing outputs

   PUBLIC :: ED_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: ED_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: ED_UpdateDiscState                ! Tight coupling routine for updating discrete states

   PUBLIC :: ED_JacobianPInput                 ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
                                               !   (Xd), and constraint-state (Z) equations all with respect to the inputs (u)
   PUBLIC :: ED_JacobianPContState             ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
                                               !   (Xd), and constraint-state (Z) equations all with respect to the continuous
                                               !   states (x)
   PUBLIC :: ED_JacobianPDiscState             ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
                                               !   (Xd), and constraint-state (Z) equations all with respect to the discrete
                                               !   states (xd)
   PUBLIC :: ED_JacobianPConstrState           ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
                                               !   (Xd), and constraint-state (Z) equations all with respect to the constraint
                                               !   states (z)

   PUBLIC :: ED_GetOP                          ! Routine to pack the operating point values (for linearization) into arrays
   
CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the start of the simulation to perform initialization steps.
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
SUBROUTINE ED_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(ED_InitInputType),       INTENT(IN   )  :: InitInp     !< Input data for initialization routine
   TYPE(ED_InputType),           INTENT(  OUT)  :: u           !< An initial guess for the input; input mesh must be defined
   TYPE(ED_ParameterType),       INTENT(  OUT)  :: p           !< Parameters
   TYPE(ED_ContinuousStateType), INTENT(  OUT)  :: x           !< Initial continuous states
   TYPE(ED_DiscreteStateType),   INTENT(  OUT)  :: xd          !< Initial discrete states
   TYPE(ED_ConstraintStateType), INTENT(  OUT)  :: z           !< Initial guess of the constraint states
   TYPE(ED_OtherStateType),      INTENT(  OUT)  :: OtherState  !< Initial other states
   TYPE(ED_OutputType),          INTENT(  OUT)  :: y           !< Initial system outputs (outputs are not calculated;
                                                               !!   only the output mesh is initialized)
   TYPE(ED_MiscVarType),         INTENT(  OUT)  :: m           !< Misc variables for optimization (not copied in glue code)
   REAL(DbKi),                   INTENT(INOUT)  :: Interval    !< Coupling interval in seconds: the rate that
                                                               !!   (1) ED_UpdateStates() is called in loose coupling &
                                                               !!   (2) ED_UpdateDiscState() is called in tight coupling.
                                                               !!   Input is the suggested time from the glue code;
                                                               !!   Output is the actual coupling interval that will be used
                                                               !!   by the glue code.
   TYPE(ED_InitOutputType),      INTENT(  OUT)  :: InitOut     !< Output for initialization routine
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


      ! Local variables

   TYPE(ED_InputFile)                           :: InputFileData           ! Data stored in the module's input file
   INTEGER(IntKi)                               :: ErrStat2                ! temporary Error status of the operation
   INTEGER(IntKi)                               :: i, K                    ! loop counters
   CHARACTER(ErrMsgLen)                         :: ErrMsg2                 ! temporary Error message if ErrStat /= ErrID_None
   REAL(R8Ki)                                   :: TransMat(3,3)            ! Initial rotation matrix at Platform Refz


      ! Initialize variables for this routine

   ErrStat = ErrID_None
   ErrMsg  = ""


      ! Initialize the NWTC Subroutine Library

   CALL NWTC_Init( EchoLibVer=.FALSE. )

      ! Display the module information

   CALL DispNVD( ED_Ver )

      !............................................................................................
      ! Read the input file and validate the data
      !............................................................................................
   p%BD4Blades = .NOT. InitInp%CompElast           ! if we're not using ElastoDyn for the blades, use BeamDyn
   p%UseAD14   = LEN_TRIM(InitInp%ADInputFile) > 0 ! if we're using AD14, we need to use the AD14 input files

   p%RootName = InitInp%RootName ! FAST already adds '.ED' to the root name

   p%Gravity = InitInp%Gravity
   
   CALL ED_ReadInput( InitInp%InputFile, InitInp%ADInputFile, InputFileData, p%BD4Blades, Interval, p%RootName, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   IF ( p%BD4Blades ) THEN
   
         ! Set DOFs to FALSE for whatever values you don't want on for BeamDyn
      InputFileData%FlapDOF1 = .FALSE.
      InputFileData%FlapDOF2 = .FALSE.
      InputFileData%EdgeDOF  = .FALSE.
      
         ! Set other values not used for BeamDyn      
      InputFileData%OoPDefl  = 0.0_ReKi
      InputFileData%IPDefl   = 0.0_ReKi
      InputFileData%TipMass  = 0.0_ReKi
      InputFileData%TipRad   = 0.0_ReKi
      InputFileData%NBlGages = 0
      InputFileData%BldGagNd = 0
      InputFileData%BldNodes = 0
       
   END IF


   CALL ED_ValidateInput( InputFileData, p%BD4Blades, InitInp%Linearize, InitInp%MHK, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN

      !............................................................................................
      ! Define parameters here:
      !............................................................................................
   CALL ED_SetParameters( InitInp, InputFileData, p, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN

      !............................................................................................
      ! Define initial system states here:
      !............................................................................................
   xd%DummyDiscState          = 0.                                             ! we don't have discrete states
   z%DummyConstrState         = 0.                                             ! we don't have constraint states

      ! initialize the continuous states:
   CALL Init_ContStates( x, p, InputFileData, OtherState, ErrStat2, ErrMsg2 )     ! initialize the continuous states
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN
   
   
      ! Initialize other states and misc variables:
   CALL Init_MiscOtherStates( m, OtherState, p, x, InputFileData, ErrStat2, ErrMsg2 )    ! initialize the other states (must do this after ED_SetParameters)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN
      

      !............................................................................................
      ! Define initial guess for the system inputs here:
      !............................................................................................

         ! allocate all the arrays that store data in the input type and initialize the values:
   CALL Init_u( u, p, x, InputFileData, m, ErrStat2, ErrMsg2 )      
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN                  
   

      !............................................................................................
      ! Define system output initializations (set up meshes) here:
      !............................................................................................

   CALL ED_AllocOutput(p, m, u, y, ErrStat2, ErrMsg2) ! u is sent so we can create sibling meshes
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN
      
   IF (p%BD4Blades) THEN
         ! we need the initial outputs to send BeamDyn for its states. I should probably put them in
         ! the InitOutput type, but I'm going to use the ED_Outputs instead
      ! 
      !   ! set the coordinate system variables:
      !CALL SetCoordSy( t, m%CoordSys, m%RtHS, u%BlPitchCom, p, x, ErrStat, ErrMsg )
      !   IF (ErrStat >= AbortErrLev) RETURN
      !
      !CALL CalculatePositions(        p, x, m%CoordSys,    m%RtHS ) ! calculate positions
      !CALL CalculateAngularPosVelPAcc(p, x, m%CoordSys,    m%RtHS ) ! calculate angular positions, velocities, and partial accelerations, including partial angular quantities
      !CALL CalculateLinearVelPAcc(    p, x, m%CoordSys,    m%RtHS ) ! calculate linear velocities and partial accelerations

      
      
      CALL ED_CalcOutput( 0.0_DbKi, u, p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2 )
         CALL CheckError( ErrStat2, ErrMsg2 )
   END IF
      
      
      !............................................................................................
      ! Define initialization-routine output here:
      !............................................................................................
   CALL AllocAry( InitOut%WriteOutputHdr, p%numOuts + p%BldNd_TotNumOuts, 'WriteOutputHdr', errStat2, errMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN

   CALL AllocAry( InitOut%WriteOutputUnt, p%numOuts + p%BldNd_TotNumOuts, 'WriteOutputUnt', errStat2, errMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN

   do i=1,p%NumOuts
      InitOut%WriteOutputHdr(i) = p%OutParam(i)%Name
      InitOut%WriteOutputUnt(i) = p%OutParam(i)%Units
   end do

      ! Set the info in WriteOutputHdr and WriteOutputUnt
   CALL AllBldNdOuts_InitOut( InitOut, p, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN
      
   InitOut%Ver         = ED_Ver
   InitOut%NumBl       = p%NumBl
   InitOut%BladeLength = p%TipRad - p%HubRad
   InitOut%TowerHeight = p%TwrFlexL
   InitOut%TowerBaseHeight = p%TowerBsHt

   ! Platform reference point wrt to global origin (0,0,0)
   InitOut%PlatformPos = x%QT(1:6)
   CALL SmllRotTrans('initial platform rotation', x%QT(4), x%QT(5), x%QT(6), TransMat, '', ErrStat2, ErrMsg2)
   InitOut%PlatformPos(1) = InitOut%PlatformPos(1) - TransMat(3,1)*p%PtfmRefzt
   InitOut%PlatformPos(2) = InitOut%PlatformPos(2) - TransMat(3,2)*p%PtfmRefzt
   InitOut%PlatformPos(3) = InitOut%PlatformPos(3) - TransMat(3,3)*p%PtfmRefzt + p%PtfmRefzt

   InitOut%HubHt            = p%HubHt
   InitOut%TwrBaseRefPos    = y%TowerLn2Mesh%Position(:,p%TwrNodes + 2)
   InitOut%TwrBaseTransDisp = y%TowerLn2Mesh%TranslationDisp(:,p%TwrNodes + 2)
   InitOut%TwrBaseRefOrient = y%TowerLn2Mesh%RefOrientation(:,:,p%TwrNodes + 2)
   InitOut%TwrBaseOrient    = y%TowerLn2Mesh%Orientation(:,:,p%TwrNodes + 2)
   InitOut%HubRad      = p%HubRad
   InitOut%RotSpeed    = p%RotSpeed
   InitOut%isFixed_GenDOF = .not. InputFileData%GenDOF
   

   if (.not. p%BD4Blades) then
      ALLOCATE(InitOut%BldRNodes(p%BldNodes),  STAT=ErrStat2)
      IF (ErrStat2 /= 0) THEN
         call CheckError( ErrStat2, ErrMsg2 )
         if (ErrStat2 >= AbortErrLev) return
      END IF
      InitOut%BldRNodes(:) = p%RNodes(:)
   else
      !Deal with BeamDyn case later
   end if

   ALLOCATE(InitOut%TwrHNodes(p%TwrNodes),  STAT=ErrStat2)
      IF (ErrStat2 /= 0) THEN
         call CheckError( ErrStat2, ErrMsg2 )
         if (ErrStat2 >= AbortErrLev) return
      END IF
   InitOut%TwrHNodes(:) = p%HNodes(:)

   CALL AllocAry(InitOut%BlPitch, p%NumBl, 'BlPitch', ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN
   InitOut%BlPitch = InputFileData%BlPitch(1:p%NumBl)

      !............................................................................................
      ! set up data needed for linearization analysis
      !............................................................................................
   
   if (InitInp%Linearize) then
      call ED_Init_Jacobian(p, u, y, InitOut, ErrStat2, ErrMsg2)
         call CheckError( ErrStat2, ErrMsg2 )
         if (ErrStat >= AbortErrLev) return
   end if
   
   
      !............................................................................................
      ! If you want to choose your own rate instead of using what the glue code suggests, tell the glue code the rate at which
      !   this module must be called here:
      !............................................................................................

   Interval = p%DT


       ! Print the summary file if requested:
   IF (InputFileData%SumPrint) THEN
      CALL ED_PrintSum( p, OtherState, ErrStat2, ErrMsg2 )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF (ErrStat >= AbortErrLev) RETURN
   END IF
       
       ! Destroy the InputFileData structure (deallocate arrays)

   CALL ED_DestroyInputFile(InputFileData, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN
      
         

CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)

      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)

      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN         
         
         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'ED_Init:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL ED_DestroyInputFile(InputFileData, ErrStat3, ErrMsg3 )
         END IF

      END IF


   END SUBROUTINE CheckError

END SUBROUTINE ED_Init
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
SUBROUTINE ED_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(ED_InputType),           INTENT(INOUT)  :: u           !< System inputs
   TYPE(ED_ParameterType),       INTENT(INOUT)  :: p           !< Parameters
   TYPE(ED_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states
   TYPE(ED_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Discrete states
   TYPE(ED_ConstraintStateType), INTENT(INOUT)  :: z           !< Constraint states
   TYPE(ED_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states
   TYPE(ED_OutputType),          INTENT(INOUT)  :: y           !< System outputs
   TYPE(ED_MiscVarType),         INTENT(INOUT)  :: m           !< Misc variables for optimization (not copied in glue code)
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Place any last minute operations or calculations here:


         ! Close files here:



         ! Destroy the input data:

      CALL ED_DestroyInput( u, ErrStat, ErrMsg )


         ! Destroy the parameter data:

      CALL ED_DestroyParam( p, ErrStat, ErrMsg )


         ! Destroy the state data:

      CALL ED_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL ED_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL ED_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL ED_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )


         ! Destroy the output data:

      CALL ED_DestroyOutput( y, ErrStat, ErrMsg )

      CALL ED_DestroyMisc( m, ErrStat, ErrMsg )



END SUBROUTINE ED_End
!----------------------------------------------------------------------------------------------------------------------------------
!> This is a loose coupling routine for solving constraint states, integrating continuous states, and updating discrete and other 
!! states. Continuous, constraint, discrete, and other states are updated to values at t + Interval.
SUBROUTINE ED_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   ) :: t          !< Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   ) :: n          !< Current simulation time step n = 0,1,...
      TYPE(ED_InputType),                 INTENT(INOUT) :: u(:)       !< Inputs at utimes (out only for mesh record-keeping in ExtrapInterp routine)
      REAL(DbKi),                         INTENT(IN   ) :: utimes(:)  !< Times associated with u(:), in seconds
      TYPE(ED_ParameterType),             INTENT(IN   ) :: p          !< Parameters
      TYPE(ED_ContinuousStateType),       INTENT(INOUT) :: x          !< Input: Continuous states at t;
                                                                      !!   Output: Continuous states at t + Interval
      TYPE(ED_DiscreteStateType),         INTENT(INOUT) :: xd         !< Input: Discrete states at t;
                                                                      !!   Output: Discrete states at t + Interval
      TYPE(ED_ConstraintStateType),       INTENT(INOUT) :: z          !< Input: Constraint states at t;
                                                                      !!   Output: Constraint states at t + Interval
      TYPE(ED_OtherStateType),            INTENT(INOUT) :: OtherState !< Other states: Other states at t;
                                                                      !!   Output: Other states at t + Interval
      TYPE(ED_MiscVarType),               INTENT(INOUT) :: m          !< Misc variables for optimization (not copied in glue code)
      INTEGER(IntKi),                     INTENT(  OUT) :: ErrStat    !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT) :: ErrMsg     !< Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""            
      

      SELECT CASE ( p%method )
         
      CASE (Method_RK4)
      
         CALL ED_RK4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
         
      CASE (Method_AB4)
      
         CALL ED_AB4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
      
      CASE (Method_ABM4)
      
         CALL ED_ABM4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
         
      CASE DEFAULT  !bjj: we already checked this at initialization, but for completeness:
         
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in ElastoDyn_UpdateStates: p%method must be 1 (RK4), 2 (AB4), or 3 (ABM4)'
         RETURN
         
      END SELECT
      
         ! Make sure the rotor azimuth is not greater or equal to 360 degrees: 
         
      ! bjj: per jmj, the subtraction of TwoPi here is so that we don't run into numerical issues with large GeAz (in large simulations)
      !   this subtraction is okay because we use x%QT(DOF_GeAz) only in equations with SIN() and/or COS() so it doesn't matter
      !   if there is a discontinunity in the channel.
      ! bjj: why don't we just do a modulo on x%QT(DOF_GeAz) instead of using x%QT(DOF_DrTr) with it?   
      
      IF ( ( x%QT(DOF_GeAz) + x%QT(DOF_DrTr) ) >= TwoPi_D )  x%QT(DOF_GeAz) = x%QT(DOF_GeAz) - TwoPi_D
            
      
END SUBROUTINE ED_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
!! This SUBROUTINE is used to compute the output channels (motions and loads) and place them in the WriteOutput() array.
!! NOTE: the descriptions of the output channels are not given here. Please see the included OutListParameters.xlsx sheet for
!! for a complete description of each output parameter.
!! NOTE: no matter how many channels are selected for output, all of the outputs are calculated
!! All of the calculated output channels are placed into the m\%AllOuts(:), while the channels selected for outputs are
!! placed in the y\%WriteOutput(:) array.
SUBROUTINE ED_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

   REAL(DbKi),                   INTENT(IN   )  :: t           !< Current simulation time in seconds
   TYPE(ED_InputType),           INTENT(IN   )  :: u           !< Inputs at Time t
   TYPE(ED_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(ED_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
   TYPE(ED_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
   TYPE(ED_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t
   TYPE(ED_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states
   TYPE(ED_OutputType),          INTENT(INOUT)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                               !!   nectivity information does not have to be recalculated)
   TYPE(ED_MiscVarType),         INTENT(INOUT)  :: m           !< Misc variables for optimization (not copied in glue code)
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


      ! Local variables:
   INTEGER, PARAMETER           :: NDims = 3

   REAL(ReKi)                   :: AngAccEB  (NDims)                               ! Angular acceleration of the base plate                                                (body B) in the inertia frame (body E for earth).
   REAL(ReKi)                   :: AngAccEN  (NDims)                               ! Angular acceleration of the nacelle                                                   (body N) in the inertia frame (body E for earth).
   REAL(ReKi)                   :: AngAccEH  (NDims)                               ! Angular acceleration of the hub                                                   (body N) in the inertia frame (body E for earth).
   REAL(ReKi)                   :: AngAccER  (NDims)                               ! Angular acceleration of the structure that furls with the rotor (not including rotor) (body R) in the inertia frame (body E for earth).
   REAL(ReKi)                   :: AngAccEX  (NDims)                               ! Angular acceleration of the platform                                                  (body X) in the inertia frame (body E for earth).
!   REAL(ReKi)                   :: ComDenom                                        ! Common denominator used in several expressions.
!   REAL(ReKi)                   :: CThrstys                                        ! Estimate of the ys-location of the center of thrust.
!   REAL(ReKi)                   :: CThrstzs                                        ! Estimate of the zs-location of the center of thrust.
   REAL(R8Ki)                   :: FrcMGagB  (NDims)                               ! Total force at the blade element   (body M) / blade strain gage location            (point S) due to the blade above the strain gage.
   REAL(ReKi)                   :: FrcFGagT  (NDims)                               ! Total force at the tower element   (body F) / tower strain gage location            (point T) due to the nacelle and rotor and tower above the strain gage.
   REAL(ReKi)                   :: FrcONcRt  (NDims)                               ! Total force at the yaw bearing (point O  ) due to the nacelle, generator, and rotor
   REAL(ReKi)                   :: FrcPRot   (NDims)                               ! Total force at the teeter pin  (point P  ) due to the rotor
   REAL(ReKi)                   :: FrcT0Trb  (NDims)                               ! Total force at the base of flexible portion of the tower (point T(0)) due to the entire wind turbine
   REAL(ReKi)                   :: FZHydro   (NDims)                               ! Total platform hydrodynamic force at the platform reference (point Z)
!   REAL(ReKi)                   :: HHWndVec  (NDims)                               ! Hub-height wind vector in the AeroDyn coordinate system
   REAL(ReKi)                   :: LinAccEIMU(NDims)                               ! Total linear acceleration of the nacelle IMU (point IMU) in the inertia frame (body E for earth)
   REAL(ReKi)                   :: LinAccEO  (NDims)                               ! Total linear acceleration of the base plate (point O) in the inertia frame (body E for earth)
   REAL(ReKi)                   :: LinAccEZ  (NDims)                               ! Total linear acceleration of the platform refernce (point Z) in the inertia frame (body E for earth)
   REAL(ReKi)                   :: MomBNcRt  (NDims)                               ! Total moment at the base plate      (body B) / yaw bearing                           (point O) due to the nacelle, generator, and rotor.
   REAL(ReKi)                   :: MomFGagT  (NDims)                               ! Total moment at the tower element   (body F) / tower strain gage location            (point T) due to the nacelle and rotor and tower above the strain gage.
   REAL(ReKi)                   :: MomLPRot  (NDims)                               ! Total moment at the low-speed shaft (body L) / teeter pin                            (point P) due to the rotor.
   REAL(ReKi)                   :: MomMGagB  (NDims)                               ! Total moment at the blade element   (body M) / blade strain gage location            (point S) due to the blade above the strain gage.
   REAL(ReKi)                   :: MomNGnRt  (NDims)                               ! Total moment at the nacelle         (body N) / specified point on rotor-furl axis    (point V) due to the structure that furls with the rotor, generator, and rotor.
   REAL(ReKi)                   :: MomNTail  (NDims)                               ! Total moment at the nacelle         (body N) / specified point on  tail-furl axis    (point W) due to the tail.
   REAL(ReKi)                   :: MomX0Trb  (NDims)                               ! Total moment at the tower base      (body X) / base of flexible portion of the tower (point T(0)) due to the entire wind turbine.
   REAL(ReKi)                   :: MXHydro   (NDims)                               ! Total platform hydrodynamic moment acting at the platform (body X) / platform reference (point Z).
   REAL(R8Ki)                   :: rOPO      (NDims)                               ! Position vector from the undeflected tower top (point O prime) to the deflected tower top (point O).
   REAL(R8Ki)                   :: rOSTip    (NDims)                               ! Position vector from the deflected tower top (point O) to the deflected blade tip (point S tip).
   REAL(R8Ki)                   :: rOSTipxn                                        ! Component of rOSTip directed along the xn-axis.
   REAL(R8Ki)                   :: rOSTipyn                                        ! Component of rOSTip directed along the yn-axis.
   REAL(R8Ki)                   :: rOSTipzn                                        ! Component of rOSTip directed along the zn-axis.
   REAL(R8Ki)                   :: rTPT      (NDims)                               ! Position vector from the undeflected tower node (point T prime) to the deflected node (point T)
   REAL(R8Ki)                   :: rSPS      (NDims)                               ! Position vector from the undeflected blade node (point S prime) to the deflected node (point S)
   REAL(R8Ki)                   :: rSTipPSTip(NDims)                               ! Position vector from the undeflected blade tip (point S tip prime) to the deflected blade tip (point S tip).
   REAL(R8Ki)                   :: TmpVec    (NDims)                               ! A temporary vector used in various computations.
   REAL(R8Ki)                   :: TmpVec2   (NDims)                               ! A temporary vector.

   REAL(ReKi)                   :: LinAccES (NDims,0:p%TipNode,p%NumBl)            ! Total linear acceleration of a point on a   blade (point S) in the inertia frame (body E for earth).
   REAL(ReKi)                   :: AngAccEK (NDims,0:p%TipNode,p%NumBl)            ! Total rotational acceleration of a point on a blade (point S) in the inertia frame (body E for earth).
   REAL(ReKi)                   :: LinAccET (NDims,0:p%TwrNodes)                   ! Total linear acceleration of a point on the tower (point T) in the inertia frame (body E for earth).
   REAL(ReKi)                   :: AngAccEF (NDims,0:p%TwrNodes)                   ! Total angular acceleration of tower element J (body F) in the inertia frame (body E for earth).
   REAL(ReKi)                   :: FrcS0B   (NDims,p%NumBl)                        ! Total force at the blade root (point S(0)) due to the blade.
   REAL(ReKi)                   :: FTTower  (NDims,p%TwrNodes)                     ! Total hydrodynamic + aerodynamic force per unit length acting on the tower at point T.
   REAL(ReKi)                   :: MFHydro  (NDims,p%TwrNodes)                     ! Total hydrodynamic + aerodynamic moment per unit length acting on a tower element (body F) at point T.
   REAL(ReKi)                   :: MomH0B   (NDims,p%NumBl)                        ! Total moment at the hub (body H) / blade root (point S(0)) due to the blade.


   INTEGER(IntKi)               :: I                                               ! Generic index
   INTEGER(IntKi)               :: J, J2                                           ! Loops through nodes / elements
   INTEGER(IntKi)               :: K                                               ! Loops through blades
   INTEGER(IntKi)               :: NodeNum                                         ! Mesh node number for given blade/node
   
   INTEGER(IntKi)               :: ErrStat2                                        ! Temporary Error code
   CHARACTER(ErrMsgLen)         :: ErrMsg2                                         ! Temporary error message
   

   LOGICAL, PARAMETER           :: UpdateValues  = .TRUE.                          ! determines if the OtherState values need to be updated
   TYPE(ED_ContinuousStateType) :: dxdt                                            ! Continuous state derivs at t

         ! Initialize some output values
      ErrStat = ErrID_None
      ErrMsg  = ""

      
      ! SEE IF THESE NEED TO BE CALLED (i.e., if UpdateStates was called, these values are already calculated)
   IF ( UpdateValues ) THEN    
         ! Update the OtherState data by calculating the derivative...
      !OtherState%HSSBrTrqC = SIGN( u%HSSBrTrqC, x%QDT(DOF_GeAz) )
      CALL ED_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg ) ! sets m%QD2T = dxdt%QDT
      CALL ED_DestroyContState( dxdt, ErrStat2, ErrMsg2 )  
      IF (ErrStat >= AbortErrLev) RETURN
   END IF      


   
      !..............................
      ! Outputs for HydroDyn, UsrPtfm and UsrTwr
      !..............................
            
      !u_UsrPtfm%X  = x_ED%QT(1:6)
      !u_UsrPtfm%XD = x_ED%QDT(1:6)

      !u_UsrTwr%X(:,J) = (/ OtherSt_ED%RtHS%rT(1,J),       -OtherSt_ED%RtHS%rT(3,J),       OtherSt_ED%RtHS%rT( 2,J)+ p%PtfmRefzt,&
      !                     OtherSt_ED%RtHS%AngPosEF(1,J), -OtherSt_ED%RtHS%AngPosEF(3,J), OtherSt_ED%RtHS%AngPosEF(2,J)         /) 
      !u_UsrTwr%XD(:,J) = (/ OtherSt_ED%RtHS%LinVelET(1,J), -OtherSt_ED%RtHS%LinVelET(3,J), OtherSt_ED%RtHS%LinVelET(2,J),&
      !                      OtherSt_ED%RtHS%AngVelEF(1,J), -OtherSt_ED%RtHS%AngVelEF(3,J), OtherSt_ED%RtHS%AngVelEF(2,J) /)                     
   
   
   
   
   ! Array m%AllOuts() is initialized to 0.0 in initialization, so we are not going to reinitialize it here.

   !...............................................................................................................................
   ! Calculate all of the total forces and moments using all of the partial forces and moments calculated in RtHS().  Also,
   !   calculate all of the total angular and linear accelerations using all of the partial accelerations calculated in RtHS().
   !   To do this, first initialize the variables using the portions not associated with the accelerations.  Then add the portions
   !   associated with the accelerations one by one:
   !...............................................................................................................................

   AngAccEB   = m%RtHS%AngAccEBt
   AngAccEH   = m%RtHS%AngAccEHt
   AngAccEN   = m%RtHS%AngAccENt
   AngAccER   = m%RtHS%AngAccERt
   AngAccEX   = m%RtHS%AngAccEXt
   LinAccEIMU = m%RtHS%LinAccEIMUt
   LinAccEO   = m%RtHS%LinAccEOt
   LinAccEZ   = m%RtHS%LinAccEZt
   FrcONcRt   = m%RtHS%FrcONcRtt
   FrcPRot    = m%RtHS%FrcPRott
   FrcT0Trb   = m%RtHS%FrcT0Trbt

   ! was FZHydro    = m%RtHS%FZHydrot
   FZHydro    = u%PlatformPtMesh%Force(DOF_Sg,1)*m%RtHS%PLinVelEZ(DOF_Sg,0,:) &
              + u%PlatformPtMesh%Force(DOF_Sw,1)*m%RtHS%PLinVelEZ(DOF_Sw,0,:) &
              + u%PlatformPtMesh%Force(DOF_Hv,1)*m%RtHS%PLinVelEZ(DOF_Hv,0,:)

   MomBNcRt   = m%RtHS%MomBNcRtt
   MomLPRot   = m%RtHS%MomLPRott
   MomNGnRt   = m%RtHS%MomNGnRtt
   MomNTail   = m%RtHS%MomNTailt
   MomX0Trb   = m%RtHS%MomX0Trbt

   ! was MXHydro = m%RtHS%MXHydrot
   MXHydro    = u%PlatformPtMesh%Moment(DOF_R-3,1)*m%RtHS%PAngVelEX(DOF_R ,0,:) &
              + u%PlatformPtMesh%Moment(DOF_P-3,1)*m%RtHS%PAngVelEX(DOF_P ,0,:) &
              + u%PlatformPtMesh%Moment(DOF_Y-3,1)*m%RtHS%PAngVelEX(DOF_Y ,0,:)

   DO I = 1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs
      AngAccEB   = AngAccEB   + m%RtHS%PAngVelEB  (p%DOFs%SrtPS(I),0,:)*m%QD2T(p%DOFs%SrtPS(I))
      AngAccEH   = AngAccEH   + m%RtHS%PAngVelEH  (p%DOFs%SrtPS(I),0,:)*m%QD2T(p%DOFs%SrtPS(I))      
      AngAccEN   = AngAccEN   + m%RtHS%PAngVelEN  (p%DOFs%SrtPS(I),0,:)*m%QD2T(p%DOFs%SrtPS(I))      
      AngAccER   = AngAccER   + m%RtHS%PAngVelER  (p%DOFs%SrtPS(I),0,:)*m%QD2T(p%DOFs%SrtPS(I))
      LinAccEIMU = LinAccEIMU + m%RtHS%PLinVelEIMU(p%DOFs%SrtPS(I),0,:)*m%QD2T(p%DOFs%SrtPS(I))
      LinAccEO   = LinAccEO   + m%RtHS%PLinVelEO  (p%DOFs%SrtPS(I),0,:)*m%QD2T(p%DOFs%SrtPS(I))
      FrcONcRt   = FrcONcRt   + m%RtHS%PFrcONcRt  (:,p%DOFs%SrtPS(I)  )*m%QD2T(p%DOFs%SrtPS(I))
      FrcPRot    = FrcPRot    + m%RtHS%PFrcPRot   (:,p%DOFs%SrtPS(I)  )*m%QD2T(p%DOFs%SrtPS(I))
      FrcT0Trb   = FrcT0Trb   + m%RtHS%PFrcT0Trb  (:,p%DOFs%SrtPS(I)  )*m%QD2T(p%DOFs%SrtPS(I))
      MomBNcRt   = MomBNcRt   + m%RtHS%PMomBNcRt  (:,p%DOFs%SrtPS(I)  )*m%QD2T(p%DOFs%SrtPS(I))
      MomLPRot   = MomLPRot   + m%RtHS%PMomLPRot  (:,p%DOFs%SrtPS(I)  )*m%QD2T(p%DOFs%SrtPS(I))
      MomNGnRt   = MomNGnRt   + m%RtHS%PMomNGnRt  (:,p%DOFs%SrtPS(I)  )*m%QD2T(p%DOFs%SrtPS(I))
      MomNTail   = MomNTail   + m%RtHS%PMomNTail  (:,p%DOFs%SrtPS(I)  )*m%QD2T(p%DOFs%SrtPS(I))
      MomX0Trb   = MomX0Trb   + m%RtHS%PMomX0Trb  (:,p%DOFs%SrtPS(I)  )*m%QD2T(p%DOFs%SrtPS(I))
   ENDDO             ! I - All active (enabled) DOFs
   DO I = 1,p%DOFs%NPYE     ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the platform center of mass (point Y)
      AngAccEX   = AngAccEX   + m%RtHS%PAngVelEX  (p%DOFs%PYE  (I),0,:)*m%QD2T(p%DOFs%PYE  (I))
      LinAccEZ   = LinAccEZ   + m%RtHS%PLinVelEZ  (p%DOFs%PYE  (I),0,:)*m%QD2T(p%DOFs%PYE  (I))


      FZHydro    = FZHydro    + (- u%PtfmAddedMass(DOF_Sg,p%DOFs%PYE(I))*m%RtHS%PLinVelEZ(DOF_Sg,0,:)   &  ! was  FZHydro = FZHydro + m%RtHS%PFZHydro(p%DOFs%PYE(I),:)*m%QD2T(p%DOFs%PYE  (I))
                                 - u%PtfmAddedMass(DOF_Sw,p%DOFs%PYE(I))*m%RtHS%PLinVelEZ(DOF_Sw,0,:)   &
                                 - u%PtfmAddedMass(DOF_Hv,p%DOFs%PYE(I))*m%RtHS%PLinVelEZ(DOF_Hv,0,:) ) &
                                *m%QD2T(p%DOFs%PYE  (I))
      ! was MXHydro = MXHydro    +  m%RtHS%PMXHydro   (p%DOFs%PYE  (I),  :)*m%QD2T(p%DOFs%PYE  (I))
      MXHydro    = MXHydro    +  (- u%PtfmAddedMass(DOF_R ,p%DOFs%PYE(I))*m%RtHS%PAngVelEX(DOF_R ,0,:)   &
                                  - u%PtfmAddedMass(DOF_P ,p%DOFs%PYE(I))*m%RtHS%PAngVelEX(DOF_P ,0,:)   &
                                  - u%PtfmAddedMass(DOF_Y ,p%DOFs%PYE(I))*m%RtHS%PAngVelEX(DOF_Y ,0,:) ) &
                                *m%QD2T(p%DOFs%PYE  (I))

   ENDDO             ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the platform center of mass (point Y)



   DO K = 1,p%NumBl ! Loop through all blades

      FrcS0B  (:          ,K) = m%RtHS%FrcS0Bt  (:,K          )
      MomH0B  (:          ,K) = m%RtHS%MomH0Bt  (:,K          )

      DO I = 1,p%DOFs%NPSE(K)  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of blade K
         FrcS0B  (:          ,K) = FrcS0B  (:          ,K) + m%RtHS%PFrcS0B  (:,K,          p%DOFs%PSE(K,I)  )*m%QD2T(p%DOFs%PSE(K,I))
         MomH0B  (:          ,K) = MomH0B  (:          ,K) + m%RtHS%PMomH0B  (:,K,          p%DOFs%PSE(K,I)  )*m%QD2T(p%DOFs%PSE(K,I))
      ENDDO             ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of blade K

      DO J = 0,p%TipNode ! Loop through the blade nodes / elements

         LinAccES(:,J,K) = m%RtHS%LinAccESt(:,K,J)
         AngAccEK(:,J,K) = m%RtHs%AngAccEKt(:,J,K)
         DO I = 1,p%DOFs%NPSE(K)  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of blade K
            LinAccES(:,J,K) = LinAccES(:,J,K) + m%RtHS%PLinVelES(K,J,p%DOFs%PSE(K,I),0,:)*m%QD2T(p%DOFs%PSE(K,I))
            AngAccEK(:,J,K) = AngAccEK(:,J,K) + m%RtHS%PAngVelEM(K,J,p%DOFs%PSE(K,I),0,:)*m%QD2T(p%DOFs%PSE(K,I))
         ENDDO             ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of blade K

      ENDDO             ! J - Blade nodes / elements

   ENDDO          ! K - All blades

   DO J = 0,p%TwrNodes  ! Loop through the tower nodes / elements, starting at the tower base (0)

      LinAccET(:,J) = m%RtHS%LinAccETt(:,J)
      AngAccEF(:,J) = m%RtHS%AngAccEFt(:,J)

      DO I = 1,p%DOFs%NPTE  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the yaw bearing center of mass (point O)
         LinAccET(:,J) = LinAccET(:,J) + m%RtHS%PLinVelET(J,p%DOFs%PTE(I),0,:)*m%QD2T(p%DOFs%PTE(I))
         AngAccEF(:,J) = AngAccEF(:,J) + m%RtHS%PAngVelEF(J,p%DOFs%PTE(I),0,:)*m%QD2T(p%DOFs%PTE(I))
      ENDDO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the yaw bearing center of mass (point O)

   ENDDO ! J - Tower nodes / elements


   DO J = 1,p%TwrNodes  ! Loop through the tower nodes / elements

      FTTower (:,J) = m%RtHS%FTHydrot (:,J)
      MFHydro (:,J) = m%RtHS%MFHydrot (:,J)

      DO I = 1,p%DOFs%NPTE  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the yaw bearing center of mass (point O)
         FTTower (:,J) = FTTower (:,J) + m%RtHS%PFTHydro (:,J,p%DOFs%PTE(I)  )*m%QD2T(p%DOFs%PTE(I))
         MFHydro (:,J) = MFHydro (:,J) + m%RtHS%PMFHydro (:,J,p%DOFs%PTE(I)  )*m%QD2T(p%DOFs%PTE(I))
      ENDDO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the yaw bearing center of mass (point O)

   ENDDO ! J - Tower nodes / elements


      ! Convert the units of the forces and moments from N and N-m
      !    to kN and kN-m:

   FrcONcRt = 0.001*FrcONcRt
   FrcPRot  = 0.001*FrcPRot
   FrcT0Trb = 0.001*FrcT0Trb
   FZHydro  = 0.001*FZHydro
   MomBNcRt = 0.001*MomBNcRt
   MomLPRot = 0.001*MomLPRot
   MomNGnRt = 0.001*MomNGnRt
   MomNTail = 0.001*MomNTail
   MomX0Trb = 0.001*MomX0Trb
   MXHydro  = 0.001*MXHydro
   FrcS0B   = 0.001*FrcS0B
   MomH0B   = 0.001*MomH0B


   !...............................................................................................................................
   ! set the values in the AllOuts array:
   !...............................................................................................................................
      ! Define the output channel specifying the current simulation time:

   m%AllOuts(  Time) = REAL( t, ReKi )


      ! Blade (1-3) Tip Motions:

   DO K = 1,p%NumBl
      rSTipPSTip = m%RtHS%rS0S(:,K,p%TipNode) - p%BldFlexL*m%CoordSys%j3(K,:)  ! Position vector from the undeflected blade tip (point S tip prime) to the deflected blade tip (point S tip) of blade 1.
      rOSTip     = m%RtHS%rS  (:,K,p%TipNode) - m%RtHS%rO                ! Position vector from the deflected tower top (point O) to the deflected blade tip (point S tip) of blade 1.
      rOSTipxn   =      DOT_PRODUCT( rOSTip, m%CoordSys%d1 )                ! Component of rOSTip directed along the xn-axis.
      rOSTipyn   = -1.0*DOT_PRODUCT( rOSTip, m%CoordSys%d3 )                ! Component of rOSTip directed along the yn-axis.
      rOSTipzn   =      DOT_PRODUCT( rOSTip, m%CoordSys%d2 )                ! Component of rOSTip directed along the zn-axis.

      IF (.NOT. p%BD4Blades) THEN
         m%AllOuts(  TipDxc(K) ) = DOT_PRODUCT(            rSTipPSTip, m%CoordSys%i1(K,         :) )
         m%AllOuts(  TipDyc(K) ) = DOT_PRODUCT(            rSTipPSTip, m%CoordSys%i2(K,         :) )
         m%AllOuts(  TipDzc(K) ) = DOT_PRODUCT(            rSTipPSTip, m%CoordSys%i3(K,         :) )
         m%AllOuts(  TipDxb(K) ) = DOT_PRODUCT(            rSTipPSTip, m%CoordSys%j1(K,         :) )
         m%AllOuts(  TipDyb(K) ) = DOT_PRODUCT(            rSTipPSTip, m%CoordSys%j2(K,         :) )
      !JASON: USE TipNode HERE INSTEAD OF BldNodes IF YOU ALLOCATE AND DEFINE n1, n2, n3, m1, m2, AND m3 TO USE TipNode.  THIS WILL REQUIRE THAT THE AERODYNAMIC AND STRUCTURAL TWISTS, AeroTwst() AND ThetaS(), BE KNOWN AT THE TIP!!!
         m%AllOuts( TipALxb(K) ) = DOT_PRODUCT( LinAccES(:,p%TipNode,K), m%CoordSys%n1(K,p%BldNodes,:) )
         m%AllOuts( TipALyb(K) ) = DOT_PRODUCT( LinAccES(:,p%TipNode,K), m%CoordSys%n2(K,p%BldNodes,:) )
         m%AllOuts( TipALzb(K) ) = DOT_PRODUCT( LinAccES(:,p%TipNode,K), m%CoordSys%n3(K,p%BldNodes,:) )
         m%AllOuts( TipRDxb(K) ) = DOT_PRODUCT( m%RtHS%AngPosHM(:,K,p%TipNode), m%CoordSys%j1(K,         :) )*R2D
         m%AllOuts( TipRDyb(K) ) = DOT_PRODUCT( m%RtHS%AngPosHM(:,K,p%TipNode), m%CoordSys%j2(K,         :) )*R2D
         ! There is no sense computing AllOuts( TipRDzc(K) ) here since it is always zero for FAST simulation results.
         IF ( p%MHK == 2 ) THEN
            IF ( rOSTipzn < 0.0 )  THEN   ! Tip of blade K is above the yaw bearing.
               m%AllOuts(TipClrnc(K) ) = SQRT( rOSTipxn*rOSTipxn + rOSTipyn*rOSTipyn + rOSTipzn*rOSTipzn ) ! Absolute distance from the tower top / yaw bearing to the tip of blade 1.
            ELSE                          ! Tip of blade K is below the yaw bearing.
               m%AllOuts(TipClrnc(K) ) = SQRT( rOSTipxn*rOSTipxn + rOSTipyn*rOSTipyn                     ) ! Perpendicular distance from the yaw axis / tower centerline to the tip of blade 1.
            ENDIF
         ELSE
            IF ( rOSTipzn > 0.0 )  THEN   ! Tip of blade K is above the yaw bearing.
               m%AllOuts(TipClrnc(K) ) = SQRT( rOSTipxn*rOSTipxn + rOSTipyn*rOSTipyn + rOSTipzn*rOSTipzn ) ! Absolute distance from the tower top / yaw bearing to the tip of blade 1.
            ELSE                          ! Tip of blade K is below the yaw bearing.
               m%AllOuts(TipClrnc(K) ) = SQRT( rOSTipxn*rOSTipxn + rOSTipyn*rOSTipyn                     ) ! Perpendicular distance from the yaw axis / tower centerline to the tip of blade 1.
            ENDIF
         ENDIF
      END IF      

   END DO !K

      ! Blade (1-3) Local Span Motions:

   DO K = 1,p%NumBl
      DO I = 1, p%NBlGages

         m%AllOuts( SpnALxb(I,K) ) = DOT_PRODUCT( LinAccES(:,p%BldGagNd(I),K), m%CoordSys%n1(K,p%BldGagNd(I),:) )
         m%AllOuts( SpnALyb(I,K) ) = DOT_PRODUCT( LinAccES(:,p%BldGagNd(I),K), m%CoordSys%n2(K,p%BldGagNd(I),:) )
         m%AllOuts( SpnALzb(I,K) ) = DOT_PRODUCT( LinAccES(:,p%BldGagNd(I),K), m%CoordSys%n3(K,p%BldGagNd(I),:) )

         rSPS                      = m%RtHS%rS0S(:,K,p%BldGagNd(I)) - p%RNodes(p%BldGagNd(I))*m%CoordSys%j3(K,:)

         m%AllOuts( SpnTDxb(I,K) ) = DOT_PRODUCT( rSPS, m%CoordSys%j1(K,:) )
         m%AllOuts( SpnTDyb(I,K) ) = DOT_PRODUCT( rSPS, m%CoordSys%j2(K,:) )
         m%AllOuts( SpnTDzb(I,K) ) = DOT_PRODUCT( rSPS, m%CoordSys%j3(K,:) )

         m%AllOuts( SpnRDxb(I,K) ) = DOT_PRODUCT( m%RtHS%AngPosHM(:,K,p%BldGagNd(I)), m%CoordSys%j1(K,:) )*R2D
         m%AllOuts( SpnRDyb(I,K) ) = DOT_PRODUCT( m%RtHS%AngPosHM(:,K,p%BldGagNd(I)), m%CoordSys%j2(K,:) )*R2D
        !m%AllOuts( SpnRDzb(I,K) ) = DOT_PRODUCT( m%RtHS%AngPosHM(:,K,p%BldGagNd(I)), m%CoordSys%j3(K,:) )*R2D           ! this is always zero for FAST

      END DO !I
   END DO !K



      ! Blade Pitch Motions:

   m%AllOuts(PtchPMzc1) = u%BlPitchCom(1)*R2D
IF ( p%NumBl > 1 ) THEN
   m%AllOuts(PtchPMzc2) = u%BlPitchCom(2)*R2D
   IF ( p%NumBl > 2 )  THEN ! 3-blader

      m%AllOuts(PtchPMzc3) = u%BlPitchCom(3)*R2D

   ELSE  ! 2-blader


      ! Teeter Motions:

      m%AllOuts(  TeetPya) = x%QT  (DOF_Teet)*R2D
      m%AllOuts(  TeetVya) = x%QDT (DOF_Teet)*R2D
      m%AllOuts(  TeetAya) = m%QD2T(DOF_Teet)*R2D

   ENDIF
END IF

      ! Shaft Motions:

   y%LSSTipPxa = x%QT (DOF_GeAz) + x%QT  (DOF_DrTr) + p%AzimB1Up + PiBy2
   if (.not. m%IgnoreMod) CALL Zero2TwoPi(y%LSSTipPxa)  ! Return value between 0 and 2pi (LSSTipPxa is used only in calculations of SIN and COS, so it's okay to take MOD/MODULO here; this wouldn't be oaky for linearization)
   m%AllOuts(LSSTipPxa) = y%LSSTipPxa*R2D
   
   m%AllOuts(LSSGagPxa) = x%QT (DOF_GeAz) + p%AzimB1Up + PiBy2 
   if (.not. m%IgnoreMod) CALL Zero2TwoPi(m%AllOuts(LSSGagPxa))  ! Return value between 0 and 2pi 
   m%AllOuts(LSSGagPxa) = m%AllOuts(LSSGagPxa)*R2D ! convert to degrees
   
   m%AllOuts(   LSSTipVxa) =      (     x%QDT (DOF_GeAz) +          x%QDT (DOF_DrTr) )*RPS2RPM
   m%AllOuts(   LSSTipAxa) = ( m%QD2T(DOF_GeAz) + m%QD2T(DOF_DrTr) )*R2D
   m%AllOuts(   LSSGagVxa) =            x%QDT (DOF_GeAz)                              *RPS2RPM
   m%AllOuts(   LSSGagAxa) =   m%QD2T(DOF_GeAz)                              *R2D
   m%AllOuts(     HSShftV) = ABS(p%GBRatio)*m%AllOuts(LSSGagVxa)
   m%AllOuts(     HSShftA) = ABS(p%GBRatio)*m%AllOuts(LSSGagAxa)

   !IF ( .NOT. EqualRealNos( m%AllOuts(WindVxi), 0.0_ReKi ) )  THEN  ! .TRUE. if the denominator in the following equation is not zero.
   !   m%AllOuts(TipSpdRat) =      ( x%QDT (DOF_GeAz) + x%QDT (DOF_DrTr) )*p%AvgNrmTpRd / m%AllOuts(  WindVxi)
   !ELSE
   !   m%AllOuts(TipSpdRat) = 0.0
   !ENDIF


      ! Nacelle IMU Motions:

   m%AllOuts(NcIMUTVxs) =      DOT_PRODUCT( m%RtHS%LinVelEIMU, m%CoordSys%c1 )
   m%AllOuts(NcIMUTVys) = -1.0*DOT_PRODUCT( m%RtHS%LinVelEIMU, m%CoordSys%c3 )
   m%AllOuts(NcIMUTVzs) =      DOT_PRODUCT( m%RtHS%LinVelEIMU, m%CoordSys%c2 )
   m%AllOuts(NcIMUTAxs) =      DOT_PRODUCT(                 LinAccEIMU, m%CoordSys%c1 )
   m%AllOuts(NcIMUTAys) = -1.0*DOT_PRODUCT(                 LinAccEIMU, m%CoordSys%c3 )
   m%AllOuts(NcIMUTAzs) =      DOT_PRODUCT(                 LinAccEIMU, m%CoordSys%c2 )
   m%AllOuts(NcIMURVxs) =      DOT_PRODUCT( m%RtHS%AngVelER  , m%CoordSys%c1 )*R2D
   m%AllOuts(NcIMURVys) = -1.0*DOT_PRODUCT( m%RtHS%AngVelER  , m%CoordSys%c3 )*R2D
   m%AllOuts(NcIMURVzs) =      DOT_PRODUCT( m%RtHS%AngVelER  , m%CoordSys%c2 )*R2D
   m%AllOuts(NcIMURAxs) =      DOT_PRODUCT(                 AngAccER  , m%CoordSys%c1 )*R2D
   m%AllOuts(NcIMURAys) = -1.0*DOT_PRODUCT(                 AngAccER  , m%CoordSys%c3 )*R2D
   m%AllOuts(NcIMURAzs) =      DOT_PRODUCT(                 AngAccER  , m%CoordSys%c2 )*R2D


      ! Rotor-Furl Motions:

   m%AllOuts( RotFurlP) = x%QT  (DOF_RFrl)*R2D
   m%AllOuts( RotFurlV) = x%QDT (DOF_RFrl)*R2D
   m%AllOuts( RotFurlA) = m%QD2T(DOF_RFrl)*R2D


      ! Tail-Furl Motions:

   m%AllOuts(TailFurlP) = x%QT  (DOF_TFrl)*R2D
   m%AllOuts(TailFurlV) = x%QDT (DOF_TFrl)*R2D
   m%AllOuts(TailFurlA) = m%QD2T(DOF_TFrl)*R2D


      ! Yaw Motions:

   m%AllOuts(   YawPzn) = x%QT  (DOF_Yaw )*R2D
   m%AllOuts(   YawVzn) = x%QDT (DOF_Yaw )*R2D
   m%AllOuts(   YawAzn) = m%QD2T(DOF_Yaw )*R2D


   ! Tower-Top / Yaw Bearing Motions:

   rOPO     = m%RtHS%rT0O - p%TwrFlexL*m%CoordSys%a2 ! Position vector from the undeflected tower top (point O prime) to the deflected tower top (point O).

   ! p%TwrNodes+1 is the tower top:
   J = p%TwrNodes+1
   m%AllOuts(TwrTpTDxi) =     m%RtHS%rO(1) - y%TowerLn2Mesh%Position(1,J)
   m%AllOuts(TwrTpTDyi) = -1.*m%RtHS%rO(3) - y%TowerLn2Mesh%Position(2,J)
   m%AllOuts(TwrTpTDzi) =     m%RtHS%rO(2) - y%TowerLn2Mesh%Position(3,J) + p%PtfmRefzt
   m%AllOuts(YawBrTDxp) =  DOT_PRODUCT(     rOPO, m%CoordSys%b1 )
   m%AllOuts(YawBrTDyp) = -DOT_PRODUCT(     rOPO, m%CoordSys%b3 )
   m%AllOuts(YawBrTDzp) =  DOT_PRODUCT(     rOPO, m%CoordSys%b2 )
   m%AllOuts(YawBrTDxt) =  DOT_PRODUCT(     rOPO, m%CoordSys%a1 )
   m%AllOuts(YawBrTDyt) = -DOT_PRODUCT(     rOPO, m%CoordSys%a3 )
   m%AllOuts(YawBrTDzt) =  DOT_PRODUCT(     rOPO, m%CoordSys%a2 )
   m%AllOuts(YawBrTVxp) =  DOT_PRODUCT( m%RtHS%LinVelEO, m%CoordSys%b1 )
   m%AllOuts(YawBrTVyp) = -DOT_PRODUCT( m%RtHS%LinVelEO, m%CoordSys%b3 )
   m%AllOuts(YawBrTVzp) =  DOT_PRODUCT( m%RtHS%LinVelEO, m%CoordSys%b2 )  
   m%AllOuts(YawBrTAxp) =  DOT_PRODUCT( LinAccEO, m%CoordSys%b1 )
   m%AllOuts(YawBrTAyp) = -DOT_PRODUCT( LinAccEO, m%CoordSys%b3 )
   m%AllOuts(YawBrTAzp) =  DOT_PRODUCT( LinAccEO, m%CoordSys%b2 )
   m%AllOuts(YawBrRDxt) =  DOT_PRODUCT( m%RtHS%AngPosXB, m%CoordSys%a1 )*R2D
   m%AllOuts(YawBrRDyt) = -DOT_PRODUCT( m%RtHS%AngPosXB, m%CoordSys%a3 )*R2D
   ! There is no sense computing m%AllOuts(YawBrRDzt) here since it is always zero for FAST simulation results.
   m%AllOuts(YawBrRVxp) =  DOT_PRODUCT( m%RtHS%AngVelEB, m%CoordSys%b1 )*R2D
   m%AllOuts(YawBrRVyp) = -DOT_PRODUCT( m%RtHS%AngVelEB, m%CoordSys%b3 )*R2D
   m%AllOuts(YawBrRVzp) =  DOT_PRODUCT( m%RtHS%AngVelEB, m%CoordSys%b2 )*R2D
   m%AllOuts(YawBrRAxp) =  DOT_PRODUCT( AngAccEB, m%CoordSys%b1 )*R2D
   m%AllOuts(YawBrRAyp) = -DOT_PRODUCT( AngAccEB, m%CoordSys%b3 )*R2D
   m%AllOuts(YawBrRAzp) =  DOT_PRODUCT( AngAccEB, m%CoordSys%b2 )*R2D


      ! Local Tower Motions:

   DO I = 1, p%NTwGages

      m%AllOuts( TwHtALxt(I) ) =      DOT_PRODUCT( LinAccET(:,p%TwrGagNd(I)), m%CoordSys%t1(p%TwrGagNd(I),:) )
      m%AllOuts( TwHtALyt(I) ) = -1.0*DOT_PRODUCT( LinAccET(:,p%TwrGagNd(I)), m%CoordSys%t3(p%TwrGagNd(I),:) )
      m%AllOuts( TwHtALzt(I) ) =      DOT_PRODUCT( LinAccET(:,p%TwrGagNd(I)), m%CoordSys%t2(p%TwrGagNd(I),:) )

      rTPT                   = m%RtHS%rT0T(:,p%TwrGagNd(I)) - p%HNodes(p%TwrGagNd(I))*m%CoordSys%a2(:)

      m%AllOuts( TwHtTDxt(I) ) =      DOT_PRODUCT( rTPT,     m%CoordSys%a1 )
      m%AllOuts( TwHtTDyt(I) ) = -1.0*DOT_PRODUCT( rTPT,     m%CoordSys%a3 )
      m%AllOuts( TwHtTDzt(I) ) =      DOT_PRODUCT( rTPT,     m%CoordSys%a2 )

      m%AllOuts( TwHtRDxt(I) ) =      DOT_PRODUCT( m%RtHS%AngPosXF(:,p%TwrGagNd(I)), m%CoordSys%a1 )*R2D  !why is this zero???
      m%AllOuts( TwHtRDyt(I) ) = -1.0*DOT_PRODUCT( m%RtHS%AngPosXF(:,p%TwrGagNd(I)), m%CoordSys%a3 )*R2D
   !   m%AllOuts( TwHtRDzt(I) ) =     DOT_PRODUCT( m%RtHS%AngPosXF(:,p%TwrGagNd(I)), m%CoordSys%a2 )*R2D  !this will always be 0 in FAST, so no need to calculate


      m%AllOuts( TwHtTPxi(I) ) =      m%RtHS%rT(1,p%TwrGagNd(I))
      m%AllOuts( TwHtTPyi(I) ) = -1.0*m%RtHS%rT(3,p%TwrGagNd(I))
      m%AllOuts( TwHtTPzi(I) ) =      m%RtHS%rT(2,p%TwrGagNd(I)) + p%PtfmRefzt

      m%AllOuts( TwHtRPxi(I) ) =  m%RtHS%AngPosEF(1,p%TwrGagNd(I))*R2D
      m%AllOuts( TwHtRPyi(I) ) = -m%RtHS%AngPosEF(3,p%TwrGagNd(I))*R2D
      m%AllOuts( TwHtRPzi(I) ) =  m%RtHS%AngPosEF(2,p%TwrGagNd(I))*R2D

   END DO !I

      ! Platform Motions:

   m%AllOuts( PtfmTDxt) =  DOT_PRODUCT(       m%RtHS%rZ, m%CoordSys%a1 )
   m%AllOuts( PtfmTDyt) = -DOT_PRODUCT(       m%RtHS%rZ, m%CoordSys%a3 )
   m%AllOuts( PtfmTDzt) =  DOT_PRODUCT(       m%RtHS%rZ, m%CoordSys%a2 )
   m%AllOuts( PtfmTDxi) = x%QT  (DOF_Sg )
   m%AllOuts( PtfmTDyi) = x%QT  (DOF_Sw )
   m%AllOuts( PtfmTDzi) = x%QT  (DOF_Hv )
   m%AllOuts( PtfmTVxt) =  DOT_PRODUCT( m%RtHS%LinVelEZ, m%CoordSys%a1 )
   m%AllOuts( PtfmTVyt) = -DOT_PRODUCT( m%RtHS%LinVelEZ, m%CoordSys%a3 )
   m%AllOuts( PtfmTVzt) =  DOT_PRODUCT( m%RtHS%LinVelEZ, m%CoordSys%a2 )
   m%AllOuts( PtfmTVxi) = x%QDT (DOF_Sg )
   m%AllOuts( PtfmTVyi) = x%QDT (DOF_Sw )
   m%AllOuts( PtfmTVzi) = x%QDT (DOF_Hv )
   m%AllOuts( PtfmTAxt) =  DOT_PRODUCT(                 LinAccEZ, m%CoordSys%a1 )
   m%AllOuts( PtfmTAyt) = -DOT_PRODUCT(                 LinAccEZ, m%CoordSys%a3 )
   m%AllOuts( PtfmTAzt) =  DOT_PRODUCT(                 LinAccEZ, m%CoordSys%a2 )
   m%AllOuts( PtfmTAxi) = m%QD2T(DOF_Sg  )
   m%AllOuts( PtfmTAyi) = m%QD2T(DOF_Sw  )
   m%AllOuts( PtfmTAzi) = m%QD2T(DOF_Hv  )
   m%AllOuts( PtfmRDxi) = x%QT  (DOF_R )*R2D
   m%AllOuts( PtfmRDyi) = x%QT  (DOF_P )*R2D
   m%AllOuts( PtfmRDzi) = x%QT  (DOF_Y )*R2D
   m%AllOuts( PtfmRVxt) =  DOT_PRODUCT( m%RtHS%AngVelEX, m%CoordSys%a1 )*R2D
   m%AllOuts( PtfmRVyt) = -DOT_PRODUCT( m%RtHS%AngVelEX, m%CoordSys%a3 )*R2D
   m%AllOuts( PtfmRVzt) =  DOT_PRODUCT( m%RtHS%AngVelEX, m%CoordSys%a2 )*R2D
   m%AllOuts( PtfmRVxi) = x%QDT (DOF_R )*R2D
   m%AllOuts( PtfmRVyi) = x%QDT (DOF_P )*R2D
   m%AllOuts( PtfmRVzi) = x%QDT (DOF_Y )*R2D
   m%AllOuts( PtfmRAxt) =  DOT_PRODUCT(                 AngAccEX, m%CoordSys%a1 )*R2D
   m%AllOuts( PtfmRAyt) = -DOT_PRODUCT(                 AngAccEX, m%CoordSys%a3 )*R2D
   m%AllOuts( PtfmRAzt) =  DOT_PRODUCT(                 AngAccEX, m%CoordSys%a2 )*R2D
   m%AllOuts( PtfmRAxi) = m%QD2T(DOF_R )*R2D
   m%AllOuts( PtfmRAyi) = m%QD2T(DOF_P )*R2D
   m%AllOuts( PtfmRAzi) = m%QD2T(DOF_Y )*R2D



      ! Blade Root Loads:

   DO K=1,p%NumBl
      m%AllOuts( RootFxc(K) ) = DOT_PRODUCT( FrcS0B(:,K), m%CoordSys%i1(K,:) )
      m%AllOuts( RootFyc(K) ) = DOT_PRODUCT( FrcS0B(:,K), m%CoordSys%i2(K,:) )
      m%AllOuts( RootFzc(K) ) = DOT_PRODUCT( FrcS0B(:,K), m%CoordSys%i3(K,:) )
      m%AllOuts( RootFxb(K) ) = DOT_PRODUCT( FrcS0B(:,K), m%CoordSys%j1(K,:) )
      m%AllOuts( RootFyb(K) ) = DOT_PRODUCT( FrcS0B(:,K), m%CoordSys%j2(K,:) )
      m%AllOuts( RootMxc(K) ) = DOT_PRODUCT( MomH0B(:,K), m%CoordSys%i1(K,:) )
      m%AllOuts( RootMyc(K) ) = DOT_PRODUCT( MomH0B(:,K), m%CoordSys%i2(K,:) )
      m%AllOuts( RootMzc(K) ) = DOT_PRODUCT( MomH0B(:,K), m%CoordSys%i3(K,:) )
      m%AllOuts( RootMxb(K) ) = DOT_PRODUCT( MomH0B(:,K), m%CoordSys%j1(K,:) )
      m%AllOuts( RootMyb(K) ) = DOT_PRODUCT( MomH0B(:,K), m%CoordSys%j2(K,:) )
   END DO !K


      ! Blade Local Span Loads:

   DO K = 1,p%NumBl
      DO I = 1,p%NBlGages

      ! Initialize FrcMGagB and MomMGagB using the tip brake effects:

         FrcMGagB = m%RtHS%FSTipDrag(:,K) - p%TipMass(K)*( p%Gravity*m%CoordSys%z2 + LinAccES(:,p%TipNode,K) )
         MomMGagB = CROSS_PRODUCT( m%RtHS%rS0S(:,K,p%TipNode) - m%RtHS%rS0S(:,K,p%BldGagNd(I)), FrcMGagB )

      ! Integrate to find FrcMGagB and MomMGagB using all of the nodes / elements above the current strain gage location:
         DO J = ( p%BldGagNd(I) + 1 ),p%BldNodes ! Loop through blade nodes / elements above strain gage node

            TmpVec2  = m%RtHS%FSAero(:,K,J) - p%MassB(K,J)*( p%Gravity*m%CoordSys%z2 + LinAccES(:,J,K) )  ! Portion of FrcMGagB associated with element J
            FrcMGagB = FrcMGagB + TmpVec2*p%DRNodes(J)

            TmpVec = CROSS_PRODUCT( m%RtHS%rS0S(:,K,J) - m%RtHS%rS0S(:,K,p%BldGagNd(I)), TmpVec2 )           ! Portion of MomMGagB associated with element J
            MomMGagB = MomMGagB + ( TmpVec + m%RtHS%MMAero(:,K,J) )*p%DRNodes(J)

         ENDDO ! J - Blade nodes / elements above strain gage node

      ! Add the effects of 1/2 the strain gage element:
      ! NOTE: for the radius in this calculation, assume that there is no
      !   shortening effect (due to blade bending) within the element.  Thus,
      !   the moment arm for the force is 1/4 of p%DRNodes() and the element
      !   length is 1/2 of p%DRNodes().

         TmpVec2  = m%RtHS%FSAero(:,K,p%BldGagNd(I)) - p%MassB(K,p%BldGagNd(I))* ( p%Gravity*m%CoordSys%z2 + LinAccES(:,p%BldGagNd(I),K) ) ! Portion of FrcMGagB associated with 1/2 of the strain gage element
         FrcMGagB = FrcMGagB + TmpVec2 * 0.5 * p%DRNodes(p%BldGagNd(I))                                                    ! Portion of FrcMGagB associated with 1/2 of the strain gage element
         FrcMGagB = 0.001*FrcMGagB           ! Convert the local force to kN


         TmpVec = CROSS_PRODUCT( ( 0.25_R8Ki*p%DRNodes(p%BldGagNd(I)) )*m%CoordSys%j3(K,:), TmpVec2 )                              ! Portion of MomMGagB associated with 1/2 of the strain gage element

         MomMGagB = MomMGagB + ( TmpVec + m%RtHS%MMAero(:,K,p%BldGagNd(I)) )* ( 0.5 *p%DRNodes(p%BldGagNd(I)) )
         MomMGagB = 0.001*MomMGagB           ! Convert the local moment to kN-m


         m%AllOuts(SpnFLxb(I,K)) = DOT_PRODUCT( FrcMGagB, m%CoordSys%n1(K,p%BldGagNd(I),:) )
         m%AllOuts(SpnFLyb(I,K)) = DOT_PRODUCT( FrcMGagB, m%CoordSys%n2(K,p%BldGagNd(I),:) )
         m%AllOuts(SpnFLzb(I,K)) = DOT_PRODUCT( FrcMGagB, m%CoordSys%n3(K,p%BldGagNd(I),:) )

         m%AllOuts(SpnMLxb(I,K)) = DOT_PRODUCT( MomMGagB, m%CoordSys%n1(K,p%BldGagNd(I),:) )
         m%AllOuts(SpnMLyb(I,K)) = DOT_PRODUCT( MomMGagB, m%CoordSys%n2(K,p%BldGagNd(I),:) )
         m%AllOuts(SpnMLzb(I,K)) = DOT_PRODUCT( MomMGagB, m%CoordSys%n3(K,p%BldGagNd(I),:) )
      END DO ! I
   END DO ! K



      ! Hub and Rotor Loads:

   !ComDenom = 0.5*p%AirDens*p%ProjArea*m%AllOuts(  WindVxi)*m%AllOuts(  WindVxi)   ! Common denominator used in several expressions

   m%AllOuts(LSShftFxa) =  DOT_PRODUCT(  FrcPRot, m%CoordSys%e1 )
   m%AllOuts(LSShftFya) =  DOT_PRODUCT(  FrcPRot, m%CoordSys%e2 )
   m%AllOuts(LSShftFza) =  DOT_PRODUCT(  FrcPRot, m%CoordSys%e3 )
   m%AllOuts(LSShftFys) = -DOT_PRODUCT(  FrcPRot, m%CoordSys%c3 )
   m%AllOuts(LSShftFzs) =  DOT_PRODUCT(  FrcPRot, m%CoordSys%c2 )
   m%AllOuts(LSShftMxa) =  DOT_PRODUCT( MomLPRot, m%CoordSys%e1 )
   m%AllOuts(LSSTipMya) =  DOT_PRODUCT( MomLPRot, m%CoordSys%e2 )
   m%AllOuts(LSSTipMza) =  DOT_PRODUCT( MomLPRot, m%CoordSys%e3 )
   m%AllOuts(LSSTipMys) = -DOT_PRODUCT( MomLPRot, m%CoordSys%c3 )
   m%AllOuts(LSSTipMzs) =  DOT_PRODUCT( MomLPRot, m%CoordSys%c2 )

!   IF ( .NOT. EqualRealNos( m%AllOuts(LSShftFxa), 0.0_ReKi ) )   THEN ! .TRUE. if the denominator in the following equations is not zero.
!
!      CThrstys = -m%AllOuts(LSSTipMzs)/m%AllOuts(LSShftFxa)  ! Estimate of the ys-location of the center of thrust
!      CThrstzs =  m%AllOuts(LSSTipMys)/m%AllOuts(LSShftFxa)  ! Estimate of the zs-location of the center of thrust
!
!!      m%AllOuts(CThrstAzm) = MOD( ( ATAN2( -CThrstzs, -CThrstys ) + p%AzimB1Up )*R2D + 360.0 + 90.0, 360.0 )  !bjj: IgnoreMod was used for linearization... perhaps these outputs should not use the MOD function; only WriteOutputs should have that...
!!      m%AllOuts(CThrstRad) = MIN( 1.0, SQRT( CThrstys*CThrstys + CThrstzs*CThrstzs )/p%AvgNrmTpRd )
!
!   ELSE
!
!      !m%AllOuts(CThrstAzm) = 0.0
!      !m%AllOuts(CThrstRad) = 0.0
!
!   ENDIF

   m%AllOuts(   RotPwr) = ( x%QDT(DOF_GeAz) + x%QDT(DOF_DrTr) )*m%AllOuts(LSShftMxa)

   !IF ( .NOT. EqualRealNos( ComDenom, 0.0_ReKi ) )  THEN   ! .TRUE. if the denominator in the following equations is not zero.
   !
   !   m%AllOuts( RotCq) = 1000.0*m%AllOuts(LSShftMxa) / ( ComDenom*p%TipRad )
   !   m%AllOuts( RotCp) = 1000.0*m%AllOuts(   RotPwr) / ( ComDenom*m%AllOuts(  WindVxi) )
   !   m%AllOuts( RotCt) = 1000.0*m%AllOuts(LSShftFxa) /   ComDenom
   !
   !ELSE
   !
   !   m%AllOuts( RotCq) = 0.0
   !   m%AllOuts( RotCp) = 0.0
   !   m%AllOuts( RotCt) = 0.0
   !
   !ENDIF


      ! Shaft Strain Gage Loads:

   m%AllOuts(LSSGagMya) = m%AllOuts(LSSTipMya) + p%ShftGagL*m%AllOuts(LSShftFza)
   m%AllOuts(LSSGagMza) = m%AllOuts(LSSTipMza) - p%ShftGagL*m%AllOuts(LSShftFya)
   m%AllOuts(LSSGagMys) = m%AllOuts(LSSTipMys) + p%ShftGagL*m%AllOuts(LSShftFzs)
   m%AllOuts(LSSGagMzs) = m%AllOuts(LSSTipMzs) - p%ShftGagL*m%AllOuts(LSShftFys)


      ! Generator and High-Speed Shaft Loads:

   m%AllOuts( HSShftTq)  = m%AllOuts(LSShftMxa)*m%RtHS%GBoxEffFac/ABS(p%GBRatio)
   m%AllOuts(HSShftPwr)  = m%AllOuts( HSShftTq)*ABS(p%GBRatio)*x%QDT(DOF_GeAz)
   m%AllOuts(HSSBrTq)    = OtherState%HSSBrTrq*0.001_ReKi


   !IF ( .NOT. EqualRealNos( ComDenom, 0.0_ReKi ) )  THEN  ! .TRUE. if the denominator in the following equations is not zero (ComDenom is the same as it is calculated above).
   !
   !   m%AllOuts( HSShftCq) = 1000.0*m%AllOuts( HSShftTq) / ( ComDenom*p%TipRad )
   !   m%AllOuts( HSShftCp) = 1000.0*m%AllOuts(HSShftPwr) / ( ComDenom*m%AllOuts(  WindVxi) )
   !   m%AllOuts(    GenCq) = 1000.0*m%AllOuts(    GenTq) / ( ComDenom*p%TipRad )
   !   m%AllOuts(    GenCp) = 1000.0*m%AllOuts(   GenPwr) / ( ComDenom*m%AllOuts(  WindVxi) )
   !
   !ELSE
   !
   !   m%AllOuts( HSShftCq) = 0.0
   !   m%AllOuts( HSShftCp) = 0.0
   !   m%AllOuts(    GenCq) = 0.0
   !   m%AllOuts(    GenCp) = 0.0
   !
   !ENDIF


      ! Rotor-Furl Axis Loads:

   m%AllOuts(RFrlBrM  ) =  DOT_PRODUCT( MomNGnRt, m%CoordSys%rfa )


      ! Tail-Furl Axis Loads:

   m%AllOuts(TFrlBrM  ) =  DOT_PRODUCT( MomNTail, m%CoordSys%tfa )


      ! Tower-Top / Yaw Bearing Loads:

   m%AllOuts( YawBrFxn) =  DOT_PRODUCT( FrcONcRt, m%CoordSys%d1 )
   m%AllOuts( YawBrFyn) = -DOT_PRODUCT( FrcONcRt, m%CoordSys%d3 )
   m%AllOuts( YawBrFzn) =  DOT_PRODUCT( FrcONcRt, m%CoordSys%d2 )
   m%AllOuts( YawBrFxp) =  DOT_PRODUCT( FrcONcRt, m%CoordSys%b1 )
   m%AllOuts( YawBrFyp) = -DOT_PRODUCT( FrcONcRt, m%CoordSys%b3 )
   m%AllOuts( YawBrMxn) =  DOT_PRODUCT( MomBNcRt, m%CoordSys%d1 )
   m%AllOuts( YawBrMyn) = -DOT_PRODUCT( MomBNcRt, m%CoordSys%d3 )
   m%AllOuts( YawBrMzn) =  DOT_PRODUCT( MomBNcRt, m%CoordSys%d2 )
   m%AllOuts( YawBrMxp) =  DOT_PRODUCT( MomBNcRt, m%CoordSys%b1 )
   m%AllOuts( YawBrMyp) = -DOT_PRODUCT( MomBNcRt, m%CoordSys%b3 )


      ! Tower Base Loads:

   m%AllOuts( TwrBsFxt) =  DOT_PRODUCT( FrcT0Trb, m%CoordSys%a1 )
   m%AllOuts( TwrBsFyt) = -DOT_PRODUCT( FrcT0Trb, m%CoordSys%a3 )
   m%AllOuts( TwrBsFzt) =  DOT_PRODUCT( FrcT0Trb, m%CoordSys%a2 )
   m%AllOuts( TwrBsMxt) =  DOT_PRODUCT( MomX0Trb, m%CoordSys%a1 )
   m%AllOuts( TwrBsMyt) = -DOT_PRODUCT( MomX0Trb, m%CoordSys%a3 )
   m%AllOuts( TwrBsMzt) =  DOT_PRODUCT( MomX0Trb, m%CoordSys%a2 )


      ! Local Tower Loads:

   FrcONcRt = 1000.0*FrcONcRt ! Convert the units of these forces and moments
   MomBNcRt = 1000.0*MomBNcRt ! from kN and kN-m back to N and N-m, respectively.

   DO I=1,p%NTwGages

      ! Initialize FrcFGagT and MomFGagT using the tower-top and yaw bearing mass effects:
      FrcFGagT = FrcONcRt - p%YawBrMass*( p%Gravity*m%CoordSys%z2 + LinAccEO )
      MomFGagT = CROSS_PRODUCT( m%RtHS%rZO - m%RtHS%rZT(:,p%TwrGagNd(I)), FrcFGagT )
      MomFGagT = MomFGagT + MomBNcRt

      ! Integrate to find FrcFGagT and MomFGagT using all of the nodes / elements above the current strain gage location:
      DO J = ( p%TwrGagNd(I) + 1 ),p%TwrNodes ! Loop through tower nodes / elements above strain gage node
         TmpVec2  = FTTower(:,J) - p%MassT(J)*( p%Gravity*m%CoordSys%z2 + LinAccET(:,J) )           ! Portion of FrcFGagT associated with element J
         FrcFGagT = FrcFGagT + TmpVec2*abs(p%DHNodes(J))

         TmpVec = CROSS_PRODUCT( m%RtHS%rZT(:,J) - m%RtHS%rZT(:,p%TwrGagNd(I)), TmpVec2 )                          ! Portion of MomFGagT associated with element J
         MomFGagT = MomFGagT + ( TmpVec + MFHydro(:,J) )*abs(p%DHNodes(J))
      ENDDO ! J -Tower nodes / elements above strain gage node

      ! Add the effects of 1/2 the strain gage element:
      ! NOTE: for the radius in this calculation, assume that there is no shortening
      !   effect (due to tower bending) within the element.  Thus, the moment arm
      !   for the force is 1/4 of DHNodes() and the element length is 1/2 of DHNodes().

      TmpVec2  = FTTower(:,p%TwrGagNd(I)) - p%MassT(p%TwrGagNd(I))*( p%Gravity*m%CoordSys%z2 + LinAccET(:,p%TwrGagNd(I)))

      FrcFGagT = FrcFGagT + TmpVec2 * 0.5 * abs(p%DHNodes(p%TwrGagNd(I)))
      FrcFGagT = 0.001*FrcFGagT  ! Convert the local force to kN

      TmpVec = CROSS_PRODUCT( ( 0.25_R8Ki*p%DHNodes( p%TwrGagNd(I)) )*m%CoordSys%a2, TmpVec2 )              ! Portion of MomFGagT associated with 1/2 of the strain gage element
      TmpVec   = TmpVec   + MFHydro(:,p%TwrGagNd(I))
      MomFGagT = MomFGagT + TmpVec * 0.5 * abs(p%DHNodes(p%TwrGagNd(I)))
      MomFGagT = 0.001*MomFGagT  ! Convert the local moment to kN-m

      m%AllOuts( TwHtFLxt(I) ) =     DOT_PRODUCT( FrcFGagT, m%CoordSys%t1(p%TwrGagNd(I),:) )
      m%AllOuts( TwHtFLyt(I) ) = -1.*DOT_PRODUCT( FrcFGagT, m%CoordSys%t3(p%TwrGagNd(I),:) )
      m%AllOuts( TwHtFLzt(I) ) =     DOT_PRODUCT( FrcFGagT, m%CoordSys%t2(p%TwrGagNd(I),:) )

      m%AllOuts( TwHtMLxt(I) ) =     DOT_PRODUCT( MomFGagT, m%CoordSys%t1(p%TwrGagNd(I),:) )
      m%AllOuts( TwHtMLyt(I) ) = -1.*DOT_PRODUCT( MomFGagT, m%CoordSys%t3(p%TwrGagNd(I),:) )
      m%AllOuts( TwHtMLzt(I) ) =     DOT_PRODUCT( MomFGagT, m%CoordSys%t2(p%TwrGagNd(I),:) )

   END DO


   !   ! Platform Loads:
   !
   !m%AllOuts(  PtfmFxt) =  DOT_PRODUCT( FZHydro, m%CoordSys%a1 )
   !m%AllOuts(  PtfmFyt) = -DOT_PRODUCT( FZHydro, m%CoordSys%a3 )
   !m%AllOuts(  PtfmFzt) =  DOT_PRODUCT( FZHydro, m%CoordSys%a2 )
   !m%AllOuts(  PtfmFxi) =  DOT_PRODUCT( FZHydro, m%CoordSys%z1 )
   !m%AllOuts(  PtfmFyi) = -DOT_PRODUCT( FZHydro, m%CoordSys%z3 )
   !m%AllOuts(  PtfmFzi) =  DOT_PRODUCT( FZHydro, m%CoordSys%z2 )
   !m%AllOuts(  PtfmMxt) =  DOT_PRODUCT( MXHydro, m%CoordSys%a1 )
   !m%AllOuts(  PtfmMyt) = -DOT_PRODUCT( MXHydro, m%CoordSys%a3 )
   !m%AllOuts(  PtfmMzt) =  DOT_PRODUCT( MXHydro, m%CoordSys%a2 )
   !m%AllOuts(  PtfmMxi) =  DOT_PRODUCT( MXHydro, m%CoordSys%z1 )
   !m%AllOuts(  PtfmMyi) = -DOT_PRODUCT( MXHydro, m%CoordSys%z3 )
   !m%AllOuts(  PtfmMzi) =  DOT_PRODUCT( MXHydro, m%CoordSys%z2 )
   !

      ! Internal p%DOFs outputs:

   m%AllOuts( Q_B1E1   ) = x%QT(   DOF_BE(1,1) )
   m%AllOuts( Q_B1F1   ) = x%QT(   DOF_BF(1,1) )
   m%AllOuts( Q_B1F2   ) = x%QT(   DOF_BF(1,2) )
   m%AllOuts( Q_DrTr   ) = x%QT(   DOF_DrTr    )
   m%AllOuts( Q_GeAz   ) = x%QT(   DOF_GeAz    )
   m%AllOuts( Q_RFrl   ) = x%QT(   DOF_RFrl    )
   m%AllOuts( Q_TFrl   ) = x%QT(   DOF_TFrl    )
   m%AllOuts( Q_Yaw    ) = x%QT(   DOF_Yaw     )
   m%AllOuts( Q_TFA1   ) = x%QT(   DOF_TFA1    )
   m%AllOuts( Q_TSS1   ) = x%QT(   DOF_TSS1    )
   m%AllOuts( Q_TFA2   ) = x%QT(   DOF_TFA2    )
   m%AllOuts( Q_TSS2   ) = x%QT(   DOF_TSS2    )
   m%AllOuts( Q_Sg     ) = x%QT(   DOF_Sg      )
   m%AllOuts( Q_Sw     ) = x%QT(   DOF_Sw      )
   m%AllOuts( Q_Hv     ) = x%QT(   DOF_Hv      )
   m%AllOuts( Q_R      ) = x%QT(   DOF_R       )
   m%AllOuts( Q_P      ) = x%QT(   DOF_P       )
   m%AllOuts( Q_Y      ) = x%QT(   DOF_Y       )

   m%AllOuts( QD_B1E1  ) = x%QDT(  DOF_BE(1,1) )
   m%AllOuts( QD_B1F1  ) = x%QDT(  DOF_BF(1,1) )
   m%AllOuts( QD_B1F2  ) = x%QDT(  DOF_BF(1,2) )
   m%AllOuts( QD_DrTr  ) = x%QDT(  DOF_DrTr    )
   m%AllOuts( QD_GeAz  ) = x%QDT(  DOF_GeAz    )
   m%AllOuts( QD_RFrl  ) = x%QDT(  DOF_RFrl    )
   m%AllOuts( QD_TFrl  ) = x%QDT(  DOF_TFrl    )
   m%AllOuts( QD_Yaw   ) = x%QDT(  DOF_Yaw     )
   m%AllOuts( QD_TFA1  ) = x%QDT(  DOF_TFA1    )
   m%AllOuts( QD_TSS1  ) = x%QDT(  DOF_TSS1    )
   m%AllOuts( QD_TFA2  ) = x%QDT(  DOF_TFA2    )
   m%AllOuts( QD_TSS2  ) = x%QDT(  DOF_TSS2    )
   m%AllOuts( QD_Sg    ) = x%QDT(  DOF_Sg      )
   m%AllOuts( QD_Sw    ) = x%QDT(  DOF_Sw      )
   m%AllOuts( QD_Hv    ) = x%QDT(  DOF_Hv      )
   m%AllOuts( QD_R     ) = x%QDT(  DOF_R       )
   m%AllOuts( QD_P     ) = x%QDT(  DOF_P       )
   m%AllOuts( QD_Y     ) = x%QDT(  DOF_Y       )

   m%AllOuts( QD2_B1E1 ) = m%QD2T( DOF_BE(1,1) )
   m%AllOuts( QD2_B1F1 ) = m%QD2T( DOF_BF(1,1) )
   m%AllOuts( QD2_B1F2 ) = m%QD2T( DOF_BF(1,2) )
   m%AllOuts( QD2_DrTr ) = m%QD2T( DOF_DrTr    )
   m%AllOuts( QD2_GeAz ) = m%QD2T( DOF_GeAz    )
   m%AllOuts( QD2_RFrl ) = m%QD2T( DOF_RFrl    )
   m%AllOuts( QD2_TFrl ) = m%QD2T( DOF_TFrl    )
   m%AllOuts( QD2_Yaw  ) = m%QD2T( DOF_Yaw     )
   m%AllOuts( QD2_TFA1 ) = m%QD2T( DOF_TFA1    )
   m%AllOuts( QD2_TSS1 ) = m%QD2T( DOF_TSS1    )
   m%AllOuts( QD2_TFA2 ) = m%QD2T( DOF_TFA2    )
   m%AllOuts( QD2_TSS2 ) = m%QD2T( DOF_TSS2    )
   m%AllOuts( QD2_Sg   ) = m%QD2T( DOF_Sg      )
   m%AllOuts( QD2_Sw   ) = m%QD2T( DOF_Sw      )
   m%AllOuts( QD2_Hv   ) = m%QD2T( DOF_Hv      )
   m%AllOuts( QD2_R    ) = m%QD2T( DOF_R       )
   m%AllOuts( QD2_P    ) = m%QD2T( DOF_P       )
   m%AllOuts( QD2_Y    ) = m%QD2T( DOF_Y       )

IF ( p%NumBl > 1 ) THEN

   m%AllOuts( Q_B2E1   ) = x%QT(   DOF_BE(2,1) )
   m%AllOuts( Q_B2F1   ) = x%QT(   DOF_BF(2,1) )
   m%AllOuts( Q_B2F2   ) = x%QT(   DOF_BF(2,2) )
      
   m%AllOuts( QD_B2E1  ) = x%QDT(  DOF_BE(2,1) )
   m%AllOuts( QD_B2F1  ) = x%QDT(  DOF_BF(2,1) )
   m%AllOuts( QD_B2F2  ) = x%QDT(  DOF_BF(2,2) )

   m%AllOuts( QD2_B2E1 ) = m%QD2T( DOF_BE(2,1) )
   m%AllOuts( QD2_B2F1 ) = m%QD2T( DOF_BF(2,1) )
   m%AllOuts( QD2_B2F2 ) = m%QD2T( DOF_BF(2,2) )
   
   IF ( p%NumBl > 2 ) THEN
      m%AllOuts( Q_B3E1   ) = x%QT(   DOF_BE(3,1) )
      m%AllOuts( Q_B3F1   ) = x%QT(   DOF_BF(3,1) )
      m%AllOuts( Q_B3F2   ) = x%QT(   DOF_BF(3,2) )

      m%AllOuts( QD_B3E1  ) = x%QDT(  DOF_BE(3,1) )
      m%AllOuts( QD_B3F1  ) = x%QDT(  DOF_BF(3,1) )
      m%AllOuts( QD_B3F2  ) = x%QDT(  DOF_BF(3,2) )

      m%AllOuts( QD2_B3E1 ) = m%QD2T( DOF_BE(3,1) )
      m%AllOuts( QD2_B3F1 ) = m%QD2T( DOF_BF(3,1) )
      m%AllOuts( QD2_B3F2 ) = m%QD2T( DOF_BF(3,2) )
   ELSE
      m%AllOuts( Q_Teet   ) = x%QT(   DOF_Teet    )
      m%AllOuts( QD_Teet  ) = x%QDT(  DOF_Teet    )
      m%AllOuts( QD2_Teet ) = m%QD2T( DOF_Teet    )
   END IF
      
END IF

   !...............................................................................................................................
   ! Place the selected output channels into the WriteOutput(:) array with the proper sign:
   !...............................................................................................................................

   DO I = 1,p%NumOuts  ! Loop through all selected output channels

      y%WriteOutput(I) = p%OutParam(I)%SignM * m%AllOuts( p%OutParam(I)%Indx )

   ENDDO             ! I - All selected output channels

   IF ( .NOT. p%BD4Blades ) THEN
      y%WriteOutput(p%NumOuts+1:) = 0.0_ReKi

         ! Now we need to populate the blade node outputs here
      call Calc_WriteAllBldNdOutput( p, u, m, y, LinAccES, ErrStat2, ErrMsg2 )   ! Call after normal writeoutput.  Will just postpend data on here.
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ED_CalcOutput')
   ENDIF


   !...............................................................................................................................
   ! Outputs required for AeroDyn
   !...............................................................................................................................
   
   !JASON: WE SHOULD REALLY BE PASSING TO AERODYN THE LINEAR VELOCITIES OF THE AERODYNAMIC CENTER IN THE INERTIA FRAME, NOT SIMPLY THE LINEAR VELOCITIES OF POINT S.  
   !       IS THERE ANY WAY OF GETTING THIS VELOCITY?<--DO THIS, WHEN YOU ADD THE COUPLED MODE SHAPES!!!!
   
   !...........
   ! Blade elements:
   !...........
   IF ( ALLOCATED(y%BladeLn2Mesh) ) THEN
      DO K = 1,p%NumBl ! Loop through all blades
         DO J = 0,p%TipNode ! Loop through the blade nodes / elements
            
            J2 = J            
            if (j==0) then
                  ! blade root
               NodeNum = p%BldNodes + 2
               if (p%UseAD14) j2 = 1                  
            elseif (j==p%TipNode) then
               ! blade tip
               NodeNum = p%BldNodes + 1
               if (p%UseAD14) j2 = p%BldNodes
            else
               NodeNum = J
            end if
                                                                                                        
            if (p%UseAD14) then                  
                  ! Translational Displacement (first calculate absolute position)
               y%BladeLn2Mesh(K)%TranslationDisp(1,NodeNum) =     m%RtHS%rS (1,K,J2) + m%RtHS%rSAerCen(1,J2,K)               ! = the distance from the undeflected tower centerline                                     to the current blade aerodynamic center in the xi ( z1) direction
               y%BladeLn2Mesh(K)%TranslationDisp(2,NodeNum) = -1.*m%RtHS%rS (3,K,J2) - m%RtHS%rSAerCen(3,J2,K)               ! = the distance from the undeflected tower centerline                                     to the current blade aerodynamic center in the yi (-z3) direction
               y%BladeLn2Mesh(K)%TranslationDisp(3,NodeNum) =     m%RtHS%rS (2,K,J2) + m%RtHS%rSAerCen(2,J2,K) + p%PtfmRefzt ! = the distance from the nominal tower base position (i.e., the undeflected position of the tower base) to the current blade aerodynamic center in the zi ( z2) direction
               
                  ! Orientation 
               y%BladeLn2Mesh(K)%Orientation(1,1,NodeNum) =     m%CoordSys%te1(K,J2,1)
               y%BladeLn2Mesh(K)%Orientation(2,1,NodeNum) =     m%CoordSys%te2(K,J2,1)
               y%BladeLn2Mesh(K)%Orientation(3,1,NodeNum) =     m%CoordSys%te3(K,J2,1)
               y%BladeLn2Mesh(K)%Orientation(1,2,NodeNum) = -1.*m%CoordSys%te1(K,J2,3)
               y%BladeLn2Mesh(K)%Orientation(2,2,NodeNum) = -1.*m%CoordSys%te2(K,J2,3)
               y%BladeLn2Mesh(K)%Orientation(3,2,NodeNum) = -1.*m%CoordSys%te3(K,J2,3)
               y%BladeLn2Mesh(K)%Orientation(1,3,NodeNum) =     m%CoordSys%te1(K,J2,2)
               y%BladeLn2Mesh(K)%Orientation(2,3,NodeNum) =     m%CoordSys%te2(K,J2,2)
               y%BladeLn2Mesh(K)%Orientation(3,3,NodeNum) =     m%CoordSys%te3(K,J2,2)
               
            else         
                  ! Translational Displacement (first calculate absolute position)
               y%BladeLn2Mesh(K)%TranslationDisp(1,NodeNum) =     m%RtHS%rS (1,K,J2)                ! = the distance from the undeflected tower centerline to the current blade node in the xi ( z1) direction
               y%BladeLn2Mesh(K)%TranslationDisp(2,NodeNum) = -1.*m%RtHS%rS (3,K,J2)                ! = the distance from the undeflected tower centerline to the current blade node in the yi (-z3) direction
               y%BladeLn2Mesh(K)%TranslationDisp(3,NodeNum) =     m%RtHS%rS (2,K,J2)  + p%PtfmRefzt ! = the distance from the nominal tower base position (i.e., the undeflected position of the tower base) to the current blade node in the zi ( z2) direction
               
                  ! Orientation
               y%BladeLn2Mesh(K)%Orientation(1,1,NodeNum) =     m%CoordSys%n1(K,J2,1)
               y%BladeLn2Mesh(K)%Orientation(2,1,NodeNum) =     m%CoordSys%n2(K,J2,1)
               y%BladeLn2Mesh(K)%Orientation(3,1,NodeNum) =     m%CoordSys%n3(K,J2,1)
               y%BladeLn2Mesh(K)%Orientation(1,2,NodeNum) = -1.*m%CoordSys%n1(K,J2,3)
               y%BladeLn2Mesh(K)%Orientation(2,2,NodeNum) = -1.*m%CoordSys%n2(K,J2,3)
               y%BladeLn2Mesh(K)%Orientation(3,2,NodeNum) = -1.*m%CoordSys%n3(K,J2,3)
               y%BladeLn2Mesh(K)%Orientation(1,3,NodeNum) =     m%CoordSys%n1(K,J2,2)
               y%BladeLn2Mesh(K)%Orientation(2,3,NodeNum) =     m%CoordSys%n2(K,J2,2)
               y%BladeLn2Mesh(K)%Orientation(3,3,NodeNum) =     m%CoordSys%n3(K,J2,2)
            end if
            
               ! Translational Displacement (get displacement, not absolute position):
            y%BladeLn2Mesh(K)%TranslationDisp(:,NodeNum) = y%BladeLn2Mesh(K)%TranslationDisp(:,NodeNum) - y%BladeLn2Mesh(K)%Position(:,NodeNum)
            
               ! Translational Velocity
            y%BladeLn2Mesh(K)%TranslationVel(1,NodeNum) =     m%RtHS%LinVelES(1,J2,K)
            y%BladeLn2Mesh(K)%TranslationVel(2,NodeNum) = -1.*m%RtHS%LinVelES(3,J2,K)
            y%BladeLn2Mesh(K)%TranslationVel(3,NodeNum) =     m%RtHS%LinVelES(2,J2,K)  
            
               ! Rotational Velocity
            y%BladeLn2Mesh(K)%RotationVel(1,NodeNum) =     m%RtHS%AngVelEM(1,J2,K)
            y%BladeLn2Mesh(K)%RotationVel(2,NodeNum) = -1.*m%RtHS%AngVelEM(3,J2,K)
            y%BladeLn2Mesh(K)%RotationVel(3,NodeNum) =     m%RtHS%AngVelEM(2,J2,K)  

               ! Translational Acceleration
            y%BladeLn2Mesh(K)%TranslationAcc(1,NodeNum) =     LinAccES(1,J2,K)
            y%BladeLn2Mesh(K)%TranslationAcc(2,NodeNum) = -1.*LinAccES(3,J2,K)
            y%BladeLn2Mesh(K)%TranslationAcc(3,NodeNum) =     LinAccES(2,J2,K)  

               ! Rotational Acceleration
            y%BladeLn2Mesh(K)%RotationAcc(1,NodeNum)     =     AngAccEK(1,J2,K)
            y%BladeLn2Mesh(K)%RotationAcc(2,NodeNum)     = -1.*AngAccEK(3,J2,K)
            y%BladeLn2Mesh(K)%RotationAcc(3,NodeNum)     =     AngAccEK(2,J2,K) 
               
            
         END DO !J = 1,p%BldNodes ! Loop through the blade nodes / elements
                  
      END DO !K = 1,p%NumBl
   END IF
   
      
   !...........
   ! Hub (for Lidar and AeroDyn15):
   !...........   
   
         ! Translation (absolute position - starting position):
   y%HubPtMotion%TranslationDisp(1,1)  =     m%RtHS%rQ(1)
   y%HubPtMotion%TranslationDisp(2,1)  = -1.*m%RtHS%rQ(3)
   y%HubPtMotion%TranslationDisp(3,1)  =     m%RtHS%rQ(2) + p%PtfmRefzt
   y%HubPtMotion%TranslationDisp       = y%HubPtMotion%TranslationDisp - y%HubPtMotion%Position   ! relative position
   
      ! Orientation:        
   y%HubPtMotion%Orientation(1,1,1)    =     m%CoordSys%g1(1) 
   y%HubPtMotion%Orientation(2,1,1)    =     m%CoordSys%g2(1)
   y%HubPtMotion%Orientation(3,1,1)    =     m%CoordSys%g3(1)   
   y%HubPtMotion%Orientation(1,2,1)    = -1.*m%CoordSys%g1(3)
   y%HubPtMotion%Orientation(2,2,1)    = -1.*m%CoordSys%g2(3) 
   y%HubPtMotion%Orientation(3,2,1)    = -1.*m%CoordSys%g3(3) 
   y%HubPtMotion%Orientation(1,3,1)    =     m%CoordSys%g1(2)
   y%HubPtMotion%Orientation(2,3,1)    =     m%CoordSys%g2(2)
   y%HubPtMotion%Orientation(3,3,1)    =     m%CoordSys%g3(2)
   
      ! Rotational velocity:
   y%HubPtMotion%RotationVel(1,1)      =     m%RtHS%AngVelEH(1)
   y%HubPtMotion%RotationVel(2,1)      = -1.*m%RtHS%AngVelEH(3)
   y%HubPtMotion%RotationVel(3,1)      =     m%RtHS%AngVelEH(2)   
   
   !...........
   ! Blade roots (BeamDyn/AeroDyn v15):
   !...........   
         
   DO K=1,p%NumBl
         
      ! Translation displacement  ! rS at the root      
      y%BladeRootMotion(K)%TranslationDisp(1,1) =            m%RtHS%rS (1,K,0)                ! = the distance from the undeflected tower centerline to the current blade node in the xi ( z1) direction
      y%BladeRootMotion(K)%TranslationDisp(2,1) =        -1.*m%RtHS%rS (3,K,0)                ! = the distance from the undeflected tower centerline to the current blade node in the yi (-z3) direction
      y%BladeRootMotion(K)%TranslationDisp(3,1) =            m%RtHS%rS (2,K,0)  + p%PtfmRefzt ! = the distance from the nominal tower base position (i.e., the undeflected position of the tower base) to the current blade node in the zi ( z2) direction
      y%BladeRootMotion(K)%TranslationDisp      = y%BladeRootMotion(K)%TranslationDisp - y%BladeRootMotion(K)%Position ! make it relative
      
      
      ! Orientation 
      y%BladeRootMotion(K)%Orientation(1,1,1)   =     m%CoordSys%j1(K,1)
      y%BladeRootMotion(K)%Orientation(2,1,1)   =     m%CoordSys%j2(K,1)
      y%BladeRootMotion(K)%Orientation(3,1,1)   =     m%CoordSys%j3(K,1)
      y%BladeRootMotion(K)%Orientation(1,2,1)   = -1.*m%CoordSys%j1(K,3)
      y%BladeRootMotion(K)%Orientation(2,2,1)   = -1.*m%CoordSys%j2(K,3)
      y%BladeRootMotion(K)%Orientation(3,2,1)   = -1.*m%CoordSys%j3(K,3)
      y%BladeRootMotion(K)%Orientation(1,3,1)   =     m%CoordSys%j1(K,2)
      y%BladeRootMotion(K)%Orientation(2,3,1)   =     m%CoordSys%j2(K,2)
      y%BladeRootMotion(K)%Orientation(3,3,1)   =     m%CoordSys%j3(K,2)

      ! Translation velocity 
      y%BladeRootMotion(K)%TranslationVel(1,1)  =     m%RtHS%LinVelES(1,0,K)
      y%BladeRootMotion(K)%TranslationVel(2,1)  = -1.*m%RtHS%LinVelES(3,0,K)
      y%BladeRootMotion(K)%TranslationVel(3,1)  =     m%RtHS%LinVelES(2,0,K)

      ! Rotation velocity  
      y%BladeRootMotion(K)%RotationVel(1,1)     =      m%RtHS%AngVelEH(1)
      y%BladeRootMotion(K)%RotationVel(2,1)     =  -1.*m%RtHS%AngVelEH(3)
      y%BladeRootMotion(K)%RotationVel(3,1)     =      m%RtHS%AngVelEH(2)
      
      ! Translation acceleration
      y%BladeRootMotion(K)%TranslationAcc(1,1)  =      LinAccES(1,0,K)
      y%BladeRootMotion(K)%TranslationAcc(2,1)  =  -1.*LinAccES(3,0,K)
      y%BladeRootMotion(K)%TranslationAcc(3,1)  =      LinAccES(2,0,K)
      
      ! Rotation acceleration  
      y%BladeRootMotion(K)%RotationAcc(1,1)     =      AngAccEH(1) 
      y%BladeRootMotion(K)%RotationAcc(2,1)     =  -1.*AngAccEH(3) 
      y%BladeRootMotion(K)%RotationAcc(3,1)     =      AngAccEH(2)
      
   END DO   
   
   !...........
   ! Hub (for AeroDyn v14):
   !...........   

      ! the hub position should use rQ instead of rP, but AeroDyn 14 treats
      ! teeter deflections like blade deflections:
   
   y%HubPtMotion14%TranslationDisp(1,1)  =     m%RtHS%rP(1)
   y%HubPtMotion14%TranslationDisp(2,1)  = -1.*m%RtHS%rP(3)
   y%HubPtMotion14%TranslationDisp(3,1)  =     m%RtHS%rP(2) + p%PtfmRefzt
   
   y%HubPtMotion14%TranslationDisp  = y%HubPtMotion14%TranslationDisp - y%HubPtMotion14%Position   
   
        ! Hub orientation should use the g instead of e system, but the current version
        ! of AeroDyn calculates forces normal and tangential to the cone of rotation
         
   y%HubPtMotion14%Orientation(1,1,1)  =     m%CoordSys%e1(1) 
   y%HubPtMotion14%Orientation(2,1,1)  =     m%CoordSys%e2(1)
   y%HubPtMotion14%Orientation(3,1,1)  =     m%CoordSys%e3(1)   
   y%HubPtMotion14%Orientation(1,2,1)  = -1.*m%CoordSys%e1(3)
   y%HubPtMotion14%Orientation(2,2,1)  = -1.*m%CoordSys%e2(3) 
   y%HubPtMotion14%Orientation(3,2,1)  = -1.*m%CoordSys%e3(3) 
   y%HubPtMotion14%Orientation(1,3,1)  =     m%CoordSys%e1(2)
   y%HubPtMotion14%Orientation(2,3,1)  =     m%CoordSys%e2(2)
   y%HubPtMotion14%Orientation(3,3,1)  =     m%CoordSys%e3(2)
   
      ! Note the hub rotational velocity should be AngVelEH instead AngVelEL, but AeroDyn (13.00.00)
      ! treats teeter deflections like blade deflections:
   
   y%HubPtMotion14%RotationVel(1,1) =     m%RtHS%AngVelEL(1)
   y%HubPtMotion14%RotationVel(2,1) = -1.*m%RtHS%AngVelEL(3)
   y%HubPtMotion14%RotationVel(3,1) =     m%RtHS%AngVelEL(2)
      
   !...........
   ! Blade roots (AeroDyn v14):
   !...........   
   
        ! Blade root orientations should use the j instead of i system, but the current version
        ! of AeroDyn calculates forces normal and tangential to the cone of rotation
      
   DO K=1,p%NumBl
         
      y%BladeRootMotion14%Orientation(1,1,K) =     m%CoordSys%i1(K,1)
      y%BladeRootMotion14%Orientation(2,1,K) =     m%CoordSys%i2(K,1)
      y%BladeRootMotion14%Orientation(3,1,K) =     m%CoordSys%i3(K,1)
      y%BladeRootMotion14%Orientation(1,2,K) = -1.*m%CoordSys%i1(K,3)
      y%BladeRootMotion14%Orientation(2,2,K) = -1.*m%CoordSys%i2(K,3)
      y%BladeRootMotion14%Orientation(3,2,K) = -1.*m%CoordSys%i3(K,3)
      y%BladeRootMotion14%Orientation(1,3,K) =     m%CoordSys%i1(K,2)
      y%BladeRootMotion14%Orientation(2,3,K) =     m%CoordSys%i2(K,2)
      y%BladeRootMotion14%Orientation(3,3,K) =     m%CoordSys%i3(K,2)
            
   END DO
   
    
   !...........
   ! Rotor furl:
   !...........   
   
      ! Rotor furl position should be rP instead of rV, but AeroDyn needs this for the HubVDue2Yaw calculation:
   
   y%RotorFurlMotion14%TranslationDisp(1,1) =     m%RtHS%rV(1)
   y%RotorFurlMotion14%TranslationDisp(2,1) = -1.*m%RtHS%rV(3)
   y%RotorFurlMotion14%TranslationDisp(3,1) =     m%RtHS%rV(2) + p%PtfmRefzt
   
   y%RotorFurlMotion14%TranslationDisp      = y%RotorFurlMotion14%TranslationDisp - y%RotorFurlMotion14%Position   
         
        ! Rotor furl orientation (note the different order than hub and blade root!)
   
   y%RotorFurlMotion14%Orientation(1,1,1) =     m%CoordSys%c1(1)
   y%RotorFurlMotion14%Orientation(2,1,1) = -1.*m%CoordSys%c3(1)
   y%RotorFurlMotion14%Orientation(3,1,1) =     m%CoordSys%c2(1)
   y%RotorFurlMotion14%Orientation(1,2,1) = -1.*m%CoordSys%c1(3)
   y%RotorFurlMotion14%Orientation(2,2,1) =     m%CoordSys%c3(3)
   y%RotorFurlMotion14%Orientation(3,2,1) = -1.*m%CoordSys%c2(3)
   y%RotorFurlMotion14%Orientation(1,3,1) =     m%CoordSys%c1(2)
   y%RotorFurlMotion14%Orientation(2,3,1) = -1.*m%CoordSys%c3(2)
   y%RotorFurlMotion14%Orientation(3,3,1) =     m%CoordSys%c2(2) 
   
      ! rotaional velocity:
   y%RotorFurlMotion14%RotationVel(1,1) =     m%RtHS%AngVelER(1)
   y%RotorFurlMotion14%RotationVel(2,1) = -1.*m%RtHS%AngVelER(3)
   y%RotorFurlMotion14%RotationVel(3,1) =     m%RtHS%AngVelER(2)

   !...........
   ! TailFin :
   !...........   
   ! Translation (absolute position - starting position):
   y%TFinCMMotion%TranslationDisp(1,1) =     m%RtHS%rJ(1)
   y%TFinCMMotion%TranslationDisp(2,1) = -1.*m%RtHS%rJ(3)
   y%TFinCMMotion%TranslationDisp(3,1) =     m%RtHS%rJ(2) + p%PtfmRefzt
   y%TFinCMMotion%TranslationDisp      = y%TFinCMMotion%TranslationDisp - y%TFinCMMotion%Position
   ! Orientation:        
   y%TFinCMMotion%Orientation(1,1,1)   =     m%CoordSys%tf1(1) 
   y%TFinCMMotion%Orientation(2,1,1)   =     m%CoordSys%tf2(1)
   y%TFinCMMotion%Orientation(3,1,1)   =     m%CoordSys%tf3(1)   
   y%TFinCMMotion%Orientation(1,2,1)   = -1.*m%CoordSys%tf1(3)
   y%TFinCMMotion%Orientation(2,2,1)   = -1.*m%CoordSys%tf2(3) 
   y%TFinCMMotion%Orientation(3,2,1)   = -1.*m%CoordSys%tf3(3) 
   y%TFinCMMotion%Orientation(1,3,1)   =     m%CoordSys%tf1(2)
   y%TFinCMMotion%Orientation(2,3,1)   =     m%CoordSys%tf2(2)
   y%TFinCMMotion%Orientation(3,3,1)   =     m%CoordSys%tf3(2)
   ! Rotational velocity:
   y%TFinCMMotion%RotationVel(1,1)     =     m%RtHS%AngVelEA(1)
   y%TFinCMMotion%RotationVel(2,1)     = -1.*m%RtHS%AngVelEA(3)
   y%TFinCMMotion%RotationVel(3,1)     =     m%RtHS%AngVelEA(2)   
   ! Linear velocity:
   y%TFinCMMotion%TranslationVel(1,1)  =     m%RtHS%LinVelEJ(1)
   y%TFinCMMotion%TranslationVel(2,1)  = -1.*m%RtHS%LinVelEJ(3)
   y%TFinCMMotion%TranslationVel(3,1)  =     m%RtHS%LinVelEJ(2)


      
   !...........
   ! Nacelle :
   !...........   
      
   y%NacelleMotion%TranslationDisp(1,1) =     m%RtHS%rO(1)
   y%NacelleMotion%TranslationDisp(2,1) = -1.*m%RtHS%rO(3)
   y%NacelleMotion%TranslationDisp(3,1) =     m%RtHS%rO(2) + p%PtfmRefzt
               
   y%NacelleMotion%TranslationDisp      = y%NacelleMotion%TranslationDisp - y%NacelleMotion%Position   
   
      ! Nacelle orientation (note the different order than hub and blade root!)
   
   y%NacelleMotion%Orientation(1,1,1) =     m%CoordSys%d1(1)
   y%NacelleMotion%Orientation(2,1,1) = -1.*m%CoordSys%d3(1)
   y%NacelleMotion%Orientation(3,1,1) =     m%CoordSys%d2(1)
   y%NacelleMotion%Orientation(1,2,1) = -1.*m%CoordSys%d1(3)
   y%NacelleMotion%Orientation(2,2,1) =     m%CoordSys%d3(3)
   y%NacelleMotion%Orientation(3,2,1) = -1.*m%CoordSys%d2(3)
   y%NacelleMotion%Orientation(1,3,1) =     m%CoordSys%d1(2)
   y%NacelleMotion%Orientation(2,3,1) = -1.*m%CoordSys%d3(2)
   y%NacelleMotion%Orientation(3,3,1) =     m%CoordSys%d2(2) 
   
   y%NacelleMotion%RotationVel(1,1)   =     m%RtHS%AngVelEN(1)
   y%NacelleMotion%RotationVel(2,1)   = -1.*m%RtHS%AngVelEN(3)
   y%NacelleMotion%RotationVel(3,1)   =     m%RtHS%AngVelEN(2) 
      
   y%NacelleMotion%TranslationVel(1,1)  =     m%RtHS%LinVelEO(1)
   y%NacelleMotion%TranslationVel(2,1)  = -1.*m%RtHS%LinVelEO(3)
   y%NacelleMotion%TranslationVel(3,1)  =     m%RtHS%LinVelEO(2)
      
   y%NacelleMotion%RotationAcc(   1,1)  =      AngAccEN(1) 
   y%NacelleMotion%RotationAcc(   2,1)  =  -1.*AngAccEN(3) 
   y%NacelleMotion%RotationAcc(   3,1)  =      AngAccEN(2)
   
   y%NacelleMotion%TranslationAcc(1,1)  =      LinAccEO(1)
   y%NacelleMotion%TranslationAcc(2,1)  =  -1.*LinAccEO(3)
   y%NacelleMotion%TranslationAcc(3,1)  =      LinAccEO(2)
   
   
   !...........
   ! Tower :
   !...........         
   
      ! Tower base position should be rT(0) instead of rZ, but AeroDyn needs this for  the HubVDue2Yaw calculation:
   y%TowerBaseMotion14%TranslationDisp(1,1) =     m%RtHS%rZ(1)
   y%TowerBaseMotion14%TranslationDisp(2,1) = -1.*m%RtHS%rZ(3)
   y%TowerBaseMotion14%TranslationDisp(3,1) =     m%RtHS%rZ(2) + p%PtfmRefzt
   
   y%TowerBaseMotion14%TranslationDisp      = y%TowerBaseMotion14%TranslationDisp - y%TowerBaseMotion14%Position   
      
   y%TowerBaseMotion14%RotationVel(1,1)     =     m%RtHS%AngVelEX(1)
   y%TowerBaseMotion14%RotationVel(2,1)     = -1.*m%RtHS%AngVelEX(3)
   y%TowerBaseMotion14%RotationVel(3,1)     =     m%RtHS%AngVelEX(2) 
   
   !...............................................................................................................................
   ! Outputs required for HydroDyn
   !...............................................................................................................................
   
   y%PlatformPtMesh%TranslationDisp(1,1) = x%QT(DOF_Sg)
   y%PlatformPtMesh%TranslationDisp(2,1) = x%QT(DOF_Sw)
   y%PlatformPtMesh%TranslationDisp(3,1) = x%QT(DOF_Hv)
   
   y%PlatformPtMesh%RotationVel(1,1)    = x%QDT(DOF_R )
   y%PlatformPtMesh%RotationVel(2,1)    = x%QDT(DOF_P )
   y%PlatformPtMesh%RotationVel(3,1)    = x%QDT(DOF_Y )
   
   y%PlatformPtMesh%TranslationVel(1,1) = x%QDT(DOF_Sg)
   y%PlatformPtMesh%TranslationVel(2,1) = x%QDT(DOF_Sw)
   y%PlatformPtMesh%TranslationVel(3,1) = x%QDT(DOF_Hv) 
   

   CALL SmllRotTrans( 'platform displacement (ED_CalcOutput)', x%QT(DOF_R ),x%QT(DOF_P ),x%QT(DOF_Y ), &
          y%PlatformPtMesh%Orientation(:,:,1), errstat=ErrStat, errmsg=ErrMsg )
      IF (ErrStat /= ErrID_None)    ErrMsg = TRIM(ErrMsg)//' (occurred at '//TRIM(Num2LStr(t))//' s)'
     !IF (ErrStat >= AbortErrLev) RETURN

   y%PlatformPtMesh%RotationAcc(1,1) = m%QD2T(DOF_R )     
   y%PlatformPtMesh%RotationAcc(2,1) = m%QD2T(DOF_P )     
   y%PlatformPtMesh%RotationAcc(3,1) = m%QD2T(DOF_Y )     

   y%PlatformPtMesh%TranslationAcc(1,1) = m%QD2T(DOF_Sg)     
   y%PlatformPtMesh%TranslationAcc(2,1) = m%QD2T(DOF_Sw)    
   y%PlatformPtMesh%TranslationAcc(3,1) = m%QD2T(DOF_Hv)    
      
   !...............................................................................................................................
   ! Outputs required for external tower loads
   !...............................................................................................................................
      
   DO J=1,p%TwrNodes
      y%TowerLn2Mesh%TranslationDisp(1,J) =     m%RtHS%rT( 1,J) - y%TowerLn2Mesh%Position(1,J)
      y%TowerLn2Mesh%TranslationDisp(2,J) = -1.*m%RtHS%rT( 3,J) - y%TowerLn2Mesh%Position(2,J)
      y%TowerLn2Mesh%TranslationDisp(3,J) =     m%RtHS%rT( 2,J) - y%TowerLn2Mesh%Position(3,J) + p%PtfmRefzt
            
      y%TowerLn2Mesh%Orientation(1,1,J)   =     m%CoordSys%t1(J,1)
      y%TowerLn2Mesh%Orientation(3,1,J)   =     m%CoordSys%t2(J,1)
      y%TowerLn2Mesh%Orientation(2,1,J)   = -1.*m%CoordSys%t3(J,1)
      y%TowerLn2Mesh%Orientation(1,2,J)   = -1.*m%CoordSys%t1(J,3)
      y%TowerLn2Mesh%Orientation(3,2,J)   = -1.*m%CoordSys%t2(J,3)
      y%TowerLn2Mesh%Orientation(2,2,J)   =     m%CoordSys%t3(J,3)
      y%TowerLn2Mesh%Orientation(1,3,J)   =     m%CoordSys%t1(J,2)
      y%TowerLn2Mesh%Orientation(3,3,J)   =     m%CoordSys%t2(J,2)
      y%TowerLn2Mesh%Orientation(2,3,J)   = -1.*m%CoordSys%t3(J,2)     
      
      y%TowerLn2Mesh%TranslationVel(1,J)  =     m%RtHS%LinVelET(1,J)
      y%TowerLn2Mesh%TranslationVel(2,J)  = -1.*m%RtHS%LinVelET(3,J)
      y%TowerLn2Mesh%TranslationVel(3,J)  =     m%RtHS%LinVelET(2,J)
            
      y%TowerLn2Mesh%RotationVel(1,J)     =     m%RtHS%AngVelEF(1,J)
      y%TowerLn2Mesh%RotationVel(2,J)     = -1.*m%RtHS%AngVelEF(3,J)
      y%TowerLn2Mesh%RotationVel(3,J)     =     m%RtHS%AngVelEF(2,J) 
            
      y%TowerLn2Mesh%TranslationAcc(1,J)  =     LinAccET(1,J)
      y%TowerLn2Mesh%TranslationAcc(2,J)  = -1.*LinAccET(3,J)
      y%TowerLn2Mesh%TranslationAcc(3,J)  =     LinAccET(2,J)
            
      y%TowerLn2Mesh%RotationAcc(1,J)     =     AngAccEF(1,J)
      y%TowerLn2Mesh%RotationAcc(2,J)     = -1.*AngAccEF(3,J)
      y%TowerLn2Mesh%RotationAcc(3,J)     =     AngAccEF(2,J) 
      
   END DO
               
   
   ! p%TwrNodes+1 is the tower top:
   J = p%TwrNodes+1
   
   y%TowerLn2Mesh%TranslationDisp(1,J) =     m%RtHS%rO(1) - y%TowerLn2Mesh%Position(1,J)
   y%TowerLn2Mesh%TranslationDisp(2,J) = -1.*m%RtHS%rO(3) - y%TowerLn2Mesh%Position(2,J)
   y%TowerLn2Mesh%TranslationDisp(3,J) =     m%RtHS%rO(2) - y%TowerLn2Mesh%Position(3,J) + p%PtfmRefzt
   
   y%TowerLn2Mesh%Orientation(1,1,J)   =     m%CoordSys%b1(1)
   y%TowerLn2Mesh%Orientation(3,1,J)   =     m%CoordSys%b2(1)
   y%TowerLn2Mesh%Orientation(2,1,J)   = -1.*m%CoordSys%b3(1)
   y%TowerLn2Mesh%Orientation(1,2,J)   = -1.*m%CoordSys%b1(3)
   y%TowerLn2Mesh%Orientation(3,2,J)   = -1.*m%CoordSys%b2(3)
   y%TowerLn2Mesh%Orientation(2,2,J)   =     m%CoordSys%b3(3)
   y%TowerLn2Mesh%Orientation(1,3,J)   =     m%CoordSys%b1(2)
   y%TowerLn2Mesh%Orientation(3,3,J)   =     m%CoordSys%b2(2)
   y%TowerLn2Mesh%Orientation(2,3,J)   = -1.*m%CoordSys%b3(2)         
          
   y%TowerLn2Mesh%TranslationVel(1,J)  =     m%RtHS%LinVelEO(1)         
   y%TowerLn2Mesh%TranslationVel(2,J)  = -1.*m%RtHS%LinVelEO(3)   
   y%TowerLn2Mesh%TranslationVel(3,J)  =     m%RtHS%LinVelEO(2)        

   y%TowerLn2Mesh%RotationVel(1,J)     =     m%RtHS%AngVelEB(1)
   y%TowerLn2Mesh%RotationVel(2,J)     = -1.*m%RtHS%AngVelEB(3)
   y%TowerLn2Mesh%RotationVel(3,J)     =     m%RtHS%AngVelEB(2) 

   y%TowerLn2Mesh%TranslationAcc(1,J)  =     LinAccEO(1)
   y%TowerLn2Mesh%TranslationAcc(2,J)  = -1.*LinAccEO(3)
   y%TowerLn2Mesh%TranslationAcc(3,J)  =     LinAccEO(2)
   
   y%TowerLn2Mesh%RotationAcc(1,J)     =     AngAccEB(1)
   y%TowerLn2Mesh%RotationAcc(2,J)     = -1.*AngAccEB(3)
   y%TowerLn2Mesh%RotationAcc(3,J)     =     AngAccEB(2) 

   
   ! p%TwrNodes+2 is the tower base:
   J = p%TwrNodes+2

   y%TowerLn2Mesh%TranslationDisp(1,J) =     m%RtHS%rZ(1) + m%RtHS%rZT0(1) - y%TowerLn2Mesh%Position(1,J)
   y%TowerLn2Mesh%TranslationDisp(2,J) = -1.*m%RtHS%rZ(3) - m%RtHS%rZT0(3) - y%TowerLn2Mesh%Position(2,J)
   y%TowerLn2Mesh%TranslationDisp(3,J) =     m%RtHS%rZ(2) + m%RtHS%rZT0(2) - y%TowerLn2Mesh%Position(3,J) + p%PtfmRefzt
      
   y%TowerLn2Mesh%Orientation(1,1,J)   =     m%CoordSys%a1(1)
   y%TowerLn2Mesh%Orientation(3,1,J)   =     m%CoordSys%a2(1)
   y%TowerLn2Mesh%Orientation(2,1,J)   = -1.*m%CoordSys%a3(1)
   y%TowerLn2Mesh%Orientation(1,2,J)   = -1.*m%CoordSys%a1(3)
   y%TowerLn2Mesh%Orientation(3,2,J)   = -1.*m%CoordSys%a2(3)
   y%TowerLn2Mesh%Orientation(2,2,J)   =     m%CoordSys%a3(3)
   y%TowerLn2Mesh%Orientation(1,3,J)   =     m%CoordSys%a1(2)
   y%TowerLn2Mesh%Orientation(3,3,J)   =     m%CoordSys%a2(2)
   y%TowerLn2Mesh%Orientation(2,3,J)   = -1.*m%CoordSys%a3(2)

   y%TowerLn2Mesh%TranslationVel(1,J)  =     m%RtHS%LinVelET(1,0)       
   y%TowerLn2Mesh%TranslationVel(2,J)  = -1.*m%RtHS%LinVelET(3,0) 
   y%TowerLn2Mesh%TranslationVel(3,J)  =     m%RtHS%LinVelET(2,0)   
   
   y%TowerLn2Mesh%RotationVel(1,J)     =     m%RtHS%AngVelEF(1,0)
   y%TowerLn2Mesh%RotationVel(2,J)     = -1.*m%RtHS%AngVelEF(3,0)
   y%TowerLn2Mesh%RotationVel(3,J)     =     m%RtHS%AngVelEF(2,0) 
   
   y%TowerLn2Mesh%TranslationAcc(1,J)  =     LinAccET(1,0)
   y%TowerLn2Mesh%TranslationAcc(2,J)  = -1.*LinAccET(3,0)
   y%TowerLn2Mesh%TranslationAcc(3,J)  =     LinAccET(2,0)
   
   y%TowerLn2Mesh%RotationAcc(1,J)     =     AngAccEF(1,0)
   y%TowerLn2Mesh%RotationAcc(2,J)     = -1.*AngAccEF(3,0)
   y%TowerLn2Mesh%RotationAcc(3,J)     =     AngAccEF(2,0) 
   
   !...............................................................................................................................
   ! Outputs required for ServoDyn
   !...............................................................................................................................
   
   y%Yaw      = x%QT( DOF_Yaw)
   y%YawRate  = x%QDT(DOF_Yaw)
   y%YawAngle = x%QT( DOF_Yaw) + x%QT(DOF_Y)  !crude approximation for yaw error... (without subtracting it from the wind direction)   
   y%BlPitch  = u%BlPitchCom !OtherState%BlPitch
   y%LSS_Spd  = x%QDT(DOF_GeAz)
   y%HSS_Spd  = ABS(p%GBRatio)*x%QDT(DOF_GeAz)
   y%RotSpeed = x%QDT(DOF_GeAz) + x%QDT(DOF_DrTr)
   
   IF ( t > 0.0_DbKi  )  THEN

      ! Calculate tower-top acceleration (fore-aft mode only) in the tower-top system:

      LinAccEO = m%RtHS%LinAccEOt
      DO I = 1,p%DOFs%NPTE  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the yaw bearing center of mass (point O)
         LinAccEO = LinAccEO + m%RtHS%PLinVelEO(p%DOFs%PTE(I),0,:)*m%QD2T(p%DOFs%PTE(I))
      ENDDO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the yaw bearing center of mass (point O)

      y%TwrAccel = DOT_PRODUCT( LinAccEO, m%CoordSys%b1 )
   ELSE
      y%TwrAccel = 0
   END IF      
   
   !Control outputs for Bladed DLL:
   y%RotPwr    = m%AllOuts(  RotPwr)*1000.
   DO K=1,p%NumBl
      y%RootMxc(K)   = m%AllOuts( RootMxc(K) )*1000.
      y%RootMyc(K)   = m%AllOuts( RootMyc(K) )*1000.
   END DO
   y%YawBrTAxp = m%AllOuts( YawBrTAxp)
   y%YawBrTAyp = m%AllOuts( YawBrTAyp)
   !y%LSSTipPxa = m%AllOuts( LSSTipPxa)*D2R ! bjj: did this above already

   y%LSSTipMxa = m%AllOuts(LSShftMxa)*1000.
   y%LSSTipMya = m%AllOuts(LSSTipMya)*1000.                ! Rotating hub My (GL co-ords) (Nm)
   y%LSSTipMza = m%AllOuts(LSSTipMza)*1000.                ! Rotating hub Mz (GL co-ords) (Nm)
   y%LSSTipMys = m%AllOuts(LSSTipMys)*1000.                ! Fixed hub My (GL co-ords) (Nm)
   y%LSSTipMzs = m%AllOuts(LSSTipMzs)*1000.                ! Fixed hub Mz (GL co-ords) (Nm)
   y%YawBrMyn  = m%AllOuts( YawBrMyn)*1000.                ! Yaw bearing My (GL co-ords) (Nm) !tower accel
   y%YawBrMzn  = m%AllOuts( YawBrMzn)*1000.                ! Yaw bearing Mz (GL co-ords) (Nm)
   y%NcIMURAxs = m%AllOuts(NcIMURAxs)*D2R                  ! Nacelle roll    acceleration (rad/s^2) -- this is in the shaft (tilted) coordinate system, instead of the nacelle (nontilted) coordinate system
   y%NcIMURAys = m%AllOuts(NcIMURAys)*D2R                  ! Nacelle nodding acceleration (rad/s^2)
   y%NcIMURAzs = m%AllOuts(NcIMURAzs)*D2R                  ! Nacelle yaw     acceleration (rad/s^2) -- this is in the shaft (tilted) coordinate system, instead of the nacelle (nontilted) coordinate system
   y%LSShftFxa = m%AllOuts(LSShftFxa)*1000.                ! Rotating low-speed shaft force x (GL co-ords) (N)
   y%LSShftFys = m%AllOuts(LSShftFys)*1000.                ! Nonrotating low-speed shaft force y (GL co-ords) (N)
   y%LSShftFzs = m%AllOuts(LSShftFzs)*1000.                ! Nonrotating low-speed shaft force z (GL co-ords) (N)

               
   RETURN
   

END SUBROUTINE ED_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for computing derivatives of continuous states.
SUBROUTINE ED_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )
!..................................................................................................................................

   REAL(DbKi),                   INTENT(IN   )  :: t           !< Current simulation time in seconds
   TYPE(ED_InputType),           INTENT(IN   )  :: u           !< Inputs at t
   TYPE(ED_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(ED_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
   TYPE(ED_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
   TYPE(ED_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t
   TYPE(ED_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states
   TYPE(ED_MiscVarType),         INTENT(INOUT)  :: m           !< Misc variables for optimization (not copied in glue code)
   TYPE(ED_ContinuousStateType), INTENT(INOUT)  :: dxdt        !< Continuous state derivatives at t [intent in so we don't need to allocate/deallocate constantly]
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! LOCAL variables
   LOGICAL, PARAMETER           :: UpdateValues  = .TRUE.      ! determines if the OtherState values need to be updated
      
   INTEGER(IntKi)                         :: I                 ! Loops through some or all of the DOFs.
   INTEGER(IntKi)                         :: ErrStat2          ! The error status code
   CHARACTER(ErrMsgLen)                   :: ErrMsg2           ! The error message, if an error occurred
   CHARACTER(*), PARAMETER                :: RoutineName = 'ED_CalcContStateDeriv'
   
      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""
   
   !OtherState%HSSBrTrqC = SIGN( u%HSSBrTrqC, x%QDT(DOF_GeAz) ) !need correct value of x%QDT(DOF_GeAz) here

         ! Compute the first time derivatives of the continuous states here:

   ! See if the values stored in m%RtHS and m%CoordSys need to be updated:
   IF ( UpdateValues ) THEN       
      
       !OtherState%BlPitch = u%BlPitchCom
       
         ! set the coordinate system variables:
      CALL SetCoordSy( t, m%CoordSys, m%RtHS, u%BlPitchCom, p, x, ErrStat2, ErrMsg2 )
         call setErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN
   
      CALL CalculatePositions(        p, x, m%CoordSys,    m%RtHS ) ! calculate positions
      CALL CalculateAngularPosVelPAcc(p, x, m%CoordSys,    m%RtHS ) ! calculate angular positions, velocities, and partial accelerations, including partial angular quantities
      CALL CalculateLinearVelPAcc(    p, x, m%CoordSys,    m%RtHS ) ! calculate linear velocities and partial accelerations
      CALL CalculateForcesMoments(    p, x, m%CoordSys, u, m%RtHS ) ! calculate the forces and moments (requires AeroBladeForces and AeroBladeMoments)            
      
   END IF
      
   !.....................................
   !  TeetMom,  RFrlMom, TFrlMom (possibly from user routines)
   ! bjj: we will want to revisit these routines:
   !.....................................
   
      ! Compute the moments from teeter springs and dampers, rotor-furl springs and dampers, tail-furl springs and dampers

   CALL Teeter  ( t, p, m%RtHS%TeetAng, m%RtHS%TeetAngVel, m%RtHS%TeetMom ) ! Compute moment from teeter     springs and dampers, TeetMom; NOTE: TeetMom will be zero for a 3-blader since TeetAng = TeetAngVel = 0
   CALL RFurling( t, p, x%QT(DOF_RFrl),          x%QDT(DOF_RFrl),            m%RtHS%RFrlMom ) ! Compute moment from rotor-furl springs and dampers, RFrlMom
   CALL TFurling( t, p, x%QT(DOF_TFrl),          x%QDT(DOF_TFrl),            m%RtHS%TFrlMom ) ! Compute moment from tail-furl  springs and dampers, TFrlMom
   
   !bjj: note m%RtHS%GBoxEffFac needed in OtherState only to fix HSSBrTrq (and used in FillAugMat)
   m%RtHS%GBoxEffFac  = p%GBoxEff**OtherState%SgnPrvLSTQ      ! = GBoxEff if SgnPrvLSTQ = 1 OR 1/GBoxEff if SgnPrvLSTQ = -1
   
   CALL FillAugMat( p, x, m%CoordSys, u, OtherState%HSSBrTrq, m%RtHS, m%AugMat )
   

   ! Invert the matrix to solve for the accelerations. The accelerations are returned by Gauss() in the first NActvDOF elements
   !   of the solution vector, SolnVec(). These are transfered to the proper index locations of the acceleration vector QD2T()
   !   using the vector subscript array SrtPS(), after Gauss() has been called:
   ! NOTE: QD2T( SrtPS(1:NActvDOF) ) cannot be sent directly because arrays sections with vector subscripts must not be used 
   !   in INTENT(OUT) arguments.

   IF ( p%DOFs%NActvDOF > 0 ) THEN
      m%AugMat_factor = m%AugMat( p%DOFs%SrtPS( 1:p%DOFs%NActvDOF ), p%DOFs%SrtPSNAUG(1:p%DOFs%NActvDOF) )
      m%SolnVec       = m%AugMat( p%DOFs%SrtPS( 1:p%DOFs%NActvDOF ), p%DOFs%SrtPSNAUG(1+p%DOFs%NActvDOF) )
   

      CALL LAPACK_getrf( M=p%DOFs%NActvDOF, N=p%DOFs%NActvDOF, A=m%AugMat_factor, IPIV=m%AugMat_pivot, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)   
         IF ( ErrStat >= AbortErrLev ) RETURN
      
      CALL LAPACK_getrs( TRANS='N',N=p%DOFs%NActvDOF, A=m%AugMat_factor,IPIV=m%AugMat_pivot, B=m%SolnVec, ErrStat=ErrStat2, ErrMsg=ErrMsg2)
   
      !CALL GaussElim( m%AugMat( p%DOFs%SrtPS(    1: p%DOFs%NActvDOF   ),     &
      !                                   p%DOFs%SrtPSNAUG(1:(p%DOFs%NActvDOF+1)) ),   &
      !                                   p%DOFs%NActvDOF,  SolnVec, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN
   END IF
   

   !bjj: if the deriv is INTENT(OUT), this is reallocated each time:
IF (.NOT. ALLOCATED(dxdt%qt) ) THEN
   CALL AllocAry( dxdt%qt,  SIZE(x%qt),  'dxdt%qt',  ErrStat2, ErrMsg2 )
   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   IF ( ErrStat >= AbortErrLev ) RETURN
END IF

IF (.NOT. ALLOCATED(dxdt%qdt) ) THEN
   CALL AllocAry( dxdt%qdt, SIZE(x%qdt), 'dxdt%qdt', ErrStat2, ErrMsg2 )
   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   IF ( ErrStat >= AbortErrLev ) RETURN
END IF

   dxdt%QT = x%QDT
   
   dxdt%QDT = 0.0
   DO I = 1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs
      dxdt%QDT(p%DOFs%SrtPS(I)) = m%SolnVec(I)    ! dxdt%QDT = m%QD2T
   ENDDO             ! I - All active (enabled) DOFs

   m%QD2T = dxdt%QDT
      
   
      ! Let's calculate the sign (+/-1) of the low-speed shaft torque for this time step and store it in SgnPrvLSTQ.
      !  This will be used during the next call to RtHS (bjj: currently violates framework, but DOE wants a hack for HSS brake).
      ! need m%QD2T set before calling this
      
   !OtherState%SgnPrvLSTQ = SignLSSTrq(p, m)
   
   
END SUBROUTINE ED_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for updating discrete states
SUBROUTINE ED_UpdateDiscState( t, n, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                   INTENT(IN   )  :: t           !< Current simulation time in seconds
      INTEGER(IntKi),               INTENT(IN   )  :: n           !< Current step of the simulation: t = n*Interval
      TYPE(ED_InputType),           INTENT(IN   )  :: u           !< Inputs at t
      TYPE(ED_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(ED_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
      TYPE(ED_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Input: Discrete states at t;
                                                                  !!   Output: Discrete states at t + Interval
      TYPE(ED_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t
      TYPE(ED_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states
      TYPE(ED_MiscVarType),         INTENT(INOUT)  :: m           !< Misc variables for optimization (not copied in glue code)
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Update discrete states here:

      ! StateData%DiscState =

END SUBROUTINE ED_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for solving for the residual of the constraint state equations
SUBROUTINE ED_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, m, z_residual, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                   INTENT(IN   )  :: Time        !< Current simulation time in seconds
      TYPE(ED_InputType),           INTENT(IN   )  :: u           !< Inputs at Time
      TYPE(ED_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(ED_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(ED_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at Time
      TYPE(ED_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time (possibly a guess)
      TYPE(ED_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states
      TYPE(ED_MiscVarType),         INTENT(INOUT)  :: m           !< Misc variables for optimization (not copied in glue code)
      TYPE(ED_ConstraintStateType), INTENT(  OUT)  :: z_residual  !< Residual of the constraint state equations using
                                                                  !!     the input values described above
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Solve for the constraint states here:

      z_residual%DummyConstrState = 0.

END SUBROUTINE ED_CalcConstrStateResidual
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine sets the parameters, based on the data stored in InputFileData
SUBROUTINE ED_SetParameters( InitInp, InputFileData, p, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(ED_InitInputType),   INTENT(IN   )    :: InitInp        !< Input data for initialization routine
   TYPE(ED_InputFile),       INTENT(IN)       :: InputFileData  !< Data stored in the module's input file
   TYPE(ED_ParameterType),   INTENT(INOUT)    :: p              !< The module's parameter data
   INTEGER(IntKi),           INTENT(OUT)      :: ErrStat        !< The error status code
   CHARACTER(*),             INTENT(OUT)      :: ErrMsg         !< The error message, if an error occurred

      ! Local variables
!   INTEGER(IntKi)                             :: K              ! Loop counter (for blades)
   INTEGER(IntKi)                             :: ErrStat2       ! Temporary error ID
   CHARACTER(ErrMsgLen)                       :: ErrMsg2        ! Temporary message describing error

      ! Initialize variables

   ErrStat = ErrID_None
   ErrMsg  = ''



      ! Set parameters from primary input file
   CALL SetPrimaryParameters( InitInp, p, InputFileData, ErrStat2, ErrMsg2  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   p%DT24 = p%DT/24.0_DbKi    ! Time-step parameter needed for Solver().


      ! Set furling parameters
   CALL SetFurlParameters( p, InputFileData, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! Set blade parameters
   CALL SetBladeParameters( p, InputFileData%InpBl, InputFileData%InpBlMesh, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! Set tower parameters
   CALL SetTowerParameters( p, InputFileData, ErrStat2, ErrMsg2 ) ! It requires p%TwrFlexL, and p%TwrNodes to be set first.
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! Set the remaining (calculuated) parameters (basically the former Coeff routine)
   CALL SetOtherParameters( p, InputFileData, ErrStat2, ErrMsg2 ) ! requires MANY things to be set first
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   CALL AllBldNdOuts_SetParameters( p, InputFileData, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
    !p%BldNd_NumOuts = 0_IntKi
    !p%BldNd_TotNumOuts = 0_IntKi

CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'ED_SetParameters:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
         END IF

      END IF


   END SUBROUTINE CheckError

END SUBROUTINE ED_SetParameters
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine initializes the ActiveDOF data type as well as the variables related to DOFs, including p%NAug and p%NDOF.
!! It assumes that p\%NumBl is set.
SUBROUTINE Init_DOFparameters( InputFileData, p, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(ED_InputFile),       INTENT(IN)       :: InputFileData  !< Data stored in the module's input file
   TYPE(ED_ParameterType),   INTENT(INOUT)    :: p              !< The module's parameter data
   INTEGER(IntKi),           INTENT(OUT)      :: ErrStat        !< The error status code
   CHARACTER(*),             INTENT(OUT)      :: ErrMsg         !< The error message, if an error occurred

      ! Local variables
   INTEGER(IntKi)                             :: K              ! Loop counter (for blades)

      ! Initialize variables

   ErrStat = ErrID_None
   ErrMsg  = ''

   IF ( p%NumBl == 1 )  THEN
      p%NDOF = 18 
   ELSEIF ( p%NumBl == 2 )  THEN
      p%NDOF = 22
   ELSE
      p%NDOF = ED_MaxDOFs
   ENDIF

   p%NAug = p%NDOF + 1

   ! ...........................................................................................................................
   ! allocate and set DOF_Flag and DOF_Desc
   ! ...........................................................................................................................
   CALL AllocAry( p%DOF_Flag, p%NDOF,   'DOF_Flag',  ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%DOF_Desc, p%NDOF,   'DOF_Desc',  ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN


   p%DOF_Flag = .false.
   p%DOF_Desc = ''
   IF ( p%NumBl == 2 )  THEN ! the 3rd blade overwrites the DOF_Teet position of the array, so don't use an "ELSE" for this statement
      p%DOF_Flag(DOF_Teet) = InputFileData%TeetDOF
      p%DOF_Desc(DOF_Teet) = 'Hub teetering DOF (internal DOF index = DOF_Teet), rad'
   END IF !


   DO K = 1,p%NumBl
      p%DOF_Flag( DOF_BF(K,1) ) = InputFileData%FlapDOF1
      p%DOF_Desc( DOF_BF(K,1) ) = '1st flapwise bending-mode DOF of blade '//TRIM(Num2LStr( K ))// &
                                  ' (internal DOF index = DOF_BF('         //TRIM(Num2LStr( K ))//',1)), m'

      p%DOF_Flag( DOF_BE(K,1) ) = InputFileData%EdgeDOF
      p%DOF_Desc( DOF_BE(K,1) ) = '1st edgewise bending-mode DOF of blade '//TRIM(Num2LStr( K ))// &
                                  ' (internal DOF index = DOF_BE('         //TRIM(Num2LStr( K ))//',1)), m'

      p%DOF_Flag( DOF_BF(K,2) ) = InputFileData%FlapDOF2
      p%DOF_Desc( DOF_BF(K,2) ) = '2nd flapwise bending-mode DOF of blade '//TRIM(Num2LStr( K ))// &
                                  ' (internal DOF index = DOF_BF('         //TRIM(Num2LStr( K ))//',2)), m'
   ENDDO          ! K - All blades

   p%DOF_Flag(DOF_DrTr) = InputFileData%DrTrDOF
   p%DOF_Desc(DOF_DrTr) = 'Drivetrain rotational-flexibility DOF (internal DOF index = DOF_DrTr), rad'
   p%DOF_Flag(DOF_GeAz) = InputFileData%GenDOF
   p%DOF_Desc(DOF_GeAz) = 'Variable speed generator DOF (internal DOF index = DOF_GeAz), rad'
   p%DOF_Flag(DOF_RFrl) = InputFileData%RFrlDOF
   p%DOF_Desc(DOF_RFrl) = 'Rotor-furl DOF (internal DOF index = DOF_RFrl), rad'
   p%DOF_Flag(DOF_TFrl) = InputFileData%TFrlDOF
   p%DOF_Desc(DOF_TFrl) = 'Tail-furl DOF (internal DOF index = DOF_TFrl), rad'
   p%DOF_Flag(DOF_Yaw ) = InputFileData%YawDOF
   p%DOF_Desc(DOF_Yaw ) = 'Nacelle yaw DOF (internal DOF index = DOF_Yaw), rad'
   p%DOF_Flag(DOF_TFA1) = InputFileData%TwFADOF1
   p%DOF_Desc(DOF_TFA1) = '1st tower fore-aft bending mode DOF (internal DOF index = DOF_TFA1), m'
   p%DOF_Flag(DOF_TSS1) = InputFileData%TwSSDOF1
   p%DOF_Desc(DOF_TSS1) = '1st tower side-to-side bending mode DOF (internal DOF index = DOF_TSS1), m'
   p%DOF_Flag(DOF_TFA2) = InputFileData%TwFADOF2
   p%DOF_Desc(DOF_TFA2) = '2nd tower fore-aft bending mode DOF (internal DOF index = DOF_TFA2), m'
   p%DOF_Flag(DOF_TSS2) = InputFileData%TwSSDOF2
   p%DOF_Desc(DOF_TSS2) = '2nd tower side-to-side bending mode DOF (internal DOF index = DOF_TSS2), m'
   p%DOF_Flag(DOF_Sg  ) = InputFileData%PtfmSgDOF
   p%DOF_Desc(DOF_Sg  ) = 'Platform horizontal surge translation DOF (internal DOF index = DOF_Sg), m'
   p%DOF_Flag(DOF_Sw  ) = InputFileData%PtfmSwDOF
   p%DOF_Desc(DOF_Sw  ) = 'Platform horizontal sway translation DOF (internal DOF index = DOF_Sw), m'
   p%DOF_Flag(DOF_Hv  ) = InputFileData%PtfmHvDOF
   p%DOF_Desc(DOF_Hv  ) = 'Platform vertical heave translation DOF (internal DOF index = DOF_Hv), m'
   p%DOF_Flag(DOF_R   ) = InputFileData%PtfmRDOF
   p%DOF_Desc(DOF_R   ) = 'Platform roll tilt rotation DOF (internal DOF index = DOF_R), rad'
   p%DOF_Flag(DOF_P   ) = InputFileData%PtfmPDOF
   p%DOF_Desc(DOF_P   ) = 'Platform pitch tilt rotation DOF (internal DOF index = DOF_P), rad'
   p%DOF_Flag(DOF_Y   ) = InputFileData%PtfmYDOF
   p%DOF_Desc(DOF_Y   ) = 'Platform yaw rotation DOF (internal DOF index = DOF_Y), rad'

   ! ...........................................................................................................................
   ! allocate the arrays stored in the p%DOFs structure:
   ! ...........................................................................................................................

      ! BJJ: note that this method will cause an error if allocating data that has already been allocated...

   ALLOCATE ( p%DOFs%NPSBE(p%NumBl), p%DOFs%NPSE(p%NumBl),  STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ExitThisRoutine( ErrID_Fatal, ' Could not allocate memory for the ActiveAOFs NPSBE and NPSE arrays.' )
      RETURN
   ENDIF


   ALLOCATE ( p%DOFs%PCE(p%NDOF), p%DOFs%PDE(p%NDOF), p%DOFs%PIE(p%NDOF), STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ExitThisRoutine( ErrID_Fatal, ' Could not allocate memory for the ActiveAOFs PCE, PDE, and PIE arrays.' )
      RETURN
   ENDIF


   ALLOCATE (  p%DOFs%PTTE(p%NDOF), p%DOFs%PTE(p%NDOF), p%DOFs%PS(p%NDOF), STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ExitThisRoutine( ErrID_Fatal, ' Could not allocate memory for the ActiveAOFs PTTE, PTE, and PS arrays.' )
      RETURN
   ENDIF


   ALLOCATE ( p%DOFs%PUE(p%NDOF), p%DOFs%PYE(p%NDOF),  STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ExitThisRoutine( ErrID_Fatal, ' Could not allocate memory for the ActiveAOFs PUE and PYE arrays.' )
      RETURN
   ENDIF


!bjj was   ALLOCATE ( p%DOFs%PSBE(p%NumBl,3), p%DOFs%PSE(p%NumBl,p%NDOF),  STAT=ErrStat )
   ALLOCATE ( p%DOFs%PSBE(p%NumBl,(NumBE+NumBF)), p%DOFs%PSE(p%NumBl,p%NDOF),  STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ExitThisRoutine( ErrID_Fatal, ' Could not allocate memory for the ActiveAOFs PSBE and PSE arrays.' )
      RETURN
   ENDIF


   ALLOCATE ( p%DOFs%SrtPS(p%NDOF), p%DOFs%SrtPSNAUG(p%NAug),  p%DOFs%Diag(p%NDOF), STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ExitThisRoutine( ErrID_Fatal, ' Could not allocate memory for the ActiveAOFs SrtPS, SrtPSNAUG, and Diag arrays.' )
      RETURN
   ENDIF


   !...............................................................................................................................
   ! Allocate and Initialize arrays for DOFS that contribute to the angular velocity of the hub and blade elements
   !...............................................................................................................................
   ! Define arrays of DOF indices (pointers) that contribute to the angular
   !   velocities of each rigid body of the wind turbine in the inertia frame:
   ! NOTE: We must include ALL of the appropriate DOF indices in these arrays,
   !       not just the indices of the enabled DOFs, since disabling a DOF only
   !       implies that each DOF acceleration is zero--it does not imply
   !       that each DOF velocity is zero (for example, consider disabled
   !       generator DOF, which still spins at constant speed).



   IF ( p%NumBl == 2 )  THEN ! 2-blader
      p%NPH = 12                         ! Number of DOFs that contribute to the angular velocity of the hub            (body H) in the inertia frame.
      p%NPM = 15                         ! Number of DOFs that contribute to the angular velocity of the blade elements (body M) in the inertia frame.
   ELSE                    ! 3-blader
      p%NPH = 11                         ! Number of DOFs that contribute to the angular velocity of the hub            (body H) in the inertia frame.
      p%NPM = 14                         ! Number of DOFs that contribute to the angular velocity of the blade elements (body M) in the inertia frame.
   ENDIF


   ALLOCATE ( p%PH(p%NPH),  p%PM(p%NumBl,p%NPM), STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ExitThisRoutine( ErrID_Fatal, ' Could not allocate memory for the ActiveDOFs PH and PM arrays.' )
      RETURN
   ENDIF

      ! Array of DOF indices (pointers) that contribute to the angular velocity of the hub (body H) in the inertia frame:
   p%PH(1:11) = (/ DOF_R, DOF_P, DOF_Y, DOF_TFA1, DOF_TSS1, DOF_TFA2, DOF_TSS2, DOF_Yaw, DOF_RFrl, DOF_GeAz, DOF_DrTr /)

   IF ( p%NumBl == 2 )  THEN ! 2-blader (add DOF_Teet to the arrays)

      p%PH(12) = DOF_Teet

         ! Array of DOF indices (pointers) that contribute to the angular velocity of the blade elements (body M) in the inertia frame:
      DO K = 1,p%NumBl ! Loop through all blades
         p%PM(K,:) = (/ DOF_R, DOF_P, DOF_Y, DOF_TFA1, DOF_TSS1, DOF_TFA2, DOF_TSS2, DOF_Yaw, DOF_RFrl, DOF_GeAz, DOF_DrTr, &
                        DOF_Teet,  DOF_BF(K,1) , DOF_BE(K,1)    , DOF_BF(K,2)          /)
      ENDDO          ! K - All blades

   ELSE  ! 3-blader

         ! Array of DOF indices (pointers) that contribute to the angular velocity of the blade elements (body M) in the inertia frame:
      DO K = 1,p%NumBl ! Loop through all blades
         p%PM(K,:) = (/ DOF_R, DOF_P, DOF_Y, DOF_TFA1, DOF_TSS1, DOF_TFA2, DOF_TSS2, DOF_Yaw, DOF_RFrl, DOF_GeAz, DOF_DrTr, &
                                   DOF_BF(K,1) , DOF_BE(K,1)    , DOF_BF(K,2)         /)
      ENDDO          ! K - All blades

   ENDIF


   !...............................................................................................................................
   ! Calculate the number of active (enabled) DOFs in the model, p%DOFs%NActvDOF:
   !...............................................................................................................................
   CALL SetEnabledDOFIndexArrays( p )

   RETURN


CONTAINS
   !............................................................................................................................
   SUBROUTINE ExitThisRoutine(ErrID,Msg)
   ! This subroutine cleans up all the allocatable arrays, closes the file, and sets the error status/message
   !............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error ID (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)

         ! Set error status/message

      ErrStat = ErrID
      ErrMsg  = Msg
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg = 'Error in Init_DOFparameters: '//TRIM(ErrMsg)
      END IF


   END SUBROUTINE ExitThisRoutine

END SUBROUTINE Init_DOFparameters
!----------------------------------------------------------------------------------------------------------------------------------
!> SHP calculates the Derive-derivative of the shape function ModShpAry at Fract.
!! NOTE: This function only works for Deriv = 0, 1, or 2.
FUNCTION SHP(Fract, FlexL, ModShpAry, Deriv, ErrStat, ErrMsg)
!..................................................................................................................................

      ! Passed variables:

   REAL(ReKi),     INTENT(IN )    :: FlexL                     !< Length of flexible beam, (m)
   REAL(ReKi),     INTENT(IN )    :: Fract                     !< Fractional distance along flexible beam, 0<=Frac<=1
   REAL(ReKi),     INTENT(IN )    :: ModShpAry(:)              !< Array holding mode shape coefficients (2:PolyOrd)
   REAL(ReKi)                     :: SHP                       !< The shape function returned by this function.

   INTEGER(IntKi), INTENT(IN )    :: Deriv                     !< Which derivative to compute Deriv = 0 (regular function SHP), 1 (D(SHP)/DZ), 2 (D2(SHP)/DZ2)
   INTEGER(IntKi), INTENT(OUT)    :: ErrStat                   !< A error level that indicates if/what error occurred
   CHARACTER(*),   INTENT(OUT)    :: ErrMsg                    !< A message indicating the error if one occurred


      ! Local variables:

   INTEGER(IntKi)                 :: CoefTmp                   ! Temporary coefficient
   INTEGER(IntKi)                 :: I                         ! Counts through polynomial array.
   INTEGER(IntKi)                 :: J                         ! I+1
   INTEGER(IntKi)                 :: Swtch(0:2)                ! Corresponds to which derivative to compute.  Sets all portions of the coefficient = 0 except those that are relevant.


   IF ( Deriv < 0 .OR. Deriv > 2 ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Function SHP input Deriv='//TRIM(Num2LStr(Deriv))//' is invalid. Deriv must be 0, 1, or 2.'
      RETURN
   ELSEIF ( Fract < 0.0_ReKi .OR. Fract > 1.0_ReKi ) THEN
      ErrStat = ErrID_Warn
      ErrMsg  = 'Function SHP input Fract='//TRIM(Num2LStr(Fract))//' does not meet the condition 0<=Fract<=1.'
   ELSE
      ErrStat = ErrID_None
   END IF

   Swtch        = 0 ! Initialize Swtch(:) to 0
   Swtch(Deriv) = 1
   SHP          = 0.0

   DO I = 1,SIZE(ModShpAry,DIM=1,KIND=IntKi) ! =2,PolyOrd
      J = I + 1
      CoefTmp = Swtch(0) + Swtch(1)*J + Swtch(2)*I*J

      IF ( (J == 2) .AND. (Deriv == 2) ) THEN !bjj this could be removed as Fract**0 = 1 (0**0 = 1 in Fortran)
         SHP =       ModShpAry(I)*CoefTmp                         /( FlexL**Deriv )
      ELSE
         SHP = SHP + ModShpAry(I)*CoefTmp*( Fract**( J - Deriv ) )/( FlexL**Deriv )
      ENDIF
   ENDDO !I

   RETURN

END FUNCTION SHP
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine allocates the coordinate systems in the ED_CoordSys type.
SUBROUTINE Alloc_CoordSys( CoordSys, p, ErrStat, ErrMsg )
!..................................................................................................................................

IMPLICIT NONE

   ! passed arguments

TYPE(ED_CoordSys),        INTENT(OUT) :: CoordSys       !< The coordinate systems, with arrays to be allocated
TYPE(ED_ParameterType),   INTENT(IN)  :: p              !< Parameters of the structural dynamics module

INTEGER(IntKi),           INTENT(OUT) :: ErrStat        !< Error status
CHARACTER(*),             INTENT(OUT) :: ErrMsg         !< Err msg


   ! local variables

CHARACTER(200), PARAMETER        :: ErrTxt = 'coordinate system arrays in SUBROUTINE Alloc_CoordSys.'


   ! Initialize ErrStat and ErrMsg

ErrStat = ErrID_None
ErrMsg  = ""


  ! Allocate coordinate system arrays:

ALLOCATE ( CoordSys%i1(p%NumBl,3), CoordSys%i2(p%NumBl,3), CoordSys%i3(p%NumBl,3), STAT=ErrStat ) !this argument doesn't work in IVF 10.1: , ERRMSG=ErrMsg
IF ( ErrStat /= 0 )  THEN
   ErrStat = ErrID_Fatal
   ErrMsg  = 'Error allocating the i1, i2, and i3 '//TRIM(ErrTxt)//' '//TRIM(ErrMsg)
   RETURN
END IF


ALLOCATE ( CoordSys%j1(p%NumBl,3), CoordSys%j2(p%NumBl,3), CoordSys%j3(p%NumBl,3), STAT=ErrStat ) !this argument doesn't work in IVF 10.1: , ERRMSG=ErrMsg
IF ( ErrStat /= 0 )  THEN
   ErrStat = ErrID_Fatal
   ErrMsg  = 'Error allocating the j1, j2, and j3 '//TRIM(ErrTxt)//' '//TRIM(ErrMsg)
   RETURN
END IF


ALLOCATE ( CoordSys%m1(p%NumBl,p%BldNodes,3), CoordSys%m2(p%NumBl,p%BldNodes,3), &
           CoordSys%m3(p%NumBl,p%BldNodes,3), STAT=ErrStat ) !this argument doesn't work in IVF 10.1: , ERRMSG=ErrMsg
IF ( ErrStat /= 0 )  THEN
   ErrStat = ErrID_Fatal
   ErrMsg  = 'Error allocating the m1, m2, and m3 '//TRIM(ErrTxt)//' '//TRIM(ErrMsg)
   RETURN
END IF


ALLOCATE ( CoordSys%n1(p%NumBl,0:p%TipNode,3), CoordSys%n2(p%NumBl,0:p%TipNode,3), &
           CoordSys%n3(p%NumBl,0:p%TipNode,3), STAT=ErrStat ) !this argument doesn't work in IVF 10.1: , ERRMSG=ErrMsg
IF ( ErrStat /= 0 )  THEN
   ErrStat = ErrID_Fatal
   ErrMsg  = 'Error allocating the n1, n2, and n3 '//TRIM(ErrTxt)//' '//TRIM(ErrMsg)
   RETURN
END IF


ALLOCATE ( CoordSys%t1(p%TwrNodes,3), CoordSys%t2(p%TwrNodes,3), CoordSys%t3(p%TwrNodes,3), STAT=ErrStat ) !this argument doesn't work in IVF 10.1: , ERRMSG=ErrMsg
IF ( ErrStat /= 0 )  THEN
   ErrStat = ErrID_Fatal
   ErrMsg  = 'Error allocating the t1, t2, and t3 '//TRIM(ErrTxt)//' '//TRIM(ErrMsg)
   RETURN
END IF


ALLOCATE ( CoordSys%te1(p%NumBl,p%BldNodes,3), CoordSys%te2(p%NumBl,p%BldNodes,3), &
           CoordSys%te3(p%NumBl,p%BldNodes,3), STAT=ErrStat ) !this argument doesn't work in IVF 10.1: , ERRMSG=ErrMsg
IF ( ErrStat /= 0 )  THEN
   ErrStat = ErrID_Fatal
   ErrMsg  = 'Error allocating the te1, te2, and te3 '//TRIM(ErrTxt)//' '//TRIM(ErrMsg)
   RETURN
END IF


RETURN
END SUBROUTINE Alloc_CoordSys
!----------------------------------------------------------------------------------------------------------------------------------
!> This takes the blade input file data and sets the corresponding blade parameters, performing linear interpolation of the
!! input data to the specified blade mesh.
!! This routine assumes p\%HubRad and p\%BldFlexL are already set.
SUBROUTINE SetBladeParameters( p, BladeInData, BladeMeshData, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(ED_ParameterType),                      INTENT(INOUT)  :: p                 !< The parameters of the structural dynamics module
   TYPE(BladeInputData),         ALLOCATABLE,   INTENT(IN)     :: BladeInData(:)    !< Program input data for all blades
   TYPE(ED_BladeMeshInputData),  ALLOCATABLE,   INTENT(IN)     :: BladeMeshData(:)  !< Program input mesh data for all blades
   INTEGER(IntKi),                              INTENT(OUT)    :: ErrStat           !< Error status
   CHARACTER(*),                                INTENT(OUT)    :: ErrMsg            !< Error message

      ! Local variables:
   REAL(ReKi)                                                  :: x                 ! Fractional location between two points in linear interpolation
   INTEGER(IntKi )                                             :: K                 ! Blade number
   INTEGER(IntKi )                                             :: J                 ! Index for the node arrays
   INTEGER(IntKi)                                              :: InterpInd         ! Index for the interpolation routine

      ! initialize variables
   ErrStat = ErrID_None
   ErrMsg  = ''
   

   ! ..............................................................................................................................
   ! Set the blade discretization information here:
   ! ..............................................................................................................................

   DO K=1,1 ! we're going to assume the discretization is the same for all blades

      IF (p%BD4Blades) THEN
         p%BldNodes = 0
      ELSE         
         p%BldNodes = BladeMeshData(K)%BldNodes
      END IF

      p%TipNode  = p%BldNodes + 1    ! The index for the blade tip and tower top nodes

   END DO

      ! .......... Allocate arrays for the blade parameters being set in this routine ..........:

   CALL Alloc_BladeParameters( p, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN

   
   IF ( .not. p%BD4Blades) then
      
      DO K=1,1 ! we're going to assume the discretization is the same for all blades

         IF ( allocated( BladeMeshData(K)%Chord ) ) THEN      
      
            p%RNodes   = BladeMeshData(K)%RNodes - p%HubRad   ! Radius to blade analysis nodes relative to root ( 0 < RNodes(:) < p%BldFlexL ) (Convert RNodes to be relative to the hub)

            p%DRNodes(1) = 2.0*p%RNodes(1)
            DO J = 2,p%BldNodes
               p%DRNodes(J) = 2.0*( p%RNodes(J) - p%RNodes(J-1) ) - p%DRNodes(J-1)
            END DO

            p%Chord     = BladeMeshData(K)%Chord
            p%AeroTwst  = BladeMeshData(K)%AeroTwst
            p%CAeroTwst = COS(p%AeroTwst)
            p%SAeroTwst = SIN(p%AeroTwst)
         
         ELSE
         
         
               ! DRNodes (Let's use constant-spaced nodes for now, but the rest of the code is written to handle variable-spaced nodes--
               !          this will be a future input!):
            p%DRNodes = p%BldFlexL/p%BldNodes !array

               ! RNodes:
            p%RNodes(1) = 0.5*p%DRNodes(1)
            DO J=2,p%BldNodes
               p%RNodes(J) = p%RNodes( J - 1 ) + 0.5*( p%DRNodes(J) + p%DRNodes( J - 1 ) )
            END DO
         
               ! these values aren't used (at least they shouldn't be):
            p%Chord     = 0.0_ReKi
            p%AeroTwst  = 0.0_ReKi
            p%CAeroTwst = 1.0_ReKi
            p%SAeroTwst = 0.0_ReKi
                        
         END IF
         

      END DO


      ! ..............................................................................................................................
      ! Interpolate the blade properties to this discretization:
      ! ..............................................................................................................................

      ! Array definitions:

      !    Input      Interp    Description
      !    -----      ------    -----------
      !    BlFract    RNodesNorm Fractional radius (0 at root, 1 at tip)
      !    PitchAx    PitchAxis  Pitch axis (0 at LE, 1 at TE)
      !    StrcTwst   ThetaS     Structural twist
      !    BMassDen   MassB      Lineal mass density
      !    FlpStff    StiffBF    Flapwise stiffness
      !    EdgStff    StiffBE    Edgewise stiffness


         ! Define RNodesNorm() which is common to all the blades:

      p%RNodesNorm = p%RNodes/p%BldFlexL  ! Normalized radius to analysis nodes relative to hub ( 0 < RNodesNorm(:) < 1 )
      
      

         ! Perform a linear interpolation of the input data to map to the meshed data for simulation:

      DO K=1,p%NumBl
         InterpInd = 1

         p%ThetaS  (K,0)         = BladeInData(K)%StrcTwst(1)
         p%ThetaS  (K,p%TipNode) = BladeInData(K)%StrcTwst(BladeInData(K)%NBlInpSt)
      
      
         DO J=1,p%BldNodes

               ! Get the index into BlFract for all of the arrays, using the NWTC Subroutine Library
            !p%ThetaS  (K,J) = InterpStp( p%RNodesNorm(J), BladeInData(K)%BlFract, BladeInData(K)%StrcTwst, &
            !                             InterpInd, BladeInData(K)%NBlInpSt )
            p%PitchAxis(K,J) = InterpStp( p%RNodesNorm(J), BladeInData(K)%BlFract, BladeInData(K)%PitchAx, &
                                         InterpInd, BladeInData(K)%NBlInpSt )


               ! The remaining arrays will have the same x value for the linear interpolation,
               ! so we'll do it manually (with a local subroutine) instead of calling the InterpStp routine again
            IF ( BladeInData(K)%NBlInpSt < 2_IntKi ) THEN
               x         = 1.0
               InterpInd = 0
            ELSE
               x = ( p%RNodesNorm(J)                     - BladeInData(K)%BlFract(InterpInd) ) / &
                   ( BladeInData(K)%BlFract(InterpInd+1) - BladeInData(K)%BlFract(InterpInd) )
            END IF

            p%ThetaS  (K,J) = InterpAry( x, BladeInData(K)%StrcTwst, InterpInd )
            p%MassB   (K,J) = InterpAry( x, BladeInData(K)%BMassDen, InterpInd )
            p%StiffBF (K,J) = InterpAry( x, BladeInData(K)%FlpStff , InterpInd )
            p%StiffBE (K,J) = InterpAry( x, BladeInData(K)%EdgStff , InterpInd )

         END DO ! J (Blade nodes)


            ! Set the blade damping and stiffness tuner
         p%BldFDamp(K,:) = BladeInData(K)%BldFlDmp
         p%BldEDamp(K,:) = BladeInData(K)%BldEdDmp
         p%FStTunr (K,:) = BladeInData(K)%FlStTunr



            ! Set the mode shape arrays
         p%BldEdgSh(:,K) = BladeInData(K)%BldEdgSh
         p%BldFl1Sh(:,K) = BladeInData(K)%BldFl1Sh
         p%BldFl2Sh(:,K) = BladeInData(K)%BldFl2Sh


      END DO ! ( Blades )

            
   else
      
      p%ThetaS  = 0.0_ReKi
      
         ! Set the blade damping and stiffness tuner
      p%BldFDamp = 0.0_ReKi
      p%BldEDamp = 0.0_ReKi
      p%FStTunr  = 0.0_ReKi

         ! Set the mode shape arrays
      p%BldEdgSh = 0.0_ReKi
      p%BldFl1Sh = 0.0_ReKi
      p%BldFl2Sh = 0.0_ReKi      
      
   end if
   
   p%CThetaS = COS(REAL(p%ThetaS,R8Ki))
   p%SThetaS = SIN(REAL(p%ThetaS,R8Ki))
   

RETURN


CONTAINS
!..................................................................................................................................
   FUNCTION InterpAry( x, YAry, Ind )
      ! This subroutine is used to interpolate the arrays more efficiently (all arrays have the same X value)
      ! See InterpStpReal() for comparison. This assumes we already know Ind and that
      ! x = ( XVal - XAry(Ind) )/( XAry(Ind+1) - XAry(Ind) )


      REAL(ReKi),      INTENT(IN) :: x                ! the relative distance between Ind and Ind+ 1
      REAL(ReKi),      INTENT(IN) :: YAry (:)         ! Array of Y values to be interpolated.
      INTEGER(IntKi) , INTENT(IN) :: Ind              ! the index into the array

      REAL(ReKi)                  :: InterpAry        ! the value calculated in this function

      IF ( Ind >= SIZE(YAry) ) THEN
         InterpAry = YAry( SIZE(YAry) )
      ELSE
         InterpAry = ( YAry(Ind+1) - YAry(Ind) ) * x  + YAry(Ind)
      END IF

   END FUNCTION InterpAry
!..................................................................................................................................
END SUBROUTINE SetBladeParameters
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine allocates arrays for the blade parameters.
SUBROUTINE Alloc_BladeParameters( p, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(ED_ParameterType),   INTENT(INOUT)  :: p                                   !< The parameters of the structural dynamics module
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                             !< Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                              !< Err msg


      ! Allocate arrays to hold the blade analysis nodes.
   CALL AllocAry  ( p%RNodes,             p%BldNodes, 'RNodes'   , ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%DRNodes,            p%BldNodes, 'DRNodes'  , ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%Chord,              p%BldNodes, 'Chord'    , ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%AeroTwst,           p%BldNodes, 'AeroTwst' , ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%CAeroTwst,          p%BldNodes, 'CAeroTwst', ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%SAeroTwst,          p%BldNodes, 'SAeroTwst', ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN


      ! Allocate arrays to hold blade data at the analysis nodes.
   CALL AllocAry  ( p%RNodesNorm,              p%BldNodes, 'RNodesNorm' , ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%PitchAxis,   p%NumBl,    p%BldNodes, 'PitchAxis'  , ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   
   ALLOCATE( p%ThetaS( p%NumBl,0:P%TipNode) &
           , p%CThetaS(p%NumBl,0:P%TipNode) &
           , p%SThetaS(p%NumBl,0:P%TipNode), STAT=ErrStat ) 
   IF (ErrStat /= 0) then
      ErrStat = ErrID_Fatal
      ErrMsg = 'Error allocating ThetaS, CThetaS, and SThetaS'
   END IF
   
      
   CALL AllocAry  ( p%MassB,       p%NumBl,    p%BldNodes, 'MassB'      , ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%StiffBF,     p%NumBl,    p%BldNodes, 'StiffBF'    , ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%StiffBE,     p%NumBl,    p%BldNodes, 'StiffBE'    , ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN


   CALL AllocAry  ( p%BldEDamp,    p%NumBl,    NumBE,      'BldEDamp'   , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%BldFDamp,    p%NumBl,    NumBF,      'BldFDamp'   , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%FStTunr,     p%NumBl,    NumBF,      'FStTunr'    , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN


         ! Allocate space for the mode shape arrays:

   ALLOCATE( p%BldEdgSh(2:PolyOrd,p%NumBl), p%BldFl1Sh(2:PolyOrd,p%NumBl), p%BldFl2Sh(2:PolyOrd,p%NumBl), STAT = ErrStat )
   IF ( ErrStat /= 0 ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = ' Error allocating BldEdgSh, BldFl1Sh, and BldFl2Sh arrays.'
      RETURN
   END IF


END SUBROUTINE Alloc_BladeParameters
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine allocates arrays for the tower parameters.
SUBROUTINE Alloc_TowerParameters( p, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(ED_ParameterType),   INTENT(INOUT)  :: p                                   !< The parameters of the structural dynamics module
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                             !< Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                              !< Err msg




      ! Allocate arrays to hold tower data at the analysis nodes.
   CALL AllocAry  ( p%HNodesNorm,    p%TwrNodes, 'HNodesNorm', ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%HNodes,        p%TwrNodes, 'HNodes'    , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%DHNodes,       p%TwrNodes, 'DHNodes'   , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%MassT,         p%TwrNodes, 'MassT'     , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%StiffTFA,      p%TwrNodes, 'StiffTFA'  , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%StiffTSS,      p%TwrNodes, 'StiffTSS'  , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN


   !   ! these are for HydroDyn?
   !CALL AllocAry  ( p%DiamT,         p%TwrNodes, 'DiamT'     , ErrStat, ErrMsg )
   !IF ( ErrStat /= ErrID_None ) RETURN
   !CALL AllocAry  ( p%CAT,           p%TwrNodes, 'CAT'       , ErrStat, ErrMsg )
   !IF ( ErrStat /= ErrID_None ) RETURN
   !CALL AllocAry  ( p%CDT,           p%TwrNodes, 'CDT'       , ErrStat, ErrMsg )
   !IF ( ErrStat /= ErrID_None ) RETURN



   !      ! Allocate space for the mode shape arrays:
   !
   !ALLOCATE( p%BldEdgSh(2:PolyOrd,p%NumBl), p%BldFl1Sh(2:PolyOrd,p%NumBl), p%BldFl2Sh(2:PolyOrd,p%NumBl), STAT = ErrStat )
   !IF ( ErrStat /= 0 ) THEN
   !   ErrStat = ErrID_Fatal
   !   ErrMsg  = ' Error allocating BldEdgSh, BldFl1Sh, and BldFl2Sh arrays.'
   !   RETURN
   !END IF



END SUBROUTINE Alloc_TowerParameters
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the remaining parameters (replacing the former FAST Initialize routine), first allocating necessary arrays.
!! It requires p\%NDOF, p\%NumBl, p\%TTopNode, p\%TipNode to be set before calling this routine.
SUBROUTINE SetOtherParameters( p, InputFileData, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(ED_InputFile),       INTENT(IN)     :: InputFileData                !< Data stored in the module's input file
   TYPE(ED_ParameterType),   INTENT(INOUT)  :: p                            !< The parameters of the structural dynamics module
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                      !< Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                       !< Err msg




      ! Allocate the arrays needed in the Coeff routine:

   !CALL AllocAry( p%AxRedTFA, 2,       2_IntKi, p%TTopNode,         'AxRedTFA',  ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   !CALL AllocAry( p%AxRedTSS, 2,       2_IntKi, p%TTopNode,         'AxRedTSS',  ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   !CALL AllocAry( p%AxRedBld, p%NumBl, 3_IntKi, 3_IntKi, p%TipNode, 'AxRedBld',  ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%BldCG,    p%NumBl,                              'BldCG',     ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%KBF,      p%NumBl, 2_IntKi, 2_IntKi,            'KBF',       ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%KBE,      p%NumBl, 1_IntKi, 1_IntKi,            'KBE',       ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%CBF,      p%NumBl, 2_IntKi, 2_IntKi,            'CBF',       ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%CBE,      p%NumBl, 1_IntKi, 1_IntKi,            'CBE',       ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%SecondMom,p%NumBl,                              'SecondMom', ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%FirstMom, p%NumBl,                              'FirstMom',  ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%FreqBE,   p%NumBl, NumBE, 3_IntKi,              'FreqBE',    ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%FreqBF,   p%NumBl, NumBF, 3_IntKi,              'FreqBF',    ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%BldMass,  p%NumBl,                              'BldMass',   ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%rSAerCenn1,p%NumBl,p%BldNodes,  'rSAerCenn1',  ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%rSAerCenn2,p%NumBl,p%BldNodes,  'rSAerCenn2',  ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry(p%BElmntMass, p%BldNodes, p%NumBl, 'BElmntMass', ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry(p%TElmntMass, p%TwrNodes,          'TElmntMass', ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN

   !CALL AllocAry( p%AxRedBld, p%NumBl, 3_IntKi, 3_IntKi, p%TipNode, 'AxRedBld',  ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   ALLOCATE ( p%AxRedBld(p%NumBl, 3_IntKi, 3_IntKi, 0:p%TipNode) , STAT=ErrStat )
   IF ( ErrStat /= 0 ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating AxRedBld array.'
      RETURN
   END IF

   
   ALLOCATE ( p%TwrFASF(2,0:p%TTopNode,0:2) , &
              p%TwrSSSF(2,0:p%TTopNode,0:2) , & 
              p%AxRedTFA(2,2,0:p%TTopNode)  , &
              p%AxRedTSS(2,2,0:p%TTopNode)  , STAT=ErrStat )
   IF ( ErrStat /= 0 ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating TwrFASF, TwrSSSF, AxRedTFA, and p%AxRedTSS arrays.'
      RETURN
   END IF

   ALLOCATE ( p%TwistedSF(p%NumBl,2,3,0:p%TipNode,0:2) , STAT=ErrStat )
   IF ( ErrStat /= 0 ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating TwistedSF array.'
      RETURN
   END IF


   CALL Coeff(p, InputFileData, ErrStat, ErrMsg)
   IF ( ErrStat /= ErrID_None ) RETURN
   
   
END SUBROUTINE SetOtherParameters
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine allocates arrays in the RtHndSide data structure.
!! It requires p\%TwrNodes, p\%NumBl, p\%TipNode, p\%NDOF, p\%BldNodes to be set before calling this routine.
SUBROUTINE Alloc_RtHS( RtHS, p, ErrStat, ErrMsg  )
!..................................................................................................................................

   TYPE(ED_RtHndSide),       INTENT(INOUT)  :: RtHS                         !< RtHndSide data type
   TYPE(ED_ParameterType),   INTENT(IN)     :: p                            !< Parameters of the structural dynamics module
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                      !< Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                       !< Error message

   ! local variables:
   INTEGER(IntKi),   PARAMETER              :: Dims = 3                     ! The position arrays all must be allocated with a dimension for X,Y,and Z
   CHARACTER(*),     PARAMETER              :: RoutineName = 'Alloc_RtHS'

      ! positions:
  !CALL AllocAry( RtHS%rZT,       Dims, p%TwrNodes,        'rZT',       ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%rT,        Dims, p%TwrNodes,        'rT',        ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%rT0T,      Dims, p%TwrNodes,        'rT0T',      ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
  !CALL AllocAry( RtHS%rQS,       Dims, p%NumBl,p%TipNode, 'rQS',       ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
  !CALL AllocAry( RtHS%rS,        Dims, p%NumBl,p%TipNode, 'rS',        ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%rS0S,      Dims, p%NumBl,p%TipNode, 'rS0S',      ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%rPS0,      Dims, p%NumBl,           'rPS0',      ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%rSAerCen,  Dims, p%TipNode, p%NumBl,'rSAerCen',  ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
  
      ! tower
   allocate(RtHS%rZT(      Dims,  0:p%TwrNodes), &
            RtHS%AngPosEF( Dims,  0:p%TwrNodes), &
            RtHS%AngPosXF( Dims,  0:p%TwrNodes), &
            RtHS%AngVelEF( Dims,  0:p%TwrNodes), &
            RtHS%LinVelET( Dims,  0:p%TwrNodes), &
            RtHS%AngAccEFt(Dims,  0:p%TwrNodes), &
            RtHS%LinAccETt(Dims,  0:p%TwrNodes), &
      STAT=ErrStat)
      if (ErrStat /= 0) then
         ErrStat = ErrID_Fatal
         ErrMsg  = "Error allocating rZT, AngPosEF, AngPosXF, LinVelET, AngVelEF, LinAccETt, and AngAccEFt arrays."
         RETURN
      end if
               
      
      ! blades
   allocate(RtHS%rS( Dims, p%NumBl,0:p%TipNode), &
            RtHS%rQS(Dims, p%NumBl,0:p%TipNode), STAT=ErrStat)
      if (ErrStat /= 0) then
         ErrStat = ErrID_Fatal
         ErrMsg  = "Error allocating rS and rQS."
         RETURN
      end if
   

      ! angular velocities (including partial angular velocities):
   !CALL AllocAry( RtHS%AngVelEF,  Dims, p%TwrNodes,        'AngVelEF',  ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   !CALL AllocAry( RtHS%AngPosEF,  Dims, p%TwrNodes,        'AngPosEF',  ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   !CALL AllocAry( RtHS%AngPosXF,  Dims, p%TwrNodes,        'AngPosXF',  ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN


         ! These angular velocities are allocated to start numbering a dimension with 0 instead of 1:
   ALLOCATE ( RtHS%PAngVelEB(p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PAngVelEB array.'
      RETURN
   ENDIF

   ALLOCATE ( RtHS%PAngVelER(p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PAngVelER array.'
      RETURN
   ENDIF

   ALLOCATE ( RtHS%PAngVelEX(p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PAngVelEX array.'
      RETURN
   ENDIF

   ALLOCATE ( RtHS%PAngVelEA(p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PAngVelEA array.'
      RETURN
   ENDIF

   ALLOCATE ( RtHS%PAngVelEF(0:p%TwrNodes, p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PAngVelEF array.'
      RETURN
   ENDIF
   ALLOCATE ( RtHS%PAngVelEG(                  p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PAngVelEG array.'
      RETURN
   ENDIF
   ALLOCATE ( RtHS%PAngVelEH(                  p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PAngVelEH array.'
      RETURN
   ENDIF
   ALLOCATE ( RtHS%PAngVelEL(                  p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PAngVelEL array.'
      RETURN
   ENDIF
   ALLOCATE ( RtHS%PAngVelEM(p%NumBl,0:p%TipNode,p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PAngVelEM array.'
      RETURN
   ENDIF
   ALLOCATE ( RtHS%PAngVelEN(                  p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PAngVelEN array.'
      RETURN
   ENDIF

      ! angular accelerations:
   !CALL AllocAry( RtHS%AngAccEFt, Dims, p%TwrNodes,         'AngAccEFt', ErrStat, ErrMsg );  IF ( ErrStat /= ErrID_None ) RETURN

      ! linear velocities (including partial linear velocities):
   !CALL AllocAry( RtHS%LinVelET,  Dims, p%TwrNodes,         'LinVelET',  ErrStat, ErrMsg );  IF ( ErrStat /= ErrID_None ) RETURN         

   !CALL AllocAry( RtHS%LinVelESm2,                 p%NumBl, 'LinVelESm2',ErrStat, ErrMsg );  IF ( ErrStat /= ErrID_None ) RETURN ! The m2-component (closest to tip) of LinVelES
   ALLOCATE( RtHS%LinVelES( Dims, 0:p%TipNode, p%NumBl ), &
             RtHS%AngVelEM( Dims, 0:p%TipNode, p%NumBl ), STAT=ErrStat )
   IF (ErrStat /= 0 ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg = RoutineName//":Error allocating LinVelES and AngVelEM."
      RETURN
   END IF
   
            ! These linear velocities are allocated to start numbering a dimension with 0 instead of 1:

   ALLOCATE ( RtHS%PLinVelEIMU(p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PLinVelEIMU array.'
      RETURN
   ENDIF

   ALLOCATE ( RtHS%PLinVelEO(p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PLinVelEO array.'
      RETURN
   ENDIF

   ALLOCATE ( RtHS%PLinVelES(p%NumBl,0:p%TipNode,p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PLinVelES array.'
      RETURN
   ENDIF

   ALLOCATE ( RtHS%PLinVelET(0:p%TwrNodes,p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PLinVelET array.'
      RETURN
   ENDIF

   ALLOCATE ( RtHS%PLinVelEZ(p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PLinVelEZ array.'
      RETURN
   ENDIF

   ALLOCATE ( RtHS%PLinVelEC(p%NDOF,0:1,3) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PLinVelEC array.'
      RETURN
   ENDIF
   ALLOCATE ( RtHS%PLinVelED(p%NDOF,0:1,3) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PLinVelED array.'
      RETURN
   ENDIF
   ALLOCATE ( RtHS%PLinVelEI(p%NDOF,0:1,3) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PLinVelEI array.'
      RETURN
   ENDIF
   ALLOCATE ( RtHS%PLinVelEJ(p%NDOF,0:1,3) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PLinVelEJ array.'
      RETURN
   ENDIF
   ALLOCATE ( RtHS%PLinVelEP(p%NDOF,0:1,3) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PLinVelEP array.'
      RETURN
   ENDIF
   ALLOCATE ( RtHS%PLinVelEQ(p%NDOF,0:1,3) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PLinVelEQ array.'
      RETURN
   ENDIF
   
   ALLOCATE ( RtHS%PLinVelEU(p%NDOF,0:1,3) , &
              RtHS%PLinVelEV(p%NDOF,0:1,3) , &
              RtHS%PLinVelEW(p%NDOF,0:1,3) , &
              RtHS%PLinVelEY(p%NDOF,0:1,3) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PLinVelEU, PLinVelEV, PLinVelEW and PLinVelEY arrays.'
      RETURN
   ENDIF

   
   ALLOCATE( RtHS%LinAccESt( Dims, p%NumBl, 0:p%TipNode ), STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for LinAccESt.'
      RETURN
   ENDIF
   ALLOCATE ( RtHS%AngAccEKt( Dims, 0:p%TipNode, p%NumBl ) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for AngAccEKt.'
      RETURN
   ENDIF

   ALLOCATE(RtHS%AngPosHM(Dims, p%NumBl, 0:p%TipNode), STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the AngPosHM arrays.'
      RETURN
   ENDIF

   !CALL AllocAry( RtHS%LinAccESt, Dims, p%NumBl, p%TipNode,'LinAccESt', ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   !CALL AllocAry( RtHS%LinAccETt, Dims, p%TwrNodes,        'LinAccETt', ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PFrcS0B,   Dims, p%NumBl,p%NDOF,    'PFrcS0B',   ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%FrcS0Bt,   Dims, p%NumBl,           'FrcS0Bt',   ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PMomH0B,   Dims, p%NumBl, p%NDOF,   'PMomH0B',   ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%MomH0Bt,   Dims, p%NumBl,           'MomH0Bt',   ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PFrcPRot,  Dims, p%NDOF,            'PFrcPRot',  ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PMomLPRot, Dims, p%NDOF,            'PMomLPRot', ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PMomNGnRt, Dims, p%NDOF,            'PMomNGnRt', ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PMomNTail, Dims, p%NDOF,            'PMomNTail', ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PFrcONcRt, Dims, p%NDOF,            'PFrcONcRt', ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PMomBNcRt, Dims, p%NDOF,            'PMomBNcRt', ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PFrcT0Trb, Dims, p%NDOF,            'PFrcT0Trb', ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PMomX0Trb, Dims, p%NDOF,            'PMomX0Trb', ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%FSAero,    Dims, p%NumBl,p%BldNodes,'FSAero',    ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%MMAero,    Dims, p%NumBl,p%BldNodes,'MMAero',    ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%FSTipDrag, Dims, p%NumBl,           'FSTipDrag', ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PFTHydro,  Dims, p%TwrNodes, p%NDOF,'PFTHydro',  ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PMFHydro,  Dims, p%TwrNodes, p%NDOF,'PMFHydro',  ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%FTHydrot,  Dims, p%TwrNodes,        'FTHydrot',  ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%MFHydrot,  Dims, p%TwrNodes,        'MFHydrot',  ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN

   CALL AllocAry( RtHS%PFrcVGnRt, Dims, p%NDOF,            'PFrcVGnRt', ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PFrcWTail, Dims, p%NDOF,            'PFrcWTail', ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PFrcZAll,  Dims, p%NDOF,            'PFrcZAll',  ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PMomXAll,  Dims, p%NDOF,            'PMomXAll',  ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN

END SUBROUTINE Alloc_RtHS
!----------------------------------------------------------------------------------------------------------------------------------
!> This takes the tower input file data and sets the corresponding tower parameters, performing linear interpolation of the
!! input data to the specified tower mesh.
!! It requires p\%TwrFlexL, and p\%TwrNodes to be set first.
SUBROUTINE SetTowerParameters( p, InputFileData, ErrStat, ErrMsg  )
!..................................................................................................................................

      ! Passed variables

   TYPE(ED_ParameterType),   INTENT(INOUT)  :: p                            !< Parameters of the structural dynamics module
   TYPE(ED_InputFile),       INTENT(IN)     :: InputFileData                !< Data stored in the module's input file
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                      !< Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                       !< Error message

      ! Local variables:

   REAL(ReKi)                               :: x                            ! Fractional location between two points in linear interpolation
   INTEGER(IntKi )                          :: J                            ! Index for the node arrays
   INTEGER(IntKi)                           :: InterpInd                    ! Index for the interpolation routine


      ! Initialize data
   ErrStat   = ErrID_None
   ErrMsg    = ''

   CALL Alloc_TowerParameters( p, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN


   !...............................................................................................................................
   ! Define the tower discretization arrays:
   !...............................................................................................................................

      ! DHNodes (Let's use constant-spaced nodes for now, but the rest of the code is written to handle variable-spaced nodes--
      !          this will be a future input!):
   p%DHNodes = p%TwrFlexL/p%TwrNodes

      ! HNodes:
   p%HNodes(1) = 0.5*p%DHNodes(1)
   DO J=2,p%TwrNodes
      p%HNodes(J) = p%HNodes( J - 1 ) + 0.5*( p%DHNodes(J) + p%DHNodes( J - 1 ) )
   END DO

      ! HNodesNorm:
   p%HNodesNorm = p%HNodes/p%TwrFlexL


   !...............................................................................................................................
   ! Interpolate the input data to the tower discretization
   !...............................................................................................................................
   ! Array definitions:

   !    Input      Interp    Description
   !    -----      ------    -----------
   !    HtFract    HNodesNorm Fractional height (0 at top of rigid section, 1 at tower top)
   !    TMassDen   MassT      Lineal mass density
   !    TwFAStif   StiffTFA   Tower fore-aft stiffness
   !    TwSSStif   StiffTSS   Tower side-to-side stiffness

   InterpInd = 1


   DO J=1,p%TwrNodes

         ! Get the index into HtFract for all of the arrays, using the NWTC Subroutine Library
      p%MassT     (J) = InterpStp( p%HNodesNorm(J), InputFileData%HtFract, InputFileData%TMassDen, InterpInd, InputFileData%NTwInpSt )
      p%StiffTFA  (J) = InterpStp( p%HNodesNorm(J), InputFileData%HtFract, InputFileData%TwFAStif, InterpInd, InputFileData%NTwInpSt )
      p%StiffTSS  (J) = InterpStp( p%HNodesNorm(J), InputFileData%HtFract, InputFileData%TwSSStif, InterpInd, InputFileData%NTwInpSt )
   END DO ! J
   p%MassT = abs(p%MassT)
   p%StiffTFA = abs(p%StiffTFA)
   p%StiffTSS = abs(p%StiffTSS)


   !...............................................................................................................................
   ! Set other tower parameters:
   !...............................................................................................................................

   p%TTopNode = p%TwrNodes + 1

   !   ! these are for HydroDyn ?
   !p%DiamT(:) = InputFileData%TwrDiam
   !p%CAT(:)   = InputFileData%TwrCA
   !p%CDT(:)   = InputFileData%TwrCD
   !

RETURN


CONTAINS
!..................................................................................................................................
   !> This subroutine is used to interpolate the arrays more efficiently (all arrays have the same X value)
   !! See InterpStpReal() for comparison. This assumes we already know Ind and that
   !! x = ( XVal - XAry(Ind) )/( XAry(Ind+1) - XAry(Ind) )
   FUNCTION InterpAry( x, YAry, Ind )

      REAL(ReKi),      INTENT(IN) :: x                !< the relative distance between Ind and Ind+ 1
      REAL(ReKi),      INTENT(IN) :: YAry (:)         !< Array of Y values to be interpolated.
      INTEGER(IntKi) , INTENT(IN) :: Ind              !< the index into the array

      REAL(ReKi)                  :: InterpAry        !< the value calculated in this function

      InterpAry = ( YAry(Ind+1) - YAry(Ind) ) * x  + YAry(Ind)

   END FUNCTION InterpAry

END SUBROUTINE SetTowerParameters
!----------------------------------------------------------------------------------------------------------------------------------
!> This takes the furling input file data and sets the corresponding furling parameters.
SUBROUTINE SetFurlParameters( p, InputFileData, ErrStat, ErrMsg  )
!..................................................................................................................................

      ! Passed variables

   TYPE(ED_ParameterType),   INTENT(INOUT)  :: p                            !< Parameters of the structural dynamics module
   TYPE(ED_InputFile),       INTENT(IN)     :: InputFileData                !< Data stored in the module's input file
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                      !< Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                       !< Error message

      ! Local variables:

   REAL(ReKi)                               :: x                            ! Fractional location between two points in linear interpolation
!   INTEGER(IntKi )                          :: J                            ! Index for the node arrays
!   INTEGER(IntKi)                           :: InterpInd                    ! Index for the interpolation routine


      ! Initialize error data
   ErrStat = ErrID_None
   ErrMsg  = ''


      ! Direct copy of InputFileData to parameters

   p%RFrlMass = InputFileData%RFrlMass
   p%BoomMass = InputFileData%BoomMass
   p%TFinMass = InputFileData%TFinMass
   p%TFrlIner = InputFileData%TFrlIner

   p%RFrlMod  = InputFileData%RFrlMod
   p%TFrlMod  = InputFileData%TFrlMod

   p%RFrlSpr  = InputFileData%RFrlSpr
   p%RFrlDmp  = InputFileData%RFrlDmp
   p%RFrlUSSP = InputFileData%RFrlUSSP
   p%RFrlDSSP = InputFileData%RFrlDSSP
   p%RFrlDSSpr= InputFileData%RFrlDSSpr
   p%RFrlUSSpr= InputFileData%RFrlUSSpr
   p%RFrlUSDP = InputFileData%RFrlUSDP
   p%RFrlDSDP = InputFileData%RFrlDSDP
   p%RFrlUSDmp= InputFileData%RFrlUSDmp
   p%RFrlDSDmp= InputFileData%RFrlDSDmp

   p%TFrlSpr  = InputFileData%TFrlSpr
   p%TFrlDmp  = InputFileData%TFrlDmp
   p%TFrlUSSP = InputFileData%TFrlUSSP
   p%TFrlDSSP = InputFileData%TFrlDSSP
   p%TFrlUSSpr= InputFileData%TFrlUSSpr
   p%TFrlDSSpr= InputFileData%TFrlDSSpr
   p%TFrlUSDP = InputFileData%TFrlUSDP
   p%TFrlDSDP = InputFileData%TFrlDSDP
   p%TFrlUSDmp= InputFileData%TFrlUSDmp
   p%TFrlDSDmp= InputFileData%TFrlDSDmp

   p%RFrlPnt_n = InputFileData%RFrlPnt_n

   p%TFrlPnt_n = InputFileData%TFrlPnt_n


      ! Store sine/cosine values instead of some input angles:

   p%CShftSkew = COS( REAL(InputFileData%ShftSkew,R8Ki) )
   p%SShftSkew = SIN( REAL(InputFileData%ShftSkew,R8Ki) )

   p%CRFrlSkew = COS( REAL(InputFileData%RFrlSkew, R8Ki) )
   p%SRFrlSkew = SIN( REAL(InputFileData%RFrlSkew, R8Ki) )
   p%CRFrlTilt = COS( REAL(InputFileData%RFrlTilt, R8Ki) )
   p%SRFrlTilt = SIN( REAL(InputFileData%RFrlTilt, R8Ki) )

   p%CTFrlSkew = COS( REAL(InputFileData%TFrlSkew, R8Ki) )
   p%STFrlSkew = SIN( REAL(InputFileData%TFrlSkew, R8Ki) )
   p%CTFrlTilt = COS( REAL(InputFileData%TFrlTilt, R8Ki) )
   p%STFrlTilt = SIN( REAL(InputFileData%TFrlTilt, R8Ki) )


      ! Common multiplications of sines and cosines:

   p%CRFrlSkw2 = p%CRFrlSkew**2
   p%SRFrlSkw2 = p%SRFrlSkew**2
   p%CSRFrlSkw = p%CRFrlSkew*p%SRFrlSkew
   p%CRFrlTlt2 = p%CRFrlTilt**2
   p%SRFrlTlt2 = p%SRFrlTilt**2
   p%CSRFrlTlt = p%CRFrlTilt*p%SRFrlTilt

   p%CTFrlSkw2 = p%CTFrlSkew**2
   p%STFrlSkw2 = p%STFrlSkew**2
   p%CSTFrlSkw = p%CTFrlSkew*p%STfrlSkew
   p%CTFrlTlt2 = p%CTFrlTilt**2
   p%STFrlTlt2 = p%STFrlTilt**2
   p%CSTFrlTlt = p%CTFrlTilt*p%STFrlTilt


      ! Calculate some positions:

   p%rWIxn     = InputFileData%BoomCM_n(1) - p%TFrlPnt_n(1)
   p%rWIyn     = InputFileData%BoomCM_n(2) - p%TFrlPnt_n(2)
   p%rWIzn     = InputFileData%BoomCM_n(3) - p%TFrlPnt_n(3)

   p%rWJxn     = InputFileData%TFinCM_n(1) - p%TFrlPnt_n(1)
   p%rWJyn     = InputFileData%TFinCM_n(2) - p%TFrlPnt_n(2)
   p%rWJzn     = InputFileData%TFinCM_n(3) - p%TFrlPnt_n(3)

   p%rVDxn     = InputFileData%RFrlCM_n(1) - p%RFrlPnt_n(1)
   p%rVDyn     = InputFileData%RFrlCM_n(2) - p%RFrlPnt_n(2)
   p%rVDzn     = InputFileData%RFrlCM_n(3) - p%RFrlPnt_n(3)

   p%rVPxn     =        0.0_ReKi        - p%RFrlPnt_n(1)
   p%rVPyn     = InputFileData%Yaw2Shft - p%RFrlPnt_n(2)


      ! Note: These positions are also used for non-furling machines:

   p%rVPzn     = InputFileData%Twr2Shft - p%RFrlPnt_n(3)
   p%rVIMUxn   = InputFileData%NcIMUxn  - p%RFrlPnt_n(1)
   p%rVIMUyn   = InputFileData%NcIMUyn  - p%RFrlPnt_n(2)
   p%rVIMUzn   = InputFileData%NcIMUzn  - p%RFrlPnt_n(3)

END SUBROUTINE SetFurlParameters
!----------------------------------------------------------------------------------------------------------------------------------
!> This takes the primary input file data and sets the corresponding parameters.
SUBROUTINE SetPrimaryParameters( InitInp, p, InputFileData, ErrStat, ErrMsg  )
!..................................................................................................................................

      ! Passed variables

   TYPE(ED_InitInputType),   INTENT(IN   )  :: InitInp                      !< Input data for initialization routine
   TYPE(ED_ParameterType),   INTENT(INOUT)  :: p                            !< Parameters of the structural dynamics module
   TYPE(ED_InputFile),       INTENT(IN)     :: InputFileData                !< Data stored in the module's input file
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                      !< Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                       !< Error message

!bjj: ERROR CHECKING!!!

      ! Initialize error data
   ErrStat = ErrID_None
   ErrMsg  = ''

   !p%Twr2Shft  = InputFileData%Twr2Shft
   !p%HubIner   = InputFileData%HubIner
   !p%NacYIner  = InputFileData%NacYIner


   !...............................................................................................................................
   ! Direct copy of variables:
   !...............................................................................................................................
   p%NumBl     = InputFileData%NumBl
   p%TipRad    = InputFileData%TipRad
   p%HubRad    = InputFileData%HubRad
   p%method    = InputFileData%method
   p%TwrNodes  = InputFileData%TwrNodes
   p%MHK       = InitInp%MHK

   p%PtfmCMxt = InputFileData%PtfmCMxt
   p%PtfmCMyt = InputFileData%PtfmCMyt   
   
   p%DT        = InputFileData%DT
   p%OverHang  = InputFileData%OverHang
   p%ShftGagL  = InputFileData%ShftGagL
   IF ( InitInp%MHK == 1 ) THEN
      p%TowerHt   = InputFileData%TowerHt - InitInp%WtrDpth
      p%TowerBsHt = InputFileData%TowerBsHt - InitInp%WtrDpth
      p%PtfmRefzt = InputFileData%PtfmRefzt - InitInp%WtrDpth
   ELSE
      p%TowerHt   = InputFileData%TowerHt
      p%TowerBsHt = InputFileData%TowerBsHt
      p%PtfmRefzt = InputFileData%PtfmRefzt
   END IF
   
   p%HubMass   = InputFileData%HubMass
   p%GenIner   = InputFileData%GenIner
   p%NacMass   = InputFileData%NacMass
   p%YawBrMass = InputFileData%YawBrMass
   p%PtfmMass  = InputFileData%PtfmMass
   p%PtfmRIner = InputFileData%PtfmRIner
   p%PtfmPIner = InputFileData%PtfmPIner
   p%PtfmYIner = InputFileData%PtfmYIner
   p%GBoxEff   = InputFileData%GBoxEff
   p%GBRatio   = InputFileData%GBRatio
   p%DTTorSpr  = InputFileData%DTTorSpr
   p%DTTorDmp  = InputFileData%DTTorDmp


   p%NTwGages  = InputFileData%NTwGages
   p%TwrGagNd  = InputFileData%TwrGagNd
   p%NBlGages  = InputFileData%NBlGages
   p%BldGagNd  = InputFileData%BldGagNd
   !p%OutFile   = InputFileData%OutFile
   !p%OutFileFmt= InputFileData%OutFileFmt !wrbinoutput, wrtxtoutput???
   p%OutFmt    = InputFileData%OutFmt
   p%Tstart    = InputFileData%Tstart
   !p%DecFact   = InputFileData%DecFact
   p%NumOuts   = InputFileData%NumOuts

   IF ( p%NumBl == 2 ) THEN
      p%UndSling = InputFileData%UndSling
      p%TeetMod  = InputFileData%TeetMod
      p%TeetDmpP = InputFileData%TeetDmpP
      p%TeetDmp  = InputFileData%TeetDmp
      p%TeetCDmp = InputFileData%TeetCDmp
      p%TeetSStP = InputFileData%TeetSStP
      p%TeetHStP = InputFileData%TeetHStP
      p%TeetSSSp = InputFileData%TeetSSSp
      p%TeetHSSp = InputFileData%TeetHSSp
   ELSE ! Three-bladed turbines don't use these parameters, so set them to zero.
      p%UndSling = 0.0
      p%TeetMod  = 0
      p%TeetDmpP = 0.0
      p%TeetDmp  = 0.0
      p%TeetCDmp = 0.0
      p%TeetSStP = 0.0
      p%TeetHStP = 0.0
      p%TeetSSSp = 0.0
      p%TeetHSSp = 0.0
   END IF

   
   CALL AllocAry( p%TipMass, p%NumBl, 'TipMass', ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN
   p%TipMass   = InputFileData%TipMass

      ! initialize all of the DOF parameters:
   CALL Init_DOFparameters( InputFileData, p, ErrStat, ErrMsg ) !sets p%NDOF and p%NAug
      IF (ErrStat >= AbortErrLev) RETURN

      ! Set parameters for output channels:
   CALL SetOutParam(InputFileData%OutList, p, ErrStat, ErrMsg ) ! requires: p%NumOuts, p%NumBl, p%NBlGages, p%NTwGages; sets: p%OutParam.
      IF (ErrStat >= AbortErrLev) RETURN

   IF ( InputFileData%TabDelim ) THEN
      p%Delim = TAB
   ELSE
      p%Delim = ' '
   !ELSE
   !   p%Delim = ','
   END IF

   !...............................................................................................................................
   ! Calculate some indirect inputs:
   !...............................................................................................................................
   p%TwoPiNB   = TwoPi_D/p%NumBl                                                   ! 2*Pi/NumBl is used in RtHS().

   p%rZT0zt    = p%TowerBsHt - p%PtfmRefzt                                         ! zt-component of position vector rZT0.
   p%RefTwrHt  = p%TowerHt   - p%PtfmRefzt                                         ! Vertical distance between ElastoDyn's undisplaced tower height (variable TowerHt) and ElastoDyn's inertia frame reference point (variable PtfmRef).
   p%TwrFlexL  = p%TowerHt   - p%TowerBsHt                                         ! Height / length of the flexible portion of the tower.
   p%BldFlexL  = p%TipRad    - p%HubRad                                            ! Length of the flexible portion of the blade.
   if (p%BD4Blades) p%BldFlexL = 0.0_ReKi
   
   IF ( InitInp%MHK == 1 ) THEN
      p%rZYzt     = InputFileData%PtfmCMzt - InitInp%WtrDpth - p%PtfmRefzt
   ELSE
      p%rZYzt     = InputFileData%PtfmCMzt - p%PtfmRefzt
   END IF

   !...............................................................................................................................
   ! set cosine and sine of Precone and Delta3 angles:
   !...............................................................................................................................
   CALL AllocAry( p%CosPreC,  p%NumBl,                              'CosPreC',   ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%SinPreC,  p%NumBl,                              'SinPreC',   ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN

   p%CosPreC  = COS( REAL(InputFileData%Precone(1:p%NumBl),R8Ki) )
   p%SinPreC  = SIN( REAL(InputFileData%Precone(1:p%NumBl),R8Ki) )
 
   IF ( p%NumBl == 2 ) THEN
      p%CosDel3  = COS( REAL(InputFileData%Delta3,R8Ki) )
      p%SinDel3  = SIN( REAL(InputFileData%Delta3,R8Ki) )
   ELSE
      p%CosDel3  = 1.0_R8Ki
      p%SinDel3  = 0.0_R8Ki
   END IF
   
   !...............................................................................................................................

      ! Calculate the average tip radius normal to the shaft (AvgNrmTpRd)
      !   and the swept area of the rotor (ProjArea):

   p%AvgNrmTpRd = p%TipRad*SUM(p%CosPreC)/p%NumBl     ! Average tip radius normal to the saft.
   p%ProjArea   = pi*( p%AvgNrmTpRd**2 )              ! Swept area of the rotor projected onto the rotor plane (the plane normal to the low-speed shaft).

   p%RotSpeed  = InputFileData%RotSpeed               ! Rotor speed in rad/sec.
   p%CShftTilt = COS( REAL(InputFileData%ShftTilt,R8Ki) )
   p%SShftTilt = SIN( REAL(InputFileData%ShftTilt,R8Ki) )

   p%HubHt     = p%TowerHt + InputFileData%Twr2Shft + p%OverHang*p%SShftTilt


      ! Direct copy of InputFileData to parameters

   !p%FlapDOF1  = InputFileData%FlapDOF1
   !p%FlapDOF2  = InputFileData%FlapDOF2
   !p%EdgeDOF   = InputFileData%EdgeDOF
   !p%TeetDOF   = InputFileData%TeetDOF
   !p%DrTrDOF   = InputFileData%DrTrDOF
   !p%GenDOF    = InputFileData%GenDOF
   !p%YawDOF    = InputFileData%YawDOF
   !p%TwFADOF1  = InputFileData%TwFADOF1
   !p%TwFADOF2  = InputFileData%TwFADOF2
   !p%TwSSDOF1  = InputFileData%TwSSDOF1
   !p%TwSSDOF2  = InputFileData%TwSSDOF2
   !p%PtfmSgDOF = InputFileData%PtfmSgDOF
   !p%PtfmSwDOF = InputFileData%PtfmSwDOF
   !p%PtfmHvDOF = InputFileData%PtfmHvDOF
   !p%PtfmRDOF  = InputFileData%PtfmRDOF
   !p%PtfmPDOF  = InputFileData%PtfmPDOF
   !p%PtfmYDOF  = InputFileData%PtfmYDOF
   !p%Azimuth   = InputFileData%Azimuth
   p%RotSpeed  = InputFileData%RotSpeed
   !p%TTDspFA   = InputFileData%TTDspFA
   !p%TTDspSS   = InputFileData%TTDspSS
   !p%PtfmSurge = InputFileData%PtfmSurge
   !p%PtfmSway  = InputFileData%PtfmSway
   !p%PtfmHeave = InputFileData%PtfmHeave
   !p%PtfmRoll  = InputFileData%PtfmRoll
   !p%PtfmPitch = InputFileData%PtfmPitch
   !p%PtfmYaw   = InputFileData%PtfmYaw
   p%HubCM     = InputFileData%HubCM
   p%AzimB1Up  = InputFileData%AzimB1Up

   p%NacCMxn   = InputFileData%NacCMxn
   p%NacCMyn   = InputFileData%NacCMyn
   p%NacCMzn   = InputFileData%NacCMzn
   !p%NcIMUxn   = InputFileData%NcIMUxn
   !p%NcIMUyn   = InputFileData%NcIMUyn
   !p%NcIMUzn   = InputFileData%NcIMUzn


   ! plus everything else from FAST_Initialize



END SUBROUTINE SetPrimaryParameters
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine initializes the continuous states of the module.
!! It assumes the parameters are set and that InputFileData contains initial conditions for the continuous states.
SUBROUTINE Init_ContStates( x, p, InputFileData, OtherState, ErrStat, ErrMsg  )
!..................................................................................................................................
   TYPE(ED_ContinuousStateType), INTENT(OUT)    :: x                 !< Initial continuous states
   TYPE(ED_ParameterType),       INTENT(IN)     :: p                 !< Parameters of the structural dynamics module
   TYPE(ED_InputFile),           INTENT(IN)     :: InputFileData     !< Data stored in the module's input file
   TYPE(ED_OtherStateType),      INTENT(IN)     :: OtherState        !< Initial other states
   INTEGER(IntKi),               INTENT(OUT)    :: ErrStat           !< Error status
   CHARACTER(*),                 INTENT(OUT)    :: ErrMsg            !< Error message

      ! local variables
   REAL(ReKi)                                   :: InitQE1(p%NumBl)  ! Initial value of the 1st blade edge DOF
   REAL(ReKi)                                   :: InitQF1(p%NumBl)  ! Initial value of the 1st blade flap DOF
   REAL(ReKi)                                   :: InitQF2(p%NumBl)  ! Initial value of the 2nd blade flap DOF
!   INTEGER(IntKi)                               :: I                 ! loop counter

      
      ! First allocate the arrays stored here:

   CALL AllocAry( x%QT, p%NDOF,   'QT',   ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN

   CALL AllocAry( x%QDT, p%NDOF,  'QDT',  ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   
   
      ! Calculate/apply the initial blade DOF values to the corresponding DOFs.
   IF (.NOT. p%BD4Blades) THEN  !Skipping subroutine if BeamDyn = TRUE
      CALL InitBlDefl ( p, InputFileData, InitQF1, InitQF2, InitQE1, ErrStat, ErrMsg  )
   ELSE
      InitQF1 = 0.0_ReKi
      InitQF2 = 0.0_ReKi
      InitQE1 = 0.0_ReKi
   END IF
   
      
   x%QT ( DOF_BF(1:p%NumBl,1) ) = InitQF1   ! These come from InitBlDefl().
   x%QT ( DOF_BF(1:p%NumBl,2) ) = InitQF2   ! These come from InitBlDefl().
   x%QT ( DOF_BE(1:p%NumBl,1) ) = InitQE1   ! These come from InitBlDefl().
   x%QDT( DOF_BF(1:p%NumBl,1) ) = 0.0
   x%QDT( DOF_BF(1:p%NumBl,2) ) = 0.0
   x%QDT( DOF_BE(1:p%NumBl,1) ) = 0.0

      ! Teeter Motion

   IF ( p%NumBl == 2 )  THEN !note, DOF_Teet doesn't exist for 3-bladed turbine, so don't include an ELSE here

      ! Set initial teeter angle to TeetDefl and initial teeter angular velocity to 0.

      x%QT (DOF_Teet) = InputFileData%TeetDefl
      x%QDT(DOF_Teet) = 0.0
   ENDIF

      ! Generator azimuth

      ! Set initial generator azimuth angle.  Turn rotor on, whether it is
      !   fixed or variable speed.  If it is fixed speed, set up the
      !   fixed rpm.

   !JASON: CHANGE THESE MOD() FUNCTIONS INTO MODULO() FUNCTIONS SO THAT YOU CAN ELIMINATE ADDING 360:
!   x%QT (DOF_GeAz) = MOD( (InputFileData%Azimuth - p%AzimB1Up)*R2D + 270.0 + 360.0, 360.0 )*D2R   ! Internal position of blade 1
   
   x%QT (DOF_GeAz) = REAL(InputFileData%Azimuth, R8Ki) - p%AzimB1Up - REAL(Piby2_D, R8Ki)
   CALL Zero2TwoPi( x%QT (DOF_GeAz) )
   x%QDT(DOF_GeAz) = p%RotSpeed                                               ! Rotor speed in rad/sec.


      ! Shaft compliance

   ! The initial shaft compliance displacements and velocities are all zero.
   !   They will remain zero if the drivetrain DOF is disabled:

   x%QT (DOF_DrTr) = 0.0
   x%QDT(DOF_DrTr) = 0.0




      ! Rotor-furl motion

      ! Set initial rotor-furl angle to RotFurl.  If rotor-furl is off, this
      !   becomes a fixed rotor-furl angle.

   x%QT (DOF_RFrl) = InputFileData%RotFurl
   x%QDT(DOF_RFrl) = 0.0



      ! Tail-furl motion

      ! Set initial tail-furl angle to TailFurl.  If tail-furl is off, this becomes a fixed tail-furl angle.

   x%QT (DOF_TFrl) = InputFileData%TailFurl
   x%QDT(DOF_TFrl) = 0.0



      ! Yaw Motion

      ! Set initial yaw angle to NacYaw.  If yaw is off, this becomes a fixed yaw angle.

   x%QT (DOF_Yaw) = InputFileData%NacYaw
   x%QDT(DOF_Yaw) = 0.0



      ! Tower motion

      ! Assign all the displacements to mode 1 unless it is disabled.  If mode 1
      !   is disabled and mode 2 is enabled, assign all displacements to mode 2.
      ! If both modes are disabled, set the displacements to zero.

   x%QT   (DOF_TFA1) =  0.0
   x%QT   (DOF_TSS1) =  0.0
   x%QT   (DOF_TFA2) =  0.0
   x%QT   (DOF_TSS2) =  0.0

   IF (    InputFileData%TwFADOF1 )  THEN   ! First fore-aft tower mode is enabled.
      x%QT(DOF_TFA1) =  InputFileData%TTDspFA
   ELSEIF( InputFileData%TwFADOF2 )  THEN   ! Second fore-aft tower mode is enabled, but first is not.
      x%QT(DOF_TFA2) =  InputFileData%TTDspFA
   ENDIF

   IF (    InputFileData%TwSSDOF1 )  THEN   ! First side-to-side tower mode is enabled.
      x%QT(DOF_TSS1) = -InputFileData%TTDspSS
   ELSEIF( InputFileData%TwSSDOF2 )  THEN   ! Second side-to-side tower mode is enabled, but first is not.
      x%QT(DOF_TSS2) = -InputFileData%TTDspSS
   ENDIF

   x%QDT  (DOF_TFA1) =  0.0
   x%QDT  (DOF_TSS1) =  0.0
   x%QDT  (DOF_TFA2) =  0.0
   x%QDT  (DOF_TSS2) =  0.0



      ! Platform Motion

      ! Set initial platform displacements.  If platform DOFs are off, these
      !   become fixed platform displacements.

   x%QT (DOF_Sg) = InputFileData%PtfmSurge
   x%QT (DOF_Sw) = InputFileData%PtfmSway
   x%QT (DOF_Hv) = InputFileData%PtfmHeave
   x%QT (DOF_R ) = InputFileData%PtfmRoll
   x%QT (DOF_P ) = InputFileData%PtfmPitch
   x%QT (DOF_Y ) = InputFileData%PtfmYaw
   x%QDT(DOF_Sg) = 0.0
   x%QDT(DOF_Sw) = 0.0
   x%QDT(DOF_Hv) = 0.0
   x%QDT(DOF_R ) = 0.0
   x%QDT(DOF_P ) = 0.0
   x%QDT(DOF_Y ) = 0.0

   


END SUBROUTINE Init_ContStates
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine initializes the other states of the module.
!! It assumes the parameters are set and that InputFileData contains initial conditions for the continuous states. and p%NumBl is set
SUBROUTINE Init_MiscOtherStates( m, OtherState, p, x, InputFileData, ErrStat, ErrMsg  )
!..................................................................................................................................
   TYPE(ED_MiscVarType),         INTENT(OUT)    :: m                 !< Initial misc variables
   TYPE(ED_OtherStateType),      INTENT(OUT)    :: OtherState        !< Initial other states
   TYPE(ED_ParameterType),       INTENT(IN)     :: p                 !< Parameters of the structural dynamics module
   TYPE(ED_ContinuousStateType), INTENT(IN)     :: x                 !< Initial continuous states
   TYPE(ED_InputFile),           INTENT(IN)     :: InputFileData     !< Data stored in the module's input file
   INTEGER(IntKi),               INTENT(OUT)    :: ErrStat           !< Error status
   CHARACTER(*),                 INTENT(OUT)    :: ErrMsg            !< Error message

   INTEGER(IntKi)               :: I                                 ! Generic loop counter.

      ! First allocate the arrays stored here:

   CALL Alloc_RtHS( m%RtHS, p, ErrStat, ErrMsg  )
      IF ( ErrStat >= AbortErrLev ) RETURN

   CALL Alloc_CoordSys( m%CoordSys, p, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) RETURN
   
   CALL AllocAry( m%QD2T, p%NDOF,   'm%QD2T',  ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) RETURN

   
   ALLOCATE ( m%AllOuts(0:MaxOutPts) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error allocating memory for the AllOuts array.'
         RETURN
      ENDIF   
   m%AllOuts = 0.0_ReKi
   
   m%IgnoreMod = .false. ! for general time steps, we don't ignore the modulos in ED_CalcOutput
   
      ! for loose coupling:
   CALL AllocAry( OtherState%IC,  ED_NMX,   'IC',   ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) RETURN
   
   
   
   !CALL AllocAry(OtherState%BlPitch, p%NumBl, 'BlPitch', ErrStat, ErrMsg )
   !      IF ( ErrStat >= AbortErrLev ) RETURN
   !   OtherState%BlPitch = InputFileData%BlPitch(1:p%NumBl)

   CALL AllocAry( m%AugMat,       p%NDOF,          p%NAug,          'AugMat',       ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) RETURN
   CALL AllocAry( m%OgnlGeAzRo,                    p%NAug,          'OgnlGeAzRo',   ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) RETURN
   CALL AllocAry( m%SolnVec,      p%DOFs%NActvDOF,                  'SolnVec',      ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) RETURN
   CALL AllocAry( m%AugMat_pivot, p%DOFs%NActvDOF,                  'AugMat_pivot', ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) RETURN
   CALL AllocAry( m%AugMat_factor,p%DOFs%NActvDOF, p%DOFs%NActvDOF, 'AugMat_factor',ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) RETURN

   
      ! Now initialize the IC array = (/NMX, NMX-1, ... , 1 /)
      ! this keeps track of the position in the array of continuous states (stored in other states)

   OtherState%IC(1) = ED_NMX
   DO I = 2,ED_NMX
      OtherState%IC(I) = OtherState%IC(I-1) - 1
   ENDDO


   !   ! Initialize the accelerations to zero.
   !
   !OtherState%QD2 = 0.0_ReKi


   OtherState%n   = -1  ! we haven't updated OtherState%xdot, yet
   
   DO i = LBOUND(OtherState%xdot,1), UBOUND(OtherState%xdot,1)
      CALL ED_CopyContState( x, OtherState%xdot(i), MESH_NEWCOPY, ErrStat, ErrMsg)
         IF ( ErrStat >= AbortErrLev ) RETURN 
   ENDDO
   
      ! hacks for HSS brake function:
   
   OtherState%HSSBrTrq   = 0.0_ReKi
   OtherState%HSSBrTrqC  = 0.0_ReKi
   OtherState%SgnPrvLSTQ = 1
   OtherState%SgnLSTQ    = 1
   
   
END SUBROUTINE Init_MiscOtherStates

!----------------------------------------------------------------------------------------------------------------------------------

!**********************************************************************************************************************************
! NOTE: The following lines of code were generated by a Matlab script called "Write_ChckOutLst.m"
!      using the parameters listed in the "OutListParameters.xlsx" Excel file. Any changes to these 
!      lines should be modified in the Matlab script and/or Excel worksheet as necessary. 
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine checks to see if any requested output channel names (stored in the OutList(:)) are invalid. It returns a 
!! warning if any of the channels are not available outputs from the module.
!!  It assigns the settings for OutParam(:) (i.e, the index, name, and units of the output channels, WriteOutput(:)).
!!  the sign is set to 0 if the channel is invalid.
!! It sets assumes the value p%NumOuts has been set before this routine has been called, and it sets the values of p%OutParam here.
!! 
!! This routine was generated by Write_ChckOutLst.m using the parameters listed in OutListParameters.xlsx at 25-Jan-2021 13:23:51.
SUBROUTINE SetOutParam(OutList, p, ErrStat, ErrMsg )
!..................................................................................................................................

   IMPLICIT                        NONE

      ! Passed variables

   CHARACTER(ChanLen),        INTENT(IN)     :: OutList(:)                        !< The list out user-requested outputs
   TYPE(ED_ParameterType),    INTENT(INOUT)  :: p                                 !< The module parameters
   INTEGER(IntKi),            INTENT(OUT)    :: ErrStat                           !< The error status code
   CHARACTER(*),              INTENT(OUT)    :: ErrMsg                            !< The error message, if an error occurred

      ! Local variables

   INTEGER                      :: ErrStat2                                        ! temporary (local) error status
   INTEGER                      :: I                                               ! Generic loop-counting index
   INTEGER                      :: J                                               ! Generic loop-counting index
   INTEGER                      :: INDX                                            ! Index for valid arrays
   INTEGER                      :: startIndx                                       ! Index for BeamDyn

   LOGICAL                      :: CheckOutListAgain                               ! Flag used to determine if output parameter starting with "M" is valid (or the negative of another parameter)
   LOGICAL                      :: InvalidOutput(0:MaxOutPts)                      ! This array determines if the output channel is valid for this configuration
   CHARACTER(ChanLen)           :: OutListTmp                                      ! A string to temporarily hold OutList(I)
   CHARACTER(*), PARAMETER      :: RoutineName = "SetOutParam"

      CHARACTER(OutStrLenM1), PARAMETER  :: ValidParamAry(981) =  (/  &   ! This lists the names of the allowed parameters, which must be sorted alphabetically
                                  "AZIMUTH  ","BLDPITCH1","BLDPITCH2","BLDPITCH3","BLPITCH1 ","BLPITCH2 ","BLPITCH3 ","GENACCEL ", &
                                  "GENSPEED ","HSSBRTQ  ","HSSHFTA  ","HSSHFTPWR","HSSHFTTQ ","HSSHFTV  ","IPDEFL1  ","IPDEFL2  ", &
                                  "IPDEFL3  ","LSSGAGA  ","LSSGAGAXA","LSSGAGAXS","LSSGAGFXA","LSSGAGFXS","LSSGAGFYA","LSSGAGFYS", &
                                  "LSSGAGFZA","LSSGAGFZS","LSSGAGMXA","LSSGAGMXS","LSSGAGMYA","LSSGAGMYS","LSSGAGMZA","LSSGAGMZS", &
                                  "LSSGAGP  ","LSSGAGPXA","LSSGAGPXS","LSSGAGV  ","LSSGAGVXA","LSSGAGVXS","LSSHFTFXA","LSSHFTFXS", &
                                  "LSSHFTFYA","LSSHFTFYS","LSSHFTFZA","LSSHFTFZS","LSSHFTMXA","LSSHFTMXS","LSSHFTPWR","LSSHFTTQ ", &
                                  "LSSTIPA  ","LSSTIPAXA","LSSTIPAXS","LSSTIPMYA","LSSTIPMYS","LSSTIPMZA","LSSTIPMZS","LSSTIPP  ", &
                                  "LSSTIPPXA","LSSTIPPXS","LSSTIPV  ","LSSTIPVXA","LSSTIPVXS","NACYAW   ","NACYAWA  ","NACYAWP  ", &
                                  "NACYAWV  ","NCIMURAXS","NCIMURAYS","NCIMURAZS","NCIMURVXS","NCIMURVYS","NCIMURVZS","NCIMUTAXS", &
                                  "NCIMUTAYS","NCIMUTAZS","NCIMUTVXS","NCIMUTVYS","NCIMUTVZS","OOPDEFL1 ","OOPDEFL2 ","OOPDEFL3 ", &
                                  "PTCHDEFL1","PTCHDEFL2","PTCHDEFL3","PTCHPMZB1","PTCHPMZB2","PTCHPMZB3","PTCHPMZC1","PTCHPMZC2", &
                                  "PTCHPMZC3","PTFMHEAVE","PTFMPITCH","PTFMRAXI ","PTFMRAXT ","PTFMRAYI ","PTFMRAYT ","PTFMRAZI ", &
                                  "PTFMRAZT ","PTFMRDXI ","PTFMRDYI ","PTFMRDZI ","PTFMROLL ","PTFMRVXI ","PTFMRVXT ","PTFMRVYI ", &
                                  "PTFMRVYT ","PTFMRVZI ","PTFMRVZT ","PTFMSURGE","PTFMSWAY ","PTFMTAXI ","PTFMTAXT ","PTFMTAYI ", &
                                  "PTFMTAYT ","PTFMTAZI ","PTFMTAZT ","PTFMTDXI ","PTFMTDXT ","PTFMTDYI ","PTFMTDYT ","PTFMTDZI ", &
                                  "PTFMTDZT ","PTFMTVXI ","PTFMTVXT ","PTFMTVYI ","PTFMTVYT ","PTFMTVZI ","PTFMTVZT ","PTFMYAW  ", &
                                  "QD2_B1E1 ","QD2_B1F1 ","QD2_B1F2 ","QD2_B2E1 ","QD2_B2F1 ","QD2_B2F2 ","QD2_B3E1 ","QD2_B3F1 ", &
                                  "QD2_B3F2 ","QD2_DRTR ","QD2_GEAZ ","QD2_HV   ","QD2_P    ","QD2_R    ","QD2_RFRL ","QD2_SG   ", &
                                  "QD2_SW   ","QD2_TEET ","QD2_TFA1 ","QD2_TFA2 ","QD2_TFRL ","QD2_TSS1 ","QD2_TSS2 ","QD2_Y    ", &
                                  "QD2_YAW  ","QD_B1E1  ","QD_B1F1  ","QD_B1F2  ","QD_B2E1  ","QD_B2F1  ","QD_B2F2  ","QD_B3E1  ", &
                                  "QD_B3F1  ","QD_B3F2  ","QD_DRTR  ","QD_GEAZ  ","QD_HV    ","QD_P     ","QD_R     ","QD_RFRL  ", &
                                  "QD_SG    ","QD_SW    ","QD_TEET  ","QD_TFA1  ","QD_TFA2  ","QD_TFRL  ","QD_TSS1  ","QD_TSS2  ", &
                                  "QD_Y     ","QD_YAW   ","Q_B1E1   ","Q_B1F1   ","Q_B1F2   ","Q_B2E1   ","Q_B2F1   ","Q_B2F2   ", &
                                  "Q_B3E1   ","Q_B3F1   ","Q_B3F2   ","Q_DRTR   ","Q_GEAZ   ","Q_HV     ","Q_P      ","Q_R      ", &
                                  "Q_RFRL   ","Q_SG     ","Q_SW     ","Q_TEET   ","Q_TFA1   ","Q_TFA2   ","Q_TFRL   ","Q_TSS1   ", &
                                  "Q_TSS2   ","Q_Y      ","Q_YAW    ","RFRLBRM  ","ROLLDEFL1","ROLLDEFL2","ROLLDEFL3","ROOTFXB1 ", &
                                  "ROOTFXB2 ","ROOTFXB3 ","ROOTFXC1 ","ROOTFXC2 ","ROOTFXC3 ","ROOTFYB1 ","ROOTFYB2 ","ROOTFYB3 ", &
                                  "ROOTFYC1 ","ROOTFYC2 ","ROOTFYC3 ","ROOTFZB1 ","ROOTFZB2 ","ROOTFZB3 ","ROOTFZC1 ","ROOTFZC2 ", &
                                  "ROOTFZC3 ","ROOTMEDG1","ROOTMEDG2","ROOTMEDG3","ROOTMFLP1","ROOTMFLP2","ROOTMFLP3","ROOTMIP1 ", &
                                  "ROOTMIP2 ","ROOTMIP3 ","ROOTMOOP1","ROOTMOOP2","ROOTMOOP3","ROOTMXB1 ","ROOTMXB2 ","ROOTMXB3 ", &
                                  "ROOTMXC1 ","ROOTMXC2 ","ROOTMXC3 ","ROOTMYB1 ","ROOTMYB2 ","ROOTMYB3 ","ROOTMYC1 ","ROOTMYC2 ", &
                                  "ROOTMYC3 ","ROOTMZB1 ","ROOTMZB2 ","ROOTMZB3 ","ROOTMZC1 ","ROOTMZC2 ","ROOTMZC3 ","ROTACCEL ", &
                                  "ROTFURL  ","ROTFURLA ","ROTFURLP ","ROTFURLV ","ROTPWR   ","ROTSPEED ","ROTTEETA ","ROTTEETP ", &
                                  "ROTTEETV ","ROTTHRUST","ROTTORQ  ","SPN1ALXB1","SPN1ALXB2","SPN1ALXB3","SPN1ALYB1","SPN1ALYB2", &
                                  "SPN1ALYB3","SPN1ALZB1","SPN1ALZB2","SPN1ALZB3","SPN1FLXB1","SPN1FLXB2","SPN1FLXB3","SPN1FLYB1", &
                                  "SPN1FLYB2","SPN1FLYB3","SPN1FLZB1","SPN1FLZB2","SPN1FLZB3","SPN1MLXB1","SPN1MLXB2","SPN1MLXB3", &
                                  "SPN1MLYB1","SPN1MLYB2","SPN1MLYB3","SPN1MLZB1","SPN1MLZB2","SPN1MLZB3","SPN1RDXB1","SPN1RDXB2", &
                                  "SPN1RDXB3","SPN1RDYB1","SPN1RDYB2","SPN1RDYB3","SPN1RDZB1","SPN1RDZB2","SPN1RDZB3","SPN1TDXB1", &
                                  "SPN1TDXB2","SPN1TDXB3","SPN1TDYB1","SPN1TDYB2","SPN1TDYB3","SPN1TDZB1","SPN1TDZB2","SPN1TDZB3", &
                                  "SPN2ALXB1","SPN2ALXB2","SPN2ALXB3","SPN2ALYB1","SPN2ALYB2","SPN2ALYB3","SPN2ALZB1","SPN2ALZB2", &
                                  "SPN2ALZB3","SPN2FLXB1","SPN2FLXB2","SPN2FLXB3","SPN2FLYB1","SPN2FLYB2","SPN2FLYB3","SPN2FLZB1", &
                                  "SPN2FLZB2","SPN2FLZB3","SPN2MLXB1","SPN2MLXB2","SPN2MLXB3","SPN2MLYB1","SPN2MLYB2","SPN2MLYB3", &
                                  "SPN2MLZB1","SPN2MLZB2","SPN2MLZB3","SPN2RDXB1","SPN2RDXB2","SPN2RDXB3","SPN2RDYB1","SPN2RDYB2", &
                                  "SPN2RDYB3","SPN2RDZB1","SPN2RDZB2","SPN2RDZB3","SPN2TDXB1","SPN2TDXB2","SPN2TDXB3","SPN2TDYB1", &
                                  "SPN2TDYB2","SPN2TDYB3","SPN2TDZB1","SPN2TDZB2","SPN2TDZB3","SPN3ALXB1","SPN3ALXB2","SPN3ALXB3", &
                                  "SPN3ALYB1","SPN3ALYB2","SPN3ALYB3","SPN3ALZB1","SPN3ALZB2","SPN3ALZB3","SPN3FLXB1","SPN3FLXB2", &
                                  "SPN3FLXB3","SPN3FLYB1","SPN3FLYB2","SPN3FLYB3","SPN3FLZB1","SPN3FLZB2","SPN3FLZB3","SPN3MLXB1", &
                                  "SPN3MLXB2","SPN3MLXB3","SPN3MLYB1","SPN3MLYB2","SPN3MLYB3","SPN3MLZB1","SPN3MLZB2","SPN3MLZB3", &
                                  "SPN3RDXB1","SPN3RDXB2","SPN3RDXB3","SPN3RDYB1","SPN3RDYB2","SPN3RDYB3","SPN3RDZB1","SPN3RDZB2", &
                                  "SPN3RDZB3","SPN3TDXB1","SPN3TDXB2","SPN3TDXB3","SPN3TDYB1","SPN3TDYB2","SPN3TDYB3","SPN3TDZB1", &
                                  "SPN3TDZB2","SPN3TDZB3","SPN4ALXB1","SPN4ALXB2","SPN4ALXB3","SPN4ALYB1","SPN4ALYB2","SPN4ALYB3", &
                                  "SPN4ALZB1","SPN4ALZB2","SPN4ALZB3","SPN4FLXB1","SPN4FLXB2","SPN4FLXB3","SPN4FLYB1","SPN4FLYB2", &
                                  "SPN4FLYB3","SPN4FLZB1","SPN4FLZB2","SPN4FLZB3","SPN4MLXB1","SPN4MLXB2","SPN4MLXB3","SPN4MLYB1", &
                                  "SPN4MLYB2","SPN4MLYB3","SPN4MLZB1","SPN4MLZB2","SPN4MLZB3","SPN4RDXB1","SPN4RDXB2","SPN4RDXB3", &
                                  "SPN4RDYB1","SPN4RDYB2","SPN4RDYB3","SPN4RDZB1","SPN4RDZB2","SPN4RDZB3","SPN4TDXB1","SPN4TDXB2", &
                                  "SPN4TDXB3","SPN4TDYB1","SPN4TDYB2","SPN4TDYB3","SPN4TDZB1","SPN4TDZB2","SPN4TDZB3","SPN5ALXB1", &
                                  "SPN5ALXB2","SPN5ALXB3","SPN5ALYB1","SPN5ALYB2","SPN5ALYB3","SPN5ALZB1","SPN5ALZB2","SPN5ALZB3", &
                                  "SPN5FLXB1","SPN5FLXB2","SPN5FLXB3","SPN5FLYB1","SPN5FLYB2","SPN5FLYB3","SPN5FLZB1","SPN5FLZB2", &
                                  "SPN5FLZB3","SPN5MLXB1","SPN5MLXB2","SPN5MLXB3","SPN5MLYB1","SPN5MLYB2","SPN5MLYB3","SPN5MLZB1", &
                                  "SPN5MLZB2","SPN5MLZB3","SPN5RDXB1","SPN5RDXB2","SPN5RDXB3","SPN5RDYB1","SPN5RDYB2","SPN5RDYB3", &
                                  "SPN5RDZB1","SPN5RDZB2","SPN5RDZB3","SPN5TDXB1","SPN5TDXB2","SPN5TDXB3","SPN5TDYB1","SPN5TDYB2", &
                                  "SPN5TDYB3","SPN5TDZB1","SPN5TDZB2","SPN5TDZB3","SPN6ALXB1","SPN6ALXB2","SPN6ALXB3","SPN6ALYB1", &
                                  "SPN6ALYB2","SPN6ALYB3","SPN6ALZB1","SPN6ALZB2","SPN6ALZB3","SPN6FLXB1","SPN6FLXB2","SPN6FLXB3", &
                                  "SPN6FLYB1","SPN6FLYB2","SPN6FLYB3","SPN6FLZB1","SPN6FLZB2","SPN6FLZB3","SPN6MLXB1","SPN6MLXB2", &
                                  "SPN6MLXB3","SPN6MLYB1","SPN6MLYB2","SPN6MLYB3","SPN6MLZB1","SPN6MLZB2","SPN6MLZB3","SPN6RDXB1", &
                                  "SPN6RDXB2","SPN6RDXB3","SPN6RDYB1","SPN6RDYB2","SPN6RDYB3","SPN6RDZB1","SPN6RDZB2","SPN6RDZB3", &
                                  "SPN6TDXB1","SPN6TDXB2","SPN6TDXB3","SPN6TDYB1","SPN6TDYB2","SPN6TDYB3","SPN6TDZB1","SPN6TDZB2", &
                                  "SPN6TDZB3","SPN7ALXB1","SPN7ALXB2","SPN7ALXB3","SPN7ALYB1","SPN7ALYB2","SPN7ALYB3","SPN7ALZB1", &
                                  "SPN7ALZB2","SPN7ALZB3","SPN7FLXB1","SPN7FLXB2","SPN7FLXB3","SPN7FLYB1","SPN7FLYB2","SPN7FLYB3", &
                                  "SPN7FLZB1","SPN7FLZB2","SPN7FLZB3","SPN7MLXB1","SPN7MLXB2","SPN7MLXB3","SPN7MLYB1","SPN7MLYB2", &
                                  "SPN7MLYB3","SPN7MLZB1","SPN7MLZB2","SPN7MLZB3","SPN7RDXB1","SPN7RDXB2","SPN7RDXB3","SPN7RDYB1", &
                                  "SPN7RDYB2","SPN7RDYB3","SPN7RDZB1","SPN7RDZB2","SPN7RDZB3","SPN7TDXB1","SPN7TDXB2","SPN7TDXB3", &
                                  "SPN7TDYB1","SPN7TDYB2","SPN7TDYB3","SPN7TDZB1","SPN7TDZB2","SPN7TDZB3","SPN8ALXB1","SPN8ALXB2", &
                                  "SPN8ALXB3","SPN8ALYB1","SPN8ALYB2","SPN8ALYB3","SPN8ALZB1","SPN8ALZB2","SPN8ALZB3","SPN8FLXB1", &
                                  "SPN8FLXB2","SPN8FLXB3","SPN8FLYB1","SPN8FLYB2","SPN8FLYB3","SPN8FLZB1","SPN8FLZB2","SPN8FLZB3", &
                                  "SPN8MLXB1","SPN8MLXB2","SPN8MLXB3","SPN8MLYB1","SPN8MLYB2","SPN8MLYB3","SPN8MLZB1","SPN8MLZB2", &
                                  "SPN8MLZB3","SPN8RDXB1","SPN8RDXB2","SPN8RDXB3","SPN8RDYB1","SPN8RDYB2","SPN8RDYB3","SPN8RDZB1", &
                                  "SPN8RDZB2","SPN8RDZB3","SPN8TDXB1","SPN8TDXB2","SPN8TDXB3","SPN8TDYB1","SPN8TDYB2","SPN8TDYB3", &
                                  "SPN8TDZB1","SPN8TDZB2","SPN8TDZB3","SPN9ALXB1","SPN9ALXB2","SPN9ALXB3","SPN9ALYB1","SPN9ALYB2", &
                                  "SPN9ALYB3","SPN9ALZB1","SPN9ALZB2","SPN9ALZB3","SPN9FLXB1","SPN9FLXB2","SPN9FLXB3","SPN9FLYB1", &
                                  "SPN9FLYB2","SPN9FLYB3","SPN9FLZB1","SPN9FLZB2","SPN9FLZB3","SPN9MLXB1","SPN9MLXB2","SPN9MLXB3", &
                                  "SPN9MLYB1","SPN9MLYB2","SPN9MLYB3","SPN9MLZB1","SPN9MLZB2","SPN9MLZB3","SPN9RDXB1","SPN9RDXB2", &
                                  "SPN9RDXB3","SPN9RDYB1","SPN9RDYB2","SPN9RDYB3","SPN9RDZB1","SPN9RDZB2","SPN9RDZB3","SPN9TDXB1", &
                                  "SPN9TDXB2","SPN9TDXB3","SPN9TDYB1","SPN9TDYB2","SPN9TDYB3","SPN9TDZB1","SPN9TDZB2","SPN9TDZB3", &
                                  "TAILFURL ","TAILFURLA","TAILFURLP","TAILFURLV","TEETAYA  ","TEETDEFL ","TEETPYA  ","TEETVYA  ", &
                                  "TFRLBRM  ","TIP2TWR1 ","TIP2TWR2 ","TIP2TWR3 ","TIPALXB1 ","TIPALXB2 ","TIPALXB3 ","TIPALYB1 ", &
                                  "TIPALYB2 ","TIPALYB3 ","TIPALZB1 ","TIPALZB2 ","TIPALZB3 ","TIPCLRNC1","TIPCLRNC2","TIPCLRNC3", &
                                  "TIPDXB1  ","TIPDXB2  ","TIPDXB3  ","TIPDXC1  ","TIPDXC2  ","TIPDXC3  ","TIPDYB1  ","TIPDYB2  ", &
                                  "TIPDYB3  ","TIPDYC1  ","TIPDYC2  ","TIPDYC3  ","TIPDZB1  ","TIPDZB2  ","TIPDZB3  ","TIPDZC1  ", &
                                  "TIPDZC2  ","TIPDZC3  ","TIPRDXB1 ","TIPRDXB2 ","TIPRDXB3 ","TIPRDYB1 ","TIPRDYB2 ","TIPRDYB3 ", &
                                  "TIPRDZB1 ","TIPRDZB2 ","TIPRDZB3 ","TIPRDZC1 ","TIPRDZC2 ","TIPRDZC3 ","TTDSPAX  ","TTDSPFA  ", &
                                  "TTDSPPTCH","TTDSPROLL","TTDSPSS  ","TTDSPTWST","TWHT1ALXT","TWHT1ALYT","TWHT1ALZT","TWHT1FLXT", &
                                  "TWHT1FLYT","TWHT1FLZT","TWHT1MLXT","TWHT1MLYT","TWHT1MLZT","TWHT1RDXT","TWHT1RDYT","TWHT1RDZT", &
                                  "TWHT1RPXI","TWHT1RPYI","TWHT1RPZI","TWHT1TDXT","TWHT1TDYT","TWHT1TDZT","TWHT1TPXI","TWHT1TPYI", &
                                  "TWHT1TPZI","TWHT2ALXT","TWHT2ALYT","TWHT2ALZT","TWHT2FLXT","TWHT2FLYT","TWHT2FLZT","TWHT2MLXT", &
                                  "TWHT2MLYT","TWHT2MLZT","TWHT2RDXT","TWHT2RDYT","TWHT2RDZT","TWHT2RPXI","TWHT2RPYI","TWHT2RPZI", &
                                  "TWHT2TDXT","TWHT2TDYT","TWHT2TDZT","TWHT2TPXI","TWHT2TPYI","TWHT2TPZI","TWHT3ALXT","TWHT3ALYT", &
                                  "TWHT3ALZT","TWHT3FLXT","TWHT3FLYT","TWHT3FLZT","TWHT3MLXT","TWHT3MLYT","TWHT3MLZT","TWHT3RDXT", &
                                  "TWHT3RDYT","TWHT3RDZT","TWHT3RPXI","TWHT3RPYI","TWHT3RPZI","TWHT3TDXT","TWHT3TDYT","TWHT3TDZT", &
                                  "TWHT3TPXI","TWHT3TPYI","TWHT3TPZI","TWHT4ALXT","TWHT4ALYT","TWHT4ALZT","TWHT4FLXT","TWHT4FLYT", &
                                  "TWHT4FLZT","TWHT4MLXT","TWHT4MLYT","TWHT4MLZT","TWHT4RDXT","TWHT4RDYT","TWHT4RDZT","TWHT4RPXI", &
                                  "TWHT4RPYI","TWHT4RPZI","TWHT4TDXT","TWHT4TDYT","TWHT4TDZT","TWHT4TPXI","TWHT4TPYI","TWHT4TPZI", &
                                  "TWHT5ALXT","TWHT5ALYT","TWHT5ALZT","TWHT5FLXT","TWHT5FLYT","TWHT5FLZT","TWHT5MLXT","TWHT5MLYT", &
                                  "TWHT5MLZT","TWHT5RDXT","TWHT5RDYT","TWHT5RDZT","TWHT5RPXI","TWHT5RPYI","TWHT5RPZI","TWHT5TDXT", &
                                  "TWHT5TDYT","TWHT5TDZT","TWHT5TPXI","TWHT5TPYI","TWHT5TPZI","TWHT6ALXT","TWHT6ALYT","TWHT6ALZT", &
                                  "TWHT6FLXT","TWHT6FLYT","TWHT6FLZT","TWHT6MLXT","TWHT6MLYT","TWHT6MLZT","TWHT6RDXT","TWHT6RDYT", &
                                  "TWHT6RDZT","TWHT6RPXI","TWHT6RPYI","TWHT6RPZI","TWHT6TDXT","TWHT6TDYT","TWHT6TDZT","TWHT6TPXI", &
                                  "TWHT6TPYI","TWHT6TPZI","TWHT7ALXT","TWHT7ALYT","TWHT7ALZT","TWHT7FLXT","TWHT7FLYT","TWHT7FLZT", &
                                  "TWHT7MLXT","TWHT7MLYT","TWHT7MLZT","TWHT7RDXT","TWHT7RDYT","TWHT7RDZT","TWHT7RPXI","TWHT7RPYI", &
                                  "TWHT7RPZI","TWHT7TDXT","TWHT7TDYT","TWHT7TDZT","TWHT7TPXI","TWHT7TPYI","TWHT7TPZI","TWHT8ALXT", &
                                  "TWHT8ALYT","TWHT8ALZT","TWHT8FLXT","TWHT8FLYT","TWHT8FLZT","TWHT8MLXT","TWHT8MLYT","TWHT8MLZT", &
                                  "TWHT8RDXT","TWHT8RDYT","TWHT8RDZT","TWHT8RPXI","TWHT8RPYI","TWHT8RPZI","TWHT8TDXT","TWHT8TDYT", &
                                  "TWHT8TDZT","TWHT8TPXI","TWHT8TPYI","TWHT8TPZI","TWHT9ALXT","TWHT9ALYT","TWHT9ALZT","TWHT9FLXT", &
                                  "TWHT9FLYT","TWHT9FLZT","TWHT9MLXT","TWHT9MLYT","TWHT9MLZT","TWHT9RDXT","TWHT9RDYT","TWHT9RDZT", &
                                  "TWHT9RPXI","TWHT9RPYI","TWHT9RPZI","TWHT9TDXT","TWHT9TDYT","TWHT9TDZT","TWHT9TPXI","TWHT9TPYI", &
                                  "TWHT9TPZI","TWRBSFXT ","TWRBSFYT ","TWRBSFZT ","TWRBSMXT ","TWRBSMYT ","TWRBSMZT ","TWRCLRNC1", &
                                  "TWRCLRNC2","TWRCLRNC3","TWRTPTDXI","TWRTPTDYI","TWRTPTDZI","TWSTDEFL1","TWSTDEFL2","TWSTDEFL3", &
                                  "YAWACCEL ","YAWAZN   ","YAWAZP   ","YAWBRFXN ","YAWBRFXP ","YAWBRFYN ","YAWBRFYP ","YAWBRFZN ", &
                                  "YAWBRFZP ","YAWBRMXN ","YAWBRMXP ","YAWBRMYN ","YAWBRMYP ","YAWBRMZN ","YAWBRMZP ","YAWBRRAXP", &
                                  "YAWBRRAYP","YAWBRRAZP","YAWBRRDXT","YAWBRRDYT","YAWBRRDZT","YAWBRRVXP","YAWBRRVYP","YAWBRRVZP", &
                                  "YAWBRTAXP","YAWBRTAYP","YAWBRTAZP","YAWBRTDXI","YAWBRTDXP","YAWBRTDXT","YAWBRTDYI","YAWBRTDYP", &
                                  "YAWBRTDYT","YAWBRTDZI","YAWBRTDZP","YAWBRTDZT","YAWBRTVXP","YAWBRTVYP","YAWBRTVZP","YAWPOS   ", &
                               "YAWPZN   ","YAWPZP   ","YAWRATE  ","YAWVZN   ","YAWVZP   "/)
      INTEGER(IntKi), PARAMETER :: ParamIndxAry(981) =  (/ &                            ! This lists the index into AllOuts(:) of the allowed parameters ValidParamAry(:)
                                   LSSTipPxa , PtchPMzc1 , PtchPMzc2 , PtchPMzc3 , PtchPMzc1 , PtchPMzc2 , PtchPMzc3 ,   HSShftA , &
                                     HSShftV ,   HSSBrTq ,   HSShftA , HSShftPwr ,  HSShftTq ,   HSShftV ,   TipDyc1 ,   TipDyc2 , &
                                     TipDyc3 , LSSGagAxa , LSSGagAxa , LSSGagAxa , LSShftFxa , LSShftFxa , LSShftFya , LSShftFys , &
                                   LSShftFza , LSShftFzs , LSShftMxa , LSShftMxa , LSSGagMya , LSSGagMys , LSSGagMza , LSSGagMzs , &
                                   LSSGagPxa , LSSGagPxa , LSSGagPxa , LSSGagVxa , LSSGagVxa , LSSGagVxa , LSShftFxa , LSShftFxa , &
                                   LSShftFya , LSShftFys , LSShftFza , LSShftFzs , LSShftMxa , LSShftMxa ,    RotPwr , LSShftMxa , &
                                   LSSTipAxa , LSSTipAxa , LSSTipAxa , LSSTipMya , LSSTipMys , LSSTipMza , LSSTipMzs , LSSTipPxa , &
                                   LSSTipPxa , LSSTipPxa , LSSTipVxa , LSSTipVxa , LSSTipVxa ,    YawPzn ,    YawAzn ,    YawPzn , &
                                      YawVzn , NcIMURAxs , NcIMURAys , NcIMURAzs , NcIMURVxs , NcIMURVys , NcIMURVzs , NcIMUTAxs , &
                                   NcIMUTAys , NcIMUTAzs , NcIMUTVxs , NcIMUTVys , NcIMUTVzs ,   TipDxc1 ,   TipDxc2 ,   TipDxc3 , &
                                    TipRDyb1 ,  TipRDyb2 ,  TipRDyb3 , PtchPMzc1 , PtchPMzc2 , PtchPMzc3 , PtchPMzc1 , PtchPMzc2 , &
                                   PtchPMzc3 ,  PtfmTDzi ,  PtfmRDyi ,  PtfmRAxi ,  PtfmRAxt ,  PtfmRAyi ,  PtfmRAyt ,  PtfmRAzi , &
                                    PtfmRAzt ,  PtfmRDxi ,  PtfmRDyi ,  PtfmRDzi ,  PtfmRDxi ,  PtfmRVxi ,  PtfmRVxt ,  PtfmRVyi , &
                                    PtfmRVyt ,  PtfmRVzi ,  PtfmRVzt ,  PtfmTDxi ,  PtfmTDyi ,  PtfmTAxi ,  PtfmTAxt ,  PtfmTAyi , &
                                    PtfmTAyt ,  PtfmTAzi ,  PtfmTAzt ,  PtfmTDxi ,  PtfmTDxt ,  PtfmTDyi ,  PtfmTDyt ,  PtfmTDzi , &
                                    PtfmTDzt ,  PtfmTVxi ,  PtfmTVxt ,  PtfmTVyi ,  PtfmTVyt ,  PtfmTVzi ,  PtfmTVzt ,  PtfmRDzi , &
                                    QD2_B1E1 ,  QD2_B1F1 ,  QD2_B1F2 ,  QD2_B2E1 ,  QD2_B2F1 ,  QD2_B2F2 ,  QD2_B3E1 ,  QD2_B3F1 , &
                                    QD2_B3F2 ,  QD2_DrTr ,  QD2_GeAz ,    QD2_Hv ,     QD2_P ,     QD2_R ,  QD2_RFrl ,    QD2_Sg , &
                                      QD2_Sw ,  QD2_Teet ,  QD2_TFA1 ,  QD2_TFA2 ,  QD2_TFrl ,  QD2_TSS1 ,  QD2_TSS2 ,     QD2_Y , &
                                     QD2_Yaw ,   QD_B1E1 ,   QD_B1F1 ,   QD_B1F2 ,   QD_B2E1 ,   QD_B2F1 ,   QD_B2F2 ,   QD_B3E1 , &
                                     QD_B3F1 ,   QD_B3F2 ,   QD_DrTr ,   QD_GeAz ,     QD_Hv ,      QD_P ,      QD_R ,   QD_RFrl , &
                                       QD_Sg ,     QD_Sw ,   QD_Teet ,   QD_TFA1 ,   QD_TFA2 ,   QD_TFrl ,   QD_TSS1 ,   QD_TSS2 , &
                                        QD_Y ,    QD_Yaw ,    Q_B1E1 ,    Q_B1F1 ,    Q_B1F2 ,    Q_B2E1 ,    Q_B2F1 ,    Q_B2F2 , &
                                      Q_B3E1 ,    Q_B3F1 ,    Q_B3F2 ,    Q_DrTr ,    Q_GeAz ,      Q_Hv ,       Q_P ,       Q_R , &
                                      Q_RFrl ,      Q_Sg ,      Q_Sw ,    Q_Teet ,    Q_TFA1 ,    Q_TFA2 ,    Q_TFrl ,    Q_TSS1 , &
                                      Q_TSS2 ,       Q_Y ,     Q_Yaw ,   RFrlBrM ,  TipRDxb1 ,  TipRDxb2 ,  TipRDxb3 ,  RootFxb1 , &
                                    RootFxb2 ,  RootFxb3 ,  RootFxc1 ,  RootFxc2 ,  RootFxc3 ,  RootFyb1 ,  RootFyb2 ,  RootFyb3 , &
                                    RootFyc1 ,  RootFyc2 ,  RootFyc3 ,  RootFzc1 ,  RootFzc2 ,  RootFzc3 ,  RootFzc1 ,  RootFzc2 , &
                                    RootFzc3 ,  RootMxb1 ,  RootMxb2 ,  RootMxb3 ,  RootMyb1 ,  RootMyb2 ,  RootMyb3 ,  RootMxc1 , &
                                    RootMxc2 ,  RootMxc3 ,  RootMyc1 ,  RootMyc2 ,  RootMyc3 ,  RootMxb1 ,  RootMxb2 ,  RootMxb3 , &
                                    RootMxc1 ,  RootMxc2 ,  RootMxc3 ,  RootMyb1 ,  RootMyb2 ,  RootMyb3 ,  RootMyc1 ,  RootMyc2 , &
                                    RootMyc3 ,  RootMzc1 ,  RootMzc2 ,  RootMzc3 ,  RootMzc1 ,  RootMzc2 ,  RootMzc3 , LSSTipAxa , &
                                    RotFurlP ,  RotFurlA ,  RotFurlP ,  RotFurlV ,    RotPwr , LSSTipVxa ,   TeetAya ,   TeetPya , &
                                     TeetVya , LSShftFxa , LSShftMxa , Spn1ALxb1 , Spn1ALxb2 , Spn1ALxb3 , Spn1ALyb1 , Spn1ALyb2 , &
                                   Spn1ALyb3 , Spn1ALzb1 , Spn1ALzb2 , Spn1ALzb3 , Spn1FLxb1 , Spn1FLxb2 , Spn1FLxb3 , Spn1FLyb1 , &
                                   Spn1FLyb2 , Spn1FLyb3 , Spn1FLzb1 , Spn1FLzb2 , Spn1FLzb3 , Spn1MLxb1 , Spn1MLxb2 , Spn1MLxb3 , &
                                   Spn1MLyb1 , Spn1MLyb2 , Spn1MLyb3 , Spn1MLzb1 , Spn1MLzb2 , Spn1MLzb3 , Spn1RDxb1 , Spn1RDxb2 , &
                                   Spn1RDxb3 , Spn1RDyb1 , Spn1RDyb2 , Spn1RDyb3 , Spn1RDzb1 , Spn1RDzb2 , Spn1RDzb3 , Spn1TDxb1 , &
                                   Spn1TDxb2 , Spn1TDxb3 , Spn1TDyb1 , Spn1TDyb2 , Spn1TDyb3 , Spn1TDzb1 , Spn1TDzb2 , Spn1TDzb3 , &
                                   Spn2ALxb1 , Spn2ALxb2 , Spn2ALxb3 , Spn2ALyb1 , Spn2ALyb2 , Spn2ALyb3 , Spn2ALzb1 , Spn2ALzb2 , &
                                   Spn2ALzb3 , Spn2FLxb1 , Spn2FLxb2 , Spn2FLxb3 , Spn2FLyb1 , Spn2FLyb2 , Spn2FLyb3 , Spn2FLzb1 , &
                                   Spn2FLzb2 , Spn2FLzb3 , Spn2MLxb1 , Spn2MLxb2 , Spn2MLxb3 , Spn2MLyb1 , Spn2MLyb2 , Spn2MLyb3 , &
                                   Spn2MLzb1 , Spn2MLzb2 , Spn2MLzb3 , Spn2RDxb1 , Spn2RDxb2 , Spn2RDxb3 , Spn2RDyb1 , Spn2RDyb2 , &
                                   Spn2RDyb3 , Spn2RDzb1 , Spn2RDzb2 , Spn2RDzb3 , Spn2TDxb1 , Spn2TDxb2 , Spn2TDxb3 , Spn2TDyb1 , &
                                   Spn2TDyb2 , Spn2TDyb3 , Spn2TDzb1 , Spn2TDzb2 , Spn2TDzb3 , Spn3ALxb1 , Spn3ALxb2 , Spn3ALxb3 , &
                                   Spn3ALyb1 , Spn3ALyb2 , Spn3ALyb3 , Spn3ALzb1 , Spn3ALzb2 , Spn3ALzb3 , Spn3FLxb1 , Spn3FLxb2 , &
                                   Spn3FLxb3 , Spn3FLyb1 , Spn3FLyb2 , Spn3FLyb3 , Spn3FLzb1 , Spn3FLzb2 , Spn3FLzb3 , Spn3MLxb1 , &
                                   Spn3MLxb2 , Spn3MLxb3 , Spn3MLyb1 , Spn3MLyb2 , Spn3MLyb3 , Spn3MLzb1 , Spn3MLzb2 , Spn3MLzb3 , &
                                   Spn3RDxb1 , Spn3RDxb2 , Spn3RDxb3 , Spn3RDyb1 , Spn3RDyb2 , Spn3RDyb3 , Spn3RDzb1 , Spn3RDzb2 , &
                                   Spn3RDzb3 , Spn3TDxb1 , Spn3TDxb2 , Spn3TDxb3 , Spn3TDyb1 , Spn3TDyb2 , Spn3TDyb3 , Spn3TDzb1 , &
                                   Spn3TDzb2 , Spn3TDzb3 , Spn4ALxb1 , Spn4ALxb2 , Spn4ALxb3 , Spn4ALyb1 , Spn4ALyb2 , Spn4ALyb3 , &
                                   Spn4ALzb1 , Spn4ALzb2 , Spn4ALzb3 , Spn4FLxb1 , Spn4FLxb2 , Spn4FLxb3 , Spn4FLyb1 , Spn4FLyb2 , &
                                   Spn4FLyb3 , Spn4FLzb1 , Spn4FLzb2 , Spn4FLzb3 , Spn4MLxb1 , Spn4MLxb2 , Spn4MLxb3 , Spn4MLyb1 , &
                                   Spn4MLyb2 , Spn4MLyb3 , Spn4MLzb1 , Spn4MLzb2 , Spn4MLzb3 , Spn4RDxb1 , Spn4RDxb2 , Spn4RDxb3 , &
                                   Spn4RDyb1 , Spn4RDyb2 , Spn4RDyb3 , Spn4RDzb1 , Spn4RDzb2 , Spn4RDzb3 , Spn4TDxb1 , Spn4TDxb2 , &
                                   Spn4TDxb3 , Spn4TDyb1 , Spn4TDyb2 , Spn4TDyb3 , Spn4TDzb1 , Spn4TDzb2 , Spn4TDzb3 , Spn5ALxb1 , &
                                   Spn5ALxb2 , Spn5ALxb3 , Spn5ALyb1 , Spn5ALyb2 , Spn5ALyb3 , Spn5ALzb1 , Spn5ALzb2 , Spn5ALzb3 , &
                                   Spn5FLxb1 , Spn5FLxb2 , Spn5FLxb3 , Spn5FLyb1 , Spn5FLyb2 , Spn5FLyb3 , Spn5FLzb1 , Spn5FLzb2 , &
                                   Spn5FLzb3 , Spn5MLxb1 , Spn5MLxb2 , Spn5MLxb3 , Spn5MLyb1 , Spn5MLyb2 , Spn5MLyb3 , Spn5MLzb1 , &
                                   Spn5MLzb2 , Spn5MLzb3 , Spn5RDxb1 , Spn5RDxb2 , Spn5RDxb3 , Spn5RDyb1 , Spn5RDyb2 , Spn5RDyb3 , &
                                   Spn5RDzb1 , Spn5RDzb2 , Spn5RDzb3 , Spn5TDxb1 , Spn5TDxb2 , Spn5TDxb3 , Spn5TDyb1 , Spn5TDyb2 , &
                                   Spn5TDyb3 , Spn5TDzb1 , Spn5TDzb2 , Spn5TDzb3 , Spn6ALxb1 , Spn6ALxb2 , Spn6ALxb3 , Spn6ALyb1 , &
                                   Spn6ALyb2 , Spn6ALyb3 , Spn6ALzb1 , Spn6ALzb2 , Spn6ALzb3 , Spn6FLxb1 , Spn6FLxb2 , Spn6FLxb3 , &
                                   Spn6FLyb1 , Spn6FLyb2 , Spn6FLyb3 , Spn6FLzb1 , Spn6FLzb2 , Spn6FLzb3 , Spn6MLxb1 , Spn6MLxb2 , &
                                   Spn6MLxb3 , Spn6MLyb1 , Spn6MLyb2 , Spn6MLyb3 , Spn6MLzb1 , Spn6MLzb2 , Spn6MLzb3 , Spn6RDxb1 , &
                                   Spn6RDxb2 , Spn6RDxb3 , Spn6RDyb1 , Spn6RDyb2 , Spn6RDyb3 , Spn6RDzb1 , Spn6RDzb2 , Spn6RDzb3 , &
                                   Spn6TDxb1 , Spn6TDxb2 , Spn6TDxb3 , Spn6TDyb1 , Spn6TDyb2 , Spn6TDyb3 , Spn6TDzb1 , Spn6TDzb2 , &
                                   Spn6TDzb3 , Spn7ALxb1 , Spn7ALxb2 , Spn7ALxb3 , Spn7ALyb1 , Spn7ALyb2 , Spn7ALyb3 , Spn7ALzb1 , &
                                   Spn7ALzb2 , Spn7ALzb3 , Spn7FLxb1 , Spn7FLxb2 , Spn7FLxb3 , Spn7FLyb1 , Spn7FLyb2 , Spn7FLyb3 , &
                                   Spn7FLzb1 , Spn7FLzb2 , Spn7FLzb3 , Spn7MLxb1 , Spn7MLxb2 , Spn7MLxb3 , Spn7MLyb1 , Spn7MLyb2 , &
                                   Spn7MLyb3 , Spn7MLzb1 , Spn7MLzb2 , Spn7MLzb3 , Spn7RDxb1 , Spn7RDxb2 , Spn7RDxb3 , Spn7RDyb1 , &
                                   Spn7RDyb2 , Spn7RDyb3 , Spn7RDzb1 , Spn7RDzb2 , Spn7RDzb3 , Spn7TDxb1 , Spn7TDxb2 , Spn7TDxb3 , &
                                   Spn7TDyb1 , Spn7TDyb2 , Spn7TDyb3 , Spn7TDzb1 , Spn7TDzb2 , Spn7TDzb3 , Spn8ALxb1 , Spn8ALxb2 , &
                                   Spn8ALxb3 , Spn8ALyb1 , Spn8ALyb2 , Spn8ALyb3 , Spn8ALzb1 , Spn8ALzb2 , Spn8ALzb3 , Spn8FLxb1 , &
                                   Spn8FLxb2 , Spn8FLxb3 , Spn8FLyb1 , Spn8FLyb2 , Spn8FLyb3 , Spn8FLzb1 , Spn8FLzb2 , Spn8FLzb3 , &
                                   Spn8MLxb1 , Spn8MLxb2 , Spn8MLxb3 , Spn8MLyb1 , Spn8MLyb2 , Spn8MLyb3 , Spn8MLzb1 , Spn8MLzb2 , &
                                   Spn8MLzb3 , Spn8RDxb1 , Spn8RDxb2 , Spn8RDxb3 , Spn8RDyb1 , Spn8RDyb2 , Spn8RDyb3 , Spn8RDzb1 , &
                                   Spn8RDzb2 , Spn8RDzb3 , Spn8TDxb1 , Spn8TDxb2 , Spn8TDxb3 , Spn8TDyb1 , Spn8TDyb2 , Spn8TDyb3 , &
                                   Spn8TDzb1 , Spn8TDzb2 , Spn8TDzb3 , Spn9ALxb1 , Spn9ALxb2 , Spn9ALxb3 , Spn9ALyb1 , Spn9ALyb2 , &
                                   Spn9ALyb3 , Spn9ALzb1 , Spn9ALzb2 , Spn9ALzb3 , Spn9FLxb1 , Spn9FLxb2 , Spn9FLxb3 , Spn9FLyb1 , &
                                   Spn9FLyb2 , Spn9FLyb3 , Spn9FLzb1 , Spn9FLzb2 , Spn9FLzb3 , Spn9MLxb1 , Spn9MLxb2 , Spn9MLxb3 , &
                                   Spn9MLyb1 , Spn9MLyb2 , Spn9MLyb3 , Spn9MLzb1 , Spn9MLzb2 , Spn9MLzb3 , Spn9RDxb1 , Spn9RDxb2 , &
                                   Spn9RDxb3 , Spn9RDyb1 , Spn9RDyb2 , Spn9RDyb3 , Spn9RDzb1 , Spn9RDzb2 , Spn9RDzb3 , Spn9TDxb1 , &
                                   Spn9TDxb2 , Spn9TDxb3 , Spn9TDyb1 , Spn9TDyb2 , Spn9TDyb3 , Spn9TDzb1 , Spn9TDzb2 , Spn9TDzb3 , &
                                   TailFurlP , TailFurlA , TailFurlP , TailFurlV ,   TeetAya ,   TeetPya ,   TeetPya ,   TeetVya , &
                                     TFrlBrM , TipClrnc1 , TipClrnc2 , TipClrnc3 ,  TipALxb1 ,  TipALxb2 ,  TipALxb3 ,  TipALyb1 , &
                                    TipALyb2 ,  TipALyb3 ,  TipALzb1 ,  TipALzb2 ,  TipALzb3 , TipClrnc1 , TipClrnc2 , TipClrnc3 , &
                                     TipDxb1 ,   TipDxb2 ,   TipDxb3 ,   TipDxc1 ,   TipDxc2 ,   TipDxc3 ,   TipDyb1 ,   TipDyb2 , &
                                     TipDyb3 ,   TipDyc1 ,   TipDyc2 ,   TipDyc3 ,   TipDzc1 ,   TipDzc2 ,   TipDzc3 ,   TipDzc1 , &
                                     TipDzc2 ,   TipDzc3 ,  TipRDxb1 ,  TipRDxb2 ,  TipRDxb3 ,  TipRDyb1 ,  TipRDyb2 ,  TipRDyb3 , &
                                    TipRDzc1 ,  TipRDzc2 ,  TipRDzc3 ,  TipRDzc1 ,  TipRDzc2 ,  TipRDzc3 , YawBrTDzt , YawBrTDxt , &
                                   YawBrRDyt , YawBrRDxt , YawBrTDyt , YawBrRDzt , TwHt1ALxt , TwHt1ALyt , TwHt1ALzt , TwHt1FLxt , &
                                   TwHt1FLyt , TwHt1FLzt , TwHt1MLxt , TwHt1MLyt , TwHt1MLzt , TwHt1RDxt , TwHt1RDyt , TwHt1RDzt , &
                                   TwHt1RPxi , TwHt1RPyi , TwHt1RPzi , TwHt1TDxt , TwHt1TDyt , TwHt1TDzt , TwHt1TPxi , TwHt1TPyi , &
                                   TwHt1TPzi , TwHt2ALxt , TwHt2ALyt , TwHt2ALzt , TwHt2FLxt , TwHt2FLyt , TwHt2FLzt , TwHt2MLxt , &
                                   TwHt2MLyt , TwHt2MLzt , TwHt2RDxt , TwHt2RDyt , TwHt2RDzt , TwHt2RPxi , TwHt2RPyi , TwHt2RPzi , &
                                   TwHt2TDxt , TwHt2TDyt , TwHt2TDzt , TwHt2TPxi , TwHt2TPyi , TwHt2TPzi , TwHt3ALxt , TwHt3ALyt , &
                                   TwHt3ALzt , TwHt3FLxt , TwHt3FLyt , TwHt3FLzt , TwHt3MLxt , TwHt3MLyt , TwHt3MLzt , TwHt3RDxt , &
                                   TwHt3RDyt , TwHt3RDzt , TwHt3RPxi , TwHt3RPyi , TwHt3RPzi , TwHt3TDxt , TwHt3TDyt , TwHt3TDzt , &
                                   TwHt3TPxi , TwHt3TPyi , TwHt3TPzi , TwHt4ALxt , TwHt4ALyt , TwHt4ALzt , TwHt4FLxt , TwHt4FLyt , &
                                   TwHt4FLzt , TwHt4MLxt , TwHt4MLyt , TwHt4MLzt , TwHt4RDxt , TwHt4RDyt , TwHt4RDzt , TwHt4RPxi , &
                                   TwHt4RPyi , TwHt4RPzi , TwHt4TDxt , TwHt4TDyt , TwHt4TDzt , TwHt4TPxi , TwHt4TPyi , TwHt4TPzi , &
                                   TwHt5ALxt , TwHt5ALyt , TwHt5ALzt , TwHt5FLxt , TwHt5FLyt , TwHt5FLzt , TwHt5MLxt , TwHt5MLyt , &
                                   TwHt5MLzt , TwHt5RDxt , TwHt5RDyt , TwHt5RDzt , TwHt5RPxi , TwHt5RPyi , TwHt5RPzi , TwHt5TDxt , &
                                   TwHt5TDyt , TwHt5TDzt , TwHt5TPxi , TwHt5TPyi , TwHt5TPzi , TwHt6ALxt , TwHt6ALyt , TwHt6ALzt , &
                                   TwHt6FLxt , TwHt6FLyt , TwHt6FLzt , TwHt6MLxt , TwHt6MLyt , TwHt6MLzt , TwHt6RDxt , TwHt6RDyt , &
                                   TwHt6RDzt , TwHt6RPxi , TwHt6RPyi , TwHt6RPzi , TwHt6TDxt , TwHt6TDyt , TwHt6TDzt , TwHt6TPxi , &
                                   TwHt6TPyi , TwHt6TPzi , TwHt7ALxt , TwHt7ALyt , TwHt7ALzt , TwHt7FLxt , TwHt7FLyt , TwHt7FLzt , &
                                   TwHt7MLxt , TwHt7MLyt , TwHt7MLzt , TwHt7RDxt , TwHt7RDyt , TwHt7RDzt , TwHt7RPxi , TwHt7RPyi , &
                                   TwHt7RPzi , TwHt7TDxt , TwHt7TDyt , TwHt7TDzt , TwHt7TPxi , TwHt7TPyi , TwHt7TPzi , TwHt8ALxt , &
                                   TwHt8ALyt , TwHt8ALzt , TwHt8FLxt , TwHt8FLyt , TwHt8FLzt , TwHt8MLxt , TwHt8MLyt , TwHt8MLzt , &
                                   TwHt8RDxt , TwHt8RDyt , TwHt8RDzt , TwHt8RPxi , TwHt8RPyi , TwHt8RPzi , TwHt8TDxt , TwHt8TDyt , &
                                   TwHt8TDzt , TwHt8TPxi , TwHt8TPyi , TwHt8TPzi , TwHt9ALxt , TwHt9ALyt , TwHt9ALzt , TwHt9FLxt , &
                                   TwHt9FLyt , TwHt9FLzt , TwHt9MLxt , TwHt9MLyt , TwHt9MLzt , TwHt9RDxt , TwHt9RDyt , TwHt9RDzt , &
                                   TwHt9RPxi , TwHt9RPyi , TwHt9RPzi , TwHt9TDxt , TwHt9TDyt , TwHt9TDzt , TwHt9TPxi , TwHt9TPyi , &
                                   TwHt9TPzi ,  TwrBsFxt ,  TwrBsFyt ,  TwrBsFzt ,  TwrBsMxt ,  TwrBsMyt ,  TwrBsMzt , TipClrnc1 , &
                                   TipClrnc2 , TipClrnc3 , TwrTpTDxi , TwrTpTDyi , TwrTpTDzi ,  TipRDzc1 ,  TipRDzc2 ,  TipRDzc3 , &
                                      YawAzn ,    YawAzn ,    YawAzn ,  YawBrFxn ,  YawBrFxp ,  YawBrFyn ,  YawBrFyp ,  YawBrFzn , &
                                    YawBrFzn ,  YawBrMxn ,  YawBrMxp ,  YawBrMyn ,  YawBrMyp ,  YawBrMzn ,  YawBrMzn , YawBrRAxp , &
                                   YawBrRAyp , YawBrRAzp , YawBrRDxt , YawBrRDyt , YawBrRDzt , YawBrRVxp , YawBrRVyp , YawBrRVzp , &
                                   YawBrTAxp , YawBrTAyp , YawBrTAzp , TwrTpTDxi , YawBrTDxp , YawBrTDxt , TwrTpTDyi , YawBrTDyp , &
                                   YawBrTDyt , TwrTpTDzi , YawBrTDzp , YawBrTDzt , YawBrTVxp , YawBrTVyp , YawBrTVzp ,    YawPzn , &
                                   YawPzn ,    YawPzn ,    YawVzn ,    YawVzn ,    YawVzn /)
      CHARACTER(ChanLen), PARAMETER :: ParamUnitsAry(981) =  (/  &  ! This lists the units corresponding to the allowed parameters
                                  "(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg/s^2)", &
                                  "(rpm)    ","(kN-m)   ","(deg/s^2)","(kW)     ","(kN-m)   ","(rpm)    ","(m)      ","(m)      ", &
                                  "(m)      ","(deg/s^2)","(deg/s^2)","(deg/s^2)","(kN)     ","(kN)     ","(kN)     ","(kN)     ", &
                                  "(kN)     ","(kN)     ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ", &
                                  "(deg)    ","(deg)    ","(deg)    ","(rpm)    ","(rpm)    ","(rpm)    ","(kN)     ","(kN)     ", &
                                  "(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN-m)   ","(kN-m)   ","(kW)     ","(kN-m)   ", &
                                  "(deg/s^2)","(deg/s^2)","(deg/s^2)","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(deg)    ", &
                                  "(deg)    ","(deg)    ","(rpm)    ","(rpm)    ","(rpm)    ","(deg)    ","(deg/s^2)","(deg)    ", &
                                  "(deg/s)  ","(deg/s^2)","(deg/s^2)","(deg/s^2)","(deg/s)  ","(deg/s)  ","(deg/s)  ","(m/s^2)  ", &
                                  "(m/s^2)  ","(m/s^2)  ","(m/s)    ","(m/s)    ","(m/s)    ","(m)      ","(m)      ","(m)      ", &
                                  "(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ", &
                                  "(deg)    ","(m)      ","(deg)    ","(deg/s^2)","(deg/s^2)","(deg/s^2)","(deg/s^2)","(deg/s^2)", &
                                  "(deg/s^2)","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg/s)  ","(deg/s)  ","(deg/s)  ", &
                                  "(deg/s)  ","(deg/s)  ","(deg/s)  ","(m)      ","(m)      ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ", &
                                  "(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ", &
                                  "(m)      ","(m/s)    ","(m/s)    ","(m/s)    ","(m/s)    ","(m/s)    ","(m/s)    ","(deg)    ", &
                                  "(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ", &
                                  "(m/s^2)  ","(rad/s^2)","(rad/s^2)","(m/s^2)  ","(rad/s^2)","(rad/s^2)","(rad/s^2)","(m/s^2)  ", &
                                  "(m/s^2)  ","(rad/s^2)","(m/s^2)  ","(m/s^2)  ","(rad/s^2)","(m/s^2)  ","(m/s^2)  ","(rad/s^2)", &
                                  "(rad/s^2)","(m/s)    ","(m/s)    ","(m/s)    ","(m/s)    ","(m/s)    ","(m/s)    ","(m/s)    ", &
                                  "(m/s)    ","(m/s)    ","(rad/s)  ","(rad/s)  ","(m/s)    ","(rad/s)  ","(rad/s)  ","(rad/s)  ", &
                                  "(m/s)    ","(m/s)    ","(rad/s)  ","(m/s)    ","(m/s)    ","(rad/s)  ","(m/s)    ","(m/s)    ", &
                                  "(rad/s)  ","(rad/s)  ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ", &
                                  "(m)      ","(m)      ","(m)      ","(rad)    ","(rad)    ","(m)      ","(rad)    ","(rad)    ", &
                                  "(rad)    ","(m)      ","(m)      ","(rad)    ","(m)      ","(m)      ","(rad)    ","(m)      ", &
                                  "(m)      ","(rad)    ","(rad)    ","(kN-m)   ","(deg)    ","(deg)    ","(deg)    ","(kN)     ", &
                                  "(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ", &
                                  "(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ", &
                                  "(kN)     ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ", &
                                  "(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ", &
                                  "(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ", &
                                  "(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(deg/s^2)", &
                                  "(deg)    ","(deg/s^2)","(deg)    ","(deg/s)  ","(kW)     ","(rpm)    ","(deg/s^2)","(deg)    ", &
                                  "(deg/s)  ","(kN)     ","(kN-m)   ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ", &
                                  "(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(kN)     ","(kN)     ","(kN)     ","(kN)     ", &
                                  "(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN-m)   ","(kN-m)   ","(kN-m)   ", &
                                  "(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(deg)    ","(deg)    ", &
                                  "(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(m)      ", &
                                  "(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ", &
                                  "(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ", &
                                  "(m/s^2)  ","(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ", &
                                  "(kN)     ","(kN)     ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ", &
                                  "(kN-m)   ","(kN-m)   ","(kN-m)   ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ", &
                                  "(deg)    ","(deg)    ","(deg)    ","(deg)    ","(m)      ","(m)      ","(m)      ","(m)      ", &
                                  "(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ", &
                                  "(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(kN)     ","(kN)     ", &
                                  "(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN-m)   ", &
                                  "(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ", &
                                  "(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ", &
                                  "(deg)    ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ", &
                                  "(m)      ","(m)      ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ", &
                                  "(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ", &
                                  "(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ", &
                                  "(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(deg)    ","(deg)    ","(deg)    ", &
                                  "(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(m)      ","(m)      ", &
                                  "(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m/s^2)  ", &
                                  "(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ", &
                                  "(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ", &
                                  "(kN)     ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ", &
                                  "(kN-m)   ","(kN-m)   ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ", &
                                  "(deg)    ","(deg)    ","(deg)    ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ", &
                                  "(m)      ","(m)      ","(m)      ","(m)      ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ", &
                                  "(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(kN)     ","(kN)     ","(kN)     ", &
                                  "(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN-m)   ","(kN-m)   ", &
                                  "(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(deg)    ", &
                                  "(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ", &
                                  "(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ", &
                                  "(m)      ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ", &
                                  "(m/s^2)  ","(m/s^2)  ","(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ", &
                                  "(kN)     ","(kN)     ","(kN)     ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ", &
                                  "(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(deg)    ","(deg)    ","(deg)    ","(deg)    ", &
                                  "(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(m)      ","(m)      ","(m)      ", &
                                  "(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m/s^2)  ","(m/s^2)  ", &
                                  "(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(kN)     ", &
                                  "(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ", &
                                  "(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ", &
                                  "(kN-m)   ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ", &
                                  "(deg)    ","(deg)    ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ", &
                                  "(m)      ","(m)      ","(m)      ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ", &
                                  "(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(kN)     ","(kN)     ","(kN)     ","(kN)     ", &
                                  "(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN-m)   ","(kN-m)   ","(kN-m)   ", &
                                  "(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(deg)    ","(deg)    ", &
                                  "(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(m)      ", &
                                  "(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ", &
                                  "(deg)    ","(deg/s^2)","(deg)    ","(deg/s)  ","(deg/s^2)","(deg)    ","(deg)    ","(deg/s)  ", &
                                  "(kN-m)   ","(m)      ","(m)      ","(m)      ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ", &
                                  "(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m)      ","(m)      ","(m)      ", &
                                  "(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ", &
                                  "(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ", &
                                  "(m)      ","(m)      ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ", &
                                  "(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(m)      ","(m)      ", &
                                  "(deg)    ","(deg)    ","(m)      ","(deg)    ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(kN)     ", &
                                  "(kN)     ","(kN)     ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(deg)    ","(deg)    ","(deg)    ", &
                                  "(deg)    ","(deg)    ","(deg)    ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ", &
                                  "(m)      ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(kN)     ","(kN)     ","(kN)     ","(kN-m)   ", &
                                  "(kN-m)   ","(kN-m)   ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ", &
                                  "(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m/s^2)  ","(m/s^2)  ", &
                                  "(m/s^2)  ","(kN)     ","(kN)     ","(kN)     ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(deg)    ", &
                                  "(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(m)      ","(m)      ","(m)      ", &
                                  "(m)      ","(m)      ","(m)      ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(kN)     ","(kN)     ", &
                                  "(kN)     ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(deg)    ","(deg)    ","(deg)    ","(deg)    ", &
                                  "(deg)    ","(deg)    ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ", &
                                  "(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(kN)     ","(kN)     ","(kN)     ","(kN-m)   ","(kN-m)   ", &
                                  "(kN-m)   ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(m)      ", &
                                  "(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ", &
                                  "(kN)     ","(kN)     ","(kN)     ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(deg)    ","(deg)    ", &
                                  "(deg)    ","(deg)    ","(deg)    ","(deg)    ","(m)      ","(m)      ","(m)      ","(m)      ", &
                                  "(m)      ","(m)      ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(kN)     ","(kN)     ","(kN)     ", &
                                  "(kN-m)   ","(kN-m)   ","(kN-m)   ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ", &
                                  "(deg)    ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(m/s^2)  ", &
                                  "(m/s^2)  ","(m/s^2)  ","(kN)     ","(kN)     ","(kN)     ","(kN-m)   ","(kN-m)   ","(kN-m)   ", &
                                  "(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(m)      ","(m)      ", &
                                  "(m)      ","(m)      ","(m)      ","(m)      ","(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(kN)     ", &
                                  "(kN)     ","(kN)     ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(deg)    ","(deg)    ","(deg)    ", &
                                  "(deg)    ","(deg)    ","(deg)    ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ", &
                                  "(m)      ","(kN)     ","(kN)     ","(kN)     ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(m)      ", &
                                  "(m)      ","(m)      ","(m)      ","(m)      ","(m)      ","(deg)    ","(deg)    ","(deg)    ", &
                                  "(deg/s^2)","(deg/s^2)","(deg/s^2)","(kN)     ","(kN)     ","(kN)     ","(kN)     ","(kN)     ", &
                                  "(kN)     ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(kN-m)   ","(deg/s^2)", &
                                  "(deg/s^2)","(deg/s^2)","(deg)    ","(deg)    ","(deg)    ","(deg/s)  ","(deg/s)  ","(deg/s)  ", &
                                  "(m/s^2)  ","(m/s^2)  ","(m/s^2)  ","(m)      ","(m)      ","(m)      ","(m)      ","(m)      ", &
                                  "(m)      ","(m)      ","(m)      ","(m)      ","(m/s)    ","(m/s)    ","(m/s)    ","(deg)    ", &
                                  "(deg)    ","(deg)    ","(deg/s)  ","(deg/s)  ","(deg/s)  "/)


      ! Initialize values
   ErrStat = ErrID_None
   ErrMsg = ""
   InvalidOutput = .FALSE.


!   ..... Developer must add checking for invalid inputs here: .....
if (p%BD4Blades) then
   startIndx = 1
else
   startIndx = p%NumBl+1
end if

   DO I = startIndx,3  ! Invalid blades

         ! motions

      InvalidOutput(   TipDxc(  I) ) = .TRUE.
      InvalidOutput(   TipDyc(  I) ) = .TRUE.
      InvalidOutput(   TipDzc(  I) ) = .TRUE.
      InvalidOutput(   TipDxb(  I) ) = .TRUE.
      InvalidOutput(   TipDyb(  I) ) = .TRUE.
      InvalidOutput(  TipALxb(  I) ) = .TRUE.
      InvalidOutput(  TipALyb(  I) ) = .TRUE.
      InvalidOutput(  TipALzb(  I) ) = .TRUE.
      InvalidOutput(  TipRDxb(  I) ) = .TRUE.
      InvalidOutput(  TipRDyb(  I) ) = .TRUE.
      InvalidOutput(  TipRDzc(  I) ) = .TRUE.
      InvalidOutput( TipClrnc(  I) ) = .TRUE.

         ! loads

      InvalidOutput(  RootFxc(  I) ) = .TRUE.
      InvalidOutput(  RootFyc(  I) ) = .TRUE.
      InvalidOutput(  RootFzc(  I) ) = .TRUE.
      InvalidOutput(  RootFxb(  I) ) = .TRUE.
      InvalidOutput(  RootFyb(  I) ) = .TRUE.
      InvalidOutput(  RootMxc(  I) ) = .TRUE.
      InvalidOutput(  RootMyc(  I) ) = .TRUE.
      InvalidOutput(  RootMzc(  I) ) = .TRUE.
      InvalidOutput(  RootMxb(  I) ) = .TRUE.
      InvalidOutput(  RootMyb(  I) ) = .TRUE.

         ! Blade node motions

      InvalidOutput(  SpnALxb(:,I) ) = .TRUE.
      InvalidOutput(  SpnALyb(:,I) ) = .TRUE.
      InvalidOutput(  SpnALzb(:,I) ) = .TRUE.

      InvalidOutput(  SpnTDxb(:,I) ) = .TRUE.
      InvalidOutput(  SpnTDyb(:,I) ) = .TRUE.
      InvalidOutput(  SpnTDzb(:,I) ) = .TRUE.

      InvalidOutput(  SpnRDxb(:,I) ) = .TRUE.
      InvalidOutput(  SpnRDyb(:,I) ) = .TRUE.
      InvalidOutput(  SpnRDzb(:,I) ) = .TRUE.

         ! Blade node loads

      InvalidOutput(  SpnMLxb(:,I) ) = .TRUE.
      InvalidOutput(  SpnMLyb(:,I) ) = .TRUE.
      InvalidOutput(  SpnMLzb(:,I) ) = .TRUE.

      InvalidOutput(  SpnFLxb(:,I) ) = .TRUE.
      InvalidOutput(  SpnFLyb(:,I) ) = .TRUE.
      InvalidOutput(  SpnFLzb(:,I) ) = .TRUE.

   END DO


   DO I = 1,p%NumBl

      DO J = p%NBlGages+1,9 ! Invalid blade gages

         InvalidOutput(  SpnALxb(J,I) ) = .TRUE.
         InvalidOutput(  SpnALyb(J,I) ) = .TRUE.
         InvalidOutput(  SpnALzb(J,I) ) = .TRUE.

         InvalidOutput(  SpnTDxb(J,I) ) = .TRUE.
         InvalidOutput(  SpnTDyb(J,I) ) = .TRUE.
         InvalidOutput(  SpnTDzb(J,I) ) = .TRUE.

         InvalidOutput(  SpnRDxb(J,I) ) = .TRUE.
         InvalidOutput(  SpnRDyb(J,I) ) = .TRUE.
         InvalidOutput(  SpnRDzb(J,I) ) = .TRUE.

            ! Loads

         InvalidOutput(  SpnMLxb(J,I) ) = .TRUE.
         InvalidOutput(  SpnMLyb(J,I) ) = .TRUE.
         InvalidOutput(  SpnMLzb(J,I) ) = .TRUE.

         InvalidOutput(  SpnFLxb(J,I) ) = .TRUE.
         InvalidOutput(  SpnFLyb(J,I) ) = .TRUE.
         InvalidOutput(  SpnFLzb(J,I) ) = .TRUE.


      END DO !J

   END DO !I

   DO J = p%NTwGages+1,9 !Invalid tower gages

         ! Motions

      InvalidOutput( TwHtALxt(J) ) = .TRUE.
      InvalidOutput( TwHtALyt(J) ) = .TRUE.
      InvalidOutput( TwHtALzt(J) ) = .TRUE.

      InvalidOutput( TwHtTDxt(J) ) = .TRUE.
      InvalidOutput( TwHtTDyt(J) ) = .TRUE.
      InvalidOutput( TwHtTDzt(J) ) = .TRUE.

      InvalidOutput( TwHtRDxt(J) ) = .TRUE.
      InvalidOutput( TwHtRDyt(J) ) = .TRUE.
      InvalidOutput( TwHtRDzt(J) ) = .TRUE.

      InvalidOutput( TwHtTPxi(J) ) = .TRUE.
      InvalidOutput( TwHtTPyi(J) ) = .TRUE.
      InvalidOutput( TwHtTPzi(J) ) = .TRUE.

      InvalidOutput( TwHtRPxi(J) ) = .TRUE.
      InvalidOutput( TwHtRPyi(J) ) = .TRUE.
      InvalidOutput( TwHtRPzi(J) ) = .TRUE.

         ! Loads

      InvalidOutput( TwHtMLxt(J) ) = .TRUE.
      InvalidOutput( TwHtMLyt(J) ) = .TRUE.
      InvalidOutput( TwHtMLzt(J) ) = .TRUE.

      InvalidOutput( TwHtFLxt(J) ) = .TRUE.
      InvalidOutput( TwHtFLyt(J) ) = .TRUE.
      InvalidOutput( TwHtFLzt(J) ) = .TRUE.

   END DO


   ! Invalid outputs based on number of blades
   IF ( p%NumBl < 3 ) THEN
         InvalidOutput(PtchPMzc3) = .TRUE.
         InvalidOutput(   Q_B3E1) = .TRUE.
         InvalidOutput(   Q_B3F1) = .TRUE.
         InvalidOutput(   Q_B3F2) = .TRUE.
         InvalidOutput(  QD_B3E1) = .TRUE.
         InvalidOutput(  QD_B3F1) = .TRUE.
         InvalidOutput(  QD_B3F2) = .TRUE.
         InvalidOutput( QD2_B3E1) = .TRUE.
         InvalidOutput( QD2_B3F1) = .TRUE.
         InvalidOutput( QD2_B3F2) = .TRUE.      
   ENDIF
   IF ( p%NumBl < 2 ) THEN
         InvalidOutput(PtchPMzc2) = .TRUE.
         InvalidOutput(   Q_B2E1) = .TRUE.
         InvalidOutput(   Q_B2F1) = .TRUE.
         InvalidOutput(   Q_B2F2) = .TRUE.
         InvalidOutput(  QD_B2E1) = .TRUE.
         InvalidOutput(  QD_B2F1) = .TRUE.
         InvalidOutput(  QD_B2F2) = .TRUE.
         InvalidOutput( QD2_B2E1) = .TRUE.
         InvalidOutput( QD2_B2F1) = .TRUE.
         InvalidOutput( QD2_B2F2) = .TRUE.      
   ENDIF
   ! 1-bladed or 3-bladed, no teeter
   IF ( p%NumBl /= 2 ) THEN
      InvalidOutput(  TeetPya) = .TRUE.
      InvalidOutput(  TeetVya) = .TRUE.
      InvalidOutput(  TeetAya) = .TRUE.

      InvalidOutput(   Q_Teet) = .TRUE.
      InvalidOutput(  QD_Teet) = .TRUE.
      InvalidOutput( QD2_Teet) = .TRUE.
   END IF
   
   InvalidOutput(HSSBrTq) = p%method == Method_RK4

    IF ( p%BD4Blades ) THEN
      InvalidOutput(   Q_B1E1) = .TRUE.
      InvalidOutput(   Q_B1F1) = .TRUE.
      InvalidOutput(   Q_B1F2) = .TRUE.

      InvalidOutput(  QD_B1E1) = .TRUE.
      InvalidOutput(  QD_B1F1) = .TRUE.
      InvalidOutput(  QD_B1F2) = .TRUE.

      InvalidOutput( QD2_B1E1) = .TRUE.
      InvalidOutput( QD2_B1F1) = .TRUE.
      InvalidOutput( QD2_B1F2) = .TRUE.      
      
      InvalidOutput(   Q_B2E1) = .TRUE.
      InvalidOutput(   Q_B2F1) = .TRUE.
      InvalidOutput(   Q_B2F2) = .TRUE.

      InvalidOutput(  QD_B2E1) = .TRUE.
      InvalidOutput(  QD_B2F1) = .TRUE.
      InvalidOutput(  QD_B2F2) = .TRUE.

      InvalidOutput( QD2_B2E1) = .TRUE.
      InvalidOutput( QD2_B2F1) = .TRUE.
      InvalidOutput( QD2_B2F2) = .TRUE.
      
      InvalidOutput(   Q_B3E1) = .TRUE.
      InvalidOutput(   Q_B3F1) = .TRUE.
      InvalidOutput(   Q_B3F2) = .TRUE.

      InvalidOutput(  QD_B3E1) = .TRUE.
      InvalidOutput(  QD_B3F1) = .TRUE.
      InvalidOutput(  QD_B3F2) = .TRUE.

      InvalidOutput( QD2_B3E1) = .TRUE.
      InvalidOutput( QD2_B3F1) = .TRUE.
      InvalidOutput( QD2_B3F2) = .TRUE.      
   END IF
!   ................. End of validity checking .................

   !-------------------------------------------------------------------------------------------------
   ! Allocate and set index, name, and units for the output channels
   ! If a selected output channel is not available in this module, set error flag.
   !-------------------------------------------------------------------------------------------------

   ALLOCATE ( p%OutParam(0:p%NumOuts) , STAT=ErrStat2 )
   IF ( ErrStat2 /= 0_IntKi )  THEN
      CALL SetErrStat( ErrID_Fatal,"Error allocating memory for the ElastoDyn OutParam array.", ErrStat, ErrMsg, RoutineName )
      RETURN
   ENDIF

      ! Set index, name, and units for the time output channel:

   p%OutParam(0)%Indx  = Time
   p%OutParam(0)%Name  = "Time"    ! OutParam(0) is the time channel by default.
   p%OutParam(0)%Units = "(s)"
   p%OutParam(0)%SignM = 1


      ! Set index, name, and units for all of the output channels.
      ! If a selected output channel is not available by this module set ErrStat = ErrID_Warn.

   DO I = 1,p%NumOuts

      p%OutParam(I)%Name  = OutList(I)
      OutListTmp          = OutList(I)

      ! Reverse the sign (+/-) of the output channel if the user prefixed the
      !   channel name with a "-", "_", "m", or "M" character indicating "minus".


      CheckOutListAgain = .FALSE.

      IF      ( INDEX( "-_", OutListTmp(1:1) ) > 0 ) THEN
         p%OutParam(I)%SignM = -1                         ! ex, "-TipDxc1" causes the sign of TipDxc1 to be switched.
         OutListTmp          = OutListTmp(2:)
      ELSE IF ( INDEX( "mM", OutListTmp(1:1) ) > 0 ) THEN ! We'll assume this is a variable name for now, (if not, we will check later if OutListTmp(2:) is also a variable name)
         CheckOutListAgain   = .TRUE.
         p%OutParam(I)%SignM = 1
      ELSE
         p%OutParam(I)%SignM = 1
      END IF

      CALL Conv2UC( OutListTmp )    ! Convert OutListTmp to upper case


      Indx = IndexCharAry( OutListTmp(1:OutStrLenM1), ValidParamAry )


         ! If it started with an "M" (CheckOutListAgain) we didn't find the value in our list (Indx < 1)

      IF ( CheckOutListAgain .AND. Indx < 1 ) THEN    ! Let's assume that "M" really meant "minus" and then test again
         p%OutParam(I)%SignM = -1                     ! ex, "MTipDxc1" causes the sign of TipDxc1 to be switched.
         OutListTmp          = OutListTmp(2:)

         Indx = IndexCharAry( OutListTmp(1:OutStrLenM1), ValidParamAry )
      END IF


      IF ( Indx > 0 ) THEN ! we found the channel name
         p%OutParam(I)%Indx     = ParamIndxAry(Indx)
         IF ( InvalidOutput( ParamIndxAry(Indx) ) ) THEN  ! but, it isn't valid for these settings
            p%OutParam(I)%Units = "INVALID"
            p%OutParam(I)%SignM = 0
         ELSE
            p%OutParam(I)%Units = ParamUnitsAry(Indx) ! it's a valid output
         END IF
      ELSE ! this channel isn't valid
         p%OutParam(I)%Indx  = Time                 ! pick any valid channel (I just picked "Time" here because it's universal)
         p%OutParam(I)%Units = "INVALID"
         p%OutParam(I)%SignM = 0                    ! multiply all results by zero

         CALL SetErrStat(ErrID_Fatal, TRIM(p%OutParam(I)%Name)//" is not an available output channel.",ErrStat,ErrMsg,RoutineName)
      END IF

   END DO

   RETURN
END SUBROUTINE SetOutParam
!----------------------------------------------------------------------------------------------------------------------------------
!End of code generated by Matlab script
!**********************************************************************************************************************************
!> This routine is used to compute rotor (blade and hub) properties:
!!   KBF(), KBE(), CBF(), CBE(), FreqBF(), FreqBE(), AxRedBld(),
!!   TwistedSF(), BldMass(), FirstMom(), SecondMom(), BldCG(),
!!   RotMass, RotIner, Hubg1Iner, Hubg2Iner, rSAerCenn1(), and
!!   rSAerCenn2(), BElmtMass()
!! tower properties:
!!   KTFA(), KTSS(), CTFA(), CTSS(), FreqTFA(), FreqTSS(),
!!   AxRedTFA(), AxRedTSS(), TwrFASF(), TwrSSSF(), TwrMass, and
!!   TwrTpMass, TElmtMass()
!! structure that furls with the rotor (not including rotor) properties:
!!   RrfaIner
!! tail boom properties:
!!   AtfaIner
!! and nacelle properties:
!!   Nacd2Iner
SUBROUTINE Coeff(p,InputFileData, ErrStat, ErrMsg)
!..................................................................................................................................

      ! Passed variables

   TYPE(ED_ParameterType),        INTENT(INOUT)    :: p                             !< Parameters of the structural dynamics module
   TYPE(ED_InputFile),            INTENT(IN)       :: InputFileData                 !< all the data in the ElastoDyn input file
   INTEGER(IntKi),                INTENT(OUT)      :: ErrStat                       !< Error status
   CHARACTER(*),                  INTENT(OUT)      :: ErrMsg                        !< Error message when ErrStat =/ ErrID_None


      ! Local variables.

   REAL(ReKi)                   :: AxRdBld   (3,3)                                 ! Temporary result holding the current addition to the p%AxRedBld() array.
   REAL(ReKi)                   :: AxRdBldOld(3,3)                                 ! Previous AxRdBld (i.e., AxRdBld from the previous node)
   REAL(ReKi)                   :: AxRdTFA   (2,2)                                 ! Temporary result holding the current addition to the AxRedTFA() array.
   REAL(ReKi)                   :: AxRdTFAOld(2,2)                                 ! Previous AxRdTFA (i.e., AxRdTFA from the previous node)
   REAL(ReKi)                   :: AxRdTSS   (2,2)                                 ! Temporary result holding the current addition to the AxRedTSS() array.
   REAL(ReKi)                   :: AxRdTSSOld(2,2)                                 ! Previous AxRdTSS (i.e., AxRdTSS from the previous node)
   REAL(ReKi)                   :: TmpDist                                         ! Temporary distance used in the calculation of the aero center locations.
   REAL(ReKi)                   :: TmpDistj1                                       ! Temporary distance used in the calculation of the aero center locations.
   REAL(ReKi)                   :: TmpDistj2                                       ! Temporary distance used in the calculation of the aero center locations.
   REAL(ReKi)                   :: ElmntStff                                       ! (Temporary) stiffness of an element.
   REAL(ReKi)                   :: ElStffFA                                        ! (Temporary) tower fore-aft stiffness of an element
   REAL(ReKi)                   :: ElStffSS                                        ! (Temporary) tower side-to-side  stiffness of an element
   REAL(ReKi)                   :: FMomAbvNd (p%NumBl,p%BldNodes)                  ! FMomAbvNd(K,J) = portion of the first moment of blade K about the rotor centerline (not root, like FirstMom(K)) associated with everything above node J (including tip brake masses).
   REAL(ReKi)                   :: KBECent   (p%NumBl,1,1)                         ! Centrifugal-term of generalized edgewise stiffness of the blades.
   REAL(ReKi)                   :: KBFCent   (p%NumBl,2,2)                         ! Centrifugal-term of generalized flapwise stiffness of the blades.
   REAL(ReKi)                   :: KTFAGrav  (2,2)                                 ! Gravitational-term of generalized fore-aft stiffness of the tower.
   REAL(ReKi)                   :: KTSSGrav  (2,2)                                 ! Gravitational-term of generalized side-to-side stiffness of the tower.
   REAL(ReKi)                   :: MBE       (p%NumBl,1,1)                         ! Generalized edgewise mass of the blades.
   REAL(ReKi)                   :: MBF       (p%NumBl,2,2)                         ! Generalized flapwise mass of the blades.
   REAL(ReKi)                   :: MTFA      (2,2)                                 ! Generalized fore-aft mass of the tower.
   REAL(ReKi)                   :: MTSS      (2,2)                                 ! Generalized side-to-side mass of the tower.
   REAL(ReKi)                   :: Shape                                           ! Temporary result holding a value from the SHP function
   REAL(ReKi)                   :: Shape1                                          ! Temporary result holding a value from the SHP function
   REAL(ReKi)                   :: Shape2                                          ! Temporary result holding a value from the SHP function
   REAL(ReKi)                   :: TMssAbvNd (p%TwrNodes)                          ! Portion of the tower mass associated with everything above node J (including tower-top effects)
   REAL(ReKi)                   :: TwstdSF   (2,3,0:1)                             ! Temperory result holding the current addition to the TwistedSF() array.
   REAL(ReKi)                   :: TwstdSFOld(2,3,0:1)                             ! Previous TwstdSF (i.e., TwstdSF from the previous node)

   INTEGER(IntKi)               :: I                                               ! Generic index.
   INTEGER(IntKi)               :: J                                               ! Loops through nodes / elements.
   INTEGER(IntKi)               :: K                                               ! Loops through blades.
   INTEGER(IntKi)               :: L                                               ! Generic index


   ErrStat = ErrID_None
   ErrMsg  = ''

   !...............................................................................................................................
   ! Calculate the distances from point S on a blade to the aerodynamic center in the j1 and j2 directions:
   !...............................................................................................................................

   DO K = 1,p%NumBl          ! Loop through the blades

      DO J = 1,p%BldNodes    ! Loop through the blade nodes / elements

         TmpDist           = ( 0.25 - p%PitchAxis(K,J) )*p%Chord(J)   ! Distance along the chordline from point S (25% chord) to the aerodynamic center of the blade element J--positive towards the trailing edge.
         TmpDistj1         = TmpDist*p%SAeroTwst(J)                   ! Distance along the j1-axis   from point S (25% chord) to the aerodynamic center of the blade element J
         TmpDistj2         = TmpDist*p%CAeroTwst(J)                   ! Distance along the j2-axis   from point S (25% chord) to the aerodynamic center of the blade element J
         p%rSAerCenn1(K,J) = TmpDistj1*p%CThetaS(K,J) - TmpDistj2*p%SThetaS(K,J)
         p%rSAerCenn2(K,J) = TmpDistj1*p%SThetaS(K,J) + TmpDistj2*p%CThetaS(K,J)

      ENDDO ! J - Blade nodes / elements

   ENDDO    ! K - Blades


   !...............................................................................................................................
   ! Calculate the structure that furls with the rotor inertia term:
   !...............................................................................................................................

   p%RrfaIner  = InputFileData%RFrlIner - p%RFrlMass*(      (p%rVDxn**2    )*( 1.0 - p%CRFrlSkw2*p%CRFrlTlt2 ) &
                                     +    (p%rVDzn**2    )*                    p%CRFrlTlt2   &
                                     +    (p%rVDyn**2    )*( 1.0 - p%SRFrlSkw2*p%CRFrlTlt2 ) &
                                     - 2.0*p%rVDxn*p%rVDzn*        p%CRFrlSkew*p%CSRFrlTlt   &
                                     - 2.0*p%rVDxn*p%rVDyn*        p%CSRFrlSkw*p%CRFrlTlt2   &
                                     - 2.0*p%rVDzn*p%rVDyn*        p%SRFrlSkew*p%CSRFrlTlt     )
   IF ( p%RrfaIner < 0.0 )   THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' RFrlIner must not be less than RFrlMass*( perpendicular distance between rotor-furl'// &
               ' axis and CM of the structure that furls with the rotor [not including rotor] )^2.'
      RETURN
   END IF

   !...............................................................................................................................
   ! Calculate the tail boom inertia term:
   !...............................................................................................................................

   p%AtfaIner  = p%TFrlIner - p%BoomMass*(   p%rWIxn*p%rWIxn*( 1.0 - p%CTFrlSkw2*p%CTFrlTlt2 ) &
                                       +     p%rWIzn*p%rWIzn*                    p%CTFrlTlt2   &
                                       +     p%rWIyn*p%rWIyn*( 1.0 - p%STFrlSkw2*p%CTFrlTlt2 ) &
                                       - 2.0*p%rWIxn*p%rWIzn*        p%CTFrlSkew*p%CSTFrlTlt   &
                                       - 2.0*p%rWIxn*p%rWIyn*        p%CSTFrlSkw*p%CTFrlTlt2   &
                                       - 2.0*p%rWIzn*p%rWIyn*        p%STFrlSkew*p%CSTFrlTlt     )
   IF ( p%AtfaIner < 0.0 ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' TFrlIner must not be less than BoomMass*( perpendicular distance between tail-furl'// &
                                        ' axis and tail boom CM )^2.'
      RETURN
   ENDIF

   !...............................................................................................................................
   ! Calculate the nacelle inertia terms:
   !...............................................................................................................................

   p%Nacd2Iner = InputFileData%NacYIner - p%NacMass*( p%NacCMxn**2 + p%NacCMyn**2 ) ! Nacelle inertia about the d2-axis
   IF ( p%Nacd2Iner < 0.0 ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' NacYIner must not be less than NacMass*( NacCMxn^2 + NacCMyn^2 ).'
      RETURN
   END IF

      ! Calculate hub inertia about its centerline passing through its c.g..
      !   This calculation assumes that the hub for a 2-blader is essentially
      !   a uniform cylinder whose centerline is transverse through the cylinder
      !   passing through its c.g..  That is, for a 2-blader, Hubg1Iner =
      !   Hubg2Iner is the inertia of the hub about both the g1- and g2- axes.  For
      !   3-bladers, Hubg1Iner is simply equal to HubIner and Hubg2Iner is zero.
      ! Also, Initialize RotMass and RotIner to associated hub properties:

   IF ( p%NumBl == 2 )  THEN ! 2-blader
      p%Hubg1Iner = ( InputFileData%HubIner - p%HubMass*( ( p%UndSling - p%HubCM )**2 ) )/( p%CosDel3**2 )
      p%Hubg2Iner = p%Hubg1Iner
      IF ( p%Hubg1Iner < 0.0 ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg = ' HubIner must not be less than HubMass*( UndSling - HubCM )^2 for 2-blader.'
         RETURN
      END IF
   ELSE                    ! 3-blader
      p%Hubg1Iner = InputFileData%HubIner
      p%Hubg2Iner = 0.0
   ENDIF

   p%RotMass   = p%HubMass
   p%RotIner   = p%Hubg1Iner


   !...............................................................................................................................

      ! Initialize several variables to 0.0:

   p%KBF     = 0.0
   p%KBE     = 0.0
   KBFCent   = 0.0
   KBECent   = 0.0

   p%TwrMass = 0.0
   p%KTFA    = 0.0
   p%KTSS    = 0.0
   KTFAGrav  = 0.0
   KTSSGrav  = 0.0



   DO K = 1,p%NumBl          ! Loop through the blades


      ! Initialize BldMass(), FirstMom(), and SecondMom() using TipMass() effects:

      p%BldMass  (K) = p%TipMass(K)
      p%FirstMom (K) = p%TipMass(K)*p%BldFlexL
      p%SecondMom(K) = p%TipMass(K)*p%BldFlexL*p%BldFlexL


      DO J = p%BldNodes,1,-1 ! Loop through the blade nodes / elements in reverse


      ! Calculate the mass of the current element

         p%BElmntMass(J,K) = p%MassB(K,J)*p%DRNodes(J)                        ! Mass of blade element J


      ! Integrate to find some blade properties which will be output in .fsm

         p%BldMass  (K) = p%BldMass  (K) + p%BElmntMass(J,K)
         p%FirstMom (K) = p%FirstMom (K) + p%BElmntMass(J,K)*p%RNodes(J)
         p%SecondMom(K) = p%SecondMom(K) + p%BElmntMass(J,K)*p%RNodes(J)*p%RNodes(J)


      ! Integrate to find FMomAbvNd:

         FMomAbvNd   (K,J) = ( 0.5*p%BElmntMass(J,K) )*( p%HubRad + p%RNodes(J  ) + 0.5*p%DRNodes(J  ) )

         IF ( J == p%BldNodes )  THEN ! Outermost blade element
      ! Add the TipMass() effects:

            FMomAbvNd(K,J) = FmomAbvNd(K,J) + p%TipMass(K)*p%TipRad
         ELSE                       ! All other blade elements
      ! Add to FMomAbvNd(K,J) the effects from the (not yet used) portion of element J+1

            FMomAbvNd(K,J) = FMomAbvNd(K,J) + FMomAbvNd(K,J+1) &
                           + ( 0.5*p%BElmntMass(J+1,K) )*( p%HubRad + p%RNodes(J+1) - 0.5*p%DRNodes(J+1) )
         ENDIF


      ENDDO ! J - Blade nodes / elements in reverse

      IF (.NOT. p%BD4Blades) THEN
         ! Calculate BldCG() using FirstMom() and BldMass(); and calculate
         !   RotMass and RotIner:

         p%BldCG    (K) = p%FirstMom (K) / p%BldMass    (K)
         p%RotMass      = p%RotMass      + p%BldMass    (K)
         p%RotIner      = p%RotIner      + ( p%SecondMom(K) + p%BldMass  (K)*p%HubRad*( 2.0*p%BldCG(K) + p%HubRad ) )*( p%CosPreC(K)**2 )
      END IF

   ENDDO ! K - Blades



   DO K = 1,p%NumBl          ! Loop through the blades


      ! Initialize the generalized blade masses using tip mass effects:

      MBF(K,1,1) = p%TipMass(K)
      MBF(K,2,2) = p%TipMass(K)
      MBE(K,1,1) = p%TipMass(K)


      DO J = 1,p%BldNodes    ! Loop through the blade nodes / elements


      ! Integrate to find the generalized mass of the blade (including tip mass effects).
      !   Ignore the cross-correlation terms of MBF (i.e. MBF(i,j) where i /= j) since
      !   these terms will never be used.

         Shape1 = SHP( p%RNodesNorm(J), p%BldFlexL, p%BldFl1Sh(:,K), 0, ErrStat, ErrMsg )
         Shape2 = SHP( p%RNodesNorm(J), p%BldFlexL, p%BldFl2Sh(:,K), 0, ErrStat, ErrMsg )
         MBF    (K,1,1) = MBF    (K,1,1) + p%BElmntMass(J,K)*Shape1*Shape1
         MBF    (K,2,2) = MBF    (K,2,2) + p%BElmntMass(J,K)*Shape2*Shape2

         Shape  = SHP( p%RNodesNorm(J), p%BldFlexL, p%BldEdgSh(:,K), 0, ErrStat, ErrMsg )
         MBE    (K,1,1) = MBE    (K,1,1) + p%BElmntMass(J,K)*Shape *Shape


      ! Integrate to find the generalized stiffness of the blade (not including centrifugal
      !    effects).

         ElmntStff      = p%StiffBF(K,J)*p%DRNodes(J)                       ! Flapwise stiffness of blade element J
         Shape1 = SHP( p%RNodesNorm(J), p%BldFlexL, p%BldFl1Sh(:,K), 2, ErrStat, ErrMsg )
         Shape2 = SHP( p%RNodesNorm(J), p%BldFlexL, p%BldFl2Sh(:,K), 2, ErrStat, ErrMsg )
         p%KBF    (K,1,1) = p%KBF    (K,1,1) + ElmntStff*Shape1*Shape1
         p%KBF    (K,1,2) = p%KBF    (K,1,2) + ElmntStff*Shape1*Shape2
         p%KBF    (K,2,1) = p%KBF    (K,2,1) + ElmntStff*Shape2*Shape1
         p%KBF    (K,2,2) = p%KBF    (K,2,2) + ElmntStff*Shape2*Shape2

         ElmntStff      = p%StiffBE(K,J)*p%DRNodes(J)                       ! Edgewise stiffness of blade element J
         Shape  = SHP( p%RNodesNorm(J), p%BldFlexL, p%BldEdgSh(:,K), 2, ErrStat, ErrMsg )
         p%KBE    (K,1,1) = p%KBE    (K,1,1) + ElmntStff*Shape *Shape


      ! Integrate to find the centrifugal-term of the generalized flapwise and edgewise
      !   stiffness of the blades.  Ignore the cross-correlation terms of KBFCent (i.e.
      !   KBFCent(i,j) where i /= j) since these terms will never be used.

         ElmntStff      = FMomAbvNd(K,J)*p%DRNodes(J)*p%RotSpeed**2   ! Centrifugal stiffness of blade element J

         Shape1 = SHP( p%RNodesNorm(J), p%BldFlexL, p%BldFl1Sh(:,K), 1, ErrStat, ErrMsg )
         Shape2 = SHP( p%RNodesNorm(J), p%BldFlexL, p%BldFl2Sh(:,K), 1, ErrStat, ErrMsg )
         KBFCent(K,1,1) = KBFCent(K,1,1) + ElmntStff*Shape1*Shape1
         KBFCent(K,2,2) = KBFCent(K,2,2) + ElmntStff*Shape2*Shape2

         Shape  = SHP( p%RNodesNorm(J), p%BldFlexL, p%BldEdgSh(:,K), 1, ErrStat, ErrMsg )
         KBECent(K,1,1) = KBECent(K,1,1) + ElmntStff*Shape *Shape


      ! Calculate the 2nd derivatives of the twisted shape functions:

         Shape  = SHP( p%RNodesNorm(J), p%BldFlexL, p%BldFl1Sh(:,K), 2, ErrStat, ErrMsg )
         p%TwistedSF(K,1,1,J,2) =  Shape*p%CThetaS(K,J)                  ! 2nd deriv. of Phi1(J) for blade K
         p%TwistedSF(K,2,1,J,2) = -Shape*p%SThetaS(K,J)                  ! 2nd deriv. of Psi1(J) for blade K

         Shape  = SHP( p%RNodesNorm(J), p%BldFlexL, p%BldFl2Sh(:,K), 2, ErrStat, ErrMsg )
         p%TwistedSF(K,1,2,J,2) =  Shape*p%CThetaS(K,J)                  ! 2nd deriv. of Phi2(J) for blade K
         p%TwistedSF(K,2,2,J,2) = -Shape*p%SThetaS(K,J)                  ! 2nd deriv. of Psi2(J) for blade K

         Shape  = SHP( p%RNodesNorm(J), p%BldFlexL, p%BldEdgSh(:,K), 2, ErrStat, ErrMsg )
         p%TwistedSF(K,1,3,J,2) =  Shape*p%SThetaS(K,J)                  ! 2nd deriv. of Phi3(J) for blade K
         p%TwistedSF(K,2,3,J,2) =  Shape*p%CThetaS(K,J)                  ! 2nd deriv. of Psi3(J) for blade K


      ! Integrate to find the 1st derivatives of the twisted shape functions:

         DO I = 1,2     ! Loop through Phi and Psi
            DO L = 1,3  ! Loop through all blade DOFs
               TwstdSF     (  I,L,  1) = p%TwistedSF(K,I,L,J,2)*0.5*p%DRNodes(J)
               p%TwistedSF   (K,I,L,J,1) = TwstdSF   ( I,L,  1)
            ENDDO       ! L - All blade DOFs
         ENDDO          ! I - Phi and Psi

         IF ( J /= 1 )  THEN  ! All but the innermost blade element
      ! Add the effects from the (not yet used) portion of element J-1

            DO I = 1,2     ! Loop through Phi and Psi
               DO L = 1,3  ! Loop through all blade DOFs
                  p%TwistedSF(K,I,L,J,1) = p%TwistedSF(K,I,L,J,1) + p%TwistedSF(K,I,L,J-1,1) &
                                       + TwstdSFOld( I,L,  1)
               ENDDO       ! L - All blade DOFs
            ENDDO          ! I - Phi and Psi
         ENDIF


      ! Integrate to find the twisted shape functions themselves (i.e., their zeroeth derivative):

         DO I = 1,2     ! Loop through Phi and Psi
            DO L = 1,3  ! Loop through all blade DOFs
               TwstdSF     (  I,L,  0) = p%TwistedSF(K,I,L,J,1)*0.5*p%DRNodes(J)
               p%TwistedSF   (K,I,L,J,0) = TwstdSF   ( I,L,  0)
            ENDDO       ! L - All blade DOFs
         ENDDO          ! I - Phi and Psi

         IF ( J /= 1 )  THEN  ! All but the innermost blade element
      ! Add the effects from the (not yet used) portion of element J-1

            DO I = 1,2     ! Loop through Phi and Psi
               DO L = 1,3  ! Loop through all blade DOFs
                  p%TwistedSF(K,I,L,J,0) = p%TwistedSF(K,I,L,J,0) + p%TwistedSF(K,I,L,J-1,0) &
                                       + TwstdSFOld( I,L,  0)
               ENDDO       ! L - All blade DOFs
            ENDDO          ! I - Phi and Psi
         ENDIF


      ! Integrate to find the blade axial reduction shape functions:

         DO I = 1,3     ! Loop through all blade DOFs
            DO L = 1,3  ! Loop through all blade DOFs
               AxRdBld    (  I,L  ) = 0.5*p%DRNodes(J)*(                          &
                                      p%TwistedSF(K,1,I,J,1)*p%TwistedSF(K,1,L,J,1) &
                                    + p%TwistedSF(K,2,I,J,1)*p%TwistedSF(K,2,L,J,1) )
               p%AxRedBld   (K,I,L,J) = AxRdBld(I,L)
            ENDDO       ! L - All blade DOFs
         ENDDO          ! I - All blade DOFs

         IF ( J /= 1 )  THEN  ! All but the innermost blade element
      ! Add the effects from the (not yet used) portion of element J-1

            DO I = 1,3     ! Loop through all blade DOFs
               DO L = 1,3  ! Loop through all blade DOFs
                  p%AxRedBld(K,I,L,J) = p%AxRedBld(K,I,L,J) + p%AxRedBld(K,I,L,J-1)   &
                                    + AxRdBldOld(I,L)
               ENDDO       ! L - All blade DOFs
            ENDDO          ! I - All blade DOFs
         ENDIF


      ! Store the TwstdSF and AxRdBld terms of the current element (these will be used for the next element)

         TwstdSFOld = TwstdSF
         AxRdBldOld = AxRdBld


      ENDDO ! J - Blade nodes / elements




   IF (p%BD4Blades) THEN

      !p%KBF     ( K,:,:    ) = 0.0_ReKi
      
         ! the 1st and zeroeth derivatives of the twisted shape functions at the blade root:
      p%TwistedSF(K,:,:,:,1) = 0.0_ReKi
      p%TwistedSF(K,:,:,:,0) = 0.0_ReKi 
      p%AxRedBld( K,:,:,:  ) = 0.0_ReKi
   ELSE
            
      ! Apply the flapwise modal stiffness tuners of the blades to KBF():

      DO I = 1,2     ! Loop through flap DOFs
         DO L = 1,2  ! Loop through flap DOFs
            p%KBF(K,I,L) = SQRT( p%FStTunr(K,I)*p%FStTunr(K,L) )*p%KBF(K,I,L)
         ENDDO       ! L - Flap DOFs
      ENDDO          ! I - Flap DOFs
      
      ! Calculate the blade natural frequencies:
      
      DO I = 1,NumBF     ! Loop through flap DOFs
         p%FreqBF(K,I,1) = Inv2Pi*SQRT(   p%KBF(K,I,I)                   /( MBF(K,I,I) - p%TipMass(K) ) )   ! Natural blade I-flap frequency w/o centrifugal stiffening nor     tip mass effects
         p%FreqBF(K,I,2) = Inv2Pi*SQRT(   p%KBF(K,I,I)                   /  MBF(K,I,I)                )     ! Natural blade I-flap frequency w/o centrifugal stiffening, but w/ tip mass effects
         p%FreqBF(K,I,3) = Inv2Pi*SQRT( ( p%KBF(K,I,I) + KBFCent(K,I,I) )/  MBF(K,I,I)                )     ! Natural blade I-flap frequency w/  centrifugal stiffening and     tip mass effects
      ENDDO          ! I - Flap DOFs

      p%FreqBE   (K,1,1) = Inv2Pi*SQRT(   p%KBE(K,1,1)                   /( MBE(K,1,1) - p%TipMass(K) ) )   ! Natural blade 1-edge frequency w/o centrifugal stiffening nor      tip mass effects
      p%FreqBE   (K,1,2) = Inv2Pi*SQRT(   p%KBE(K,1,1)                   /  MBE(K,1,1)                )     ! Natural Blade 1-edge frequency w/o  centrifugal stiffening, but w/ tip mass effects
      p%FreqBE   (K,1,3) = Inv2Pi*SQRT( ( p%KBE(K,1,1) + KBECent(K,1,1) )/  MBE(K,1,1)                )     ! Natural Blade 1-edge frequency w/  centrifugal stiffening and      tip mass effects


      ! Calculate the generalized damping of the blades:

      DO I = 1,NumBF     ! Loop through flap DOFs
         DO L = 1,NumBF  ! Loop through flap DOFs
            p%CBF(K,I,L) = ( 0.01*p%BldFDamp(K,L) )*p%KBF(K,I,L)/( Pi*p%FreqBF(K,L,1) )
         ENDDO       ! L - Flap DOFs
      ENDDO          ! I - Flap DOFs

      p%CBE      (K,1,1) = ( 0.01*p%BldEDamp(K,1) )*p%KBE(K,1,1)/( Pi*p%FreqBE(K,1,1) )


      ! Calculate the 2nd derivatives of the twisted shape functions at the blade root:

      Shape  = SHP( 0.0_ReKi, p%BldFlexL, p%BldFl1Sh(:,K), 2, ErrStat, ErrMsg )
      p%TwistedSF(K,1,1,0,2) =  Shape*p%CThetaS(K,0)        ! 2nd deriv. of Phi1(0) for blade K
      p%TwistedSF(K,2,1,0,2) = -Shape*p%SThetaS(K,0)        ! 2nd deriv. of Psi1(0) for blade K

      Shape  = SHP( 0.0_ReKi, p%BldFlexL, p%BldFl2Sh(:,K), 2, ErrStat, ErrMsg )
      p%TwistedSF(K,1,2,0,2) =  Shape*p%CThetaS(K,0)        ! 2nd deriv. of Phi2(0) for blade K
      p%TwistedSF(K,2,2,0,2) = -Shape*p%SThetaS(K,0)        ! 2nd deriv. of Psi2(0) for blade K

      Shape  = SHP( 0.0_ReKi, p%BldFlexL, p%BldEdgSh(:,K), 2, ErrStat, ErrMsg )
      p%TwistedSF(K,1,3,0,2) =  Shape*p%SThetaS(K,0)        ! 2nd deriv. of Phi3(0) for blade K
      p%TwistedSF(K,2,3,0,2) =  Shape*p%CThetaS(K,0)        ! 2nd deriv. of Psi3(0) for blade K
      
      
      ! Calculate the 2nd derivatives of the twisted shape functions at the tip:

      Shape  = SHP( 1.0_ReKi, p%BldFlexL, p%BldFl1Sh(:,K), 2, ErrStat, ErrMsg )
      p%TwistedSF(K,1,1,p%TipNode,2) =  Shape*p%CThetaS(K,p%TipNode)        ! 2nd deriv. of Phi1(p%TipNode) for blade K
      p%TwistedSF(K,2,1,p%TipNode,2) = -Shape*p%SThetaS(K,p%TipNode)        ! 2nd deriv. of Psi1(p%TipNode) for blade K

      Shape  = SHP( 1.0_ReKi, p%BldFlexL, p%BldFl2Sh(:,K), 2, ErrStat, ErrMsg )
      p%TwistedSF(K,1,2,p%TipNode,2) =  Shape*p%CThetaS(K,p%TipNode)        ! 2nd deriv. of Phi2(p%TipNode) for blade K
      p%TwistedSF(K,2,2,p%TipNode,2) = -Shape*p%SThetaS(K,p%TipNode)        ! 2nd deriv. of Psi2(p%TipNode) for blade K

      Shape  = SHP( 1.0_ReKi, p%BldFlexL, p%BldEdgSh(:,K), 2, ErrStat, ErrMsg )
      p%TwistedSF(K,1,3,p%TipNode,2) =  Shape*p%SThetaS(K,p%TipNode)        ! 2nd deriv. of Phi3(p%TipNode) for blade K
      p%TwistedSF(K,2,3,p%TipNode,2) =  Shape*p%CThetaS(K,p%TipNode)        ! 2nd deriv. of Psi3(p%TipNode) for blade K


      ! Integrate to find the 1st and zeroeth derivatives of the twisted shape functions
      !   at the tip:

      DO I = 1,2     ! Loop through Phi and Psi
         DO L = 1,3  ! Loop through all blade DOFs
            p%TwistedSF(K,I,L,p%TipNode,1) = p%TwistedSF(K,I,L,p%BldNodes,1) + TwstdSFOld(I,L,1)
            p%TwistedSF(K,I,L,p%TipNode,0) = p%TwistedSF(K,I,L,p%BldNodes,0) + TwstdSFOld(I,L,0)
         ENDDO       ! L - All blade DOFs
      ENDDO          ! I - Phi and Psi

         ! the 1st and zeroeth derivatives of the twisted shape functions at the blade root:
      p%TwistedSF(K,:,:,0,1) = 0.0_ReKi
      p%TwistedSF(K,:,:,0,0) = 0.0_ReKi 
      p%AxRedBld( K,:,:,0  ) = 0.0_ReKi

      ! Integrate to find the blade axial reduction shape functions at the tip:

      DO I = 1,3     ! Loop through all blade DOFs
         DO L = 1,3  ! Loop through all blade DOFs
            p%AxRedBld(K,I,L,p%TipNode) = p%AxRedBld(K,I,L,p%BldNodes) + AxRdBldOld(I,L)
         ENDDO       ! L - All blade DOFs
      ENDDO          ! I - All blade DOFs
   END IF ! p%BD4Blades


   ENDDO ! K - Blades



      ! Calculate the tower-top mass:

   p%TwrTpMass = p%RotMass + p%RFrlMass + p%BoomMass + p%TFinMass + p%NacMass + p%YawBrMass


   DO J = p%TwrNodes,1,-1 ! Loop through the tower nodes / elements in reverse


      ! Calculate the mass of the current element

      p%TElmntMass(J)    = p%MassT(J)*abs(p%DHNodes(J))     ! Mass of tower element J


      ! Integrate to find the tower mass which will be output in .fsm

      p%TwrMass      = p%TwrMass + p%TElmntMass(J)


      ! Integrate to find TMssAbvNd:

      TMssAbvNd   (J) = 0.5*p%TElmntMass(J)

      IF ( J == p%TwrNodes )  THEN ! Uppermost tower element
      ! Add the TwrTpMass effects:

         TMssAbvNd(J) = TMssAbvNd(J) + p%TwrTpMass
      ELSE                       ! All other tower elements
      ! Add to TMssAbvNd(J) the effects from the (not yet used) portion of element J+1

         TMssAbvNd(J) = 0.5*p%TElmntMass(J+1) + TMssAbvNd(J) + TMssAbvNd(J+1)
      ENDIF


   ENDDO ! J - Tower nodes / elements in reverse



      ! Initialize the generalized tower masses using tower-top mass effects:

   DO I = 1,2  ! Loop through all tower modes in a single direction
      MTFA(I,I) = p%TwrTpMass
      MTSS(I,I) = p%TwrTpMass
   ENDDO       ! I - All tower modes in a single direction

      ! set values for tower base (note that we haven't corrctly defined the values for (:,0,2) in the arrays below):
   p%TwrFASF(   :,0,0:1) = 0.0_ReKi
   p%TwrSSSF(   :,0,0:1) = 0.0_ReKi
   p%AxRedTFA(:,:,0)     = 0.0_ReKi
   p%AxRedTSS(:,:,0)     = 0.0_ReKi
   
   DO J = 1,p%TwrNodes    ! Loop through the tower nodes / elements


      ! Calculate the tower shape functions (all derivatives):

      p%TwrFASF(1,J,2) = SHP( p%HNodesNorm(J), p%TwrFlexL, InputFileData%TwFAM1Sh(:), 2, ErrStat, ErrMsg )
      p%TwrFASF(2,J,2) = SHP( p%HNodesNorm(J), p%TwrFlexL, InputFileData%TwFAM2Sh(:), 2, ErrStat, ErrMsg )
      p%TwrFASF(1,J,1) = SHP( p%HNodesNorm(J), p%TwrFlexL, InputFileData%TwFAM1Sh(:), 1, ErrStat, ErrMsg )
      p%TwrFASF(2,J,1) = SHP( p%HNodesNorm(J), p%TwrFlexL, InputFileData%TwFAM2Sh(:), 1, ErrStat, ErrMsg )
      p%TwrFASF(1,J,0) = SHP( p%HNodesNorm(J), p%TwrFlexL, InputFileData%TwFAM1Sh(:), 0, ErrStat, ErrMsg )
      p%TwrFASF(2,J,0) = SHP( p%HNodesNorm(J), p%TwrFlexL, InputFileData%TwFAM2Sh(:), 0, ErrStat, ErrMsg )

      p%TwrSSSF(1,J,2) = SHP( p%HNodesNorm(J), p%TwrFlexL, InputFileData%TwSSM1Sh(:), 2, ErrStat, ErrMsg )
      p%TwrSSSF(2,J,2) = SHP( p%HNodesNorm(J), p%TwrFlexL, InputFileData%TwSSM2Sh(:), 2, ErrStat, ErrMsg )
      p%TwrSSSF(1,J,1) = SHP( p%HNodesNorm(J), p%TwrFlexL, InputFileData%TwSSM1Sh(:), 1, ErrStat, ErrMsg )
      p%TwrSSSF(2,J,1) = SHP( p%HNodesNorm(J), p%TwrFlexL, InputFileData%TwSSM2Sh(:), 1, ErrStat, ErrMsg )
      p%TwrSSSF(1,J,0) = SHP( p%HNodesNorm(J), p%TwrFlexL, InputFileData%TwSSM1Sh(:), 0, ErrStat, ErrMsg )
      p%TwrSSSF(2,J,0) = SHP( p%HNodesNorm(J), p%TwrFlexL, InputFileData%TwSSM2Sh(:), 0, ErrStat, ErrMsg )


      ! Integrate to find the generalized mass of the tower (including tower-top mass effects).
      !   Ignore the cross-correlation terms of MTFA (i.e. MTFA(i,j) where i /= j) and MTSS
      !   since these terms will never be used.


      DO I = 1,2     ! Loop through all tower DOFs in one direction
         MTFA  (I,I) = MTFA  (I,I) + p%TElmntMass(J)*p%TwrFASF(I,J,0)**2
         MTSS  (I,I) = MTSS  (I,I) + p%TElmntMass(J)*p%TwrSSSF(I,J,0)**2
      ENDDO          ! I - through all tower DOFs in one direction


      ! Integrate to find the generalized stiffness of the tower (not including gravitational
      !    effects).

      ElStffFA       = p%StiffTFA(J)*abs(p%DHNodes(J))                        ! Fore-aft stiffness of tower element J
      ElStffSS       = p%StiffTSS(J)*abs(p%DHNodes(J))                        ! Side-to-side stiffness of tower element J

      DO I = 1,2     ! Loop through all tower DOFs in one direction
         DO L = 1,2  ! Loop through all tower DOFs in one direction
            p%KTFA (I,L) = p%KTFA    (I,L) + ElStffFA *p%TwrFASF(I,J,2)*p%TwrFASF(L,J,2)
            p%KTSS (I,L) = p%KTSS    (I,L) + ElStffSS *p%TwrSSSF(I,J,2)*p%TwrSSSF(L,J,2)
         ENDDO       ! L - All tower DOFs in one direction
      ENDDO          ! I - through all tower DOFs in one direction


      ! Integrate to find the gravitational-term of the generalized stiffness of the tower.
      !   Ignore the cross-correlation terms of KTFAGrav (i.e. KTFAGrav(i,j) where i /= j)
      !   and KTSSGrav since these terms will never be used.

      ElmntStff      = -TMssAbvNd(J)*abs(p%DHNodes(J))*p%Gravity              ! Gravitational stiffness of tower element J

      DO I = 1,2     ! Loop through all tower DOFs in one direction
         KTFAGrav(I,I) = KTFAGrav(I,I) + ElmntStff*p%TwrFASF(I,J,1)**2
         KTSSGrav(I,I) = KTSSGrav(I,I) + ElmntStff*p%TwrSSSF(I,J,1)**2
      ENDDO


      ! Integrate to find the tower axial reduction shape functions:

      DO I = 1,2     ! Loop through all tower DOFs in one direction
         DO L = 1,2  ! Loop through all tower DOFs in one direction
            AxRdTFA (I,L) = 0.5*p%DHNodes(J)*p%TwrFASF(I,J,1)*p%TwrFASF(L,J,1)
            AxRdTSS (I,L) = 0.5*p%DHNodes(J)*p%TwrSSSF(I,J,1)*p%TwrSSSF(L,J,1)

            p%AxRedTFA(I,L,J) = AxRdTFA(I,L)
            p%AxRedTSS(I,L,J) = AxRdTSS(I,L)
         ENDDO       ! L - All tower DOFs in one direction
      ENDDO

      IF ( J /= 1 )  THEN  ! All but the lowermost tower element
      ! Add the effects from the (not yet used) portion of element J-1

         DO I = 1,2     ! Loop through all tower DOFs in one direction
            DO L = 1,2  ! Loop through all tower DOFs in one direction
               p%AxRedTFA(I,L,J) = p%AxRedTFA(I,L,J) + p%AxRedTFA(I,L,J-1)+ AxRdTFAOld(I,L)
               p%AxRedTSS(I,L,J) = p%AxRedTSS(I,L,J) + p%AxRedTSS(I,L,J-1)+ AxRdTSSOld(I,L)
            ENDDO       ! L - All tower DOFs in one direction
         ENDDO
      ENDIF


      ! Store the AxRdTFA and AxRdTSS terms of the current element (these will be used for the next element)

      AxRdTFAOld = AxRdTFA
      AxRdTSSOld = AxRdTSS


   ENDDO ! J - Tower nodes / elements


   ! Apply the modal stiffness tuners of the tower to KTFA() and KTSS():

   DO I = 1,2     ! Loop through all tower DOFs in one direction
      DO L = 1,2  ! Loop through all tower DOFs in one direction
         p%KTFA(I,L) = SQRT( InputFileData%FAStTunr(I)*InputFileData%FAStTunr(L) )*p%KTFA(I,L)

         p%KTSS(I,L) = SQRT( InputFileData%SSStTunr(I)*InputFileData%SSStTunr(L) )*p%KTSS(I,L)
      ENDDO       ! L - All tower DOFs in one direction
   ENDDO          ! I - through all tower DOFs in one direction


      ! Calculate the tower natural frequencies:

   DO I = 1,2     ! Loop through all tower DOFs in one direction
      if ( EqualRealNos(( MTFA(I,I) - p%TwrTpMass ), 0.0_ReKi) ) then
         p%FreqTFA(I,1) = NaN ! Avoid creating a divide by zero signal, but set p%FreqTFA(I,1) = NaN
      else        
         p%FreqTFA(I,1) = Inv2Pi*SQRT(   p%KTFA(I,I)/( MTFA(I,I) - p%TwrTpMass ) )  ! Natural tower I-fore-aft frequency w/o gravitational destiffening nor tower-top mass effects
      end if
      if ( EqualRealNos(( MTSS(I,I) - p%TwrTpMass ), 0.0_ReKi) ) then
         p%FreqTSS(I,1) = NaN ! Avoid creating a divide by zero signal, but set p%FreqTFS(I,1) = NaN
      else
         p%FreqTSS(I,1) = Inv2Pi*SQRT(   p%KTSS(I,I)/( MTSS(I,I) - p%TwrTpMass ) )  ! Natural tower I-side-to-side frequency w/o gravitational destiffening nor tower-top mass effects
      end if
      p%FreqTFA(I,2) = Inv2Pi*SQRT( ( p%KTFA(I,I) + KTFAGrav(I,I) )/MTFA(I,I)               )  ! Natural tower I-fore-aft frequency w/  gravitational destiffening and tower-top mass effects
      p%FreqTSS(I,2) = Inv2Pi*SQRT( ( p%KTSS(I,I) + KTSSGrav(I,I) )/MTSS(I,I)               )  ! Natural tower I-side-to-side frequency w/  gravitational destiffening and tower-top mass effects
   ENDDO          ! I - All tower DOFs in one direction


      ! Calculate the generalized damping of the tower:

   DO I = 1,2     ! Loop through all tower DOFs in one direction
      DO L = 1,2  ! Loop through all tower DOFs in one direction
         p%CTFA(I,L) = ( 0.01*InputFileData%TwrFADmp(L) )*p%KTFA(I,L)/( Pi*p%FreqTFA(L,1) )

         p%CTSS(I,L) = ( 0.01*InputFileData%TwrSSDmp(L) )*p%KTSS(I,L)/( Pi*p%FreqTSS(L,1) )
      ENDDO       ! L - All tower DOFs in one direction
   ENDDO          ! I - All tower DOFs in one direction


      ! Calculate the tower shape functions (all derivatives) at the tower-top:

   p%TwrFASF(1,p%TTopNode,2) = SHP( 1.0_ReKi, p%TwrFlexL, InputFileData%TwFAM1Sh(:), 2, ErrStat, ErrMsg )
   p%TwrFASF(2,p%TTopNode,2) = SHP( 1.0_ReKi, p%TwrFlexL, InputFileData%TwFAM2Sh(:), 2, ErrStat, ErrMsg )
   p%TwrFASF(1,p%TTopNode,1) = SHP( 1.0_ReKi, p%TwrFlexL, InputFileData%TwFAM1Sh(:), 1, ErrStat, ErrMsg )
   p%TwrFASF(2,p%TTopNode,1) = SHP( 1.0_ReKi, p%TwrFlexL, InputFileData%TwFAM2Sh(:), 1, ErrStat, ErrMsg )
   p%TwrFASF(1,p%TTopNode,0) = SHP( 1.0_ReKi, p%TwrFlexL, InputFileData%TwFAM1Sh(:), 0, ErrStat, ErrMsg )
   p%TwrFASF(2,p%TTopNode,0) = SHP( 1.0_ReKi, p%TwrFlexL, InputFileData%TwFAM2Sh(:), 0, ErrStat, ErrMsg )

   p%TwrSSSF(1,p%TTopNode,2) = SHP( 1.0_ReKi, p%TwrFlexL, InputFileData%TwSSM1Sh(:), 2, ErrStat, ErrMsg )
   p%TwrSSSF(2,p%TTopNode,2) = SHP( 1.0_ReKi, p%TwrFlexL, InputFileData%TwSSM2Sh(:), 2, ErrStat, ErrMsg )
   p%TwrSSSF(1,p%TTopNode,1) = SHP( 1.0_ReKi, p%TwrFlexL, InputFileData%TwSSM1Sh(:), 1, ErrStat, ErrMsg )
   p%TwrSSSF(2,p%TTopNode,1) = SHP( 1.0_ReKi, p%TwrFlexL, InputFileData%TwSSM2Sh(:), 1, ErrStat, ErrMsg )
   p%TwrSSSF(1,p%TTopNode,0) = SHP( 1.0_ReKi, p%TwrFlexL, InputFileData%TwSSM1Sh(:), 0, ErrStat, ErrMsg )
   p%TwrSSSF(2,p%TTopNode,0) = SHP( 1.0_ReKi, p%TwrFlexL, InputFileData%TwSSM2Sh(:), 0, ErrStat, ErrMsg )


      ! Integrate to find the tower axial reduction shape functions at the tower-top:

   DO I = 1,2     ! Loop through all tower DOFs in one direction
      DO L = 1,2  ! Loop through all tower DOFs in one direction
         p%AxRedTFA(I,L,p%TTopNode) = p%AxRedTFA(I,L,p%TwrNodes)+ AxRdTFAOld(I,L)
         p%AxRedTSS(I,L,p%TTopNode) = p%AxRedTSS(I,L,p%TwrNodes)+ AxRdTSSOld(I,L)
      ENDDO       ! L - All tower DOFs in one direction
   ENDDO


      ! Calculate the turbine mass:

   p%TurbMass  = p%TwrTpMass + p%TwrMass


   RETURN
END SUBROUTINE Coeff
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine calculates the initial blade deflections.
!! Base the intial values of the blade DOFs, INITQF1, INITQF2, and
!!   INITQE1, on OoPDefl and IPDefl.
!! Write messages to the screen if the specified initial tip displacements
!!  are incompatible with the enabled DOFs.
SUBROUTINE InitBlDefl ( p, InputFileData, InitQF1, InitQF2, InitQE1, ErrStat, ErrMsg )
!..................................................................................................................................


      ! Passed variables:
   TYPE(ED_ParameterType),  INTENT(IN)  :: p                                       !< parameters of the structural dynamics module
   TYPE(ED_InputFile),      INTENT(IN)  :: InputFileData                           !< all the data in the ElastoDyn input file

   REAL(ReKi),              INTENT(OUT) :: InitQE1(p%NumBl)                        !< Initial edge deflection (output).
   REAL(ReKi),              INTENT(OUT) :: InitQF1(p%NumBl)                        !< Initial flap deflection for mode 1 (output).
   REAL(ReKi),              INTENT(OUT) :: InitQF2(p%NumBl)                        !< Initial flap deflection for mode 2 (output).

   INTEGER(IntKi),          INTENT(OUT) :: ErrStat                                 !< Error status
   CHARACTER(ErrMsgLen),    INTENT(OUT) :: ErrMsg                                  !< Error message when ErrStat =/ ErrID_None


      ! Local variables:
   REAL(ReKi)                   :: A(2,3)                                          ! Augmented matrix for solution of initial deflections.
   REAL(ReKi)                   :: CosPitch                                        ! Cosine of the pitch for this blade.
   REAL(ReKi)                   :: Det                                             ! Determinate of right-hand side of A.
   REAL(ReKi)                   :: SinPitch                                        ! Sine of the pitch for this blade.
   REAL(ReKi)                   :: TotResid                                        ! Generator torque.

   INTEGER(IntKi)               :: K                                               ! Blade number

      ! some warning messages
   CHARACTER(*), PARAMETER      :: Approx   = ' An approximate characterization of the specified blade deflection will be made.'
   CHARACTER(*), PARAMETER      :: BadIP    = ' Initial blade in-plane tip displacement will be ignored.'
   CHARACTER(*), PARAMETER      :: BadOoP   = ' Initial blade out-of-plane tip displacement will be ignored.'
   CHARACTER(*), PARAMETER      :: Ignore   = ' All initial blade tip displacements will be ignored.'


      ! Initialize variables
   ErrStat = ErrID_None
   ErrMsg  = ''

   InitQE1 = 0.0
   InitQF1 = 0.0
   InitQF2 = 0.0
   !bjj: replace InitQF1 and InitQF2 with an array to avoid so much duplication of logic here...

   DO K=1,p%NumBl

         ! Calculate the array of deflections(???).

      CosPitch = COS( InputFileData%BlPitch(K) )
      SinPitch = SIN( InputFileData%BlPitch(K) )

      A(1,2) =  p%TwistedSF(K,1,3,p%TipNode,0)*CosPitch + p%TwistedSF(K,2,3,p%TipNode,0)*SinPitch
      A(2,2) = -p%TwistedSF(K,1,3,p%TipNode,0)*SinPitch + p%TwistedSF(K,2,3,p%TipNode,0)*CosPitch
      A(1,3) =  InputFileData%OoPDefl
      A(2,3) =  InputFileData%IPDefl

      IF ( InputFileData%FlapDOF1 )  THEN                                ! Blade flap mode 1 is enabled

         A(1,1) =  p%TwistedSF(K,1,1,p%TipNode,0)*CosPitch + p%TwistedSF(K,2,1,p%TipNode,0)*SinPitch
         A(2,1) = -p%TwistedSF(K,1,1,p%TipNode,0)*SinPitch + p%TwistedSF(K,2,1,p%TipNode,0)*CosPitch

         DET = ( A(1,1)*A(2,2) - A(1,2)*A(2,1) )

         IF ( .NOT. EqualRealNos( DET, 0.0_ReKi ) ) THEN                  ! Apply all flap deflection to mode 1

            InitQF1(K) = ( A(1,3)*A(2,2) - A(1,2)*A(2,3) )/DET
            InitQE1(K) = ( A(1,1)*A(2,3) - A(1,3)*A(2,1) )/DET

         ELSEIF ( .NOT. InputFileData%EdgeDOF )  THEN                     ! Blade edge mode 1 is not enabled which caused DET = 0.

            InitQE1(K) = 0.0

            IF ( .NOT. EqualRealNos( A(1,1), 0.0_ReKi ) )  THEN
               IF ( .NOT. EqualRealNos( A(2,1), 0.0_ReKi ) )  THEN        ! Find a solution of the 2 equations in 1 variable that
                                                                          !  minimizes the sum of the squares of the equation's residuals.

                  InitQF1(K) = ( A(1,1)*A(1,3) + A(2,1)*A(2,3) )/( A(1,1)**2 + A(2,1)**2 )

                  TotResid = SQRT( ( A(1,1)*InitQF1(K) - A(1,3) )**2 + ( A(2,1)*InitQF1(K) - A(2,3) )**2 )

                  IF ( .NOT. EqualRealNos( TotResid, 0.0_ReKi ) ) THEN
                     CALL CheckError( ErrID_Warn, Approx )
                  ENDIF

               ELSE !A(1,1) /= 0; A(2,1) == 0

                  InitQF1(K) = A(1,3)/A(1,1)

                  IF ( .NOT. EqualRealNos( InputFileData%IPDefl,  0.0_ReKi ) )  THEN
                     CALL CheckError( ErrID_Warn, BadIP )
                  ENDIF

               ENDIF

            ELSE ! A(1,1) == 0

               IF ( .NOT. EqualRealNos( InputFileData%OoPDefl, 0.0_ReKi ) ) THEN
                  CALL CheckError( ErrID_Warn, BadOoP )
               END IF

               IF ( .NOT. EqualRealNos( A(2,1), 0.0_ReKi ) )   THEN
                  InitQF1(K) = A(2,3)/A(2,1)
               ELSE
                  InitQF1(K) = 0.0

                  IF ( .NOT. EqualRealNos( InputFileData%IPDefl,  0.0_ReKi ) )  THEN
                     CALL CheckError( ErrID_Warn, BadIP )
                  ENDIF

               ENDIF
            ENDIF

         ELSE                                     ! It is impossible to find any "good" solution, so ignore the initial tip displacements

            InitQF1(K) = 0.0
            InitQE1(K) = 0.0

            IF ( ( InputFileData%OoPDefl /= 0.0 ) .OR. ( InputFileData%IPDefl /= 0.0 ) )  THEN
               CALL CheckError( ErrID_Warn, Ignore )
            ENDIF

         ENDIF

      ELSE                                        ! Blade flap mode 1 is not enabled.

         InitQF1(K) = 0.0

         IF ( InputFileData%FlapDOF2 )  THEN                    ! Blade flap mode 2 is enabled.

            A(1,1) =  p%TwistedSF(K,1,2,p%TipNode,0)*CosPitch + p%TwistedSF(K,2,2,p%TipNode,0)*SinPitch
            A(2,1) = -p%TwistedSF(K,1,2,p%TipNode,0)*SinPitch + p%TwistedSF(K,2,2,p%TipNode,0)*CosPitch

            DET = ( A(1,1)*A(2,2) - A(1,2)*A(2,1) )

            IF ( .NOT. EqualRealNos( DET, 0.0_ReKi ) ) THEN      ! Apply all flap deflection to mode 2
               InitQF2 = ( A(1,3)*A(2,2) - A(1,2)*A(2,3) )/DET
               InitQE1 = ( A(1,1)*A(2,3) - A(1,3)*A(2,1) )/DET

            ELSEIF ( .NOT. InputFileData%EdgeDOF )  THEN          ! Blade edge mode 1 is not enabled which caused DET = 0.

               InitQE1(K) = 0.0

               IF ( .NOT. EqualRealNos( A(1,1), 0.0_ReKi ) )  THEN
                  IF ( .NOT. EqualRealNos( A(2,1), 0.0_ReKi ) )   THEN      ! Find a solution of the 2 equations in 1 variable that
                                                                            !  minimizes the sum of the squares of the equation's residuals
                     InitQF2(K) = ( A(1,1)*A(1,3) + A(2,1)*A(2,3) )/( A(1,1)**2 + A(2,1)**2 )

                     TotResid = SQRT( ( A(1,1)*InitQF2(K) - A(1,3))**2 + ( A(2,1)*InitQF2(K) - A(2,3) )**2 )

                     IF ( .NOT. EqualRealNos( TotResid, 0.0_ReKi ) )  THEN
                        CALL CheckError( ErrID_Warn, Approx )
                     ENDIF
                  ELSE
                     InitQF2(K) = A(1,3)/A(1,1)

                     IF ( .NOT. EqualRealNos( InputFileData%IPDefl, 0.0_ReKi ) ) THEN
                        CALL CheckError( ErrID_Warn, BadIP )
                     ENDIF
                  ENDIF
               ELSE
                  IF ( .NOT. EqualRealNos( InputFileData%OoPDefl, 0.0_ReKi ) ) THEN
                     CALL CheckError( ErrID_Warn, BadOoP )
                  END IF

                  IF ( .NOT. EqualRealNos( A(2,1), 0.0_ReKi ) )  THEN
                     InitQF2(K) = A(2,3)/A(2,1)
                  ELSE
                     InitQF2(K) = 0.0

                     IF ( .NOT. EqualRealNos( InputFileData%IPDefl, 0.0_ReKi ) )  THEN
                        CALL CheckError( ErrID_Warn, BadIP )
                     ENDIF

                  ENDIF
               ENDIF

            ELSE                                  ! It is impossible to find any "good" solution, so ignore
                                                  ! the initial tip displacements.
               InitQF2(K) = 0.0
               InitQE1(K) = 0.0

               IF ( .NOT. EqualRealNos( InputFileData%OoPDefl,  0.0_ReKi ) .OR. &
                    .NOT. EqualRealNos( InputFileData%IPDefl,   0.0_ReKi ) )  THEN
                  CALL CheckError( ErrID_Warn, Ignore )
                ENDIF
            ENDIF

         ELSE                                     ! Blade flap mode 2 is not enabled.

            InitQF2(K) = 0.0

            IF ( .NOT. EqualRealNos( A(1,2), 0.0_ReKi ) )  THEN

               IF ( .NOT. EqualRealNos( A(2,2), 0.0_ReKi ) )  THEN         ! Find a solution of the 2 equations in 1 variable that minimizes
                                                                           !  the sum of the squares of the equation's residuals.
                  InitQE1(K) = ( A(1,2)*A(1,3) + A(2,2)*A(2,3) )/( A(1,2)**2 + A(2,2)**2 )

                  TotResid = SQRT( ( A(1,2)*InitQE1(K) - A(1,3) )**2 + ( A(2,2)*InitQE1(K) - A(2,3) )**2)

                  IF ( .NOT. EqualRealNos( TotResid, 0.0_ReKi ) )  THEN
                     CALL CheckError( ErrID_Warn, Approx )
                  ENDIF

               ELSE

                  InitQE1(K) = A(1,3)/A(1,2)

                  IF ( .NOT. EqualRealNos( InputFileData%IPDefl, 0.0_ReKi ) )  THEN
                     CALL CheckError( ErrID_Warn, BadIP )
                  ENDIF

               ENDIF

            ELSE

               IF ( .NOT. EqualRealNos( InputFileData%OoPDefl, 0.0_ReKi ) ) THEN
                  CALL CheckError( ErrID_Warn, BadOoP )
               END IF

               IF ( .NOT. EqualRealNos( A(2,2), 0.0_ReKi ) )  THEN
                  InitQE1(K) = A(2,3)/A(2,2)
               ELSE
                  InitQE1(K) = 0.0

                  IF ( .NOT. EqualRealNos( InputFileData%IPDefl,  0.0_ReKi ) )  THEN
                     CALL CheckError( ErrID_Warn, BadIP )
                  ENDIF

               ENDIF

            ENDIF

         ENDIF

      ENDIF

   END DO !K

   RETURN
CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'InitBlDefl:Blade '//TRIM(Num2LStr(K))// &
                     ' initial blade tip displacements are Incompat with enabled DOFs: '//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
         END IF

      END IF


   END SUBROUTINE CheckError

END SUBROUTINE InitBlDefl
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is used create arrays of DOF indices (pointers / (vector susbscript arrays) that contribute to the QD2T-related
!!   linear accelerations of various points within the system in the inertia frame, based on which DOFs are presently enabled.
!! NOTE: The order in which the DOFs are tested within this routine and hence the order in which the DOF indices appear in the
!!       vector subscript arrays, determines the order in which the states will appear in the linearized model created by FAST
!!       when AnalMode == 2.  This order is not necessarily sorted from smallest to largest DOF index.
!! bjj: note that this routine is now called only in the initialization routine. It is not available during time simulation.
SUBROUTINE SetEnabledDOFIndexArrays( p )
!----------------------------------------------------------------------------------------------------------------------------------

      ! passed variables
   TYPE(ED_ParameterType), INTENT(INOUT)   :: p                                  !< Parameters of the structural dynamics module

      ! Local Variables:
   INTEGER(IntKi)                   :: I                                          ! Loops through all DOFs.
   INTEGER(IntKi)                   :: K                                          ! Loops through blades.



      ! Initialize total counts to zero.

   p%DOFs%NActvDOF = 0
   p%DOFs%NPCE     = 0
   p%DOFs%NPDE     = 0
   p%DOFs%NPIE     = 0
   p%DOFs%NPTTE    = 0
   p%DOFs%NPTE     = 0
   p%DOFs%NPSBE(:) = 0
   p%DOFs%NPSE (:) = 0
   p%DOFs%NPUE     = 0
   p%DOFs%NPYE     = 0


      ! Test each DOF and include the appropriate indices in the subscript arrays
      !  and total counts:

   IF ( p%DOF_Flag(DOF_Sg  ) )  THEN  ! Platform surge.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPDE     = p%DOFs%NPDE     + 1
      p%DOFs%NPIE     = p%DOFs%NPIE     + 1
      p%DOFs%NPTE     = p%DOFs%NPTE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1
      p%DOFs%NPUE     = p%DOFs%NPUE     + 1
      p%DOFs%NPYE     = p%DOFs%NPYE     + 1

      p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_Sg
      p%DOFs%PCE     (  p%DOFs%NPCE    ) = DOF_Sg
      p%DOFs%PDE     (  p%DOFs%NPDE    ) = DOF_Sg
      p%DOFs%PIE     (  p%DOFs%NPIE    ) = DOF_Sg
      p%DOFs%PTE     (  p%DOFs%NPTE    ) = DOF_Sg
      p%DOFs%PSE     (:,p%DOFs%NPSE (:)) = DOF_Sg
      p%DOFs%PUE     (  p%DOFs%NPUE    ) = DOF_Sg
      p%DOFs%PYE     (  p%DOFs%NPYE    ) = DOF_Sg

   ENDIF


   IF ( p%DOF_Flag(DOF_Sw  ) )  THEN  ! Platform sway.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPDE     = p%DOFs%NPDE     + 1
      p%DOFs%NPIE     = p%DOFs%NPIE     + 1
      p%DOFs%NPTE     = p%DOFs%NPTE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1
      p%DOFs%NPUE     = p%DOFs%NPUE     + 1
      p%DOFs%NPYE     = p%DOFs%NPYE     + 1

       p%DOFs%PS     (  p%DOFs%NActvDOF) = DOF_Sw
       p%DOFs%PCE    (  p%DOFs%NPCE    ) = DOF_Sw
       p%DOFs%PDE    (  p%DOFs%NPDE    ) = DOF_Sw
       p%DOFs%PIE    (  p%DOFs%NPIE    ) = DOF_Sw
       p%DOFs%PTE    (  p%DOFs%NPTE    ) = DOF_Sw
       p%DOFs%PSE    (:,p%DOFs%NPSE (:)) = DOF_Sw
       p%DOFs%PUE    (  p%DOFs%NPUE    ) = DOF_Sw
       p%DOFs%PYE    (  p%DOFs%NPYE    ) = DOF_Sw

   ENDIF


   IF ( p%DOF_Flag(DOF_Hv  ) )  THEN  ! Platform heave.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPDE     = p%DOFs%NPDE     + 1
      p%DOFs%NPIE     = p%DOFs%NPIE     + 1
      p%DOFs%NPTE     = p%DOFs%NPTE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1
      p%DOFs%NPUE     = p%DOFs%NPUE     + 1
      p%DOFs%NPYE     = p%DOFs%NPYE     + 1

      p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_Hv
      p%DOFs%PCE     (  p%DOFs%NPCE    ) = DOF_Hv
      p%DOFs%PDE     (  p%DOFs%NPDE    ) = DOF_Hv
      p%DOFs%PIE     (  p%DOFs%NPIE    ) = DOF_Hv
      p%DOFs%PTE     (  p%DOFs%NPTE    ) = DOF_Hv
      p%DOFs%PSE     (:,p%DOFs%NPSE (:)) = DOF_Hv
      p%DOFs%PUE     (  p%DOFs%NPUE    ) = DOF_Hv
      p%DOFs%PYE     (  p%DOFs%NPYE    ) = DOF_Hv

   ENDIF


   IF ( p%DOF_Flag(DOF_R   ) )  THEN  ! Platform roll.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPDE     = p%DOFs%NPDE     + 1
      p%DOFs%NPIE     = p%DOFs%NPIE     + 1
      p%DOFs%NPTE     = p%DOFs%NPTE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1
      p%DOFs%NPUE     = p%DOFs%NPUE     + 1
      p%DOFs%NPYE     = p%DOFs%NPYE     + 1

      p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_R
      p%DOFs%PCE     (  p%DOFs%NPCE    ) = DOF_R
      p%DOFs%PDE     (  p%DOFs%NPDE    ) = DOF_R
      p%DOFs%PIE     (  p%DOFs%NPIE    ) = DOF_R
      p%DOFs%PTE     (  p%DOFs%NPTE    ) = DOF_R
      p%DOFs%PSE     (:,p%DOFs%NPSE (:)) = DOF_R
      p%DOFs%PUE     (  p%DOFs%NPUE    ) = DOF_R
      p%DOFs%PYE     (  p%DOFs%NPYE    ) = DOF_R

   ENDIF


   IF ( p%DOF_Flag(DOF_P   ) )  THEN  ! Platform pitch.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPDE     = p%DOFs%NPDE     + 1
      p%DOFs%NPIE     = p%DOFs%NPIE     + 1
      p%DOFs%NPTE     = p%DOFs%NPTE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1
      p%DOFs%NPUE     = p%DOFs%NPUE     + 1
      p%DOFs%NPYE     = p%DOFs%NPYE     + 1

      p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_P
      p%DOFs%PCE     (  p%DOFs%NPCE    ) = DOF_P
      p%DOFs%PDE     (  p%DOFs%NPDE    ) = DOF_P
      p%DOFs%PIE     (  p%DOFs%NPIE    ) = DOF_P
      p%DOFs%PTE     (  p%DOFs%NPTE    ) = DOF_P
      p%DOFs%PSE     (:,p%DOFs%NPSE (:)) = DOF_P
      p%DOFs%PUE     (  p%DOFs%NPUE    ) = DOF_P
      p%DOFs%PYE     (  p%DOFs%NPYE    ) = DOF_P

   ENDIF


   IF ( p%DOF_Flag(DOF_Y   ) )  THEN  ! Platform yaw.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPDE     = p%DOFs%NPDE     + 1
      p%DOFs%NPIE     = p%DOFs%NPIE     + 1
      p%DOFs%NPTE     = p%DOFs%NPTE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1
      p%DOFs%NPUE     = p%DOFs%NPUE     + 1
      p%DOFs%NPYE     = p%DOFs%NPYE     + 1

      p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_Y
      p%DOFs%PCE     (  p%DOFs%NPCE    ) = DOF_Y
      p%DOFs%PDE     (  p%DOFs%NPDE    ) = DOF_Y
      p%DOFs%PIE     (  p%DOFs%NPIE    ) = DOF_Y
      p%DOFs%PTE     (  p%DOFs%NPTE    ) = DOF_Y
      p%DOFs%PSE     (:,p%DOFs%NPSE (:)) = DOF_Y
      p%DOFs%PUE     (  p%DOFs%NPUE    ) = DOF_Y
      p%DOFs%PYE     (  p%DOFs%NPYE    ) = DOF_Y

   ENDIF


   IF ( p%DOF_Flag(DOF_TFA1) )  THEN  ! 1st tower fore-aft.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPDE     = p%DOFs%NPDE     + 1
      p%DOFs%NPIE     = p%DOFs%NPIE     + 1
      p%DOFs%NPTTE    = p%DOFs%NPTTE    + 1
      p%DOFs%NPTE     = p%DOFs%NPTE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1
      p%DOFs%NPUE     = p%DOFs%NPUE     + 1

      p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_TFA1
      p%DOFs%PCE     (  p%DOFs%NPCE    ) = DOF_TFA1
      p%DOFs%PDE     (  p%DOFs%NPDE    ) = DOF_TFA1
      p%DOFs%PIE     (  p%DOFs%NPIE    ) = DOF_TFA1
      p%DOFs%PTTE    (  p%DOFs%NPTTE   ) = DOF_TFA1
      p%DOFs%PTE     (  p%DOFs%NPTE    ) = DOF_TFA1
      p%DOFs%PSE     (:,p%DOFs%NPSE (:)) = DOF_TFA1
      p%DOFs%PUE     (  p%DOFs%NPUE    ) = DOF_TFA1

   ENDIF


   IF ( p%DOF_Flag(DOF_TSS1) )  THEN  ! 1st tower side-to-side.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPDE     = p%DOFs%NPDE     + 1
      p%DOFs%NPIE     = p%DOFs%NPIE     + 1
      p%DOFs%NPTTE    = p%DOFs%NPTTE    + 1
      p%DOFs%NPTE     = p%DOFs%NPTE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1
      p%DOFs%NPUE     = p%DOFs%NPUE     + 1

      p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_TSS1
      p%DOFs%PCE     (  p%DOFs%NPCE    ) = DOF_TSS1
      p%DOFs%PDE     (  p%DOFs%NPDE    ) = DOF_TSS1
      p%DOFs%PIE     (  p%DOFs%NPIE    ) = DOF_TSS1
      p%DOFs%PTTE    (  p%DOFs%NPTTE   ) = DOF_TSS1
      p%DOFs%PTE     (  p%DOFs%NPTE    ) = DOF_TSS1
      p%DOFs%PSE     (:,p%DOFs%NPSE (:)) = DOF_TSS1
      p%DOFs%PUE     (  p%DOFs%NPUE    ) = DOF_TSS1

   ENDIF


   IF ( p%DOF_Flag(DOF_TFA2) )  THEN  ! 2nd tower fore-aft.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPDE     = p%DOFs%NPDE     + 1
      p%DOFs%NPIE     = p%DOFs%NPIE     + 1
      p%DOFs%NPTTE    = p%DOFs%NPTTE    + 1
      p%DOFs%NPTE     = p%DOFs%NPTE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1
      p%DOFs%NPUE     = p%DOFs%NPUE     + 1

      p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_TFA2
      p%DOFs%PCE     (  p%DOFs%NPCE    ) = DOF_TFA2
      p%DOFs%PDE     (  p%DOFs%NPDE    ) = DOF_TFA2
      p%DOFs%PIE     (  p%DOFs%NPIE    ) = DOF_TFA2
      p%DOFs%PTTE    (  p%DOFs%NPTTE   ) = DOF_TFA2
      p%DOFs%PTE     (  p%DOFs%NPTE    ) = DOF_TFA2
      p%DOFs%PSE     (:,p%DOFs%NPSE (:)) = DOF_TFA2
      p%DOFs%PUE     (  p%DOFs%NPUE    ) = DOF_TFA2

   ENDIF


   IF ( p%DOF_Flag(DOF_TSS2) )  THEN  ! 2nd tower side-to-side.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPDE     = p%DOFs%NPDE     + 1
      p%DOFs%NPIE     = p%DOFs%NPIE     + 1
      p%DOFs%NPTTE    = p%DOFs%NPTTE    + 1
      p%DOFs%NPTE     = p%DOFs%NPTE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1
      p%DOFs%NPUE     = p%DOFs%NPUE     + 1

      p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_TSS2
      p%DOFs%PCE     (  p%DOFs%NPCE    ) = DOF_TSS2
      p%DOFs%PDE     (  p%DOFs%NPDE    ) = DOF_TSS2
      p%DOFs%PIE     (  p%DOFs%NPIE    ) = DOF_TSS2
      p%DOFs%PTTE    (  p%DOFs%NPTTE   ) = DOF_TSS2
      p%DOFs%PTE     (  p%DOFs%NPTE    ) = DOF_TSS2
      p%DOFs%PSE     (:,p%DOFs%NPSE (:)) = DOF_TSS2
      p%DOFs%PUE     (  p%DOFs%NPUE    ) = DOF_TSS2

   ENDIF


   IF ( p%DOF_Flag(DOF_Yaw ) )  THEN  ! Nacelle yaw.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPDE     = p%DOFs%NPDE     + 1
      p%DOFs%NPIE     = p%DOFs%NPIE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1
      p%DOFs%NPUE     = p%DOFs%NPUE     + 1

      p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_Yaw
      p%DOFs%PCE     (  p%DOFs%NPCE    ) = DOF_Yaw
      p%DOFs%PDE     (  p%DOFs%NPDE    ) = DOF_Yaw
      p%DOFs%PIE     (  p%DOFs%NPIE    ) = DOF_Yaw
      p%DOFs%PSE     (:,p%DOFs%NPSE (:)) = DOF_Yaw
      p%DOFs%PUE     (  p%DOFs%NPUE    ) = DOF_Yaw

   ENDIF


   IF ( p%DOF_Flag(DOF_TFrl) )  THEN  ! Tail-furl.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPIE     = p%DOFs%NPIE     + 1

      p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_TFrl
      p%DOFs%PIE     (  p%DOFs%NPIE    ) = DOF_TFrl

   ENDIF


   IF ( p%DOF_Flag(DOF_RFrl) )  THEN  ! Rotor-furl.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPDE     = p%DOFs%NPDE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1

      p%DOFs%PS     (  p%DOFs%NActvDOF) = DOF_RFrl
      p%DOFs%PCE    (  p%DOFs%NPCE    ) = DOF_RFrl
      p%DOFs%PDE    (  p%DOFs%NPDE    ) = DOF_RFrl
      p%DOFs%PSE    (:,p%DOFs%NPSE (:)) = DOF_RFrl

   ENDIF


   IF ( p%DOF_Flag(DOF_GeAz) )  THEN  ! Generator azimuth.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1

      p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_GeAz
      p%DOFs%PCE     (  p%DOFs%NPCE    ) = DOF_GeAz
      p%DOFs%PSE     (:,p%DOFs%NPSE (:)) = DOF_GeAz

   ENDIF


   IF ( p%DOF_Flag(DOF_DrTr) )  THEN  ! Drivetrain torsion.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1

      p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_DrTr
      p%DOFs%PCE     (  p%DOFs%NPCE    ) = DOF_DrTr
      p%DOFs%PSE     (:,p%DOFs%NPSE (:)) = DOF_DrTr

   ENDIF


   IF ( p%NumBl == 2 )  THEN
      IF ( p%DOF_Flag(DOF_Teet   ) )  THEN  ! Rotor-teeter.

         p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
         p%DOFs%NPCE     = p%DOFs%NPCE     + 1
         p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1

         p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_Teet
         p%DOFs%PCE     (  p%DOFs%NPCE    ) = DOF_Teet
         p%DOFs%PSE     (:,p%DOFs%NPSE (:)) = DOF_Teet

      ENDIF
   ENDIF


   DO K = 1,p%NumBl ! Loop through all blades
      IF ( p%DOF_Flag(DOF_BF(K,1)) )  THEN  ! 1st blade flap.

         p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
         p%DOFs%NPSBE(K) = p%DOFs%NPSBE(K) + 1
         p%DOFs%NPSE (K) = p%DOFs%NPSE (K) + 1

         p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_BF(K,1)
         p%DOFs%PSBE    (K,p%DOFs%NPSBE(K)) = DOF_BF(K,1)
         p%DOFs%PSE     (K,p%DOFs%NPSE (K)) = DOF_BF(K,1)

      ENDIF
   ENDDO          ! K - Blades


   DO K = 1,p%NumBl ! Loop through all blades
      IF ( p%DOF_Flag(DOF_BE(K,1)) )  THEN  ! 1st blade edge.

         p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
         p%DOFs%NPSBE(K) = p%DOFs%NPSBE(K) + 1
         p%DOFs%NPSE (K) = p%DOFs%NPSE (K) + 1

         p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_BE(K,1)
         p%DOFs%PSBE    (K,p%DOFs%NPSBE(K)) = DOF_BE(K,1)
         p%DOFs%PSE     (K,p%DOFs%NPSE (K)) = DOF_BE(K,1)

      ENDIF
   ENDDO          ! K - Blades


   DO K = 1,p%NumBl ! Loop through all blades
      IF ( p%DOF_Flag(DOF_BF(K,2)) )  THEN  ! 2nd blade flap.

         p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
         p%DOFs%NPSBE(K) = p%DOFs%NPSBE(K) + 1
         p%DOFs%NPSE (K) = p%DOFs%NPSE (K) + 1

         p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_BF(K,2)
         p%DOFs%PSBE    (K,p%DOFs%NPSBE(K)) = DOF_BF(K,2)
         p%DOFs%PSE     (K,p%DOFs%NPSE (K)) = DOF_BF(K,2)

      ENDIF
   ENDDO          ! K - Blades



      ! Compute the sorted (from smallest to largest p%DOFs index) version of PS(),
      !   SrtPS(), and SrtPSNAUG().  At the same time compute Diag(), which is an
      !   array containing the indices of SrtPS() associated with each enabled
      !   DOF; that is, SrtPS(Diag(I)) = I:
      ! NOTE: This calculation is recomputing NActvDOF as computed above.  This is
      !       of no concern however, since the resulting value will be the same.

   p%DOFs%NActvDOF = 0
   DO I = 1,p%NDOF  ! Loop through all DOFs
      IF ( p%DOF_Flag(I) )  THEN   ! .TRUE. if the corresponding DOF is enabled

         p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1

         p%DOFs%SrtPS     (p%DOFs%NActvDOF) = I
         p%DOFs%SrtPSNAUG (p%DOFs%NActvDOF) = I
         p%DOFs%Diag      (I           ) = p%DOFs%NActvDOF

      ENDIF
   ENDDO          ! I - All DOFs

   p%DOFs%SrtPSNAUG ( p%DOFs%NActvDOF + 1 ) = p%NAug


   RETURN
END SUBROUTINE SetEnabledDOFIndexArrays
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is used to define the internal coordinate systems for this particular time step.
!! It also sets the TeeterAng and TeetAngVel for this time step.
SUBROUTINE SetCoordSy( t, CoordSys, RtHSdat, BlPitch, p, x, ErrStat, ErrMsg )

      ! Subroutine arguments (passed variables)

   REAL(DbKi),                   INTENT(IN)    :: t                             !< Current simulation time, in seconds (used only for SmllRotTrans error messages)
   REAL(ReKi),                   INTENT(IN)    :: BlPitch (:)                   !< The current blade pitch
   TYPE(ED_CoordSys),            INTENT(INOUT) :: CoordSys                      !< The coordinate systems to be set
   TYPE(ED_RtHndSide),           INTENT(INOUT) :: RtHSdat                       !< data from the RtHndSid module
   TYPE(ED_ParameterType),       INTENT(IN)    :: p                             !< The module's parameters
   TYPE(ED_ContinuousStateType), INTENT(IN)    :: x                             !< The module's continuous states

   INTEGER(IntKi),               INTENT(OUT)    :: ErrStat                      !< Error status
   CHARACTER(*),                 INTENT(OUT)    :: ErrMsg                       !< Error message

      ! Local variables

   REAL(R8Ki)                   :: CAzimuth                                        ! COS( rotor azimuth angle ).
   REAL(R8Ki)                   :: CgRotAng                                        ! COS( gRotAng ).
   REAL(R8Ki)                   :: CNacYaw                                         ! COS( nacelle yaw angle ).
   REAL(R8Ki)                   :: CosPitch                                        ! COS( the current pitch angle ).
   REAL(R8Ki)                   :: CPitPTwstA                                      ! COS( BlPitch(K) + AeroTwst(J) ) found using the sum of angles formulae of cosine.
   REAL(R8Ki)                   :: CPitPTwstS                                      ! COS( BlPitch(K) + ThetaS(K,J) ) found using the sum of angles formulae of cosine.
   REAL(R8Ki)                   :: CRotFurl                                        ! COS( rotor-furl angle ).
   REAL(R8Ki)                   :: CTailFurl                                       ! COS( tail-furl angle ).
   REAL(R8Ki)                   :: CTeetAng                                        ! COS( TeetAng ).
   REAL(R8Ki)                   :: g1Prime   (3)                                   ! = g1.
   REAL(R8Ki)                   :: g2Prime   (3)                                   ! completes the right-handed gPrime-vector triad
   REAL(R8Ki)                   :: g3Prime   (3)                                   ! = g3 rotated about g1 so that parallel to the pitching axis of blade K (i.e., the current blade in the blade loop).
   REAL(R8Ki)                   :: gRotAng                                         ! Angle of rotation about g1 to get from the g to the gPrime system.
   REAL(R8Ki)                   :: Lj1       (3)                                   ! vector / direction Lj1 at node J for blade K.
   REAL(R8Ki)                   :: Lj2       (3)                                   ! vector / direction Lj2 at node J for blade K.
   REAL(R8Ki)                   :: Lj3       (3)                                   ! vector / direction Lj3 at node J for blade K.
   REAL(R8Ki)                   :: SAzimuth                                        ! SIN( rotor azimuth angle ).
   REAL(R8Ki)                   :: SgRotAng                                        ! SIN( gRotAng ).
   REAL(R8Ki)                   :: SinPitch                                        ! SIN( the current pitch angle ).
   REAL(R8Ki)                   :: SNacYaw                                         ! SIN( nacelle yaw angle ).
   REAL(R8Ki)                   :: SPitPTwstA                                      ! SIN( BlPitch(K) + AeroTwst(J) ) found using the sum of angles formulae of sine.
   REAL(R8Ki)                   :: SPitPTwstS                                      ! SIN( BlPitch(K) + ThetaS(K,J) ) found using the sum of angles formulae of sine.
   REAL(R8Ki)                   :: SRotFurl                                        ! SIN( rotor-furl angle ).
   REAL(R8Ki)                   :: STailFurl                                       ! SIN( tail-furl angle ).
   REAL(R8Ki)                   :: STeetAng                                        ! SIN( TeetAng ).
   REAL(R8Ki)                   :: ThetaFA                                         ! Tower fore-aft tilt deflection angle.
   REAL(R8Ki)                   :: ThetaIP                                         ! Blade in-plane deflection angle at node J for blade K.
   REAL(R8Ki)                   :: ThetaLxb                                        ! Blade deflection angle about the Lxb (n1) -axis at node J for blade K.
   REAL(R8Ki)                   :: ThetaLyb                                        ! Blade deflection angle about the Lyb (n2) -axis at node J for blade K.
   REAL(R8Ki)                   :: ThetaOoP                                        ! Blade out-of-plane deflection angle at node J for blade K.
   REAL(R8Ki)                   :: ThetaSS                                         ! Tower side-to-side tilt deflection angle.
   REAL(R8Ki)                   :: TransMat  (3,3)                                 ! The resulting transformation matrix due to three orthogonal rotations, (-).

   INTEGER(IntKi)               :: J                                               ! Loops through nodes / elements.
   INTEGER(IntKi)               :: K                                               ! Loops through blades.


   INTEGER(IntKi)               :: ErrStat2                      ! Temporary error status
   CHARACTER(ErrMsgLen)         :: ErrMsg2                       ! Temporary error message


   ErrStat = ErrID_None
   ErrMsg  = ''

      ! Inertial frame coordinate system:

   CoordSys%z1 = (/ 1.0_R8Ki, 0.0_R8Ki, 0.0_R8Ki /)   ! Vector / direction z1 (=  xi from the IEC coord. system).
   CoordSys%z2 = (/ 0.0_R8Ki, 1.0_R8Ki, 0.0_R8Ki /)   ! Vector / direction z2 (=  zi from the IEC coord. system).
   CoordSys%z3 = (/ 0.0_R8Ki, 0.0_R8Ki, 1.0_R8Ki /)   ! Vector / direction z3 (= -yi from the IEC coord. system).


      ! Tower base / platform coordinate system:

   CALL SmllRotTrans( 'platform displacement (ElastoDyn SetCoordSy)', x%QT(DOF_R), x%QT(DOF_Y), -x%QT(DOF_P), TransMat, TRIM(Num2LStr(t))//' s', ErrStat2, ErrMsg2 )  ! Get the transformation matrix, TransMat, from inertial frame to tower base / platform coordinate systems.
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN

   CoordSys%a1 = TransMat(1,1)*CoordSys%z1 + TransMat(1,2)*CoordSys%z2 + TransMat(1,3)*CoordSys%z3 ! Vector / direction a1 (=  xt from the IEC coord. system).
   CoordSys%a2 = TransMat(2,1)*CoordSys%z1 + TransMat(2,2)*CoordSys%z2 + TransMat(2,3)*CoordSys%z3 ! Vector / direction a2 (=  zt from the IEC coord. system).
   CoordSys%a3 = TransMat(3,1)*CoordSys%z1 + TransMat(3,2)*CoordSys%z2 + TransMat(3,3)*CoordSys%z3 ! Vector / direction a3 (= -yt from the IEC coord. system).


   DO J = 1,p%TwrNodes ! Loop through the tower nodes / elements


      ! Tower element-fixed coordinate system:

      ThetaFA = -p%TwrFASF(1,J       ,1)*x%QT(DOF_TFA1) - p%TwrFASF(2,J       ,1)*x%QT(DOF_TFA2)
      ThetaSS =  p%TwrSSSF(1,J       ,1)*x%QT(DOF_TSS1) + p%TwrSSSF(2,J       ,1)*x%QT(DOF_TSS2)

      CALL SmllRotTrans( 'tower deflection (ElastoDyn SetCoordSy)', ThetaSS, 0.0_R8Ki, ThetaFA, TransMat, TRIM(Num2LStr(t))//' s', ErrStat2, ErrMsg2 )   ! Get the transformation matrix, TransMat, from tower-base to tower element-fixed coordinate systems.
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF (ErrStat >= AbortErrLev) RETURN

      CoordSys%t1(J,:) = TransMat(1,1)*CoordSys%a1 + TransMat(1,2)*CoordSys%a2 + TransMat(1,3)*CoordSys%a3  ! Vector / direction t1 for tower node J (=  Lxt from the IEC coord. system).
      CoordSys%t2(J,:) = TransMat(2,1)*CoordSys%a1 + TransMat(2,2)*CoordSys%a2 + TransMat(2,3)*CoordSys%a3  ! Vector / direction t2 for tower node J (=  Lzt from the IEC coord. system).
      CoordSys%t3(J,:) = TransMat(3,1)*CoordSys%a1 + TransMat(3,2)*CoordSys%a2 + TransMat(3,3)*CoordSys%a3  ! Vector / direction t3 for tower node J (= -Lyt from the IEC coord. system).


   ENDDO ! J - Tower nodes / elements


      ! Tower-top / base plate coordinate system:

   ThetaFA    = -p%TwrFASF(1,p%TTopNode,1)*x%QT(DOF_TFA1) - p%TwrFASF(2,p%TTopNode,1)*x%QT(DOF_TFA2)
   ThetaSS    =  p%TwrSSSF(1,p%TTopNode,1)*x%QT(DOF_TSS1) + p%TwrSSSF(2,p%TTopNode,1)*x%QT(DOF_TSS2)

   CALL SmllRotTrans( 'tower deflection (ElastoDyn SetCoordSy)', ThetaSS, 0.0_R8Ki, ThetaFA, TransMat, TRIM(Num2LStr(t))//' s', ErrStat2, ErrMsg2 )   ! Get the transformation matrix, TransMat, from tower-base to tower-top/base-plate coordinate systems.
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN

   CoordSys%b1 = TransMat(1,1)*CoordSys%a1 + TransMat(1,2)*CoordSys%a2 + TransMat(1,3)*CoordSys%a3 ! Vector / direction b1 (=  xp from the IEC coord. system).
   CoordSys%b2 = TransMat(2,1)*CoordSys%a1 + TransMat(2,2)*CoordSys%a2 + TransMat(2,3)*CoordSys%a3 ! Vector / direction b2 (=  zp from the IEC coord. system).
   CoordSys%b3 = TransMat(3,1)*CoordSys%a1 + TransMat(3,2)*CoordSys%a2 + TransMat(3,3)*CoordSys%a3 ! Vector / direction b3 (= -yp from the IEC coord. system).


      ! Nacelle / yaw coordinate system:

   CNacYaw  = COS( x%QT(DOF_Yaw ) )
   SNacYaw  = SIN( x%QT(DOF_Yaw ) )

   CoordSys%d1 = CNacYaw*CoordSys%b1 - SNacYaw*CoordSys%b3     ! Vector / direction d1 (=  xn from the IEC coord. system).
   CoordSys%d2 = CoordSys%b2                                   ! Vector / direction d2 (=  zn from the IEC coord. system).
   CoordSys%d3 = SNacYaw*CoordSys%b1 + CNacYaw*CoordSys%b3     ! Vector / direction d3 (= -yn from the IEC coord. system).


      ! Rotor-furl coordinate system:

   CRotFurl = COS( x%QT(DOF_RFrl) )
   SRotFurl = SIN( x%QT(DOF_RFrl) )

   CoordSys%rf1 = ( (   1.0 - p%CRFrlSkw2*p%CRFrlTlt2 )*CRotFurl   + p%CRFrlSkw2*p%CRFrlTlt2          )*CoordSys%d1 &
                + ( p%CRFrlSkew*p%CSRFrlTlt*( 1.0 -     CRotFurl ) - p%SRFrlSkew*p%CRFrlTilt*SRotFurl )*CoordSys%d2 &
                + ( p%CSRFrlSkw*p%CRFrlTlt2*( CRotFurl - 1.0     ) -             p%SRFrlTilt*SRotFurl )*CoordSys%d3
   CoordSys%rf2 = ( p%CRFrlSkew*p%CSRFrlTlt*( 1.0 -     CRotFurl ) + p%SRFrlSkew*p%CRFrlTilt*SRotFurl )*CoordSys%d1 &
                + (             p%CRFrlTlt2*            CRotFurl   +             p%SRFrlTlt2          )*CoordSys%d2 &
                + ( p%SRFrlSkew*p%CSRFrlTlt*( CRotFurl - 1.0     ) + p%CRFrlSkew*p%CRFrlTilt*SRotFurl )*CoordSys%d3
   CoordSys%rf3 = ( p%CSRFrlSkw*p%CRFrlTlt2*( CRotFurl - 1.0     ) +             p%SRFrlTilt*SRotFurl )*CoordSys%d1 &
                + ( p%SRFrlSkew*p%CSRFrlTlt*( CRotFurl - 1.0     ) - p%CRFrlSkew*p%CRFrlTilt*SRotFurl )*CoordSys%d2 &
                + ( (   1.0 - p%SRFrlSkw2*p%CRFrlTlt2 )*CRotFurl   + p%SRFrlSkw2*p%CRFrlTlt2          )*CoordSys%d3
   CoordSys%rfa = p%CRFrlSkew*p%CRFrlTilt*CoordSys%d1 + p%SRFrlTilt*CoordSys%d2 - p%SRFrlSkew*p%CRFrlTilt*CoordSys%d3


      ! Shaft coordinate system:

   CoordSys%c1 =  p%CShftSkew*p%CShftTilt*CoordSys%rf1 + p%SShftTilt*CoordSys%rf2 - p%SShftSkew*p%CShftTilt*CoordSys%rf3  ! Vector / direction c1 (=  xs from the IEC coord. system).
   CoordSys%c2 = -p%CShftSkew*p%SShftTilt*CoordSys%rf1 + p%CShftTilt*CoordSys%rf2 + p%SShftSkew*p%SShftTilt*CoordSys%rf3  ! Vector / direction c2 (=  zs from the IEC coord. system).
   CoordSys%c3 =  p%SShftSkew*            CoordSys%rf1                            + p%CShftSkew*            CoordSys%rf3  ! Vector / direction c3 (= -ys from the IEC coord. system).


      ! Azimuth coordinate system:

   CAzimuth = COS( x%QT(DOF_DrTr) + x%QT(DOF_GeAz) )
   SAzimuth = SIN( x%QT(DOF_DrTr) + x%QT(DOF_GeAz) )

   CoordSys%e1 =  CoordSys%c1                                  ! Vector / direction e1 (=  xa from the IEC coord. system).
   CoordSys%e2 =  CAzimuth*CoordSys%c2 + SAzimuth*CoordSys%c3  ! Vector / direction e2 (=  ya from the IEC coord. system).
   CoordSys%e3 = -SAzimuth*CoordSys%c2 + CAzimuth*CoordSys%c3  ! Vector / direction e3 (=  za from the IEC coord. system).


      ! Teeter coordinate system:

      ! Lets define TeetAng, which is the current teeter angle (= QT(DOF_Teet) for
      !   2-blader or 0 for 3-blader) and is used in place of QT(DOF_Teet)
      !   throughout SUBROUTINE RtHS().  Doing it this way, we can run the same
      !   equations of motion for both the 2 and 3-blader configurations even
      !   though a 3-blader does not have a teetering DOF.

   IF ( p%NumBl == 2 )  THEN ! 2-blader
      RtHSdat%TeetAng    = x%QT (DOF_Teet)
      RtHSdat%TeetAngVel = x%QDT(DOF_Teet)
   ELSE                    ! 3-blader
      RtHSdat%TeetAng    = 0.0  ! Teeter is not an available DOF for a 3-blader
      RtHSdat%TeetAngVel = 0.0  ! Teeter is not an available DOF for a 3-blader
   ENDIF
   CTeetAng = COS( RtHSdat%TeetAng )
   STeetAng = SIN( RtHSdat%TeetAng )

   CoordSys%f1 = CTeetAng*CoordSys%e1 - STeetAng*CoordSys%e3       ! Vector / direction f1.
   CoordSys%f2 = CoordSys%e2                                       ! Vector / direction f2.
   CoordSys%f3 = STeetAng*CoordSys%e1 + CTeetAng*CoordSys%e3       ! Vector / direction f3.


      ! Hub / delta-3 coordinate system:

   CoordSys%g1 =  CoordSys%f1                                      ! Vector / direction g1 (=  xh from the IEC coord. system).
   CoordSys%g2 =  p%CosDel3*CoordSys%f2 + p%SinDel3*CoordSys%f3    ! Vector / direction g2 (=  yh from the IEC coord. system).
   CoordSys%g3 = -p%SinDel3*CoordSys%f2 + p%CosDel3*CoordSys%f3    ! Vector / direction g3 (=  zh from the IEC coord. system).


   DO K = 1,p%NumBl ! Loop through all blades


      ! Hub (Prime) coordinate system rotated to match blade K.

       gRotAng = p%TwoPiNB*(K-1)
      CgRotAng = COS( gRotAng )
      SgRotAng = SIN( gRotAng )

      g1Prime =  CoordSys%g1
      g2Prime =  CgRotAng*CoordSys%g2 + SgRotAng*CoordSys%g3
      g3Prime = -SgRotAng*CoordSys%g2 + CgRotAng*CoordSys%g3


      ! Coned coordinate system:

      CoordSys%i1(K,:) = p%CosPreC(K)*g1Prime - p%SinPreC(K)*g3Prime  ! i1(K,:) = vector / direction i1 for blade K (=  xcK from the IEC coord. system).
      CoordSys%i2(K,:) = g2Prime                                      ! i2(K,:) = vector / direction i2 for blade K (=  ycK from the IEC coord. system).
      CoordSys%i3(K,:) = p%SinPreC(K)*g1Prime + p%CosPreC(K)*g3Prime  ! i3(K,:) = vector / direction i3 for blade K (=  zcK from the IEC coord. system).


      ! Blade / pitched coordinate system:

      CosPitch = COS( REAL(BlPitch(K),R8Ki) )
      SinPitch = SIN( REAL(BlPitch(K),R8Ki) )

      CoordSys%j1(K,:) = CosPitch*CoordSys%i1(K,:) - SinPitch*CoordSys%i2(K,:)      ! j1(K,:) = vector / direction j1 for blade K (=  xbK from the IEC coord. system).
      CoordSys%j2(K,:) = SinPitch*CoordSys%i1(K,:) + CosPitch*CoordSys%i2(K,:)      ! j2(K,:) = vector / direction j2 for blade K (=  ybK from the IEC coord. system).
      CoordSys%j3(K,:) = CoordSys%i3(K,:)                                           ! j3(K,:) = vector / direction j3 for blade K (=  zbK from the IEC coord. system).


      DO J = 0,p%TipNode ! Loop through the blade nodes / elements


      ! Blade coordinate system aligned with local structural axes (not element fixed):

         Lj1 = p%CThetaS(K,J)*CoordSys%j1(K,:) - p%SThetaS(K,J)*CoordSys%j2(K,:)  ! vector / direction Lj1 at node J for blade K
         Lj2 = p%SThetaS(K,J)*CoordSys%j1(K,:) + p%CThetaS(K,J)*CoordSys%j2(K,:)  ! vector / direction Lj2 at node J for blade K
         Lj3 = CoordSys%j3(K,:)                                                   ! vector / direction Lj3 at node J for blade K


      ! Blade element-fixed coordinate system aligned with local structural axes:

         ThetaOoP =   p%TwistedSF(K,1,1,J,1)*x%QT( DOF_BF(K,1) ) &
                    + p%TwistedSF(K,1,2,J,1)*x%QT( DOF_BF(K,2) ) &
                    + p%TwistedSF(K,1,3,J,1)*x%QT( DOF_BE(K,1) )
         ThetaIP  = - p%TwistedSF(K,2,1,J,1)*x%QT( DOF_BF(K,1) ) &
                    - p%TwistedSF(K,2,2,J,1)*x%QT( DOF_BF(K,2) ) &
                    - p%TwistedSF(K,2,3,J,1)*x%QT( DOF_BE(K,1) )

         ThetaLxb = p%CThetaS(K,J)*ThetaIP - p%SThetaS(K,J)*ThetaOoP
         ThetaLyb = p%SThetaS(K,J)*ThetaIP + p%CThetaS(K,J)*ThetaOoP

         CALL SmllRotTrans( 'blade deflection (ElastoDyn SetCoordSy)', ThetaLxb, ThetaLyb, 0.0_R8Ki, TransMat, TRIM(Num2LStr(t))//' s', ErrStat2, ErrMsg2 ) ! Get the transformation matrix, TransMat, from blade coordinate system aligned with local structural axes (not element fixed) to blade element-fixed coordinate system aligned with local structural axes.
            CALL CheckError( ErrStat2, ErrMsg2 )
            IF (ErrStat >= AbortErrLev) RETURN

         CoordSys%n1(K,J,:) = TransMat(1,1)*Lj1 + TransMat(1,2)*Lj2 + TransMat(1,3)*Lj3   ! Vector / direction n1 for node J of blade K (= LxbK from the IEC coord. system).
         CoordSys%n2(K,J,:) = TransMat(2,1)*Lj1 + TransMat(2,2)*Lj2 + TransMat(2,3)*Lj3   ! Vector / direction n2 for node J of blade K (= LybK from the IEC coord. system).
         CoordSys%n3(K,J,:) = TransMat(3,1)*Lj1 + TransMat(3,2)*Lj2 + TransMat(3,3)*Lj3   ! Vector / direction n3 for node J of blade K (= LzbK from the IEC coord. system).

      ! skip these next CoordSys variables at the root and the tip; they are required only for AD14:
         
         if (j == 0 .or. j==p%TipNode) cycle  
      
         
      ! Blade element-fixed coordinate system used for calculating and returning
      !    aerodynamics loads:
      ! This coordinate system is rotated about positive n3 by the angle
      !    BlPitch(K) + ThetaS(K,J) and is coincident with the i-vector triad
      !    when the blade is undeflected.

         CPitPTwstS = CosPitch*p%CThetaS(K,J) - SinPitch*p%SThetaS(K,J)  ! = COS( BlPitch(K) + ThetaS(K,J) ) found using the sum of angles formulae of cosine.
         SPitPTwstS = CosPitch*p%SThetaS(K,J) + SinPitch*p%CThetaS(K,J)  ! = SIN( BlPitch(K) + ThetaS(K,J) ) found using the sum of angles formulae of   sine.

         CoordSys%m1(K,J,:)  =  CPitPTwstS*CoordSys%n1(K,J,:) + SPitPTwstS*CoordSys%n2(K,J,:)   ! m1(K,J,:) = vector / direction m1 for node J of blade K (used to calc. and return aerodynamic loads from AeroDyn).
         CoordSys%m2(K,J,:)  = -SPitPTwstS*CoordSys%n1(K,J,:) + CPitPTwstS*CoordSys%n2(K,J,:)   ! m2(K,J,:) = vector / direction m2 for node J of blade K (used to calc. and return aerodynamic loads from AeroDyn).
         CoordSys%m3(K,J,:)  =  CoordSys%n3(K,J,:)                                              ! m3(K,J,:) = vector / direction m3 for node J of blade K (used to calc. and return aerodynamic loads from AeroDyn).


      ! Calculate the trailing edge coordinate system used in noise calculations.
      ! This coordinate system is blade element-fixed and oriented with the local
      !   aerodynamic axes (te2 points toward trailing edge, te1 points toward
      !   suction surface):

         CPitPTwstA = CosPitch*p%CAeroTwst(J) - SinPitch*p%SAeroTwst(J)  ! = COS( BlPitch(K) + AeroTwst(J) ) found using the sum of angles formulae of cosine.
         SPitPTwstA = CosPitch*p%SAeroTwst(J) + SinPitch*p%CAeroTwst(J)  ! = SIN( BlPitch(K) + AeroTwst(J) ) found using the sum of angles formulae of   sine.

         CoordSys%te1(K,J,:) =  CPitPTwstA*CoordSys%m1(K,J,:) - SPitPTwstA*CoordSys%m2(K,J,:)   ! te1(K,J,:) = vector / direction te1 for node J of blade K (used to calc. noise and to calc. and return aerodynamic loads from AeroDyn).
         CoordSys%te2(K,J,:) =  SPitPTwstA*CoordSys%m1(K,J,:) + CPitPTwstA*CoordSys%m2(K,J,:)   ! te2(K,J,:) = vector / direction te2 for node J of blade K (used to calc. noise and to calc. and return aerodynamic loads from AeroDyn).
         CoordSys%te3(K,J,:) =  CoordSys%m3(K,J,:)                                              ! te3(K,J,:) = vector / direction te3 for node J of blade K (used to calc. noise and to calc. and return aerodynamic loads from AeroDyn).


      ENDDO ! J - Blade nodes / elements


   ENDDO ! K - Blades


      ! Tail-furl coordinate system:

   CTailFurl = COS( x%QT(DOF_TFrl) )
   STailFurl = SIN( x%QT(DOF_TFrl) )

   CoordSys%tf1 = ( ( 1.0 - p%CTFrlSkw2*p%CTFrlTlt2 )*CTailFurl  + p%CTFrlSkw2*p%CTFrlTlt2           )*CoordSys%d1 &
                + ( p%CTFrlSkew*p%CSTFrlTlt*(  1.0 - CTailFurl ) - p%STFrlSkew*p%CTFrlTilt*STailFurl )*CoordSys%d2 &
                + ( p%CSTFrlSkw*p%CTFrlTlt2*( CTailFurl - 1.0  ) -             p%STFrlTilt*STailFurl )*CoordSys%d3
   CoordSys%tf2 = ( p%CTFrlSkew*p%CSTFrlTlt*(  1.0 - CTailFurl ) + p%STFrlSkew*p%CTFrlTilt*STailFurl )*CoordSys%d1 &
                + (             p%CTFrlTlt2*         CTailFurl +               p%STFrlTlt2           )*CoordSys%d2 &
                + ( p%STFrlSkew*p%CSTFrlTlt*( CTailFurl - 1.0  ) + p%CTFrlSkew*p%CTFrlTilt*STailFurl )*CoordSys%d3
   CoordSys%tf3 = ( p%CSTFrlSkw*p%CTFrlTlt2*( CTailFurl - 1.0  ) +             p%STFrlTilt*STailFurl )*CoordSys%d1 &
                + ( p%STFrlSkew*p%CSTFrlTlt*( CTailFurl - 1.0  ) - p%CTFrlSkew*p%CTFrlTilt*STailFurl )*CoordSys%d2 &
                + ( ( 1.0 - p%STFrlSkw2*p%CTFrlTlt2 )*CTailFurl  + p%STFrlSkw2*p%CTFrlTlt2           )*CoordSys%d3
   CoordSys%tfa = p%CTFrlSkew*p%CTFrlTilt*CoordSys%d1 + p%STFrlTilt*CoordSys%d2 - p%STFrlSkew*p%CTFrlTilt*CoordSys%d3


   RETURN
CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'SetCoordSy:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
         END IF

      END IF


   END SUBROUTINE CheckError
!----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE SetCoordSy
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine computes the rotor-furl moment due to rotor-furl deflection and rate.
SUBROUTINE RFurling( t, p, RFrlDef, RFrlRate, RFrlMom )
      ! Passed Variables:
   REAL(DbKi), INTENT(IN)              :: t                                   !< simulation time
   TYPE(ED_ParameterType), INTENT(IN)  :: p                                   !< parameters from the structural dynamics module
   REAL(R8Ki), INTENT(IN )             :: RFrlDef                             !< The rotor-furl deflection, x%QT(DOF_RFrl)
   REAL(ReKi), INTENT(OUT)             :: RFrlMom                             !< The total moment supplied by the springs, and dampers
   REAL(R8Ki), INTENT(IN )             :: RFrlRate                            !< The rotor-furl rate, x%QDT(DOF_RFrl)
      ! Local variables:
   REAL(ReKi)                   :: RFrlDMom                                   ! The moment supplied by the rotor-furl dampers
   REAL(ReKi)                   :: RFrlSMom                                   ! The moment supplied by the rotor-furl springs

   SELECT CASE ( p%RFrlMod ) ! Which rotor-furl model are we using?

      CASE ( 0_IntKi )       ! None!

         RFrlMom = 0.0

      CASE ( 1_IntKi )        ! Standard (using inputs from the FAST furling input file).

         ! Linear spring:
         RFrlSMom = -p%RFrlSpr*RFrlDef

         ! Add spring-stops:
         IF ( RFrlDef > p%RFrlUSSP )  THEN       ! Up-stop
            RFrlSMom = RFrlSMom - p%RFrlUSSpr*( RFrlDef - p%RFrlUSSP )
         ELSEIF ( RFrlDef < p%RFrlDSSP )  THEN   ! Down-stop
            RFrlSMom = RFrlSMom - p%RFrlDSSpr*( RFrlDef - p%RFrlDSSP )
         ENDIF

         ! Linear damper:
         RFrlDMom = -p%RFrlDmp*RFrlRate

         ! Add damper-stops:
         IF ( RFrlDef > p%RFrlUSDP )  THEN       ! Up-stop
            RFrlDMom = RFrlDMom - p%RFrlUSDmp*RFrlRate
         ELSEIF ( RFrlDef < p%RFrlDSDP )  THEN   ! Down-stop
            RFrlDMom = RFrlDMom - p%RFrlDSDmp*RFrlRate
         ENDIF

         ! Total up all the moments.
         RFrlMom = RFrlSMom + RFrlDMom

      CASE ( 2_IntKi )              ! User-defined rotor-furl spring/damper model.

         CALL UserRFrl ( RFrlDef, RFrlRate, t, p%RootName, RFrlMom )

   END   SELECT
END SUBROUTINE RFurling
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine computes the teeter moment due to teeter deflection and rate.
SUBROUTINE Teeter( t, p, TeetDef, TeetRate, TeetMom )
!..................................................................................................................................

      ! Passed Variables:
   REAL(DbKi), INTENT(IN)             :: t                                       !< simulation time
   TYPE(ED_ParameterType), INTENT(IN) :: p                                       !< parameters from the structural dynamics module
   REAL(R8Ki), INTENT(IN )            :: TeetDef                                 !< The teeter deflection, x%QT(DOF_Teet).
   REAL(ReKi), INTENT(OUT)            :: TeetMom                                 !< The total moment supplied by the stop, spring, and damper.
   REAL(R8Ki), INTENT(IN )            :: TeetRate                                !< The teeter rate, x%QDT(DOF_Teet).


      ! Local variables:
   REAL(ReKi)                         :: AbsDef                                   ! Absolute value of the teeter deflection.
   REAL(ReKi)                         :: SprgDef                                  ! Deflection past the spring.
   REAL(ReKi)                         :: StopDef                                  ! Deflection past the stop.
   REAL(ReKi)                         :: TeetDMom                                 ! The moment supplied by the damper.
   REAL(ReKi)                         :: TeetFMom                                 ! The moment supplied by Coulomb-friction damping.
   REAL(ReKi)                         :: TeetKMom                                 ! The moment supplied by the spring.
   REAL(ReKi)                         :: TeetSMom                                 ! The moment supplied by the stop.



   SELECT CASE ( p%TeetMod ) ! Which teeter model are we using?

   CASE ( 0_IntKi )              ! None!


      TeetMom = 0.0_ReKi


   CASE ( 1_IntKi )              ! Standard (using inputs from the primary FAST input file).


      ! Compute the absulute value of the deflection.

      AbsDef  = ABS( TeetDef )


      ! Linear teeter spring.

      SprgDef = AbsDef - p%TeetSStP

      IF ( SprgDef > 0.0_ReKi )  THEN
         TeetKMom = -SIGN( SprgDef*p%TeetSSSp, real(TeetDef,ReKi) )
      ELSE
         TeetKMom = 0.0_ReKi
      ENDIF


      ! Compute teeter-stop moment if hard stop has been contacted.

      StopDef = AbsDef - p%TeetHStP

      IF ( StopDef > 0.0_ReKi )  THEN
         TeetSMom = -p%TeetHSSp*SIGN( StopDef, real(TeetDef,reKi) )
      ELSE
         TeetSMom = 0.0_ReKi
      ENDIF


      ! Compute linear teeter-damper moment.

      IF ( ABS(TeetDef) > p%TeetDmpP ) THEN
         TeetDMom = -p%TeetDmp*TeetRate
      ELSE
         TeetDMom = 0.0_ReKi
      END IF
      
      

      ! Add coulomb friction to the teeter hinge.

      IF ( .NOT. EqualRealNos( TeetRate, 0.0_R8Ki ) )  THEN
         TeetFMom = 0.0_ReKi
      ELSE
         TeetFMom = -SIGN( p%TeetCDmp, real(TeetRate,reKi) )
      ENDIF


      ! Total up all the moments.

      TeetMom = TeetSMom + TeetDMom + TeetKMom + TeetFMom


   CASE ( 2_IntKi )              ! User-defined teeter spring/damper model.


      CALL UserTeet ( TeetDef, TeetRate, t, p%RootName, TeetMom )


   END SELECT


   RETURN
END SUBROUTINE Teeter
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine computes the tail-furl moment due to tail-furl deflection and rate.
SUBROUTINE TFurling( t, p, TFrlDef, TFrlRate, TFrlMom )
      ! Passed Variables:
   REAL(DbKi), INTENT(IN)             :: t                                       !< simulation time
   TYPE(ED_ParameterType), INTENT(IN) :: p                                       !< parameters from the structural dynamics module
   REAL(R8Ki), INTENT(IN )            :: TFrlDef                                 !< The tail-furl deflection, QT(DOF_TFrl).
   REAL(ReKi), INTENT(OUT)            :: TFrlMom                                 !< The total moment supplied by the springs, and dampers.
   REAL(R8Ki), INTENT(IN )            :: TFrlRate                                !< The tail-furl rate, QDT(DOF_TFrl).
      ! Local variables:
   REAL(ReKi)                         :: TFrlDMom                                ! The moment supplied by the tail-furl dampers.
   REAL(ReKi)                         :: TFrlSMom                                ! The moment supplied by the tail-furl springs.

   SELECT CASE ( p%TFrlMod ) ! Which tail-furl model are we using?

      CASE ( 0_IntKi )              ! None!

         TFrlMom = 0.0

      CASE ( 1_IntKi )              ! Standard (using inputs from the FAST furling input file).

         ! Linear spring:
         TFrlSMom = -p%TFrlSpr*TFrlDef

         ! Add spring-stops:
         IF ( TFrlDef > p%TFrlUSSP )  THEN      ! Up-stop
            TFrlSMom = TFrlSMom - p%TFrlUSSpr*( TFrlDef - p%TFrlUSSP )
         ELSEIF ( TFrlDef < p%TFrlDSSP )  THEN  ! Down-stop
            TFrlSMom = TFrlSMom - p%TFrlDSSpr*( TFrlDef - p%TFrlDSSP )
         ENDIF

         ! Linear damper:
         TFrlDMom = -p%TFrlDmp*TFrlRate

         ! Add damper-stops:
         IF ( TFrlDef > p%TFrlUSDP )  THEN      ! Up-stop
            TFrlDMom = TFrlDMom - p%TFrlUSDmp*TFrlRate
         ELSEIF ( TFrlDef < p%TFrlDSDP )  THEN  ! Down-stop
            TFrlDMom = TFrlDMom - p%TFrlDSDmp*TFrlRate
         ENDIF

         ! Total up all the moments.
         TFrlMom = TFrlSMom + TFrlDMom

      CASE ( 2 )              ! User-defined tail-furl spring/damper model.

         CALL UserTFrl ( TFrlDef, TFrlRate, t, p%RootName, TFrlMom )

   END SELECT
END SUBROUTINE TFurling
!----------------------------------------------------------------------------------------------------------------------------------
!> This function calculates the sign (+/-1) of the low-speed shaft torque for
!!   this time step.  MomLPRot is the moment on the
!!   low-speed shaft at the teeter pin caused by the rotor.
FUNCTION SignLSSTrq( p, m )

      ! Passed variables

   TYPE(ED_ParameterType),  INTENT(IN)  :: p                 !< Parameters
   TYPE(ED_MiscVarType),    INTENT(IN)  :: m                 !< Misc variables

   INTEGER(IntKi)                       :: SignLSSTrq        !< The sign of the LSS_Trq, output from this function

      ! Local variables

   REAL(ReKi)                           :: MomLPRot  (3)     ! The total moment on the low-speed shaft at point P caused by the rotor.
   INTEGER(IntKi)                       :: I                 ! loop counter


   MomLPRot = m%RtHS%MomLPRott ! Initialize MomLPRot using MomLPRott
   DO I = 1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs

      MomLPRot = MomLPRot + m%RtHS%PMomLPRot(:,p%DOFs%SrtPS(I))*m%QD2T(p%DOFs%SrtPS(I))  ! Add the moments associated with the accelerations of the DOFs

   ENDDO             ! I - All active (enabled) DOFs

      ! MomLProt has now been found.  Now dot this with e1 to get the
      !   low-speed shaft torque and take the SIGN of the result:

   SignLSSTrq = NINT( SIGN( 1.0_R8Ki, DOT_PRODUCT( MomLPRot, m%CoordSys%e1 ) ) )

END FUNCTION SignLSSTrq
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is used to calculate the positions stored in other states that are used in both the
!! CalcOutput and CalcContStateDeriv routines.
SUBROUTINE CalculatePositions( p, x, CoordSys, RtHSdat )
!..................................................................................................................................

      ! Passed variables
   TYPE(ED_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(ED_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
   TYPE(ED_CoordSys),            INTENT(IN   )  :: CoordSys    !< The coordinate systems that have been set for these states/time
   TYPE(ED_RtHndSide),           INTENT(INOUT)  :: RtHSdat     !< data from the RtHndSid module (contains positions to be set)

      !Local variables
   !REAL(R8Ki)                   :: rQ        (3)                                   ! Position vector from inertial frame origin to apex of rotation (point Q).

   INTEGER(IntKi)               :: J                                               ! Counter for elements
   INTEGER(IntKi)               :: K                                               ! Counter for blades

      !-------------------------------------------------------------------------------------------------
      ! Positions
      !-------------------------------------------------------------------------------------------------

      ! Define the position vectors between the various points on the wind turbine
      !   that are not dependent on the distributed tower or blade parameters:

   RtHSdat%rZ    = x%QT(DOF_Sg)* CoordSys%z1 + x%QT(DOF_Hv)* CoordSys%z2 - x%QT(DOF_Sw)* CoordSys%z3                          ! Position vector from inertia frame origin to platform reference (point Z).
   RtHSdat%rZY   = p%rZYzt*  CoordSys%a2 + p%PtfmCMxt*CoordSys%a1 - p%PtfmCMyt*CoordSys%a3                                    ! Position vector from platform reference (point Z) to platform mass center (point Y).      
   RtHSdat%rZT0  = p%rZT0zt* CoordSys%a2                                                                                      ! Position vector from platform reference (point Z) to tower base (point T(0))
   RtHSdat%rZO   = ( x%QT(DOF_TFA1) + x%QT(DOF_TFA2)                                                        )*CoordSys%a1 &   ! Position vector from platform reference (point Z) to tower-top / base plate (point O).
                    + ( p%RefTwrHt - 0.5*(      p%AxRedTFA(1,1,p%TTopNode)*x%QT(DOF_TFA1)*x%QT(DOF_TFA1) &
                                          +     p%AxRedTFA(2,2,p%TTopNode)*x%QT(DOF_TFA2)*x%QT(DOF_TFA2) &
                                          + 2.0*p%AxRedTFA(1,2,p%TTopNode)*x%QT(DOF_TFA1)*x%QT(DOF_TFA2) &
                                          +     p%AxRedTSS(1,1,p%TTopNode)*x%QT(DOF_TSS1)*x%QT(DOF_TSS1) &
                                          +     p%AxRedTSS(2,2,p%TTopNode)*x%QT(DOF_TSS2)*x%QT(DOF_TSS2) &
                                          + 2.0*p%AxRedTSS(1,2,p%TTopNode)*x%QT(DOF_TSS1)*x%QT(DOF_TSS2)   ) )*CoordSys%a2 &
                    + ( x%QT(DOF_TSS1) + x%QT(DOF_TSS2)                                                      )*CoordSys%a3
   RtHSdat%rOU   =   p%NacCMxn*CoordSys%d1  +  p%NacCMzn  *CoordSys%d2  -  p%NacCMyn  *CoordSys%d3                            ! Position vector from tower-top / base plate (point O) to nacelle center of mass (point U).
   RtHSdat%rOV   = p%RFrlPnt_n(1)*CoordSys%d1  +  p%RFrlPnt_n(3)*CoordSys%d2  -  p%RFrlPnt_n(2)*CoordSys%d3                            ! Position vector from tower-top / base plate (point O) to specified point on rotor-furl axis (point V).
   RtHSdat%rVIMU =   p%rVIMUxn*CoordSys%rf1 +  p%rVIMUzn  *CoordSys%rf2 -   p%rVIMUyn *CoordSys%rf3                           ! Position vector from specified point on rotor-furl axis (point V) to nacelle IMU (point IMU).
   RtHSdat%rVD   =     p%rVDxn*CoordSys%rf1 +    p%rVDzn  *CoordSys%rf2 -     p%rVDyn *CoordSys%rf3                           ! Position vector from specified point on rotor-furl axis (point V) to center of mass of structure that furls with the rotor (not including rotor) (point D).
   RtHSdat%rVP   =     p%rVPxn*CoordSys%rf1 +    p%rVPzn  *CoordSys%rf2 -     p%rVPyn *CoordSys%rf3 + p%OverHang*CoordSys%c1  ! Position vector from specified point on rotor-furl axis (point V) to teeter pin (point P).
   RtHSdat%rPQ   = -p%UndSling*CoordSys%g1                                                                                    ! Position vector from teeter pin (point P) to apex of rotation (point Q).
   RtHSdat%rQC   =     p%HubCM*CoordSys%g1                                                                                    ! Position vector from apex of rotation (point Q) to hub center of mass (point C).
   RtHSdat%rOW   = p%TFrlPnt_n(1)*CoordSys%d1  + p%TFrlPnt_n(3) *CoordSys%d2 -  p%TFrlPnt_n(2)*CoordSys%d3                             ! Position vector from tower-top / base plate (point O) to specified point on  tail-furl axis (point W).
   RtHSdat%rWI   =     p%rWIxn*CoordSys%tf1 +      p%rWIzn*CoordSys%tf2 -     p%rWIyn*CoordSys%tf3                            ! Position vector from specified point on  tail-furl axis (point W) to tail boom center of mass     (point I).
   RtHSdat%rWJ   =     p%rWJxn*CoordSys%tf1 +      p%rWJzn*CoordSys%tf2 -     p%rWJyn*CoordSys%tf3                            ! Position vector from specified point on  tail-furl axis (point W) to tail fin  center of mass     (point J).
   RtHSdat%rPC   = RtHSdat%rPQ + RtHSdat%rQC                                                                                  ! Position vector from teeter pin (point P) to hub center of mass (point C).
   RtHSdat%rT0O  = RtHSdat%rZO - RtHSdat%rZT0                                                                                 ! Position vector from the tower base (point T(0)) to tower-top / base plate (point O).
   RtHSdat%rO    = RtHSdat%rZ  + RtHSdat%rZO                                                                                  ! Position vector from inertial frame origin to tower-top / base plate (point O).
   RtHSdat%rV    = RtHSdat%rO  + RtHSdat%rOV                                                                                  ! Position vector from inertial frame origin to specified point on rotor-furl axis (point V)
   !RtHSdat%rP    = RtHSdat%rO  + RtHSdat%rOV + RtHSdat%rVP                                                                   ! Position vector from inertial frame origin to teeter pin (point P).
   RtHSdat%rP    = RtHSdat%rV  + RtHSdat%rVP                                                                                  ! Position vector from inertial frame origin to teeter pin (point P).
   RtHSdat%rQ    = RtHSdat%rP  + RtHSdat%rPQ                                                                                  ! Position vector from inertial frame origin to apex of rotation (point Q).
   RtHSdat%rJ    = RtHSdat%rO  + RtHSdat%rOW + RtHSdat%rWJ                                                                    ! Position vector from inertial frame origin to tail fin center of mass (point J).


   DO K = 1,p%NumBl ! Loop through all blades

      ! Calculate the position vector of the tip:
      RtHSdat%rS0S(:,K,p%TipNode) = ( p%TwistedSF(K,1,1,p%TipNode,0)*x%QT( DOF_BF(K,1) ) &                                       ! Position vector from the blade root (point S(0)) to the blade tip (point S(p%BldFlexL)).
                                    + p%TwistedSF(K,1,2,p%TipNode,0)*x%QT( DOF_BF(K,2) ) &
                                    + p%TwistedSF(K,1,3,p%TipNode,0)*x%QT( DOF_BE(K,1) )                     )*CoordSys%j1(K,:) &
                                  + ( p%TwistedSF(K,2,1,p%TipNode,0)*x%QT( DOF_BF(K,1) ) &
                                    + p%TwistedSF(K,2,2,p%TipNode,0)*x%QT( DOF_BF(K,2) ) &
                                    + p%TwistedSF(K,2,3,p%TipNode,0)*x%QT( DOF_BE(K,1) )                     )*CoordSys%j2(K,:) &
                                  + ( p%BldFlexL - 0.5* &
                                  (      p%AxRedBld(K,1,1,p%TipNode)*x%QT( DOF_BF(K,1) )*x%QT( DOF_BF(K,1) ) &
                                    +    p%AxRedBld(K,2,2,p%TipNode)*x%QT( DOF_BF(K,2) )*x%QT( DOF_BF(K,2) ) &
                                    +    p%AxRedBld(K,3,3,p%TipNode)*x%QT( DOF_BE(K,1) )*x%QT( DOF_BE(K,1) ) &
                                    + 2.*p%AxRedBld(K,1,2,p%TipNode)*x%QT( DOF_BF(K,1) )*x%QT( DOF_BF(K,2) ) &
                                    + 2.*p%AxRedBld(K,2,3,p%TipNode)*x%QT( DOF_BF(K,2) )*x%QT( DOF_BE(K,1) ) &
                                    + 2.*p%AxRedBld(K,1,3,p%TipNode)*x%QT( DOF_BF(K,1) )*x%QT( DOF_BE(K,1) ) ) )*CoordSys%j3(K,:)
      RtHSdat%rQS (:,K,p%TipNode) = RtHSdat%rS0S(:,K,p%TipNode) + p%HubRad*CoordSys%j3(K,:)                                      ! Position vector from apex of rotation (point Q) to the blade tip (point S(p%BldFlexL)).
      RtHSdat%rS  (:,K,p%TipNode) = RtHSdat%rQS (:,K,p%TipNode) + RtHSdat%rQ                                                     ! Position vector from inertial frame origin      to the blade tip (point S(p%BldFlexL)).
      
      ! position vectors for blade root node:
      RtHSdat%rQS (:,K,0) = p%HubRad*CoordSys%j3(K,:)    
      RtHSdat%rS  (:,K,0) = p%HubRad*CoordSys%j3(K,:) + RtHSdat%rQ
      
      
         ! Calculate the position vector from the teeter pin to the blade root:
   
      RtHSdat%rPS0(:,K) = RtHSdat%rPQ + p%HubRad*CoordSys%j3(K,:)   ! Position vector from teeter pin (point P) to blade root (point S(0)).

      
      DO J = 1,p%BldNodes ! Loop through the blade nodes / elements


      ! Calculate the position vector of the current node:

         RtHSdat%rS0S(:,K,J) = (  p%TwistedSF(K,1,1,J,0)*x%QT( DOF_BF(K,1) ) &                                                   ! Position vector from the blade root (point S(0)) to the current node (point S(RNodes(J)).
                                + p%TwistedSF(K,1,2,J,0)*x%QT( DOF_BF(K,2) ) &
                                + p%TwistedSF(K,1,3,J,0)*x%QT( DOF_BE(K,1) )                          )*CoordSys%j1(K,:) &
                            + (   p%TwistedSF(K,2,1,J,0)*x%QT( DOF_BF(K,1) ) &
                                + p%TwistedSF(K,2,2,J,0)*x%QT( DOF_BF(K,2) ) &
                                + p%TwistedSF(K,2,3,J,0)*x%QT( DOF_BE(K,1) )                          )*CoordSys%j2(K,:) &
                            + (  p%RNodes(J) - 0.5* &
                              (      p%AxRedBld(K,1,1,J)*x%QT( DOF_BF(K,1) )*x%QT( DOF_BF(K,1) ) &
                               +     p%AxRedBld(K,2,2,J)*x%QT( DOF_BF(K,2) )*x%QT( DOF_BF(K,2) ) &
                               +     p%AxRedBld(K,3,3,J)*x%QT( DOF_BE(K,1) )*x%QT( DOF_BE(K,1) ) &
                               + 2.0*p%AxRedBld(K,1,2,J)*x%QT( DOF_BF(K,1) )*x%QT( DOF_BF(K,2) ) &
                               + 2.0*p%AxRedBld(K,2,3,J)*x%QT( DOF_BF(K,2) )*x%QT( DOF_BE(K,1) ) &
                               + 2.0*p%AxRedBld(K,1,3,J)*x%QT( DOF_BF(K,1) )*x%QT( DOF_BE(K,1) )    ) )*CoordSys%j3(K,:)
         RtHSdat%rQS (:,K,J) = RtHSdat%rS0S(:,K,J) + p%HubRad*CoordSys%j3(K,:)                                                ! Position vector from apex of rotation (point Q) to the current node (point S(RNodes(J)).
         RtHSdat%rS  (:,K,J) = RtHSdat%rQS (:,K,J) + RtHSdat%rQ                                                               ! Position vector from inertial frame origin      to the current node (point S(RNodes(J)).


      END DO !J = 1,p%BldNodes ! Loop through the blade nodes / elements



   END DO !K = 1,p%NumBl
 

   

   !----------------------------------------------------------------------------------------------------
   ! Get the tower element positions
   !----------------------------------------------------------------------------------------------------
   RtHSdat%rZT (:,0) = RtHSdat%rZT0
   DO J = 1,p%TwrNodes  ! Loop through the tower nodes / elements


      ! Calculate the position vector of the current node:

      RtHSdat%rT0T(:,J) = ( p%TwrFASF(1,J,0)*x%QT(DOF_TFA1) + p%TwrFASF(2,J,0)*x%QT(DOF_TFA2)           )*CoordSys%a1 &       ! Position vector from base of flexible portion of tower (point T(0)) to current node (point T(J)).
                        + ( p%HNodes(J) - 0.5*(     p%AxRedTFA(1,1,J)*x%QT(DOF_TFA1)*x%QT(DOF_TFA1) &
                                              +     p%AxRedTFA(2,2,J)*x%QT(DOF_TFA2)*x%QT(DOF_TFA2) &
                                              + 2.0*p%AxRedTFA(1,2,J)*x%QT(DOF_TFA1)*x%QT(DOF_TFA2) &
                                              +     p%AxRedTSS(1,1,J)*x%QT(DOF_TSS1)*x%QT(DOF_TSS1) &
                                              +     p%AxRedTSS(2,2,J)*x%QT(DOF_TSS2)*x%QT(DOF_TSS2) &
                                              + 2.0*p%AxRedTSS(1,2,J)*x%QT(DOF_TSS1)*x%QT(DOF_TSS2)   ) )*CoordSys%a2 &
                        + ( p%TwrSSSF(1,J,0)*x%QT(DOF_TSS1) + p%TwrSSSF(2,J,0)*x%QT(DOF_TSS2)           )*CoordSys%a3
      RtHSdat%rZT (:,J) = RtHSdat%rZT0 + RtHSdat%rT0T(:,J)                                                                    ! Position vector from platform reference (point Z) to the current node (point T(HNodes(J)).


      RtHSdat%rT(:,J)   = RtHSdat%rZ   + RtHSdat%rZT (:,J)                                                                    ! Position vector from inertial frame origin        to the current node (point T(HNodes(J)).

   END DO


END SUBROUTINE CalculatePositions
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is used to calculate the angular positions, velocities, and partial accelerations stored in other states that are used in
!! both the CalcOutput and CalcContStateDeriv routines.
SUBROUTINE CalculateAngularPosVelPAcc( p, x, CoordSys, RtHSdat )
!..................................................................................................................................

      ! Passed variables
   TYPE(ED_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(ED_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
   TYPE(ED_CoordSys),            INTENT(IN   )  :: CoordSys    !< The coordinate systems that have been set for these states/time
   TYPE(ED_RtHndSide),           INTENT(INOUT)  :: RtHSdat     !< data from the RtHndSid module (contains positions to be set)

      !Local variables
   
   REAL(ReKi)                   :: AngVelHM  (3)                                   ! Angular velocity of eleMent J of blade K (body M) in the hub (body H).
!   REAL(ReKi)                   :: AngVelEN  (3)                                   ! Angular velocity of the nacelle (body N) in the inertia frame (body E for earth).
   REAL(ReKi)                   :: AngAccELt (3)                                   ! Portion of the angular acceleration of the low-speed shaft (body L) in the inertia frame (body E for earth) associated with everything but the QD2T()'s.
   INTEGER(IntKi)               :: J                                               ! Counter for elements
   INTEGER(IntKi)               :: K                                               ! Counter for blades

   !-------------------------------------------------------------------------------------------------
   ! Angular and partial angular velocities
   !-------------------------------------------------------------------------------------------------

   ! Define the angular and partial angular velocities of all of the rigid bodies in the inertia frame:
   ! NOTE: PAngVelEN(I,D,:) = the Dth-derivative of the partial angular velocity of DOF I for body N in body E.

   RtHSdat%PAngVelEX(       :,0,:) = 0.0
   RtHSdat%PAngVelEX(DOF_R   ,0,:) =  CoordSys%z1
   RtHSdat%PAngVelEX(DOF_P   ,0,:) = -CoordSys%z3
   RtHSdat%PAngVelEX(DOF_Y   ,0,:) =  CoordSys%z2
   RtHSdat%AngVelEX                =                     x%QDT(DOF_R   )*RtHSdat%PAngVelEX(DOF_R   ,0,:) &
                                                       + x%QDT(DOF_P   )*RtHSdat%PAngVelEX(DOF_P   ,0,:) &
                                                       + x%QDT(DOF_Y   )*RtHSdat%PAngVelEX(DOF_Y   ,0,:)
   RtHSdat%AngPosEX                =                     x%QT (DOF_R   )*RtHSdat%PAngVelEX(DOF_R   ,0,:) &
                                                       + x%QT (DOF_P   )*RtHSdat%PAngVelEX(DOF_P   ,0,:) &
                                                       + x%QT (DOF_Y   )*RtHSdat%PAngVelEX(DOF_Y   ,0,:)

   RtHSdat%PAngVelEB(       :,0,:) =  RtHSdat%PAngVelEX(:,0,:)
   RtHSdat%PAngVelEB(DOF_TFA1,0,:) = -p%TwrFASF(1,p%TTopNode,1)*CoordSys%a3
   RtHSdat%PAngVelEB(DOF_TSS1,0,:) =  p%TwrSSSF(1,p%TTopNode,1)*CoordSys%a1
   RtHSdat%PAngVelEB(DOF_TFA2,0,:) = -p%TwrFASF(2,p%TTopNode,1)*CoordSys%a3
   RtHSdat%PAngVelEB(DOF_TSS2,0,:) =  p%TwrSSSF(2,p%TTopNode,1)*CoordSys%a1
   RtHSdat%AngVelEB                =  RtHSdat%AngVelEX + x%QDT(DOF_TFA1)*RtHSdat%PAngVelEB(DOF_TFA1,0,:) &
                                                       + x%QDT(DOF_TSS1)*RtHSdat%PAngVelEB(DOF_TSS1,0,:) &
                                                       + x%QDT(DOF_TFA2)*RtHSdat%PAngVelEB(DOF_TFA2,0,:) &
                                                       + x%QDT(DOF_TSS2)*RtHSdat%PAngVelEB(DOF_TSS2,0,:)
   RtHSdat%AngPosXB                =                     x%QT (DOF_TFA1)*RtHSdat%PAngVelEB(DOF_TFA1,0,:) &
                                                       + x%QT (DOF_TSS1)*RtHSdat%PAngVelEB(DOF_TSS1,0,:) &
                                                       + x%QT (DOF_TFA2)*RtHSdat%PAngVelEB(DOF_TFA2,0,:) &
                                                       + x%QT (DOF_TSS2)*RtHSdat%PAngVelEB(DOF_TSS2,0,:)

   RtHSdat%PAngVelEN(       :,0,:)= RtHSdat%PAngVelEB(:,0,:)
   RtHSdat%PAngVelEN(DOF_Yaw ,0,:)= CoordSys%d2
   RtHSdat%AngVelEN               = RtHSdat%AngVelEB + x%QDT(DOF_Yaw )*RtHSdat%PAngVelEN(DOF_Yaw ,0,:)

   RtHSdat%PAngVelER(       :,0,:)= RtHSdat%PAngVelEN(:,0,:)
   RtHSdat%PAngVelER(DOF_RFrl,0,:)= CoordSys%rfa
   RtHSdat%AngVelER               = RtHSdat%AngVelEN + x%QDT(DOF_RFrl)*RtHSdat%PAngVelER(DOF_RFrl,0,:)

   RtHSdat%PAngVelEL(       :,0,:)= RtHSdat%PAngVelER(:,0,:)
   RtHSdat%PAngVelEL(DOF_GeAz,0,:)= CoordSys%c1
   RtHSdat%PAngVelEL(DOF_DrTr,0,:)= CoordSys%c1
   RtHSdat%AngVelEL               = RtHSdat%AngVelER + x%QDT(DOF_GeAz)*RtHSdat%PAngVelEL(DOF_GeAz,0,:) &
                                                           + x%QDT(DOF_DrTr)*RtHSdat%PAngVelEL(DOF_DrTr,0,:)

   RtHSdat%PAngVelEH(       :,0,:)= RtHSdat%PAngVelEL(:,0,:)
   RtHSdat%AngVelEH               = RtHSdat%AngVelEL
IF ( p%NumBl == 2 )  THEN ! 2-blader
   RtHSdat%PAngVelEH(DOF_Teet,0,:)= CoordSys%f2
   RtHSdat%AngVelEH               = RtHSdat%AngVelEH + x%QDT(DOF_Teet)*RtHSdat%PAngVelEH(DOF_Teet,0,:)
ENDIF

   RtHSdat%PAngVelEG(       :,0,:) = RtHSdat%PAngVelER(:,0,:)
   RtHSdat%PAngVelEG(DOF_GeAz,0,:) = p%GBRatio*CoordSys%c1
   RtHSdat%AngVelEG                = RtHSdat%AngVelER + x%QDT(DOF_GeAz)*RtHSdat%PAngVelEG(DOF_GeAz,0,:)

   RtHSdat%PAngVelEA(       :,0,:) = RtHSdat%PAngVelEN(:,0,:)
   RtHSdat%PAngVelEA(DOF_TFrl,0,:) = CoordSys%tfa
   RtHSdat%AngVelEA                = RtHSdat%AngVelEN + x%QDT(DOF_TFrl)*RtHSdat%PAngVelEA(DOF_TFrl,0,:)



   ! Define the 1st derivatives of the partial angular velocities of all
   !   of the rigid bodies in the inertia frame and the portion of the angular
   !   acceleration of the rigid bodies in the inertia frame associated with
   !   everything but the QD2T()'s:

   RtHSdat%PAngVelEX(       :,1,:) = 0.0
   RtHSdat%AngAccEXt               = 0.0

   RtHSdat%PAngVelEB(       :,1,:) =                  RtHSdat%PAngVelEX(:,1,:)
   RtHSdat%PAngVelEB(DOF_TFA1,1,:) = CROSS_PRODUCT(   RtHSdat%AngVelEX,                   RtHSdat%PAngVelEB(DOF_TFA1,0,:) )
   RtHSdat%PAngVelEB(DOF_TSS1,1,:) = CROSS_PRODUCT(   RtHSdat%AngVelEX,                   RtHSdat%PAngVelEB(DOF_TSS1,0,:) )
   RtHSdat%PAngVelEB(DOF_TFA2,1,:) = CROSS_PRODUCT(   RtHSdat%AngVelEX,                   RtHSdat%PAngVelEB(DOF_TFA2,0,:) )
   RtHSdat%PAngVelEB(DOF_TSS2,1,:) = CROSS_PRODUCT(   RtHSdat%AngVelEX,                   RtHSdat%PAngVelEB(DOF_TSS2,0,:) )
   RtHSdat%AngAccEBt               =                  RtHSdat%AngAccEXt + x%QDT(DOF_TFA1)*RtHSdat%PAngVelEB(DOF_TFA1,1,:) &
                                                                        + x%QDT(DOF_TSS1)*RtHSdat%PAngVelEB(DOF_TSS1,1,:) &
                                                                        + x%QDT(DOF_TFA2)*RtHSdat%PAngVelEB(DOF_TFA2,1,:) &
                                                                        + x%QDT(DOF_TSS2)*RtHSdat%PAngVelEB(DOF_TSS2,1,:)

   RtHSdat%PAngVelEN(       :,1,:) =                 RtHSdat%PAngVelEB(:,1,:)
   RtHSdat%PAngVelEN(DOF_Yaw ,1,:) = CROSS_PRODUCT(  RtHSdat%AngVelEB,                    RtHSdat%PAngVelEN(DOF_Yaw ,0,:) )
   RtHSdat%AngAccENt               =                 RtHSdat%AngAccEBt  + x%QDT(DOF_Yaw )*RtHSdat%PAngVelEN(DOF_Yaw ,1,:)

   RtHSdat%PAngVelER(       :,1,:) =                 RtHSdat%PAngVelEN(:,1,:)
   RtHSdat%PAngVelER(DOF_RFrl,1,:) = CROSS_PRODUCT(  RtHSdat%AngVelEN,                    RtHSdat%PAngVelER(DOF_RFrl,0,:) )
   RtHSdat%AngAccERt               =                 RtHSdat%AngAccENt  + x%QDT(DOF_RFrl)*RtHSdat%PAngVelER(DOF_RFrl,1,:)

   RtHSdat%PAngVelEL(       :,1,:) =                 RtHSdat%PAngVelER(:,1,:)
   RtHSdat%PAngVelEL(DOF_GeAz,1,:) = CROSS_PRODUCT(  RtHSdat%AngVelER,                    RtHSdat%PAngVelEL(DOF_GeAz,0,:) )
   RtHSdat%PAngVelEL(DOF_DrTr,1,:) = CROSS_PRODUCT(  RtHSdat%AngVelER,                    RtHSdat%PAngVelEL(DOF_DrTr,0,:) )
           AngAccELt               =                 RtHSdat%AngAccERt  + x%QDT(DOF_GeAz)*RtHSdat%PAngVelEL(DOF_GeAz,1,:) &
                                                                        + x%QDT(DOF_DrTr)*RtHSdat%PAngVelEL(DOF_DrTr,1,:)

   RtHSdat%PAngVelEH(       :,1,:) = RtHSdat%PAngVelEL(:,1,:)
   RtHSdat%AngAccEHt               =                  AngAccELt
IF ( p%NumBl == 2 )  THEN ! 2-blader
   RtHSdat%PAngVelEH(DOF_Teet,1,:) = CROSS_PRODUCT(  RtHSdat%AngVelEH,                    RtHSdat%PAngVelEH(DOF_Teet,0,:) )
   RtHSdat%AngAccEHt               =                 RtHSdat%AngAccEHt   + x%QDT(DOF_Teet)*RtHSdat%PAngVelEH(DOF_Teet,1,:)
ENDIF

   RtHSdat%PAngVelEG(       :,1,:) = RtHSdat%PAngVelER(:,1,:)
   RtHSdat%PAngVelEG(DOF_GeAz,1,:) = CROSS_PRODUCT(  RtHSdat%AngVelER,                    RtHSdat%PAngVelEG(DOF_GeAz,0,:) )
   RtHSdat%AngAccEGt               =                 RtHSdat%AngAccERt  + x%QDT(DOF_GeAz)*RtHSdat%PAngVelEG(DOF_GeAz,1,:)

   RtHSdat%PAngVelEA(       :,1,:) = RtHSdat%PAngVelEN(:,1,:)
   RtHSdat%PAngVelEA(DOF_TFrl,1,:) = CROSS_PRODUCT(  RtHSdat%AngVelEN,                    RtHSdat%PAngVelEA(DOF_TFrl,0,:) )
   RtHSdat%AngAccEAt               =                 RtHSdat%AngAccENt  + x%QDT(DOF_TFrl)*RtHSdat%PAngVelEA(DOF_TFrl,1,:)



   DO K = 1,p%NumBl ! Loop through all blades

      DO J = 0,p%TipNode ! Loop through the blade nodes / elements
      ! Define the partial angular velocities of the current node (body M(RNodes(J))) in the inertia frame:
      ! NOTE: PAngVelEM(K,J,I,D,:) = the Dth-derivative of the partial angular velocity
      !   of DOF I for body M of blade K, element J in body E.

         RtHSdat%PAngVelEM(K,J,          :,0,:) = RtHSdat%PAngVelEH(:,0,:)
         RtHSdat%PAngVelEM(K,J,DOF_BF(K,1),0,:) = - p%TwistedSF(K,2,1,J,1)*CoordSys%j1(K,:) &
                                                  + p%TwistedSF(K,1,1,J,1)*CoordSys%j2(K,:)
         RtHSdat%PAngVelEM(K,J,DOF_BF(K,2),0,:) = - p%TwistedSF(K,2,2,J,1)*CoordSys%j1(K,:) &
                                                  + p%TwistedSF(K,1,2,J,1)*CoordSys%j2(K,:)
         RtHSdat%PAngVelEM(K,J,DOF_BE(K,1),0,:) = - p%TwistedSF(K,2,3,J,1)*CoordSys%j1(K,:) &
                                                  + p%TwistedSF(K,1,3,J,1)*CoordSys%j2(K,:)
                                      AngVelHM  =     x%QDT(DOF_BF(K,1))*RtHSdat%PAngVelEM(K,J,DOF_BF(K,1),0,:) &
                                                    + x%QDT(DOF_BF(K,2))*RtHSdat%PAngVelEM(K,J,DOF_BF(K,2),0,:) &
                                                    + x%QDT(DOF_BE(K,1))*RtHSdat%PAngVelEM(K,J,DOF_BE(K,1),0,:)
          RtHSdat%AngVelEM(:,J,K              ) =  RtHSdat%AngVelEH + AngVelHM
          RtHSdat%AngPosHM(:,K,J              ) =     x%QT (DOF_BF(K,1))*RtHSdat%PAngVelEM(K,J,DOF_BF(K,1),0,:) &
                                                    + x%QT (DOF_BF(K,2))*RtHSdat%PAngVelEM(K,J,DOF_BF(K,2),0,:) &
                                                    + x%QT (DOF_BE(K,1))*RtHSdat%PAngVelEM(K,J,DOF_BE(K,1),0,:)
         RtHSdat%AngAccEKt(:,J              ,K) =  RtHSdat%AngAccEHt + x%QDT(DOF_BF(K,1))*RtHSdat%PAngVelEM(K,J,DOF_BF(K,1),1,:) & 
                                                                     + x%QDT(DOF_BF(K,2))*RtHSdat%PAngVelEM(K,J,DOF_BF(K,2),1,:) & 
                                                                     + x%QDT(DOF_BE(K,1))*RtHSdat%PAngVelEM(K,J,DOF_BE(K,1),1,:)   
 
      ! Define the 1st derivatives of the partial angular velocities of the current node (body M(RNodes(J))) in the inertia frame:

   ! NOTE: These are currently unused by the code, therefore, they need not
   !       be calculated.  Thus, they are currently commented out.  If it
   !       turns out that they are ever needed (i.e., if inertias of the
   !       blade elements are ever added, etc...) simply uncomment out these computations:
   !      RtHSdat%PAngVelEM(K,J,          :,1,:) = RtHSdat%PAngVelEH(:,1,:)
   !      RtHSdat%PAngVelEM(K,J,DOF_BF(K,1),1,:) = CROSS_PRODUCT(   RtHSdat%AngVelEH, PAngVelEM(K,J,DOF_BF(K,1),0,:) )
   !      RtHSdat%PAngVelEM(K,J,DOF_BF(K,2),1,:) = CROSS_PRODUCT(   RtHSdat%AngVelEH, PAngVelEM(K,J,DOF_BF(K,2),0,:) )
   !      RtHSdat%PAngVelEM(K,J,DOF_BE(K,1),1,:) = CROSS_PRODUCT(   RtHSdat%AngVelEH, PAngVelEM(K,J,DOF_BE(K,1),0,:) )


      END DO !J = 1,p%BldNodes ! Loop through the blade nodes / elements

   END DO !K = 1,p%NumBl


   !...............
   ! tower values:
   !...............

   DO J = 0,p%TwrNodes  ! Loop through the tower nodes / elements

      ! Define the partial angular velocities (and their 1st derivatives) of the
      !   current node (body F(HNodes(J))  in the inertia frame.
      ! Also define the overall angular velocity of the current node in the inertia frame.
      !   Also, define the portion of the angular acceleration of the current node
      !   in the inertia frame associated with everything but the QD2T()'s:

      ! NOTE: PAngVelEF(J,I,D,:) = the Dth-derivative of the partial angular velocity
      !   of DOF I for body F of element J in body E.

      RtHSdat%PAngVelEF (J,       :,0,:) = RtHSdat%PAngVelEX(:,0,:)
      RtHSdat%PAngVelEF (J,DOF_TFA1,0,:) = -p%TwrFASF(1,J,1)*CoordSys%a3
      RtHSdat%PAngVelEF (J,DOF_TSS1,0,:) =  p%TwrSSSF(1,J,1)*CoordSys%a1
      RtHSdat%PAngVelEF (J,DOF_TFA2,0,:) = -p%TwrFASF(2,J,1)*CoordSys%a3
      RtHSdat%PAngVelEF (J,DOF_TSS2,0,:) =  p%TwrSSSF(2,J,1)*CoordSys%a1

      RtHSdat%PAngVelEF (J,       :,1,:) = RtHSdat%PAngVelEX(:,1,:)
      RtHSdat%PAngVelEF (J,DOF_TFA1,1,:) = CROSS_PRODUCT(  RtHSdat%AngVelEX  ,  RtHSdat%PAngVelEF(J,DOF_TFA1,0,:) )
      RtHSdat%PAngVelEF (J,DOF_TSS1,1,:) = CROSS_PRODUCT(  RtHSdat%AngVelEX  ,  RtHSdat%PAngVelEF(J,DOF_TSS1,0,:) )
      RtHSdat%PAngVelEF (J,DOF_TFA2,1,:) = CROSS_PRODUCT(  RtHSdat%AngVelEX  ,  RtHSdat%PAngVelEF(J,DOF_TFA2,0,:) )
      RtHSdat%PAngVelEF (J,DOF_TSS2,1,:) = CROSS_PRODUCT(  RtHSdat%AngVelEX  ,  RtHSdat%PAngVelEF(J,DOF_TSS2,0,:) )


      RtHSdat%AngVelEF (:,J)            =  RtHSdat%AngVelEX  + x%QDT(DOF_TFA1)*RtHSdat%PAngVelEF(J,DOF_TFA1,0,:) &
                                                             + x%QDT(DOF_TSS1)*RtHSdat%PAngVelEF(J,DOF_TSS1,0,:) &
                                                             + x%QDT(DOF_TFA2)*RtHSdat%PAngVelEF(J,DOF_TFA2,0,:) &
                                                             + x%QDT(DOF_TSS2)*RtHSdat%PAngVelEF(J,DOF_TSS2,0,:)

      RtHSdat%AngPosXF (:,J)            =                      x%QT (DOF_TFA1)*RtHSdat%PAngVelEF(J,DOF_TFA1,0,:) &
                                                             + x%QT (DOF_TSS1)*RtHSdat%PAngVelEF(J,DOF_TSS1,0,:) &
                                                             + x%QT (DOF_TFA2)*RtHSdat%PAngVelEF(J,DOF_TFA2,0,:) &
                                                             + x%QT (DOF_TSS2)*RtHSdat%PAngVelEF(J,DOF_TSS2,0,:)
      RtHSdat%AngPosEF (:,J)            =  RtHSdat%AngPosEX  + RtHSdat%AngPosXF(:,J)
      RtHSdat%AngAccEFt(:,J)            =  RtHSdat%AngAccEXt + x%QDT(DOF_TFA1)*RtHSdat%PAngVelEF(J,DOF_TFA1,1,:) &
                                                             + x%QDT(DOF_TSS1)*RtHSdat%PAngVelEF(J,DOF_TSS1,1,:) &
                                                             + x%QDT(DOF_TFA2)*RtHSdat%PAngVelEF(J,DOF_TFA2,1,:) &
                                                             + x%QDT(DOF_TSS2)*RtHSdat%PAngVelEF(J,DOF_TSS2,1,:)

   END DO ! J


END SUBROUTINE CalculateAngularPosVelPAcc
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is used to calculate the linear velocities and accelerations stored in other states that are used in
!! both the CalcOutput and CalcContStateDeriv routines.
SUBROUTINE CalculateLinearVelPAcc( p, x, CoordSys, RtHSdat )
!..................................................................................................................................

      ! Passed variables
   TYPE(ED_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(ED_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
   TYPE(ED_CoordSys),            INTENT(IN   )  :: CoordSys    !< The coordinate systems that have been set for these states/time
   TYPE(ED_RtHndSide),           INTENT(INOUT)  :: RtHSdat     !< data from the RtHndSid module (contains positions to be set)

      ! Local variables
   REAL(ReKi)                   :: LinAccEPt (3)                                   ! "Portion of the linear acceleration of the teeter pin (point P) in the inertia frame (body E for earth) associated with everything but the QD2T()'s"
   REAL(ReKi)                   :: LinAccEQt (3)                                   ! "Portion of the linear acceleration of the apex of rotation (point Q) in the inertia frame (body E for earth) associated with everything but the QD2T()'s"
   REAL(ReKi)                   :: LinAccEVt (3)                                   ! "Portion of the linear acceleration of the selected point on the rotor-furl axis (point V) in the inertia frame (body E for earth) associated with everything but the QD2T()'s"
   REAL(ReKi)                   :: LinAccEWt (3)                                   ! "Portion of the linear acceleration of the selected point on the  tail-furl axis (point W) in the inertia frame (body E for earth) associated with everything but the QD2T()'s"
   REAL(ReKi)                   :: LinVelHS  (3)                                   ! "Relative linear velocity of the current point on the current blade (point S) in the hub frame (body H)"
   REAL(ReKi)                   :: LinVelXO  (3)                                   ! "Relative linear velocity of the tower-top / base plate (point O) in the platform (body X)"
   REAL(ReKi)                   :: LinVelXT  (3)                                   ! "Relative linear velocity of the current point on the tower (point T) in the platform (body X)"

   REAL(ReKi)                   :: EwAXrWI   (3)                                   ! = AngVelEA X rWI
   REAL(ReKi)                   :: EwAXrWJ   (3)                                   ! = AngVelEA X rWJ
   REAL(ReKi)                   :: EwHXrPQ   (3)                                   ! = AngVelEH X rPQ
   REAL(ReKi)                   :: EwHXrQC   (3)                                   ! = AngVelEH X rQC
   REAL(ReKi)                   :: EwHXrQS   (3)                                   ! = AngVelEH X rQS of the current blade point S.
   REAL(ReKi)                   :: EwNXrOU   (3)                                   ! = AngVelEN X rOU
   REAL(ReKi)                   :: EwNXrOV   (3)                                   ! = AngVelEN X rOV
   REAL(ReKi)                   :: EwNXrOW   (3)                                   ! = AngVelEN X rOW
   REAL(ReKi)                   :: EwRXrVD   (3)                                   ! = AngVelER X rVD
   REAL(ReKi)                   :: EwRXrVIMU (3)                                   ! = AngVelER X rVIMU
   REAL(ReKi)                   :: EwRXrVP   (3)                                   ! = AngVelER X rVP
   REAL(ReKi)                   :: EwXXrZO   (3)                                   ! = AngVelEX X rZO
   REAL(ReKi)                   :: EwXXrZT   (3)                                   ! = AngVelEX X rZT
   REAL(ReKi)                   :: EwXXrZY   (3)                                   ! = AngVelEX X rZY

   REAL(ReKi)                   :: TmpVec0   (3)                                   ! A temporary vector used in various computations.
   REAL(ReKi)                   :: TmpVec1   (3)                                   ! A temporary vector used in various computations.
   REAL(ReKi)                   :: TmpVec2   (3)                                   ! A temporary vector used in various computations.
   REAL(ReKi)                   :: TmpVec3   (3)                                   ! A temporary vector used in various computations.
   REAL(ReKi)                   :: TmpVec4   (3)                                   ! A temporary vector used in various computations.

   INTEGER(IntKi)               :: I                                               ! Loops through some or all of the DOFs
   INTEGER(IntKi)               :: J                                               ! Counter for elements
   INTEGER(IntKi)               :: K                                               ! Counter for blades


      ! Initializations:

   RtHSdat%LinAccECt   = 0.0
   RtHSdat%LinAccEDt   = 0.0
   RtHSdat%LinAccEIMUt = 0.0
   RtHSdat%LinAccEIt   = 0.0
   RtHSdat%LinAccEJt   = 0.0
   RtHSdat%LinAccEOt   = 0.0
           LinAccEPt   = 0.0
           LinAccEQt   = 0.0
   RtHSdat%LinAccESt   = 0.0
   RtHSdat%LinAccETt   = 0.0
   RtHSdat%LinAccEUt   = 0.0
           LinAccEVt   = 0.0
           LinAccEWt   = 0.0
   RtHSdat%LinAccEYt   = 0.0
   RtHSdat%LinAccEZt   = 0.0


   !-------------------------------------------------------------------------------------------------
   ! Partial linear velocities and accelerations
   !-------------------------------------------------------------------------------------------------

      ! Define the partial linear velocities (and their 1st derivatives) of all of
      !   the points on the wind turbine in the inertia frame that are not
      !   dependent on the distributed tower or blade parameters.  Also, define
      !   the portion of the linear acceleration of the points in the inertia
      !   frame associated with everything but the QD2T()'s:
      ! NOTE: PLinVelEX(I,D,:) = the Dth-derivative of the partial linear velocity
      !   of DOF I for point X in body E.

   EwXXrZY   = CROSS_PRODUCT( RtHSdat%AngVelEX, RtHSdat%rZY   ) !
   EwXXrZO   = CROSS_PRODUCT( RtHSdat%AngVelEX, RtHSdat%rZO   ) !
   EwNXrOU   = CROSS_PRODUCT( RtHSdat%AngVelEN, RtHSdat%rOU   ) !
   EwNXrOV   = CROSS_PRODUCT( RtHSdat%AngVelEN, RtHSdat%rOV   ) !
   EwRXrVD   = CROSS_PRODUCT( RtHSdat%AngVelER, RtHSdat%rVD   ) ! Cross products
   EwRXrVIMU = CROSS_PRODUCT( RtHSdat%AngVelER, RtHSdat%rVIMU ) ! that are used
   EwRXrVP   = CROSS_PRODUCT( RtHSdat%AngVelER, RtHSdat%rVP   ) ! in the following
   EwHXrPQ   = CROSS_PRODUCT( RtHSdat%AngVelEH, RtHSdat%rPQ   ) ! DO...LOOPs
   EwHXrQC   = CROSS_PRODUCT( RtHSdat%AngVelEH, RtHSdat%rQC   ) !
   EwNXrOW   = CROSS_PRODUCT( RtHSdat%AngVelEN, RtHSdat%rOW   ) !
   EwAXrWI   = CROSS_PRODUCT( RtHSdat%AngVelEA, RtHSdat%rWI   ) !
   EwAXrWJ   = CROSS_PRODUCT( RtHSdat%AngVelEA, RtHSdat%rWJ   ) !


   RtHSdat%PLinVelEZ(       :,:,:) = 0.0
   RtHSdat%PLinVelEZ(DOF_Sg  ,0,:) =  CoordSys%z1
   RtHSdat%PLinVelEZ(DOF_Sw  ,0,:) = -CoordSys%z3
   RtHSdat%PLinVelEZ(DOF_Hv  ,0,:) =  CoordSys%z2

   RtHSdat%LinVelEZ                =   x%QDT(DOF_Sg  )*RtHSdat%PLinVelEZ(DOF_Sg  ,0,:) &
                                     + x%QDT(DOF_Sw  )*RtHSdat%PLinVelEZ(DOF_Sw  ,0,:) &
                                     + x%QDT(DOF_Hv  )*RtHSdat%PLinVelEZ(DOF_Hv  ,0,:)


   RtHSdat%PLinVelEY(       :,:,:) = RtHSdat%PLinVelEZ(:,:,:)
   DO I = 1,NPX   ! Loop through all DOFs associated with the angular motion of the platform (body X)

      TmpVec0 = CROSS_PRODUCT( RtHSdat%PAngVelEX(PX(I)   ,0,:), RtHSdat%rZY  )
      TmpVec1 = CROSS_PRODUCT( RtHSdat%PAngVelEX(PX(I)   ,0,:),     EwXXrZY  )

      RtHSdat%PLinVelEY(PX(I),0,:) = TmpVec0   +                       RtHSdat%PLinVelEY(PX(I)   ,0,:)
      RtHSdat%PLinVelEY(PX(I),1,:) = TmpVec1   +                       RtHSdat%PLinVelEY(PX(I)   ,1,:)

       RtHSdat%LinAccEYt           = RtHSdat%LinAccEYt + x%QDT(PX(I) )*RtHSdat%PLinVelEY(PX(I)   ,1,:)

   ENDDO          ! I - all DOFs associated with the angular motion of the platform (body X)


   RtHSdat%PLinVelEO(       :,:,:) = RtHSdat%PLinVelEZ(:,:,:)
   RtHSdat%PLinVelEO(DOF_TFA1,0,:) = CoordSys%a1 - (   p%AxRedTFA(1,1,p%TTopNode)* x%QT(DOF_TFA1) &
                                                     + p%AxRedTFA(1,2,p%TTopNode)* x%QT(DOF_TFA2)   )*CoordSys%a2
   RtHSdat%PLinVelEO(DOF_TSS1,0,:) = CoordSys%a3 - (   p%AxRedTSS(1,1,p%TTopNode)* x%QT(DOF_TSS1) &
                                                     + p%AxRedTSS(1,2,p%TTopNode)* x%QT(DOF_TSS2)   )*CoordSys%a2
   RtHSdat%PLinVelEO(DOF_TFA2,0,:) = CoordSys%a1 - (   p%AxRedTFA(2,2,p%TTopNode)* x%QT(DOF_TFA2) &
                                                     + p%AxRedTFA(1,2,p%TTopNode)* x%QT(DOF_TFA1)   )*CoordSys%a2
   RtHSdat%PLinVelEO(DOF_TSS2,0,:) = CoordSys%a3 - (   p%AxRedTSS(2,2,p%TTopNode)* x%QT(DOF_TSS2) &
                                                     + p%AxRedTSS(1,2,p%TTopNode)* x%QT(DOF_TSS1)   )*CoordSys%a2

   TmpVec1 = CROSS_PRODUCT(   RtHSdat%AngVelEX   , RtHSdat%PLinVelEO(DOF_TFA1,0,:) )
   TmpVec2 = CROSS_PRODUCT(   RtHSdat%AngVelEX   , RtHSdat%PLinVelEO(DOF_TSS1,0,:) )
   TmpVec3 = CROSS_PRODUCT(   RtHSdat%AngVelEX   , RtHSdat%PLinVelEO(DOF_TFA2,0,:) )
   TmpVec4 = CROSS_PRODUCT(   RtHSdat%AngVelEX   , RtHSdat%PLinVelEO(DOF_TSS2,0,:) )

   RtHSdat%PLinVelEO(DOF_TFA1,1,:) = TmpVec1 - (   p%AxRedTFA(1,1,p%TTopNode)*x%QDT(DOF_TFA1) &
                                                 + p%AxRedTFA(1,2,p%TTopNode)*x%QDT(DOF_TFA2)   )*CoordSys%a2
   RtHSdat%PLinVelEO(DOF_TSS1,1,:) = TmpVec2 - (   p%AxRedTSS(1,1,p%TTopNode)*x%QDT(DOF_TSS1) &
                                                 + p%AxRedTSS(1,2,p%TTopNode)*x%QDT(DOF_TSS2)   )*CoordSys%a2
   RtHSdat%PLinVelEO(DOF_TFA2,1,:) = TmpVec3 - (   p%AxRedTFA(2,2,p%TTopNode)*x%QDT(DOF_TFA2) &
                                                 + p%AxRedTFA(1,2,p%TTopNode)*x%QDT(DOF_TFA1)   )*CoordSys%a2
   RtHSdat%PLinVelEO(DOF_TSS2,1,:) = TmpVec4 - (   p%AxRedTSS(2,2,p%TTopNode)*x%QDT(DOF_TSS2) &
                                                 + p%AxRedTSS(1,2,p%TTopNode)*x%QDT(DOF_TSS1)   )*CoordSys%a2

    LinVelXO               =              x%QDT(DOF_TFA1)*RtHSdat%PLinVelEO(DOF_TFA1,0,:) &
                                        + x%QDT(DOF_TSS1)*RtHSdat%PLinVelEO(DOF_TSS1,0,:) &
                                        + x%QDT(DOF_TFA2)*RtHSdat%PLinVelEO(DOF_TFA2,0,:) &
                                        + x%QDT(DOF_TSS2)*RtHSdat%PLinVelEO(DOF_TSS2,0,:)
    RtHSdat%LinAccEOt              =      x%QDT(DOF_TFA1)*RtHSdat%PLinVelEO(DOF_TFA1,1,:) &
                                        + x%QDT(DOF_TSS1)*RtHSdat%PLinVelEO(DOF_TSS1,1,:) &
                                        + x%QDT(DOF_TFA2)*RtHSdat%PLinVelEO(DOF_TFA2,1,:) &
                                        + x%QDT(DOF_TSS2)*RtHSdat%PLinVelEO(DOF_TSS2,1,:)
    
   RtHSdat%LinVelEO = LinVelXO + RtHSdat%LinVelEZ
   DO I = 1,NPX   ! Loop through all DOFs associated with the angular motion of the platform (body X)

      TmpVec0 = CROSS_PRODUCT( RtHSdat%PAngVelEX(PX(I)   ,0,:), RtHSdat%rZO                 )
      TmpVec1 = CROSS_PRODUCT( RtHSdat%PAngVelEX(PX(I)   ,0,:),     EwXXrZO + LinVelXO      )

      RtHSdat%PLinVelEO(PX(I),0,:) = TmpVec0    +                       RtHSdat%PLinVelEO(PX(I)   ,0,:)
      RtHSdat%PLinVelEO(PX(I),1,:) = TmpVec1    +                       RtHSdat%PLinVelEO(PX(I)   ,1,:)

      RtHSdat%LinVelEO             =  RtHSdat%LinVelEO  + x%QDT(PX(I) )*RtHSdat%PLinVelEO(PX(I)   ,0,:)
      RtHSdat%LinAccEOt            =  RtHSdat%LinAccEOt + x%QDT(PX(I) )*RtHSdat%PLinVelEO(PX(I)   ,1,:)

   ENDDO          ! I - all DOFs associated with the angular motion of the platform (body X)
                     

   RtHSdat%PLinVelEU(       :,:,:) = RtHSdat%PLinVelEO(:,:,:)
   DO I = 1,NPN   ! Loop through all DOFs associated with the angular motion of the nacelle (body N)

      TmpVec0 = CROSS_PRODUCT( RtHSdat%PAngVelEN(PN(I)   ,0,:), RtHSdat%rOU                 )
      TmpVec1 = CROSS_PRODUCT( RtHSdat%PAngVelEN(PN(I)   ,0,:),     EwNXrOU                 )
      TmpVec2 = CROSS_PRODUCT( RtHSdat%PAngVelEN(PN(I)   ,1,:), RtHSdat%rOU                 )

      RtHSdat%PLinVelEU(PN(I),0,:) = TmpVec0    +               RtHSdat%PLinVelEU(PN(I)   ,0,:)
      RtHSdat%PLinVelEU(PN(I),1,:) = TmpVec1    + TmpVec2 +     RtHSdat%PLinVelEU(PN(I)   ,1,:)

       RtHSdat%LinAccEUt           =  RtHSdat%LinAccEUt + x%QDT(PN(I) )*RtHSdat%PLinVelEU(PN(I)   ,1,:)

   ENDDO          ! I - all DOFs associated with the angular motion of the nacelle (body N)


   RtHSdat%PLinVelEV(       :,:,:) = RtHSdat%PLinVelEO(:,:,:)
   DO I = 1,NPN   ! Loop through all DOFs associated with the angular motion of the nacelle (body N)

      TmpVec0 = CROSS_PRODUCT( RtHSdat%PAngVelEN(PN(I)   ,0,:), RtHSdat%rOV                 )
      TmpVec1 = CROSS_PRODUCT( RtHSdat%PAngVelEN(PN(I)   ,0,:),     EwNXrOV                 )
      TmpVec2 = CROSS_PRODUCT( RtHSdat%PAngVelEN(PN(I)   ,1,:), RtHSdat%rOV                 )

      RtHSdat%PLinVelEV(PN(I),0,:) = TmpVec0    +               RtHSdat%PLinVelEV(PN(I)   ,0,:)
      RtHSdat%PLinVelEV(PN(I),1,:) = TmpVec1    + TmpVec2 +     RtHSdat%PLinVelEV(PN(I)   ,1,:)

       LinAccEVt                   =  LinAccEVt + x%QDT(PN(I) )*RtHSdat%PLinVelEV(PN(I)   ,1,:)

   ENDDO          ! I - all DOFs associated with the angular motion of the nacelle (body N)


   RtHSdat%PLinVelED(       :,:,:) = RtHSdat%PLinVelEV(:,:,:)
   DO I = 1,NPR   ! Loop through all DOFs associated with the angular motion of the structure that furls with the rotor (not including rotor) (body R)

      TmpVec0 = CROSS_PRODUCT( RtHSdat%PAngVelER(PR(I)   ,0,:), RtHSdat%rVD                 )
      TmpVec1 = CROSS_PRODUCT( RtHSdat%PAngVelER(PR(I)   ,0,:),     EwRXrVD                 )
      TmpVec2 = CROSS_PRODUCT( RtHSdat%PAngVelER(PR(I)   ,1,:), RtHSdat%rVD                 )

      RtHSdat%PLinVelED(PR(I),0,:) = TmpVec0    +                       RtHSdat%PLinVelED(PR(I)   ,0,:)
      RtHSdat%PLinVelED(PR(I),1,:) = TmpVec1    + TmpVec2 +             RtHSdat%PLinVelED(PR(I)   ,1,:)

      RtHSdat%LinAccEDt            =  RtHSdat%LinAccEDt + x%QDT(PR(I) )*RtHSdat%PLinVelED(PR(I)   ,1,:)

   ENDDO          ! I - all DOFs associated with the angular motion of the structure that furls with the rotor (not including rotor) (body R)


   RtHSdat%PLinVelEIMU(     :,:,:) = RtHSdat%PLinVelEV(:,:,:)
    RtHSdat%LinVelEIMU             =  RtHSdat%LinVelEZ
   DO I = 1,NPR   ! Loop through all DOFs associated with the angular motion of the structure that furls with the rotor (not including rotor) (body R)

      TmpVec0 = CROSS_PRODUCT( RtHSdat%PAngVelER(PR(I)   ,0,:), RtHSdat%rVIMU               )
      TmpVec1 = CROSS_PRODUCT( RtHSdat%PAngVelER(PR(I)   ,0,:),     EwRXrVIMU               )
      TmpVec2 = CROSS_PRODUCT( RtHSdat%PAngVelER(PR(I)   ,1,:), RtHSdat%rVIMU               )

      RtHSdat%PLinVelEIMU(PR(I),0,:) = TmpVec0    +                         RtHSdat%PLinVelEIMU(PR(I) ,0,:)
      RtHSdat%PLinVelEIMU(PR(I),1,:) = TmpVec1    + TmpVec2 +               RtHSdat%PLinVelEIMU(PR(I) ,1,:)

      RtHSdat%LinVelEIMU             =  RtHSdat%LinVelEIMU  + x%QDT(PR(I) )*RtHSdat%PLinVelEIMU(PR(I) ,0,:)
      RtHSdat%LinAccEIMUt            =  RtHSdat%LinAccEIMUt + x%QDT(PR(I) )*RtHSdat%PLinVelEIMU(PR(I) ,1,:)

   ENDDO          ! I - all DOFs associated with the angular motion of the structure that furls with the rotor (not including rotor) (body R)


   RtHSdat%PLinVelEP(       :,:,:) = RtHSdat%PLinVelEV(:,:,:)
   DO I = 1,NPR   ! Loop through all DOFs associated with the angular motion of the structure that furls with the rotor (not including rotor) (body R)

      TmpVec0 = CROSS_PRODUCT(             RtHSdat%PAngVelER(PR(I)   ,0,:),     RtHSdat%rVP                 )
      TmpVec1 = CROSS_PRODUCT(             RtHSdat%PAngVelER(PR(I)   ,0,:), EwRXrVP                 )
      TmpVec2 = CROSS_PRODUCT(             RtHSdat%PAngVelER(PR(I)   ,1,:),     RtHSdat%rVP                 )

      RtHSdat%PLinVelEP(PR(I),0,:) = TmpVec0    +               RtHSdat%PLinVelEP(PR(I)   ,0,:)
      RtHSdat%PLinVelEP(PR(I),1,:) = TmpVec1    + TmpVec2 +     RtHSdat%PLinVelEP(PR(I)   ,1,:)

       LinAccEPt           =  LinAccEPt + x%QDT(PR(I) )*RtHSdat%PLinVelEP(PR(I)   ,1,:)

   ENDDO          ! I - all DOFs associated with the angular motion of the structure that furls with the rotor (not including rotor) (body R)


   RtHSdat%PLinVelEQ(       :,:,:) = RtHSdat%PLinVelEP(:,:,:)
    RtHSdat%LinVelEQ               =  RtHSdat%LinVelEZ
   DO I = 1,p%NPH   ! Loop through all DOFs associated with the angular motion of the hub (body H)

      TmpVec0 = CROSS_PRODUCT( RtHSdat%PAngVelEH(p%PH(I)   ,0,:),   RtHSdat%rPQ  )
      TmpVec1 = CROSS_PRODUCT( RtHSdat%PAngVelEH(p%PH(I)   ,0,:),       EwHXrPQ  )
      TmpVec2 = CROSS_PRODUCT( RtHSdat%PAngVelEH(p%PH(I)   ,1,:),   RtHSdat%rPQ  )

      RtHSdat%PLinVelEQ(p%PH(I),0,:) = TmpVec0    +                 RtHSdat%PLinVelEQ(p%PH(I)   ,0,:)
      RtHSdat%PLinVelEQ(p%PH(I),1,:) = TmpVec1    + TmpVec2 +       RtHSdat%PLinVelEQ(p%PH(I)   ,1,:)

      RtHSdat%LinVelEQ               =  RtHSdat%LinVelEQ  + x%QDT(p%PH(I) )*RtHSdat%PLinVelEQ(p%PH(I)   ,0,:)
      LinAccEQt                      =          LinAccEQt + x%QDT(p%PH(I) )*RtHSdat%PLinVelEQ(p%PH(I)   ,1,:)
      
   ENDDO          ! I - all DOFs associated with the angular motion of the hub (body H)


   RtHSdat%PLinVelEC(       :,:,:) = RtHSdat%PLinVelEQ(:,:,:)
   DO I = 1,p%NPH   ! Loop through all DOFs associated with the angular motion of the hub (body H)

      TmpVec0 = CROSS_PRODUCT( RtHSdat%PAngVelEH(p%PH(I)   ,0,:), RtHSdat%rQC )
      TmpVec1 = CROSS_PRODUCT( RtHSdat%PAngVelEH(p%PH(I)   ,0,:),     EwHXrQC )
      TmpVec2 = CROSS_PRODUCT( RtHSdat%PAngVelEH(p%PH(I)   ,1,:), RtHSdat%rQC )

      RtHSdat%PLinVelEC(p%PH(I),0,:) = TmpVec0    +                         RtHSdat%PLinVelEC(p%PH(I)   ,0,:)
      RtHSdat%PLinVelEC(p%PH(I),1,:) = TmpVec1    + TmpVec2 +               RtHSdat%PLinVelEC(p%PH(I)   ,1,:)

      RtHSdat%LinAccECt              =  RtHSdat%LinAccECt + x%QDT(p%PH(I) )*RtHSdat%PLinVelEC(p%PH(I)   ,1,:)

   ENDDO          ! I - all DOFs associated with the angular motion of the hub (body H)




   DO K = 1,p%NumBl ! Loop through all blades

      DO J = 0,p%TipNode ! Loop through the blade nodes / elements

      ! Define the partial linear velocities (and their 1st derivatives) of the
      !   current node (point S(RNodes(J))) in the inertia frame.  Also define
      !   the overall linear velocity of the current node in the inertia frame.
      !   Also, define the portion of the linear acceleration of the current node
      !   in the inertia frame associated with everything but the QD2T()'s:

         EwHXrQS = CROSS_PRODUCT(  RtHSdat%AngVelEH, RtHSdat%rQS(:,K,J) )

         RtHSdat%PLinVelES(K,J,          :,:,:) = RtHSdat%PLinVelEQ(:,:,:)
         RtHSdat%PLinVelES(K,J,DOF_BF(K,1),0,:) = p%TwistedSF(K,1,1,J,0)                          *CoordSys%j1(K,:) &  !bjj: this line can be optimized
                                                + p%TwistedSF(K,2,1,J,0)                          *CoordSys%j2(K,:) &
                                                - (   p%AxRedBld(K,1,1,J)*x%QT ( DOF_BF(K,1) ) &
                                                    + p%AxRedBld(K,1,2,J)*x%QT ( DOF_BF(K,2) ) &
                                                    + p%AxRedBld(K,1,3,J)*x%QT ( DOF_BE(K,1) )   )*CoordSys%j3(K,:)
         RtHSdat%PLinVelES(K,J,DOF_BE(K,1),0,:) = p%TwistedSF(K,1,3,J,0)                          *CoordSys%j1(K,:) &
                                                + p%TwistedSF(K,2,3,J,0)                          *CoordSys%j2(K,:) &
                                                - (   p%AxRedBld(K,3,3,J)*x%QT ( DOF_BE(K,1) ) &
                                                    + p%AxRedBld(K,2,3,J)*x%QT ( DOF_BF(K,2) ) &
                                                    + p%AxRedBld(K,1,3,J)*x%QT ( DOF_BF(K,1) )   )*CoordSys%j3(K,:)
         RtHSdat%PLinVelES(K,J,DOF_BF(K,2),0,:) = p%TwistedSF(K,1,2,J,0)                          *CoordSys%j1(K,:) &
                                                + p%TwistedSF(K,2,2,J,0)                          *CoordSys%j2(K,:) &
                                                - (   p%AxRedBld(K,2,2,J)*x%QT ( DOF_BF(K,2) ) &
                                                    + p%AxRedBld(K,1,2,J)*x%QT ( DOF_BF(K,1) ) &
                                                    + p%AxRedBld(K,2,3,J)*x%QT ( DOF_BE(K,1) )   )*CoordSys%j3(K,:)

         TmpVec1 = CROSS_PRODUCT( RtHSdat%AngVelEH, RtHSdat%PLinVelES(K,J,DOF_BF(K,1),0,:) )
         TmpVec2 = CROSS_PRODUCT( RtHSdat%AngVelEH, RtHSdat%PLinVelES(K,J,DOF_BE(K,1),0,:) )
         TmpVec3 = CROSS_PRODUCT( RtHSdat%AngVelEH, RtHSdat%PLinVelES(K,J,DOF_BF(K,2),0,:) )

         RtHSdat%PLinVelES(K,J,DOF_BF(K,1),1,:) = TmpVec1 &
                                                - (   p%AxRedBld(K,1,1,J)*x%QDT( DOF_BF(K,1) ) &
                                                    + p%AxRedBld(K,1,2,J)*x%QDT( DOF_BF(K,2) ) &
                                                    + p%AxRedBld(K,1,3,J)*x%QDT( DOF_BE(K,1) )   )*CoordSys%j3(K,:)
         RtHSdat%PLinVelES(K,J,DOF_BE(K,1),1,:) = TmpVec2 &
                                                - (   p%AxRedBld(K,3,3,J)*x%QDT( DOF_BE(K,1) ) &
                                                    + p%AxRedBld(K,2,3,J)*x%QDT( DOF_BF(K,2) ) &
                                                    + p%AxRedBld(K,1,3,J)*x%QDT( DOF_BF(K,1) )   )*CoordSys%j3(K,:)
         RtHSdat%PLinVelES(K,J,DOF_BF(K,2),1,:) = TmpVec3 &
                                                - (   p%AxRedBld(K,2,2,J)*x%QDT( DOF_BF(K,2) ) &
                                                    + p%AxRedBld(K,1,2,J)*x%QDT( DOF_BF(K,1) ) &
                                                    + p%AxRedBld(K,2,3,J)*x%QDT( DOF_BE(K,1) )   )*CoordSys%j3(K,:)

         LinVelHS                 = x%QDT( DOF_BF(K,1) )*RtHSdat%PLinVelES(K,J,DOF_BF(K,1),0,:) &
                                  + x%QDT( DOF_BE(K,1) )*RtHSdat%PLinVelES(K,J,DOF_BE(K,1),0,:) &
                                  + x%QDT( DOF_BF(K,2) )*RtHSdat%PLinVelES(K,J,DOF_BF(K,2),0,:)
         RtHSdat%LinAccESt(:,K,J) = x%QDT( DOF_BF(K,1) )*RtHSdat%PLinVelES(K,J,DOF_BF(K,1),1,:) &
                                  + x%QDT( DOF_BE(K,1) )*RtHSdat%PLinVelES(K,J,DOF_BE(K,1),1,:) &
                                  + x%QDT( DOF_BF(K,2) )*RtHSdat%PLinVelES(K,J,DOF_BF(K,2),1,:)

         RtHSdat%LinVelES(:,J,K)  = LinVelHS + RtHSdat%LinVelEZ
         DO I = 1,p%NPH   ! Loop through all DOFs associated with the angular motion of the hub (body H)

            TmpVec0 = CROSS_PRODUCT(   RtHSdat%PAngVelEH(p%PH(I),0,:), RtHSdat%rQS(:,K,J)            )  !bjj: this line can be optimized
            TmpVec1 = CROSS_PRODUCT(   RtHSdat%PAngVelEH(p%PH(I),0,:),     EwHXrQS        + LinVelHS )  !bjj: this line can be optimized
            TmpVec2 = CROSS_PRODUCT(   RtHSdat%PAngVelEH(p%PH(I),1,:), RtHSdat%rQS(:,K,J)            )  !bjj: this line can be optimized

            RtHSdat%PLinVelES(K,J,p%PH(I),0,:) = RtHSdat%PLinVelES(K,J,p%PH(I),0,:) + TmpVec0            !bjj: this line can be optimized
            RtHSdat%PLinVelES(K,J,p%PH(I),1,:) = RtHSdat%PLinVelES(K,J,p%PH(I),1,:) + TmpVec1 + TmpVec2  !bjj: this line can be optimized

            RtHSdat%LinVelES(:,J,K)          = RtHSdat%LinVelES(:,J,K)   + x%QDT(p%PH(I))*RtHSdat%PLinVelES(K,J,p%PH(I),0,:)  !bjj: this line can be optimized
            RtHSdat%LinAccESt(:,K,J)         = RtHSdat%LinAccESt(:,K,J)  + x%QDT(p%PH(I))*RtHSdat%PLinVelES(K,J,p%PH(I),1,:)  !bjj: this line can be optimized

         END DO ! I - all DOFs associated with the angular motion of the hub (body H)

      END DO !J = 0,p%TipNodes ! Loop through the blade nodes / elements
      
      
   !JASON: USE TipNode HERE INSTEAD OF BldNodes IF YOU ALLOCATE AND DEFINE n1, n2, n3, m1, m2, AND m3 TO USE TipNode.  THIS WILL REQUIRE THAT THE AERODYNAMIC AND STRUCTURAL TWISTS, AeroTwst() AND ThetaS(), BE KNOWN AT THE TIP!!!
      !IF (.NOT. p%BD4Blades) THEN
      !   RtHSdat%LinVelESm2(K) = DOT_PRODUCT( RtHSdat%LinVelES(:,p%TipNode,K), CoordSys%m2(K,p%BldNodes,:) )
      !END IF
            
   END DO !K = 1,p%NumBl


   RtHSdat%PLinVelEW(       :,:,:) = RtHSdat%PLinVelEO(:,:,:)
   DO I = 1,NPN   ! Loop through all DOFs associated with the angular motion of the nacelle (body N)

      TmpVec0 = CROSS_PRODUCT( RtHSdat%PAngVelEN(PN(I)   ,0,:), RtHSdat%rOW                 )
      TmpVec1 = CROSS_PRODUCT( RtHSdat%PAngVelEN(PN(I)   ,0,:),     EwNXrOW                 )
      TmpVec2 = CROSS_PRODUCT( RtHSdat%PAngVelEN(PN(I)   ,1,:), RtHSdat%rOW                 )

      RtHSdat%PLinVelEW(PN(I),0,:) = TmpVec0    +               RtHSdat%PLinVelEW(PN(I)   ,0,:)
      RtHSdat%PLinVelEW(PN(I),1,:) = TmpVec1    + TmpVec2 +     RtHSdat%PLinVelEW(PN(I)   ,1,:)

       LinAccEWt                   =  LinAccEWt + x%QDT(PN(I) )*RtHSdat%PLinVelEW(PN(I)   ,1,:)

   ENDDO          ! I - all DOFs associated with the angular motion of the nacelle (body N)


   ! Velocities of point I (tail boom center of mass)
   RtHSdat%PLinVelEI(       :,:,:) = RtHSdat%PLinVelEW(:,:,:)
   DO I = 1,NPA   ! Loop through all DOFs associated with the angular motion of the tail (body A)

      TmpVec0 = CROSS_PRODUCT( RtHSdat%PAngVelEA(PA(I)   ,0,:), RtHSdat%rWI                 )
      TmpVec1 = CROSS_PRODUCT( RtHSdat%PAngVelEA(PA(I)   ,0,:),     EwAXrWI                 )
      TmpVec2 = CROSS_PRODUCT( RtHSdat%PAngVelEA(PA(I)   ,1,:), RtHSdat%rWI                 )

      RtHSdat%PLinVelEI(PA(I),0,:) = TmpVec0    +                       RtHSdat%PLinVelEI(PA(I)   ,0,:)
      RtHSdat%PLinVelEI(PA(I),1,:) = TmpVec1    + TmpVec2 +             RtHSdat%PLinVelEI(PA(I)   ,1,:)

      RtHSdat%LinAccEIt            =  RtHSdat%LinAccEIt + x%QDT(PA(I) )*RtHSdat%PLinVelEI(PA(I)   ,1,:)

   ENDDO          ! I - all DOFs associated with the angular motion of the tail (body A)


   ! Velocities of point J (tail fin center of mass)
   RtHSdat%PLinVelEJ(       :,:,:) = RtHSdat%PLinVelEW(:,:,:)
   RtHSdat%LinVelEJ                = RtHSdat%LinVelEZ
   DO I = 1,NPA   ! Loop through all DOFs associated with the angular motion of the tail (body A)

      TmpVec0 = CROSS_PRODUCT( RtHSdat%PAngVelEA(PA(I)   ,0,:), RtHSdat%rWJ                 )
      TmpVec1 = CROSS_PRODUCT( RtHSdat%PAngVelEA(PA(I)   ,0,:),     EwAXrWJ                 )
      TmpVec2 = CROSS_PRODUCT( RtHSdat%PAngVelEA(PA(I)   ,1,:), RtHSdat%rWJ                 )

      RtHSdat%PLinVelEJ(PA(I),0,:) = TmpVec0    +               RtHSdat%PLinVelEJ(PA(I)   ,0,:)
      RtHSdat%PLinVelEJ(PA(I),1,:) = TmpVec1    + TmpVec2 +     RtHSdat%PLinVelEJ(PA(I)   ,1,:)

       RtHSdat%LinVelEJ            =  RtHSdat%LinVelEJ  + x%QDT(PA(I) )*RtHSdat%PLinVelEJ(PA(I)   ,0,:)
       RtHSdat%LinAccEJt           =  RtHSdat%LinAccEJt + x%QDT(PA(I) )*RtHSdat%PLinVelEJ(PA(I)   ,1,:)

   ENDDO          ! I - all DOFs associated with the angular motion of the tail (body A)


   DO J = 0,p%TwrNodes  ! Loop through the tower nodes / elements


      ! Define the partial linear velocities (and their 1st derivatives) of the current node (point T(HNodes(J))) in the inertia frame.
      !  Also define the overall linear velocity of the current node in the inertia frame.
      !  Also, define the portion of the linear acceleration of the current node in the inertia frame associated with
      !    everything but the QD2T()'s:

      EwXXrZT                   = CROSS_PRODUCT(  RtHSdat%AngVelEX, RtHSdat%rZT(:,J) )

      RtHSdat%PLinVelET(J,       :,:,:) = RtHSdat%PLinVelEZ(:,:,:)  !bjj: can this line be optimized
      RtHSdat%PLinVelET(J,DOF_TFA1,0,:) = p%TwrFASF(1,J,0)*CoordSys%a1 - (   p%AxRedTFA(1,1,J)* x%QT(DOF_TFA1) &
                                                                           + p%AxRedTFA(1,2,J)* x%QT(DOF_TFA2)   )*CoordSys%a2  
      RtHSdat%PLinVelET(J,DOF_TSS1,0,:) = p%TwrSSSF(1,J,0)*CoordSys%a3 - (   p%AxRedTSS(1,1,J)* x%QT(DOF_TSS1) &
                                                                           + p%AxRedTSS(1,2,J)* x%QT(DOF_TSS2)   )*CoordSys%a2
      RtHSdat%PLinVelET(J,DOF_TFA2,0,:) = p%TwrFASF(2,J,0)*CoordSys%a1 - (   p%AxRedTFA(2,2,J)* x%QT(DOF_TFA2) &
                                                                           + p%AxRedTFA(1,2,J)* x%QT(DOF_TFA1)   )*CoordSys%a2
      RtHSdat%PLinVelET(J,DOF_TSS2,0,:) = p%TwrSSSF(2,J,0)*CoordSys%a3 - (   p%AxRedTSS(2,2,J)* x%QT(DOF_TSS2) &
                                                                           + p%AxRedTSS(1,2,J)* x%QT(DOF_TSS1)   )*CoordSys%a2

      TmpVec1 = CROSS_PRODUCT( RtHSdat%AngVelEX, RtHSdat%PLinVelET(J,DOF_TFA1,0,:) )
      TmpVec2 = CROSS_PRODUCT( RtHSdat%AngVelEX, RtHSdat%PLinVelET(J,DOF_TSS1,0,:) )
      TmpVec3 = CROSS_PRODUCT( RtHSdat%AngVelEX, RtHSdat%PLinVelET(J,DOF_TFA2,0,:) )
      TmpVec4 = CROSS_PRODUCT( RtHSdat%AngVelEX, RtHSdat%PLinVelET(J,DOF_TSS2,0,:) )

      RtHSdat%PLinVelET(J,DOF_TFA1,1,:) = TmpVec1 - (   p%AxRedTFA(1,1,J)*x%QDT(DOF_TFA1) &
                                                      + p%AxRedTFA(1,2,J)*x%QDT(DOF_TFA2)   )*CoordSys%a2
      RtHSdat%PLinVelET(J,DOF_TSS1,1,:) = TmpVec2 - (   p%AxRedTSS(1,1,J)*x%QDT(DOF_TSS1) &
                                                      + p%AxRedTSS(1,2,J)*x%QDT(DOF_TSS2)   )*CoordSys%a2
      RtHSdat%PLinVelET(J,DOF_TFA2,1,:) = TmpVec3 - (   p%AxRedTFA(2,2,J)*x%QDT(DOF_TFA2) &
                                                      + p%AxRedTFA(1,2,J)*x%QDT(DOF_TFA1)   )*CoordSys%a2
      RtHSdat%PLinVelET(J,DOF_TSS2,1,:) = TmpVec4 - (   p%AxRedTSS(2,2,J)*x%QDT(DOF_TSS2) &
                                                      + p%AxRedTSS(1,2,J)*x%QDT(DOF_TSS1)   )*CoordSys%a2

              LinVelXT       = x%QDT(DOF_TFA1)*RtHSdat%PLinVelET(J,DOF_TFA1,0,:) &
                             + x%QDT(DOF_TSS1)*RtHSdat%PLinVelET(J,DOF_TSS1,0,:) &
                             + x%QDT(DOF_TFA2)*RtHSdat%PLinVelET(J,DOF_TFA2,0,:) &
                             + x%QDT(DOF_TSS2)*RtHSdat%PLinVelET(J,DOF_TSS2,0,:)
      RtHSdat%LinAccETt(:,J) = x%QDT(DOF_TFA1)*RtHSdat%PLinVelET(J,DOF_TFA1,1,:) &
                             + x%QDT(DOF_TSS1)*RtHSdat%PLinVelET(J,DOF_TSS1,1,:) &
                             + x%QDT(DOF_TFA2)*RtHSdat%PLinVelET(J,DOF_TFA2,1,:) &
                             + x%QDT(DOF_TSS2)*RtHSdat%PLinVelET(J,DOF_TSS2,1,:)

      RtHSdat%LinVelET(:,J)  = LinVelXT + RtHSdat%LinVelEZ
      DO I = 1,NPX   ! Loop through all DOFs associated with the angular motion of the platform (body X)

         TmpVec0   = CROSS_PRODUCT( RtHSdat%PAngVelEX(PX(I),0,:), RtHSdat%rZT(:,J)            )
         TmpVec1   = CROSS_PRODUCT( RtHSdat%PAngVelEX(PX(I),0,:), EwXXrZT      + LinVelXT )

         RtHSdat%PLinVelET(J,PX(I),0,:) = RtHSdat%PLinVelET(J,PX(I),0,:) + TmpVec0
         RtHSdat%PLinVelET(J,PX(I),1,:) = RtHSdat%PLinVelET(J,PX(I),1,:) + TmpVec1

         RtHSdat%LinVelET( :,        J) = RtHSdat%LinVelET( :,        J) + x%QDT(PX(I))*RtHSdat%PLinVelET(J,PX(I),0,:)
         RtHSdat%LinAccETt(:,        J) = RtHSdat%LinAccETt(:,        J) + x%QDT(PX(I))*RtHSdat%PLinVelET(J,PX(I),1,:)

      ENDDO          ! I - all DOFs associated with the angular motion of the platform (body X)


   END DO ! J


END SUBROUTINE CalculateLinearVelPAcc
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is used to calculate the forces and moments stored in other states that are used in
!! both the CalcOutput and CalcContStateDeriv routines.
SUBROUTINE CalculateForcesMoments( p, x, CoordSys, u, RtHSdat )
!..................................................................................................................................

      ! Passed variables
   TYPE(ED_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(ED_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
   TYPE(ED_CoordSys),            INTENT(IN   )  :: CoordSys    !< The coordinate systems that have been set for these states/time
   TYPE(ED_InputType),           INTENT(IN   )  :: u           !< The aero (blade) & nacelle forces/moments
   TYPE(ED_RtHndSide),           INTENT(INOUT)  :: RtHSdat     !< data from the RtHndSid module (contains positions to be set)

      ! Local variables
   REAL(ReKi)                   :: TmpVec    (3)                                   ! A temporary vector used in various computations.
   REAL(ReKi)                   :: TmpVec1   (3)                                   ! A temporary vector used in various computations.
   REAL(ReKi)                   :: TmpVec2   (3)                                   ! A temporary vector used in various computations.
   REAL(ReKi)                   :: TmpVec3   (3)                                   ! A temporary vector used in various computations.
   REAL(ReKi)                   :: TmpVec4   (3)                                   ! A temporary vector used in various computations.
   REAL(ReKi)                   :: TmpVec5   (3)                                   ! A temporary vector used in various computations.
   REAL(ReKi)                   :: Force(3)  ! External force  (e.g. from AeroDyn)
   REAL(ReKi)                   :: Moment(3) ! External moment (e.g. from AeroDyn)
   INTEGER(IntKi)               :: I                                               ! Loops through some or all of the DOFs
   INTEGER(IntKi)               :: J                                               ! Counter for elements
   INTEGER(IntKi)               :: K                                               ! Counter for blades
   INTEGER(IntKi)               :: NodeNum                                         ! Node number for blade element (on a single mesh)
      
!.....................................
! Compute forces and moments from properties related to Aero inputs   
!FSTipDrag
!FSAero and MMAero
!.....................................
   DO K = 1,p%NumBl ! Loop through all blades
      
         ! Calculate the tip drag forces if necessary:
      !bjj: add this back when we've figured out how to handle the tip brakes:
      !RtHSdat%FSTipDrag(:,K) = m%CoordSys%m2(K,p%BldNodes,:)*SIGN( 0.5*p%AirDens*(RtHSdat%LinVelESm2(K)**2)*u%TBDrCon(K), -1.*RtHSdat%LinVelESm2(K) )
      RtHSdat%FSTipDrag = 0.0_ReKi         ! Calculate the tip drag forces if necessary

      
   ! Calculate the normal and tangential aerodynamic forces and the aerodynamic
   !   pitching moment at the current element per unit span by calling AeroDyn,
   !   if necessary:
      
      DO J = 1,p%BldNodes ! Loop through the blade nodes / elements


   ! Calculate the aerodynamic pitching moment arm (i.e., the position vector
   !   from point S on the blade to the aerodynamic center of the element):

         RtHSdat%rSAerCen(:,J,K) = p%rSAerCenn1(K,J)*CoordSys%n1(K,J,:) + p%rSAerCenn2(K,J)*CoordSys%n2(K,J,:)   

!        rPAerCen     = m%RtHS%rPQ + m%RtHS%rQS(:,K,J) + m%RtHS%rSAerCen(:,J,K)     ! Position vector from teeter pin (point P)  to blade analysis node aerodynamic center.
!        rAerCen      =                       m%RtHS%rS (:,K,J) + m%RtHS%rSAerCen(:,J,K)     ! Position vector from inertial frame origin to blade analysis node aerodynamic center.
         

   ! fill FSAero() and MMAero() with the forces resulting from inputs u%BladeLn2Mesh(K)%Force(1:2,:) and u%BladeLn2Mesh(K)%Moment(3,:):
   ! [except, we're ignoring the additional nodes we added on the mesh end points]
   
         NodeNum = J ! we're ignoring the root and tip
         
         if (p%UseAD14) then
            RtHSdat%FSAero(:,K,J) = ( u%BladePtLoads(K)%Force(1,NodeNum) * CoordSys%te1(K,J,:) &
                                    + u%BladePtLoads(K)%Force(2,NodeNum) * CoordSys%te2(K,J,:) ) / p%DRNodes(J)

            RtHSdat%MMAero(:,K,J) = CROSS_PRODUCT( RtHSdat%rSAerCen(:,J,K), RtHSdat%FSAero(:,K,J) )&
                                  + u%BladePtLoads(K)%Moment(3,NodeNum)/p%DRNodes(J) * CoordSys%te3(K,J,:)        
         else
            RtHSdat%FSAero(1,K,J) =  u%BladePtLoads(K)%Force(1,NodeNum) / p%DRNodes(J)
            RtHSdat%FSAero(2,K,J) =  u%BladePtLoads(K)%Force(3,NodeNum) / p%DRNodes(J) 
            RtHSdat%FSAero(3,K,J) = -u%BladePtLoads(K)%Force(2,NodeNum) / p%DRNodes(J)

            RtHSdat%MMAero(1,K,J) =  u%BladePtLoads(K)%Moment(1,NodeNum) / p%DRNodes(J)
            RtHSdat%MMAero(2,K,J) =  u%BladePtLoads(K)%Moment(3,NodeNum) / p%DRNodes(J)
            RtHSdat%MMAero(3,K,J) = -u%BladePtLoads(K)%Moment(2,NodeNum) / p%DRNodes(J)
         end if
                     
         
      END DO !J
   END DO  ! K 
   
   
!.....................................
! PFrcS0B and PMomH0B  
!.....................................
DO K = 1,p%NumBl ! Loop through all blades

      ! Initialize the partial forces and moments (including those associated
      !   with the QD2T()'s and those that are not) at the blade root (point S(0))
      !   using the tip brake effects:

      RtHSdat%PFrcS0B(:,K,:) = 0.0 ! Initialize these partial
      RtHSdat%PMomH0B(:,K,:) = 0.0 ! forces and moments to zero
      DO I = 1,p%DOFs%NPSE(K)  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of blade K

         TmpVec1 = -p%TipMass(K)*RtHSdat%PLinVelES(K,p%TipNode,p%DOFs%PSE(K,I),0,:)                            ! The portion of PFrcS0B associated with the tip brake

         RtHSdat%PFrcS0B(:,K,p%DOFs%PSE(K,I)) = TmpVec1
         RtHSdat%PMomH0B(:,K,p%DOFs%PSE(K,I)) = CROSS_PRODUCT( RtHSdat%rS0S(:,K,p%TipNode), TmpVec1 )          ! The portion of PMomH0B associated with the tip brake

      ENDDO             ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of blade K  
   
   
      DO J = 1,p%BldNodes ! Loop through the blade nodes / elements

      ! Integrate to find the partial forces and moments (including those associated
      !   with the QD2T()'s and those that are not) at the blade root (point S(0)):

         DO I = 1,p%DOFs%NPSE(K)  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of blade K

            TmpVec1 = -p%BElmntMass(J,K)*RtHSdat%PLinVelES(K,J,p%DOFs%PSE(K,I),0,:)   ! The portion of PFrcS0B associated with blade element J

            RtHSdat%PFrcS0B(:,K,p%DOFs%PSE(K,I)) = RtHSdat%PFrcS0B(:,K,p%DOFs%PSE(K,I)) + TmpVec1
            RtHSdat%PMomH0B(:,K,p%DOFs%PSE(K,I)) = RtHSdat%PMomH0B(:,K,p%DOFs%PSE(K,I)) + &
                                                    CROSS_PRODUCT( RtHSdat%rS0S(:,K,J), TmpVec1 )                   ! The portion of PMomH0B associated with blade element J

         ENDDO             ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of blade K
      END DO
      
      
   END DO     
   
 
!.....................................
! FrcS0Bt and MomH0Bt
!.....................................
   DO K = 1,p%NumBl ! Loop through all blades
   
      TmpVec1 = RtHSdat%FSTipDrag(:,K) - p%TipMass(K)*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccESt(:,K,p%TipNode) ) ! The portion of FrcS0Bt associated with the tip brake
      RtHSdat%FrcS0Bt(:,K) = TmpVec1
      RtHSdat%MomH0Bt(:,K) = CROSS_PRODUCT(  RtHSdat%rS0S(:,K,p%TipNode), TmpVec1 )                                 ! The portion of MomH0Bt associated with the tip brake

      DO J = 1,p%BldNodes ! Loop through the blade nodes / elements      
      
         TmpVec1 = RtHSdat%FSAero(:,K,J)*p%DRNodes(J) - p%BElmntMass(J,K)*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccESt(:,K,J) ) ! The portion of FrcS0Bt associated with blade element J
         TmpVec2 = CROSS_PRODUCT( RtHSdat%rS0S(:,K,J), TmpVec1 )                                    ! The portion of MomH0Bt associated with blade element J
         TmpVec3 = RtHSdat%MMAero(:,K,J)*p%DRNodes(J)                                               ! The total external moment applied to blade element J

         RtHSdat%FrcS0Bt(:,K) = RtHSdat%FrcS0Bt(:,K) + TmpVec1
         RtHSdat%MomH0Bt(:,K) = RtHSdat%MomH0Bt(:,K) + TmpVec2 + TmpVec3
      
      END DO !J
      
   END DO !K   
         
      
!.....................................
! PFrcPRot AND PMomLPRot:  
!   ( requires PFrcS0B and PMomH0B)
!.....................................
      ! Initialize the partial forces and moments (including those associated
      !   with the QD2T()'s and those that are not) at the teeter pin (point P) using the hub mass effects:    

   RtHSdat%PFrcPRot  = 0.0   ! Initialize these partial
   RtHSdat%PMomLPRot = 0.0   ! forces and moments to zero
   DO I = 1,p%DOFs%NPCE  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the hub center of mass (point C)

      TmpVec1 = -p%HubMass*RtHSdat%PLinVelEC(p%DOFs%PCE(I),0,:)     ! The portion of PFrcPRot  associated with the HubMass
      TmpVec2 = CROSS_PRODUCT( RtHSdat%rPC, TmpVec1 )      ! The portion of PMomLPRot associated with the HubMass

      RtHSdat%PFrcPRot (:,p%DOFs%PCE(I)) = TmpVec1
      RtHSdat%PMomLPRot(:,p%DOFs%PCE(I)) = TmpVec2 - p%Hubg1Iner*CoordSys%g1*DOT_PRODUCT( CoordSys%g1, RtHSdat%PAngVelEH(p%DOFs%PCE(I),0,:) ) &
                                                   - p%Hubg2Iner*CoordSys%g2*DOT_PRODUCT( CoordSys%g2, RtHSdat%PAngVelEH(p%DOFs%PCE(I),0,:) )

   ENDDO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the hub center of mass (point C)

   
   DO K = 1,p%NumBl ! Loop through all blades
   
         ! Calculate the position vector from the teeter pin to the blade root:
   
      !rPS0(:,K) = RtHSdat%rPQ + p%HubRad*CoordSys%j3(K,:)   ! Position vector from teeter pin (point P) to blade root (point S(0)).
            
      ! Add the blade effects to the partial forces and moments (including those associated with the QD2T()'s and those that are 
      !   not) at the teeter pin (point P):

      DO I = 1,p%DOFs%NPSE(K)  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of blade K

         TmpVec = CROSS_PRODUCT( RtHSdat%rPS0(:,K), RtHSdat%PFrcS0B(:,K,p%DOFs%PSE(K,I)) ) ! The portion of PMomLPRot associated with PFrcS0B.

         RtHSdat%PFrcPRot (:,p%DOFs%PSE(K,I)) = RtHSdat%PFrcPRot (:,p%DOFs%PSE(K,I)) + RtHSdat%PFrcS0B(:,K,p%DOFs%PSE(K,I))
         RtHSdat%PMomLPRot(:,p%DOFs%PSE(K,I)) = RtHSdat%PMomLPRot(:,p%DOFs%PSE(K,I)) + RtHSdat%PMomH0B(:,K,p%DOFs%PSE(K,I))+TmpVec

      ENDDO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of blade K


   END DO   ! K   

!.....................................
! FrcPRott and MomLPRott:
!   (requires FrcS0Bt and MomH0Bt)
!.....................................

   TmpVec1 = -p%HubMass*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccECt )                     ! The portion of FrcPRott  associated with the HubMass
   TmpVec2 = CROSS_PRODUCT( RtHSdat%rPC, TmpVec1 )                                        ! The portion of MomLPRott associated with the HubMass
   TmpVec  = p%Hubg1Iner*CoordSys%g1*DOT_PRODUCT( CoordSys%g1, RtHSdat%AngVelEH ) &       ! = ( Hub inertia dyadic ) dot ( angular velocity of hub in the inertia frame )
           + p%Hubg2Iner*CoordSys%g2*DOT_PRODUCT( CoordSys%g2, RtHSdat%AngVelEH )
   TmpVec3 = CROSS_PRODUCT( -RtHSdat%AngVelEH, TmpVec )                                   ! = ( -angular velocity of hub in the inertia frame ) cross ( TmpVec )

   RtHSdat%FrcPRott(1)  = TmpVec1(1) + u%HubPtLoad%Force(1,1)
   RtHSdat%FrcPRott(2)  = TmpVec1(2) + u%HubPtLoad%Force(3,1)
   RtHSdat%FrcPRott(3)  = TmpVec1(3) - u%HubPtLoad%Force(2,1)
   
   RtHSdat%MomLPRott    = TmpVec2 + TmpVec3 - p%Hubg1Iner*CoordSys%g1*DOT_PRODUCT( CoordSys%g1, RtHSdat%AngAccEHt ) &
                                            - p%Hubg2Iner*CoordSys%g2*DOT_PRODUCT( CoordSys%g2, RtHSdat%AngAccEHt )                                          
      
   RtHSdat%MomLPRott(1) = RtHSdat%MomLPRott(1) + u%HubPtLoad%Moment(1,1)
   RtHSdat%MomLPRott(2) = RtHSdat%MomLPRott(2) + u%HubPtLoad%Moment(3,1)
   RtHSdat%MomLPRott(3) = RtHSdat%MomLPRott(3) - u%HubPtLoad%Moment(2,1)
   
   DO K = 1,p%NumBl ! Loop through all blades
   
         ! Calculate the position vector from the teeter pin to the blade root:
   
      !rPS0 = RtHSdat%rPQ + p%HubRad*m%CoordSys%j3(K,:)   ! Position vector from teeter pin (point P) to blade root (point S(0)).
            
      TmpVec = CROSS_PRODUCT( RtHSdat%rPS0(:,K), RtHSdat%FrcS0Bt(:,K) )       ! The portion of MomLPRott associated with FrcS0Bt.

      RtHSdat%FrcPRott  = RtHSdat%FrcPRott  + RtHSdat%FrcS0Bt(:,K)
      RtHSdat%MomLPRott = RtHSdat%MomLPRott + RtHSdat%MomH0Bt(:,K) + TmpVec

   END DO   ! K

   
   ! Define the partial forces and moments (including those associated with
   !   the QD2T()'s and those that are not) at the specified point on the
   !   rotor-furl axis (point V) / nacelle (body N) using the structure that
   !   furls with the rotor, generator, and rotor effects.
   
!.....................................
! PMomNGnRt and PFrcVGnRt
!  (requires PMomLPRot and PFrcPRot)
!..................................... 
   RtHSdat%PFrcVGnRt = RtHSdat%PFrcPRot    ! Initialize these partial forces and
   RtHSdat%PMomNGnRt = RtHSdat%PMomLPRot   ! moments using the rotor effects
   DO I = 1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs

      TmpVec = CROSS_PRODUCT( RtHSdat%rVP, RtHSdat%PFrcPRot(:,p%DOFs%SrtPS(I)) )  ! The portion of PMomNGnRt associated with the PFrcPRot

      RtHSdat%PMomNGnRt(:,p%DOFs%SrtPS(I)) = RtHSdat%PMomNGnRt(:,p%DOFs%SrtPS(I)) + TmpVec

   ENDDO             ! I - All active (enabled) DOFs
   DO I = 1,p%DOFs%NPDE  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the center of mass of the structure that furls with the rotor (not including rotor) (point D)

      TmpVec1 = -p%RFrlMass*RtHSdat%PLinVelED(p%DOFs%PDE(I)  ,0,:)           ! The portion of PFrcVGnRt associated with the RFrlMass
      TmpVec2 = CROSS_PRODUCT( RtHSdat%rVD,              TmpVec1 )  ! The portion of PMomNGnRt associated with the RFrlMass

      RtHSdat%PFrcVGnRt(:,p%DOFs%PDE(I)  ) = RtHSdat%PFrcVGnRt(:,p%DOFs%PDE(I)  ) + TmpVec1

      RtHSdat%PMomNGnRt(:,p%DOFs%PDE(I)  ) = RtHSdat%PMomNGnRt(:,p%DOFs%PDE(I)  ) + TmpVec2                                   &
                                           - p%RrfaIner*CoordSys%rfa*DOT_PRODUCT( CoordSys%rfa, RtHSdat%PAngVelER(p%DOFs%PDE(I) ,0,:) ) &
                                           - p%GenIner*CoordSys%c1 *DOT_PRODUCT(  CoordSys%c1 , RtHSdat%PAngVelEG(p%DOFs%PDE(I) ,0,:) )

   ENDDO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the center of mass of the structure that furls with the rotor (not including rotor) (point D)
   IF ( p%DOF_Flag(DOF_GeAz) )  THEN

      RtHSdat%PMomNGnRt(:,DOF_GeAz) = RtHSdat%PMomNGnRt(:,DOF_GeAz)                                             &     ! The previous loop (DO I = 1,NPDE) misses the DOF_GeAz-contribution to: ( Generator inertia dyadic ) dot ( partial angular velocity of the generator in the inertia frame )
                            -  p%GenIner*CoordSys%c1 *DOT_PRODUCT( CoordSys%c1, RtHSdat%PAngVelEG(DOF_GeAz,0,:) )     ! Thus, add this contribution if necessary.

   ENDIF   
   
!.....................................
! FrcVGnRtt and MomNGnRtt
!  (requires FrcPRott and MomLPRott)
!.....................................
   TmpVec1 = -p%RFrlMass*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccEDt )                ! The portion of FrcVGnRtt associated with the RFrlMass
   TmpVec2 = CROSS_PRODUCT( RtHSdat%rVD      ,  TmpVec1 )                             ! The portion of MomNGnRtt associated with the RFrlMass
   TmpVec3 = CROSS_PRODUCT( RtHSdat%rVP      , RtHSdat%FrcPRott )                     ! The portion of MomNGnRtt associated with the FrcPRott
   TmpVec  = p%RrfaIner*CoordSys%rfa*DOT_PRODUCT( CoordSys%rfa, RtHSdat%AngVelER )    ! = ( R inertia dyadic ) dot ( angular velocity of structure that furls with the rotor in the inertia frame )
   TmpVec4 = CROSS_PRODUCT( -RtHSdat%AngVelER, TmpVec )                               ! = ( -angular velocity of structure that furls with the rotor in the inertia frame ) cross ( TmpVec )
   TmpVec  =  p%GenIner*CoordSys%c1* DOT_PRODUCT( CoordSys%c1 , RtHSdat%AngVelEG )    ! = ( Generator inertia dyadic ) dot ( angular velocity of generator in the inertia frame )
   TmpVec5 = CROSS_PRODUCT( -RtHSdat%AngVelEG, TmpVec )                               ! = ( -angular velocity of generator in the inertia frame ) cross ( TmpVec )

   RtHSdat%FrcVGnRtt = RtHSdat%FrcPRott  + TmpVec1
   RtHSdat%MomNGnRtt = RtHSdat%MomLPRott + TmpVec2 + TmpVec3 + TmpVec4 + TmpVec5            &
                     - p%RrfaIner*CoordSys%rfa*DOT_PRODUCT( CoordSys%rfa, RtHSdat%AngAccERt ) &
                     -  p%GenIner*CoordSys%c1 *DOT_PRODUCT( CoordSys%c1 , RtHSdat%AngAccEGt )


!.....................................
! PFrcWTail and PMomNTail
!.....................................
      ! Define the partial forces and moments (including those associated with the QD2T()'s and 
      !   those that are not) at the specified point on the tail-furl axis (point W) / nacelle (body N) using the tail effects.

   RtHSdat%PFrcWTail = 0.0   ! Initialize these partial
   RtHSdat%PMomNTail = 0.0   ! forces and moments to zero
   DO I = 1,p%DOFs%NPIE  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the tail boom center of mass (point I)

      TmpVec1 = -p%BoomMass*RtHSdat%PLinVelEI(p%DOFs%PIE(I),0,:)    ! The portion of PFrcWTail associated with the BoomMass
      TmpVec2 = -p%TFinMass*RtHSdat%PLinVelEJ(p%DOFs%PIE(I),0,:)    ! The portion of PFrcWTail associated with the TFinMass
      TmpVec3 = CROSS_PRODUCT( RtHSdat%rWI, TmpVec1 )                      ! The portion of PMomNTail associated with the BoomMass
      TmpVec4 = CROSS_PRODUCT( RtHSdat%rWJ, TmpVec2 )                      ! The portion of PMomNTail associated with the TFinMass

      RtHSdat%PFrcWTail(:,p%DOFs%PIE(I)) = TmpVec1 + TmpVec2
      RtHSdat%PMomNTail(:,p%DOFs%PIE(I)) = TmpVec3 + TmpVec4 - p%AtfaIner*CoordSys%tfa* &
                                                             DOT_PRODUCT( CoordSys%tfa, RtHSdat%PAngVelEA(p%DOFs%PIE(I),0,:) )

   ENDDO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the tail boom center of mass (point I)

!.....................................
! FrcWTailt and MomNTailt - Forces on the tailfin
!.....................................
   ! Aerodynamic loads on TailFin CM (point K), with change of coordinate system
   Force(1:3)  = (/ u%TFinCMLoads%Force (1,1), u%TFinCMLoads%Force (3,1), -u%TFinCMLoads%Force (2,1) /)
   Moment(1:3) = (/ u%TFinCMLoads%Moment(1,1), u%TFinCMLoads%Moment(3,1), -u%TFinCMLoads%Moment(2,1) /)

   TmpVec1 = -p%BoomMass*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccEIt )                 ! The portion of FrcWTailt associated with the BoomMass
   TmpVec2 = -p%TFinMass*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccEJt )                 ! The portion of FrcWTailt associated with the TFinMass
   TmpVec3 = CROSS_PRODUCT( RtHSdat%rWI      , TmpVec1 )                           ! The portion of MomNTailt associated with the BoomMass
   TmpVec4 = CROSS_PRODUCT( RtHSdat%rWJ      , TmpVec2 )                           ! The portion of MomNTailt associated with the TFinMass
   TmpVec  = p%AtfaIner*CoordSys%tfa*DOT_PRODUCT( CoordSys%tfa, RtHSdat%AngVelEA )   ! = ( A inertia dyadic ) dot ( angular velocity of the tail in the inertia frame )
   TmpVec5 = CROSS_PRODUCT( -RtHSdat%AngVelEA, TmpVec  )                           ! = ( -angular velocity of the tail in the inertia frame ) cross ( TmpVec )

   RtHSdat%FrcWTailt = Force + TmpVec1 + TmpVec2
   RtHSdat%MomNTailt = Moment + TmpVec3 + TmpVec4 + TmpVec5         &
                     + CROSS_PRODUCT( RtHSdat%rWJ      , Force  )  &                         ! The portion of MomNTailt associated with Force with lever arm WK
                     - p%AtfaIner*CoordSys%tfa*DOT_PRODUCT( CoordSys%tfa, RtHSdat%AngAccEAt )   
   
!.....................................
! PFrcONcRt and PMomBNcRt
!  (requires PFrcVGnRt, PMomNGnRt, PFrcWTail, PMomNTail, )
!.....................................

   ! Define the partial forces and moments (including those associated with
   !   the QD2T()'s and those that are not) at the yaw bearing (point O) /
   !   base plate (body B) using the nacelle, generator, rotor, and tail effects.

   RtHSdat%PFrcONcRt = RtHSdat%PFrcVGnRt + RtHSdat%PFrcWTail   ! Initialize these partial forces and moments using
   RtHSdat%PMomBNcRt = RtHSdat%PMomNGnRt + RtHSdat%PMomNTail   ! the rotor, rotor-furl, generator, and tail effects
   DO I = 1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs

      TmpVec = CROSS_PRODUCT( RtHSdat%rOV, RtHSdat%PFrcVGnRt(:,p%DOFs%SrtPS(I)) ) ! The portion of PMomBNcRt associated with the PFrcVGnRt

      RtHSdat%PMomBNcRt(:,p%DOFs%SrtPS(I)) = RtHSdat%PMomBNcRt(:,p%DOFs%SrtPS(I)) + TmpVec

   ENDDO             ! I - All active (enabled) DOFs
   DO I = 1,p%DOFs%NPIE  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the tail boom center of mass (point I)

      TmpVec = CROSS_PRODUCT( RtHSdat%rOW, RtHSdat%PFrcWTail(:,p%DOFs%PIE(I)  ) ) ! The portion of PMomBNcRt associated with the PFrcWTail

      RtHSdat%PMomBNcRt(:,p%DOFs%PIE(I) ) = RtHSdat%PMomBNcRt(:,p%DOFs%PIE(I) ) + TmpVec

   ENDDO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the tail boom center of mass (point I)
   DO I = 1,p%DOFs%NPUE  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the nacelle center of mass (point U)

      TmpVec1 = -p%NacMass*RtHSdat%PLinVelEU(p%DOFs%PUE(I),0,:)              ! The portion of PFrcONcRt associated with the NacMass
      TmpVec2 = CROSS_PRODUCT( RtHSdat%rOU,               TmpVec1 ) ! The portion of PMomBNcRt associated with the NacMass

      RtHSdat%PFrcONcRt(:,p%DOFs%PUE(I)  ) = RtHSdat%PFrcONcRt(:,p%DOFs%PUE(I) ) + TmpVec1
      RtHSdat%PMomBNcRt(:,p%DOFs%PUE(I)  ) = RtHSdat%PMomBNcRt(:,p%DOFs%PUE(I) ) + TmpVec2 - p%Nacd2Iner*CoordSys%d2* &
                                             DOT_PRODUCT( CoordSys%d2, RtHSdat%PAngVelEN(p%DOFs%PUE(I),0,:) )

   ENDDO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the nacelle center of mass (point U)


!.....................................
! FrcONcRtt and MomBNcRtt
!  (requires FrcVGnRtt, MomNGnRtt, FrcWTailt, MomNTailt)
!.....................................

   TmpVec1 = -p%NacMass*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccEUt )                ! The portion of FrcONcRtt associated with the NacMass
   TmpVec2 = CROSS_PRODUCT( RtHSdat%rOU,           TmpVec1 )                         ! The portion of MomBNcRtt associated with the NacMass
   TmpVec3 = CROSS_PRODUCT( RtHSdat%rOV, RtHSdat%FrcVGnRtt )                         ! The portion of MomBNcRtt associated with the FrcVGnRtt
   TmpVec4 = CROSS_PRODUCT( RtHSdat%rOW, RtHSdat%FrcWTailt )                         ! The portion of MomBNcRtt associated with the FrcWTailt
   TmpVec  = p%Nacd2Iner*CoordSys%d2*DOT_PRODUCT( CoordSys%d2, RtHSdat%AngVelEN )    ! = ( Nacelle inertia dyadic ) dot ( angular velocity of nacelle in the inertia frame )
    
   RtHSdat%FrcONcRtt = RtHSdat%FrcVGnRtt + RtHSdat%FrcWTailt + TmpVec1 + (/ u%NacelleLoads%Force(1,1), u%NacelleLoads%Force(3,1), -u%NacelleLoads%Force(2,1) /)
   
   RtHSdat%MomBNcRtt = RtHSdat%MomNGnRtt + RtHSdat%MomNTailt + TmpVec2 + TmpVec3 + TmpVec4  &
                        + CROSS_PRODUCT( -RtHSdat%AngVelEN, TmpVec    )                     &    ! = ( -angular velocity of nacelle in the inertia frame ) cross ( TmpVec ) &
                        - p%Nacd2Iner*CoordSys%d2*DOT_PRODUCT( CoordSys%d2, RtHSdat%AngAccENt ) &
                        + (/ u%NacelleLoads%Moment(1,1), u%NacelleLoads%Moment(3,1), -u%NacelleLoads%Moment(2,1) /)

!.....................................
! PFTHydro and PMFHydro   
!  (requires TwrAddedMass)   
!.....................................

   ! Compute the partial hydrodynamic forces and moments per unit length
   !   (including those associated with the QD2T()'s and those that are not) at the current tower element (point T) / (body F):

   ! NOTE: These forces are named PFTHydro, PMFHydro, FTHydrot, and MFHydrot. However, the names should not imply that the 
   !       forces are a result of hydrodynamic contributions only.  These tower forces contain contributions from any external 
   !       load acting on the tower other than loads transmitted from aerodynamics.  For example, these tower forces contain 
   !       contributions from foundation stiffness and damping [not floating] or mooring line restoring and damping,
   !       as well as hydrostatic and hydrodynamic contributions [offshore].

   DO J=1,p%TwrNodes
   
      RtHSdat%PFTHydro(:,J,:) = 0.0
      RtHSdat%PMFHydro(:,J,:) = 0.0
      DO I = 1,p%DOFs%NPTE  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the tower

         RtHSdat%PFTHydro(:,J,p%DOFs%PTE(I)) = &
                                CoordSys%z1*( - u%TwrAddedMass(DOF_Sg,DOF_Sg,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,1) &
                                              + u%TwrAddedMass(DOF_Sg,DOF_Sw,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,3) &
                                              - u%TwrAddedMass(DOF_Sg,DOF_Hv,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,2) &
                                              - u%TwrAddedMass(DOF_Sg,DOF_R ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,1) &
                                              + u%TwrAddedMass(DOF_Sg,DOF_P ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,3) &
                                              - u%TwrAddedMass(DOF_Sg,DOF_Y ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,2)   ) &
                              - CoordSys%z3*( - u%TwrAddedMass(DOF_Sw,DOF_Sg,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,1) &
                                              + u%TwrAddedMass(DOF_Sw,DOF_Sw,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,3) &
                                              - u%TwrAddedMass(DOF_Sw,DOF_Hv,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,2) &
                                              - u%TwrAddedMass(DOF_Sw,DOF_R ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,1) &
                                              + u%TwrAddedMass(DOF_Sw,DOF_P ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,3) &
                                              - u%TwrAddedMass(DOF_Sw,DOF_Y ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,2)   ) &
                              + CoordSys%z2*( - u%TwrAddedMass(DOF_Hv,DOF_Sg,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,1) &
                                              + u%TwrAddedMass(DOF_Hv,DOF_Sw,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,3) &
                                              - u%TwrAddedMass(DOF_Hv,DOF_Hv,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,2) &
                                              - u%TwrAddedMass(DOF_Hv,DOF_R ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,1) &
                                              + u%TwrAddedMass(DOF_Hv,DOF_P ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,3) &
                                              - u%TwrAddedMass(DOF_Hv,DOF_Y ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,2)   )
         RtHSdat%PMFHydro(:,J,p%DOFs%PTE(I)) = &
                                CoordSys%z1*( - u%TwrAddedMass(DOF_R ,DOF_Sg,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,1) &
                                              + u%TwrAddedMass(DOF_R ,DOF_Sw,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,3) &
                                              - u%TwrAddedMass(DOF_R ,DOF_Hv,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,2) &
                                              - u%TwrAddedMass(DOF_R ,DOF_R ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,1) &
                                              + u%TwrAddedMass(DOF_R ,DOF_P ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,3) &
                                              - u%TwrAddedMass(DOF_R ,DOF_Y ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,2)   ) &
                              - CoordSys%z3*( - u%TwrAddedMass(DOF_P ,DOF_Sg,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,1) &
                                              + u%TwrAddedMass(DOF_P ,DOF_Sw,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,3) &
                                              - u%TwrAddedMass(DOF_P ,DOF_Hv,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,2) &
                                              - u%TwrAddedMass(DOF_P ,DOF_R ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,1) &
                                              + u%TwrAddedMass(DOF_P ,DOF_P ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,3) &
                                              - u%TwrAddedMass(DOF_P ,DOF_Y ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,2)   ) &
                              + CoordSys%z2*( - u%TwrAddedMass(DOF_Y ,DOF_Sg,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,1) &
                                              + u%TwrAddedMass(DOF_Y ,DOF_Sw,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,3) &
                                              - u%TwrAddedMass(DOF_Y ,DOF_Hv,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,2) &
                                              - u%TwrAddedMass(DOF_Y ,DOF_R ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,1) &
                                              + u%TwrAddedMass(DOF_Y ,DOF_P ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,3) &
                                              - u%TwrAddedMass(DOF_Y ,DOF_Y ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,2)   )

      END DO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the tower
   END DO !J

!.....................................
! FTHydrot and MFHydrot   
!  (requires TwrAddedMass)   
!.....................................
   
   DO J=1,p%TwrNodes
      RtHSdat%FTHydrot(:,J) = CoordSys%z1*( u%TowerPtLoads%Force(DOF_Sg,J)/abs(p%DHNodes(J)) &
                                                  - u%TwrAddedMass(DOF_Sg,DOF_Sg,J)*RtHSdat%LinAccETt(1,J) &
                                                  + u%TwrAddedMass(DOF_Sg,DOF_Sw,J)*RtHSdat%LinAccETt(3,J) &
                                                  - u%TwrAddedMass(DOF_Sg,DOF_Hv,J)*RtHSdat%LinAccETt(2,J) &
                                                  - u%TwrAddedMass(DOF_Sg,DOF_R ,J)*RtHSdat%AngAccEFt(1,J) &
                                                  + u%TwrAddedMass(DOF_Sg,DOF_P ,J)*RtHSdat%AngAccEFt(3,J) &
                                                  - u%TwrAddedMass(DOF_Sg,DOF_Y ,J)*RtHSdat%AngAccEFt(2,J)   ) &
                            - CoordSys%z3*( u%TowerPtLoads%Force(DOF_Sw,J)/abs(p%DHNodes(J)) &
                                                  - u%TwrAddedMass(DOF_Sw,DOF_Sg,J)*RtHSdat%LinAccETt(1,J) &
                                                  + u%TwrAddedMass(DOF_Sw,DOF_Sw,J)*RtHSdat%LinAccETt(3,J) &
                                                  - u%TwrAddedMass(DOF_Sw,DOF_Hv,J)*RtHSdat%LinAccETt(2,J) &
                                                  - u%TwrAddedMass(DOF_Sw,DOF_R ,J)*RtHSdat%AngAccEFt(1,J) &
                                                  + u%TwrAddedMass(DOF_Sw,DOF_P ,J)*RtHSdat%AngAccEFt(3,J) &
                                                  - u%TwrAddedMass(DOF_Sw,DOF_Y ,J)*RtHSdat%AngAccEFt(2,J)   ) &
                             + CoordSys%z2*( u%TowerPtLoads%Force(DOF_Hv,J)/abs(p%DHNodes(J)) &
                                                  - u%TwrAddedMass(DOF_Hv,DOF_Sg,J)*RtHSdat%LinAccETt(1,J) &
                                                  + u%TwrAddedMass(DOF_Hv,DOF_Sw,J)*RtHSdat%LinAccETt(3,J) &
                                                  - u%TwrAddedMass(DOF_Hv,DOF_Hv,J)*RtHSdat%LinAccETt(2,J) &
                                                  - u%TwrAddedMass(DOF_Hv,DOF_R ,J)*RtHSdat%AngAccEFt(1,J) &
                                                  + u%TwrAddedMass(DOF_Hv,DOF_P ,J)*RtHSdat%AngAccEFt(3,J) &
                                                  - u%TwrAddedMass(DOF_Hv,DOF_Y ,J)*RtHSdat%AngAccEFt(2,J)   )
      RtHSdat%MFHydrot(:,J) = CoordSys%z1*( u%TowerPtLoads%Moment(DOF_R-3,J)/abs(p%DHNodes(J)) &
                                                  - u%TwrAddedMass(DOF_R ,DOF_Sg,J)*RtHSdat%LinAccETt(1,J) &
                                                  + u%TwrAddedMass(DOF_R ,DOF_Sw,J)*RtHSdat%LinAccETt(3,J) &
                                                  - u%TwrAddedMass(DOF_R ,DOF_Hv,J)*RtHSdat%LinAccETt(2,J) &
                                                  - u%TwrAddedMass(DOF_R ,DOF_R ,J)*RtHSdat%AngAccEFt(1,J) &
                                                  + u%TwrAddedMass(DOF_R ,DOF_P ,J)*RtHSdat%AngAccEFt(3,J) &
                                                  - u%TwrAddedMass(DOF_R ,DOF_Y ,J)*RtHSdat%AngAccEFt(2,J)   ) &
                            - CoordSys%z3*( u%TowerPtLoads%Moment(DOF_P-3 ,J)/abs(p%DHNodes(J)) &
                                                  - u%TwrAddedMass(DOF_P ,DOF_Sg,J)*RtHSdat%LinAccETt(1,J) &
                                                  + u%TwrAddedMass(DOF_P ,DOF_Sw,J)*RtHSdat%LinAccETt(3,J) &
                                                  - u%TwrAddedMass(DOF_P ,DOF_Hv,J)*RtHSdat%LinAccETt(2,J) &
                                                  - u%TwrAddedMass(DOF_P ,DOF_R ,J)*RtHSdat%AngAccEFt(1,J) &
                                                  + u%TwrAddedMass(DOF_P ,DOF_P ,J)*RtHSdat%AngAccEFt(3,J) &
                                                  - u%TwrAddedMass(DOF_P ,DOF_Y ,J)*RtHSdat%AngAccEFt(2,J)   ) &
                            + CoordSys%z2*( u%TowerPtLoads%Moment(DOF_Y-3 ,J)/abs(p%DHNodes(J)) &
                                                  - u%TwrAddedMass(DOF_Y ,DOF_Sg,J)*RtHSdat%LinAccETt(1,J) &
                                                  + u%TwrAddedMass(DOF_Y ,DOF_Sw,J)*RtHSdat%LinAccETt(3,J) &
                                                  - u%TwrAddedMass(DOF_Y ,DOF_Hv,J)*RtHSdat%LinAccETt(2,J) &
                                                  - u%TwrAddedMass(DOF_Y ,DOF_R ,J)*RtHSdat%AngAccEFt(1,J) &
                                                  + u%TwrAddedMass(DOF_Y ,DOF_P ,J)*RtHSdat%AngAccEFt(3,J) &
                                                  - u%TwrAddedMass(DOF_Y ,DOF_Y ,J)*RtHSdat%AngAccEFt(2,J)   )

   END DO !J
   
!.....................................
! PFrcT0Trb and PMomX0Trb
!  (requires PFrcONcRt, PMomBNcRt, PFrcT0Trb, PMomX0Trb, PFTHydro, PMFHydro)
!.....................................

      ! Initialize the partial forces and moments (including those associated
      !   with the QD2T()'s and those that are not) at the tower base (point T(0))  using everything but the tower:

   RtHSdat%PFrcT0Trb = RtHSdat%PFrcONcRt   ! Initialize these partial forces and moments
   RtHSdat%PMomX0Trb = RtHSdat%PMomBNcRt   ! using all of the effects above the yaw bearing
   DO I = 1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs

      TmpVec  = CROSS_PRODUCT(  RtHSdat%rT0O, RtHSdat%PFrcONcRt(:,p%DOFs%SrtPS(I)) )   ! The portion of PMomX0Trb associated with the PFrcONcRt

      RtHSdat%PMomX0Trb(:,p%DOFs%SrtPS(I)) = RtHSdat%PMomX0Trb(:,p%DOFs%SrtPS(I)) + TmpVec

   ENDDO             ! I - All active (enabled) DOFs
   DO I = 1,p%DOFs%NPTE  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the yaw bearing center of mass (point O)

      TmpVec1 = -p%YawBrMass*RtHSdat%PLinVelEO(p%DOFs%PTE(I),0,:)               ! The portion of PFrcT0Trb associated with the YawBrMass
      TmpVec2 = CROSS_PRODUCT( RtHSdat%rT0O,               TmpVec1 )   ! The portion of PMomX0Trb associated with the YawBrMass

      RtHSdat%PFrcT0Trb(:,p%DOFs%PTE(I)  ) = RtHSdat%PFrcT0Trb(:,p%DOFs%PTE(I)  ) + TmpVec1
      RtHSdat%PMomX0Trb(:,p%DOFs%PTE(I)  ) = RtHSdat%PMomX0Trb(:,p%DOFs%PTE(I)  ) + TmpVec2

   ENDDO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the yaw bearing center of mass (point O)

   TmpVec1 = -p%YawBrMass*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccEOt ) ! The portion of FrcT0Trbt associated with the YawBrMass
   TmpVec2 = CROSS_PRODUCT( RtHSdat%rT0O,   TmpVec1 )               ! The portion of MomX0Trbt associated with the YawBrMass
   TmpVec3 = CROSS_PRODUCT( RtHSdat%rT0O, RtHSdat%FrcONcRtt )               ! The portion of MomX0Trbt associated with the FrcONcRtt

   RtHSdat%FrcT0Trbt = RtHSdat%FrcONcRtt + TmpVec1
   RtHSdat%MomX0Trbt = RtHSdat%MomBNcRtt + TmpVec2 + TmpVec3   
   
   ! Integrate to find the total partial forces and moments (including those
   !   associated with the QD2T()'s and those that are not) at the tower base (point T(0)):
   DO J=1,p%TwrNodes

      DO I = 1,p%DOFs%NPTE  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the tower

         TmpVec1 = RtHSdat%PFTHydro(:,J,p%DOFs%PTE(I))*abs(p%DHNodes(J)) - p%TElmntMass(J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,:)           ! The portion of PFrcT0Trb associated with tower element J
         TmpVec2 = CROSS_PRODUCT( RtHSdat%rT0T(:,J), TmpVec1 )                 ! The portion of PMomX0Trb associated with tower element J
         TmpVec3 = RtHSdat%PMFHydro(:,J,p%DOFs%PTE(I))*abs(p%DHNodes(J))             ! The added moment applied at tower element J

         RtHSdat%PFrcT0Trb(:,p%DOFs%PTE(I)) = RtHSdat%PFrcT0Trb(:,p%DOFs%PTE(I)) + TmpVec1
         RtHSdat%PMomX0Trb(:,p%DOFs%PTE(I)) = RtHSdat%PMomX0Trb(:,p%DOFs%PTE(I)) + TmpVec2 + TmpVec3

      ENDDO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the tower

      TmpVec1 = ( RtHSdat%FTHydrot(:,J) )*abs(p%DHNodes(J)) &
              - p%TElmntMass(J)*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccETt(:,J) )          ! The portion of FrcT0Trbt associated with tower element J
      TmpVec2 = CROSS_PRODUCT( RtHSdat%rT0T(:,J), TmpVec1 )                                 ! The portion of MomX0Trbt associated with tower element J
      TmpVec3 = ( RtHSdat%MFHydrot(:,J) )*abs(p%DHNodes(J))                                      ! The external moment applied to tower element J

      RtHSdat%FrcT0Trbt = RtHSdat%FrcT0Trbt + TmpVec1

      RtHSdat%MomX0Trbt = RtHSdat%MomX0Trbt + TmpVec2 + TmpVec3

   END DO !J

!.....................................
! PFZHydro and  PMXHydro  
!  ( requires PtfmAddedMass )
!.....................................   
   
   !..................................................................................................................................
   ! Compute the partial platform forces and moments (including those associated with the QD2T()'s and those that are not) at the
   ! platform reference (point Z) / (body X).
   !
   ! NOTE: These forces are named PFZHydro, PMXHydro, FZHydrot, and MXHydrot. However, the names should not imply that the forces
   !   are a result of hydrodynamic contributions only. These platform forces contain contributions from any external load acting
   !   on the platform other than loads transmitted from the wind turbine. For example, these platform forces contain contributions
   !   from foundation stiffness and damping [not floating] or mooring line restoring and damping [floating], as well as hydrostatic
   !   and hydrodynamic contributions [offshore].
   !bjj: m%RtHS%PFZHydro, %PMXHydro, %FZHydrot, and %MXHydrot are not used in the output routine anymore
   !      (because of their dependence on inputs, u)

   RtHSdat%PFZHydro = 0.0
   RtHSdat%PMXHydro = 0.0
   DO I = 1,p%DOFs%NPYE  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the platform center of mass (point Y)

      RtHSdat%PFZHydro(p%DOFs%PYE(I),:) = - u%PtfmAddedMass(DOF_Sg,p%DOFs%PYE(I))*RtHSdat%PLinVelEZ(DOF_Sg,0,:) &
                                          - u%PtfmAddedMass(DOF_Sw,p%DOFs%PYE(I))*RtHSdat%PLinVelEZ(DOF_Sw,0,:) &
                                          - u%PtfmAddedMass(DOF_Hv,p%DOFs%PYE(I))*RtHSdat%PLinVelEZ(DOF_Hv,0,:)
      RtHSdat%PMXHydro(p%DOFs%PYE(I),:) = - u%PtfmAddedMass(DOF_R ,p%DOFs%PYE(I))*RtHSdat%PAngVelEX(DOF_R ,0,:) &
                                          - u%PtfmAddedMass(DOF_P ,p%DOFs%PYE(I))*RtHSdat%PAngVelEX(DOF_P ,0,:) &
                                          - u%PtfmAddedMass(DOF_Y ,p%DOFs%PYE(I))*RtHSdat%PAngVelEX(DOF_Y ,0,:)

   ENDDO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the platform center of mass (point Y)

   RtHSdat%FZHydrot = u%PlatformPtMesh%Force(DOF_Sg,1)*RtHSdat%PLinVelEZ(DOF_Sg,0,:) &
                    + u%PlatformPtMesh%Force(DOF_Sw,1)*RtHSdat%PLinVelEZ(DOF_Sw,0,:) &
                    + u%PlatformPtMesh%Force(DOF_Hv,1)*RtHSdat%PLinVelEZ(DOF_Hv,0,:)
   RtHSdat%MXHydrot = u%PlatformPtMesh%Moment(DOF_R-3,1)*RtHSdat%PAngVelEX(DOF_R ,0,:) &
                    + u%PlatformPtMesh%Moment(DOF_P-3,1)*RtHSdat%PAngVelEX(DOF_P ,0,:) &
                    + u%PlatformPtMesh%Moment(DOF_Y-3,1)*RtHSdat%PAngVelEX(DOF_Y ,0,:)
   
!.....................................
! PFrcZAll and PMomXAll  
!  (requires PFrcT0Trb, PMomX0Trb, PFZHydro, PMXHydro )
!.....................................   

   ! Define the partial forces and moments (including those associated with the QD2T()'s and those that are not) at the
   !   platform reference (point Z) / (body X) using the turbine and platform effects:

   RtHSdat%PFrcZAll = RtHSdat%PFrcT0Trb ! Initialize these partial forces and moments
   RtHSdat%PMomXAll = RtHSdat%PMomX0Trb ! using the effects from the wind turbine
   DO I = 1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs

      TmpVec = CROSS_PRODUCT( RtHSdat%rZT0, RtHSdat%PFrcT0Trb(:,p%DOFs%SrtPS(I)) )   ! The portion of PMomXAll associated with the PFrcT0Trb

      RtHSdat%PMomXAll(:,p%DOFs%SrtPS(I)) = RtHSdat%PMomXAll(:,p%DOFs%SrtPS(I)) + TmpVec

   ENDDO             ! I - All active (enabled) DOFs
   DO I = 1,p%DOFs%NPYE  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the platform center of mass (point Y)

      TmpVec1 = -p%PtfmMass*RtHSdat%PLinVelEY(p%DOFs%PYE(I),0,:)                ! The portion of PFrcZAll associated with the PtfmMass
      TmpVec2 = CROSS_PRODUCT( RtHSdat%rZY ,               TmpVec1 )   ! The portion of PMomXAll associated with the PtfmMass

      RtHSdat%PFrcZAll(:,p%DOFs%PYE(I)) = RtHSdat%PFrcZAll(:,p%DOFs%PYE(I)  )        + RtHSdat%PFZHydro(p%DOFs%PYE(I),:) + TmpVec1
      RtHSdat%PMomXAll(:,p%DOFs%PYE(I)) = RtHSdat%PMomXAll(:,p%DOFs%PYE(I)  )        + RtHSdat%PMXHydro(p%DOFs%PYE(I),:) + TmpVec2 &
                                    - p%PtfmRIner*CoordSys%a1*DOT_PRODUCT( CoordSys%a1, RtHSdat%PAngVelEX(p%DOFs%PYE(I),0,:) )   &
                                    - p%PtfmYIner*CoordSys%a2*DOT_PRODUCT( CoordSys%a2, RtHSdat%PAngVelEX(p%DOFs%PYE(I),0,:) )   &
                                    - p%PtfmPIner*CoordSys%a3*DOT_PRODUCT( CoordSys%a3, RtHSdat%PAngVelEX(p%DOFs%PYE(I),0,:) )

   ENDDO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the platform center of mass (point Y)

!.....................................
! FrcZAllt and MomXAllt
!  (requires FrcT0Trbt, MomX0Trbt)
!.....................................

   TmpVec1 = -p%PtfmMass*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccEYt  )                                              ! The portion of FrcZAllt associated with the PtfmMass
   TmpVec2 = CROSS_PRODUCT( RtHSdat%rZY      ,   TmpVec1 )                                                                      ! The portion of MomXAllt associated with the PtfmMass
   TmpVec3 = CROSS_PRODUCT( RtHSdat%rZT0     , RtHSdat%FrcT0Trbt )                                                      ! The portion of MomXAllt associated with the FrcT0Trbt
   TmpVec  = p%PtfmRIner*CoordSys%a1*DOT_PRODUCT( CoordSys%a1, RtHSdat%AngVelEX  ) &      ! = ( Platform inertia dyadic ) dot ( angular velocity of platform in the inertia frame )
           + p%PtfmYIner*CoordSys%a2*DOT_PRODUCT( CoordSys%a2, RtHSdat%AngVelEX  ) &
           + p%PtfmPIner*CoordSys%a3*DOT_PRODUCT( CoordSys%a3, RtHSdat%AngVelEX  )
   TmpVec4 = CROSS_PRODUCT( -RtHSdat%AngVelEX,   TmpVec  )                                                      ! = ( -angular velocity of platform in the inertia frame ) cross ( TmpVec )

   RtHSdat%FrcZAllt = RtHSdat%FrcT0Trbt + RtHSdat%FZHydrot + TmpVec1
   RtHSdat%MomXAllt = RtHSdat%MomX0Trbt + RtHSdat%MXHydrot + TmpVec2 + TmpVec3 + TmpVec4   
   
   
END SUBROUTINE CalculateForcesMoments
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is used to populate the AugMat matrix for RtHS (CalcContStateDeriv)
SUBROUTINE FillAugMat( p, x, CoordSys, u, HSSBrTrq, RtHSdat, AugMat )
!..................................................................................................................................

      ! Passed variables
   TYPE(ED_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(ED_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
   TYPE(ED_CoordSys),            INTENT(IN   )  :: CoordSys    !< The coordinate systems that have been set for these states/time
   TYPE(ED_InputType),           INTENT(IN   )  :: u           !< The aero blade forces/moments
   TYPE(ED_RtHndSide),           INTENT(IN   )  :: RtHSdat     !< data from the RtHndSid module (contains positions to be set)
   REAL(ReKi),                   INTENT(IN )    :: HSSBrTrq    !<  SIGN( u%HSSBrTrqC, x%QDT(DOF_GeAz) ) or corrected value from FixHSS
   REAL(R8Ki),                   INTENT(OUT)    :: AugMat(:,:) !< the return matrix 
   
      ! Local variables
   REAL(ReKi)                   :: TmpVec    (3)                                   ! A temporary vector used in various computations.
   REAL(ReKi)                   :: TmpVec1   (3)                                   ! A temporary vector used in various computations.
   REAL(ReKi)                   :: TmpVec3   (3)                                   ! A temporary vector used in various computations.
   REAL(ReKi)                   :: GBoxTrq                                         ! Gearbox torque on the LSS side in N-m (calculated from inputs and parameters).
   REAL(ReKi)                   :: GBoxEffFac2                                     ! A second gearbox efficiency factor = ( 1 / GBoxEff^SgnPrvLSTQ - 1 )

   INTEGER(IntKi)               :: I                                               ! Loops through some or all of the DOFs
   INTEGER(IntKi)               :: J                                               ! Counter for elements
   INTEGER(IntKi)               :: K                                               ! Counter for blades
   INTEGER(IntKi)               :: L                                               ! Generic index

   
      ! Initialize the matrix:
      
   AugMat      = 0.0
   GBoxTrq    = ( u%GenTrq + HSSBrTrq )*ABS(p%GBRatio) ! bjj: do we use HSSBrTrqC or HSSBrTrq?
   
   DO K = 1,p%NumBl ! Loop through all blades
   

      ! Initialize the portions of the mass matrix on and below the diagonal associated with purely blade DOFs (these portions can't
      !   be calculated using partial loads) using the tip mass effects. 
      ! Also, initialize the portions of the forcing vector associated with purely blade DOFs (these portions can't be calculated 
      !   using partial loads) using the tip mass effects:
      ! NOTE: The vector subscript array, PSBE(), used in the following loops must be sorted from smallest to largest DOF index in 
      !       order for the loops to work to enter values only on and below the diagonal of AugMat():

   
      DO L = 1,p%DOFs%NPSBE(K)    ! Loop through all active (enabled) blade DOFs that contribute to the QD2T-related linear accelerations of the tip of blade K (point S(p%BldFlexL))
         DO I = L,p%DOFs%NPSBE(K) ! Loop through all active (enabled) blade DOFs greater than or equal to L
            AugMat(p%DOFs%PSBE(K,I),p%DOFs%PSBE(K,L)) = p%TipMass(K)*&
                                        DOT_PRODUCT( RtHSdat%PLinVelES(K, p%TipNode, p%DOFs%PSBE(K,I),0,:), &   ! [C(q,t)]B
                                                     RtHSdat%PLinVelES(K, p%TipNode, p%DOFs%PSBE(K,L),0,:)    )
         ENDDO             ! I - All active (enabled) blade DOFs greater than or equal to L
      ENDDO                ! L - All active (enabled) blade DOFs that contribute to the QD2T-related linear accelerations of the tip of blade K (point S(p%BldFlexL))

      TmpVec1 = RtHSdat%FSTipDrag(:,K) - p%TipMass(K)*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccESt(:,K,p%TipNode) ) ! The portion of FrcS0Bt associated with the tip brake
      DO I = 1,p%DOFs%NPSBE(K)    ! Loop through all active (enabled) blade DOFs that contribute to the QD2T-related linear accelerations of the tip of blade K (point S(p%BldFlexL))
            AugMat(p%DOFs%PSBE(K,I), p%NAug) = DOT_PRODUCT( RtHSdat%PLinVelES(K,p%TipNode,p%DOFs%PSBE(K,I),0,:), &   ! {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB
                                                              TmpVec1                               ) ! NOTE: TmpVec1 is still the portion of FrcS0Bt associated with the tip brake
      ENDDO                ! I - All active (enabled) blade DOFs that contribute to the QD2T-related linear accelerations of the tip of blade K (point S(p%BldFlexL))
   
      

      DO J = 1,p%BldNodes ! Loop through the blade nodes / elements


      ! Integrate to find the portions of the mass matrix on and below the diagonal associated with purely blade DOFs (these portions
      !   can't be calculated using partial loads).  Also, integrate to find the portions of the forcing vector associated with
      !    purely blade DOFs (these portions can't be calculated using partial loads):
      ! NOTE: The vector subscript array, PSBE(), used in the following loops must
      !       be sorted from smallest to largest DOF index in order for the loops
      !       to work to enter values only on and below the diagonal of AugMat():
  
         DO L = 1,p%DOFs%NPSBE(K)    ! Loop through all active (enabled) blade DOFs that contribute to the QD2T-related linear accelerations of the blade
            DO I = L,p%DOFs%NPSBE(K) ! Loop through all active (enabled) blade DOFs greater than or equal to L
               AugMat(p%DOFs%PSBE(K,I),p%DOFs%PSBE(K,L)) = AugMat(p%DOFs%PSBE(K,I),p%DOFs%PSBE(K,L)) + p%BElmntMass(J,K)*&
                                             DOT_PRODUCT( RtHSdat%PLinVelES(K,J,p%DOFs%PSBE(K,I),0,:), &           ! [C(q,t)]B
                                                          RtHSdat%PLinVelES(K,J,p%DOFs%PSBE(K,L),0,:)   )
            ENDDO             ! I - All active (enabled) blade DOFs greater than or equal to L
         ENDDO                ! L - All active (enabled) blade DOFs that contribute to the QD2T-related linear accelerations of the blade
      
         TmpVec1 = RtHSdat%FSAero(:,K,J)*p%DRNodes(J) - p%BElmntMass(J,K)*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccESt(:,K,J) ) ! The portion of FrcS0Bt associated with blade element J
         TmpVec3 = RtHSdat%MMAero(:,K,J)*p%DRNodes(J)                                               ! The total external moment applied to blade element J
         DO I = 1,p%DOFs%NPSBE(K)    ! Loop through all active (enabled) blade DOFs that contribute to the QD2T-related linear accelerations of the blade
               AugMat(p%DOFs%PSBE(K,I), p%NAug) = AugMat(p%DOFs%PSBE(K,I),     p%NAug)                      & ! {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB
                                           + DOT_PRODUCT( RtHSdat%PLinVelES(K,J,p%DOFs%PSBE(K,I),0,:), TmpVec1 ) & ! NOTE: TmpVec1 is still the portion of FrcS0Bt associated with blade element J
                                           + DOT_PRODUCT( RtHSdat%PAngVelEM(K,J,p%DOFs%PSBE(K,I),0,:), TmpVec3 )   !       and TmpVec3 is still the total external moment applied to blade element J
         ENDDO                ! I - All active (enabled) blade DOFs that contribute to the QD2T-related linear accelerations of the blade


      ENDDO ! J - Blade nodes / elements      
      
      
      

      ! Initialize the portions of the mass matrix below the diagonal associated
      !   with the teeter and pure blade DOFs using the partial loads at the teeter pin; only do this if necessary:

      IF ( ( p%NumBl == 2 ) ) THEN
         IF ( p%DOF_Flag(DOF_Teet) ) THEN ! NOTE: two "ifs" since DOF_Teet might be out of bound
            DO L = 1,p%DOFs%NPSBE(K) ! Loop through all active (enabled) blade DOFs that contribute to the QD2T-related linear accelerations of the blade
               AugMat(DOF_Teet,p%DOFs%PSBE(K,L)) = -DOT_PRODUCT( RtHSdat%PAngVelEH(DOF_Teet,0,:), &
                                                                 RtHSdat%PMomLPRot(:,p%DOFs%PSBE(K,L)) )  ! [C(q,t)]B
            ENDDO             ! L - All active (enabled) blade DOFs that contribute to the QD2T-related linear accelerations of the blade
         ENDIF
      ENDIF



      ! If the associated DOFs are enabled, add the blade elasticity and damping
      !   forces to the forcing vector (these portions can't be calculated using
      !   partial loads):

      IF ( p%DOF_Flag(DOF_BF(K,1)) )  THEN
         AugMat(    DOF_BF(K,1),p%NAug) = AugMat(DOF_BF(K,1),p%NAug)      & !
                                        - p%KBF(K,1,1)*x%QT( DOF_BF(K,1)) &
                                        - p%KBF(K,1,2)*x%QT( DOF_BF(K,2)) &
                                        - p%CBF(K,1,1)*x%QDT(DOF_BF(K,1)) &
                                        - p%CBF(K,1,2)*x%QDT(DOF_BF(K,2))
      ENDIF
      IF ( p%DOF_Flag(DOF_BF(K,2)) )  THEN
         AugMat(    DOF_BF(K,2),p%NAug) = AugMat(DOF_BF(K,2),p%NAug)      & ! {-f(qd,q,t)}ElasticB + {-f(qd,q,t)}DampB
                                        - p%KBF(K,2,1)*x%QT( DOF_BF(K,1)) &
                                        - p%KBF(K,2,2)*x%QT( DOF_BF(K,2)) &
                                        - p%CBF(K,2,1)*x%QDT(DOF_BF(K,1)) &
                                        - p%CBF(K,2,2)*x%QDT(DOF_BF(K,2))
      ENDIF
      IF ( p%DOF_Flag(DOF_BE(K,1)) )  THEN
         AugMat(    DOF_BE(K,1),p%NAug) = AugMat(DOF_BE(K,1),p%NAug)      & !
                                        - p%KBE(K,1,1)*x%QT( DOF_BE(K,1)) &
                                        - p%CBE(K,1,1)*x%QDT(DOF_BE(K,1))
      ENDIF
                  
      
   END DO !k
         
   
      ! Initialize the portions of the mass matrix on and below the diagonal
      !   associated with purely tower DOFs (these portions can't be calculated
      !   using partial loads) using the yaw bearing mass effects.
      !   Also, initialize the portions of the forcing vector associated with
      !   purely blade DOFs (these portions can't be calculated using partial
      !   loads) using the yaw bearing mass effects:
      ! NOTE: The vector subscript array, PTTE(), used in the following loops must
      !       be sorted from smallest to largest DOF index in order for the loops
      !       to work to enter values only on and below the diagonal of AugMat():

   DO L = 1,p%DOFs%NPTTE    ! Loop through all active (enabled) tower DOFs that contribute to the QD2T-related linear accelerations of the yaw bearing (point O)
      DO I = L,p%DOFs%NPTTE ! Loop through all active (enabled) tower DOFs greater than or equal to L
         AugMat(p%DOFs%PTTE(I),p%DOFs%PTTE(L)) = p%YawBrMass*DOT_PRODUCT( RtHSdat%PLinVelEO(p%DOFs%PTTE(I),0,:), &     ! [C(q,t)]T of YawBrMass
                                                      RtHSdat%PLinVelEO(p%DOFs%PTTE(L),0,:)    )
      ENDDO          ! I - All active (enabled) tower DOFs greater than or equal to L
   ENDDO             ! L - All active (enabled) tower DOFs that contribute to the QD2T-related linear accelerations of the yaw bearing (point O)

   TmpVec1 = -p%YawBrMass*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccEOt ) ! The portion of FrcT0Trbt associated with the YawBrMass
   DO I = 1,p%DOFs%NPTTE    ! Loop through all active (enabled) tower DOFs that contribute to the QD2T-related linear accelerations of the yaw bearing (point O)
         AugMat(p%DOFs%PTTE(I),   p%NAug) =           DOT_PRODUCT( RtHSdat%PLinVelEO(p%DOFs%PTTE(I),0,:), &     ! {-f(qd,q,t)}T + {-f(qd,q,t)}GravT of YawBrMass
                                                      TmpVec1                   )   ! NOTE: TmpVec1 is still the portion of FrcT0Trbt associated with YawBrMass
   ENDDO             ! I - All active (enabled) tower DOFs that contribute to the QD2T-related linear accelerations of the yaw bearing (point O)
   
   
   
   

   DO J = 1,p%TwrNodes

   !..................................................................................................................................
   ! Integrate to find the portions of the mass matrix on and below the diagonal associated with purely tower DOFs (these portions
   !   can't be calculated using partial loads).  Also, integrate to find the portions of the forcing vector associated with purely
   !   tower DOFs (these portions can't be calculated using partial loads).
   ! NOTE: The vector subscript array, PTTE(), used in the following loops must be sorted from smallest to largest DOF index in order
   !   for the loops to work to enter values only on and below the diagonal of AugMat():
   !..................................................................................................................................

      DO L = 1,p%DOFs%NPTTE    ! Loop through all active (enabled) tower DOFs that contribute to the QD2T-related linear accelerations of the tower
         DO I = L,p%DOFs%NPTTE ! Loop through all active (enabled) tower DOFs greater than or equal to L
            AugMat(p%DOFs%PTTE(I),p%DOFs%PTTE(L)) = AugMat(p%DOFs%PTTE(I),p%DOFs%PTTE(L))  &
                                                  + p%TElmntMass(J)   *DOT_PRODUCT( RtHSdat%PLinVelET(J,p%DOFs%PTTE(I),0,:),  &
                                                                              RtHSdat%PLinVelET(J,p%DOFs%PTTE(L),0,:) ) &   ! [C(q,t)]T + [C(q,t)]HydroT
                                                  - abs(p%DHNodes(J))*DOT_PRODUCT( RtHSdat%PLinVelET(J,p%DOFs%PTTE(I),0,:),  &
                                                                              RtHSdat%PFTHydro (:,J,p%DOFs%PTTE(L)  ) ) &
                                                  - abs(p%DHNodes(J))*DOT_PRODUCT( RtHSdat%PAngVelEF(J,p%DOFs%PTTE(I),0,:),  &
                                                                              RtHSdat%PMFHydro (:,J,p%DOFs%PTTE(L)  ) )
         ENDDO                 ! I - All active (enabled) tower DOFs greater than or equal to L
      ENDDO                    ! L - All active (enabled) tower DOFs that contribute to the QD2T-related linear accelerations of the tower

      TmpVec1 = ( RtHSdat%FTHydrot(:,J) )*abs(p%DHNodes(J)) &
              - p%TElmntMass(J)*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccETt(:,J) )          ! The portion of FrcT0Trbt associated with tower element J
      TmpVec3 = ( RtHSdat%MFHydrot(:,J) )*abs(p%DHNodes(J))             ! The external moment applied to tower element J
      DO I = 1,p%DOFs%NPTTE    ! Loop through all active (enabled) tower DOFs that contribute to the QD2T-related linear accelerations of the tower
            AugMat(p%DOFs%PTTE(I),        p%NAug) = AugMat(p%DOFs%PTTE(I),   p%NAug)                         &                 ! {-f(qd,q,t)}T + {-f(qd,q,t)}GravT + {-f(qd,q,t)}AeroT + {-f(qd,q,t)}HydroT
                                                  +  DOT_PRODUCT( RtHSdat%PLinVelET(J,p%DOFs%PTTE(I),0,:), TmpVec1        ) &  ! NOTE: TmpVec1 is still the portion of FrcT0Trbt associated with tower element J
                                                  +  DOT_PRODUCT( RtHSdat%PAngVelEF(J,p%DOFs%PTTE(I),0,:), TmpVec3        )    !       and TmpVec3 is still the total external moment to tower element J
      ENDDO                    ! I - All active (enabled) tower DOFs that contribute to the QD2T-related linear accelerations of the tower

   ENDDO ! J - Tower nodes / elements   
   
   !..................................................................................................................................
   ! If the associated DOFs are enabled, add the tower elasticity and damping forces to the forcing vector (these portions can't be
   !   calculated using partial loads):
   !..................................................................................................................................

   IF ( p%DOF_Flag(DOF_TFA1) )  THEN
      AugMat(    DOF_TFA1,p%NAug) = AugMat(DOF_TFA1,p%NAug)                                   &
                                  - p%KTFA(1,1)*x%QT( DOF_TFA1) - p%KTFA(1,2)*x%QT( DOF_TFA2) &                                     !
                                  - p%CTFA(1,1)*x%QDT(DOF_TFA1) - p%CTFA(1,2)*x%QDT(DOF_TFA2)
   ENDIF
   IF ( p%DOF_Flag(DOF_TSS1) )  THEN
      AugMat(    DOF_TSS1,p%NAug) = AugMat(DOF_TSS1,p%NAug)                                   &
                                  - p%KTSS(1,1)*x%QT( DOF_TSS1) - p%KTSS(1,2)*x%QT( DOF_TSS2) &                                     ! {-f(qd,q,t)}ElasticT + {-f(qd,q,t)}DampT
                                  - p%CTSS(1,1)*x%QDT(DOF_TSS1) - p%CTSS(1,2)*x%QDT(DOF_TSS2)
   ENDIF
   IF ( p%DOF_Flag(DOF_TFA2) )  THEN
      AugMat(    DOF_TFA2,p%NAug) = AugMat(DOF_TFA2,p%NAug)                                   &
                                  - p%KTFA(2,1)*x%QT( DOF_TFA1) - p%KTFA(2,2)*x%QT( DOF_TFA2) &                                     !
                                  - p%CTFA(2,1)*x%QDT(DOF_TFA1) - p%CTFA(2,2)*x%QDT(DOF_TFA2)
   ENDIF
   IF ( p%DOF_Flag(DOF_TSS2) )  THEN
      AugMat(    DOF_TSS2,p%NAug) = AugMat(DOF_TSS2,p%NAug)                                   &
                                  - p%KTSS(2,1)*x%QT( DOF_TSS1) - p%KTSS(2,2)*x%QT( DOF_TSS2) &                                     !
                                  - p%CTSS(2,1)*x%QDT(DOF_TSS1) - p%CTSS(2,2)*x%QDT(DOF_TSS2)
   ENDIF
   
   
   
!..................................................................................................................................
! Now that all of the partial loads have been found, let's fill in the portions of the mass matrix on and below the diagonal that
! may be calculated with the help of the partial loads.
! Also, let's fill in the portions of the forcing vector that may be calculated with the help of the partial loads.
! Also let's add in additional terms to the forcing function that can't be added with the help of the partial loads.
!
! NOTE: The vector subscript array, SrtPS(), used in the following loops must be sorted from smallest to largest DOF index in order
!   for the loops to work to enter values only on and below the diagonal of AugMat():
!..................................................................................................................................

   IF ( p%DOF_Flag (DOF_Sg  ) )  THEN
      DO I = p%DOFs%Diag(DOF_Sg  ),p%DOFs%NActvDOF   ! Loop through all active (enabled) DOFs on or below the diagonal
         AugMat(p%DOFs%SrtPS(I),DOF_Sg  ) = -1.*DOT_PRODUCT( RtHSdat%PLinVelEZ(DOF_Sg ,0,:), RtHSdat%PFrcZAll (:,p%DOFs%SrtPS(I)) ) ! [C(q,t)]X + [C(q,t)]HydroX + [C(q,t)]T + [C(q,t)]HydroT + [C(q,t)]N + [C(q,t)]R + [C(q,t)]H + [C(q,t)]B + [C(q,t)]A
      ENDDO                            ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_Sg  ,         p%NAug) =     DOT_PRODUCT( RtHSdat%PLinVelEZ(DOF_Sg ,0,:), RtHSdat%FrcZAllt              )        ! {-f(qd,q,t)}X + {-f(qd,q,t)}HydroX + {-f(qd,q,t)}T + {-f(qd,q,t)}AeroT + {-f(qd,q,t)}HydroT + {-f(qd,q,t)}N + {-f(qd,q,t)}R + {-f(qd,q,t)}H + {-f(qd,q,t)}B + {-f(qd,q,t)}AeroB + {-f(qd,q,t)}A + {-f(qd,q,t)}AeroA
   ENDIF

   IF ( p%DOF_Flag (DOF_Sw  ) )  THEN
      DO I = p%DOFs%Diag(DOF_Sw  ),p%DOFs%NActvDOF   ! Loop through all active (enabled) DOFs on or below the diagonal
         AugMat(p%DOFs%SrtPS(I),DOF_Sw  ) = -1.*DOT_PRODUCT( RtHSdat%PLinVelEZ(DOF_Sw ,0,:), RtHSdat%PFrcZAll (:,p%DOFs%SrtPS(I)) ) ! [C(q,t)]X + [C(q,t)]HydroX + [C(q,t)]T + [C(q,t)]HydroT + [C(q,t)]N + [C(q,t)]R + [C(q,t)]H + [C(q,t)]B + [C(q,t)]A
      ENDDO                            ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_Sw  ,         p%NAug) =     DOT_PRODUCT( RtHSdat%PLinVelEZ(DOF_Sw ,0,:), RtHSdat%FrcZAllt              )        ! {-f(qd,q,t)}X + {-f(qd,q,t)}HydroX + {-f(qd,q,t)}T + {-f(qd,q,t)}AeroT + {-f(qd,q,t)}HydroT + {-f(qd,q,t)}N + {-f(qd,q,t)}R + {-f(qd,q,t)}H + {-f(qd,q,t)}B + {-f(qd,q,t)}AeroB + {-f(qd,q,t)}A + {-f(qd,q,t)}AeroA
   ENDIF

   IF ( p%DOF_Flag (DOF_Hv  ) )  THEN
      DO I = p%DOFs%Diag(DOF_Hv  ),p%DOFs%NActvDOF   ! Loop through all active (enabled) DOFs on or below the diagonal
         AugMat(p%DOFs%SrtPS(I),DOF_Hv  ) = -1.*DOT_PRODUCT( RtHSdat%PLinVelEZ(DOF_Hv ,0,:), RtHSdat%PFrcZAll (:,p%DOFs%SrtPS(I)) ) ! [C(q,t)]X + [C(q,t)]HydroX + [C(q,t)]T + [C(q,t)]HydroT + [C(q,t)]N + [C(q,t)]R + [C(q,t)]H + [C(q,t)]B + [C(q,t)]A
      ENDDO                            ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_Hv  ,         p%NAug) =     DOT_PRODUCT( RtHSdat%PLinVelEZ(DOF_Hv ,0,:), RtHSdat%FrcZAllt              )        ! {-f(qd,q,t)}X + {-f(qd,q,t)}GravX + {-f(qd,q,t)}HydroX + {-f(qd,q,t)}T + {-f(qd,q,t)}GravT + {-f(qd,q,t)}AeroT + {-f(qd,q,t)}HydroT + {-f(qd,q,t)}N + {-f(qd,q,t)}GravN + {-f(qd,q,t)}R + {-f(qd,q,t)}GravR + {-f(qd,q,t)}H + {-f(qd,q,t)}GravH + {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB + {-f(qd,q,t)}A + {-f(qd,q,t)}GravA + {-f(qd,q,t)}AeroA
   ENDIF

   IF ( p%DOF_Flag (DOF_R   ) )  THEN
      DO I = p%DOFs%Diag(DOF_R   ),p%DOFs%NActvDOF   ! Loop through all active (enabled) DOFs on or below the diagonal

         AugMat(p%DOFs%SrtPS(I),DOF_R   ) = -1.*DOT_PRODUCT( RtHSdat%PAngVelEX(DOF_R  ,0,:), RtHSdat%PMomXAll (:,p%DOFs%SrtPS(I)) ) ! [C(q,t)]X + [C(q,t)]HydroX + [C(q,t)]T + [C(q,t)]HydroT + [C(q,t)]N + [C(q,t)]R + [C(q,t)]G + [C(q,t)]H + [C(q,t)]B + [C(q,t)]A
      ENDDO                            ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_R   ,         p%NAug) =     DOT_PRODUCT( RtHSdat%PAngVelEX(DOF_R  ,0,:), RtHSdat%MomXAllt              )        ! {-f(qd,q,t)}X + {-f(qd,q,t)}GravX + {-f(qd,q,t)}HydroX + {-f(qd,q,t)}T + {-f(qd,q,t)}GravT + {-f(qd,q,t)}AeroT + {-f(qd,q,t)}HydroT + {-f(qd,q,t)}N + {-f(qd,q,t)}GravN + {-f(qd,q,t)}R + {-f(qd,q,t)}GravR + {-f(qd,q,t)}G + {-f(qd,q,t)}H + {-f(qd,q,t)}GravH + {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB + {-f(qd,q,t)}A + {-f(qd,q,t)}GravA + {-f(qd,q,t)}AeroA
   ENDIF

   IF ( p%DOF_Flag (DOF_P   ) )  THEN
      DO I = p%DOFs%Diag(DOF_P   ),p%DOFs%NActvDOF    ! Loop through all active (enabled) DOFs on or below the diagonal
         AugMat(p%DOFs%SrtPS(I),DOF_P   ) = -1.*DOT_PRODUCT( RtHSdat%PAngVelEX(DOF_P  ,0,:), RtHSdat%PMomXAll (:,p%DOFs%SrtPS(I)) ) ! [C(q,t)]X + [C(q,t)]HydroX + [C(q,t)]T + [C(q,t)]HydroT + [C(q,t)]N + [C(q,t)]R + [C(q,t)]G + [C(q,t)]H + [C(q,t)]B + [C(q,t)]A
      ENDDO                                                             ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_P            ,p%NAug) =     DOT_PRODUCT( RtHSdat%PAngVelEX(DOF_P  ,0,:), RtHSdat%MomXAllt              )        ! {-f(qd,q,t)}X + {-f(qd,q,t)}GravX + {-f(qd,q,t)}HydroX + {-f(qd,q,t)}T + {-f(qd,q,t)}GravT + {-f(qd,q,t)}AeroT + {-f(qd,q,t)}HydroT + {-f(qd,q,t)}N + {-f(qd,q,t)}GravN + {-f(qd,q,t)}R + {-f(qd,q,t)}GravR + {-f(qd,q,t)}G + {-f(qd,q,t)}H + {-f(qd,q,t)}GravH + {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB + {-f(qd,q,t)}A + {-f(qd,q,t)}GravA + {-f(qd,q,t)}AeroA
   END IF

   IF ( p%DOF_Flag (DOF_Y   ) )  THEN
      DO I = p%DOFs%Diag(DOF_Y   ),p%DOFs%NActvDOF    ! Loop through all active (enabled) DOFs on or below the diagonal
         AugMat(p%DOFs%SrtPS(I),DOF_Y   ) = -1.*DOT_PRODUCT( RtHSdat%PAngVelEX(DOF_Y  ,0,:), RtHSdat%PMomXAll (:,p%DOFs%SrtPS(I)) ) ! [C(q,t)]X + [C(q,t)]HydroX + [C(q,t)]T + [C(q,t)]HydroT + [C(q,t)]N + [C(q,t)]R + [C(q,t)]G + [C(q,t)]H + [C(q,t)]B + [C(q,t)]A
      ENDDO                                                             ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_Y   ,         p%NAug) =     DOT_PRODUCT( RtHSdat%PAngVelEX(DOF_Y  ,0,:), RtHSdat%MomXAllt              )        ! {-f(qd,q,t)}X + {-f(qd,q,t)}GravX + {-f(qd,q,t)}HydroX + {-f(qd,q,t)}T + {-f(qd,q,t)}GravT + {-f(qd,q,t)}AeroT + {-f(qd,q,t)}HydroT + {-f(qd,q,t)}N + {-f(qd,q,t)}GravN + {-f(qd,q,t)}R + {-f(qd,q,t)}GravR + {-f(qd,q,t)}G + {-f(qd,q,t)}H + {-f(qd,q,t)}GravH + {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB + {-f(qd,q,t)}A + {-f(qd,q,t)}GravA + {-f(qd,q,t)}AeroA
   ENDIF

   IF ( p%DOF_Flag (DOF_TFA1) )  THEN
      DO I = p%DOFs%Diag(DOF_TFA1),p%DOFs%NActvDOF   ! Loop through all active (enabled) DOFs on or below the diagonal
         AugMat(p%DOFs%SrtPS(I),DOF_TFA1) = AugMat(p%DOFs%SrtPS(I),DOF_TFA1)                             &
                                          -  DOT_PRODUCT( RtHSdat%PLinVelEO(DOF_TFA1,0,:),       &
                                                          RtHSdat%PFrcONcRt(:,p%DOFs%SrtPS(I)) ) &                          ! [C(q,t)]N + [C(q,t)]R + [C(q,t)]G + [C(q,t)]H + [C(q,t)]B + [C(q,t)]A
                                          -  DOT_PRODUCT( RtHSdat%PAngVelEB(DOF_TFA1,0,:),       &
                                                          RtHSdat%PMomBNcRt(:,p%DOFs%SrtPS(I)) )
      ENDDO                            ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_TFA1,         p%NAug) = AugMat(DOF_TFA1,    p%NAug)                                  &
                                          +  DOT_PRODUCT( RtHSdat%PLinVelEO(DOF_TFA1,0,:), RtHSdat%FrcONcRtt  ) &   ! {-f(qd,q,t)}N + {-f(qd,q,t)}GravN + {-f(qd,q,t)}R + {-f(qd,q,t)}GravR + {-f(qd,q,t)}G + {-f(qd,q,t)}H + {-f(qd,q,t)}GravH + {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB + {-f(qd,q,t)}A + {-f(qd,q,t)}GravA + {-f(qd,q,t)}AeroA
                                          +  DOT_PRODUCT( RtHSdat%PAngVelEB(DOF_TFA1,0,:), RtHSdat%MomBNcRtt  )
   ENDIF

   IF ( p%DOF_Flag (DOF_TSS1) )  THEN
      DO I = p%DOFs%Diag(DOF_TSS1),p%DOFs%NActvDOF   ! Loop through all active (enabled) DOFs on or below the diagonal
         AugMat(p%DOFs%SrtPS(I),DOF_TSS1) = AugMat(p%DOFs%SrtPS(I),DOF_TSS1)                             &
                                          -  DOT_PRODUCT( RtHSdat%PLinVelEO(DOF_TSS1,0,:),       &
                                                          RtHSdat%PFrcONcRt(:,p%DOFs%SrtPS(I)) ) &                          ! [C(q,t)]N + [C(q,t)]R + [C(q,t)]G + [C(q,t)]H + [C(q,t)]B + [C(q,t)]A
                                          -  DOT_PRODUCT( RtHSdat%PAngVelEB(DOF_TSS1,0,:),       &
                                                          RtHSdat%PMomBNcRt(:,p%DOFs%SrtPS(I)) )
      ENDDO                            ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_TSS1,         p%NAug) = AugMat(DOF_TSS1,    p%NAug)                                  &
                                          +  DOT_PRODUCT( RtHSdat%PLinVelEO(DOF_TSS1,0,:), RtHSdat%FrcONcRtt  ) &   ! {-f(qd,q,t)}N + {-f(qd,q,t)}GravN + {-f(qd,q,t)}R + {-f(qd,q,t)}GravR + {-f(qd,q,t)}G + {-f(qd,q,t)}H + {-f(qd,q,t)}GravH + {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB + {-f(qd,q,t)}A + {-f(qd,q,t)}GravA + {-f(qd,q,t)}AeroA
                                          +  DOT_PRODUCT( RtHSdat%PAngVelEB(DOF_TSS1,0,:), RtHSdat%MomBNcRtt  )
   ENDIF

   IF ( p%DOF_Flag (DOF_TFA2) )  THEN
      DO I = p%DOFs%Diag(DOF_TFA2),p%DOFs%NActvDOF   ! Loop through all active (enabled) DOFs on or below the diagonal
         AugMat(p%DOFs%SrtPS(I),DOF_TFA2) = AugMat(p%DOFs%SrtPS(I),DOF_TFA2)                             &
                                          -  DOT_PRODUCT( RtHSdat%PLinVelEO(DOF_TFA2,0,:),       &
                                                          RtHSdat%PFrcONcRt(:,p%DOFs%SrtPS(I)) ) &                          ! [C(q,t)]N + [C(q,t)]R + [C(q,t)]G + [C(q,t)]H + [C(q,t)]B + [C(q,t)]A
                                          -  DOT_PRODUCT( RtHSdat%PAngVelEB(DOF_TFA2,0,:),       &
                                                          RtHSdat%PMomBNcRt(:,p%DOFs%SrtPS(I)) )
      ENDDO                            ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_TFA2,         p%NAug) = AugMat(DOF_TFA2,    p%NAug)                                  &
                                          +  DOT_PRODUCT( RtHSdat%PLinVelEO(DOF_TFA2,0,:), RtHSdat%FrcONcRtt  ) &   ! {-f(qd,q,t)}N + {-f(qd,q,t)}GravN + {-f(qd,q,t)}R + {-f(qd,q,t)}GravR + {-f(qd,q,t)}G + {-f(qd,q,t)}H + {-f(qd,q,t)}GravH + {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB + {-f(qd,q,t)}A + {-f(qd,q,t)}GravA + {-f(qd,q,t)}AeroA
                                          +  DOT_PRODUCT( RtHSdat%PAngVelEB(DOF_TFA2,0,:), RtHSdat%MomBNcRtt  )
   ENDIF

   IF ( p%DOF_Flag (DOF_TSS2) )  THEN
      DO I = p%DOFs%Diag(DOF_TSS2),p%DOFs%NActvDOF   ! Loop through all active (enabled) DOFs on or below the diagonal
         AugMat(p%DOFs%SrtPS(I),DOF_TSS2) = AugMat(p%DOFs%SrtPS(I),DOF_TSS2)                             &
                                          -  DOT_PRODUCT( RtHSdat%PLinVelEO(DOF_TSS2,0,:),       &
                                                          RtHSdat%PFrcONcRt(:,p%DOFs%SrtPS(I)) ) &                          ! [C(q,t)]N + [C(q,t)]R + [C(q,t)]G + [C(q,t)]H + [C(q,t)]B + [C(q,t)]A
                                          -  DOT_PRODUCT( RtHSdat%PAngVelEB(DOF_TSS2,0,:),       &
                                                          RtHSdat%PMomBNcRt(:,p%DOFs%SrtPS(I)) )
      ENDDO                            ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_TSS2,         p%NAug) = AugMat(DOF_TSS2,    p%NAug)                                  &
                                          +  DOT_PRODUCT( RtHSdat%PLinVelEO(DOF_TSS2,0,:), RtHSdat%FrcONcRtt  ) &   ! {-f(qd,q,t)}N + {-f(qd,q,t)}GravN + {-f(qd,q,t)}R + {-f(qd,q,t)}GravR + {-f(qd,q,t)}G + {-f(qd,q,t)}H + {-f(qd,q,t)}GravH + {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB + {-f(qd,q,t)}A + {-f(qd,q,t)}GravA + {-f(qd,q,t)}AeroA
                                          +  DOT_PRODUCT( RtHSdat%PAngVelEB(DOF_TSS2,0,:), RtHSdat%MomBNcRtt  )
   ENDIF
   
   IF ( p%DOF_Flag (DOF_Yaw ) )  THEN
      DO I = p%DOFs%Diag(DOF_Yaw ),p%DOFs%NActvDOF   ! Loop through all active (enabled) DOFs on or below the diagonal
         AugMat(p%DOFs%SrtPS(I),DOF_Yaw ) = -DOT_PRODUCT( RtHSdat%PAngVelEN(DOF_Yaw ,0,:), RtHSdat%PMomBNcRt(:,p%DOFs%SrtPS(I)) )   ! [C(q,t)]N + [C(q,t)]R + [C(q,t)]G + [C(q,t)]H + [C(q,t)]B + [C(q,t)]A
      ENDDO                            ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_Yaw ,         p%NAug) =  DOT_PRODUCT( RtHSdat%PAngVelEN(DOF_Yaw ,0,:), RtHSdat%MomBNcRtt             ) &        ! {-f(qd,q,t)}N + {-f(qd,q,t)}GravN + {-f(qd,q,t)}R + {-f(qd,q,t)}GravR + {-f(qd,q,t)}G + {-f(qd,q,t)}H + {-f(qd,q,t)}GravH + {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB + {-f(qd,q,t)}A + {-f(qd,q,t)}GravA + {-f(qd,q,t)}AeroA
                                                              + u%YawMom                                                            ! + {-f(qd,q,t)}SpringYaw  + {-f(qd,q,t)}DampYaw; NOTE: The neutral yaw rate, YawRateNeut, defaults to zero.  It is only used for yaw control.
   ENDIF
   
   
   IF ( p%DOF_Flag (DOF_RFrl) )  THEN
      DO I = p%DOFs%Diag(DOF_RFrl),p%DOFs%NActvDOF   ! Loop through all active (enabled) DOFs on or below the diagonal
         AugMat(p%DOFs%SrtPS(I),DOF_RFrl) = -DOT_PRODUCT( RtHSdat%PAngVelER(DOF_RFrl,0,:),       &
                                                          RtHSdat%PMomNGnRt(:,p%DOFs%SrtPS(I)) )                            ! [C(q,t)]R + [C(q,t)]G + [C(q,t)]H + [C(q,t)]B
      ENDDO                            ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_RFrl,         p%NAug) =  DOT_PRODUCT( RtHSdat%PAngVelER(DOF_RFrl,0,:), RtHSdat%MomNGnRtt  ) &   ! {-f(qd,q,t)}R + {-f(qd,q,t)}GravR + {-f(qd,q,t)}G + {-f(qd,q,t)}H + {-f(qd,q,t)}GravH + {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB
                                                                              +  RtHSdat%RFrlMom                            ! + {-f(qd,q,t)}SpringRF + {-f(qd,q,t)}DampRF
   ENDIF

   TmpVec = p%GenIner*CoordSys%c1*DOT_PRODUCT( CoordSys%c1, RtHSdat%PAngVelEG(DOF_GeAz,0,:) )  ! = ( generator inertia dyadic ) Dot ( partial angular velocity of G in E for DOF_GeAz )

   IF ( p%DOF_Flag (DOF_GeAz) )  THEN
      DO I = p%DOFs%Diag(DOF_GeAz),p%DOFs%NActvDOF   ! Loop through all active (enabled) DOFs on or below the diagonal
         AugMat(p%DOFs%SrtPS(I),DOF_GeAz) = -DOT_PRODUCT( RtHSdat%PAngVelEL(DOF_GeAz,0,:), RtHSdat%PMomLPRot(:,p%DOFs%SrtPS(I)) )! [C(q,t)]H + [C(q,t)]B
      ENDDO                            ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_GeAz,         p%NAug) =  DOT_PRODUCT( RtHSdat%PAngVelEL(DOF_GeAz,0,:), RtHSdat%MomLPRott             ) &     ! {-f(qd,q,t)}H + {-f(qd,q,t)}GravH + {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB
                                                              -  GBoxTrq                                                         ! + {-f(qd,q,t)}Gen + {-f(qd,q,t)}Brake


      ! The previous loop (DO I = p%DOFs%Diag(DOF_GeAz),p%DOFs%NActvDOF) misses the
      !   generator inertia-contribution to the mass matrix and forcing function.
      !   Thus, add these in as well:


         AugMat(DOF_GeAz,       DOF_GeAz) = AugMat(DOF_GeAz,DOF_GeAz)                                    &
                                            +  DOT_PRODUCT( RtHSdat%PAngVelEG(DOF_GeAz,0,:), TmpVec                )             ! [C(q,t)]G
         AugMat(DOF_GeAz,         p%NAug) = AugMat(DOF_GeAz,  p%NAug)                                    &
                                            -  DOT_PRODUCT( RtHSdat%AngAccEGt              , TmpVec                )             ! {-f(qd,q,t)}G


   ENDIF

   IF ( p%DOF_Flag (DOF_DrTr) )  THEN
      DO I = p%DOFs%Diag(DOF_DrTr),p%DOFs%NActvDOF   ! Loop through all active (enabled) DOFs on or below the diagonal
         AugMat(p%DOFs%SrtPS(I),DOF_DrTr) = -DOT_PRODUCT( RtHSdat%PAngVelEL(DOF_DrTr,0,:), RtHSdat%PMomLPRot(:,p%DOFs%SrtPS(I)) ) ! [C(q,t)]H + [C(q,t)]B
      ENDDO                            ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_DrTr,         p%NAug) =  DOT_PRODUCT( RtHSdat%PAngVelEL(DOF_DrTr,0,:), RtHSdat%MomLPRott             ) &      ! {-f(qd,q,t)}H + {-f(qd,q,t)}GravH + {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB
                                                          -  p%DTTorSpr*x%QT (DOF_DrTr)                                    &      ! + {-f(qd,q,t)}ElasticDrive
                                                          -  p%DTTorDmp*x%QDT(DOF_DrTr)                                           ! + {-f(qd,q,t)}DampDrive
   ENDIF

   IF ( p%DOF_Flag (DOF_TFrl) )  THEN
      ! The tail-furl DOF does not affect any DOF index larger than DOF_TFrl.  Therefore, there is no need to perform the loop: DO I = Diag(DOF_TFrl),NActvDOF
         AugMat(DOF_TFrl,       DOF_TFrl) = -DOT_PRODUCT( RtHSdat%PAngVelEA(DOF_TFrl,0,:), RtHSdat%PMomNTail(:,DOF_TFrl) )        ! [C(q,t)]A
         AugMat(DOF_TFrl,         p%NAug) =  DOT_PRODUCT( RtHSdat%PAngVelEA(DOF_TFrl,0,:), RtHSdat%MomNTailt             ) &      ! {-f(qd,q,t)}A + {-f(qd,q,t)}GravA + {-f(qd,q,t)}AeroA
                                                              +  RtHSdat%TFrlMom                                                  ! + {-f(qd,q,t)}SpringTF + {-f(qd,q,t)}DampTF
   ENDIF

   IF ( ( p%NumBl == 2 ) ) THEN 
      IF ( p%DOF_Flag(DOF_Teet) )  THEN  ! NOTE: two "ifs" since DOF_Teet may be out of bound
         ! The teeter DOF does not affect any DOF index larger than DOF_Teet.  Therefore, there is no need to perform the loop: DO I = Diag(DOF_Teet),NActvDOF
         AugMat(DOF_Teet,       DOF_Teet) = -DOT_PRODUCT( RtHSdat%PAngVelEH(DOF_Teet,0,:), RtHSdat%PMomLPRot(:,DOF_Teet) )        ! [C(q,t)]H + [C(q,t)]B
         AugMat(DOF_Teet,         p%NAug) =  DOT_PRODUCT( RtHSdat%PAngVelEH(DOF_Teet,0,:), RtHSdat%MomLPRott             ) &      ! {-f(qd,q,t)}H + {-f(qd,q,t)}GravH + {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB
                                                              +  RtHSdat%TeetMom                                                  ! + {-f(qd,q,t)}SpringTeet + {-f(qd,q,t)}DampTeet
      ENDIF
   ENDIF
   !..................................................................................................................................
   ! So far, we have only filled in the portions of the mass matrix on and below the diagonal.  Because the mass matrix is symmetric
   !   up to this point, let's fill in the portion above the diagonal by mirroring the values from below:
   ! NOTE: The vector subscript array, SrtPS(), used in the following loops must be sorted from smallest to largest DOF index in order
   !   for the loops to work to enter values only on and below the diagonal of AugMat():
   !..................................................................................................................................

      DO L = 2,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs above the diagonal (columns)
         DO I = 1,L-1   ! Loop through all active (enabled) DOFs above the diagonal (rows)
            AugMat(p%DOFs%SrtPS(I),p%DOFs%SrtPS(L)) = AugMat(p%DOFs%SrtPS(L),p%DOFs%SrtPS(I))
         ENDDO          ! I - All active (enabled) DOFs above the diagonal (rows)
      ENDDO             ! L - All active (enabled) DOFs above the diagonal (columns)

   ! Let's add the gearbox friction terms to the mass matrix and forcing
   !   function.  These only effect the equation for the generator azimuth DOF.
   ! NOTE: the MASS MATRIX WILL NO LONGER BE SYMMETRIC after adding these
   !       terms, unless the gearbox efficiency, GBoxEff, was set to 100%:
   
   
   GBoxEffFac2 = ( 1.0/RtHSdat%GBoxEffFac - 1.0 ) ! = ( 1 / GBoxEff^SgnPrvLSTQ - 1 )
   !TmpVec = p%GenIner*CoordSys%c1*DOT_PRODUCT( CoordSys%c1, RtHSdat%PAngVelEG(DOF_GeAz,0,:) )  ! = ( generator inertia dyadic ) Dot ( partial angular velocity of G in E for DOF_GeAz )

   DO I = 1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs

      AugMat(DOF_GeAz,p%DOFs%SrtPS(I)) = AugMat(DOF_GeAz,p%DOFs%SrtPS(I)) &                                          ! NOTE: TmpVec is still = ( generator inertia dyadic ) Dot ( partial angular velocity of G in E for DOF_GeAz ) in the following equation
                                + GBoxEffFac2*  DOT_PRODUCT( RtHSdat%PAngVelEG(p%DOFs%SrtPS(I),0,:), TmpVec )        ! [C(q,t)]GBFric

   ENDDO             ! I - All active (enabled) DOFs

   AugMat(   DOF_GeAz,    p%NAug) = AugMat(DOF_GeAz,    p%NAug) &                                                    ! NOTE: TmpVec is still = ( generator inertia dyadic ) Dot ( partial angular velocity of G in E for DOF_GeAz ) in the following equation
                                - GBoxEffFac2*( DOT_PRODUCT( RtHSdat%AngAccEGt              , TmpVec ) + GBoxTrq )   ! {-f(qd,q,t)}GBFric

   
   
   
END SUBROUTINE FillAugMat
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine allocates the arrays and meshes stored in the ED_OutputType data structure (y), based on the parameters (p). 
!! Inputs (u) are included only so that output meshes can be siblings of the inputs.
!! The routine assumes that the arrays/meshes are not currently allocated (It will produce a fatal error otherwise.)
!! Also note that this must be called after init_u() so that the misc variables that contain the orientations are set.
SUBROUTINE ED_AllocOutput( p, m, u, y, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(ED_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(ED_MiscVarType),         INTENT(IN   )  :: m           !< Misc vars (initial positions, set in Init_Inputs())
   TYPE(ED_InputType),           INTENT(INOUT)  :: u           !< Input meshes (sibling)
   TYPE(ED_OutputType),          INTENT(INOUT)  :: y           !< Outputs to be allocated
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
      
   
   ! local variables
   REAL(R8Ki)                                   :: Orientation(3,3) 
   REAL(ReKi)                                   :: Position(3) 
   INTEGER(IntKi)                               :: NodeNum     ! node number
   INTEGER(IntKi)                               :: J, K        ! loop counters
   INTEGER(IntKi)                               :: ErrStat2    ! The error identifier (ErrStat)
   CHARACTER(ErrMsgLen)                         :: ErrMsg2     ! The error message (ErrMsg)
   
   
      ! initialize variables:
      
   ErrStat = ErrID_None
   ErrMsg  = ""
      

   CALL AllocAry( y%WriteOutput, p%numOuts + p%BldNd_TotNumOuts, 'WriteOutput', errStat2, errMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN

   CALL AllocAry( y%BlPitch, p%NumBl, 'BlPitch', ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN
      
   !.......................................................
   ! Create Line2 Mesh for motion outputs on blades:
   !.......................................................
   IF ( .NOT. p%BD4Blades) THEN
      ALLOCATE( y%BladeLn2Mesh(p%NumBl), Stat=ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN
         CALL CheckError( ErrID_Fatal, 'ED: Could not allocate space for y%BladeLn2Mesh{p%NumBl}' )
         RETURN
      END IF
   
      DO K = 1,p%NumBl
         
         CALL MeshCreate( BlankMesh          = y%BladeLn2Mesh(K)      &
                           , NNodes          = p%BldNodes+2           &
                           , IOS             = COMPONENT_OUTPUT       &
                           , TranslationDisp = .TRUE.                 &
                           , Orientation     = .TRUE.                 &
                           , RotationVel     = .TRUE.                 &
                           , TranslationVel  = .TRUE.                 &
                           , RotationAcc     = .TRUE.                 &
                           , TranslationAcc  = .TRUE.                 &
                           , ErrStat         = ErrStat2               &
                           , ErrMess         = ErrMsg2                )
            CALL CheckError( ErrStat2, ErrMsg2 )
            IF (ErrStat >= AbortErrLev) RETURN
      
         DO J = 1,p%BldNodes               
            CALL MeshPositionNode ( y%BladeLn2Mesh(K), J, u%BladePtLoads(K)%Position(:,J), ErrStat2, ErrMsg2, Orient=u%BladePtLoads(K)%RefOrientation(:,:,J) )
               CALL CheckError( ErrStat2, ErrMsg2 )
               IF (ErrStat >= AbortErrLev) RETURN
         END DO
            
            ! now add position/orientation of nodes for AD14 or AD15
         if (p%UseAD14) then     ! position/orientation of nodes for AeroDyn v14 or v15    
         
               ! Use orientation at p%BldNodes for the extra node at the blade tip
            CALL MeshPositionNode ( y%BladeLn2Mesh(K), p%BldNodes + 1, (/0.0_ReKi, 0.0_ReKi, p%BldFlexL /), ErrStat2, ErrMsg2, Orient=u%BladePtLoads(K)%RefOrientation(:,:,p%BldNodes) )
               CALL CheckError( ErrStat2, ErrMsg2 )
               IF (ErrStat >= AbortErrLev) RETURN
            
               ! Use orientation at node 1 for the blade root            
            CALL MeshPositionNode ( y%BladeLn2Mesh(K), p%BldNodes + 2, (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi /), ErrStat2, ErrMsg2, Orient=u%BladePtLoads(K)%RefOrientation(:,:,1), ref=.true. )
               CALL CheckError( ErrStat2, ErrMsg2 )
               IF (ErrStat >= AbortErrLev) RETURN
               
         else
         
            ! position the nodes on the blade root and blade tip:
            DO J = 0,p%TipNode,p%TipNode
               if (j==0) then ! blade root
                  NodeNum = p%BldNodes + 2
                  y%BladeLn2Mesh(K)%RefNode = NodeNum
               elseif (j==p%TipNode) then ! blade tip
                  NodeNum = p%BldNodes + 1
               end if
         
               Orientation(1,1) =     m%CoordSys%n1(K,J,1)
               Orientation(2,1) =     m%CoordSys%n2(K,J,1)
               Orientation(3,1) =     m%CoordSys%n3(K,J,1)
               Orientation(1,2) = -1.*m%CoordSys%n1(K,J,3)
               Orientation(2,2) = -1.*m%CoordSys%n2(K,J,3)
               Orientation(3,2) = -1.*m%CoordSys%n3(K,J,3)
               Orientation(1,3) =     m%CoordSys%n1(K,J,2)
               Orientation(2,3) =     m%CoordSys%n2(K,J,2)
               Orientation(3,3) =     m%CoordSys%n3(K,J,2) 
               
                  ! Translational Displacement 
               position(1) =     m%RtHS%rS (1,K,J)                ! = the distance from the undeflected tower centerline to the current blade node in the xi ( z1) direction
               position(2) = -1.*m%RtHS%rS (3,K,J)                ! = the distance from the undeflected tower centerline to the current blade node in the yi (-z3) direction
               position(3) =     m%RtHS%rS (2,K,J)  + p%PtfmRefzt ! = the distance from the nominal tower base position (i.e., the undeflected position of the tower base) to the current blade node in the zi ( z2) direction
               
               
               CALL MeshPositionNode ( y%BladeLn2Mesh(K), NodeNum, position, ErrStat2, ErrMsg2, Orient=Orientation )
                  CALL CheckError( ErrStat2, ErrMsg2 )
                  IF (ErrStat >= AbortErrLev) RETURN
                                    
            END DO ! nodes 
            
         end if ! position/orientation of nodes for AeroDyn v14 or v15
         
         ! create elements:      
         DO J = 2,p%TipNode !p%BldNodes + 1
            
            CALL MeshConstructElement ( Mesh      = y%BladeLn2Mesh(K)  &
                                       , Xelement = ELEMENT_LINE2      &
                                       , P1       = J-1                &   ! node1 number
                                       , P2       = J                  &   ! node2 number
                                       , ErrStat  = ErrStat2           &
                                       , ErrMess  = ErrMsg2            )
               CALL CheckError( ErrStat2, ErrMsg2 )
               IF (ErrStat >= AbortErrLev) RETURN
      
         END DO ! J (blade nodes)

            ! add the other extra element, connecting the first node on the blade:
         CALL MeshConstructElement ( Mesh      = y%BladeLn2Mesh(K)  &
                                    , Xelement = ELEMENT_LINE2      &
                                    , P1       = p%BldNodes + 2     &   ! node1 number (extra node at root)
                                    , P2       = 1                  &   ! node2 number (first node on blade)
                                    , ErrStat  = ErrStat2           &
                                    , ErrMess  = ErrMsg2            )         
            CALL CheckError( ErrStat2, ErrMsg2 )
            IF (ErrStat >= AbortErrLev) RETURN
      
         
            ! that's our entire mesh:
         CALL MeshCommit ( y%BladeLn2Mesh(K), ErrStat2, ErrMsg2 )   
            CALL CheckError( ErrStat2, ErrMsg2 )
            IF (ErrStat >= AbortErrLev) RETURN
   
      END DO
      
   END IF
   
   !.......................................................
   ! Create Point Mesh for Motions Output at Platform Reference Point:
   !.......................................................
      
   CALL MeshCopy ( SrcMesh  = u%PlatformPtMesh &
                 , DestMesh = y%PlatformPtMesh &
                 , CtrlCode = MESH_SIBLING     &
                 , IOS      = COMPONENT_OUTPUT &
                 , TranslationDisp = .TRUE.    &
                 , Orientation     = .TRUE.    &
                 , RotationVel     = .TRUE.    &
                 , TranslationVel  = .TRUE.    &
                 , RotationAcc     = .TRUE.    &
                 , TranslationAcc  = .TRUE.    &
                 , ErrStat  = ErrStat2         &
                 , ErrMess  = ErrMsg2          )  ! automatically sets    y%PlatformPtMesh%RemapFlag = .TRUE.
   
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN
   
   
   !.......................................................
   ! Create Line2 Mesh for Motions Output on Tower Line2 Mesh:
   !  first p%TwrNodes nodes are the same as the input TowerPtLoads mesh
   !.......................................................
      
   CALL MeshCreate( BlankMesh = y%TowerLn2Mesh           &
                    , IOS             = COMPONENT_OUTPUT &
                    , NNodes          = p%TwrNodes + 2   &
                    , TranslationDisp = .TRUE.           &
                    , Orientation     = .TRUE.           &
                    , RotationVel     = .TRUE.           &
                    , TranslationVel  = .TRUE.           &  
                    , RotationAcc     = .TRUE.           &  
                    , TranslationAcc  = .TRUE.           &
                    , ErrStat         = ErrStat2         &
                    , ErrMess         = ErrMsg2          )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF (ErrStat >= AbortErrLev) RETURN
   
      ! position the nodes on the tower:
      DO J = 1,p%TwrNodes      
         CALL MeshPositionNode ( y%TowerLn2Mesh, J, u%TowerPtLoads%Position(:,J), ErrStat2, ErrMsg2, &
                                 orient = u%TowerPtLoads%RefOrientation(:,:,J) )
            CALL CheckError(ErrStat2,ErrMsg2)
            IF (ErrStat >= AbortErrLev) RETURN
      END DO

   ! for now, we're going to add two nodes, one at the beginning and the other at the end
   ! they're numbered this way so that I don't have to redo all the loops in the computational part of the code
   ! I am not going to use them in my input.
      CALL MeshPositionNode ( y%TowerLn2Mesh, p%TwrNodes + 1, (/0.0_ReKi, 0.0_ReKi, p%TowerHt /), ErrStat2, ErrMsg2 ) 
         CALL CheckError(ErrStat2,ErrMsg2)
         IF (ErrStat >= AbortErrLev) RETURN
         
      CALL MeshPositionNode ( y%TowerLn2Mesh, p%TwrNodes + 2, (/0.0_ReKi, 0.0_ReKi, p%TowerBsHt  /), ErrStat2, ErrMsg2, ref=.true. )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF (ErrStat >= AbortErrLev) RETURN

            
     ! create line2 elements from the tower nodes (we've alread checked that p%TwrNodes is a positive number):
      DO J = 2,p%TwrNodes+1  !the plus 1 includes one of the end nodes
         CALL MeshConstructElement ( Mesh      = y%TowerLn2Mesh     &
                                    , Xelement = ELEMENT_LINE2      &
                                    , P1       = J-1                &   ! node1 number
                                    , P2       = J                  &   ! node2 number
                                    , ErrStat  = ErrStat2           &
                                    , ErrMess  = ErrMsg2            )
         
         CALL CheckError(ErrStat2,ErrMsg2)
         IF (ErrStat >= AbortErrLev) RETURN
      END DO
      
   ! add the other extra element, connecting the first node:
      CALL MeshConstructElement ( Mesh      = y%TowerLn2Mesh     &
                                 , Xelement = ELEMENT_LINE2      &
                                 , P1       = p%TwrNodes + 2     &   ! node1 number
                                 , P2       = 1                  &   ! node2 number
                                 , ErrStat  = ErrStat2           &
                                 , ErrMess  = ErrMsg2            )
         
      CALL CheckError(ErrStat2,ErrMsg2)
      IF (ErrStat >= AbortErrLev) RETURN
                                          
   
      ! that's our entire mesh:
   CALL MeshCommit ( y%TowerLn2Mesh, ErrStat2, ErrMsg2 )   
      CALL CheckError(ErrStat2,ErrMsg2)
      IF (ErrStat >= AbortErrLev) RETURN
                                                          
   !.......................................................
   ! Create Point Meshes for motions AeroDyn/BeamDyn needs:
   !.......................................................
   
   ! -------------- Hub -----------------------------------
      !BJJ: sibling of u%HubPtLoad
   CALL MeshCopy (     SrcMesh  = u%HubPtLoad             &
                     , DestMesh = y%HubPtMotion           &
                     , CtrlCode = MESH_SIBLING            &
                     , IOS      = COMPONENT_OUTPUT        &      
                     , TranslationDisp = .TRUE.           &
                     , Orientation     = .TRUE.           &
                     , RotationVel     = .TRUE.           &
                     ,ErrStat          = ErrStat2         &
                     ,ErrMess          = ErrMsg2          )
      CALL CheckError(ErrStat2,ErrMsg2)
      IF (ErrStat >= AbortErrLev) RETURN
      
   ! -------------- pseudo-Hub (for AD v14)  -----------------------------------
   CALL MeshCreate( BlankMesh          = y%HubPtMotion14  &
                     ,IOS              = COMPONENT_OUTPUT &
                     ,NNodes           = 1                &
                     , TranslationDisp = .TRUE.           &
                     , Orientation     = .TRUE.           &
                     , RotationVel     = .TRUE.           &
                     , ErrStat         = ErrStat2         &
                     , ErrMess         = ErrMsg2          )      
      CALL CheckError(ErrStat2,ErrMsg2)
      IF (ErrStat >= AbortErrLev) RETURN
      
      ! pseudo-Hub position and orientation (relative here as before, but should not be)
      
   CALL MeshPositionNode ( y%HubPtMotion14, 1, (/0.0_ReKi, 0.0_ReKi, p%HubHt /), ErrStat, ErrMsg ) !orientation is identity by default
      CALL CheckError(ErrStat2,ErrMsg2)
      IF (ErrStat >= AbortErrLev) RETURN
            
   CALL CommitPointMesh( y%HubPtMotion14 )
      IF (ErrStat >= AbortErrLev) RETURN
 
      
   ! -------------- Blade Roots -----------------------------------
   ALLOCATE( y%BladeRootMotion(p%NumBl), Stat=ErrStat2 )
   IF ( ErrStat2 /= 0 ) THEN
      CALL CheckError( ErrID_Fatal, 'ED: Could not allocate space for y%BladeRootMotions{p%NumBl}' )
      RETURN
   END IF
                  
      
   DO k=1,p%NumBl      
      CALL MeshCreate( BlankMesh       = y%BladeRootMotion(k)   &
                     ,IOS              = COMPONENT_OUTPUT       &
                     ,NNodes           = 1                      &
                     ,TranslationDisp  = .TRUE.                 & 
                     ,Orientation      = .TRUE.                 & 
                     ,TranslationVel   = .TRUE.                 & 
                     ,TranslationAcc   = .TRUE.                 & 
                     ,RotationVel      = .TRUE.                 & 
                     ,RotationAcc      = .TRUE.                 & 
                     ,ErrStat          = ErrStat2               &
                     ,ErrMess          = ErrMsg2                )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF (ErrStat >= AbortErrLev) RETURN            
   END DO
   
      
   CALL MeshCreate( BlankMesh          = y%BladeRootMotion14    &
                     ,IOS              = COMPONENT_OUTPUT       &
                     ,NNodes           = p%NumBl                &
                     , Orientation     = .TRUE.                 &
                     ,ErrStat          = ErrStat2               &
                     ,ErrMess          = ErrMsg2                )
      CALL CheckError(ErrStat2,ErrMsg2)
      IF (ErrStat >= AbortErrLev) RETURN

   DO K=1,p%NumBl      
      
      Orientation(1,1) =               p%CosPreC(K)
      Orientation(2,1) =  0.0_R8Ki
      Orientation(3,1) =  1.0_R8Ki *   p%SinPreC(K)

      Orientation(1,2) =  0.0_R8Ki
      Orientation(2,2) =  1.0_R8Ki
      Orientation(3,2) =  0.0_R8Ki

      Orientation(1,3) = -1.0_R8Ki *    p%SinPreC(K)
      Orientation(2,3) =  0.0_R8Ki
      Orientation(3,3) =                p%CosPreC(K)
                  
      Position(1) = p%HubRad*p%SinPreC(K)
      Position(2) = 0.0_ReKi
      Position(3) = p%HubRad*p%CosPreC(K)      
      
      CALL MeshPositionNode ( y%BladeRootMotion14, K, Position, &
                            ErrStat, ErrMsg, Orient=Orientation ) 
         CALL CheckError(ErrStat2,ErrMsg2)
         IF (ErrStat >= AbortErrLev) RETURN
                  
               
      position(1) =     m%RtHS%rS (1,K,0)                ! = the distance from the undeflected tower centerline to the current blade node in the xi ( z1) direction
      position(2) = -1.*m%RtHS%rS (3,K,0)                ! = the distance from the undeflected tower centerline to the current blade node in the yi (-z3) direction
      position(3) =     m%RtHS%rS (2,K,0)  + p%PtfmRefzt ! = the distance from the nominal tower base position (i.e., the undeflected position of the tower base) to the current blade node in the zi ( z2) direction
      
      
      Orientation(1,1) =     m%CoordSys%j1(K,1)
      Orientation(2,1) =     m%CoordSys%j2(K,1)
      Orientation(3,1) =     m%CoordSys%j3(K,1)
      Orientation(1,2) = -1.*m%CoordSys%j1(K,3)
      Orientation(2,2) = -1.*m%CoordSys%j2(K,3)
      Orientation(3,2) = -1.*m%CoordSys%j3(K,3)
      Orientation(1,3) =     m%CoordSys%j1(K,2)
      Orientation(2,3) =     m%CoordSys%j2(K,2)
      Orientation(3,3) =     m%CoordSys%j3(K,2)
      
      CALL MeshPositionNode ( y%BladeRootMotion(K), 1, Position, &
                            ErrStat, ErrMsg, Orient=Orientation ) 
         CALL CheckError(ErrStat2,ErrMsg2)
         IF (ErrStat >= AbortErrLev) RETURN
         
   END DO
                     
   CALL CommitPointMesh( y%BladeRootMotion14 )
      IF (ErrStat >= AbortErrLev) RETURN
   
   DO k=1,p%NumBl      
      CALL CommitPointMesh( y%BladeRootMotion(K) )
         IF (ErrStat >= AbortErrLev) RETURN
   END DO
   
      
   ! -------------- Rotor Furl -----------------------------------
   CALL MeshCreate( BlankMesh          = y%RotorFurlMotion14    &
                     ,IOS              = COMPONENT_OUTPUT       &
                     ,NNodes           = 1                      &
                     , TranslationDisp = .TRUE.                 &
                     , Orientation     = .TRUE.                 &
                     , RotationVel     = .TRUE.                 &
                     ,ErrStat          = ErrStat2               &
                     ,ErrMess          = ErrMsg2                )
      CALL CheckError(ErrStat2,ErrMsg2)
      IF (ErrStat >= AbortErrLev) RETURN

!bjj: FIX THIS>>>>     
!call wrscr(newline//'fix RotorFurlMotion initialization')
   CALL MeshPositionNode ( y%RotorFurlMotion14, 1, (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi /), ErrStat, ErrMsg ) !orientation is identity by default
!<<<<<FIX THIS
      CALL CheckError(ErrStat2,ErrMsg2)
      IF (ErrStat >= AbortErrLev) RETURN
      
   CALL CommitPointMesh( y%RotorFurlMotion14 )
      IF (ErrStat >= AbortErrLev) RETURN      
      
   ! -------------- Nacelle -----------------------------------      
   CALL MeshCopy ( SrcMesh  = u%NacelleLoads   &
                 , DestMesh = y%NacelleMotion  &
                 , CtrlCode = MESH_SIBLING     &
                 , IOS      = COMPONENT_OUTPUT &
                 , TranslationDisp = .TRUE.    &
                 , Orientation     = .TRUE.    &
                 , TranslationVel  = .TRUE.    &
                 , RotationVel     = .TRUE.    &
                 , TranslationAcc  = .TRUE.    &
                 , RotationAcc     = .TRUE.    &   
                 , ErrStat         = ErrStat2  &
                 , ErrMess         = ErrMsg2   )      ! automatically sets    y%NacelleMotion%RemapFlag   = .TRUE.
   
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN
      
      
   ! -------------- Tailfin -----------------------------------
   call MeshCopy ( SrcMesh  = u%TFinCMLoads    &
                 , DestMesh = y%TFinCMMotion   &
                 , CtrlCode = MESH_SIBLING     &
                 , IOS      = COMPONENT_OUTPUT &
                 , TranslationDisp = .TRUE.    &
                 , Orientation     = .TRUE.    &
                 , TranslationVel  = .TRUE.    &
                 , RotationVel     = .TRUE.    &
                 , TranslationAcc  = .TRUE.    &
                 , RotationAcc     = .TRUE.    &   
                 , ErrStat  = ErrStat2         &
                 , ErrMess  = ErrMsg2          )

   call CheckError( ErrStat2, ErrMsg2 )
   if (ErrStat >= AbortErrLev) RETURN         
     
   ! -------------- Tower Base-----------------------------------
   CALL MeshCreate( BlankMesh          = y%TowerBaseMotion14    &
                     ,IOS              = COMPONENT_OUTPUT       &
                     ,NNodes           = 1                      &
                     , TranslationDisp = .TRUE.                 &
                     , RotationVel     = .TRUE.                 &
                     ,ErrStat          = ErrStat2               &
                     ,ErrMess          = ErrMsg2                )
      CALL CheckError(ErrStat2,ErrMsg2)
      IF (ErrStat >= AbortErrLev) RETURN

!bjj: FIX THIS>>>>     
!call wrscr(newline//'fix TowerBaseMotion14 initialization')
   CALL MeshPositionNode ( y%TowerBaseMotion14, 1, (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi /), ErrStat, ErrMsg ) !orientation is identity by default
!<<<<<FIX THIS
      CALL CheckError(ErrStat2,ErrMsg2)
      IF (ErrStat >= AbortErrLev) RETURN
      
   CALL CommitPointMesh( y%TowerBaseMotion14 )
      IF (ErrStat >= AbortErrLev) RETURN      
      
      
CONTAINS
   !...............................................................................................................................
   SUBROUTINE CommitPointMesh(NewMesh)
      ! This routine makes every node a point element and then commits the mesh.
      
      TYPE(MeshType), INTENT(INOUT)  :: NewMesh  
      
      INTEGER(IntKi) :: Node
      
      DO Node = 1,NewMesh%Nnodes
         
            ! create an element from this point      
         CALL MeshConstructElement ( Mesh = NewMesh                 &
                                    , Xelement = ELEMENT_POINT      &
                                    , P1       = Node               &   ! node number
                                    , ErrStat  = ErrStat            &
                                    , ErrMess  = ErrMsg             )
            CALL CheckError(ErrStat2,ErrMsg2)
            IF (ErrStat >= AbortErrLev) RETURN

      END DO
      
         ! that's our entire mesh:
      CALL MeshCommit ( NewMesh, ErrStat2, ErrMsg2 )   
         CALL CheckError(ErrStat2,ErrMsg2)
         IF (ErrStat >= AbortErrLev) RETURN                

   END SUBROUTINE CommitPointMesh
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................
         ! Passed arguments
         
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)      
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'ED_AllocOutput:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................

      END IF

   END SUBROUTINE CheckError 
   !...............................................................................................................................
END SUBROUTINE ED_AllocOutput
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine allocates the arrays stored in the ED_InputType data structure (u), based on the parameters (p). 
!! The routine assumes that the arrays are not currently allocated (It will produce a fatal error otherwise.) It also initializes the inputs
SUBROUTINE Init_u( u, p, x, InputFileData, m, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(ED_InputType),           INTENT(INOUT)  :: u                 !< Inputs to be allocated
   TYPE(ED_ParameterType),       INTENT(IN   )  :: p                 !< Parameters
   TYPE(ED_ContinuousStateType), INTENT(IN   )  :: x                 !< Continuous states
   TYPE(ED_InputFile),           INTENT(IN   )  :: InputFileData     !< Data stored in the module's input file
   TYPE(ED_MiscVarType),         INTENT(INOUT)  :: m                 !< Misc variables (used to calculate initial position/orientation for meshes)
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat           !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None
      
   
   ! local variables
   REAL(R8Ki)                                   :: Orientation(3,3)  ! reference orientation matrix
   REAL(ReKi)                                   :: Position(3)       ! position vector
   TYPE(ED_ContinuousStateType)                 :: x_tmp             ! continuous states (set to 0)
   INTEGER(IntKi)                               :: J, K              ! loop counters
   INTEGER(IntKi)                               :: NodeNum           ! number of current blade node
   INTEGER(IntKi)                               :: ErrStat2          ! The error identifier (ErrStat)
   CHARACTER(ErrMsgLen)                         :: ErrMsg2           ! The error message (ErrMsg)
   CHARACTER(*), PARAMETER                      :: RoutineName = 'Init_u'
   
      ! initialize variables:
      
   ErrStat = ErrID_None
   ErrMsg  = ""

   !.......................................................
   ! allocate the u%BlPitchCom array    
   !.......................................................

   CALL AllocAry( u%BlPitchCom, p%NumBl, 'BlPitchCom', ErrStat2, ErrMsg2 )
   if (Failed()) return
   ! will initialize u%BlPitchCom later, after getting undisplaced positions    
   
   !.......................................................
   ! we're going to calculate the non-displaced positions of
   ! several variables so we can set up meshes properly later.
   ! want inputs and states initialized to 0 first.
   !.......................................................
   CALL ED_CopyContState( x, x_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
   if (Failed()) return
      x_tmp%qt  = 0.0_ReKi
      x_tmp%qdt = 0.0_ReKi
      x_tmp%QT (DOF_GeAz) = - p%AzimB1Up - REAL(Piby2_D, R8Ki)
         CALL Zero2TwoPi( x_tmp%QT (DOF_GeAz) )

      u%BlPitchCom = 0.0_ReKi
      
      ! set the coordinate system variables:
   CALL SetCoordSy( -p%DT, m%CoordSys, m%RtHS, u%BlPitchCom, p, x_tmp, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
   CALL CalculatePositions( p, x_tmp, m%CoordSys, m%RtHS ) ! calculate positions
   
   CALL ED_DestroyContState(x_tmp, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   !.......................................................
   ! initialize the u%BlPitchCom array    
   !.......................................................
      ! was allocated above to call SetCoordSy and CalculatePositions with undisplaced values
   u%BlPitchCom = InputFileData%BlPitch(1:p%NumBl)

   
   !.......................................................
   ! Create Line2 Meshes for loads input on blades:
   !.......................................................
      
   IF (.not. p%BD4Blades) THEN
      ALLOCATE( u%BladePtLoads(p%NumBl), STAT=ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN
         CALL SetErrStat( ErrID_Fatal, "Could not allocate u%BladePtLoads", ErrStat, ErrMsg, RoutineName )
         CALL Cleanup()
         RETURN
      END IF
   
   
      DO K=1,p%NumBl
      
         CALL MeshCreate( BlankMesh         = u%BladePtLoads(K)      &
                           ,IOS             = COMPONENT_INPUT        &
                           ,NNodes          = p%BldNodes             &
                           ,Force           = .TRUE.                 &
                           ,Moment          = .TRUE.                 &
                           ,ErrStat         = ErrStat2               &
                           ,ErrMess         = ErrMsg2                )
         if (Failed()) return
      
         if (p%UseAD14) then
            ! position the nodes on the blades:
            DO J = 1,p%BldNodes
         
               NodeNum = J
         
               Orientation(1,1) =  p%CAeroTwst(J)
               Orientation(2,1) =  p%SAeroTwst(J)
               Orientation(3,1) =  0.0_ReKi

               Orientation(1,2) = -p%SAeroTwst(J)
               Orientation(2,2) =  p%CAeroTwst(J)
               Orientation(3,2) =  0.0_ReKi

               Orientation(1,3) =  0.0_ReKi
               Orientation(2,3) =  0.0_ReKi
               Orientation(3,3) =  1.0_ReKi
                           
               CALL MeshPositionNode ( u%BladePtLoads(K), NodeNum, (/0.0_ReKi, 0.0_ReKi, p%RNodes(J) /), ErrStat2, ErrMsg2, Orient=Orientation )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  IF (ErrStat >= AbortErrLev) THEN
                     CALL Cleanup()
                     RETURN
                  END IF
                                                               
            END DO ! nodes  
         else
            ! position the nodes on the blades:
            DO J = 1,p%BldNodes
               NodeNum = J
         
               Orientation(1,1) =     m%CoordSys%n1(K,J,1)
               Orientation(2,1) =     m%CoordSys%n2(K,J,1)
               Orientation(3,1) =     m%CoordSys%n3(K,J,1)
               Orientation(1,2) = -1.*m%CoordSys%n1(K,J,3)
               Orientation(2,2) = -1.*m%CoordSys%n2(K,J,3)
               Orientation(3,2) = -1.*m%CoordSys%n3(K,J,3)
               Orientation(1,3) =     m%CoordSys%n1(K,J,2)
               Orientation(2,3) =     m%CoordSys%n2(K,J,2)
               Orientation(3,3) =     m%CoordSys%n3(K,J,2) 
               
                  ! Translational Displacement 
               position(1) =     m%RtHS%rS (1,K,J)                ! = the distance from the undeflected tower centerline to the current blade node in the xi ( z1) direction
               position(2) = -1.*m%RtHS%rS (3,K,J)                ! = the distance from the undeflected tower centerline to the current blade node in the yi (-z3) direction
               position(3) =     m%RtHS%rS (2,K,J)  + p%PtfmRefzt ! = the distance from the nominal tower base position (i.e., the undeflected position of the tower base) to the current blade node in the zi ( z2) direction
               
               
               CALL MeshPositionNode ( u%BladePtLoads(K), NodeNum, position, ErrStat2, ErrMsg2, Orient=Orientation )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  IF (ErrStat >= AbortErrLev) THEN
                     CALL Cleanup()
                     RETURN
                  END IF                           
                                    
            END DO ! nodes              
         end if ! position/orientation of nodes for AeroDyn v14 or v15
         
         ! create elements:      
         DO J = 1,p%BldNodes !p%BldNodes + 1
            
            CALL MeshConstructElement ( Mesh      = u%BladePtLoads(K)  &
                                       , Xelement = ELEMENT_POINT      &
                                       , P1       = J                  &   ! node1 number
                                       , ErrStat  = ErrStat2           &
                                       , ErrMess  = ErrMsg2            )
         
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF (ErrStat >= AbortErrLev) THEN
                  CALL Cleanup()
                  RETURN
               END IF
      
         END DO ! J (blade nodes)

            ! that's our entire mesh:
         CALL MeshCommit ( u%BladePtLoads(K), ErrStat2, ErrMsg2 )   
         if (Failed()) return

   
            ! initialize it
         u%BladePtLoads(K)%Moment   = 0.0_ReKi
         u%BladePtLoads(K)%Force    = 0.0_ReKi         
         
                     
      END DO ! blades
   END IF ! p%BD4Blades
   
                     
   !.......................................................
   ! Create Point Mesh for loads input at hub point (from BeamDyn):
   !....................................................... 
   ! place single node at hub; position affects mapping/coupling with other modules      
   Position(1)  =     m%RtHS%rQ(1)
   Position(2)  = -1.*m%RtHS%rQ(3)
   Position(3)  =     m%RtHS%rQ(2) + p%PtfmRefzt
   
   Orientation(1,1) =     m%CoordSys%g1(1)
   Orientation(2,1) =     m%CoordSys%g2(1)
   Orientation(3,1) =     m%CoordSys%g3(1)
   Orientation(1,2) = -1.*m%CoordSys%g1(3)
   Orientation(2,2) = -1.*m%CoordSys%g2(3)
   Orientation(3,2) = -1.*m%CoordSys%g3(3)
   Orientation(1,3) =     m%CoordSys%g1(2)
   Orientation(2,3) =     m%CoordSys%g2(2)
   Orientation(3,3) =     m%CoordSys%g3(2) 
   call CreatePointMesh(u%HubPtLoad, Position, Orientation, errStat2, errMsg2, hasMotion=.False., hasLoads=.True.)
   if (Failed()) return
         
                     
   !.......................................................
   ! Create Point Mesh for loads input at Platform Reference Point:
   !.......................................................
   Position = (/0.0_ReKi, 0.0_ReKi, p%PtfmRefzt /)
   call Eye(Orientation, ErrStat2, errMsg2)
   call CreatePointMesh(u%PlatformPtMesh, Position, Orientation, errStat2, errMsg2, hasMotion=.False., hasLoads=.True.)
   if (Failed()) return
      
   !.......................................................
   ! Create Point Mesh for loads input at nacelle:
   !.......................................................
   Position = (/0.0_ReKi, 0.0_ReKi, p%TowerHt /)
   call Eye(Orientation, ErrStat2, errMsg2)
   call CreatePointMesh(u%NacelleLoads, Position, Orientation, errStat2, errMsg2, hasMotion=.False., hasLoads=.True.)
   if (Failed()) return
      
   !.......................................................
   ! Create Point Mesh for loads on Rotor tailfin:
   !.......................................................
   Position(1) =     m%RtHS%rJ(1)               ! undeflected position of the tailfin CM in the xi ( z1) direction
   Position(2) = -1.*m%RtHS%rJ(3)               ! undeflected position of the tailfin CM in the yi (-z3) direction
   Position(3) =     m%RtHS%rJ(2) + p%PtfmRefzt ! undeflected position of the tailfin CM in the zi ( z2) direction
   Orientation(1,1) =     m%CoordSys%tf1(1)
   Orientation(2,1) =     m%CoordSys%tf2(1)
   Orientation(3,1) =     m%CoordSys%tf3(1)
   Orientation(1,2) = -1.*m%CoordSys%tf1(3)
   Orientation(2,2) = -1.*m%CoordSys%tf2(3)
   Orientation(3,2) = -1.*m%CoordSys%tf3(3)
   Orientation(1,3) =     m%CoordSys%tf1(2)
   Orientation(2,3) =     m%CoordSys%tf2(2)
   Orientation(3,3) =     m%CoordSys%tf3(2) 
   call CreatePointMesh(u%TFinCMLoads, Position, Orientation, errStat2, errMsg2, hasMotion=.False., hasLoads=.True.)
   if (Failed()) return


   !.......................................................
   ! Create u%TwrAddedMass for loads input on tower:
   ! SHOULD REMOVE EVENTUALLY
   !.......................................................
         
   CALL AllocAry( u%TwrAddedMass,  6_IntKi, 6_IntKi, p%TwrNodes,   'TwrAddedMass',    ErrStat2, ErrMsg2 )
   if (Failed()) return
      
      ! initialize it
   u%TwrAddedMass          = 0.0_ReKi  
      
   !.......................................................
   ! Create point Mesh for lumped load input on tower:
   ! note that this does not contain the end points
   !.......................................................
      
   CALL MeshCreate( BlankMesh      = u%TowerPtLoads         &
                     ,IOS          = COMPONENT_INPUT        &
                     ,NNodes       = p%TwrNodes             &
                     ,Force        = .TRUE.                 &
                     ,Moment       = .TRUE.                 &
                     ,ErrStat      = ErrStat2               &
                     ,ErrMess      = ErrMsg2                )
   if (Failed()) return
   
      ! position the nodes on the tower:
   DO J = 1,p%TwrNodes      
      CALL MeshPositionNode ( u%TowerPtLoads, J, (/0.0_ReKi, 0.0_ReKi, p%HNodes(J) + p%TowerBsHt /), ErrStat2, ErrMsg2 )
      if (Failed()) return
   END DO
   
      ! create elements:      
   DO J = 1,p%TwrNodes
      CALL MeshConstructElement ( Mesh      = u%TowerPtLoads     &
                                 , Xelement = ELEMENT_POINT      &
                                 , P1       = J                  &   ! node1 number
                                 , ErrStat  = ErrStat2           &
                                 , ErrMess  = ErrMsg2            )
         
      if (Failed()) return
   END DO
      
   
      ! that's our entire mesh:
   CALL MeshCommit ( u%TowerPtLoads, ErrStat2, ErrMsg2 )   
   if (Failed()) return
      
      ! initialize fields
   u%TowerPtLoads%Moment   = 0.0_ReKi
   u%TowerPtLoads%Force    = 0.0_ReKi   
         
   !.......................................................
   ! initialize all remaining inputs (non-allocatable):
   !.......................................................
         
   u%PtfmAddedMass         = 0.0_ReKi      
   u%YawMom                = 0.0_ReKi
   u%GenTrq                = 0.0_ReKi
   u%HSSBrTrqC             = 0.0_ReKi      
         
   CALL Cleanup()
   
CONTAINS
   !...............................................................................................................................
   SUBROUTINE Cleanup()
   ! This subroutine cleans up if the error for returning to calling routine
   !...............................................................................................................................      
         !.........................................................................................................................
         ! close files, deallocate local arrays
         !.........................................................................................................................
         CALL ED_DestroyContState( x_tmp, ErrStat2, ErrMsg2 )
         
   END SUBROUTINE Cleanup   
   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call Cleanup()
   end function Failed
            
END SUBROUTINE Init_u
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine implements the fourth-order Runge-Kutta Method (RK4) for numerically integrating ordinary differential equations:
!!
!!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!!   Define constants k1, k2, k3, and k4 as 
!!        k1 = dt * f(t        , x_t        )
!!        k2 = dt * f(t + dt/2 , x_t + k1/2 )
!!        k3 = dt * f(t + dt/2 , x_t + k2/2 ), and
!!        k4 = dt * f(t + dt   , x_t + k3   ).
!!   Then the continuous states at t = t + dt are
!!        x_(t+dt) = x_t + k1/6 + k2/3 + k3/3 + k4/6 + O(dt^5)
!!
!! For details, see:
!! Press, W. H.; Flannery, B. P.; Teukolsky, S. A.; and Vetterling, W. T. "Runge-Kutta Method" and "Adaptive Step Size Control for 
!!   Runge-Kutta." Sections 16.1 and 16.2 in Numerical Recipes in FORTRAN: The Art of Scientific Computing, 2nd ed. Cambridge, England: 
!!   Cambridge University Press, pp. 704-716, 1992.
SUBROUTINE ED_RK4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                   INTENT(IN   )  :: t           !< Current simulation time in seconds
      INTEGER(IntKi),               INTENT(IN   )  :: n           !< time step number
      TYPE(ED_InputType),           INTENT(INOUT)  :: u(:)        !< Inputs at t (out only for mesh record-keeping in ExtrapInterp routine)
      REAL(DbKi),                   INTENT(IN   )  :: utimes(:)   !< times of input
      TYPE(ED_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(ED_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
      TYPE(ED_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(ED_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t (possibly a guess)
      TYPE(ED_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states
      TYPE(ED_MiscVarType),         INTENT(INOUT)  :: m           !< misc/optimization variables
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! local variables
         
      TYPE(ED_ContinuousStateType)                 :: xdot        ! time derivatives of continuous states      
      TYPE(ED_ContinuousStateType)                 :: k1          ! RK4 constant; see above
      TYPE(ED_ContinuousStateType)                 :: k2          ! RK4 constant; see above 
      TYPE(ED_ContinuousStateType)                 :: k3          ! RK4 constant; see above 
      TYPE(ED_ContinuousStateType)                 :: k4          ! RK4 constant; see above 
      TYPE(ED_ContinuousStateType)                 :: x_tmp       ! Holds temporary modification to x
      TYPE(ED_InputType)                           :: u_interp    ! interpolated value of inputs 

      INTEGER(IntKi)                               :: ErrStat2    ! local error status
      CHARACTER(ErrMsgLen)                         :: ErrMsg2     ! local error message (ErrMsg)
      
      
      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      CALL ED_CopyContState( x, k1, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
      CALL ED_CopyContState( x, k2, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
      CALL ED_CopyContState( x, k3, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
      CALL ED_CopyContState( x, k4,    MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
      CALL ED_CopyContState( x, x_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN


      CALL ED_CopyInput( u(1), u_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN
                     
      ! interpolate u to find u_interp = u(t)
      CALL ED_Input_ExtrapInterp( u, utimes, u_interp, t, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN
!      HSSBrTrq_at_t = u_interp%HSSBrTrqC
!      OtherState%HSSBrTrqC = SIGN( u_interp%HSSBrTrqC, x%QDT(DOF_GeAz) )         
!      OtherState%HSSBrTrq  = OtherState%HSSBrTrqC         

      ! find xdot at t
      CALL ED_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, m, xdot, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      k1%qt  = p%dt * xdot%qt
      k1%qdt = p%dt * xdot%qdt
  
      x_tmp%qt  = x%qt  + 0.5 * k1%qt
      x_tmp%qdt = x%qdt + 0.5 * k1%qdt

      ! interpolate u to find u_interp = u(t + dt/2)
      CALL ED_Input_ExtrapInterp(u, utimes, u_interp, t+0.5*p%dt, ErrStat2, ErrMsg2)
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN
!      u_interp%HSSBrTrqC = max(0.0_ReKi, min(u_interp%HSSBrTrqC, HSSBrTrq_at_t )) ! hack for extrapolation of limits       
!      OtherState%HSSBrTrqC = SIGN( u_interp%HSSBrTrqC, x_tmp%QDT(DOF_GeAz) )         
!      OtherState%HSSBrTrq  = OtherState%HSSBrTrqC         

      ! find xdot at t + dt/2
      CALL ED_CalcContStateDeriv( t + 0.5*p%dt, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      k2%qt  = p%dt * xdot%qt
      k2%qdt = p%dt * xdot%qdt

      x_tmp%qt  = x%qt  + 0.5 * k2%qt
      x_tmp%qdt = x%qdt + 0.5 * k2%qdt

      ! find xdot at t + dt/2
!      u_interp%HSSBrTrqC = max(0.0_ReKi, min(u_interp%HSSBrTrqC, HSSBrTrq_at_t )) ! hack for extrapolation of limits       
!      OtherState%HSSBrTrqC = SIGN( u_interp%HSSBrTrqC, x_tmp%QDT(DOF_GeAz) )         
!      OtherState%HSSBrTrq  = OtherState%HSSBrTrqC         
      CALL ED_CalcContStateDeriv( t + 0.5*p%dt, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      k3%qt  = p%dt * xdot%qt
      k3%qdt = p%dt * xdot%qdt

      x_tmp%qt  = x%qt  + k3%qt
      x_tmp%qdt = x%qdt + k3%qdt

      ! interpolate u to find u_interp = u(t + dt)
      CALL ED_Input_ExtrapInterp(u, utimes, u_interp, t + p%dt, ErrStat2, ErrMsg2)
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN
!      u_interp%HSSBrTrqC = max(0.0_ReKi, min(u_interp%HSSBrTrqC, HSSBrTrq_at_t )) ! hack for extrapolation of limits       
!      OtherState%HSSBrTrqC = SIGN( u_interp%HSSBrTrqC, x_tmp%QDT(DOF_GeAz) )         
!      OtherState%HSSBrTrq  = OtherState%HSSBrTrqC         

      ! find xdot at t + dt
      CALL ED_CalcContStateDeriv( t + p%dt, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      k4%qt  = p%dt * xdot%qt
      k4%qdt = p%dt * xdot%qdt

      x%qt  = x%qt  +  ( k1%qt  + 2. * k2%qt  + 2. * k3%qt  + k4%qt  ) / 6.      
      x%qdt = x%qdt +  ( k1%qdt + 2. * k2%qdt + 2. * k3%qdt + k4%qdt ) / 6.      

         ! clean up local variables:
      CALL ExitThisRoutine(  )
         
CONTAINS      
   !...............................................................................................................................
   SUBROUTINE ExitThisRoutine()
   ! This subroutine destroys all the local variables
   !...............................................................................................................................

         ! local variables
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)
   
   
      CALL ED_DestroyContState( xdot,     ErrStat3, ErrMsg3 )
      CALL ED_DestroyContState( k1,       ErrStat3, ErrMsg3 )
      CALL ED_DestroyContState( k2,       ErrStat3, ErrMsg3 )
      CALL ED_DestroyContState( k3,       ErrStat3, ErrMsg3 )
      CALL ED_DestroyContState( k4,       ErrStat3, ErrMsg3 )
      CALL ED_DestroyContState( x_tmp,    ErrStat3, ErrMsg3 )

      CALL ED_DestroyInput(     u_interp, ErrStat3, ErrMsg3 )
         
   END SUBROUTINE ExitThisRoutine      
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)

         ! local variables
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)

      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'ED_RK4:'//TRIM(Msg)         
         ErrStat = MAX(ErrStat,ErrID)
         
         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         
         IF ( ErrStat >= AbortErrLev ) CALL ExitThisRoutine( )                  
                  
         
      END IF

   END SUBROUTINE CheckError                    
      
END SUBROUTINE ED_RK4
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine implements the fourth-order Adams-Bashforth Method (RK4) for numerically integrating ordinary differential 
!! equations:
!!
!!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!!
!!   x(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
!!
!!  See, e.g.,
!!  http://en.wikipedia.org/wiki/Linear_multistep_method
!!
!!  or
!!
!!  K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
SUBROUTINE ED_AB4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           !< time step number
      TYPE(ED_InputType),             INTENT(INOUT)  :: u(:)        !< Inputs at t (out only for mesh record-keeping in ExtrapInterp routine)
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   !< times of input
      TYPE(ED_ParameterType),         INTENT(IN   )  :: p           !< Parameters
      TYPE(ED_ContinuousStateType),   INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
      TYPE(ED_DiscreteStateType),     INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(ED_ConstraintStateType),   INTENT(IN   )  :: z           !< Constraint states at t (possibly a guess)
      TYPE(ED_OtherStateType),        INTENT(INOUT)  :: OtherState  !< Other states
      TYPE(ED_MiscVarType),           INTENT(INOUT)  :: m           !< misc/optimization variables
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


      ! local variables
      TYPE(ED_InputType)                             :: u_interp
      TYPE(ED_ContinuousStateType)                   :: xdot
         
      INTEGER(IntKi)                                 :: ErrStat2    ! local error status
      CHARACTER(ErrMsgLen)                           :: ErrMsg2     ! local error message (ErrMsg)


      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 
      
      
      if (OtherState%n .lt. n) then

         OtherState%n = n
            
         ! Update IC() index so IC(1) is the location of xdot values at n.
         ! (this allows us to shift the indices into the array, not copy all of the values)
         OtherState%IC = CSHIFT( OtherState%IC, -1 ) ! circular shift of all values to the right
            
      elseif (OtherState%n .gt. n) then
 
         CALL CheckError( ErrID_Fatal, ' Backing up in time is not supported with a multistep method.')
         RETURN

      endif        
      
      
      ! Allocate the input arrays
      CALL ED_CopyInput( u(1), u_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      
      ! need xdot at t
      CALL ED_Input_ExtrapInterp(u, utimes, u_interp, t, ErrStat2, ErrMsg2)
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN                  
      IF (EqualRealNos( x%qdt(DOF_GeAz) ,0.0_R8Ki ) ) THEN
         OtherState%HSSBrTrqC = u_interp%HSSBrTrqC
      ELSE
         OtherState%HSSBrTrqC  = SIGN( u_interp%HSSBrTrqC, real(x%qdt(DOF_GeAz),ReKi) ) ! hack for HSS brake (need correct sign)
      END IF
      OtherState%HSSBrTrq   = OtherState%HSSBrTrqC
      OtherState%SgnPrvLSTQ = OtherState%SgnLSTQ(OtherState%IC(2))
      
      CALL ED_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, m, xdot, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         
         CALL ED_CopyContState(xdot, OtherState%xdot ( OtherState%IC(1) ), MESH_NEWCOPY, ErrStat2, ErrMsg2)
            CALL CheckError(ErrStat2,ErrMsg2)
            IF ( ErrStat >= AbortErrLev ) RETURN

                                                    
      if (n .le. 2) then
                                               
         CALL ED_RK4(t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat2, ErrMsg2 )
            CALL CheckError(ErrStat2,ErrMsg2)
            IF ( ErrStat >= AbortErrLev ) RETURN

      else
         
         x%qt  = x%qt  + p%DT24 * ( 55.*OtherState%xdot(OtherState%IC(1))%qt  - 59.*OtherState%xdot(OtherState%IC(2))%qt   &
                                  + 37.*OtherState%xdot(OtherState%IC(3))%qt   - 9.*OtherState%xdot(OtherState%IC(4))%qt )

         x%qdt = x%qdt + p%DT24 * ( 55.*OtherState%xdot(OtherState%IC(1))%qdt - 59.*OtherState%xdot(OtherState%IC(2))%qdt  &
                                  + 37.*OtherState%xdot(OtherState%IC(3))%qdt  - 9.*OtherState%xdot(OtherState%IC(4))%qdt )
         
         
            ! Make sure the HSS brake will not reverse the direction of the HSS
            !   for the next time step.  Do this by computing the predicted value
            !   of x%qt(); QD(DOF_GeAz,IC(NMX)) as will be done during the next time step.
            ! Only do this after the first few time steps since it doesn't work
            !   for the Runga-Kutta integration scheme.
   
         
         CALL FixHSSBrTq ( 'P', p, x, OtherState, m, ErrStat2, ErrMsg2 )
            CALL CheckError(ErrStat2,ErrMsg2)
            IF ( ErrStat >= AbortErrLev ) RETURN
            
      endif
            
      OtherState%SgnPrvLSTQ = SignLSSTrq(p, m)   
      OtherState%SgnLSTQ(OtherState%IC(1)) = OtherState%SgnPrvLSTQ 
      
      
         ! clean up local variables:
      CALL ExitThisRoutine()
      
CONTAINS      
   !...............................................................................................................................
   SUBROUTINE ExitThisRoutine()
   ! This subroutine destroys all the local variables
   !...............................................................................................................................

         ! local variables
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)
   
   
      CALL ED_DestroyInput(     u_interp, ErrStat3, ErrMsg3 )
      CALL ED_DestroyContState( xdot,     ErrStat2, ErrMsg3 )
      
   END SUBROUTINE ExitThisRoutine    
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)

         ! local variables
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)

      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'ED_AB4:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         
         IF ( ErrStat >= AbortErrLev ) CALL ExitThisRoutine( )                  
         
      END IF

   END SUBROUTINE CheckError            
         
END SUBROUTINE ED_AB4
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine implements the fourth-order Adams-Bashforth-Moulton Method (RK4) for numerically integrating ordinary 
!! differential equations:
!!
!!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!!
!!   Adams-Bashforth Predictor: \n
!!   x^p(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
!!
!!   Adams-Moulton Corrector: \n
!!   x(t+dt) = x(t)  + (dt / 24.) * ( 9.*f(t+dt,x^p) + 19.*f(t,x) - 5.*f(t-dt,x) + 1.*f(t-2.*dt,x) )
!!
!!  See, e.g.,
!!  http://en.wikipedia.org/wiki/Linear_multistep_method
!!
!!  or
!!
!!  K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
SUBROUTINE ED_ABM4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           !< time step number
      TYPE(ED_InputType),             INTENT(INOUT)  :: u(:)        !< Inputs at t (out only for mesh record-keeping in ExtrapInterp routine)
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   !< times of input
      TYPE(ED_ParameterType),         INTENT(IN   )  :: p           !< Parameters
      TYPE(ED_ContinuousStateType),   INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
      TYPE(ED_DiscreteStateType),     INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(ED_ConstraintStateType),   INTENT(IN   )  :: z           !< Constraint states at t (possibly a guess)
      TYPE(ED_OtherStateType),        INTENT(INOUT)  :: OtherState  !< Other states
      TYPE(ED_MiscVarType),           INTENT(INOUT)  :: m           !< misc/optimization variables
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! local variables

      TYPE(ED_InputType)                             :: u_interp    ! Inputs at t
      TYPE(ED_ContinuousStateType)                   :: x_pred      ! Continuous states at t
      TYPE(ED_ContinuousStateType)                   :: xdot_pred   ! Derivative of continuous states at t

      INTEGER(IntKi)                                 :: ErrStat2    ! local error status
      CHARACTER(ErrMsgLen)                           :: ErrMsg2     ! local error message (ErrMsg)
      
      
      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 
      
         ! predict:

      CALL ED_CopyContState(x, x_pred, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      CALL ED_AB4( t, n, u, utimes, p, x_pred, xd, z, OtherState, m, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

         
      if (n .gt. 2_IntKi) then
         
            ! correct:
         
            ! allocate the arrays in u_interp
         CALL ED_CopyInput( u(1), u_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
            CALL CheckError(ErrStat2,ErrMsg2)
            IF ( ErrStat >= AbortErrLev ) RETURN
         
         CALL ED_Input_ExtrapInterp(u, utimes, u_interp, t + p%dt, ErrStat2, ErrMsg2)
            CALL CheckError(ErrStat2,ErrMsg2)
            IF ( ErrStat >= AbortErrLev ) RETURN
            
         u_interp%HSSBrTrqC = max(0.0_ReKi, min(u_interp%HSSBrTrqC, ABS( OtherState%HSSBrTrqC) )) ! hack for extrapolation of limits  (OtherState%HSSBrTrqC is HSSBrTrqC at t)     
         IF (EqualRealNos( x_pred%qdt(DOF_GeAz) ,0.0_R8Ki ) ) THEN
            OtherState%HSSBrTrqC = u_interp%HSSBrTrqC
         ELSE
            OtherState%HSSBrTrqC  = SIGN( u_interp%HSSBrTrqC, real(x_pred%qdt(DOF_GeAz),ReKi) ) ! hack for HSS brake (need correct sign)
         END IF
         OtherState%HSSBrTrq  = OtherState%HSSBrTrqC

         CALL ED_CalcContStateDeriv(t + p%dt, u_interp, p, x_pred, xd, z, OtherState, m, xdot_pred, ErrStat2, ErrMsg2 )
            CALL CheckError(ErrStat2,ErrMsg2)
            IF ( ErrStat >= AbortErrLev ) RETURN

         
         x%qt  = x%qt  + p%DT24 * ( 9. * xdot_pred%qt +  19. * OtherState%xdot(OtherState%IC(1))%qt &
                                                        - 5. * OtherState%xdot(OtherState%IC(2))%qt &
                                                        + 1. * OtherState%xdot(OtherState%IC(3))%qt )

         x%qdt = x%qdt + p%DT24 * ( 9. * xdot_pred%qdt + 19. * OtherState%xdot(OtherState%IC(1))%qdt &
                                                       -  5. * OtherState%xdot(OtherState%IC(2))%qdt &
                                                       +  1. * OtherState%xdot(OtherState%IC(3))%qdt )
         
                  
          ! Make sure the HSS brake has not reversed the direction of the HSS:
         
         CALL FixHSSBrTq ( 'C', p, x, OtherState, m, ErrStat2, ErrMsg2 )      
            CALL CheckError(ErrStat2,ErrMsg2)
            IF ( ErrStat >= AbortErrLev ) RETURN
         OtherState%SgnPrvLSTQ = SignLSSTrq(p, m)
         OtherState%SgnLSTQ(OtherState%IC(1)) = OtherState%SgnPrvLSTQ 
                                                                    
      else

         x%qt  = x_pred%qt
         x%qdt = x_pred%qdt

      endif
      
      
         ! clean up local variables:
      CALL ExitThisRoutine()
      
CONTAINS      
   !...............................................................................................................................
   SUBROUTINE ExitThisRoutine()
   ! This subroutine destroys all the local variables
   !...............................................................................................................................

         ! local variables
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)
   
   
      CALL ED_DestroyContState( xdot_pred,  ErrStat3, ErrMsg3 )
      CALL ED_DestroyContState( x_pred,     ErrStat3, ErrMsg3 )
      CALL ED_DestroyInput(     u_interp,   ErrStat3, ErrMsg3 )               
      
   END SUBROUTINE ExitThisRoutine    
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)

         ! local variables
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)

      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'ED_ABM4:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) CALL ExitThisRoutine( )                  
         
      END IF

   END SUBROUTINE CheckError                 

END SUBROUTINE ED_ABM4
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine generates the summary file, which contains a regurgitation of  the input data and interpolated flexible body data.
SUBROUTINE ED_PrintSum( p, OtherState, ErrStat, ErrMsg )

      ! passed variables
   TYPE(ED_ParameterType),    INTENT(IN   )  :: p                       !< Parameters of the structural dynamics module
   TYPE(ED_OtherStateType),   INTENT(IN   )  :: OtherState              !< Other states of the structural dynamics module 
   INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat                 !< Error status of the operation
   CHARACTER(*),              INTENT(  OUT)  :: ErrMsg                  !< Error message if ErrStat /= ErrID_None


      ! Local variables.

   INTEGER(IntKi)               :: I                                               ! Index for the nodes.
   INTEGER(IntKi)               :: K                                               ! Generic index (also for the blade number).
   INTEGER(IntKi)               :: UnSu                                            ! I/O unit number for the summary output file

   CHARACTER(*), PARAMETER      :: Fmt1      = "(34X,3(6X,'Blade',I2,:))"          ! Format for outputting blade headings.
   CHARACTER(*), PARAMETER      :: Fmt2      = "(34X,3(6X,A,:))"                   ! Format for outputting blade headings.
   CHARACTER(*), PARAMETER      :: FmtDat    = '(A,T35,3(:,F13.3))'                ! Format for outputting mass and modal data.
   CHARACTER(*), PARAMETER      :: FmtDatT   = '(A,T35,1(:,F13.8))'                ! Format for outputting time steps.
   CHARACTER(100)               :: RotorType                                       ! Text description of rotor.

   CHARACTER(30)                :: OutPFmtS                                        ! Format to print list of selected output channel names to summary file
   CHARACTER(30)                :: OutPFmt                                         ! Format to print list of selected output channels to summary file
   CHARACTER(10)                :: DOFEnabled                                      ! String to say if a DOF is enabled or disabled
   CHARACTER(ChanLen),PARAMETER :: TitleStr(2) = (/ 'Parameter', 'Units    ' /)
   CHARACTER(ChanLen),PARAMETER :: TitleStrLines(2) = (/ '---------------', '---------------' /)

   ! Open the summary file and give it a heading.
   
   CALL GetNewUnit( UnSu, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN

   CALL OpenFOutFile ( UnSu, TRIM( p%RootName )//'.sum', ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN

   
      ! Heading:
   WRITE (UnSu,'(/,A)')  'This summary information was generated by '//TRIM( GetNVD(ED_Ver) )// &
                         ' on '//CurDate()//' at '//CurTime()//'.'
   !WRITE (UnSu,'(//,1X,A)')  TRIM( p%FTitle )


   !..................................
   ! Turbine features.
   !..................................

   WRITE (UnSu,'(//,A,/)')  'Turbine features:'

   IF ( p%OverHang > 0.0 )  THEN
      RotorType = 'Downwind,'
   ELSE
      RotorType = 'Upwind,'
   ENDIF
   RotorType = TRIM(RotorType)//' '//trim(Num2LStr(p%NumBl))//'-bladed rotor'
   IF ( p%NumBl==2 )  THEN
      IF ( p%DOF_Flag(DOF_Teet) ) THEN ! NOTE: two "ifs" required since DOF_Teet might be out of bound
         RotorType = TRIM(RotorType)//' with teetering hub.'
      ELSE
         RotorType = TRIM(RotorType)//' with rigid hub.'
      ENDIF
   ELSE
      RotorType = TRIM(RotorType)//' with rigid hub.'
   ENDIF

   WRITE    (UnSu,'(A)')  '            '//TRIM(RotorType)

   WRITE    (UnSu,'(A)')  '            The model has '//TRIM(Num2LStr( p%DOFs%NActvDOF ))//' of '// &
                                        TRIM(Num2LStr( p%NDOF ))//' DOFs active (enabled) at start-up.'

   DO I=1,p%NDOF
   
      IF ( p%DOF_Flag( I ) ) THEN
         DOFEnabled = 'Enabled'
      ELSE
         DOFEnabled = 'Disabled'
      END IF
      
      K = INDEX( p%DOF_Desc(I), ' (internal DOF index' )
      IF (K == 0) K = LEN_TRIM(p%DOF_Desc(I))
      
      WRITE ( UnSu, '(1x,A10,1x,A)' ) DOFEnabled, p%DOF_Desc(I)(:K)
         
      
   END DO
   

      ! Time steps.

   WRITE (UnSu,'(//,A,/)')  'Time steps:'

   WRITE (UnSu,FmtDatT) '    Structural            (s)     ', p%DT

      ! Some calculated parameters.

   WRITE (UnSu,'(//,A,/)')  'Some calculated parameters:'

   WRITE (UnSu,FmtDat ) '    Hub-Height            (m)     ', p%HubHt
   WRITE (UnSu,FmtDat ) '    Flexible Tower Length (m)     ', p%TwrFlexL
IF (.NOT. p%BD4Blades) THEN
   WRITE (UnSu,FmtDat ) '    Flexible Blade Length (m)     ', p%BldFlexL
ELSE
   WRITE (UnSu,'(A)' )  '    Flexible Blade Length (m)               N/A'
END IF

      ! Rotor properties:

   WRITE (UnSu,'(//,A,/)')  'Rotor mass properties:'

   WRITE (UnSu,FmtDat ) '    Rotor Mass            (kg)    ', p%RotMass
   WRITE (UnSu,FmTDat ) '    Rotor Inertia         (kg-m^2)', p%RotINer

   WRITE (UnSu,Fmt1   ) ( K,         K=1,p%NumBl )
   WRITE (UnSu,Fmt2   ) ( '-------', K=1,p%NumBl )

IF (.NOT. p%BD4Blades) THEN

   WRITE (UnSu,FmtDat ) '    Mass                  (kg)    ', ( p%BldMass  (K), K=1,p%NumBl )
   WRITE (UnSu,FmtDat ) '    Second Mass Moment    (kg-m^2)', ( p%SecondMom(K), K=1,p%NumBl )
   WRITE (UnSu,FmtDat ) '    First Mass Moment     (kg-m)  ', ( p%FirstMom (K), K=1,p%NumBl )
   WRITE (UnSu,FmtDat ) '    Center of Mass        (m)     ', ( p%BldCG    (K), K=1,p%NumBl )

ELSE

   WRITE (UnSu,'(A)' ) '    Mass                  (kg)              N/A'
   WRITE (UnSu,'(A)' ) '    Second Mass Moment    (kg-m^2)          N/A'
   WRITE (UnSu,'(A)' ) '    First Mass Moment     (kg-m)            N/A'
   WRITE (UnSu,'(A)' ) '    Center of Mass        (m)               N/A'
   
END IF 

      ! Output additional masses:

   WRITE (UnSu,'(//,A,/)')  'Additional mass properties:'

   WRITE (UnSu,FmtDat ) '    Tower-top Mass        (kg)    ', p%TwrTpMass
   WRITE (UnSu,FmtDat ) '    Tower Mass            (kg)    ', p%TwrMass
   !WRITE (UnSu,FmtDat ) '    Turbine Mass          (kg)    ', p%TurbMass
   WRITE (UnSu,FmtDat ) '    Platform Mass         (kg)    ', p%PtfmMass
   WRITE (UnSu,FmtDat ) '    Mass Incl. Platform   (kg)    ', p%TurbMass + p%PtfmMass !TotalMass !bjj TotalMass not used anywhere else so removed it


      ! Interpolated tower properties.

   WRITE (UnSu,"(//,'Interpolated tower properties:',/)")

   WRITE (UnSu,'(A)')  'Node  TwFract   HNodes  DHNodes  TMassDen    FAStiff    SSStiff'
   WRITE (UnSu,'(A)')  ' (-)      (-)      (m)      (m)    (kg/m)     (Nm^2)     (Nm^2)'

   DO I=1,p%TwrNodes
      WRITE(UnSu,'(I4,3F9.3,F10.3,2ES11.3)')  I, p%HNodesNorm(I), p%HNodes(I), p%DHNodes(I), p%MassT(I), &
                                                   p%StiffTFA(I),  p%StiffTSS(I)
   ENDDO ! I


      ! Interpolated blade properties.
IF (.NOT. p%BD4Blades) THEN

   DO K=1,p%NumBl

      WRITE (UnSu,'(//,A,I1,A,/)')  'Interpolated blade ', K, ' properties:'

      WRITE (UnSu,'(A)')  'Node  BlFract   RNodes  DRNodes PitchAxis  StrcTwst  BMassDen    FlpStff    EdgStff'
      WRITE (UnSu,'(A)')  ' (-)      (-)      (m)      (m)       (-)     (deg)    (kg/m)     (Nm^2)     (Nm^2)'

      DO I=1,p%BldNodes
         WRITE(UnSu,'(I4,3F9.3,3F10.3,2ES11.3)')  I, p%RNodesNorm(I), p%RNodes(I) + p%HubRad, p%DRNodes(I), &
                                                      p%PitchAxis(K,I),p%ThetaS(K,I)*R2D, p%MassB(K,I), &
                                                      p%StiffBF(K,I), p%StiffBE(K,I)
      ENDDO ! I

   ENDDO ! K
   
END IF

   OutPFmt  = '( I4, 3X,A '//TRIM(Num2LStr(ChanLen))//',1 X, A'//TRIM(Num2LStr(ChanLen))//' )'
   OutPFmtS = '( A4, 3X,A '//TRIM(Num2LStr(ChanLen))//',1 X, A'//TRIM(Num2LStr(ChanLen))//' )'
   WRITE (UnSu,'(//,A,//)')  'Requested Outputs:'
  !WRITE (UnSu,"(/, '  Col  Parameter       Units', /, '  ---  --------------  ----------')")
   WRITE (UnSu,OutPFmtS)  "Col", TitleStr
   WRITE (UnSu,OutPFmtS)  "---", TitleStrLines
   DO I = 0,p%NumOuts
      WRITE (UnSu,OutPFmt)  I, p%OutParam(I)%Name, p%OutParam(I)%Units
   END DO             

   IF (.not. p%BD4Blades) THEN
      WRITE (UnSu,'(2x,A)')
      WRITE (UnSu,'(2x,A)')
      WRITE (UnSu,'(2x,A)')  'Requested Output Channels at each blade station:'
      WRITE (UnSu,OutPFmtS)  "Col", TitleStr
      WRITE (UnSu,OutPFmtS)  "---", TitleStrLines
      !WRITE (UnSu,'(2x,A)')  'Col   Parameter       Units'
      !WRITE (UnSu,'(2x,A)')  '----  --------------  ----------'
      DO I = 1,p%BldNd_NumOuts
         WRITE (UnSu,OutPFmt)  I, p%BldNd_OutParam(I)%Name, p%BldNd_OutParam(I)%Units
      END DO
   ENDIF

   CLOSE(UnSu)

RETURN
END SUBROUTINE ED_PrintSum
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is used to adjust the HSSBrTrq value if the absolute
!!   magnitudue of the HSS brake torque was strong enough to reverse
!!   the direction of the HSS, which is a physically impossible
!!   situation.  The problem arises since we are integrating in
!!   discrete time, not continuous time.
SUBROUTINE FixHSSBrTq ( Integrator, p, x, OtherState, m, ErrStat, ErrMsg )


   ! Passed variables:

   TYPE(ED_ParameterType),      INTENT(IN   )  :: p                       !< Parameters of the structural dynamics module
   TYPE(ED_OtherStateType),     INTENT(INOUT)  :: OtherState              !< Other states of the structural dynamics module 
   TYPE(ED_MiscVarType),        INTENT(INOUT)  :: m                       !< misc (optimization) variables
   TYPE(ED_ContinuousStateType),INTENT(INOUT)  :: x                       !< Continuous states of the structural dynamics module at n+1
   CHARACTER(1),                INTENT(IN   )  :: Integrator              !< A string holding the current integrator being used.
   INTEGER(IntKi),              INTENT(  OUT)  :: ErrStat                 !< Error status of the operation
   CHARACTER(*),                INTENT(  OUT)  :: ErrMsg                  !< Error message if ErrStat /= ErrID_None


   ! Local variables:

   REAL(ReKi)                             :: RqdFrcGeAz                           ! The force term required to produce RqdQD2GeAz.
   REAL(ReKi)                             :: RqdQD2GeAz                           ! The required QD2T(DOF_GeAz) to cause the HSS to stop rotating.

   INTEGER                                :: I                                    ! Loops through all DOFs.
   INTEGER(IntKi)                         :: ErrStat2
   CHARACTER(ErrMsgLen)                   :: ErrMsg2
   CHARACTER(*), PARAMETER                :: RoutineName = 'FixHSSBrTq'

   

   ErrStat = ErrID_None
   ErrMsg  = ""

   IF ( .NOT. p%DOF_Flag(DOF_GeAz) .OR. EqualRealNos(OtherState%HSSBrTrqC, 0.0_ReKi ) )  RETURN


      ! The absolute magnitude of the HSS brake must have been too great
      !   that the HSS direction was reversed.  What should have happened
      !   is that the HSS should have stopped rotating.  In other words,
      !   QD(DOF_GeAz,IC(NMX)) should equal zero!  Determining what
      !   QD2T(DOF_GeAz) will make QD(DOF_GeAz,IC(NMX)) = 0, depends on
      !   which integrator we are using.

   
   SELECT CASE (Integrator)

   CASE ('C')   ! Corrector

      ! Find the required QD2T(DOF_GeAz) to cause the HSS to stop rotating (RqdQD2GeAz).
      ! This is found by solving the corrector formula for QD2(DOF_GeAz,IC(NMX))
      !   when QD(DOF_GeAz,IC(NMX)) equals zero.

      RqdQD2GeAz = ( -      OtherState%xdot(OtherState%IC(1))%qt (DOF_GeAz)/ p%DT24 &
                     - 19.0*OtherState%xdot(OtherState%IC(1))%qdt(DOF_GeAz)         &
                     +  5.0*OtherState%xdot(OtherState%IC(2))%qdt(DOF_GeAz)         &
                     -      OtherState%xdot(OtherState%IC(3))%qdt(DOF_GeAz)         ) / 9.0
      
   CASE ('P')   ! Predictor

      ! Find the required QD2T(DOF_GeAz) to cause the HSS to stop rotating (RqdQD2GeAz).
      ! This is found by solving the predictor formula for QD2(DOF_GeAz,IC(1))
      !   when QD(DOF_GeAz,IC(NMX)) equals zero.

      RqdQD2GeAz = ( -      OtherState%xdot(OtherState%IC(1))%qt( DOF_GeAz)  / p%DT24 &
                     + 59.0*OtherState%xdot(OtherState%IC(2))%qdt(DOF_GeAz) &
                     - 37.0*OtherState%xdot(OtherState%IC(3))%qdt(DOF_GeAz) &
                     +  9.0*OtherState%xdot(OtherState%IC(4))%qdt(DOF_GeAz)   )/55.0
            
   END SELECT


   ! Rearrange the augmented matrix of equations of motion to account
   !   for the known acceleration of the generator azimuth DOF.  To
   !   do this, make the known inertia like an applied force to the
   !   system.  Then set force QD2T(DOF_GeAz) to equal the known
   !   acceleration in the augmented matrix of equations of motion:
   ! Here is how the new equations are derived.  First partition the
   !   augmented matrix as follows, where Qa are the unknown
   !   accelerations, Qb are the known accelerations, Fa are the
   !   known forces, and Fb are the unknown forces:
   !      [Caa Cab]{Qa}={Fa}
   !      [Cba Cbb]{Qb}={Fb}
   !   By rearranging, the equations for the unknown and known
   !   accelerations are as follows:
   !      [Caa]{Qa}={Fa}-[Cab]{Qb} and [I]{Qb}={Qb}
   !   Combining these two sets of equations into one set yields:
   !      [Caa 0]{Qa}={{Fa}-[Cab]{Qb}}
   !      [  0 I]{Qb}={          {Qb}}
   !   Once this equation is solved, the unknown force can be found from:
   !      {Fb}=[Cba]{Qa}+[Cbb]{Qb}

   m%OgnlGeAzRo    = m%AugMat(DOF_GeAz,:)  ! used for HSS Brake hack; copy this row before modifying the old matrix
   
  
   DO I = 1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs

      m%AugMat(p%DOFs%SrtPS(I),    p%NAUG) = m%AugMat(p%DOFs%SrtPS(I),p%NAUG) &
                                                    - m%AugMat(p%DOFs%SrtPS(I),DOF_GeAz)*RqdQD2GeAz  ! {{Fa}-[Cab]{Qb}}
      m%AugMat(p%DOFs%SrtPS(I),DOF_GeAz)   = 0.0                                                     ! [0]
      m%AugMat(DOF_GeAz, p%DOFs%SrtPS(I))  = 0.0                                                     ! [0]

   ENDDO             ! I - All active (enabled) DOFs

   m%AugMat(DOF_GeAz,DOF_GeAz) = 1.0                                                           ! [I]{Qb}={Qb}
   m%AugMat(DOF_GeAz,  p%NAUG) = RqdQD2GeAz                                                    !


   ! Invert the matrix to solve for the new (updated) accelerations.  Like in
   !   CalcContStateDeriv(), the accelerations are returned by Gauss() in the first NActvDOF
   !   elements of the solution vector, SolnVec().  These are transfered to the
   !   proper index locations of the acceleration vector QD2T() using the
   !   vector subscript array SrtPS(), after Gauss() has been called:

   ! Invert the matrix to solve for the accelerations. The accelerations are returned by Gauss() in the first NActvDOF elements
   !   of the solution vector, SolnVec(). These are transfered to the proper index locations of the acceleration vector QD2T()
   !   using the vector subscript array SrtPS(), after Gauss() has been called:

      m%AugMat_factor = m%AugMat( p%DOFs%SrtPS( 1:p%DOFs%NActvDOF ), p%DOFs%SrtPSNAUG(1:p%DOFs%NActvDOF) )
      m%SolnVec       = m%AugMat( p%DOFs%SrtPS( 1:p%DOFs%NActvDOF ), p%DOFs%SrtPSNAUG(1+p%DOFs%NActvDOF) )
   
      CALL LAPACK_getrf( M=p%DOFs%NActvDOF, N=p%DOFs%NActvDOF, A=m%AugMat_factor, IPIV=m%AugMat_pivot, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)   
         IF ( ErrStat >= AbortErrLev ) RETURN
      
      CALL LAPACK_getrs( TRANS='N',N=p%DOFs%NActvDOF, A=m%AugMat_factor,IPIV=m%AugMat_pivot, B=m%SolnVec, ErrStat=ErrStat2, ErrMsg=ErrMsg2)
   
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) RETURN
   
   
      ! Find the force required to produce RqdQD2GeAz from the equations of
      !   motion using the new accelerations:

   RqdFrcGeAz = 0.0
   DO I = 1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs
      ! bjj: use m%SolnVec(I) instead of m%QD2T(p%DOFs%SrtPS(I)) here; then update m%QD2T(p%DOFs%SrtPS(I))
      !      later if necessary
      !RqdFrcGeAz = RqdFrcGeAz + OgnlGeAzRo(SrtPS(I))*m%QD2T(p%DOFs%SrtPS(I))  ! {Fb}=[Cba]{Qa}+[Cbb]{Qb}
      RqdFrcGeAz = RqdFrcGeAz + m%OgnlGeAzRo(p%DOFs%SrtPS(I))*m%SolnVec(I)  ! {Fb}=[Cba]{Qa}+[Cbb]{Qb}
   ENDDO             ! I - All active (enabled) DOFs


      ! Find the HSSBrTrq necessary to bring about this force:

   OtherState%HSSBrTrq = OtherState%HSSBrTrqC & 
                       + ( ( m%OgnlGeAzRo(p%NAUG) - RqdFrcGeAz )*m%RtHS%GBoxEffFac/ABS(p%GBRatio) )


      ! Make sure this new HSSBrTrq isn't larger in absolute magnitude than
      !   the original HSSBrTrq.  Indeed, the new HSSBrTrq can't be larger than
      !   the old HSSBrTrq, since the old HSSBrTrq was found solely as a
      !   function of time--and is thus the maximum possible at the current
      !   time.  If the new HSSBrTrq is larger, then the reversal in direction
      !   was caused by factors other than the HSS brake--thus the original HSS
      !   brake torque values were OK to begin with.  Thus, restore the
      !   variables changed by this subroutine, back to their original values:

   IF ( ABS( OtherState%HSSBrTrq ) > ABS( OtherState%HSSBrTrqC ) )  THEN

      OtherState%HSSBrTrq = OtherState%HSSBrTrqC !OtherState%HSSBrTrqC = SIGN( u%HSSBrTrqC, x%QDT(DOF_GeAz) )
      !m%QD2T     = QD2TC

   ELSE

      ! overwrite QD2T with the new values
      m%QD2T = 0.0
      DO I = 1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs
         m%QD2T(p%DOFs%SrtPS(I)) = m%SolnVec(I)
      ENDDO             ! I - All active (enabled) DOFs
      
            
      ! Use the new accelerations to update the DOF values.  Again, this
      !   depends on the integrator type:

      SELECT CASE (Integrator)

      CASE ('C')  ! Corrector

      ! Update QD and QD2 with the new accelerations using the corrector.
      ! This will make QD(DOF_GeAz,IC(NMX)) equal to zero and adjust all
      !    of the other QDs as necessary.
      ! The Q's are unnaffected by this change.     
      
         x%qdt =                   OtherState%xdot(OtherState%IC(1))%qt &  ! qd at n
                 + p%DT24 * ( 9. * m%QD2T &                                ! the value we just changed
                           + 19. * OtherState%xdot(OtherState%IC(1))%qdt &
                           -  5. * OtherState%xdot(OtherState%IC(2))%qdt &
                           +  1. * OtherState%xdot(OtherState%IC(3))%qdt )
            

         
      CASE ('P')  ! Predictor

      ! Update QD and QD2 with the new accelerations using predictor.  
         
         x%qdt =                OtherState%xdot(OtherState%IC(1))%qt + &  ! qd at n
                 p%DT24 * ( 55.*m%QD2T &                                  ! the value we just changed
                          - 59.*OtherState%xdot(OtherState%IC(2))%qdt  &
                          + 37.*OtherState%xdot(OtherState%IC(3))%qdt  &
                           - 9.*OtherState%xdot(OtherState%IC(4))%qdt )
         
         OtherState%xdot ( OtherState%IC(1) )%qdt = m%QD2T        ! fix the history

         
      END SELECT
      
   ENDIF

   RETURN
END SUBROUTINE FixHSSBrTq
!----------------------------------------------------------------------------------------------------------------------------------

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ###### The following four routines are Jacobian routines for linearization capabilities #######
! If the module does not implement them, set ErrStat = ErrID_Fatal in ED_Init() when InitInp%Linearize is .true.
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the inputs (u). The partial derivatives dY/du, dX/du, dXd/du, and dZ/du are returned.
SUBROUTINE ED_JacobianPInput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdu, dXdu, dXddu, dZdu )
!..................................................................................................................................

   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(ED_InputType),                   INTENT(INOUT)           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(ED_ParameterType),               INTENT(IN   )           :: p          !< Parameters
   TYPE(ED_ContinuousStateType),         INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(ED_DiscreteStateType),           INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(ED_ConstraintStateType),         INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(ED_OtherStateType),              INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(ED_OutputType),                  INTENT(INOUT)           :: y          !< Output (change to inout if a mesh copy is required);
                                                                               !!   Output fields are not used by this routine, but type is   
                                                                               !!   available here so that mesh parameter information (i.e.,  
                                                                               !!   connectivity) does not have to be recalculated for dYdu.
   TYPE(ED_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dYdu(:,:)  !< Partial derivatives of output functions (Y) with respect 
                                                                               !!   to the inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXdu(:,:)  !< Partial derivatives of continuous state functions (X) with 
                                                                               !!   respect to the inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXddu(:,:) !< Partial derivatives of discrete state functions (Xd) with 
                                                                               !!   respect to the inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dZdu(:,:)  !< Partial derivatives of constraint state functions (Z) with 
                                                                               !!   respect to the inputs (u) [intent in to avoid deallocation]

   
      ! local variables
   TYPE(ED_OutputType)                               :: y_p
   TYPE(ED_OutputType)                               :: y_m
   TYPE(ED_ContinuousStateType)                      :: x_p
   TYPE(ED_ContinuousStateType)                      :: x_m
   TYPE(ED_InputType)                                :: u_perturb
   REAL(R8Ki)                                        :: delta        ! delta change in input or state
   INTEGER(IntKi)                                    :: i, j   
   
   INTEGER(IntKi)                                    :: ErrStat2
   CHARACTER(ErrMsgLen)                              :: ErrMsg2
   CHARACTER(*), PARAMETER                           :: RoutineName = 'ED_JacobianPInput'
   
   
      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''
   m%IgnoreMod = .true. ! to compute perturbations, we need to ignore the modulo function
   
      ! make a copy of the inputs to perturb
   call ED_CopyInput( u, u_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat>=AbortErrLev) then
         call cleanup()
         return
      end if
      
   

   IF ( PRESENT( dYdu ) ) THEN

      ! Calculate the partial derivative of the output functions (Y) with respect to the inputs (u) here:

      ! allocate dYdu if necessary
      if (.not. allocated(dYdu)) then
         call AllocAry(dYdu, p%Jac_ny, size(p%Jac_u_indx,1)+1, 'dYdu', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) then
            call cleanup()
            return
         end if
      end if
      
         ! make a copy of outputs because we will need two for the central difference computations (with orientations)
      call ED_CopyOutput( y, y_p, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call ED_CopyOutput( y, y_m, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) then
            call cleanup()
            return
         end if
         
      do i=1,size(p%Jac_u_indx,1)
         
            ! get u_op + delta u
         call ED_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
         call ED_Perturb_u( p, i, 1, u_perturb, delta )

            ! compute y at u_op + delta u
         call ED_CalcOutput( t, u_perturb, p, x, xd, z, OtherState, y_p, m, ErrStat2, ErrMsg2 ) 
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
         
            
            ! get u_op - delta u
         call ED_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later
         call ED_Perturb_u( p, i, -1, u_perturb, delta )
         
            ! compute y at u_op - delta u
         call ED_CalcOutput( t, u_perturb, p, x, xd, z, OtherState, y_m, m, ErrStat2, ErrMsg2 ) 
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
         
            
            ! get central difference:            
         call Compute_dY( p, y_p, y_m, delta, dYdu(:,i) )
         
      end do
      
      ! now do the extended input: sum the p%NumBl blade pitch columns
      dYdu(:,size(p%Jac_u_indx,1)+1) = dYdu(:,size(p%Jac_u_indx,1)-p%NumBl-1) ! last NumBl+2 columns are: GenTrq, YawMom, and BlPitchCom   
      do i=2,p%NumBl
         dYdu(:,size(p%Jac_u_indx,1)+1) = dYdu(:,size(p%Jac_u_indx,1)+1) + dYdu(:,size(p%Jac_u_indx,1)-p%NumBl-2+i) 
      end do
      
      
      if (ErrStat>=AbortErrLev) then
         call cleanup()
         return
      end if
      call ED_DestroyOutput( y_p, ErrStat2, ErrMsg2 ) ! we don't need this any more
      call ED_DestroyOutput( y_m, ErrStat2, ErrMsg2 ) ! we don't need this any more
      

   END IF
   

   IF ( PRESENT( dXdu ) ) THEN

      ! Calculate the partial derivative of the continuous state functions (X) with respect to the inputs (u) here:

      ! allocate dXdu if necessary
      if (.not. allocated(dXdu)) then
         call AllocAry(dXdu, p%DOFs%NActvDOF * 2, size(p%Jac_u_indx,1)+1, 'dXdu', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) then
            call cleanup()
            return
         end if
      end if
      
         
      do i=1,size(p%Jac_u_indx,1)
         
            ! get u_op + delta u
         call ED_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
         call ED_Perturb_u( p, i, 1, u_perturb, delta )

            ! compute x at u_op + delta u
         call ED_CalcContStateDeriv( t, u_perturb, p, x, xd, z, OtherState, m, x_p, ErrStat2, ErrMsg2 ) 
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            
                                         
            ! get u_op - delta u
         call ED_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
            
         call ED_Perturb_u( p, i, -1, u_perturb, delta )
         
            ! compute x at u_op - delta u
         call ED_CalcContStateDeriv( t, u_perturb, p, x, xd, z, OtherState, m, x_m, ErrStat2, ErrMsg2 ) 
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
            
            
            ! get central difference:            
            
            ! we may have had an error allocating memory, so we'll check
         if (ErrStat>=AbortErrLev) then 
            call cleanup()
            return
         end if         
         
         do j=1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs
            dXdu(j, i) = x_p%QT( p%DOFs%PS(j) ) - x_m%QT( p%DOFs%PS(j) )
         end do
         do j=1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs
            dXdu(j+p%DOFs%NActvDOF, i) = x_p%QDT( p%DOFs%PS(j) ) - x_m%QDT( p%DOFs%PS(j) )
         end do              
         dXdu(:,i) = dXdu(:,i) / (2*delta) 
         
      end do
      
      
      ! now do the extended input: sum the p%NumBl blade pitch columns
      dXdu(:,size(p%Jac_u_indx,1)+1) = dXdu(:,size(p%Jac_u_indx,1)-p%NumBl-1) ! last NumBl+2 columns are: GenTrq, YawMom, and BlPitchCom   
      do i=2,p%NumBl
         dXdu(:,size(p%Jac_u_indx,1)+1) = dXdu(:,size(p%Jac_u_indx,1)+1) + dXdu(:,size(p%Jac_u_indx,1)-p%NumBl-2+i) 
      end do
      
      
      call ED_DestroyContState( x_p, ErrStat2, ErrMsg2 ) ! we don't need this any more
      call ED_DestroyContState( x_m, ErrStat2, ErrMsg2 ) ! we don't need this any more
      
      
      
   END IF

   
   
   IF ( PRESENT( dXddu ) ) THEN
      if (allocated(dXddu)) deallocate(dXddu)
   END IF

   IF ( PRESENT( dZdu ) ) THEN
      if (allocated(dZdu)) deallocate(dZdu)
   END IF
   
   call cleanup()
   
contains
   subroutine cleanup()
      call ED_DestroyOutput(       y_p, ErrStat2, ErrMsg2 )
      call ED_DestroyOutput(       y_m, ErrStat2, ErrMsg2 )
      call ED_DestroyContState(    x_p, ErrStat2, ErrMsg2 )
      call ED_DestroyContState(    x_m, ErrStat2, ErrMsg2 )
      call ED_DestroyInput(  u_perturb, ErrStat2, ErrMsg2 )
      m%IgnoreMod = .false.
   end subroutine cleanup
   
END SUBROUTINE ED_JacobianPInput
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the continuous states (x). The partial derivatives dY/dx, dX/dx, dXd/dx, and dZ/dx are returned.
SUBROUTINE ED_JacobianPContState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdx, dXdx, dXddx, dZdx )
!..................................................................................................................................

   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(ED_InputType),                   INTENT(IN   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(ED_ParameterType),               INTENT(IN   )           :: p          !< Parameters
   TYPE(ED_ContinuousStateType),         INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(ED_DiscreteStateType),           INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(ED_ConstraintStateType),         INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(ED_OtherStateType),              INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(ED_OutputType),                  INTENT(INOUT)           :: y          !< Output (change to inout if a mesh copy is required);
                                                                               !!   Output fields are not used by this routine, but type is   
                                                                               !!   available here so that mesh parameter information (i.e.,  
                                                                               !!   connectivity) does not have to be recalculated for dYdu.
   TYPE(ED_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dYdx(:,:)  !< Partial derivatives of output functions (Y) with respect 
                                                                               !!   to the continuous states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXdx(:,:)  !< Partial derivatives of continuous state functions (X) with respect 
                                                                               !!   to the continuous states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXddx(:,:) !< Partial derivatives of discrete state functions (Xd) with respect 
                                                                               !!   to the continuous states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dZdx(:,:)  !< Partial derivatives of constraint state functions (Z) with respect 
                                                                               !!   to the continuous states (x) [intent in to avoid deallocation]
   
      ! local variables
   TYPE(ED_OutputType)                               :: y_p
   TYPE(ED_OutputType)                               :: y_m
   TYPE(ED_ContinuousStateType)                      :: x_p
   TYPE(ED_ContinuousStateType)                      :: x_m
   TYPE(ED_ContinuousStateType)                      :: x_perturb
   REAL(R8Ki)                                        :: delta        ! delta change in input or state
   INTEGER(IntKi)                                    :: i, j   
   
   INTEGER(IntKi)                                    :: ErrStat2
   CHARACTER(ErrMsgLen)                              :: ErrMsg2
   CHARACTER(*), PARAMETER                           :: RoutineName = 'ED_JacobianPContState'
   
   
      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''
   m%IgnoreMod = .true. ! to get true perturbations, we can't use the modulo function

      ! make a copy of the continuous states to perturb
   call ED_CopyContState( x, x_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat>=AbortErrLev) then
         call cleanup()
         return
      end if

   IF ( PRESENT( dYdx ) ) THEN

      ! Calculate the partial derivative of the output functions (Y) with respect to the continuous states (x) here:

      ! allocate dYdx if necessary
      if (.not. allocated(dYdx)) then
         call AllocAry(dYdx, p%Jac_ny, p%DOFs%NActvDOF*2, 'dYdx', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) then
            call cleanup()
            return
         end if
      end if
      
         ! make a copy of outputs because we will need two for the central difference computations (with orientations)
      call ED_CopyOutput( y, y_p, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call ED_CopyOutput( y, y_m, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) then
            call cleanup()
            return
         end if
         
         
      do i=1,p%DOFs%NActvDOF*2
         
            ! get x_op + delta x
         call ED_CopyContState( x, x_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
         call ED_Perturb_x( p, i, 1, x_perturb, delta )

            ! compute y at x_op + delta x
         call ED_CalcOutput( t, u, p, x_perturb, xd, z, OtherState, y_p, m, ErrStat2, ErrMsg2 ) 
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
         
            
            ! get x_op - delta x
         call ED_CopyContState( x, x_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later
         call ED_Perturb_x( p, i, -1, x_perturb, delta )
         
            ! compute y at x_op - delta x
         call ED_CalcOutput( t, u, p, x_perturb, xd, z, OtherState, y_m, m, ErrStat2, ErrMsg2 ) 
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
         
            
            ! get central difference:            
         call Compute_dY( p, y_p, y_m, delta, dYdx(:,i) )
         
      end do
      
      if (ErrStat>=AbortErrLev) then
         call cleanup()
         return
      end if
      call ED_DestroyOutput( y_p, ErrStat2, ErrMsg2 ) ! we don't need this any more
      call ED_DestroyOutput( y_m, ErrStat2, ErrMsg2 ) ! we don't need this any more
      
   END IF

   IF ( PRESENT( dXdx ) ) THEN

      ! Calculate the partial derivative of the continuous state functions (X) with respect to the continuous states (x) here:

      ! allocate dXdx if necessary
      if (.not. allocated(dXdx)) then
         call AllocAry(dXdx, p%DOFs%NActvDOF * 2, p%DOFs%NActvDOF * 2, 'dXdx', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) then
            call cleanup()
            return
         end if
      end if
      
      do i=1,p%DOFs%NActvDOF * 2
         
            ! get x_op + delta x
         call ED_CopyContState( x, x_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
         call ED_Perturb_x( p, i, 1, x_perturb, delta )

            ! compute x at x_op + delta x
         call ED_CalcContStateDeriv( t, u, p, x_perturb, xd, z, OtherState, m, x_p, ErrStat2, ErrMsg2 ) 
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            
                                         
            ! get x_op - delta x
         call ED_CopyContState( x, x_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
         call ED_Perturb_x( p, i, -1, x_perturb, delta )
         
            ! compute x at x_op - delta x
         call ED_CalcContStateDeriv( t, u, p, x_perturb, xd, z, OtherState, m, x_m, ErrStat2, ErrMsg2 ) 
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
                                         
            
            ! get central difference:            
            
            ! we may have had an error allocating memory, so we'll check
         if (ErrStat>=AbortErrLev) then 
            call cleanup()
            return
         end if         
         
         do j=1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs
            dXdx(j, i) = x_p%QT( p%DOFs%PS(j) ) - x_m%QT( p%DOFs%PS(j) )
         end do
         do j=1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs
            dXdx(j+p%DOFs%NActvDOF, i) = x_p%QDT( p%DOFs%PS(j) ) - x_m%QDT( p%DOFs%PS(j) )
         end do              
         dXdx(:,i) = dXdx(:,i) / (2*delta) 
         
      end do
      
      call ED_DestroyContState( x_p, ErrStat2, ErrMsg2 ) ! we don't need this any more
      call ED_DestroyContState( x_m, ErrStat2, ErrMsg2 ) ! we don't need this any more
   END IF

   IF ( PRESENT( dXddx ) ) THEN
      if (allocated(dXddx)) deallocate(dXddx)
   END IF

   IF ( PRESENT( dZdx ) ) THEN
      if (allocated(dZdx)) deallocate(dZdx)
   END IF

   call cleanup()
   
contains
   subroutine cleanup()
      call ED_DestroyOutput(         y_p, ErrStat2, ErrMsg2 )
      call ED_DestroyOutput(         y_m, ErrStat2, ErrMsg2 )
      call ED_DestroyContState(      x_p, ErrStat2, ErrMsg2 )
      call ED_DestroyContState(      x_m, ErrStat2, ErrMsg2 )
      call ED_DestroyContState(x_perturb, ErrStat2, ErrMsg2 )
      m%IgnoreMod = .false.
   end subroutine cleanup

END SUBROUTINE ED_JacobianPContState
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the discrete states (xd). The partial derivatives dY/dxd, dX/dxd, dXd/dxd, and dZ/dxd are returned.
SUBROUTINE ED_JacobianPDiscState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdxd, dXdxd, dXddxd, dZdxd )
!..................................................................................................................................

   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(ED_InputType),                   INTENT(INOUT)           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(ED_ParameterType),               INTENT(IN   )           :: p          !< Parameters
   TYPE(ED_ContinuousStateType),         INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(ED_DiscreteStateType),           INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(ED_ConstraintStateType),         INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(ED_OtherStateType),              INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(ED_OutputType),                  INTENT(INOUT)           :: y          !< Output (change to inout if a mesh copy is required);
                                                                               !!   Output fields are not used by this routine, but type is   
                                                                               !!   available here so that mesh parameter information (i.e.,  
                                                                               !!   connectivity) does not have to be recalculated for dYdu.
   TYPE(ED_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dYdxd(:,:) !< Partial derivatives of output functions
                                                                               !!  (Y) with respect to the discrete
                                                                               !!  states (xd) [intent in to avoid deallocation]
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXdxd(:,:) !< Partial derivatives of continuous state
                                                                               !!   functions (X) with respect to the
                                                                               !!   discrete states (xd) [intent in to avoid deallocation]
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXddxd(:,:)!< Partial derivatives of discrete state
                                                                               !!   functions (Xd) with respect to the
                                                                               !!   discrete states (xd) [intent in to avoid deallocation]
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dZdxd(:,:) !< Partial derivatives of constraint state
                                                                               !!   functions (Z) with respect to the
                                                                               !!   discrete states (xd) [intent in to avoid deallocation]


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''


   IF ( PRESENT( dYdxd ) ) THEN

      ! Calculate the partial derivative of the output functions (Y) with respect to the discrete states (xd) here:

      ! allocate and set dYdxd

   END IF

   IF ( PRESENT( dXdxd ) ) THEN

      ! Calculate the partial derivative of the continuous state functions (X) with respect to the discrete states (xd) here:

      ! allocate and set dXdxd

   END IF

   IF ( PRESENT( dXddxd ) ) THEN

      ! Calculate the partial derivative of the discrete state functions (Xd) with respect to the discrete states (xd) here:

      ! allocate and set dXddxd

   END IF

   IF ( PRESENT( dZdxd ) ) THEN

      ! Calculate the partial derivative of the constraint state functions (Z) with respect to the discrete states (xd) here:

      ! allocate and set dZdxd

   END IF


END SUBROUTINE ED_JacobianPDiscState
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the constraint states (z). The partial derivatives dY/dz, dX/dz, dXd/dz, and dZ/dz are returned.
SUBROUTINE ED_JacobianPConstrState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdz, dXdz, dXddz, dZdz )
!..................................................................................................................................

   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(ED_InputType),                   INTENT(INOUT)           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(ED_ParameterType),               INTENT(IN   )           :: p          !< Parameters
   TYPE(ED_ContinuousStateType),         INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(ED_DiscreteStateType),           INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(ED_ConstraintStateType),         INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(ED_OtherStateType),              INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(ED_OutputType),                  INTENT(INOUT)           :: y          !< Output (change to inout if a mesh copy is required);
                                                                               !!   Output fields are not used by this routine, but type is   
                                                                               !!   available here so that mesh parameter information (i.e.,  
                                                                               !!   connectivity) does not have to be recalculated for dYdu.
   TYPE(ED_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dYdz(:,:)  !< Partial derivatives of output functions (Y) with respect 
                                                                               !!  to the constraint states (z) [intent in to avoid deallocation]
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXdz(:,:)  !< Partial derivatives of continuous state functions (X) with respect 
                                                                               !!  to the constraint states (z) [intent in to avoid deallocation]
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXddz(:,:) !< Partial derivatives of discrete state functions (Xd) with respect 
                                                                               !!  to the constraint states (z) [intent in to avoid deallocation]
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dZdz(:,:)  !< Partial derivatives of constraint state functions (Z) with respect 
                                                                               !! to the constraint states (z) [intent in to avoid deallocation]


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''

   IF ( PRESENT( dYdz ) ) THEN

         ! Calculate the partial derivative of the output functions (Y) with respect to the constraint states (z) here:

      ! allocate and set dYdz

   END IF

   IF ( PRESENT( dXdz ) ) THEN

         ! Calculate the partial derivative of the continuous state functions (X) with respect to the constraint states (z) here:

      ! allocate and set dXdz

   END IF

   IF ( PRESENT( dXddz ) ) THEN

         ! Calculate the partial derivative of the discrete state functions (Xd) with respect to the constraint states (z) here:

      ! allocate and set dXddz

   END IF

   IF ( PRESENT( dZdz ) ) THEN

         ! Calculate the partial derivative of the constraint state functions (Z) with respect to the constraint states (z) here:

      ! allocate and set dZdz

   END IF


END SUBROUTINE ED_JacobianPConstrState
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine initializes the Jacobian parameters and initialization outputs for the linearized outputs.
SUBROUTINE ED_Init_Jacobian_y( p, y, InitOut, ErrStat, ErrMsg)

   TYPE(ED_ParameterType)            , INTENT(INOUT) :: p                     !< parameters
   TYPE(ED_OutputType)               , INTENT(IN   ) :: y                     !< outputs
   TYPE(ED_InitOutputType)           , INTENT(INOUT) :: InitOut               !< Output for initialization routine   
   
   INTEGER(IntKi)                    , INTENT(  OUT) :: ErrStat               !< Error status of the operation
   CHARACTER(*)                      , INTENT(  OUT) :: ErrMsg                !< Error message if ErrStat /= ErrID_None
   
      ! local variables:
   INTEGER(IntKi)                :: i,j,k, index_last, index_next
   INTEGER(IntKi)                                    :: ErrStat2
   CHARACTER(ErrMsgLen)                              :: ErrMsg2
   CHARACTER(*), PARAMETER                           :: RoutineName = 'ED_Init_Jacobian_y'
   LOGICAL                                           :: Mask(FIELDMASK_SIZE)   ! flags to determine if this field is part of the packing
   logical, allocatable                              :: AllOut(:)
   
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
      ! determine how many outputs there are in the Jacobians      
   p%Jac_ny = 0         
   if (allocated(y%BladeLn2Mesh)) then
      do i=1,p%NumBl
         p%Jac_ny = p%Jac_ny + y%BladeLn2Mesh(i)%NNodes * 18  ! 3 TranslationDisp, Orientation, TranslationVel, RotationVel, TranslationAcc, and RotationAcc at each node on each blade
      end do      
   end if
   
   p%Jac_ny = p%Jac_ny &
      + y%PlatformPtMesh%NNodes  * 18           & ! 3 TranslationDisp, Orientation, TranslationVel, RotationVel, TranslationAcc, and RotationAcc at each node
      + y%TowerLn2Mesh%NNodes    * 18           & ! 3 TranslationDisp, Orientation, TranslationVel, RotationVel, TranslationAcc, and RotationAcc at each node
      + y%HubPtMotion%NNodes     * 9            & ! 3 TranslationDisp, Orientation, and RotationVel at each node
      + y%NacelleMotion%NNodes   * 18           & ! 3 TranslationDisp, Orientation, TranslationVel, RotationVel, TranslationAcc, and RotationAcc at each node
      + 3                                       & ! Yaw, YawRate, and HSS_Spd
      + p%NumOuts  + p%BldNd_TotNumOuts           ! WriteOutput values 
      
   do i=1,p%NumBl
      p%Jac_ny = p%Jac_ny + y%BladeRootMotion(i)%NNodes * 18  ! 3 TranslationDisp, Orientation, TranslationVel, RotationVel, TranslationAcc, and RotationAcc at each (1) node on each blade
   end do

   
      !.................   
      ! set linearization output names:
      !.................   
   CALL AllocAry(InitOut%LinNames_y, p%Jac_ny, 'LinNames_y', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   CALL AllocAry(InitOut%RotFrame_y, p%Jac_ny, 'RotFrame_y', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return
   
   InitOut%RotFrame_y = .false. ! note that meshes are in the global, not rotating frame
   
   ! note that this Mask is for the y%HubPtMotion mesh ONLY. The others pack *all* of the motion fields
   Mask  = .false.
   Mask(MASKID_TRANSLATIONDISP) = .true.
   Mask(MASKID_ORIENTATION) = .true.
   Mask(MASKID_ROTATIONVEL) = .true.
   
   index_next = 1
   if (allocated(y%BladeLn2Mesh)) then
      index_last = index_next
      do i=1,p%NumBl
         call PackMotionMesh_Names(y%BladeLn2Mesh(i), 'Blade '//trim(num2lstr(i)), InitOut%LinNames_y, index_next)
      end do      
      !InitOut%RotFrame_y(index_last:index_next-1) = .true. ! values on the mesh are in global, not rotating frame
   end if
   call PackMotionMesh_Names(y%PlatformPtMesh, 'Platform', InitOut%LinNames_y, index_next)
   call PackMotionMesh_Names(y%TowerLn2Mesh, 'Tower', InitOut%LinNames_y, index_next)
   call PackMotionMesh_Names(y%HubPtMotion, 'Hub', InitOut%LinNames_y, index_next, FieldMask=Mask)
   index_last = index_next
   do i=1,p%NumBl
      call PackMotionMesh_Names(y%BladeRootMotion(i), 'Blade root '//trim(num2lstr(i)), InitOut%LinNames_y, index_next)
   end do   
   !InitOut%RotFrame_y(index_last:index_next-1) = .true. ! values on the mesh are in global, not rotating frame

   call PackMotionMesh_Names(y%NacelleMotion, 'Nacelle', InitOut%LinNames_y, index_next)
   InitOut%LinNames_y(index_next) = 'Yaw, rad'; index_next = index_next+1
   InitOut%LinNames_y(index_next) = 'YawRate, rad/s'; index_next = index_next+1
   InitOut%LinNames_y(index_next) = 'HSS_Spd, rad/s'
         
   do i=1,p%NumOuts + p%BldNd_TotNumOuts
      InitOut%LinNames_y(i+index_next) = trim(InitOut%WriteOutputHdr(i))//', '//trim(InitOut%WriteOutputUnt(i)) !trim(p%OutParam(i)%Name)//', '//p%OutParam(i)%Units
   end do   
   
   
   !! check for AllOuts in rotating frame
   allocate( AllOut(0:MaxOutPts), STAT=ErrStat2 ) ! allocate starting at zero to account for invalid output channels
   if (ErrStat2 /=0 ) then
      call SetErrStat(ErrID_Info, 'error allocating temporary space for AllOut',ErrStat,ErrMsg,RoutineName)
      return;
   end if
   
   AllOut = .false.
   do k=1,3
      AllOut(TipDxc(  k)) = .true.
      AllOut(TipDyc(  k)) = .true.
      AllOut(TipDzc(  k)) = .true.
      AllOut(TipDxb(  k)) = .true.
      AllOut(TipDyb(  k)) = .true.
      AllOut(TipALxb( k)) = .true.
      AllOut(TipALyb( k)) = .true.
      AllOut(TipALzb( k)) = .true.
      AllOut(TipRDxb( k)) = .true.
      AllOut(TipRDyb( k)) = .true.
      AllOut(TipRDzc( k)) = .true.
      AllOut(TipClrnc(k)) = .true.
      AllOut(PtchPMzc(k)) = .true.
      AllOut(RootFxc( k)) = .true.
      AllOut(RootFyc( k)) = .true.
      AllOut(RootFzc( k)) = .true.
      AllOut(RootFxb( k)) = .true.
      AllOut(RootFyb( k)) = .true.
      AllOut(RootMxc( k)) = .true.
      AllOut(RootMyc( k)) = .true.
      AllOut(RootMzc( k)) = .true.
      AllOut(RootMxb( k)) = .true.
      AllOut(RootMyb( k)) = .true.
      
      do j=1,9            
         AllOut(SpnALxb( j,k)) = .true.         
         AllOut(SpnALyb( j,k)) = .true.
         AllOut(SpnALzb( j,k)) = .true.
         AllOut(SpnFLxb( j,k)) = .true.
         AllOut(SpnFLyb( j,k)) = .true.
         AllOut(SpnFLzb( j,k)) = .true.
         AllOut(SpnMLxb( j,k)) = .true.
         AllOut(SpnMLyb( j,k)) = .true.
         AllOut(SpnMLzb( j,k)) = .true.
         AllOut(SpnTDxb( j,k)) = .true.
         AllOut(SpnTDyb( j,k)) = .true.
         AllOut(SpnTDzb( j,k)) = .true.
         AllOut(SpnRDxb( j,k)) = .true.
         AllOut(SpnRDyb( j,k)) = .true.
         AllOut(SpnRDzb( j,k)) = .true.
      end do
   end do
   
   do i=1,p%NumOuts
      InitOut%RotFrame_y(i+index_next) = AllOut( p%OutParam(i)%Indx )      
   end do    
   
   do i=1, p%BldNd_TotNumOuts
      InitOut%RotFrame_y(i+p%NumOuts+index_next) = .true.     
   end do
   
   deallocate(AllOut)         
   
   
END SUBROUTINE ED_Init_Jacobian_y
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine initializes the Jacobian parameters and initialization outputs for the linearized continuous states.
SUBROUTINE ED_Init_Jacobian_x( p, InitOut, ErrStat, ErrMsg)

   TYPE(ED_ParameterType)            , INTENT(INOUT) :: p                     !< parameters
   TYPE(ED_InitOutputType)           , INTENT(INOUT) :: InitOut               !< Output for initialization routine   
   
   INTEGER(IntKi)                    , INTENT(  OUT) :: ErrStat               !< Error status of the operation
   CHARACTER(*)                      , INTENT(  OUT) :: ErrMsg                !< Error message if ErrStat /= ErrID_None
   
   INTEGER(IntKi)                                    :: ErrStat2
   CHARACTER(ErrMsgLen)                              :: ErrMsg2
   CHARACTER(*), PARAMETER                           :: RoutineName = 'ED_Init_Jacobian_x'
   
      ! local variables:
   INTEGER(IntKi)                :: i
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
      ! allocate space for the row/column names and for perturbation sizes
   call allocAry(p%dx,               p%NDof,            'p%dx',       ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   CALL AllocAry(InitOut%LinNames_x, p%DOFs%NActvDOF*2, 'LinNames_x', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   CALL AllocAry(InitOut%RotFrame_x, p%DOFs%NActvDOF*2, 'RotFrame_x', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   CALL AllocAry(InitOut%DerivOrder_x, p%DOFs%NActvDOF*2, 'DerivOrder_x', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return
   
      ! All Elastodyn continuous states are max order = 2
   if ( allocated(InitOut%DerivOrder_x) ) InitOut%DerivOrder_x = 2
   
   p%dx = 0.0_R8Ki ! initialize in case we have only 1 blade
   
   ! set perturbation sizes: p%dx
   p%dx(DOF_Sg  :DOF_Hv)   = 0.2_R8Ki * D2R_D * max(p%TowerHt, 1.0_ReKi)     ! platform translational displacement states
   p%dx(DOF_R   :DOF_Y )   = 2.0_R8Ki * D2R_D                                ! platform rotational states
   p%dx(DOF_TFA1:DOF_TSS1) = 0.020_R8Ki * D2R_D * p%TwrFlexL                 ! tower deflection states: 1st tower
   p%dx(DOF_TFA2:DOF_TSS2) = 0.002_R8Ki * D2R_D * p%TwrFlexL                 ! tower deflection states: 2nd tower
   p%dx(DOF_Yaw :DOF_TFrl) = 2.0_R8Ki * D2R_D                                ! nacelle-yaw, rotor-furl, generator azimuth, drivetrain, and tail-furl rotational states

   do i=1,p%NumBl
      p%dx(DOF_BF(i,1))= 0.20_R8Ki * D2R_D * p%BldFlexL ! blade-deflection states: 1st blade flap mode 
      p%dx(DOF_BF(i,2))= 0.02_R8Ki * D2R_D * p%BldFlexL ! blade-deflection states: 2nd blade flap mode for blades (1/10 of the other perturbations)
      p%dx(DOF_BE(i,1))= 0.20_R8Ki * D2R_D * p%BldFlexL ! blade-deflection states: 1st blade edge mode
   end do
         
   if ( p%NumBl == 2 ) then
      p%dx(DOF_Teet)       = 2.0_R8Ki * D2R_D              ! rotor-teeter rotational state
   end if
   
   !Set some limits in case perturbation is very small
   do i=1,p%NDof
      p%dx(i) = max(p%dx(i), MinPerturb)
   end do
      
   InitOut%RotFrame_x   = .false.
   do i=1,p%DOFs%NActvDOF
      if (  p%DOFs%PS(i) >=  DOF_BF(1,1) ) then
         if ( p%NumBl == 2 ) then
            InitOut%RotFrame_x(i) = p%DOFs%PS(i) < DOF_Teet
         else
            InitOut%RotFrame_x(i) = .true. ! = p%DOFs%PS(i) <= DOF_BF (MaxBl,NumBF)
         end if
      end if      
   end do
   
      ! set linearization output names:
   do i=1,p%DOFs%NActvDOF
      InitOut%LinNames_x(i) = p%DOF_Desc( p%DOFs%PS(i) )
   end do
   
   do i=1,p%DOFs%NActvDOF
      InitOut%LinNames_x(i+p%DOFs%NActvDOF) = 'First time derivative of '//trim(InitOut%LinNames_x(i))//'/s'
      InitOut%RotFrame_x(i+p%DOFs%NActvDOF) = InitOut%RotFrame_x(i)
   end do      
   
END SUBROUTINE ED_Init_Jacobian_x
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine initializes the array that maps rows/columns of the Jacobian to specific mesh fields.
!! Do not change the order of this packing without changing corresponding linearization routines !
SUBROUTINE ED_Init_Jacobian( p, u, y, InitOut, ErrStat, ErrMsg)

   TYPE(ED_ParameterType)            , INTENT(INOUT) :: p                     !< parameters
   TYPE(ED_InputType)                , INTENT(IN   ) :: u                     !< inputs
   TYPE(ED_OutputType)               , INTENT(IN   ) :: y                     !< outputs
   TYPE(ED_InitOutputType)           , INTENT(INOUT) :: InitOut               !< Output for initialization routine   
   INTEGER(IntKi)                    , INTENT(  OUT) :: ErrStat               !< Error status of the operation
   CHARACTER(*)                      , INTENT(  OUT) :: ErrMsg                !< Error message if ErrStat /= ErrID_None
   
   INTEGER(IntKi)                                    :: ErrStat2
   CHARACTER(ErrMsgLen)                              :: ErrMsg2
   CHARACTER(*), PARAMETER                           :: RoutineName = 'ED_Init_Jacobian'
   
      ! local variables:
   INTEGER(IntKi)                :: i, j, k, index, index_last, nu, i_meshField, m
   REAL(R8Ki)                    :: MaxThrust, MaxTorque
   REAL(R8Ki)                    :: ScaleLength
   
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
         
   
   call ED_Init_Jacobian_y( p, y, InitOut, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         
   call ED_Init_Jacobian_x( p, InitOut, ErrStat2, ErrMsg2)      
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   
      
      
      ! determine how many inputs there are in the Jacobians
   nu = 0;
   if (allocated(u%BladePtLoads)) then
      do i=1,p%NumBl
         nu = nu + u%BladePtLoads(i)%NNodes * 6  ! 3 forces + 3 moments at each node on each blade
      end do      
   end if
   nu = nu &
      + u%PlatformPtMesh%NNodes * 6            & ! 3 forces + 3 moments at each node
      + u%TowerPtLoads%NNodes   * 6            & ! 3 forces + 3 moments at each node
      + u%HubPtLoad%NNodes      * 6            & ! 3 forces + 3 moments at each node
      + u%NacelleLoads%NNodes   * 6            & ! 3 forces + 3 moments at each node
      + p%NumBl                                & ! blade pitch command (BlPitchCom)    
      + 2                                        ! YawMom and GenTrq
         
   ! note: all other inputs are ignored
      
   !....................                        
   ! fill matrix to store index to help us figure out what the ith value of the u vector really means
   ! (see elastodyn::ed_perturb_u ... these MUST match )
   ! column 1 indicates module's mesh and field
   ! column 2 indicates the first index of the acceleration/load field
   ! column 3 is the node
   !....................
      
   !...............
   ! ED input mappings stored in p%Jac_u_indx:   
   !...............
   call AllocAry(p%Jac_u_indx, nu, 3, 'p%Jac_u_indx', ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)   
   if (ErrStat >= AbortErrLev) return
     
      
   index = 1
   if (allocated(u%BladePtLoads)) then
      !Module/Mesh/Field: u%BladePtLoads(1)%Force  = 1;
      !Module/Mesh/Field: u%BladePtLoads(1)%Moment = 2;
      !Module/Mesh/Field: u%BladePtLoads(2)%Force  = 3;
      !Module/Mesh/Field: u%BladePtLoads(2)%Moment = 4;
      !Module/Mesh/Field: u%BladePtLoads(3)%Force  = 5;
      !Module/Mesh/Field: u%BladePtLoads(3)%Moment = 6;      
      do k=1,p%NumBl
         
         do i_meshField = 1,2
            do i=1,u%BladePtLoads(k)%NNodes
               do j=1,3
                  p%Jac_u_indx(index,1) =  i_meshField + (k-1)*2 !Module/Mesh/Field: u%BladePtLoads(k)%{Force/Moment} = m
                  p%Jac_u_indx(index,2) =  j !index:  j
                  p%Jac_u_indx(index,3) =  i !Node:   i
                  index = index + 1
               end do !j      
            end do !i
            
         end do !i_meshField                            
      end do !k
                        
   end if
   
   !if MaxBl ever changes (i.e., MaxBl /=3), we need to modify this accordingly:
   do i_meshField = 7,8
      do i=1,u%PlatformPtMesh%NNodes
         do j=1,3
            p%Jac_u_indx(index,1) =  i_meshField !Module/Mesh/Field: u%PlatformPtMesh%Force = 7; u%PlatformPtMesh%Moment = 8;
            p%Jac_u_indx(index,2) =  j !index:  j
            p%Jac_u_indx(index,3) =  i !Node:   i
            index = index + 1
         end do !j      
      end do !i
   end do
  
   do i_meshField = 9,10
      do i=1,u%TowerPtLoads%NNodes
         do j=1,3
            p%Jac_u_indx(index,1) =  i_meshField !Module/Mesh/Field: u%TowerPtLoads%Force = 9; u%TowerPtLoads%Moment = 10;
            p%Jac_u_indx(index,2) =  j !index:  j
            p%Jac_u_indx(index,3) =  i !Node:   i
            index = index + 1
         end do !j      
      end do !i
   end do
         
   do i_meshField = 11,12
      do i=1,u%HubPtLoad%NNodes
         do j=1,3
            p%Jac_u_indx(index,1) =  i_meshField !Module/Mesh/Field: u%HubPtLoad%Force = 11; u%HubPtLoad%Moment = 12;
            p%Jac_u_indx(index,2) =  j !index:  j
            p%Jac_u_indx(index,3) =  i !Node:   i
            index = index + 1
         end do !j      
      end do !i
   end do   
   
   do i_meshField = 13,14
      do i=1,u%NacelleLoads%NNodes
         do j=1,3
            p%Jac_u_indx(index,1) =  i_meshField !Module/Mesh/Field: u%NacelleLoads%Force = 13; u%NacelleLoads%Moment = 14;
            p%Jac_u_indx(index,2) =  j !index:  j
            p%Jac_u_indx(index,3) =  i !Node:   i
            index = index + 1
         end do !j      
      end do !i
   end do
   
   do i_meshField = 1,p%NumBl ! scalars   
      p%Jac_u_indx(index,1) =  15 !Module/Mesh/Field: u%BlPitchCom = 15;
      p%Jac_u_indx(index,2) =  1 !index:  n/a
      p%Jac_u_indx(index,3) =  i_meshField !Node:   blade
      index = index + 1      
   end do
   
   do i_meshField = 16,17 ! scalars   
      p%Jac_u_indx(index,1) =  i_meshField !Module/Mesh/Field: u%YawMom = 16; u%GenTrq = 17;
      p%Jac_u_indx(index,2) =  1 !index:  j
      p%Jac_u_indx(index,3) =  1 !Node:   i
      index = index + 1
   end do
   
   !................
   ! input perturbations, du:
   !................
   call AllocAry(p%du, 17, 'p%du', ErrStat2, ErrMsg2) ! 17 = number of unique values in p%Jac_u_indx(:,1) 
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)      
   if (ErrStat >= AbortErrLev) return
   
      ! p%TipRad is set to 0 for BeamDyn simulations, so we're using a copy of the value from the input file here
   ScaleLength = max(p%TipRad, p%TowerHt, 1.0_ReKi)
   MaxThrust = 490.0_R8Ki * pi_D /  9.0_R8Ki * ScaleLength**2
   MaxTorque = 122.5_R8Ki * pi_D / 27.0_R8Ki * ScaleLength**3
   
   if (allocated(u%BladePtLoads)) then
      do k=1,p%NumBl
         p%du(2*k-1) = MaxThrust / real(100*p%NumBl*u%BladePtLoads(k)%NNodes,R8Ki) ! u%BladePtLoads(k)%Force  = 2*k-1
         p%du(2*k  ) = MaxTorque / real(100*p%NumBl*u%BladePtLoads(k)%NNodes,R8Ki) ! u%BladePtLoads(k)%Moment = 2*k
      end do !k
   else
      p%du(1:6) = 0.0_R8Ki
   end if
   
   p%du( 7) = MaxThrust / 100.0_R8Ki                           ! u%PlatformPtMesh%Force = 7
   p%du( 8) = MaxTorque / 100.0_R8Ki                           ! u%PlatformPtMesh%Moment = 8
   p%du( 9) = MaxThrust / real(100*u%TowerPtLoads%NNodes,R8Ki) ! u%TowerPtLoads%Force = 9
   p%du(10) = MaxTorque / real(100*u%TowerPtLoads%NNodes,R8Ki) ! u%TowerPtLoads%Moment = 10
   p%du(11) = MaxThrust / 100.0_R8Ki                           ! u%HubPtLoad%Force = 11
   p%du(12) = MaxTorque / 100.0_R8Ki                           ! u%HubPtLoad%Moment = 12
   p%du(13) = MaxThrust / 100.0_R8Ki                           ! u%NacelleLoads%Force = 13
   p%du(14) = MaxTorque / 100.0_R8Ki                           ! u%NacelleLoads%Moment = 14   
   p%du(15) = 2.0_R8Ki * D2R_D                                 ! u%BlPitchCom = 15 
   p%du(16) = MaxTorque / 100.0_R8Ki                           ! u%YawMom = 16
   p%du(17) = MaxTorque / (100.0_R8Ki*p%GBRatio)               ! u%GenTrq = 17
      
   !Set some limits in case perturbation is very small
   do i=1,size(p%du)
      p%du(i) = max(p%du(i), MinPerturb)
   end do

   !................
   ! names of the columns, InitOut%LinNames_u:
   !................
   call AllocAry(InitOut%LinNames_u, nu+1, 'LinNames_u', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   call AllocAry(InitOut%RotFrame_u, nu+1, 'RotFrame_u', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   call AllocAry(InitOut%IsLoad_u,   nu+1, 'IsLoad_u',   ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return
      
   InitOut%IsLoad_u   = .true.  ! most of ED's inputs are loads; we will override the non-load inputs below.
   InitOut%RotFrame_u = .false.
   index = 1
   if (allocated(u%BladePtLoads)) then
      index_last = index
      do k=1,p%NumBl
         call PackLoadMesh_Names(u%BladePtLoads(k), 'Blade '//trim(num2lstr(k)), InitOut%LinNames_u, index)   
      end do
      !InitOut%RotFrame_u(index_last:index-1) = .true. ! values on the mesh are in global, not rotating frame
   end if
   call PackLoadMesh_Names(u%PlatformPtMesh, 'Platform', InitOut%LinNames_u, index)   
   call PackLoadMesh_Names(u%TowerPtLoads, 'Tower', InitOut%LinNames_u, index)   
   call PackLoadMesh_Names(u%HubPtLoad, 'Hub', InitOut%LinNames_u, index)   
   call PackLoadMesh_Names(u%NacelleLoads, 'Nacelle', InitOut%LinNames_u, index)   
      
   do k = 1,p%NumBl ! scalars
      InitOut%LinNames_u(index) = 'Blade '//trim(num2lstr(k))//' pitch command, rad'
      InitOut%IsLoad_u(  index) = .false.
      InitOut%RotFrame_u(index) = .true.
      index = index + 1
   end do

   InitOut%LinNames_u(index) = 'Yaw moment, Nm' ; index = index + 1
   InitOut%LinNames_u(index) = 'Generator torque, Nm' ; index = index + 1
   InitOut%LinNames_u(index) = 'Extended input: collective blade-pitch command, rad'
   InitOut%IsLoad_u(  index) = .false.
   
END SUBROUTINE ED_Init_Jacobian
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine perturbs the nth element of the u array (and mesh/field it corresponds to)
!! Do not change this without making sure subroutine elastodyn::ed_init_jacobian is consistant with this routine!
SUBROUTINE ED_Perturb_u( p, n, perturb_sign, u, du )

   TYPE(ED_ParameterType)              , INTENT(IN   ) :: p                      !< parameters
   INTEGER( IntKi )                    , INTENT(IN   ) :: n                      !< number of array element to use 
   INTEGER( IntKi )                    , INTENT(IN   ) :: perturb_sign           !< +1 or -1 (value to multiply perturbation by; positive or negative difference)
   TYPE(ED_InputType)                  , INTENT(INOUT) :: u                      !< perturbed ED inputs
   REAL( R8Ki )                        , INTENT(  OUT) :: du                     !< amount that specific input was perturbed
   

   ! local variables
   INTEGER                                             :: fieldIndx
   INTEGER                                             :: node
   
      
   fieldIndx = p%Jac_u_indx(n,2) 
   node      = p%Jac_u_indx(n,3) 
   
   du = p%du(  p%Jac_u_indx(n,1) )
   
      ! determine which mesh we're trying to perturb and perturb the input:
   SELECT CASE( p%Jac_u_indx(n,1) )
      
   CASE ( 1) !Module/Mesh/Field: u%BladePtLoads(1)%Force = 1      
      u%BladePtLoads(1)%Force( fieldIndx,node) = u%BladePtLoads(1)%Force( fieldIndx,node) + du * perturb_sign       
   CASE ( 2) !Module/Mesh/Field: u%BladePtLoads(1)%Moment = 2
      u%BladePtLoads(1)%Moment(fieldIndx,node) = u%BladePtLoads(1)%Moment(fieldIndx,node) + du * perturb_sign            
   CASE ( 3) !Module/Mesh/Field: u%BladePtLoads(2)%Force = 3
      u%BladePtLoads(2)%Force( fieldIndx,node) = u%BladePtLoads(2)%Force( fieldIndx,node) + du * perturb_sign       
   CASE ( 4) !Module/Mesh/Field: u%BladePtLoads(2)%Moment = 4
      u%BladePtLoads(2)%Moment(fieldIndx,node) = u%BladePtLoads(2)%Moment(fieldIndx,node) + du * perturb_sign            
   CASE ( 5) !Module/Mesh/Field: u%BladePtLoads(2)%Force = 5
      u%BladePtLoads(3)%Force( fieldIndx,node) = u%BladePtLoads(3)%Force( fieldIndx,node) + du * perturb_sign       
   CASE ( 6) !Module/Mesh/Field: u%BladePtLoads(2)%Moment = 6
      u%BladePtLoads(3)%Moment(fieldIndx,node) = u%BladePtLoads(3)%Moment(fieldIndx,node) + du * perturb_sign            
               
   CASE ( 7) !Module/Mesh/Field: u%PlatformPtMesh%Force = 7
      u%PlatformPtMesh%Force( fieldIndx,node) = u%PlatformPtMesh%Force( fieldIndx,node) + du * perturb_sign       
   CASE ( 8) !Module/Mesh/Field: u%PlatformPtMesh%Moment = 8
      u%PlatformPtMesh%Moment(fieldIndx,node) = u%PlatformPtMesh%Moment(fieldIndx,node) + du * perturb_sign            
                     
   CASE ( 9) !Module/Mesh/Field: u%TowerPtLoads%Force = 9
      u%TowerPtLoads%Force( fieldIndx,node) = u%TowerPtLoads%Force( fieldIndx,node) + du * perturb_sign       
   CASE (10) !Module/Mesh/Field: u%TowerPtLoads%Moment = 10
      u%TowerPtLoads%Moment(fieldIndx,node) = u%TowerPtLoads%Moment(fieldIndx,node) + du * perturb_sign            

   CASE (11) !Module/Mesh/Field: u%HubPtLoad%Force = 11
      u%HubPtLoad%Force( fieldIndx,node) = u%HubPtLoad%Force( fieldIndx,node) + du * perturb_sign       
   CASE (12) !Module/Mesh/Field: u%HubPtLoad%Moment = 12
      u%HubPtLoad%Moment(fieldIndx,node) = u%HubPtLoad%Moment(fieldIndx,node) + du * perturb_sign            
  
   CASE (13) !Module/Mesh/Field: u%NacelleLoads%Force = 13
      u%NacelleLoads%Force( fieldIndx,node) = u%NacelleLoads%Force( fieldIndx,node) + du * perturb_sign       
   CASE (14) !Module/Mesh/Field: u%NacelleLoads%Moment = 14
      u%NacelleLoads%Moment(fieldIndx,node) = u%NacelleLoads%Moment(fieldIndx,node) + du * perturb_sign            
   
   CASE (15) !Module/Mesh/Field: u%BlPitchCom = 15
      u%BlPitchCom(node) = u%BlPitchCom(node) + du * perturb_sign
   CASE (16) !Module/Mesh/Field: u%YawMom = 16
      u%YawMom = u%YawMom + du * perturb_sign
   CASE (17) !Module/Mesh/Field: u%GenTrq = 17
      u%GenTrq = u%GenTrq + du * perturb_sign
      
   END SELECT
                                             
END SUBROUTINE ED_Perturb_u
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine perturbs the nth element of the continuous state array.
!! Do not change this without making sure subroutine elastodyn::ed_init_jacobian is consistant with this routine!
SUBROUTINE ED_Perturb_x( p, n, perturb_sign, x, dx )

   TYPE(ED_ParameterType)              , INTENT(IN   ) :: p                      !< parameters
   INTEGER( IntKi )                    , INTENT(IN   ) :: n                      !< number of array element to use 
   INTEGER( IntKi )                    , INTENT(IN   ) :: perturb_sign           !< +1 or -1 (value to multiply perturbation by; positive or negative difference)
   TYPE(ED_ContinuousStateType)        , INTENT(INOUT) :: x                      !< perturbed ED states
   REAL( R8Ki )                        , INTENT(  OUT) :: dx                     !< amount that specific state was perturbed
   

   ! local variables
   integer(intKi)                                      :: indx
   
   
   if (n > p%DOFs%NActvDOF) then
      indx = p%DOFs%PS(n-p%DOFs%NActvDOF)
      dx   = p%dx( indx )

      x%QDT( indx ) = x%QDT( indx ) + dx * perturb_sign 
   else
      indx = p%DOFs%PS(n)
      dx   = p%dx( indx )
      
      x%QT(  indx ) = x%QT(  indx ) + dx * perturb_sign 
   end if
                                                
END SUBROUTINE ED_Perturb_x
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine uses values of two output types to compute an array of differences.
!! Do not change this packing without making sure subroutine elastodyn::ed_init_jacobian is consistant with this routine!
SUBROUTINE Compute_dY(p, y_p, y_m, delta, dY)
   
   TYPE(ED_ParameterType)            , INTENT(IN   ) :: p         !< parameters
   TYPE(ED_OutputType)               , INTENT(IN   ) :: y_p       !< ED outputs at \f$ u + \Delta u \f$ or \f$ x + \Delta x \f$ (p=plus)
   TYPE(ED_OutputType)               , INTENT(IN   ) :: y_m       !< ED outputs at \f$ u - \Delta u \f$ or \f$ x - \Delta x \f$ (m=minus)   
   REAL(R8Ki)                        , INTENT(IN   ) :: delta     !< difference in inputs or states \f$ delta = \Delta u \f$ or \f$ delta = \Delta x \f$
   REAL(R8Ki)                        , INTENT(INOUT) :: dY(:)     !< column of dYdu or dYdx: \f$ \frac{\partial Y}{\partial u_i} = \frac{y_p - y_m}{2 \, \Delta u}\f$ or \f$ \frac{\partial Y}{\partial x_i} = \frac{y_p - y_m}{2 \, \Delta x}\f$
   
      ! local variables:
   INTEGER(IntKi)                                    :: k                      ! loop over blades
   INTEGER(IntKi)                                    :: indx_first             ! index indicating next value of dY to be filled 
   LOGICAL                                           :: Mask(FIELDMASK_SIZE)   ! flags to determine if this field is part of the packing

   
   Mask  = .false.
   Mask(MASKID_TRANSLATIONDISP) = .true.
   Mask(MASKID_ORIENTATION) = .true.
   Mask(MASKID_ROTATIONVEL) = .true.   
   
   
   indx_first = 1            
   if (allocated(y_p%BladeLn2Mesh)) then
      do k=1,p%NumBl
         call PackMotionMesh_dY(y_p%BladeLn2Mesh(k), y_m%BladeLn2Mesh(k), dY, indx_first)                  
      end do      
   end if
   
   call PackMotionMesh_dY(y_p%PlatformPtMesh, y_m%PlatformPtMesh, dY, indx_first, UseSmlAngle=.true.)
   call PackMotionMesh_dY(y_p%TowerLn2Mesh,   y_m%TowerLn2Mesh,   dY, indx_first, UseSmlAngle=.true.)
   call PackMotionMesh_dY(y_p%HubPtMotion,    y_m%HubPtMotion,    dY, indx_first, FieldMask=Mask)
   do k=1,p%NumBl
      call PackMotionMesh_dY(y_p%BladeRootMotion(k),   y_m%BladeRootMotion(k),   dY, indx_first)                  
   end do
   call PackMotionMesh_dY(y_p%NacelleMotion,  y_m%NacelleMotion,  dY, indx_first)                  
                     
   dY(indx_first) = y_p%Yaw     - y_m%Yaw;       indx_first = indx_first + 1    
   dY(indx_first) = y_p%YawRate - y_m%YawRate;   indx_first = indx_first + 1    
   dY(indx_first) = y_p%HSS_Spd - y_m%HSS_Spd;   indx_first = indx_first + 1
   
   !indx_last = indx_first + p%NumOuts - 1
   do k=1,p%NumOuts + p%BldNd_TotNumOuts
      dY(k+indx_first-1) = y_p%WriteOutput(k) - y_m%WriteOutput(k)
   end do   
   
   dY = dY / (2.0_R8Ki*delta)
   
END SUBROUTINE Compute_dY
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to pack the data structures representing the operating points into arrays for linearization.
SUBROUTINE ED_GetOP( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, u_op, y_op, x_op, dx_op, xd_op, z_op, NeedTrimOP )

   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(ED_InputType),                   INTENT(IN   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(ED_ParameterType),               INTENT(IN   )           :: p          !< Parameters
   TYPE(ED_ContinuousStateType),         INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(ED_DiscreteStateType),           INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(ED_ConstraintStateType),         INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(ED_OtherStateType),              INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(ED_OutputType),                  INTENT(IN   )           :: y          !< Output at operating point
   TYPE(ED_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: u_op(:)    !< values of linearized inputs
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: y_op(:)    !< values of linearized outputs
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: x_op(:)    !< values of linearized continuous states
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dx_op(:)   !< values of first time derivatives of linearized continuous states
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: xd_op(:)   !< values of linearized discrete states
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: z_op(:)    !< values of linearized constraint states
   LOGICAL,                 OPTIONAL,    INTENT(IN   )           :: NeedTrimOP !< whether a y_op values should contain values for trim solution (3-value representation instead of full orientation matrices, no rotation acc)



   INTEGER(IntKi)                                    :: i, k, index
   INTEGER(IntKi)                                    :: ny
   INTEGER(IntKi)                                    :: ErrStat2
   CHARACTER(ErrMsgLen)                              :: ErrMsg2
   CHARACTER(*), PARAMETER                           :: RoutineName = 'ED_GetOP'
   LOGICAL                                           :: ReturnTrimOP
   TYPE(ED_ContinuousStateType)                      :: dx          !< derivative of continuous states at operating point
   LOGICAL                                           :: Mask(FIELDMASK_SIZE)               !< flags to determine if this field is part of the packing
   
   
      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''

   !..................................
   IF ( PRESENT( u_op ) ) THEN
      if (.not. allocated(u_op)) then         
         call AllocAry(u_op, size(p%Jac_u_indx,1)+1,'u_op',ErrStat2,ErrMsg2) ! +1 for extended input here
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) return
      end if
            
      index = 1
      if (allocated(u%BladePtLoads)) then
         do k=1,p%NumBl
            call PackLoadMesh(u%BladePtLoads(k), u_op, index)   
         end do
      end if
      call PackLoadMesh(u%PlatformPtMesh, u_op, index)   
      call PackLoadMesh(u%TowerPtLoads, u_op, index)   
      call PackLoadMesh(u%HubPtLoad, u_op, index)   
      call PackLoadMesh(u%NacelleLoads, u_op, index)   
      
      do k = 1,p%NumBl ! scalars
         u_op(index) = u%BlPitchCom(k)
         index = index + 1
      end do
      u_op(index) = u%YawMom ; index = index + 1
      u_op(index) = u%GenTrq ; index = index + 1
      
         ! extended input:
      u_op(index) = u%BlPitchCom(1)
      
      do k = 2,p%NumBl
         if (.not. EqualRealNos( u%BlPitchCom(1), u%BlPitchCom(k) ) ) then
            call SetErrStat(ErrID_Info,"Operating point of collective pitch extended input is invalid because "// &
                     "the commanded blade pitch angles are not the same for each blade.", ErrStat, ErrMsg, RoutineName)
            exit
         end if      
      end do      
      
   END IF

   !..................................
   IF ( PRESENT( y_op ) ) THEN
      if (present(NeedTrimOP)) then
         ReturnTrimOP = NeedTrimOP
      else
         ReturnTrimOP = .false.
      end if
      
      if (.not. allocated(y_op)) then 
            ! our operating point includes DCM (orientation) matrices, not just small angles like the perturbation matrices do
         ny = p%Jac_ny + y%PlatformPtMesh%NNodes * 6 & ! Jac_ny has 3 for Orientation, but we need 9 at each node
                       + y%TowerLn2Mesh%NNodes   * 6 & ! Jac_ny has 3 for Orientation, but we need 9 at each node 
                       + y%HubPtMotion%NNodes    * 6 & ! Jac_ny has 3 for Orientation, but we need 9 at each node
                       + y%NacelleMotion%NNodes  * 6   ! Jac_ny has 3 for Orientation, but we need 9 at each node
            
         if (allocated(y%BladeLn2Mesh)) then
            do k=1,p%NumBl
               ny = ny + y%BladeLn2Mesh(k)%NNodes * 6  ! Jac_ny has 3 for Orientation, but we need 9 (at each node on each blade)
            end do      
         end if
         do k=1,p%NumBl
            ny = ny + y%BladeRootMotion(k)%NNodes * 6  ! Jac_ny has 3 for Orientation, but we need 9 at each node on each blade
         end do
                                    
         call AllocAry(y_op, ny,'y_op',ErrStat2,ErrMsg2)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) return
      end if
            
      if (ReturnTrimOP) y_op = 0.0_ReKi ! initialize in case we are returning packed orientations and don't fill the entire array

      
      Mask  = .false.
      Mask(MASKID_TRANSLATIONDISP) = .true.
      Mask(MASKID_ORIENTATION) = .true.
      Mask(MASKID_ROTATIONVEL) = .true.
   
      index = 1
      if (allocated(y%BladeLn2Mesh)) then
         do k=1,p%NumBl
            call PackMotionMesh(y%BladeLn2Mesh(k), y_op, index, TrimOP=ReturnTrimOP)
         end do      
      end if
      call PackMotionMesh(y%PlatformPtMesh, y_op, index, TrimOP=ReturnTrimOP)
      call PackMotionMesh(y%TowerLn2Mesh, y_op, index, TrimOP=ReturnTrimOP)
      call PackMotionMesh(y%HubPtMotion, y_op, index, FieldMask=Mask, TrimOP=ReturnTrimOP)
      
      do k=1,p%NumBl
         call PackMotionMesh(y%BladeRootMotion(k), y_op, index, TrimOP=ReturnTrimOP)
      end do
      call PackMotionMesh(y%NacelleMotion, y_op, index, TrimOP=ReturnTrimOP)
      
      y_op(index) = y%Yaw     ; index = index + 1    
      y_op(index) = y%YawRate ; index = index + 1    
      y_op(index) = y%HSS_Spd 
   
      do i=1,p%NumOuts + p%BldNd_TotNumOuts
         y_op(i+index) = y%WriteOutput(i)
      end do
                        
   END IF

   !..................................
   IF ( PRESENT( x_op ) ) THEN

      if (.not. allocated(x_op)) then                           
         call AllocAry(x_op, p%DOFs%NActvDOF * 2,'x_op',ErrStat2,ErrMsg2)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) return
      end if
      
      do i=1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs
         x_op(i) = x%QT( p%DOFs%PS(i) )
      end do
      do i=1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs
         x_op(i+p%DOFs%NActvDOF) = x%QDT( p%DOFs%PS(i) )
      end do                                
      
   END IF

   !..................................
   IF ( PRESENT( dx_op ) ) THEN

      if (.not. allocated(dx_op)) then                           
         call AllocAry(dx_op, p%DOFs%NActvDOF * 2,'dx_op',ErrStat2,ErrMsg2)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) return
      end if
      
      call ED_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dx, ErrStat2, ErrMsg2 ) 
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
         if (ErrStat>=AbortErrLev) then
            call ED_DestroyContState( dx, ErrStat2, ErrMsg2)
            return
         end if
                     
      do i=1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs
         dx_op(i) = dx%QT( p%DOFs%PS(i) )
      end do
      do i=1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs
         dx_op(i+p%DOFs%NActvDOF) = dx%QDT( p%DOFs%PS(i) )
      end do                                
      
      call ED_DestroyContState( dx, ErrStat2, ErrMsg2)
            
   END IF

   !..................................
   IF ( PRESENT( xd_op ) ) THEN
   END IF
   
   !..................................
   IF ( PRESENT( z_op ) ) THEN
   END IF

END SUBROUTINE ED_GetOP
!----------------------------------------------------------------------------------------------------------------------------------


END MODULE ElastoDyn
!**********************************************************************************************************************************
