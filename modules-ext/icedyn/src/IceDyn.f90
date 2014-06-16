!**********************************************************************************************************************************
! File last committed: $Date$
! (File) Revision #: $Rev$
! URL: $HeadURL$
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2014  National Renewable Energy Laboratory
!
!    This file is part of module IceDyn.
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
!**********************************************************************************************************************************
!    IceDyn is a module describing ice load on offshore wind turbine supporting structures.
!
!    The module is given the module name ModuleName = IceDyn and the abbreviated name ModName = ID. The mathematical
!    formulation of this module is a subset of the most general form permitted by the FAST modularization framework in tight
!    coupling, thus, the module is developed to support both loose and tight coupling (tight coupling for both time marching and
!    linearization).
!
!
!    References:
!
!    Ice Module Manual, by Bingbin Yu, Dale Karr.
!
!**********************************************************************************************************************************
MODULE IceDyn

   USE IceDyn_Types
   USE NWTC_Library

   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER  :: ID_Ver = ProgDesc( 'IceDyn', 'v1.00.01', '27-Oct-2013' )

   ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: ID_Init                           ! Initialization routine
   PUBLIC :: ID_End                            ! Ending routine (includes clean up)

   PUBLIC :: ID_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                               !   continuous states, and updating discrete states
   PUBLIC :: ID_CalcOutput                     ! Routine for computing outputs

   PUBLIC :: ID_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: ID_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: ID_UpdateDiscState                ! Tight coupling routine for updating discrete states

CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_Init( InitInp, u, p, x, xd, z, OtherState, y, Interval, InitOut, ErrStat, ErrMsg )
! This routine is called at the start of the simulation to perform initialization steps.
! The parameters are set here and not changed during the simulation.
! The initial states and initial guess for the input are defined.
!..................................................................................................................................

   TYPE(ID_InitInputType),       INTENT(IN   )  :: InitInp     ! Input data for initialization routine
   TYPE(ID_InputType),           INTENT(  OUT)  :: u           ! An initial guess for the input; input mesh must be defined
   TYPE(ID_ParameterType),       INTENT(  OUT)  :: p           ! Parameters
   TYPE(ID_ContinuousStateType), INTENT(  OUT)  :: x           ! Initial continuous states
   TYPE(ID_DiscreteStateType),   INTENT(  OUT)  :: xd          ! Initial discrete states
   TYPE(ID_ConstraintStateType), INTENT(  OUT)  :: z           ! Initial guess of the constraint states
   TYPE(ID_OtherStateType),      INTENT(  OUT)  :: OtherState  ! Initial other/optimization states
   TYPE(ID_OutputType),          INTENT(  OUT)  :: y           ! Initial system outputs (outputs are not calculated;
                                                               !   only the output mesh is initialized)
   REAL(DbKi),                   INTENT(INOUT)  :: Interval    ! Coupling interval in seconds: the rate that
                                                               !   (1) ID_UpdateStates() is called in loose coupling &
                                                               !   (2) ID_UpdateDiscState() is called in tight coupling.
                                                               !   Input is the suggested time from the glue code;
                                                               !   Output is the actual coupling interval that will be used
                                                               !   by the glue code.
   TYPE(ID_InitOutputType),      INTENT(  OUT)  :: InitOut     ! Output for initialization routine
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


      ! Local variables

   TYPE(ID_InputFile)                           :: InputFileData           ! Data stored in the module's input file
   INTEGER(IntKi)                               :: ErrStat2                ! temporary Error status of the operation
   CHARACTER(LEN(ErrMsg))                       :: ErrMsg2                 ! temporary Error message if ErrStat /= ErrID_None
   REAL(ReKi)                                   :: TmpPos(3)


      ! Initialize variables for this routine

   ErrStat = ErrID_None
   ErrMsg  = ""


      ! Initialize the NWTC Subroutine Library

   CALL NWTC_Init( EchoLibVer=.FALSE. )

      ! Display the module information

   CALL DispNVD( ID_Ver )

      !............................................................................................
      ! Read the input file and validate the data
      !............................................................................................
   p%RootName = TRIM(InitInp%RootName)//'.ID' ! all of the output file names from this module will end with '.ID'

   CALL ID_ReadInput( InitInp, InputFileData, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   ! CALL ID_ValidateInput( InputFileData, ErrStat2, ErrMsg2 )
   !   CALL CheckError( ErrStat2, ErrMsg2 )
   !   IF (ErrStat >= AbortErrLev) RETURN

      !............................................................................................
      ! Define parameters here:
      !............................................................................................
   CALL ID_SetParameters( InputFileData, p, Interval, InitInp%Tmax, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN

      p%verif = 1
      p%method = 1
      p%dt = Interval
      p%tolerance = 1e-6
      

   !p%DT  = Interval !bjj: this is set in ID_SetParameters


      !............................................................................................
      ! Define initial system states here:
      !............................................................................................

   z%DummyConstrState         = 0                                             ! we don't have constraint states


      ! initialize the continuous states:
   ! CALL Init_ContStates( x, p, InputFileData, OtherState, ErrStat2, ErrMsg2 )     ! initialize the continuous states
   !   CALL CheckError( ErrStat2, ErrMsg2 )
   !   IF (ErrStat >= AbortErrLev) RETURN

    IF (p%ModNo /= 6) THEN

       x%q              = 0                                             ! we don't have continuous states for ice model 1-5
       x%dqdt           = 0                                             

    ELSE

       x%q    = p%InitLoc    ! Initial ice floe location
       x%dqdt = p%v          ! Initial ice velocity

    ENDIF

      ! initialize the discrete states:
   !CALL ID_Init_DiscrtStates( xd, p, InputFileData, ErrStat2, ErrMsg2 )     ! initialize the continuous states
   !   CALL CheckError( ErrStat2, ErrMsg2 )
   !   IF (ErrStat >= AbortErrLev) RETURN
   
   xd%DummyDiscState      = 0

      ! Initialize other states:
    CALL ID_Init_OtherStates( OtherState, p, x, InputFileData, ErrStat2, ErrMsg2 )    ! initialize the other states
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN
      
      if ( p%method .eq. 2) then

         Allocate( OtherState%xdot(4), STAT=ErrStat )
         IF (ErrStat /= 0) THEN
            ErrStat = ErrID_Fatal
            ErrMsg = ' Error in IceDyn: could not allocate OtherStat%xdot.'
            RETURN
         END IF

      elseif ( p%method .eq. 3) then

         Allocate( OtherState%xdot(4), STAT=ErrStat )
         IF (ErrStat /= 0) THEN
            ErrStat = ErrID_Fatal
            ErrMsg = ' Error in IceDyn: could not allocate OtherStat%xdot.'
            RETURN
         END IF

      endif

      !............................................................................................
      ! Define initial guess for the system inputs and output (set up meshes) here:
      !............................................................................................


    ! Define initial guess for the system inputs here:

      u%q    = 0.
      u%dqdt = 0.
      
      y%fice = 0.
      
      ! Define system output initializations (set up mesh) here:
       CALL MeshCreate( BlankMesh       = u%PointMesh            &
                      ,IOS             = COMPONENT_INPUT        &
                      ,NNodes          = 1                      &
                      ,TranslationDisp = .TRUE.                 &
                      ,TranslationVel  = .TRUE.                 &
                      ,TranslationAcc  = .TRUE.                 &
                      ,nScalars        = 0                      &
                      ,ErrStat         = ErrStat                &
                      ,ErrMess         = ErrMsg                 )

       CALL MeshConstructElement ( Mesh = u%PointMesh            &
                                , Xelement = ELEMENT_POINT      &
                                , P1       = 1                  &
                                , ErrStat  = ErrStat            &
                                , ErrMess  = ErrMsg             )

      ! place single node at origin; position affects mapping/coupling with other modules
      TmpPos(1) = 0.
      TmpPos(2) = 0.
      TmpPos(3) = 0.

      CALL MeshPositionNode ( Mesh  = u%PointMesh   &
                            , INode = 1             &
                            , Pos   = TmpPos        &
                            , ErrStat   = ErrStat   &
                            , ErrMess   = ErrMsg    )

      CALL MeshCommit ( Mesh    = u%PointMesh     &
                       ,ErrStat = ErrStat         &
                       ,ErrMess = ErrMsg          )

      CALL MeshCopy ( SrcMesh  = u%PointMesh      &
                    , DestMesh = y%PointMesh      &
                    , CtrlCode = MESH_SIBLING     &
                    , Force    = .TRUE.           &
                    , ErrStat  = ErrStat          &
                    , ErrMess  = ErrMsg           )

      ! Define initial guess for the system inputs here:

      y%PointMesh%Force(1,1)   = 0.
      y%PointMesh%Force(2,1)   = 0.
      y%PointMesh%Force(3,1)   = 0.

      u%PointMesh%TranslationDisp(1,1) = 0.
      u%PointMesh%TranslationDisp(2,1) = 0.
      u%PointMesh%TranslationDisp(3,1) = 0.

      u%PointMesh%TranslationVel(1,1) = 0.
      u%PointMesh%TranslationVel(2,1) = 0.
      u%PointMesh%TranslationVel(3,1) = 0.

      u%PointMesh%TranslationAcc(1,1) = 0.
      u%PointMesh%TranslationAcc(2,1) = 0.
      u%PointMesh%TranslationAcc(3,1) = 0.

      ! set remap flags to true
      y%PointMesh%RemapFlag = .True.
      u%PointMesh%RemapFlag = .True.



      !............................................................................................
      ! Define initialization-routine output here:
      !............................................................................................
   CALL AllocAry( InitOut%WriteOutputHdr, p%NumOuts, 'WriteOutputHdr', ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN
   CALL AllocAry( InitOut%WriteOutputUnt, p%NumOuts, 'WriteOutputUnt', ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN

   InitOut%WriteOutputHdr = p%OutName
   InitOut%WriteOutputUnt = p%OutUnit

      !............................................................................................
      ! If you want to choose your own rate instead of using what the glue code suggests, tell the glue code the rate at which
      !   this module must be called here:
      !............................................................................................

  ! Interval = p%DT


  !     ! Print the summary file if requested:
  ! IF (InputFileData%SumPrint) THEN
  !    CALL ID_PrintSum( p, OtherState, GetAdamsVals, ErrStat2, ErrMsg2 )
  !       CALL CheckError( ErrStat2, ErrMsg2 )
  !       IF (ErrStat >= AbortErrLev) RETURN
  ! END IF

       ! Destroy the InputFileData structure (deallocate arrays)

   CALL ID_DestroyInputFile(InputFileData, ErrStat2, ErrMsg2 )
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
      CHARACTER(1024)            :: ErrMsg3     ! The error message (ErrMsg)

      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF ( LEN_TRIM(ErrMsg) > 0 ) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//' '//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL ID_DestroyInputFile(InputFileData, ErrStat3, ErrMsg3 )
         END IF

      END IF


   END SUBROUTINE CheckError

END SUBROUTINE ID_Init
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
!
! This routine is called at the end of the simulation.
!..................................................................................................................................

      TYPE(ID_InputType),           INTENT(INOUT)  :: u           ! System inputs
      TYPE(ID_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
      TYPE(ID_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states
      TYPE(ID_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Discrete states
      TYPE(ID_ConstraintStateType), INTENT(INOUT)  :: z           ! Constraint states
      TYPE(ID_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(ID_OutputType),          INTENT(INOUT)  :: y           ! System outputs
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""

      ! Place any last minute operations or calculations here:

      ! Close files here:

      ! Destroy the input data:

      CALL ID_DestroyInput( u, ErrStat, ErrMsg )

      ! Destroy the parameter data:

      CALL ID_DestroyParam( p, ErrStat, ErrMsg )

      ! Destroy the state data:

      CALL ID_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL ID_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL ID_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL ID_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )

      ! Destroy the output data:

      CALL ID_DestroyOutput( y, ErrStat, ErrMsg )


END SUBROUTINE ID_End
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! Routine for solving for constraint states, integrating continuous states, and updating discrete states
! Constraint states are solved for input t; Continuous and discrete states are updated for t + p%dt
! (stepsize dt assumed to be in ModName parameter)
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   ) :: t          ! Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   ) :: n          ! Current simulation time step n = 0,1,...
      TYPE(ID_InputType),                 INTENT(INOUT) :: u(:)       ! Inputs at utimes
      REAL(DbKi),                         INTENT(IN   ) :: utimes(:)  ! Times associated with u(:), in seconds
      TYPE(ID_ParameterType),             INTENT(IN   ) :: p          ! Parameters
      TYPE(ID_ContinuousStateType),       INTENT(INOUT) :: x          ! Input: Continuous states at t;
                                                                      !   Output: Continuous states at t + Interval
      TYPE(ID_DiscreteStateType),         INTENT(INOUT) :: xd         ! Input: Discrete states at t;
                                                                      !   Output: Discrete states at t  + Interval
      TYPE(ID_ConstraintStateType),       INTENT(INOUT) :: z          ! Input: Initial guess of constraint states at t+dt;
                                                                      !   Output: Constraint states at t+dt
      TYPE(ID_OtherStateType),            INTENT(INOUT) :: OtherState ! Other/optimization states
      INTEGER(IntKi),                     INTENT(  OUT) :: ErrStat    ! Error status of the operation
      CHARACTER(*),                       INTENT(  OUT) :: ErrMsg     ! Error message if ErrStat /= ErrID_None

      ! local variables

      !TYPE(ID_InputType)                               :: u_interp      ! input interpolated from given u at utimes
      !TYPE(ID_ContinuousStateType)                     :: xdot          ! continuous state time derivative
      INTEGER(IntKi)                                    :: I             ! Loop count
      REAL(ReKi)                                        :: Del2          ! Deflection of the current ice tooth, for model 2,3
!      REAL(ReKi)                                       :: Del(p%Zn)     ! Deflection of ice tooth in each zone, for model 4
      REAL(ReKi)                                        :: StrRt         ! Strain rate (s^-1)
      REAL(ReKi)                                        :: SigCrp        ! Creep stress (Pa)
      
      INTEGER(IntKi)                                    :: nt            ! Current time step
      INTEGER(IntKi)                                    :: Nc            ! Current ice tooth number
      REAL(ReKi)                                        :: Psumc         ! Current sum of pitch of broken ice teeth (m)
      REAL(ReKi)                                        :: Dmc           ! Delmax of the current tooth (m)
      REAL(ReKi)                                        :: Pc            ! Pitch of the current tooth (m)
      REAL(ReKi)                                        :: Pnxt          ! Pitch of the next ice tooth (m)
       ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""
      
      nt = INT( t / p%dt ) + 1  !bjj: should be the same as n+1, no?

      ! Update Other states here
      
      IF (p%ModNo == 2) THEN

        Del2 = p%InitLoc + p%v * t - p%Pitch * (OtherState%IceTthNo2 -1) - u(1)%q

            IF ( Del2 >= (p%Delmax2 - p%tolerance) ) THEN

                OtherState%IceTthNo2 = OtherState%IceTthNo2 + 1
                !CALL WrScr( 'IceTthNo=' // Num2LStr(xd%IceTthNo2) )

            ENDIF

      ENDIF
      
      ! Update Other states here
     
      IF (p%ModNo == 3) THEN
         Nc = OtherState%Nc (nt)
         Psumc = OtherState%Psum (nt)
         Pc = p%rdmP (Nc)
      
         IF (p%SubModNo == 3) THEN
      
            Del2 = p%InitLoc + p%v * t - OtherState%Psum (nt) - u(1)%q ! Deflection of the current ice tooth
            Dmc  = p%RdmDm ( OtherState%Nc(nt) )
            
            IF ( Del2 >= ( p%RdmDm ( OtherState%Nc(nt) ) - p%tolerance) ) THEN ! Current ice tooth breaks
                
                DO I= nt+1 , p%TmStep
                   OtherState%Nc (I) = Nc + 1
                   OtherState%Psum (I) = Psumc + Pc
                   Pnxt = OtherState%Psum (I)
                ENDDO
         
            ENDIF 
            
         ENDIF
          
      ENDIF
      

      ! Update Continuous states here
      
      IF (p%ModNo == 6) THEN

           if (p%method .eq. 1) then

               CALL ID_RK4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )

            elseif (p%method .eq. 2) then

               CALL ID_AB4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )

            elseif (p%method .eq. 3) then

               CALL ID_ABM4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )

            else

               ErrStat = ErrID_Fatal
               ErrMsg  = ' Error in ID_UpdateStates: p%method must be 1 (RK4), 2 (AB4), or 3 (ABM4)'
               RETURN

            endif

      END IF
           

      IF ( ErrStat >= AbortErrLev ) RETURN

END SUBROUTINE ID_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_CalcOutput( t, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
!
! Routine for computing outputs, used in both loose and tight coupling.
!..................................................................................................................................

      REAL(DbKi),                    INTENT(IN   )  :: t           ! Current simulation time in seconds
      TYPE(ID_InputType),            INTENT(IN   )  :: u           ! Inputs at t
      TYPE(ID_ParameterType),        INTENT(IN   )  :: p           ! Parameters
      TYPE(ID_ContinuousStateType),  INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(ID_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(ID_ConstraintStateType),  INTENT(IN   )  :: z           ! Constraint states at t
      TYPE(ID_OtherStateType),       INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(ID_OutputType),           INTENT(INOUT)  :: y           ! Outputs computed at t (Input only so that mesh con-
                                                                    !   nectivity information does not have to be recalculated)
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
      INTEGER(IntKi)                    :: I                        ! Loop count
      INTEGER(IntKi)                    :: nt                       ! Current time step
      REAL(ReKi)                        :: Del2                     ! Deflection of the current ice tooth, for model 2
!      REAL(ReKi)                        :: Del(p%Zn)                ! Deflection of each ice tooth, for model 4
!      REAL(ReKi)                        :: ZoneF(p%Zn)              ! Ice force of each ice tooth, for model 4
      REAL(ReKi)                        :: R                        ! Cone radius
      REAL(ReKi)                        :: Pn1                      ! Normal force from cone to ice piece 1, for model 5
      REAL(ReKi)                        :: Pn2                      ! Normal force from cone to ice piece 2, for model 5
      REAL(ReKi)                        :: Xb                       ! Moment arm for buoyancy force, for model 5 
      REAL(ReKi)                        :: Fb                       ! Buoyancy force, for model5
      REAL(ReKi)                        :: pbeta
      REAL(ReKi)                        :: qbeta
      REAL(ReKi)                        :: Rh                       ! Horizontal force, for model5
      REAL(ReKi)                        :: Rv                       ! Vertical force, for model5
      REAL(ReKi)                        :: Wr                       ! Ride-up ice weight, for model5 
      REAL(ReKi)                        :: CrntT0

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""

      nt = INT( t / p%dt ) + 1
      
      ! Compute outputs here:

      ! Ice Model 1 --------------------------------
      IF (p%ModNo == 1) THEN

          IF (p%SubModNo == 1) THEN

              IF (t < p%t0 + p%tolerance) THEN

                  y%fice = 0

              ELSEIF (t< p%t0 + p%tm1a + p%tolerance) THEN

                  y%fice = (t-p%t0) / p%tm1a * p%Fmax1a

              ELSE
                  
                  y%fice = p%Fmax1a

              ENDIF

          ELSEIF (p%SubModNo == 2) THEN

              IF (t < p%t0 + p%tolerance) THEN

                  y%fice = 0

              ELSEIF (t< p%t0 + p%tm1b + p%tolerance) THEN

                  y%fice = (t-p%t0) / p%tm1b * p%Fmax1b

              ELSE
                  
                  y%fice = 0

              ENDIF

          ELSEIF (p%SubModNo == 3) THEN

              IF (t < p%t0 + p%tolerance) THEN

                  y%fice = 0

              ELSEIF (t< p%t0 + p%tm1c + p%tolerance) THEN

                  y%fice = (t-p%t0) / p%tm1c * p%Fmax1c

              ELSE

                  y%fice = p%Fmax1c

              ENDIF

          ENDIF

      ENDIF

      ! Ice Model 2 --------------------------------

      IF (p%ModNo == 2) THEN

          Del2  = p%InitLoc + p%v * t - p%Pitch * (OtherState%IceTthNo2 -1) - u%q

            IF (p%Delmax2 <= p%Pitch) THEN !Sub-model 1 
            
                IF ( Del2 <= 0) THEN

                    y%fice = 0.0

                ELSEIF (Del2 < p%Delmax2 + p%tolerance) THEN

                    y%fice = Del2 * p%Kice2

                ELSE

                    ErrStat = ErrID_Fatal
                    ErrMsg  = ' Error in IceDyn Model 2: Ice tooth does not break when Del>Delmax'

                ENDIF
                
            ELSEIF (p%Delmax2 <= 2*p%Pitch) THEN !Sub-model 2, two ice teeth bending
            
                IF ( Del2 <= 0) THEN

                    y%fice = 0.0

                ELSEIF (Del2 < p%Pitch) THEN

                    y%fice = Del2 * p%Kice2

                ELSEIF (Del2 < p%Delmax2 + p%tolerance) THEN 
                
                    y%fice = Del2 * p%Kice2 + (Del2 - p%Pitch) * p%Kice2
                
                ELSE

                    ErrStat = ErrID_Fatal
                    ErrMsg  = ' Error in IceDyn Model 2: Ice tooth does not break when Del>Delmax'

                ENDIF
            
            ENDIF

      END IF
      
      IF (p%ModNo == 3) THEN

          IF (p%SubModNo == 1) THEN

              IF ( t <= p%rdmt0 (nt) ) THEN

                  y%fice = 0

              ELSEIF (t <= p%rdmt0 (nt) + p%rdmtm (nt)) THEN

                  y%fice = (t-p%rdmt0 (nt)) /  p%rdmtm(nt) * p%rdmFm (nt)

              ELSE

                  y%fice = p%rdmFm (nt)

              ENDIF

          ELSEIF (p%SubModNo == 2) THEN
              
              CrntT0 = p%rdmt0 (nt)

              IF (t <= p%rdmt0 (nt)) THEN

                  y%fice = 0

              ELSEIF (t <= p%rdmt0 (nt) + p%rdmtm (nt)) THEN

                  y%fice = (t-p%rdmt0 (nt)) / p%rdmtm(nt) * p%rdmFm (nt)

              ELSE

                  y%fice = 0

              ENDIF

          ELSEIF (p%SubModNo == 3) THEN

              Del2  = p%InitLoc + p%v * t - OtherState%Psum (nt) - u%q     ! Determine the contact state between ice sheet and the tower

                !IF (Del2 >= xd%Dmaxn) THEN
                !    ErrStat = ErrID_Fatal
                !    ErrMsg  = ' Error in IceDyn Model 3c: two ice teeth break at once'
                !ENDIF

              IF (Del2 <= 0) THEN

                y%fice = 0

                ELSE IF (Del2 > 0 .AND. Del2 <= p%rdmP (OtherState%Nc(nt)) ) THEN

                   y%fice = p%rdmKi(OtherState%Nc(nt))  * Del2

                ELSE

                   y%fice = p%rdmKi(OtherState%Nc(nt)) * Del2 + p%rdmKi(OtherState%Nc(nt+1)) * (Del2-p%rdmP (OtherState%Nc(nt))) ! Two teeth in contact

                ENDIF

          ENDIF
          
      ENDIF
      
      y%PointMesh%Force(1,1) = y%fice
      y%PointMesh%Force(2,1) = 0.
      y%PointMesh%Force(3,1) = 0.
            
END SUBROUTINE ID_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, xdot, ErrStat, ErrMsg )
!
! Routine for computing derivatives of continuous states.
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      TYPE(ID_InputType),             INTENT(IN   )  :: u           ! Inputs at t
      TYPE(ID_ParameterType),         INTENT(IN   )  :: p           ! Parameters
      TYPE(ID_ContinuousStateType),   INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(ID_DiscreteStateType),     INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(ID_ConstraintStateType),   INTENT(IN   )  :: z           ! Constraint states at t
      TYPE(ID_OtherStateType),        INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(ID_ContinuousStateType),   INTENT(  OUT)  :: xdot        ! Continuous state derivatives at t
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
      ! REAL(ReKi) :: mass    ! Mass of ice feature (kg)
      REAL(ReKi) :: force     ! Ice force (N)
      REAL(ReKi) :: R       ! Structure radius
      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


      

END SUBROUTINE ID_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_UpdateDiscState( t, n, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! Routine for updating discrete states
!..................................................................................................................................

      REAL(DbKi),                    INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                INTENT(IN   )  :: n           ! Current step of the simulation: t = n*Interval
      TYPE(ID_InputType),            INTENT(IN   )  :: u           ! Inputs at t
      TYPE(ID_ParameterType),        INTENT(IN   )  :: p           ! Parameters
      TYPE(ID_ContinuousStateType),  INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(ID_DiscreteStateType),    INTENT(INOUT)  :: xd          ! Input: Discrete states at t;
                                                                    !   Output: Discrete states at t + Interval
      TYPE(ID_ConstraintStateType),  INTENT(IN   )  :: z           ! Constraint states at t
      TYPE(ID_OtherStateType),       INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                  INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""

      ! Update discrete states here:
            xd%DummyDiscState = 0.0

END SUBROUTINE ID_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_CalcConstrStateResidual( t, u, p, x, xd, z, OtherState, Z_residual, ErrStat, ErrMsg )
!
! Routine for solving for the residual of the constraint state functions
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      TYPE(ID_InputType),            INTENT(IN   )  :: u           ! Inputs at t
      TYPE(ID_ParameterType),        INTENT(IN   )  :: p           ! Parameters
      TYPE(ID_ContinuousStateType),  INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(ID_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(ID_ConstraintStateType),  INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
      TYPE(ID_OtherStateType),       INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(ID_ConstraintStateType),  INTENT(  OUT)  :: Z_residual  ! Residual of the constraint state functions using
                                                                    !     the input values described above
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


      ! Solve for the constraint states here:

      Z_residual%DummyConstrState = 0

END SUBROUTINE ID_CalcConstrStateResidual
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_ReadInput( InitInp, InputFileData, ErrStat, ErrMsg )
!     This public subroutine reads the input required for IceDyn from the file whose name is an
!     input parameter.
!----------------------------------------------------------------------------------------------------


      ! Passed variables

   TYPE(ID_InitInputType),        INTENT( IN    )   :: InitInp              ! the IceDyn initial input data
   TYPE(ID_InputFile),            INTENT(   OUT )   :: InputFileData        ! Data stored in the IceDyn's input file
   INTEGER,                       INTENT(   OUT )   :: ErrStat              ! returns a non-zero value when an error occurs
   CHARACTER(*),                  INTENT(   OUT )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None


      ! Local variables

   INTEGER                                          :: UnIn                 ! Unit number for the input file
   CHARACTER(1024)                                  :: FileName             ! Name of HydroDyn input file



      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""


   !-------------------------------------------------------------------------------------------------
   ! Open the file
   !-------------------------------------------------------------------------------------------------
   FileName = TRIM(InitInp%InputFile)

   CALL GetNewUnit( UnIn )
   CALL OpenFInpFile( UnIn, FileName, ErrStat )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to open IceDyn input file: '//FileName
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF


   !CALL WrScr( 'Opening HydroDyn input file:  '//FileName )


   !-------------------------------------------------------------------------------------------------
   ! File header
   !-------------------------------------------------------------------------------------------------

   CALL ReadCom( UnIn, FileName, 'IceDyn input file header line 1', ErrStat )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read HydroDyn input file header line 1.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF


   CALL ReadCom( UnIn, FileName, 'IceDyn input file header line 2', ErrStat )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read HydroDyn input file header line 2.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   CALL ReadCom( UnIn, FileName, 'IceDyn input file header line 3', ErrStat )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read HydroDyn input file header line 3.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Ice Model Number
   !-------------------------------------------------------------------------------------------------

      ! Header

   CALL ReadCom( UnIn, FileName, 'Ice Models header', ErrStat, ErrMsg)

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Comment line - Ice Models'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF


      ! IceModel - Number that represents different ice models.

   CALL ReadVar ( UnIn, FileName, InputFileData%IceModel, 'IceModel', 'Number that represents different ice models', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read IceModel parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF


    ! IceSubModel - Number that represents different ice sub models.

   CALL ReadVar ( UnIn, FileName, InputFileData%IceSubModel, 'IceSubModel', 'Number that represents different ice sub-models', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read IceSubModel parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Ice properties - General
   !-------------------------------------------------------------------------------------------------

      ! Header

   CALL ReadCom( UnIn, FileName, 'Ice properties - General header', ErrStat, ErrMsg)

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Ice properties - General header comment line.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

     ! IceVel - Velocity of ice sheet movement (m/s)

   CALL ReadVar ( UnIn, FileName, InputFileData%v, 'IceVel', 'Velocity of ice sheet movement', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read IceVel parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   ! IceThks  - Thickness of the ice sheet (m)

   CALL ReadVar ( UnIn, FileName, InputFileData%h, 'IceThks', 'Thickness of the ice sheet (m)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read IceThks parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   ! StWidth  - Width of the structure in contact with the ice (m)

   CALL ReadVar ( UnIn, FileName, InputFileData%StrWd, 'StWidth', 'Width of the structure in contact with the ice (m)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read StWidth parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

      ! WtDen     - Mass density of water (kg/m3)

   CALL ReadVar ( UnIn, FileName, InputFileData%rhow, 'WtDen', 'Mass density of water (kg/m3)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read WtDen parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

      ! IceDen  - Mass density of ice (kg/m3)

   CALL ReadVar ( UnIn, FileName, InputFileData%rhoi, 'IceDen', 'Mass density of ice (kg/m3)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read IceDen parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

          ! InitLoc - Ice sheet initial location (m)

   CALL ReadVar ( UnIn, FileName, InputFileData%InitLoc, 'InitLoc', 'Ice sheet initial location (m)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read InitLoc parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

          ! InitTm - Ice load starting time (s)

   CALL ReadVar ( UnIn, FileName, InputFileData%t0, 'InitTm', 'Ice load starting time (s)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read InitTm parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
          ! Seed1 - Random seed 1

   CALL ReadVar ( UnIn, FileName, InputFileData%Seed1, 'Seed1', 'Random seed 1', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Seed1 parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
           ! Seed2 - Random seed 2

   CALL ReadVar ( UnIn, FileName, InputFileData%Seed2, 'Seed2', 'Random seed 2', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Seed2 parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Ice properties - Ice Model 1
   !-------------------------------------------------------------------------------------------------

        ! Header

   CALL ReadCom( UnIn, FileName, 'Ice properties - Ice Model 1 SubModel 1', ErrStat, ErrMsg)

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Ice properties - Ice Model 1 SubModel 1 header comment line.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

    ! Ikm - Indentation factor

   CALL ReadVar ( UnIn, FileName, InputFileData%Ikm, 'Ikm', 'Indentation factor', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Ikm parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

    ! Ag - Constant depends only on ice crystal type, used in calculating uniaxial stress (MPa-3s-1)

   CALL ReadVar ( UnIn, FileName, InputFileData%Ag, 'Ag', 'Constant depends only on ice crystal type, used in calculating uniaxial stress (MPa-3s-1)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Ag parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

    ! Qg - Activation Energy (kJmol^-1)

   CALL ReadVar ( UnIn, FileName, InputFileData%Qg, 'Qg', 'Activation Energy (kJmol^-1)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Qg parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

    ! Rg - Universal gas constant (Jmol-1K-1)

   CALL ReadVar ( UnIn, FileName, InputFileData%Rg, 'Rg', 'Universal gas constant (Jmol-1K-1)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Rg parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

    ! Tice - Ice temperature (K)

   CALL ReadVar ( UnIn, FileName, InputFileData%Tice, 'Tice', 'Ice temperature (K)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Tice parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF


    ! Header

   CALL ReadCom( UnIn, FileName, 'Ice properties - Ice Model 1 SubModel 2', ErrStat, ErrMsg)

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Ice properties - Ice Model 1 SubModel 2 header comment line.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

    ! Poisson  - Poisson's ratio of ice

   CALL ReadVar ( UnIn, FileName, InputFileData%nu, 'Poisson', 'Poisson ratio of ice', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Poisson parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

    ! WgAngle - Wedge Angel, degree. Default 45 Degrees.

   CALL ReadVar ( UnIn, FileName, InputFileData%phi, 'WgAngle', ' Wedge Angel, degree. Default 45 Degrees', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read WgAngle parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

    ! EIce - Young's modulus of ice (GPa)

   CALL ReadVar ( UnIn, FileName, InputFileData%Eice, 'EIce', 'Youngs modulus of ice (GPa)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read EIce parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

    ! Header

   CALL ReadCom( UnIn, FileName, 'Ice properties - Ice Model 1 SubModel 3', ErrStat, ErrMsg)

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Ice properties - Ice Model 1 SubModel 3 header comment line.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

    ! SigN - Nominal ice stress (MPa)

   CALL ReadVar ( UnIn, FileName, InputFileData%SigNm, 'SigNm', 'Nominal ice stress (MPa)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read SigN parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Ice properties - Ice Model 2
   !-------------------------------------------------------------------------------------------------

   ! Header

   CALL ReadCom( UnIn, FileName, 'Ice properties - Ice Model 2 SubModel 1,2', ErrStat, ErrMsg)

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Ice properties - Ice Model 2 SubModel 1,2 header comment line.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

    ! Pitch - Distance between sequential ice teeth (m)

   CALL ReadVar ( UnIn, FileName, InputFileData%Pitch, 'Pitch', 'Distance between sequential ice teeth (m)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Pitch parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

    ! IceStr2 - Ice failure stress (MPa)

   CALL ReadVar ( UnIn, FileName, InputFileData%IceStr2, 'IceStr2', 'Ice failure stress (MPa)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read IceStr2 parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

    ! Delmax2 - Ice tooth maximum elastic deformation (m)

   CALL ReadVar ( UnIn, FileName, InputFileData%Delmax2, 'Delmax2', 'Ice tooth maximum elastic deformation (m)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Delmax2 parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Ice properties - Ice Model 3
   !-------------------------------------------------------------------------------------------------

     ! Header

   CALL ReadCom( UnIn, FileName, 'Ice PROPERTIES -Ice Model 3, SubModel 1,2', ErrStat, ErrMsg)

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Ice PROPERTIES -Ice Model 3, SubModel 1,2 header comment line.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   ! ThkMean - Mean value of ice thickness (m)

   CALL ReadVar ( UnIn, FileName, InputFileData%miuh, 'ThkMean', 'Mean value of ice thickness (m)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read ThkMean parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   ! ThkVar - Variance of ice thickness (m^2)

   CALL ReadVar ( UnIn, FileName, InputFileData%varh, 'ThkVar', 'Variance of ice thickness (m^2)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read ThkVar parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   ! VelMean - Mean value of ice velocity (m/s)

   CALL ReadVar ( UnIn, FileName, InputFileData%miuv, 'VelMean', 'Mean value of ice velocity (m/s)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read VelMean parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   ! VelVar - Variance of ice velocity (m^2/s^2)

   CALL ReadVar ( UnIn, FileName, InputFileData%varv, 'VelVar', 'Variance of ice velocity (m^2/s^2)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read VelVar parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   ! TeMean - Mean value of ice loading event duration (s)

   CALL ReadVar ( UnIn, FileName, InputFileData%miut, 'TeMean', 'Mean value of ice loading event duration (s)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read TeMean parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

    ! Header

   CALL ReadCom( UnIn, FileName, 'Ice PROPERTIES -Ice Model 3, SubModel 2,3', ErrStat, ErrMsg)

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Ice PROPERTIES -Ice Model 3, SubModel 2,3 header comment line.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   ! StrMean - Mean value of ice thickness (m)

   CALL ReadVar ( UnIn, FileName, InputFileData%miubr, 'StrMean', 'Mean value of ice strength (MPa)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read StrMean parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   ! StrVar - Variance of ice thickness (m^2)

   CALL ReadVar ( UnIn, FileName, InputFileData%varbr, 'StrVar', 'Variance of ice strength (MPa)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read StrVar parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   ! Header

   CALL ReadCom( UnIn, FileName, 'Ice PROPERTIES -Ice Model 3, SubModel 3', ErrStat, ErrMsg)

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Ice PROPERTIES -Ice Model 3, SubModel 3 header comment line.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   ! DelMean - Mean value of maximum ice tooth tip displacement (m)

   CALL ReadVar ( UnIn, FileName, InputFileData%miuDelm, 'DelMean', 'Mean value of maximum ice tooth tip displacement (m)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read DelMean parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   ! DelVar - Variance of maximum ice tooth tip displacement (m^2)

   CALL ReadVar ( UnIn, FileName, InputFileData%varDelm, 'DelVar', 'Variance of maximum ice tooth tip displacement (m^2)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read DelVar parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   ! PMean - Mean value of the distance between sequential ice teeth (m)

   CALL ReadVar ( UnIn, FileName, InputFileData%miuP, 'PMean', 'Mean value of the distance between sequential ice teeth (m)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read PMean parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   ! PVar - Variance of the distance between sequential ice teeth (m^2)

   CALL ReadVar ( UnIn, FileName, InputFileData%varP, 'PVar', 'Variance of the distance between sequential ice teeth (m^2)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read PVar parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF


   !-------------------------------------------------------------------------------------------------
   ! This is the end of the input file
   !-------------------------------------------------------------------------------------------------
   CLOSE ( UnIn )


   RETURN


END SUBROUTINE ID_ReadInput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_SetParameters( InputFileData, p, Interval, Tmax, ErrStat, ErrMsg  )
! This takes the primary input file data and sets the corresponding parameters.
!..................................................................................................................................

   IMPLICIT                        NONE


      ! Passed variables

   TYPE(ID_ParameterType),   INTENT(INOUT)  :: p                            ! Parameters of the IceDyn module
   TYPE(ID_InputFile),       INTENT(IN)     :: InputFileData                ! Data stored in the module's input file
   REAL(DbKi),               INTENT(IN)     :: Interval                     ! Coupling interval in seconds
   REAL(DbKi),               INTENT(IN   )  :: Tmax                         ! Total simulation time 
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                      ! Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                       ! Error message

     ! Local variables
   INTEGER (IntKi)                          :: I
     ! Ice Model 1

   REAL(ReKi)                               :: StrRt                        ! Strain rate (s^-1)
   REAL(ReKi)                               :: SigCrp                       ! Creep stress (Pa)
   REAL(ReKi)                               :: Bf                           ! Flexural rigidity of ice plate (Nm)
   REAL(ReKi)                               :: kappa                        ! Constants in Ice Model 1b
   REAL(ReKi)                               :: PhiR                         ! Phi in radius
   
     ! Ice Model 3
     
   REAL(ReKi), allocatable                  :: rdmh(:)                      ! Random ice thickness time series (m)
   REAL(ReKi), allocatable                  :: rdmv(:)                      ! Random ice velocity time series (m/s)
   REAL(ReKi), allocatable                  :: rdmte(:)                     ! Random ice loading event time (s)
   REAL(ReKi), allocatable                  :: rdmsig(:)                    ! Random ice strength time series (Pa)
   REAL(ReKi), allocatable                  :: rdmTstr(:)                   ! Random ice strength of ice teeth (Pa)
   REAL(ReKi)                               :: t                            ! Time for generating random parameters time series (s)
   REAL(ReKi)                               :: rdmScrp                      ! Random ice creeping strength (MPa)
   INTEGER(IntKi)                           :: Nthmax                       ! Approximate maximum ice teeth number
   INTEGER(IntKi)                           :: nSeeds                       ! Number of seeds needed to generate random number
   INTEGER(IntKi), allocatable              :: TmpIceSeeds(:)               ! Random seeds needed for initialization random number generater for IceDyn 
   INTEGER(IntKi)                           :: J                            ! Loop count 
   REAL(ReKi)                               :: h_dmmy                       ! Dummy random number for h (m)
   REAL(ReKi)                               :: v_dmmy                       ! Dummy random number for v (m/s)
   REAL(ReKi)                               :: t_dmmy                       ! Dummy random number for t (s)
   REAL(ReKi)                               :: s_dmmy                       ! Dummy random number for sig (Pa)
   REAL(ReKi)                               :: Dm_dmmy                      ! Dummy random number for Dmax (m)
   REAL(ReKi)                               :: P_dmmy                       ! Dummy random number for P (m)
   REAL(ReKi)                               :: CrntT0

    ! Initialize error data
   ErrStat = ErrID_None
   ErrMsg  = ''

   !...............................................................................................................................
   ! Direct copy of variables:
   !...............................................................................................................................
   p%ModNo     = InputFileData%IceModel
   p%SubModNo  = InputFileData%IceSubModel
   p%v         = InputFileData%v
   p%h         = InputFileData%h
   p%InitLoc   = InputFileData%InitLoc
   p%t0        = InputFileData%t0
   p%StrWd     = InputFileData%StrWd
   p%dt        = Interval
   p%Tmax      = Tmax   

    ! Ice Model 1
   p%Ikm       = InputFileData%Ikm
   
    ! Ice Model 2
   p%Delmax2   = InputFileData%Delmax2
   p%Pitch     = InputFileData%Pitch
   
   IF (p%Delmax2 >= 2*p%Pitch) THEN
   
       ErrStat = ErrID_Fatal
       ErrMsg  = ' Error in IceDyn Model 2: Input Delmax2 should not be larger than 2*Pitch'
       
   ENDIF
      
   !...............................................................................................................................
   ! Calculate some indirect inputs:
   !...............................................................................................................................

   StrRt       = p%v / 4.0 / p%StrWd
   p%EiPa      = InputFileData%EIce * 1.0e9

   ! Ice Model 1a -------------------------------------------------------------------------

   SigCrp      = ( 1.0/InputFileData%Ag * exp( InputFileData%Qg / InputFileData%Rg / InputFileData%Tice ) * StrRt )**(1.0/3.0) * 1e6
   p%Cstr      = ( 1.0/InputFileData%Ag * exp( InputFileData%Qg / InputFileData%Rg / InputFileData%Tice ) )**(1.0/3.0) * 1e6
   p%Fmax1a    = InputFileData%Ikm * p%StrWd * p%h * SigCrp
   p%tm1a      = InputFileData%Ikm * SigCrp / p%EiPa / StrRt

   ! Ice Model 1b

   Bf          = p%EiPa * p%h**3 / 12.0 / (1.0-InputFileData%nu**2.0)
   kappa       = ( InputFileData%rhow * 9.81 / 4.0 / Bf ) ** 0.25
   PhiR        = InputFileData%phi / 180.0 * 3.1415927
   p%Fmax1b    = 5.3 * Bf * kappa * ( kappa * p%StrWd + 2 * tan(PhiR/2.0) )
   p%tm1b      = p%Fmax1b / p%StrWd / p%h / p%EiPa / StrRt

   ! Ice Model 1c

   p%Fmax1c     = p%StrWd * p%h * InputFileData%SigNm *1e6
   p%tm1c       = InputFileData%SigNm *1e6 / p%EiPa / StrRt

   ! Ice Model 2 -------------------------------------------------------------------------

   p%Kice2      = InputFileData%IceStr2 *1e6 * p%StrWd * p%h / p%Delmax2

   ! Ice Model 3 -------------------------------------------------------------------------
   
   p%TmStep     = INT( p%Tmax / p%dt ) + 1
   Nthmax       = p%v * p%Tmax / InputFileData%miuP * 2
   
       ! Random number initialization
       
   CALL RANDOM_SEED ( SIZE = nSeeds )
      
   IF ( nSeeds /= 2 ) THEN
      CALL ProgWarn( ' The random number generator in use differs from the original code provided by NREL. This pRNG uses ' &
                                  //TRIM(Int2LStr(nSeeds))//' seeds instead of the 2 in the IceDyn input file.')
      ErrStat = ErrID_Warn
   END IF

   ALLOCATE ( TmpIceSeeds ( nSeeds ), STAT=ErrStat )
   IF (ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for TmpIceSeeds array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF   

      ! We'll just populate this with odd seeds = Seed(1) and even seeds = Seed(2)
      
   DO I = 1,nSeeds,2
      TmpIceSeeds(I) = InputFileData%Seed1
   END DO
   DO I = 2,nSeeds,2
      TmpIceSeeds(I) = InputFileData%Seed2
   END DO                   
                  
   CALL RANDOM_SEED ( PUT=TmpIceSeeds )
   DEALLOCATE(TmpIceSeeds, STAT=ErrStat)
   IF (ErrStat /= ErrID_None ) THEN
      CALL ProgWarn( ' Error deallocating space for TmpIceSeeds array.' )
      ErrStat = ErrID_Warn
   END IF 
   
       ! Submodel 1 and 2
       
   CALL AllocAry( p%rdmFm, p%TmStep, 'rdmFm', ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN
   
   CALL AllocAry( p%rdmt0, p%TmStep, 'rdmt0', ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN
   
   CALL AllocAry( p%rdmtm, p%TmStep, 'rdmtm', ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN
   
   CALL AllocAry( rdmh, p%TmStep, 'rdmh', ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN
   
   CALL AllocAry( rdmv, p%TmStep, 'rdmv', ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN
   
   CALL AllocAry( rdmte, p%TmStep, 'rdmte', ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN
   
   CALL AllocAry( rdmsig, p%TmStep, 'rdmsig', ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN
   
       ! Submodel 3
       
   CALL AllocAry( p%rdmDm, Nthmax, 'rdmDm', ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN
   
   CALL AllocAry( p%rdmP, Nthmax, 'rdmP', ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN
   
   CALL AllocAry( p%rdmKi, Nthmax, 'rdmKi', ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN
   
   !CALL AllocAry( p%rdmTstr, Nthmax, 'rdmTstr', ErrStat, ErrMsg )
   !IF ( ErrStat >= AbortErrLev ) RETURN
   

   
   IF (p%SubModNo == 1) THEN
       
       DO J = 1,p%TmStep
     
          t = J*Interval
          IF( J == 1 ) THEN
              p%rdmt0(J)    = 0
              CALL ID_Generate_RandomNum ( rdmh(J), rdmv(J), rdmte(J), s_dmmy, Dm_dmmy, P_dmmy, p, InputFileData, ErrStat, ErrMsg)
              rdmScrp       = p%Cstr * ( rdmv(J) / 4.0 / p%StrWd )**(1.0/3.0)
              p%rdmFm(J)    = InputFileData%Ikm * p%StrWd * rdmh(J) * rdmScrp
              p%rdmtm(J)    = InputFileData%Ikm * rdmScrp / p%EiPa / ( rdmv(J) / 4.0 / p%StrWd )               
          ELSEIF ( t < p%rdmt0 (J-1) + rdmte (J-1) ) THEN
              rdmh(J)       = rdmh(J-1)
              rdmv(J)       = rdmv(J-1)
              rdmte(J)      = rdmte(J-1)
              rdmsig(J)     = rdmsig(J-1)
              p%rdmFm(J)      = p%rdmFm(J-1)
              p%rdmt0(J)      = p%rdmt0(J-1)
              p%rdmtm(J)      = p%rdmtm(J-1)
          ELSE              
              p%rdmt0(J)    = p%rdmt0 (J-1) + rdmte (J-1)
              CALL ID_Generate_RandomNum ( rdmh(J), rdmv(J), rdmte(J), s_dmmy, Dm_dmmy, P_dmmy, p, InputFileData, ErrStat, ErrMsg)
              rdmScrp       = p%Cstr * ( rdmv(J) / 4.0 / p%StrWd )**(1.0/3.0)
              p%rdmFm(J)    = InputFileData%Ikm * p%StrWd * rdmh(J) * rdmScrp
              p%rdmtm(J)    = InputFileData%Ikm * rdmScrp / p%EiPa / ( rdmv(J) / 4.0 / p%StrWd )              
          ENDIF
          
       END DO
       
    ELSEIF (p%SubModNo == 2) THEN
    
       DO J = 1,p%TmStep
     
          t = J*Interval
          IF( J == 1 ) THEN
              p%rdmt0(J)       = 0
              CALL ID_Generate_RandomNum ( rdmh(J), rdmv(J), rdmte(J), rdmsig(J), Dm_dmmy, P_dmmy, p, InputFileData, ErrStat, ErrMsg)
              p%rdmFm(J)    = p%StrWd * rdmh(J) * rdmsig(J)
              p%rdmtm(J)    = rdmsig(J) / p%EiPa / ( rdmv(J) / 4.0 / p%StrWd )               
          ELSEIF ( t < p%rdmt0 (J-1) + rdmte (J-1) ) THEN
              rdmh(J)       = rdmh(J-1)
              rdmv(J)       = rdmv(J-1)
              rdmte(J)      = rdmte(J-1)
              rdmsig(J)     = rdmsig(J-1)
              p%rdmFm(J)    = p%rdmFm(J-1)
              p%rdmt0(J)    = p%rdmt0(J-1)
              p%rdmtm(J)    = p%rdmtm(J-1)
          ELSE              
              p%rdmt0(J)   = p%rdmt0 (J-1) + rdmte (J-1)
              CALL ID_Generate_RandomNum ( rdmh(J), rdmv(J), rdmte(J), rdmsig(J), Dm_dmmy, P_dmmy, p, InputFileData, ErrStat, ErrMsg)
              p%rdmFm(J)    = p%StrWd * rdmh(J) * rdmsig(J)
              p%rdmtm(J)    = rdmsig(J) / p%EiPa / ( rdmv(J) / 4.0 / p%StrWd )              
          ENDIF
              CrntT0 = p%rdmt0 (J) !bjj: doesn't appear to be used...
       END DO
       
    ELSEIF (p%SubModNo == 3) THEN
    
        DO J = 1, Nthmax        
            CALL ID_Generate_RandomNum ( h_dmmy, v_dmmy, t_dmmy, rdmsig(J), p%rdmDm(J), p%rdmP(J), p, InputFileData, ErrStat, ErrMsg)
            p%rdmKi(J) = rdmsig(J) * p%StrWd * p%h / p%rdmDm(J)
        END DO
        
    ENDIF

    ! Deallocate local variables 
   
   IF ( ALLOCATED(TmpIceSeeds) )THEN
      DEALLOCATE(TmpIceSeeds, STAT=ErrStat)
      IF (ErrStat /= ErrID_None ) THEN
         CALL ProgWarn( ' Error deallocating space for TmpIceSeeds array.' )
         ErrStat = ErrID_Warn
      END IF
   END IF
   
   
   IF ( ALLOCATED(TmpIceSeeds) )THEN
      DEALLOCATE(TmpIceSeeds, STAT=ErrStat)
      IF (ErrStat /= ErrID_None ) THEN
         CALL ProgWarn( ' Error deallocating space for TmpIceSeeds array.' )
         ErrStat = ErrID_Warn
      END IF
   END IF
   

   !...............................................................................................................................
   ! Calculate Output variables:
   !...............................................................................................................................

   p%NumOuts    = 3

      CALL AllocAry( p%OutName, p%NumOuts, 'OutName', ErrStat, ErrMsg )
      IF (ErrStat >= AbortErrLev) RETURN

      p%OutName (1) = 'Time'
      p%OutName (2) = 'IceDisp'
      p%OutName (3) = 'IceForce'

      CALL AllocAry( p%OutUnit, p%NumOuts, 'OutUnit', ErrStat, ErrMsg )
      IF (ErrStat >= AbortErrLev) RETURN

      p%OutUnit (1) = '(s)'
      p%OutUnit (2) = '(m)'
      p%OutUnit (3) = '(kN)'
      
      ! Test parameter assignments
      
END SUBROUTINE ID_SetParameters
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_Init_OtherStates( OtherState, p, x, InputFileData, ErrStat, ErrMsg  )
! This routine initializes the other states of the module.
! It assumes the parameters are set and that InputFileData contains initial conditions for the other states.
!..................................................................................................................................
   IMPLICIT                        NONE

   TYPE(ID_OtherStateType),      INTENT(  OUT)  :: OtherState        ! Initial other states
   TYPE(ID_ParameterType),       INTENT(IN)     :: p                 ! Parameters of the IceDyn module
   TYPE(ID_ContinuousStateType), INTENT(IN   )  :: x                 ! Initial continuous states
   TYPE(ID_InputFile),           INTENT(IN)     :: InputFileData     ! Data stored in the module's input file
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat           ! Error status
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg            ! Error message

      ! local variables
   INTEGER(IntKi)                               :: I                 ! loop counter
   REAL(ReKi)                                   :: StrRt             ! Strain rate (s^-1)
   REAL(ReKi)                                   :: SigCrp            ! Creep stress (Pa)
      ! Initialize error data
   ErrStat = ErrID_None
   ErrMsg  = ''

   IF  ( p%ModNo == 2 ) THEN

       OtherState%IceTthNo2 = 1 ! Initialize first ice tooth number

   ENDIF
   
   IF  ( p%ModNo == 3 ) THEN

       CALL AllocAry( OtherState%Nc, p%TmStep, 'OtherState%Nc', ErrStat, ErrMsg )
       IF ( ErrStat >= AbortErrLev ) RETURN
   
       CALL AllocAry( OtherState%Psum, p%TmStep, 'OtherState%Psum', ErrStat, ErrMsg )
       IF ( ErrStat >= AbortErrLev ) RETURN
       
       DO I = 1,p%TmStep
          OtherState%Nc (I)      = 1. ! Initialize first ice tooth number
          OtherState%Psum (I)    = 0. ! Initialize sum of pitches of broken ice teeth 
       ENDDO

   ENDIF

   
END SUBROUTINE ID_Init_OtherStates
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_Generate_RandomNum ( h, v, t, s, Dm, Pch, p, InputFileData, ErrStat, ErrMsg)
! This routine generate random numbers for the module.
! It assumes the parameters are set and that InputFileData contains input for generating random numbers.
!..................................................................................................................................
   IMPLICIT                        NONE

   REAL(ReKi),                   INTENT(OUT)    :: h                 ! Random ice thickness (m)
   REAL(ReKi),                   INTENT(OUT)    :: v                 ! Random ice velocity (m/s)
   REAL(ReKi),                   INTENT(OUT)    :: t                 ! Random ice loading event time (s)
   REAL(ReKi),                   INTENT(OUT)    :: s                 ! Random ice strength (Pa)
   REAL(ReKi),                   INTENT(OUT)    :: Dm                ! Random ice tooth maximum displacement (m)
   REAL(ReKi),                   INTENT(OUT)    :: Pch               ! Random ice tooth spacing (m)
   TYPE(ID_ParameterType),       INTENT(IN)     :: p                 ! Parameters of the IceDyn module
   TYPE(ID_InputFile),           INTENT(IN)     :: InputFileData     ! Data stored in the module's input file
   INTEGER(IntKi),               INTENT(OUT)    :: ErrStat           ! Error status
   CHARACTER(*),                 INTENT(OUT)    :: ErrMsg            ! Error message   

      ! local variables
   INTEGER(IntKi)                               :: I                 ! loop counter
   REAL(ReKi)                                   :: SigLogh           ! sigma_log(h), standard deviation of log(h)
   REAL(ReKi)                                   :: MiuLogh           ! miu_log(h), mean value of log(h)
   REAL(ReKi)                                   :: VelSig            ! parameter for a Rayleigh distribution
   REAL(ReKi)                                   :: TeLamb            ! parameter for a exponential distribution
   REAL, PARAMETER                              :: Pi = 3.1415927

      ! Initialize error data
   ErrStat = ErrID_None
   ErrMsg  = ''

      !Ice thickness has a lognormal distribution
   SigLogh = SQRT( LOG ( InputFileData%varh / InputFileData%miuh + 1) )
   MiuLogh = LOG ( InputFileData%miuh ) - 0.5 * SigLogh **2
   h       = EXP( random_normal() * SigLogh + MiuLogh )

      !Ice velocity has a Rayleigh distribution
   VelSig = InputFileData%miuv / SQRT(Pi/2)
   v      = random_rayleigh (VelSig)

      !Iceloading event time has a exponential distribution
   TeLamb = 1 / InputFileData%miut
   t      = random_exponential(TeLamb)

      !Ice strength has a Weibull distribution
   s      = random_weibull (InputFileData%miubr, InputFileData%varbr) * 1e6

      !Ice teeth Delmax and pitch have normal distributions
    
    Dm    = InputFileData%miuDelm + InputFileData%varDelm ** 0.5 * random_normal()
    Pch   = InputFileData%miuP + InputFileData%varP ** 0.5 * random_normal()


    CONTAINS

!Functions that generate random number with respect to certain distributions

         FUNCTION random_normal() RESULT(fn_val)

         ! Adapted from the following Fortran 77 code
         !      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
         !      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
         !      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

         !  The function random_normal() returns a normally distributed pseudo-random
         !  number with zero mean and unit variance.

         !  The algorithm uses the ratio of uniforms method of A.J. Kinderman
         !  and J.F. Monahan augmented with quadratic bounding curves.

         REAL :: fn_val

         !     Local variables
         REAL,parameter  :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
                            r1 = 0.27597, r2 = 0.27846 
         REAL     :: u, v, x, y, q

         !     Generate P = (u,v) uniform in rectangle enclosing acceptance region

         DO
           CALL RANDOM_NUMBER(u)
           CALL RANDOM_NUMBER(v)
           v = 1.7156 * (v - 0.5)

         !     Evaluate the quadratic form
           x = u - s
           y = ABS(v) - t
           q = x**2 + y*(a*y - b*x)

         !     Accept P if inside inner ellipse
           IF (q < r1) EXIT
         !     Reject P if outside outer ellipse
           IF (q > r2) CYCLE
         !     Reject P if outside acceptance region
           IF (v**2 < -4.0*LOG(u)*u**2) EXIT
         END DO

         !     Return ratio of P's coordinates as the normal deviate
         fn_val = v/u
         RETURN

      END FUNCTION random_normal


      FUNCTION random_exponential(Lambda) RESULT(fn_val)

         ! Adapted from Fortran 77 code from the book:
         !     Dagpunar, J. 'Principles of random variate generation'
         !     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

         ! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
         ! A NEGATIVE EXPONENTIAL DlSTRIBUTION WlTH DENSITY PROPORTIONAL
         ! TO EXP(-random_exponential), USING INVERSION.

         REAL  :: fn_val
         REAL  :: Lambda

         !     Local variable
         REAL  :: r

         DO
           CALL RANDOM_NUMBER(r)
           IF (r > 0.0) EXIT
         END DO

         fn_val = -LOG(1-r)/lambda
         RETURN

      END FUNCTION random_exponential


      FUNCTION random_rayleigh(Sigma) RESULT(fn_val)

         ! Adapted from Fortran 77 code from the book:
         !     Dagpunar, J. 'Principles of random variate generation'
         !     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

         ! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
         ! A NEGATIVE EXPONENTIAL DlSTRIBUTION WlTH DENSITY PROPORTIONAL
         ! TO EXP(-random_exponential), USING INVERSION.

         REAL  :: fn_val
         REAL  :: Sigma

         !     Local variable
         REAL  :: r

         DO
           CALL RANDOM_NUMBER(r)
           IF (r > 0.0) EXIT
         END DO

         fn_val = SQRT(-LOG(1-r)*2*sigma**2)
         RETURN

      END FUNCTION random_rayleigh


      FUNCTION random_weibull (mean,var) RESULT(fn_val)

         !Function generates a random variate in (0,infinity) from Weibull
         !distribution with input mean value and variance

         IMPLICIT NONE

         REAL :: mean
         REAL :: var
         REAL :: fn_val

         !Local variables
         REAL :: k
         REAL :: Lambda
         REAL :: u

         k = wbpar (mean, var)
         lambda = mean / NWTC_gamma(1+1/k)

         CALL RANDOM_NUMBER(u)
         fn_val = lambda * (-log(1-u)) ** (1/k)

      END FUNCTION random_weibull


      FUNCTION wbpar (mean, var) Result (k1)

         !Calculate Weibull distribution parameters due to mean value and variance of the data
         IMPLICIT NONE

         REAL :: mean
         REAL :: var
         REAL :: k 
         REAL :: k1, F1, dFdk
         REAL :: error 

         INTEGER :: I
         
         k = 10
         error = 1e-6

         DO i = 1,10000

            F1 = (NWTC_gamma(1+1/k))**2 / NWTC_gamma(1+2/k) - mean**2/(mean**2+var);

            IF (abs(F1) < error) EXIT

            dFdk = 2* (NWTC_gamma(1+1/k))**2 * (-1/k**2) / NWTC_gamma(1+2/k) * (digamma(1+1/k) -digamma(1+2/k));
            k = k - F1/dFdk;

            END DO

            !IF (abs(F1) >= error)THEN

            !    WrScr('Weibull parameters never found')

            !ENDIF


          k1 = k

        END FUNCTION wbpar

        FUNCTION digamma(z) RESULT(phy)

         !Calculate the value of digamma function of z
         REAL, INTENT(IN) :: z
         REAL             :: phy

         phy = log(z) - 1/2/z - 1/12/z**2 + 1/120/z**4 - 1/252/z**6 + 1/240/z**8 - 5/660/z**10;

      END FUNCTION digamma

END SUBROUTINE ID_Generate_RandomNum
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_RK4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! This subroutine implements the fourth-order Runge-Kutta Method (RK4) for numerically integrating ordinary differential equations:
!
!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x).
!   Define constants k1, k2, k3, and k4 as
!        k1 = dt * f(t        , x_t        )
!        k2 = dt * f(t + dt/2 , x_t + k1/2 )
!        k3 = dt * f(t + dt/2 , x_t + k2/2 ), and
!        k4 = dt * f(t + dt   , x_t + k3   ).
!   Then the continuous states at t = t + dt are
!        x_(t+dt) = x_t + k1/6 + k2/3 + k3/3 + k4/6 + O(dt^5)
!
! For details, see:
! Press, W. H.; Flannery, B. P.; Teukolsky, S. A.; and Vetterling, W. T. "Runge-Kutta Method" and "Adaptive Step Size Control for
!   Runge-Kutta." §16.1 and 16.2 in Numerical Recipes in FORTRAN: The Art of Scientific Computing, 2nd ed. Cambridge, England:
!   Cambridge University Press, pp. 704-716, 1992.
!
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           ! time step number
      TYPE(ID_InputType),             INTENT(INOUT)  :: u(:)        ! Inputs at t
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   ! times of input
      TYPE(ID_ParameterType),         INTENT(IN   )  :: p           ! Parameters
      TYPE(ID_ContinuousStateType),   INTENT(INOUT)  :: x           ! Continuous states at t on input at t + dt on output
      TYPE(ID_DiscreteStateType),     INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(ID_ConstraintStateType),   INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
      TYPE(ID_OtherStateType),        INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables

      TYPE(ID_ContinuousStateType)                 :: xdot        ! time derivatives of continuous states
      TYPE(ID_ContinuousStateType)                 :: k1          ! RK4 constant; see above
      TYPE(ID_ContinuousStateType)                 :: k2          ! RK4 constant; see above
      TYPE(ID_ContinuousStateType)                 :: k3          ! RK4 constant; see above
      TYPE(ID_ContinuousStateType)                 :: k4          ! RK4 constant; see above
      TYPE(ID_ContinuousStateType)                 :: x_tmp       ! Holds temporary modification to x
      TYPE(ID_InputType)                           :: u_interp    ! interpolated value of inputs

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


      ! interpolate u to find u_interp = u(t)
      CALL ID_Input_ExtrapInterp( u, utimes, u_interp, t, ErrStat, ErrMsg )

      ! find xdot at t
      CALL ID_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, xdot, ErrStat, ErrMsg )

      k1%q    = p%dt * xdot%q
      k1%dqdt = p%dt * xdot%dqdt

      x_tmp%q    = x%q    + 0.5 * k1%q
      x_tmp%dqdt = x%dqdt + 0.5 * k1%dqdt

      ! interpolate u to find u_interp = u(t + dt/2)
      CALL ID_Input_ExtrapInterp(u, utimes, u_interp, t+0.5*p%dt, ErrStat, ErrMsg)

      ! find xdot at t + dt/2
      CALL ID_CalcContStateDeriv( t + 0.5*p%dt, u_interp, p, x_tmp, xd, z, OtherState, xdot, ErrStat, ErrMsg )

      k2%q    = p%dt * xdot%q
      k2%dqdt = p%dt * xdot%dqdt

      x_tmp%q    = x%q    + 0.5 * k2%q
      x_tmp%dqdt = x%dqdt + 0.5 * k2%dqdt

      ! find xdot at t + dt/2
      CALL ID_CalcContStateDeriv( t + 0.5*p%dt, u_interp, p, x_tmp, xd, z, OtherState, xdot, ErrStat, ErrMsg )

      k3%q    = p%dt * xdot%q
      k3%dqdt = p%dt * xdot%dqdt

      x_tmp%q    = x%q    + k3%q
      x_tmp%dqdt = x%dqdt + k3%dqdt

      ! interpolate u to find u_interp = u(t + dt)
      CALL ID_Input_ExtrapInterp(u, utimes, u_interp, t + p%dt, ErrStat, ErrMsg)

      ! find xdot at t + dt
      CALL ID_CalcContStateDeriv( t + p%dt, u_interp, p, x_tmp, xd, z, OtherState, xdot, ErrStat, ErrMsg )

      k4%q    = p%dt * xdot%q
      k4%dqdt = p%dt * xdot%dqdt

      x%q    = x%q    +  ( k1%q    + 2. * k2%q    + 2. * k3%q    + k4%q    ) / 6.
      x%dqdt = x%dqdt +  ( k1%dqdt + 2. * k2%dqdt + 2. * k3%dqdt + k4%dqdt ) / 6.

END SUBROUTINE ID_RK4
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_AB4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! This subroutine implements the fourth-order Adams-Bashforth Method (RK4) for numerically integrating ordinary differential
! equations:
!
!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x).
!
!   x(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
!
!  See, e.g.,
!  http://en.wikipedia.org/wiki/Linear_multistep_method
!
!  or
!
!  K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
!
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           ! time step number
      TYPE(ID_InputType),             INTENT(INOUT)  :: u(:)        ! Inputs at t
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   ! times of input
      TYPE(ID_ParameterType),         INTENT(IN   )  :: p           ! Parameters
      TYPE(ID_ContinuousStateType),   INTENT(INOUT)  :: x           ! Continuous states at t on input at t + dt on output
      TYPE(ID_DiscreteStateType),     INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(ID_ConstraintStateType),   INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
      TYPE(ID_OtherStateType),        INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


      ! local variables
      TYPE(ID_ContinuousStateType) :: xdot       ! Continuous state derivs at t
      TYPE(ID_InputType)           :: u_interp


      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""

      ! need xdot at t
      CALL ID_Input_ExtrapInterp(u, utimes, u_interp, t, ErrStat, ErrMsg)
      CALL ID_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, xdot, ErrStat, ErrMsg )

      if (n .le. 2) then

         OtherState%n = n

         OtherState%xdot ( 3 - n ) = xdot

         CALL ID_RK4(t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )

      else

         if (OtherState%n .lt. n) then

            OtherState%n = n
            OtherState%xdot(4)    = OtherState%xdot(3)
            OtherState%xdot(3)    = OtherState%xdot(2)
            OtherState%xdot(2)    = OtherState%xdot(1)

         elseif (OtherState%n .gt. n) then

            ErrStat = ErrID_Fatal
            ErrMsg = ' Backing up in time is not supported with a multistep method '
            RETURN

         endif

         OtherState%xdot ( 1 )     = xdot  ! make sure this is most up to date

         x%q    = x%q    + (p%dt / 24.) * ( 55.*OtherState%xdot(1)%q - 59.*OtherState%xdot(2)%q    + 37.*OtherState%xdot(3)%q  &
                                       - 9. * OtherState%xdot(4)%q )

         x%dqdt = x%dqdt + (p%dt / 24.) * ( 55.*OtherState%xdot(1)%dqdt - 59.*OtherState%xdot(2)%dqdt  &
                                          + 37.*OtherState%xdot(3)%dqdt  - 9.*OtherState%xdot(4)%dqdt )

      endif

END SUBROUTINE ID_AB4
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_ABM4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! This subroutine implements the fourth-order Adams-Bashforth-Moulton Method (RK4) for numerically integrating ordinary
! differential equations:
!
!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x).
!
!   Adams-Bashforth Predictor:
!   x^p(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
!
!   Adams-Moulton Corrector:
!   x(t+dt) = x(t)  + (dt / 24.) * ( 9.*f(t+dt,x^p) + 19.*f(t,x) - 5.*f(t-dt,x) + 1.*f(t-2.*dt,x) )
!
!  See, e.g.,
!  http://en.wikipedia.org/wiki/Linear_multistep_method
!
!  or
!
!  K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
!
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           ! time step number
      TYPE(ID_InputType),             INTENT(INOUT)  :: u(:)        ! Inputs at t
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   ! times of input
      TYPE(ID_ParameterType),         INTENT(IN   )  :: p           ! Parameters
      TYPE(ID_ContinuousStateType),   INTENT(INOUT)  :: x           ! Continuous states at t on input at t + dt on output
      TYPE(ID_DiscreteStateType),     INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(ID_ConstraintStateType),   INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
      TYPE(ID_OtherStateType),        INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables

      TYPE(ID_InputType)            :: u_interp        ! Continuous states at t
      TYPE(ID_ContinuousStateType)  :: x_pred          ! Continuous states at t
      TYPE(ID_ContinuousStateType)  :: xdot_pred       ! Continuous states at t

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL ID_CopyContState(x, x_pred, 0, ErrStat, ErrMsg)

      CALL ID_AB4( t, n, u, utimes, p, x_pred, xd, z, OtherState, ErrStat, ErrMsg )

      if (n .gt. 2) then

         CALL ID_Input_ExtrapInterp(u, utimes, u_interp, t + p%dt, ErrStat, ErrMsg)

         CALL ID_CalcContStateDeriv(t + p%dt, u_interp, p, x_pred, xd, z, OtherState, xdot_pred, ErrStat, ErrMsg )

         x%q    = x%q    + (p%dt / 24.) * ( 9. * xdot_pred%q +  19. * OtherState%xdot(1)%q - 5. * OtherState%xdot(2)%q &
                                          + 1. * OtherState%xdot(3)%q )

         x%dqdt = x%dqdt + (p%dt / 24.) * ( 9. * xdot_pred%dqdt + 19. * OtherState%xdot(1)%dqdt - 5. * OtherState%xdot(2)%dqdt &
                                          + 1. * OtherState%xdot(3)%dqdt )

      else

         x%q    = x_pred%q
         x%dqdt = x_pred%dqdt

       endif

END SUBROUTINE ID_ABM4
!----------------------------------------------------------------------------------------------------------------------------------
!..................................................................................................................................
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! WE ARE NOT YET IMPLEMENTING THE JACOBIANS...
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
END MODULE IceDyn
!**********************************************************************************************************************************
