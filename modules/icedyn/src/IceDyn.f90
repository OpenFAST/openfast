!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013-2016  Univeristy of Michigan and National Renewable Energy Laboratory 
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
!>    IceDyn is a module describing ice load on offshore wind turbine supporting structures.
!!
!!    The module is given the module name ModuleName = IceDyn and the abbreviated name ModName = IceD. The mathematical
!!    formulation of this module is a subset of the most general form permitted by the FAST modularization framework in tight
!!    coupling, thus, the module is developed to support both loose and tight coupling (tight coupling for both time marching and
!!    linearization).
!!   
!!    References:
!!
!!    Ice Module Manual, by Bingbin Yu, Dale Karr.
!!
!**********************************************************************************************************************************
MODULE IceDyn

   USE IceDyn_Types
   USE NWTC_Library

   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER  :: IceD_Ver = ProgDesc( 'IceDyn', '', '' )

   ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: IceD_Init                           ! Initialization routine
   PUBLIC :: IceD_End                            ! Ending routine (includes clean up)

   PUBLIC :: IceD_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                                 !   continuous states, and updating discrete states
   PUBLIC :: IceD_CalcOutput                     ! Routine for computing outputs

   PUBLIC :: IceD_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: IceD_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: IceD_UpdateDiscState                ! Tight coupling routine for updating discrete states

CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the start of the simulation to perform initialization steps.
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
SUBROUTINE IceD_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(IceD_InitInputType),       INTENT(IN   )  :: InitInp     !< Input data for initialization routine
   TYPE(IceD_InputType),           INTENT(  OUT)  :: u           !< An initial guess for the input; input mesh must be defined
   TYPE(IceD_ParameterType),       INTENT(  OUT)  :: p           !< Parameters
   TYPE(IceD_ContinuousStateType), INTENT(  OUT)  :: x           !< Initial continuous states
   TYPE(IceD_DiscreteStateType),   INTENT(  OUT)  :: xd          !< Initial discrete states
   TYPE(IceD_ConstraintStateType), INTENT(  OUT)  :: z           !< Initial guess of the constraint states
   TYPE(IceD_OtherStateType),      INTENT(  OUT)  :: OtherState  !< Initial other states
   TYPE(IceD_OutputType),          INTENT(  OUT)  :: y           !< Initial system outputs (outputs are not calculated;
                                                                 !!   only the output mesh is initialized)
   TYPE(IceD_MiscVarType),         INTENT(  OUT)  :: m           !< Initial misc/optimization variables
   REAL(DbKi),                     INTENT(INOUT)  :: Interval    !! Coupling interval in seconds: the rate that
                                                                 !!   (1) IceD_UpdateStates() is called in loose coupling &
                                                                 !!   (2) IceD_UpdateDiscState() is called in tight coupling.
                                                                 !!   Input is the suggested time from the glue code;
                                                                 !!   Output is the actual coupling interval that will be used
                                                                 !!   by the glue code.
   TYPE(IceD_InitOutputType),      INTENT(  OUT)  :: InitOut     !< Output for initialization routine
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


      ! Local variables

   TYPE(IceD_InputFile)                           :: InputFileData           ! Data stored in the module's input file
   INTEGER(IntKi)                                 :: ErrStat2                ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                           :: ErrMsg2                 ! temporary Error message if ErrStat /= ErrID_None
   REAL(ReKi)                                     :: TmpPos(3)


      ! Initialize variables for this routine

   ErrStat = ErrID_None
   ErrMsg  = ""


      ! Initialize the NWTC Subroutine Library

   CALL NWTC_Init( EchoLibVer=.FALSE. )


      ! Display the module information

   IF (InitInp%LegNum == 1) CALL DispNVD( IceD_Ver )  ! We don't need to see the version info more than one time

      !............................................................................................
      ! Read the input file and validate the data
      !............................................................................................
   p%RootName = TRIM(InitInp%RootName) ! FAST makes all of the output file names from this module end contain '.IceD' before the file extension

   CALL IceD_ReadInput( InitInp, InputFileData, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

    CALL IceD_ValidateInput( InputFileData, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN

      ! make sure the InitInp%LegNum makes sense:
   IF ( InputFileData%NumLegs < InitInp%LegNum .OR. InitInp%LegNum < 1 ) THEN
      CALL CheckError( ErrID_Fatal, 'IceDyn input file specified '//trim(num2lstr(InputFileData%NumLegs))// &
                        ' legs. The glue/driver code requested data for leg '//trim(num2lstr(InitInp%LegNum))//'.')
      RETURN
   END IF


      !............................................................................................
      ! Define parameters here:
      !............................................................................................
   CALL IceD_SetParameters( InputFileData, p, Interval, InitInp%Tmax, InitInp%LegNum, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN



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
   !CALL IceD_Init_DiscrtStates( xd, p, InputFileData, ErrStat2, ErrMsg2 )     ! initialize the continuous states
   !   CALL CheckError( ErrStat2, ErrMsg2 )
   !   IF (ErrStat >= AbortErrLev) RETURN

   xd%DummyDiscState      = 0

      ! Initialize other states:
    CALL IceD_Init_OtherStates( OtherState, p, x, InputFileData, ErrStat2, ErrMsg2 )    ! initialize the other states
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN

    m%DummyMiscVar = 0
    
      !............................................................................................
      ! Define initial guess for the system inputs and output (set up meshes) here:
      !............................................................................................

      ! set up meshes for inputs and outputs:

      ! Define system output initializations (set up mesh) here:
       CALL MeshCreate( BlankMesh      = u%PointMesh            &
                      ,IOS             = COMPONENT_INPUT        &
                      ,NNodes          = 1                      &
                      ,TranslationDisp = .TRUE.                 &
                      ,TranslationVel  = .TRUE.                 &
                      ,nScalars        = 0                      &
                      ,ErrStat         = ErrStat                &
                      ,ErrMess         = ErrMsg                 )

       CALL MeshConstructElement ( Mesh = u%PointMesh           &
                                , Xelement = ELEMENT_POINT      &
                                , P1       = 1                  &
                                , ErrStat  = ErrStat            &
                                , ErrMess  = ErrMsg             )

      ! place single node at water level; position affects mapping/coupling with other modules
      TmpPos(1) = InputFileData%LegPosX( InitInp%LegNum )
      TmpPos(2) = InputFileData%LegPosY( InitInp%LegNum )
      TmpPos(3) = InitInp%MSL2SWL

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
      u%PointMesh%TranslationDisp = 0.0_ReKi    ! initialize all components of all (1) points
      u%PointMesh%TranslationVel  = 0.0_ReKi    ! initialize all components of all (1) points  [bjj: does not appear to be used...]

      y%PointMesh%Force = 0.0_ReKi              ! initialize all components of all (1) points



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

   CALL AllocAry( y%WriteOutput, p%NumOuts, 'WriteOutput', ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN
            
      
   InitOut%WriteOutputHdr = p%OutName
   InitOut%WriteOutputUnt = p%OutUnit

   InitOut%Ver     = IceD_Ver
   InitOut%numLegs = InputFileData%NumLegs
   
   !bjj todo: check that water density is the same
   
   IF ( .NOT. EqualRealNos( InputFileData%rhow, InitInp%WtrDens ) ) THEN
      CALL CheckError( ErrID_Warn, 'IceD_Init: water density from IceDyn input file ('//trim(num2Lstr(InputFileData%rhow))//& 
                                ' kg/m^3) differs from water density in glue code ('//trim(num2Lstr(InitInp%WtrDens))//' kg/m^3).')
   END IF
      
   IF ( .NOT. EqualRealNos( 9.81_ReKi, InitInp%gravity ) ) THEN
      CALL CheckError( ErrID_Warn, 'IceD_Init: gravity hard-coded in IceDyn ('//trim(num2Lstr(9.81))//& 
                                       ' m/s^2) differs from gravity in glue code ('//trim(num2Lstr(InitInp%gravity))//' m/s^2).')
   END IF
   
       

  !     ! Print the summary file if requested:
  ! IF (InputFileData%SumPrint) THEN
  !    CALL IceD_PrintSum( p, OtherState, ErrStat2, ErrMsg2 )
  !       CALL CheckError( ErrStat2, ErrMsg2 )
  !       IF (ErrStat >= AbortErrLev) RETURN
  ! END IF

       ! Destroy the InputFileData structure (deallocate arrays)

   CALL IceD_DestroyInputFile(InputFileData, ErrStat2, ErrMsg2 )
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

         IF ( LEN_TRIM(ErrMsg) > 0 ) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//' '//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL IceD_DestroyInputFile(InputFileData, ErrStat3, ErrMsg3 )
         END IF

      END IF


   END SUBROUTINE CheckError

END SUBROUTINE IceD_Init
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
SUBROUTINE IceD_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(IceD_InputType),           INTENT(INOUT)  :: u           !< System inputs
      TYPE(IceD_ParameterType),       INTENT(INOUT)  :: p           !< Parameters
      TYPE(IceD_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states
      TYPE(IceD_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Discrete states
      TYPE(IceD_ConstraintStateType), INTENT(INOUT)  :: z           !< Constraint states
      TYPE(IceD_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states
      TYPE(IceD_OutputType),          INTENT(INOUT)  :: y           !< System outputs
      TYPE(IceD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""

      ! Place any last minute operations or calculations here:

      ! Close files here:

      ! Destroy the input data:

      CALL IceD_DestroyInput( u, ErrStat, ErrMsg )

      ! Destroy the parameter data:

      CALL IceD_DestroyParam( p, ErrStat, ErrMsg )

      ! Destroy the state data:

      CALL IceD_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL IceD_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL IceD_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL IceD_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )

      CALL IceD_DestroyMisc( m, ErrStat, ErrMsg )
      
      ! Destroy the output data:

      CALL IceD_DestroyOutput( y, ErrStat, ErrMsg )


END SUBROUTINE IceD_End
!----------------------------------------------------------------------------------------------------------------------------------
!> This is a loose coupling routine for solving constraint states, integrating continuous states, and updating discrete and other 
!! states. Continuous, constraint, discrete, and other states are updated to values at t + Interval.
SUBROUTINE IceD_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   ) :: t          !< Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   ) :: n          !< Current simulation time step n = 0,1,...
      TYPE(IceD_InputType),               INTENT(INOUT) :: u(:)       !< Inputs at utimes
      REAL(DbKi),                         INTENT(IN   ) :: utimes(:)  !< Times associated with u(:), in seconds
      TYPE(IceD_ParameterType),           INTENT(IN   ) :: p          !< Parameters
      TYPE(IceD_ContinuousStateType),     INTENT(INOUT) :: x          !< Input: Continuous states at t;
                                                                      !!   Output: Continuous states at t + Interval
      TYPE(IceD_DiscreteStateType),       INTENT(INOUT) :: xd         !< Input: Discrete states at t;
                                                                      !!   Output: Discrete states at t  + Interval
      TYPE(IceD_ConstraintStateType),     INTENT(INOUT) :: z          !< Input: Constraint states at t;
                                                                      !!   Output: Constraint states at t+dt
      TYPE(IceD_OtherStateType),          INTENT(INOUT) :: OtherState !< Input: Other states at t;
                                                                      !!   Output: Other states at t+dt
      TYPE(IceD_MiscVarType),             INTENT(INOUT) :: m          !< Misc/optimization variables
      INTEGER(IntKi),                     INTENT(  OUT) :: ErrStat    !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT) :: ErrMsg     !< Error message if ErrStat /= ErrID_None

      ! local variables

      TYPE(IceD_InputType)                              :: u_interp      ! input interpolated from given u at utimes
      !TYPE(IceD_ContinuousStateType)                   :: xdot          ! continuous state time derivative
      INTEGER(IntKi)                                    :: I             ! Loop count
      REAL(ReKi)                                        :: Del2          ! Deflection of the current ice tooth, for model 2,3
      REAL(ReKi)                                        :: Del(p%Zn)     ! Deflection of ice tooth in each zone, for model 4
!      REAL(ReKi)                                        :: StrRt         ! Strain rate (s^-1)
!      REAL(ReKi)                                        :: SigCrp        ! Creep stress (Pa)
      
      INTEGER(IntKi)                                    :: nt            ! Current time step
      INTEGER(IntKi)                                    :: Nc            ! Current ice tooth number
      REAL(ReKi)                                        :: Psumc         ! Current sum of pitch of broken ice teeth (m)
      REAL(ReKi)                                        :: Dmc           ! Delmax of the current tooth (m)
      REAL(ReKi)                                        :: Pc            ! Pitch of the current tooth (m)
      REAL(ReKi)                                        :: Pnxt          ! Pitch of the next ice tooth (m)
      
      INTEGER(IntKi)                                    :: ErrStat2
      CHARACTER(ErrMsgLen)                              :: ErrMsg2
      
       ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""
      
      
      CALL IceD_CopyInput( u(1), u_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg, 'IceD_UpdateStates')
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF
                     
      ! interpolate u to find u_interp = u(t)
      CALL IceD_Input_ExtrapInterp( u, utimes, u_interp, t, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg, 'IceD_UpdateStates')
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF
      
      
      
      nt = n+1   !bjj: was this: INT( t / p%dt ) + 1

      ! Update Other states here for ice model 2
      
      IF (p%ModNo == 2) THEN

        Del2 = p%InitLoc + p%v * t - p%Pitch * (OtherState%IceTthNo2 -1) - u_interp%PointMesh%TranslationDisp(1,1)

            IF ( Del2 >= (p%Delmax2 - p%tolerance) ) THEN

                OtherState%IceTthNo2 = OtherState%IceTthNo2 + 1
                !CALL WrScr( 'IceTthNo=' // Num2LStr(xd%IceTthNo2) )

            ENDIF

      ENDIF
      
      ! Update Other states here for ice model 3
     
      IF (p%ModNo == 3) THEN
         Nc = OtherState%Nc (nt)
         Psumc = OtherState%Psum (nt)
         Pc = p%rdmP (Nc)
      
         IF (p%SubModNo == 3) THEN
      
            Del2 = p%InitLoc + p%v * t - OtherState%Psum (nt) - u_interp%PointMesh%TranslationDisp(1,1) ! Deflection of the current ice tooth
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
      
      ! Update Other States here for ice model 4
      
      IF (p%ModNo == 4) THEN
      
        DO I = 1, p%Zn
            
            Del (I) = p%Y0 (I) + p%v * t - p%ZonePitch * (OtherState%IceTthNo (I)-1) - u_interp%PointMesh%TranslationDisp(1,1)
            
            IF ( Del(I) >= (p%Delmax - p%tolerance) ) THEN 
                
                OtherState%IceTthNo (I) = OtherState%IceTthNo (I)+1
                
            ENDIF
            
        END DO
      
      END IF
      
      ! Update Other States here for ice model 5
      
      IF (p%ModNo == 5) THEN
      
            OtherState%Beta   = SolveBeta( p%alphaR, p%v, t - OtherState%Tinit, p%Lbr)
            
            IF (OtherState%Beta >= p%alphaR) THEN
               OtherState%Tinit = t
            END IF
      
      ENDIF

      ! Update Continuous states here for ice model 6
      
      IF (p%ModNo == 6) THEN

           if (p%method .eq. 1) then

               CALL IceD_RK4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )

            elseif (p%method .eq. 2) then

               CALL IceD_AB4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )

            elseif (p%method .eq. 3) then

               CALL IceD_ABM4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )

            else

               ErrStat = ErrID_Fatal
               ErrMsg  = ' Error in IceD_UpdateStates: p%method must be 1 (RK4), 2 (AB4), or 3 (ABM4)'
               RETURN

            endif
            
            !IF ((x%q - u_interp%PointMesh%TranslationDisp(1,1)) >= OtherState%dxc) THEN
            !     OtherState%dxc = x%q - u_interp%PointMesh%TranslationDisp(1,1)
            !ENDIF
            
            ! Determine whether the splitting failure happens
            IF (OtherState%Splitf == 0) THEN
                CALL Ice_Split (t, u_interp, p, x, OtherState, ErrStat, ErrMsg)
            ENDIF

      END IF
           

      IF ( ErrStat >= AbortErrLev ) RETURN
      
      RETURN
CONTAINS      
   SUBROUTINE Cleanup()
   
   CALL IceD_DestroyInput( u_interp, ErrStat2, ErrMsg2 )

   END SUBROUTINE Cleanup
   
   FUNCTION SolveBeta (alpha, vice, t, l) Result (beta1)
         
         !SOLVEBETA Solve for ice wedge uprising angle
         
         IMPLICIT NONE
            
         ! Input values
         REAL(ReKi), intent(in) :: alpha     ! Cone angle (rad)
         REAL(ReKi), intent(in) :: vice      ! Ice velocity (m/s)
         REAL(DbKi), intent(in) :: t         ! Ice time (s)
         REAL(ReKi), intent(in) :: l         ! Ice breaking length (m)
            
         REAL(ReKi)             :: beta      ! Initial Ice wedge uprising angle
         REAL(ReKi)             :: beta1     ! Ice wedge uprising angle
                                   
         REAL(ReKi)             :: Equ    
         REAL(ReKi)             :: Derv
            
         beta = 0
         DO i = 1,100
            
            Equ  = sin(beta) - tan(alpha) * cos(beta) + tan(alpha) * (1-vice*t/l);
            Derv = cos(beta) + tan(alpha) * sin(beta); 
                 
            IF ( abs(Equ) <= p%tolerance) EXIT   
                 
            beta = beta - Equ / Derv
                  
         END DO
         
         beta1 = beta
               
    END FUNCTION SolveBeta
    
    SUBROUTINE Ice_Split (t, u, p, x, OtherState, ErrStat, ErrMsg)

      REAL(DbKi),                       INTENT(IN   )  :: t           ! Current simulation time in seconds
      TYPE(IceD_InputType),             INTENT(IN   )  :: u           ! Inputs at t
      TYPE(IceD_ParameterType),         INTENT(IN   )  :: p           ! Parameters
      TYPE(IceD_ContinuousStateType),   INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(IceD_OtherStateType),        INTENT(INOUT)  :: OtherState  ! Other states at t
      INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
      
      ! Local variable
      REAL(ReKi)                                     :: IceForce 
     REAL(ReKi)                            :: R
      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Compute outputs here:
      
      R = p%StrWd/2
      
      IF ( OtherState%Splitf == 0 ) THEN
          
       IF ((x%q - u%PointMesh%TranslationDisp(1,1)) < OtherState%dxc) THEN
          
           IceForce = 0

       ELSE IF ((x%q - u%PointMesh%TranslationDisp(1,1)) >= OtherState%dxc) THEN

           IF (x%q - u%PointMesh%TranslationDisp(1,1) < R) THEN
      
              IceForce =  p%Cpa * ( 2 * p%h * ( R**2 - (R - x%q + u%PointMesh%TranslationDisp(1,1))**2 )**0.5 )**( p%dpa + 1 ) * 1.0e6
           
           ELSE

              IceForce = p%Cpa * ( 2 * R *  p%h )**( p%dpa + 1 ) * 1.0e6 

           ENDIF

       ENDIF

    ELSE

      IceForce = 0 

    ENDIF
      
      IF ( IceForce >= p%Fsp ) THEN
          OtherState%Splitf = 1
      ENDIF 

    END SUBROUTINE Ice_Split

      
END SUBROUTINE IceD_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
SUBROUTINE IceD_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                      INTENT(IN   )  :: t           !< Current simulation time in seconds
      TYPE(IceD_InputType),            INTENT(IN   )  :: u           !< Inputs at t
      TYPE(IceD_ParameterType),        INTENT(IN   )  :: p           !< Parameters
      TYPE(IceD_ContinuousStateType),  INTENT(IN   )  :: x           !< Continuous states at t
      TYPE(IceD_DiscreteStateType),    INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(IceD_ConstraintStateType),  INTENT(IN   )  :: z           !< Constraint states at t
      TYPE(IceD_OtherStateType),       INTENT(IN   )  :: OtherState  !< Other states at t
      TYPE(IceD_OutputType),           INTENT(INOUT)  :: y           !< Outputs computed at t (Input only so that mesh 
                                                                     !! connectivity information does not have to be recalculated)
      TYPE(IceD_MiscVarType),          INTENT(INOUT)  :: m           !< Misc/optimization variables
      INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! local variables
      INTEGER(IntKi)                    :: I                        ! Loop count
      INTEGER(IntKi)                    :: nt                       ! Current time step
      REAL(ReKi)                        :: Del2                     ! Deflection of the current ice tooth, for model 2
      REAL(ReKi)                        :: Del(p%Zn)                ! Deflection of each ice tooth, for model 4
      REAL(ReKi)                        :: ZoneF(p%Zn)              ! Ice force of each ice tooth, for model 4
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
      REAL(ReKi)                        :: gamma
      REAL(ReKi)                        :: CrntT0
      REAL(IntKi)                       :: dummyN

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""

      nt = INT( t / p%dt ) + 1
      dummyN = 1
      
      ! Compute outputs here:

      ! Ice Model 1 --------------------------------
      IF (p%ModNo == 1) THEN

          IF (p%SubModNo == 1) THEN

              IF (t < p%t0 + p%tolerance) THEN

                  y%PointMesh%Force(1,1) = 0

              ELSEIF (t< p%t0 + p%tm1a + p%tolerance) THEN

                  y%PointMesh%Force(1,1) = (t-p%t0) / p%tm1a * p%Fmax1a

              ELSE
                  
                  y%PointMesh%Force(1,1) = p%Fmax1a

              ENDIF

          ELSEIF (p%SubModNo == 2) THEN

              IF (t < p%t0 + p%tolerance) THEN

                  y%PointMesh%Force(1,1) = 0.0_ReKi

              ELSEIF (t< p%t0 + p%tm1b + p%tolerance) THEN

                  y%PointMesh%Force(1,1) = (t-p%t0) / p%tm1b * p%Fmax1b

              ELSE
                  
                  y%PointMesh%Force(1,1) = 0.0_ReKi

              ENDIF

          ELSEIF (p%SubModNo == 3) THEN

              IF (t < p%t0 + p%tolerance) THEN

                  y%PointMesh%Force(1,1) = 0.0_ReKi

              ELSEIF (t< p%t0 + p%tm1c + p%tolerance) THEN

                  y%PointMesh%Force(1,1) = (t-p%t0) / p%tm1c * p%Fmax1c

              ELSE

                  y%PointMesh%Force(1,1) = p%Fmax1c

              ENDIF

          ENDIF

      ENDIF

      ! Ice Model 2 --------------------------------

      IF (p%ModNo == 2) THEN

          Del2  = p%InitLoc + p%v * t - p%Pitch * (OtherState%IceTthNo2 -1) - u%PointMesh%TranslationDisp(1,1)

            IF (p%Delmax2 <= p%Pitch) THEN !Sub-model 1 
            
                IF ( Del2 <= 0) THEN

                    y%PointMesh%Force(1,1) = 0.0_ReKi

                ELSEIF (Del2 < p%Delmax2 + p%tolerance) THEN

                    y%PointMesh%Force(1,1) = Del2 * p%Kice2

                ELSE

                    !ErrStat = ErrID_Fatal
                    y%PointMesh%Force(1,1) = Del2 * p%Kice2
                    CALL WrScr ( ' Warning in IceDyn Model 2: Ice tooth does not break when Del>Delmax')

                ENDIF
                
            ELSEIF (p%Delmax2 <= 2*p%Pitch) THEN !Sub-model 2, two ice teeth bending
            
                IF ( Del2 <= 0) THEN

                    y%PointMesh%Force(1,1) = 0.0_ReKi

                ELSEIF (Del2 < p%Pitch) THEN

                    y%PointMesh%Force(1,1) = Del2 * p%Kice2

                ELSEIF (Del2 < p%Delmax2 + p%tolerance) THEN 
                
                    y%PointMesh%Force(1,1) = Del2 * p%Kice2 + (Del2 - p%Pitch) * p%Kice2
                
                ELSE

                    !ErrStat = ErrID_Fatal
                    y%PointMesh%Force(1,1) = Del2 * p%Kice2 + (Del2 - p%Pitch) * p%Kice2
                    CALL WrScr ( ' Warning in IceDyn Model 2: Ice tooth does not break when Del>Delmax' )

                ENDIF
            
            ELSE
                ErrStat = ErrID_Fatal
                ErrMsg  = ' Error in IceDyn Model 2.2: input delmax > 2*pitch'
            ENDIF

      END IF
      
      IF (p%ModNo == 3) THEN

          IF (p%SubModNo == 1) THEN

              IF ( t <= p%rdmt0 (nt) ) THEN

                  y%PointMesh%Force(1,1) = 0.0_ReKi

              ELSEIF (t <= p%rdmt0 (nt) + p%rdmtm (nt)) THEN

                  y%PointMesh%Force(1,1) = (t-p%rdmt0 (nt)) /  p%rdmtm(nt) * p%rdmFm (nt)

              ELSE

                  y%PointMesh%Force(1,1) = p%rdmFm (nt)

              ENDIF

          ELSEIF (p%SubModNo == 2) THEN
              
              CrntT0 = p%rdmt0 (nt)

              IF (t <= p%rdmt0 (nt)) THEN

                  y%PointMesh%Force(1,1) = 0.0_ReKi

              ELSEIF (t <= p%rdmt0 (nt) + p%rdmtm (nt)) THEN

                  y%PointMesh%Force(1,1) = (t-p%rdmt0 (nt)) / p%rdmtm(nt) * p%rdmFm (nt)

              ELSE

                  y%PointMesh%Force(1,1) = 0.0_ReKi

              ENDIF

          ELSEIF (p%SubModNo == 3) THEN

              Del2  = p%InitLoc + p%v * t - OtherState%Psum (nt) - u%PointMesh%TranslationDisp(1,1)     ! Determine the contact state between ice sheet and the tower

                !IF (Del2 >= xd%Dmaxn) THEN
                !    ErrStat = ErrID_Fatal
                !    ErrMsg  = ' Error in IceDyn Model 3c: two ice teeth break at once'
                !ENDIF

              IF (Del2 <= 0) THEN

                y%PointMesh%Force(1,1) = 0.0_ReKi

                ELSE IF (Del2 > 0 .AND. Del2 <= p%rdmP (OtherState%Nc(nt)) ) THEN

                   y%PointMesh%Force(1,1) = p%rdmKi(OtherState%Nc(nt))  * Del2

                ELSE

                   y%PointMesh%Force(1,1) = p%rdmKi(OtherState%Nc(nt)) * Del2 + p%rdmKi(OtherState%Nc(nt+1)) * (Del2-p%rdmP (OtherState%Nc(nt))) ! Two teeth in contact

                ENDIF

          ENDIF
          
      ENDIF
      
      ! Ice Model 4 --------------------------------
      
      IF (p%ModNo == 4) THEN
      
        DO I = 1, p%Zn
            
            Del (I) = p%Y0 (I) + p%v * t - p%ZonePitch * (OtherState%IceTthNo (I)-1) - u%PointMesh%TranslationDisp(1,1)
            
            IF ( Del (I) <= 0) THEN 
                
                ZoneF (I) = 0.0
                
            ELSEIF (Del (I) < p%Delmax + p%tolerance) THEN
                
                ZoneF (I) = Del (I) * p%Kice
                
            ELSE 
                
                !ErrStat = ErrID_Fatal
                ZoneF (I) = Del (I) * p%Kice
                CALL WrScr (' Warning in IceDyn Model 4: ZoneDel > Delmax')
                
            ENDIF          
            
        END DO
        
        IF (sum(ZoneF) <0) THEN
        
           dummyN = dummyN + 1
           
        ENDIF
        
        y%PointMesh%Force(1,1) = sum(ZoneF)
      
      END IF
      
      ! Ice Model 5 --------------------------------
      
      IF (p%ModNo == 5) THEN
      
           IF (t <= OtherState%Tinit) THEN
           
            y%PointMesh%Force(1,1) = 0
           
           ELSE IF (t <= OtherState%Tinit + p%dt) THEN  ! Ice load at breakage
           
            y%PointMesh%Force(1,1) = p%RHbr
           
           ELSE ! Ice load after breakage
           
            Wr = p%Wri * ( p%Zr - p%Lbr * sin(OtherState%Beta) ) / p%Zr
            Pn1 = Wr * cos(p%alphaR)
            gamma = p%rhoi / p%rhow
            
            IF (OtherState%Beta < gamma * p%h / p%Lbr) THEN 
        
                Xb = p%Lbr /3.0 *( ( 3.0*gamma*p%h - p%Lbr*tan(OtherState%Beta) ) / ( 2.0 *gamma*p%h - p%Lbr*tan(OtherState%Beta) ) )
                Fb = p%rhow * 9.81 * p%Dwl * (0.5 * p%Lbr**2 * tan(OtherState%Beta) + p%Lbr*(gamma*p%h - p%Lbr*tan(OtherState%Beta)) )
        
            ELSE
        
                Xb = 1.0/3.0 * gamma * p%h / sin(OtherState%Beta)
                Fb = p%rhow * 9.81 * p%Dwl * (0.5 * (gamma*p%h)**2 / tan(OtherState%Beta) )
        
            END IF
    
            pbeta = sin(OtherState%Beta) * ( sin(p%alphaR) + p%mu*cos(p%alphaR) ) + cos(OtherState%Beta) * ( cos(p%alphaR) - p%mu*sin(p%alphaR) )
            qbeta = ( sin(p%alphaR) + p%mu*cos(p%alphaR) ) * sin( p%alphaR - OtherState%Beta )
    
            Pn2 = ( p%WL * p%Lbr /2.0 * cos(OtherState%Beta) + Wr * p%Lbr * qbeta - Fb*Xb) / pbeta / p%Lbr
    
            Rh = ( Pn1 + Pn2 ) * ( sin(p%alphaR) + p%mu*cos(p%alphaR) )
            Rv = ( Pn1 + Pn2 ) * ( cos(p%alphaR) - p%mu*sin(p%alphaR) )
            
            IF (Rh < 0) THEN
                y%PointMesh%Force(1,1) = 0                
            ELSE                 
                y%PointMesh%Force(1,1) = Rh
            ENDIF 
           
           ENDIF
      
      ENDIF
      
      ! Ice Model 6 --------------------------------
      IF (p%ModNo == 6) THEN
      
           R = p%StrWd/2
      
          IF ( OtherState%Splitf == 0 ) THEN
          
              IF ((x%q - u%PointMesh%TranslationDisp(1,1)) < OtherState%dxc) THEN
          
                 y%PointMesh%Force(1,1) = 0

              ELSE IF ((x%q - u%PointMesh%TranslationDisp(1,1)) >= OtherState%dxc) THEN

                   IF (x%q - u%PointMesh%TranslationDisp(1,1) < R) THEN
      
                        y%PointMesh%Force(1,1) =  p%Cpa * ( 2 * p%h * ( R**2 - (R - x%q + u%PointMesh%TranslationDisp(1,1))**2 )**0.5 )**( p%dpa + 1 ) * 1.0e6
           
                   ELSE

                        y%PointMesh%Force(1,1) = p%Cpa * ( 2 * R *  p%h )**( p%dpa + 1 ) * 1.0e6 

                   ENDIF

              ENDIF

          ELSE

                 y%PointMesh%Force(1,1) = 0 

          ENDIF
      
      ENDIF
      
      !y%PointMesh%Force(1,1) = y%fice
      !y%PointMesh%Force(2,1) = 0.
      !y%PointMesh%Force(3,1) = 0.
            
      
      ! values to write to a file:
      y%WriteOutput(1) = x%q                     ! IceDisp  !bjj: todo: do we need to recalculate this???
      y%WriteOutput(2) = y%PointMesh%Force(1,1)  ! IceForce
      
      
END SUBROUTINE IceD_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing derivatives of continuous states.
SUBROUTINE IceD_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, xdot, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                       INTENT(IN   )  :: t           !< Current simulation time in seconds
      TYPE(IceD_InputType),             INTENT(IN   )  :: u           !< Inputs at t
      TYPE(IceD_ParameterType),         INTENT(IN   )  :: p           !< Parameters
      TYPE(IceD_ContinuousStateType),   INTENT(IN   )  :: x           !< Continuous states at t
      TYPE(IceD_DiscreteStateType),     INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(IceD_ConstraintStateType),   INTENT(IN   )  :: z           !< Constraint states at t
      TYPE(IceD_OtherStateType),        INTENT(IN   )  :: OtherState  !< Other states at t
      TYPE(IceD_MiscVarType),           INTENT(INOUT)  :: m           !< Misc/optimization variables
      TYPE(IceD_ContinuousStateType),   INTENT(  OUT)  :: xdot        !< Continuous state derivatives at t
      INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! local variables
      ! REAL(ReKi) :: mass    ! Mass of ice feature (kg)
      REAL(ReKi) :: force     ! Ice force (N)
      REAL(ReKi) :: R       ! Structure radius
      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""
      
      IF (p%ModNo == 6) THEN

      ! Compute the first time derivatives of the continuous states here:

     ! When the ice and the structure is in contact, there is ice force.     

         R = p%StrWd/2 
        
         IF ( OtherState%Splitf == 0 ) THEN
          
           IF ((x%q - u%PointMesh%TranslationDisp(1,1)) < OtherState%dxc) THEN
          
               force = 0 + p%FdrN

           ELSE IF ((x%q - u%PointMesh%TranslationDisp(1,1)) >= OtherState%dxc) THEN

                IF (x%q - u%PointMesh%TranslationDisp(1,1) < R) THEN
      
                    force = -p%Cpa * ( 2 * p%h * ( R**2 - (R - x%q + u%PointMesh%TranslationDisp(1,1))**2 )**0.5 )**( p%dpa + 1 ) * 1.0e6 + p%FdrN
           
                ELSE

                    force = -p%Cpa * ( 2 * R *  p%h )**( p%dpa + 1 ) * 1.0e6 + p%FdrN

                ENDIF
 
           ENDIF

         ELSE

              force = 0 
 
         ENDIF
        
         xdot%q = x%dqdt
   
         xdot%dqdt = force / p%Mice
      
      ENDIF
      

END SUBROUTINE IceD_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for updating discrete states
SUBROUTINE IceD_UpdateDiscState( t, n, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                      INTENT(IN   )  :: t           !< Current simulation time in seconds
      INTEGER(IntKi),                  INTENT(IN   )  :: n           !< Current step of the simulation: t = n*Interval
      TYPE(IceD_InputType),            INTENT(IN   )  :: u           !< Inputs at t
      TYPE(IceD_ParameterType),        INTENT(IN   )  :: p           !< Parameters
      TYPE(IceD_ContinuousStateType),  INTENT(IN   )  :: x           !< Continuous states at t
      TYPE(IceD_DiscreteStateType),    INTENT(INOUT)  :: xd          !< Input: Discrete states at t;
                                                                     !<   Output: Discrete states at t + Interval
      TYPE(IceD_ConstraintStateType),  INTENT(IN   )  :: z           !< Constraint states at t
      TYPE(IceD_OtherStateType),       INTENT(IN   )  :: OtherState  !< Other states at t
      TYPE(IceD_MiscVarType),          INTENT(INOUT)  :: m           !< Misc/optimization variables
      INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! local variables

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""

      ! Update discrete states here:
            xd%DummyDiscState = 0.0

END SUBROUTINE IceD_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for solving for the residual of the constraint state functions
SUBROUTINE IceD_CalcConstrStateResidual( t, u, p, x, xd, z, OtherState, m, Z_residual, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                      INTENT(IN   )  :: t           !< Current simulation time in seconds
      TYPE(IceD_InputType),            INTENT(IN   )  :: u           !< Inputs at t
      TYPE(IceD_ParameterType),        INTENT(IN   )  :: p           !< Parameters
      TYPE(IceD_ContinuousStateType),  INTENT(IN   )  :: x           !< Continuous states at t
      TYPE(IceD_DiscreteStateType),    INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(IceD_ConstraintStateType),  INTENT(IN   )  :: z           !< Constraint states at t
      TYPE(IceD_OtherStateType),       INTENT(IN   )  :: OtherState  !< Other states at t
      TYPE(IceD_ConstraintStateType),  INTENT(  OUT)  :: Z_residual  !< Residual of the constraint state functions using
                                                                     !!    the input values described above
      TYPE(IceD_MiscVarType),          INTENT(INOUT)  :: m           !< Misc/optimization variables
      INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


      ! Solve for the constraint states here:

      Z_residual%DummyConstrState = 0

END SUBROUTINE IceD_CalcConstrStateResidual
!----------------------------------------------------------------------------------------------------------------------------------
!>     This public subroutine reads the input required for IceDyn from the file whose name is an
!!     input parameter.
SUBROUTINE IceD_ReadInput( InitInp, InputFileData, ErrStat, ErrMsg )
!----------------------------------------------------------------------------------------------------

      ! Passed variables

   TYPE(IceD_InitInputType),      INTENT( IN    )   :: InitInp              !< the IceDyn initial input data
   TYPE(IceD_InputFile),          INTENT(   OUT )   :: InputFileData        !< Data stored in the IceDyn's input file
   INTEGER,                       INTENT(   OUT )   :: ErrStat              !< returns a non-zero value when an error occurs
   CHARACTER(*),                  INTENT(   OUT )   :: ErrMsg               !< Error message if ErrStat /= ErrID_None


      ! Local variables

   INTEGER                                          :: UnIn                 ! Unit number for the input file
   CHARACTER(ErrMsgLen)                             :: FileName             ! Name of HydroDyn input file
   
   INTEGER                                          :: UnEc                 ! Unit number for the echo file
   LOGICAL, PARAMETER                               :: Echo = .FALSE.       ! echo file for debugging
!  LOGICAL, PARAMETER                               :: Echo = .TRUE.        ! echo file for debugging (bjj: would like to add this feature to the input file)



      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""
   UnEc    = -1
   UnIn    = -1

   !-------------------------------------------------------------------------------------------------
   ! Open the file
   !-------------------------------------------------------------------------------------------------
   IF ( Echo ) THEN
      CALL GetNewUnit( UnEc, ErrStat, ErrMsg )      
      CALL OpenFOutFile( UnEc, TRIM(InitInp%RootName)//'.IceD.ech', ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN
         CALL WrScr( ' Error opening echo file: "'//TRIM(ErrMsg)//'". Simulation will continue with no IceDyn echo file.' )
         CLOSE( UnEc )
         UnEc = -1
      END IF
   END IF   
   
   
   
   FileName = TRIM(InitInp%InputFile)

   CALL GetNewUnit( UnIn, ErrStat, ErrMsg )
   CALL OpenFInpFile( UnIn, FileName, ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF


   !CALL WrScr( 'Opening IceDyn input file:  '//FileName )


   !-------------------------------------------------------------------------------------------------
   ! File header
   !-------------------------------------------------------------------------------------------------

   CALL ReadCom( UnIn, FileName, 'IceDyn input file header line 1', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

   CALL ReadCom( UnIn, FileName, 'IceDyn input file header line 2', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Structure properties
   !-------------------------------------------------------------------------------------------------

   CALL ReadCom( UnIn, FileName, 'IceDyn structure properties header line', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF
   
      ! NumLegs - Number of support-structure legs in contact with ice

   CALL ReadVar ( UnIn, FileName, InputFileData%NumLegs, 'NumLegs', 'Number of support-structure legs in contact with ice', ErrStat, ErrMsg, UnEc )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

   IF ( InputFileData%NumLegs < 0 ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg = 'IceD_ReadInput: NumLegs must be a positive number.'
      CALL Cleanup()
      RETURN
   END IF
   
   CALL AllocAry( InputFileData%LegPosX,  InputFileData%NumLegs, 'LegPosX', ErrStat, ErrMsg ) 
   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

   CALL AllocAry( InputFileData%LegPosY,  InputFileData%NumLegs, 'LegPosY', ErrStat, ErrMsg ) 
   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

   CALL AllocAry( InputFileData%StrWd,  InputFileData%NumLegs, 'StrWd', ErrStat, ErrMsg ) 
   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF
   
   
   ! LegPosX  - Global X position of legs 1-NumLegs (m)
   
   CALL ReadAry ( UnIn, FileName, InputFileData%LegPosX, InputFileData%NumLegs, 'LegPosX', 'Global X position of legs, 1-NumLegs (m)', ErrStat, ErrMsg, UnEc )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

   ! LegPosY  - Global Y position of legs 1-NumLegs (m)
   
   CALL ReadAry ( UnIn, FileName, InputFileData%LegPosY, InputFileData%NumLegs, 'LegPosY', 'Global Y position of legs, 1-NumLegs (m)', ErrStat, ErrMsg, UnEc )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

   
   ! StWidth  - Width of the structure in contact with the ice (m)

   CALL ReadAry ( UnIn, FileName, InputFileData%StrWd, InputFileData%NumLegs, 'StWidth', 'Width of the structure in contact with the ice (m)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF
   
   
   !-------------------------------------------------------------------------------------------------
   ! Ice Model Number
   !-------------------------------------------------------------------------------------------------

      ! Header

   CALL ReadCom( UnIn, FileName, 'Ice Models header', ErrStat, ErrMsg, UnEc)

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF


      ! IceModel - Number that represents different ice models.

   CALL ReadVar ( UnIn, FileName, InputFileData%IceModel, 'IceModel', 'Number that represents different ice models', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF


    ! IceSubModel - Number that represents different ice sub models.

   CALL ReadVar ( UnIn, FileName, InputFileData%IceSubModel, 'IceSubModel', 'Number that represents different ice sub-models', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Ice properties - General
   !-------------------------------------------------------------------------------------------------

      ! Header

   CALL ReadCom( UnIn, FileName, 'Ice properties - General header', ErrStat, ErrMsg, UnEc)

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

     ! IceVel - Velocity of ice sheet movement (m/s)

   CALL ReadVar ( UnIn, FileName, InputFileData%v, 'IceVel', 'Velocity of ice sheet movement', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

   ! IceThks  - Thickness of the ice sheet (m)

   CALL ReadVar ( UnIn, FileName, InputFileData%h, 'IceThks', 'Thickness of the ice sheet (m)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF


      ! WtDen     - Mass density of water (kg/m3)

   CALL ReadVar ( UnIn, FileName, InputFileData%rhow, 'WtDen', 'Mass density of water (kg/m3)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read WtDen parameter.'
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

      ! IceDen  - Mass density of ice (kg/m3)

   CALL ReadVar ( UnIn, FileName, InputFileData%rhoi, 'IceDen', 'Mass density of ice (kg/m3)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

          ! InitLoc - Ice sheet initial location (m)

   CALL ReadVar ( UnIn, FileName, InputFileData%InitLoc, 'InitLoc', 'Ice sheet initial location (m)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

          ! InitTm - Ice load starting time (s)

   CALL ReadVar ( UnIn, FileName, InputFileData%t0, 'InitTm', 'Ice load starting time (s)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF
   
          ! Seed1 - Random seed 1

   CALL ReadVar ( UnIn, FileName, InputFileData%Seed1, 'Seed1', 'Random seed 1', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF
   
           ! Seed2 - Random seed 2

   CALL ReadVar ( UnIn, FileName, InputFileData%Seed2, 'Seed2', 'Random seed 2', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Ice properties - Ice Model 1
   !-------------------------------------------------------------------------------------------------

        ! Header

   CALL ReadCom( UnIn, FileName, 'Ice properties - Ice Model 1 SubModel 1', ErrStat, ErrMsg, UnEc)

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

    ! Ikm - Indentation factor

   CALL ReadVar ( UnIn, FileName, InputFileData%Ikm, 'Ikm', 'Indentation factor', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

    ! Ag - Constant depends only on ice crystal type, used in calculating uniaxial stress (MPa-3s-1)

   CALL ReadVar ( UnIn, FileName, InputFileData%Ag, 'Ag', 'Constant depends only on ice crystal type, used in calculating uniaxial stress (MPa-3s-1)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

    ! Qg - Activation Energy (kJmol^-1)

   CALL ReadVar ( UnIn, FileName, InputFileData%Qg, 'Qg', 'Activation Energy (kJmol^-1)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

    ! Rg - Universal gas constant (Jmol-1K-1)

   CALL ReadVar ( UnIn, FileName, InputFileData%Rg, 'Rg', 'Universal gas constant (Jmol-1K-1)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

    ! Tice - Ice temperature (K)

   CALL ReadVar ( UnIn, FileName, InputFileData%Tice, 'Tice', 'Ice temperature (K)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF


    ! Header

   CALL ReadCom( UnIn, FileName, 'Ice properties - Ice Model 1 SubModel 2', ErrStat, ErrMsg, UnEc)

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

    ! Poisson  - Poisson's ratio of ice

   CALL ReadVar ( UnIn, FileName, InputFileData%nu, 'Poisson', 'Poisson ratio of ice', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

    ! WgAngle - Wedge Angel, degree. Default 45 Degrees.

   CALL ReadVar ( UnIn, FileName, InputFileData%phi, 'WgAngle', ' Wedge Angel, degree. Default 45 Degrees', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

    ! EIce - Young's modulus of ice (GPa)

   CALL ReadVar ( UnIn, FileName, InputFileData%Eice, 'EIce', 'Youngs modulus of ice (GPa)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

    ! Header

   CALL ReadCom( UnIn, FileName, 'Ice properties - Ice Model 1 SubModel 3', ErrStat, ErrMsg, UnEc)

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

    ! SigN - Nominal ice stress (MPa)

   CALL ReadVar ( UnIn, FileName, InputFileData%SigNm, 'SigNm', 'Nominal ice stress (MPa)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Ice properties - Ice Model 2
   !-------------------------------------------------------------------------------------------------

   ! Header

   CALL ReadCom( UnIn, FileName, 'Ice properties - Ice Model 2 SubModel 1,2', ErrStat, ErrMsg, UnEc)

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

    ! Pitch - Distance between sequential ice teeth (m)

   CALL ReadVar ( UnIn, FileName, InputFileData%Pitch, 'Pitch', 'Distance between sequential ice teeth (m)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

    ! IceStr2 - Ice failure stress (MPa)

   CALL ReadVar ( UnIn, FileName, InputFileData%IceStr2, 'IceStr2', 'Ice failure stress (MPa)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

    ! Delmax2 - Ice tooth maximum elastic deformation (m)

   CALL ReadVar ( UnIn, FileName, InputFileData%Delmax2, 'Delmax2', 'Ice tooth maximum elastic deformation (m)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Ice properties - Ice Model 3
   !-------------------------------------------------------------------------------------------------

     ! Header

   CALL ReadCom( UnIn, FileName, 'Ice PROPERTIES -Ice Model 3, SubModel 1,2', ErrStat, ErrMsg, UnEc)

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

   ! ThkMean - Mean value of ice thickness (m)

   CALL ReadVar ( UnIn, FileName, InputFileData%miuh, 'ThkMean', 'Mean value of ice thickness (m)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

   ! ThkVar - Variance of ice thickness (m^2)

   CALL ReadVar ( UnIn, FileName, InputFileData%varh, 'ThkVar', 'Variance of ice thickness (m^2)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

   ! VelMean - Mean value of ice velocity (m/s)

   CALL ReadVar ( UnIn, FileName, InputFileData%miuv, 'VelMean', 'Mean value of ice velocity (m/s)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

   ! VelVar - Variance of ice velocity (m^2/s^2)

   CALL ReadVar ( UnIn, FileName, InputFileData%varv, 'VelVar', 'Variance of ice velocity (m^2/s^2)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

   ! TeMean - Mean value of ice loading event duration (s)

   CALL ReadVar ( UnIn, FileName, InputFileData%miut, 'TeMean', 'Mean value of ice loading event duration (s)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

    ! Header

   CALL ReadCom( UnIn, FileName, 'Ice PROPERTIES -Ice Model 3, SubModel 2,3', ErrStat, ErrMsg, UnEc)

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

   ! StrMean - Mean value of ice thickness (m)

   CALL ReadVar ( UnIn, FileName, InputFileData%miubr, 'StrMean', 'Mean value of ice strength (MPa)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

   ! StrVar - Variance of ice thickness (m^2)

   CALL ReadVar ( UnIn, FileName, InputFileData%varbr, 'StrVar', 'Variance of ice strength (MPa)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

   ! Header

   CALL ReadCom( UnIn, FileName, 'Ice PROPERTIES -Ice Model 3, SubModel 3', ErrStat, ErrMsg, UnEc)

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

   ! DelMean - Mean value of maximum ice tooth tip displacement (m)

   CALL ReadVar ( UnIn, FileName, InputFileData%miuDelm, 'DelMean', 'Mean value of maximum ice tooth tip displacement (m)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

   ! DelVar - Variance of maximum ice tooth tip displacement (m^2)

   CALL ReadVar ( UnIn, FileName, InputFileData%varDelm, 'DelVar', 'Variance of maximum ice tooth tip displacement (m^2)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

   ! PMean - Mean value of the distance between sequential ice teeth (m)

   CALL ReadVar ( UnIn, FileName, InputFileData%miuP, 'PMean', 'Mean value of the distance between sequential ice teeth (m)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

   ! PVar - Variance of the distance between sequential ice teeth (m^2)

   CALL ReadVar ( UnIn, FileName, InputFileData%varP, 'PVar', 'Variance of the distance between sequential ice teeth (m^2)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CALL Cleanup()
      RETURN
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Ice properties - Ice Model 4
   !-------------------------------------------------------------------------------------------------
      
      ! Header
      
   CALL ReadCom( UnIn, FileName, 'Ice properties - Ice Model 4', ErrStat, ErrMsg, UnEc)
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

     ! PrflMean - Mean value of ice contact face position (m)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%PrflMean, 'PrflMean', 'Mean value of ice contact face position (m)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
      ! PrflSig - Standard deviation of ice contact face position (m)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%PrflSig, 'PrflSig', 'Standard deviation of ice contact face position (m)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
      ! ZoneNo1 - Number of failure zones along contact width 
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Zn1, 'ZoneNo1', 'Number of failure zones along contact width', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
      
     ! ZoneNo2 - Number of failure zones along contact height/thickness
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Zn2, 'ZoneNo2', 'Number of failure zones along contact height/thickness', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
      ! ZonePitch - Distance between sequential ice teeth (m)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%ZonePitch, 'ZonePitch', 'Distance between sequential ice teeth (m)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   
      ! IceStr - Ice failure stress within each failure region (MPa)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%IceStr, 'IceStr', 'Ice failure stress within each failure region (MPa)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   ! Delmax - Ice teeth maximum elastic deformation (m)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Delmax, 'Delmax', 'Ice teeth maximum elastic deformation (m)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   !-------------------------------------------------------------------------------------------------
   ! Ice properties - Ice Model 5
   !-------------------------------------------------------------------------------------------------
   
   ! Header
      
   CALL ReadCom( UnIn, FileName, 'Ice properties - Ice Model 5, Submodel 1,2', ErrStat, ErrMsg, UnEc)
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   
     ! ConeAgl - Slope angle of the cone (degree)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%alpha, 'ConeAgl', 'Slope angle of the cone (degree)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   
     ! ConeDwl - Cone waterline diameter (m)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Dwl, 'ConeDwl', 'Cone waterline diameter (m)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   
     ! ConeDtp - Cone top diameter (m)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Dtp, 'ConeDtp', 'Cone top diameter (m)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   
     ! RdupThk - Ride-up ice thickness (m)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%hr, 'RdupThk', 'Ride-up ice thickness (m)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   
     ! mu - Friction coefficient between structure and ice (-)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%mu, 'mu', 'Friction coefficient between structure and ice (-)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   
     ! FlxStr - Flexural strength of ice (MPa)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%sigf, 'FlxStr', 'Flexural strength of ice (MPa)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   ! StrLim - Limit strain for ice fracture failure (-)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%StrLim, 'StrLim', 'Limit strain for ice fracture failure (-)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   ! StrRtLim - Limit strain rate for ice brittle behavior (s^-1)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%StrRtLim, 'StrRtLim', 'Limit strain rate for ice brittle behavior (s^-1)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   !-------------------------------------------------------------------------------------------------
   ! Ice properties - Ice Model 6
   !-------------------------------------------------------------------------------------------------
   
   ! Header
      
   CALL ReadCom( UnIn, FileName, 'Ice properties - Ice Model 6', ErrStat, ErrMsg, UnEc)
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
     ! FloeLth - Ice floe length (m)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Ll, 'FloeLth', 'Ice floe length (m)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
 
   ! FloeWth - Ice floe width (m)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Lw, 'FloeWth', 'Ice floe width (m)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
 
   ! CPrAr - Ice crushing strength pressure-area relation constant
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Cpa, 'CPrAr', 'Ice crushing strength pressure-area relation constant', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
 
   ! dPrAr - Ice crushing strength pressure-area relation order
      
   CALL ReadVar ( UnIn, FileName, InputFileData%dpa, 'dPrAr', 'Ice crushing strength pressure-area relation order', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
 
   ! Fdr - Constant external driving force (MN)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Fdr, 'Fdr', 'Constant external driving force (MN)', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
 
   ! Kic - Fracture toughness of ice (kNm^(-3/2))
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Kic, 'Kic', 'Fracture toughness of ice (kNm^(-3/2))', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
 
   ! FspN - Non-dimensional splitting load
         
   CALL ReadVar ( UnIn, FileName, InputFileData%FspN, 'FspN', 'Non-dimensional splitting load', ErrStat, ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
 

   !-------------------------------------------------------------------------------------------------
   ! This is the end of the input file
   !-------------------------------------------------------------------------------------------------
   CALL Cleanup()


   RETURN
CONTAINS
   SUBROUTINE Cleanup()
   
      CLOSE( UnIn )
      IF (UnEc > 0) CLOSE(UnEc)
      
   END SUBROUTINE

END SUBROUTINE IceD_ReadInput
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine performs validity checks on data from the input file.
SUBROUTINE IceD_ValidateInput( InputFileData, ErrStat, ErrMsg )
   TYPE(IceD_InputFile),     INTENT(IN)     :: InputFileData                !< Data stored in the module's input file
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                      !< Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                       !< Error message

   
   ErrStat = ErrID_None
   ErrMsg  = ''

   IF ( InputFileData%IceModel < 1 .OR. InputFileData%IceModel > 6 ) CALL SetErrStat( ErrID_Fatal, 'IceModel must be a number 1-6.',   ErrStat, ErrMsg, 'IceD_ValidateInput')
   IF ( InputFileData%IceSubModel < 1  )                             CALL SetErrStat( ErrID_Fatal, 'Invalid IceSubModel value', ErrStat, ErrMsg, 'IceD_ValidateInput')
   
   
   
END SUBROUTINE IceD_ValidateInput
!----------------------------------------------------------------------------------------------------------------------------------
!> This takes the primary input file data and sets the corresponding parameters.
SUBROUTINE IceD_SetParameters( InputFileData, p, Interval, Tmax, LegNum, ErrStat, ErrMsg  )
!..................................................................................................................................

   IMPLICIT                        NONE


      ! Passed variables

   TYPE(IceD_ParameterType), INTENT(INOUT)  :: p                            !< Parameters of the IceDyn module
   TYPE(IceD_InputFile),     INTENT(IN)     :: InputFileData                !< Data stored in the module's input file
   REAL(DbKi),               INTENT(IN)     :: Interval                     !< Coupling interval in seconds
   REAL(DbKi),               INTENT(IN   )  :: Tmax                         !< Total simulation time 
   INTEGER(IntKi),           INTENT(IN   )  :: LegNum                       !< Leg number of this IceDyn instance
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                      !< Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                       !< Error message

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
!   REAL(ReKi), allocatable                  :: rdmTstr(:)                   ! Random ice strength of ice teeth (Pa)
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
   
    ! Ice Model 4
   REAL(ReKi)                               :: ZoneWd                       ! Width of a single failure zone
   REAL(ReKi)                               :: ZoneHt                       ! Height of a single failuer zone
 
    ! Ice Model 5
   REAL(ReKi)                       :: flexStrPa               ! Ice flexural strength (Pa)  
   REAL(ReKi)                       :: A(6)                    ! Coefficients when calculating ice breaking force using sub-model 1
   REAL(ReKi)                       :: Pn1
   REAL(ReKi)                       :: Pn2
   REAL(ReKi)                       :: F1
   REAL(ReKi)                       :: Lxlim1
   REAL(ReKi)                       :: Lxlim2
   REAL(ReKi)                               :: Pbr


    ! Initialize error data
   ErrStat = ErrID_None
   ErrMsg  = ''

   
   p%verif     = 1
   p%method    = 1            ! integration method (RK4)
   p%dt        = Interval
   p%tolerance = 1e-6      
   
   
   !...............................................................................................................................
   ! Direct copy of variables:
   !...............................................................................................................................
   p%ModNo     = InputFileData%IceModel
   p%SubModNo  = InputFileData%IceSubModel
   p%v         = InputFileData%v
   p%h         = InputFileData%h
   p%InitLoc   = InputFileData%InitLoc
   p%t0        = InputFileData%t0
   p%StrWd     = InputFileData%StrWd( LegNum )
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
   
     ! Ice Model 4
   p%Delmax    = InputFileData%Delmax
   p%ZonePitch = InputFileData%ZonePitch
   
   ! Ice Model 5
   p%rhoi      = InputFileData%rhoi
   p%rhow      = InputFileData%rhow
   p%Dwl    = InputFileData%Dwl
   p%mu        = InputFileData%mu
      
   ! Ice Model 6
   p%Cpa    = InputFileData%Cpa
   p%dpa    = InputFileData%dpa

      
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
   kappa       = ( InputFileData%rhow * 9.81 / 4.0 / Bf ) ** 0.25 !bjj: can you use the gravitational constant defined in FAST? now stored in InitInput%gravity
   PhiR        = InputFileData%phi / 180.0 * 3.1415927  !bjj todo: you can use "D2R" from NWTC_Library to convert degrees to radians
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
              CALL IceD_Generate_RandomNum ( rdmh(J), rdmv(J), rdmte(J), s_dmmy, Dm_dmmy, P_dmmy, p, InputFileData, ErrStat, ErrMsg)
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
              CALL IceD_Generate_RandomNum ( rdmh(J), rdmv(J), rdmte(J), s_dmmy, Dm_dmmy, P_dmmy, p, InputFileData, ErrStat, ErrMsg)
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
              CALL IceD_Generate_RandomNum ( rdmh(J), rdmv(J), rdmte(J), rdmsig(J), Dm_dmmy, P_dmmy, p, InputFileData, ErrStat, ErrMsg)
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
              CALL IceD_Generate_RandomNum ( rdmh(J), rdmv(J), rdmte(J), rdmsig(J), Dm_dmmy, P_dmmy, p, InputFileData, ErrStat, ErrMsg)
              p%rdmFm(J)    = p%StrWd * rdmh(J) * rdmsig(J)
              p%rdmtm(J)    = rdmsig(J) / p%EiPa / ( rdmv(J) / 4.0 / p%StrWd )              
          ENDIF
              CrntT0 = p%rdmt0 (J) !bjj: doesn't appear to be used...
       END DO
       
    ELSEIF (p%SubModNo == 3) THEN
    
        DO J = 1, Nthmax        
            CALL IceD_Generate_RandomNum ( h_dmmy, v_dmmy, t_dmmy, rdmsig(J), p%rdmDm(J), p%rdmP(J), p, InputFileData, ErrStat, ErrMsg)
            p%rdmKi(J) = rdmsig(J) * p%StrWd * p%h / p%rdmDm(J)
        END DO
        
    ENDIF

    ! Ice Model 4
   ZoneWd      = p%StrWd / REAL(InputFileData%Zn1)
   ZoneHt      = InputFileData%h     / REAL(InputFileData%Zn2)
   p%Kice      = InputFileData%IceStr *1e6 * ZoneWd * ZoneHt / InputFileData%Delmax
   
   p%Zn        = InputFileData%Zn1   * InputFileData%Zn2                           ! Total number of failure zones
   
   CALL AllocAry( p%ContPrfl, p%Zn, 'ContPrfl', ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN
   
    CALL AllocAry( p%Y0, p%Zn, 'ContPrfl', ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN
   
   DO I = 1,p%Zn
      p%ContPrfl (I) = InputFileData%PrflMean + InputFileData%PrflSig * random_normal()
      !CALL WrScr ( Num2LStr( p%ContPrfl (I) ) )
   END DO
   
   p%Y0        = InputFileData%InitLoc + p%ContPrfl - MAXVAL(p%ContPrfl)
   
   ! Ice Model 5
   flexStrPa   = InputFileData%sigf * 1e6
   
   p%alphaR    = InputFileData%alpha / 180.0 * 3.1415927
   p%Zr        = (InputFileData%Dwl - InputFileData%Dtp) / 2 * tan(p%alphaR)
   p%Wri       = p%rhoi * 9.81 * p%Dwl * p%h * p%Zr / sin(p%alphaR)
   
   IF (p%SubModNo == 1) THEN
         
         p%LovR = SolveLambda ( p%rhoi, p%h, p%Dwl, flexStrPa )
         p%Lbr  = p%LovR * p%Dwl / 2
         A     = BrkLdPar (p%alphaR, p%LovR, InputFileData%mu)
         
         p%RHbr = ( A(1) * flexStrPa * p%h**2 + A(2) * p%rhoi * 9.81 * p%h * p%Dwl**2 + A(3) * p%rhoi * 9.81 * p%h * (p%Dwl**2 - InputFileData%Dtp**2) ) * A(4)
         p%RVbr = A(5) * p%RHbr + A(6) * p%rhoi * 9.81 * p%h * (p%Dwl**2 - InputFileData%Dtp**2)
   
        
   ELSEIF (p%SubModNo == 2) THEN
   
         Pbr      = 8.0 * sqrt(2.0) * ( ( flexStrPa * p%h**2) / 4 )
         Pn1      = p%Wri * cos(p%alphaR)
         F1     = p%Wri * ( sin(p%alphaR) + p%mu * cos(p%alphaR) );
        Pn2    = ( Pbr + F1*sin(p%alphaR) ) / ( cos(p%alphaR) - p%mu * sin(p%alphaR) )
         
         p%RHbr = (Pn1 + Pn2) * ( sin(p%alphaR) + p%mu * cos(p%alphaR) )
         p%RVbr = (Pn1 + Pn2) * ( cos(p%alphaR) - p%mu * sin(p%alphaR) )
         
         Lxlim1 = ( 3.0 * sqrt(6.0) ) / 8.0 * ( p%v * tan(p%alphaR) ) / InputFileData%StrRtLim   !Limit strain rate criteria
         Lxlim2 = sqrt(6.0) *  ( ( flexStrPa * p%h**2) / 4 / (p%rhow * 9.81) / InputFileData%StrLim )**(1.0/3.0)
         
         IF (Lxlim1 <= Lxlim2) THEN
         
            p%Lbr = ( 3.0 * sqrt(2.0) ) / 8.0 * ( p%v * tan(p%alphaR) ) / InputFileData%StrRtLim
         
         ELSE
         
            p%Lbr = 2.0 *  ( ( flexStrPa * p%h**2) / 4 / (p%rhow * 9.81) / InputFileData%StrLim )**(1.0/3.0)
         
         ENDIF
   
   ELSE
         
         ErrMsg   = 'Sub-model number for model 5 should be 1 or 2'
         ErrStat = ErrID_Fatal
   
   ENDIF 
   
   p%WL        = p%rhoi * 9.81 * p%Dwl * p%h * p%Lbr
   
   ! Ice Model 6
   p%FdrN      = InputFileData%Fdr * 1e6
   p%Mice      = InputFileData%rhoi * p%h * InputFileData%Lw * InputFileData%Ll
   p%Fsp       = InputFileData%FspN * p%h * InputFileData%Kic * 1e3 * sqrt(InputFileData%Ll) 
    
    
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

   p%NumOuts    = 2

      CALL AllocAry( p%OutName, p%NumOuts, 'OutName', ErrStat, ErrMsg )
      IF (ErrStat >= AbortErrLev) RETURN

      !p%OutName (1) = 'Time'
      p%OutName (1) = 'IceDisp'
      p%OutName (2) = 'IceForce'

      CALL AllocAry( p%OutUnit, p%NumOuts, 'OutUnit', ErrStat, ErrMsg )
      IF (ErrStat >= AbortErrLev) RETURN

      !p%OutUnit (1) = '(s)'
      p%OutUnit (1) = '(m)'
      p%OutUnit (2) = '(N)'
      
      ! Test parameter assignments
      
      CONTAINS
      
      FUNCTION random_normal() RESULT(fn_val)
      
         ! Adapted from the following Fortran 77 code
         !      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
         !      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
         !      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.
         
         !  The function random_normal() returns a normally distributed pseudo-random
         !  number with zero mean and unit variance.
         
         !  The algorithm uses the ratio of uniforms method of A.J. Kinderman
         !  and J.F. Monahan augmented with quadratic bounding curves.
         
         REAL(ReKi) :: fn_val
         
         !     Local variables
         REAL(ReKi)     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
                     r1 = 0.27597, r2 = 0.27846, u, v, x, y, q
         
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
        
          
        FUNCTION  SolveLambda(rhoi, t, D, sigf) Result (rho)
      
      !  SOLVERHO Solve for rho according to Ralston model (Ralston 1978)
      !   Rho = A/R, A is the first circumferential crack radius, R is the cone
      !   structure waterline radius. According to first equation on Ralston 
      !   paper (P301), calculate rho.
         
         IMPLICIT NONE
      
         ! Input values
         REAL(ReKi) :: rhoi      ! Mass density of ice, (kg/m^3)
         REAL(ReKi) :: t         ! Ice thickness (m)
         REAL(ReKi) :: D         ! Cone waterline diameter (m)
         REAL(ReKi) :: sigf      ! Ice flextural strength (Pa)
         
         REAL(ReKi) :: rho       ! Rho = A/R
         REAL(ReKi) :: x = 1.01    ! Initial value of rho
         REAL(ReKi) :: Equ    
         REAL(ReKi) :: Derv
         
         DO i = 1,100
         
            Equ = x - log(x) + 0.0922 * rhoi * 9.81 * t * D**2 / sigf / t**2 * (2*x+1) * (x-1)**2 - 1.369
            Derv = 1 - 1/x + 0.0922 * rhoi * 9.81 * t * D**2 / sigf / t**2 * (2*(x-1)**2 + (2*x+1)*2*(x-1))
            
            IF ( abs(Equ) <= 1e-6) THEN
            
               rho = x
               EXIT
            
            END IF 
            
            x = x - Equ / Derv
         
         END DO
                
      END FUNCTION SolveLambda
      
      
      FUNCTION BrkLdPar (alpha, lambda, mu) Result (A)
      
      !BRKLDPAR Calculates Ralston's horizontal force paramters A1, A2, A3, A4 and B1, B2.
      !   Detailed explanation in Ralston's paper: Ice Force Desgin Consideration 
      !  for Conical Offshore Structure and Ice Module Manual
      
         IMPLICIT NONE
      
         ! Input values
         REAL(ReKi) :: alpha        ! Cone angle, (rad)
         REAL(ReKi) :: lambda       ! Ratio of breaking length over cone waterline radius
         REAL(ReKi) :: mu        ! Friction coefficient between structure and ice
         
         REAL(ReKi) :: A(6)         ! Coefficients when calculating ice breaking force
         
         ! Local variables
         REAL(ReKi) :: f
         REAL(ReKi) :: g
         REAL(ReKi) :: h
         REAL(ReKi) :: pi = 3.1415927
         
         A(1) = 1.0/3.0 * ( lambda/(lambda-1) + (1-lambda+lambda*log(lambda))/(lambda-1) + 2.422* (lambda*log(lambda))/(lambda-1) )

         A(2) = ( lambda**2 + lambda -2.0 )/12.0

         f = pi/2.0 + pi/8.0 * (sin(alpha))**2 / (1-(sin(alpha))**2) - pi/16.0 * (sin(alpha))**4 / (1-(sin(alpha))**4)
         g = ( 1.0/2.0 + alpha/sin(2*alpha) ) / ( pi/4.0*sin(alpha) + mu*alpha*cos(alpha)/sin(alpha) )
         A(3) = 1.0/4.0 * ( 1/cos(alpha) + mu*Esina(alpha,10)/sin(alpha) - mu*f*g/tan(alpha) )

         A(4) = tan(alpha) / ( 1 - mu * g)

         h = cos(alpha) - mu/sin(alpha) * ( Esina(alpha,10) - cos(alpha)**2 * Fsina(alpha) )

         A(5) = h / ( pi/4.0 * sin(alpha) + mu * alpha / tan(alpha) )
         A(6) = 1.0/4.0 * (pi/2.0*cos(alpha) - mu*alpha - f*h/ ( pi/4.0 * sin(alpha) + mu * alpha / tan(alpha) ))
            
            !CALL WrScr(Num2LStr(Esina(alpha,10)))
            !CALL WrScr(Num2LStr(Fsina(alpha)))
               
      END FUNCTION BrkLdPar
      
      
      FUNCTION Esina (alpha, n) Result (Esin)
      !ESINA calculates E(sin(alpha)). Detailed explanation in Ice Module Manual, Model 5, Submodel 1
      
         IMPLICIT NONE
         
         !Input variable
         REAL(ReKi)     :: alpha       ! Cone angle, (rad)
         INTEGER(IntKi) :: n     
         
         !Output
         REAL(ReKi)     :: Esin 
            
         
         !Local variable 
         INTEGER(IntKi) :: i  
         REAL(ReKi)     :: E = 0
            REAL(ReKi)     :: pi = 3.1415927
         
         DO i = 1,n
         
            E = E + pi/2.0*( factorial(2*(i-1)) / 2**(2*(i-1)) / (factorial(i-1))**2 )**2 * (sin(alpha))**(2*(i-1)) / (1-2*(i-1))
         
            END DO
            
            Esin = E
          
      END FUNCTION Esina
      
      
      FUNCTION Fsina (alpha) Result (F)
      !ESINA calculates F(sin(alpha)). Detailed explanation in Ice Module Manual, Model 5, Submodel 1
      
         IMPLICIT NONE
         
         !Input variable
         REAL(ReKi)     :: alpha       ! Cone angle, (rad)
         !Output
         REAL(ReKi)     :: F
         !Local variable 
         REAL(ReKi)     :: pi = 3.1415927
         
         F = pi/2.0 + pi/8.0 * sin(alpha)**2 / (1-sin(alpha)**2) - pi/16.0 * sin(alpha)**4 / (1-sin(alpha)**4)
      
      END FUNCTION Fsina
      

        FUNCTION factorial (n) Result (fac)
        ! FACTORIAL calculates the factorial of n 
        
         IMPLICIT NONE
         
         !Input variable
         INTEGER(IntKi) :: n
         
         !Output
         REAL(ReKi)     :: fac
         
         !Local variables
         INTEGER(IntKi) :: i  
         INTEGER(IntKi) :: M = 1
         
         IF (n == 0) THEN
             
             fac = 1.0
             
         ELSE
         
             DO i = 1,n
         
                M = M * i  
         
             ENDDO
             
         ENDIF
         
         fac = REAL(M)
         M = 1
        
        END FUNCTION factorial
        
END SUBROUTINE IceD_SetParameters
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine initializes the other states of the module.
!! It assumes the parameters are set and that InputFileData contains initial conditions for the other states.
SUBROUTINE IceD_Init_OtherStates( OtherState, p, x, InputFileData, ErrStat, ErrMsg  )
!..................................................................................................................................
   IMPLICIT                        NONE

   TYPE(IceD_OtherStateType),      INTENT(  OUT)  :: OtherState        !< Initial other states
   TYPE(IceD_ParameterType),       INTENT(IN)     :: p                 !< Parameters of the IceDyn module
   TYPE(IceD_ContinuousStateType), INTENT(IN   )  :: x                 !< Initial continuous states
   TYPE(IceD_InputFile),           INTENT(IN)     :: InputFileData     !< Data stored in the module's input file
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat           !< Error status
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg            !< Error message

      ! local variables
   INTEGER(IntKi)                                 :: I                 ! loop counter
!   REAL(ReKi)                                     :: StrRt             ! Strain rate (s^-1)
!   REAL(ReKi)                                     :: SigCrp            ! Creep stress (Pa)
      ! Initialize error data
   ErrStat = ErrID_None
   ErrMsg  = ''

   IF  ( p%ModNo == 2 ) THEN

       OtherState%IceTthNo2 = 1 ! Initialize first ice tooth number

   ELSEIF ( p%ModNo == 3 ) THEN

       CALL AllocAry( OtherState%Nc, p%TmStep, 'OtherState%Nc', ErrStat, ErrMsg )
       IF ( ErrStat >= AbortErrLev ) RETURN
   
       CALL AllocAry( OtherState%Psum, p%TmStep, 'OtherState%Psum', ErrStat, ErrMsg )
       IF ( ErrStat >= AbortErrLev ) RETURN
       
       DO I = 1,p%TmStep
          OtherState%Nc (I)      = 1. ! Initialize first ice tooth number
          OtherState%Psum (I)    = 0. ! Initialize sum of pitches of broken ice teeth 
       ENDDO

   ELSEIF ( p%ModNo == 4 ) THEN
       
       CALL AllocAry( OtherState%IceTthNo, p%Zn,   'IceTthNo',   ErrStat, ErrMsg )
       IF ( ErrStat /= ErrID_None ) RETURN ! Initialize first ice tooth number for each zone
       DO I = 1,p%Zn
          OtherState%IceTthNo (I) = 1
       END DO
    
    ELSEIF ( p%ModNo == 5 ) THEN
    
       OtherState%beta = 0.    ! ice plate lifted angle
       OtherState%Tinit = p%t0
       
    ELSEIF ( p%ModNo == 6 ) THEN
    
       OtherState%dxc = 0.    ! ice crushed depth
       OtherState%Splitf = 0. ! flag to indicate if the ice floe has splitted (0 not splitted, 1 splitted)
        
   ENDIF

   
   if ( p%method .eq. 2) then  !integration methods

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
   
   
   
END SUBROUTINE IceD_Init_OtherStates
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine generate random numbers for the module.
!! It assumes the parameters are set and that InputFileData contains input for generating random numbers.
SUBROUTINE IceD_Generate_RandomNum ( h, v, t, s, Dm, Pch, p, InputFileData, ErrStat, ErrMsg)
!..................................................................................................................................
   IMPLICIT                        NONE

   REAL(ReKi),                   INTENT(OUT)    :: h                 !< Random ice thickness (m)
   REAL(ReKi),                   INTENT(OUT)    :: v                 !< Random ice velocity (m/s)
   REAL(ReKi),                   INTENT(OUT)    :: t                 !< Random ice loading event time (s)
   REAL(ReKi),                   INTENT(OUT)    :: s                 !< Random ice strength (Pa)
   REAL(ReKi),                   INTENT(OUT)    :: Dm                !< Random ice tooth maximum displacement (m)
   REAL(ReKi),                   INTENT(OUT)    :: Pch               !< Random ice tooth spacing (m)
   TYPE(IceD_ParameterType),     INTENT(IN)     :: p                 !< Parameters of the IceDyn module
   TYPE(IceD_InputFile),         INTENT(IN)     :: InputFileData     !< Data stored in the module's input file
   INTEGER(IntKi),               INTENT(OUT)    :: ErrStat           !< Error status
   CHARACTER(*),                 INTENT(OUT)    :: ErrMsg            !< Error message   

      ! local variables
!   INTEGER(IntKi)                               :: I                 ! loop counter
   REAL(ReKi)                                   :: SigLogh           ! sigma_log(h), standard deviation of log(h)
   REAL(ReKi)                                   :: MiuLogh           ! miu_log(h), mean value of log(h)
   REAL(ReKi)                                   :: VelSig            ! parameter for a Rayleigh distribution
   REAL(ReKi)                                   :: TeLamb            ! parameter for a exponential distribution
   !REAL, PARAMETER                              :: Pi = 3.1415927 !bjj: comes from NWTC_Library

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

         REAL(ReKi)  :: fn_val
         REAL(ReKi)  :: Lambda

         !     Local variable
         REAL(ReKi)  :: r

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

         REAL(ReKi)  :: fn_val
         REAL(ReKi)  :: Sigma

         !     Local variable
         REAL(ReKi)  :: r

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

         REAL(ReKi) :: mean
         REAL(ReKi) :: var
         REAL(ReKi) :: fn_val

         !Local variables                                                                                                                                                                     
         REAL(ReKi) :: k
         REAL(ReKi) :: Lambda
         REAL(ReKi) :: u

         k = wbpar (mean, var)
         lambda = mean / NWTC_gamma(1+1/k)
       ! lambda = mean /      gamma(1+1/k)

         CALL RANDOM_NUMBER(u)
         fn_val = lambda * (-log(1-u)) ** (1/k)

      END FUNCTION random_weibull


      FUNCTION wbpar (mean, var) Result (k1)

         !Calculate Weibull distribution parameters due to mean value and variance of the data
         IMPLICIT NONE
                                                                                                                                                                                       
         REAL(ReKi) :: mean
         REAL(ReKi) :: var
         REAL(ReKi) :: k 
         REAL(ReKi) :: k1, F1, dFdk
         REAL(ReKi) :: error 

         INTEGER :: I
         
         k = 10
         error = 1e-6

         DO i = 1,10000

            F1 = (NWTC_gamma(1+1/k))**2 / NWTC_gamma(1+2/k) - mean**2/(mean**2+var);
           !F1 = (     gamma(1+1/k))**2 /      gamma(1+2/k) - mean**2/(mean**2+var)

            IF (abs(F1) < error) EXIT

            dFdk = 2* (NWTC_gamma(1+1/k))**2 * (-1/k**2) / NWTC_gamma(1+2/k) * (digamma(1+1/k) -digamma(1+2/k));
           !dFdk = 2* (     gamma(1+1/k))**2 * (-1/k**2) /      gamma(1+2/k) * (digamma(1+1/k) -digamma(1+2/k))
            k = k - F1/dFdk

            END DO

            !IF (abs(F1) >= error)THEN

            !    WrScr('Weibull parameters never found')

            !ENDIF


          k1 = k

        END FUNCTION wbpar

      FUNCTION digamma(z) RESULT(phy)

         !Calculate the value of digamma function of z
         REAL(ReKi), INTENT(IN) :: z
         REAL(ReKi)             :: phy

         phy = log(z) - 1./2./z - 1./12./z**2 + 1./120./z**4 - 1./252./z**6 + 1./240./z**8 - 5./660./z**10;

      END FUNCTION digamma

END SUBROUTINE IceD_Generate_RandomNum
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
!!   Runge-Kutta." §16.1 and 16.2 in Numerical Recipes in FORTRAN: The Art of Scientific Computing, 2nd ed. Cambridge, England:
!!   Cambridge University Press, pp. 704-716, 1992.
SUBROUTINE IceD_RK4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           !< time step number
      TYPE(IceD_InputType),           INTENT(INOUT)  :: u(:)        !< Inputs at t
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   !< times of input
      TYPE(IceD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(IceD_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
      TYPE(IceD_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(IceD_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t
      TYPE(IceD_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states at t
      TYPE(IceD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! local variables

      TYPE(IceD_ContinuousStateType)                 :: xdot        ! time derivatives of continuous states
      TYPE(IceD_ContinuousStateType)                 :: k1          ! RK4 constant; see above
      TYPE(IceD_ContinuousStateType)                 :: k2          ! RK4 constant; see above
      TYPE(IceD_ContinuousStateType)                 :: k3          ! RK4 constant; see above
      TYPE(IceD_ContinuousStateType)                 :: k4          ! RK4 constant; see above
      TYPE(IceD_ContinuousStateType)                 :: x_tmp       ! Holds temporary modification to x
      TYPE(IceD_InputType)                           :: u_interp    ! interpolated value of inputs

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL IceD_CopyInput( u(1), u_interp, MESH_NEWCOPY, ErrStat, ErrMsg) ! bjj: need to allocate space for u_interp so that IceD_Input_ExtrapInterp works

      ! interpolate u to find u_interp = u(t)
      CALL IceD_Input_ExtrapInterp( u, utimes, u_interp, t, ErrStat, ErrMsg )

      ! find xdot at t
      CALL IceD_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, m, xdot, ErrStat, ErrMsg )

      k1%q    = p%dt * xdot%q
      k1%dqdt = p%dt * xdot%dqdt

      x_tmp%q    = x%q    + 0.5 * k1%q
      x_tmp%dqdt = x%dqdt + 0.5 * k1%dqdt

      ! interpolate u to find u_interp = u(t + dt/2)
      CALL IceD_Input_ExtrapInterp(u, utimes, u_interp, t+0.5*p%dt, ErrStat, ErrMsg)

      ! find xdot at t + dt/2
      CALL IceD_CalcContStateDeriv( t + 0.5*p%dt, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, ErrStat, ErrMsg )

      k2%q    = p%dt * xdot%q
      k2%dqdt = p%dt * xdot%dqdt

      x_tmp%q    = x%q    + 0.5 * k2%q
      x_tmp%dqdt = x%dqdt + 0.5 * k2%dqdt

      ! find xdot at t + dt/2
      CALL IceD_CalcContStateDeriv( t + 0.5*p%dt, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, ErrStat, ErrMsg )

      k3%q    = p%dt * xdot%q
      k3%dqdt = p%dt * xdot%dqdt

      x_tmp%q    = x%q    + k3%q
      x_tmp%dqdt = x%dqdt + k3%dqdt

      ! interpolate u to find u_interp = u(t + dt)
      CALL IceD_Input_ExtrapInterp(u, utimes, u_interp, t + p%dt, ErrStat, ErrMsg)

      ! find xdot at t + dt
      CALL IceD_CalcContStateDeriv( t + p%dt, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, ErrStat, ErrMsg )

      k4%q    = p%dt * xdot%q
      k4%dqdt = p%dt * xdot%dqdt

      x%q    = x%q    +  ( k1%q    + 2. * k2%q    + 2. * k3%q    + k4%q    ) / 6.
      x%dqdt = x%dqdt +  ( k1%dqdt + 2. * k2%dqdt + 2. * k3%dqdt + k4%dqdt ) / 6.
      
      CALL Cleanup()

CONTAINS
SUBROUTINE Cleanup()
   
   CALL IceD_DestroyInput( u_interp, ErrStat, ErrMsg) 
   CALL IceD_DestroyContState( xdot , ErrStat, ErrMsg) 
   CALL IceD_DestroyContState( k1   , ErrStat, ErrMsg) 
   CALL IceD_DestroyContState( k2   , ErrStat, ErrMsg) 
   CALL IceD_DestroyContState( k3   , ErrStat, ErrMsg) 
   CALL IceD_DestroyContState( k4   , ErrStat, ErrMsg) 
   CALL IceD_DestroyContState( x_tmp, ErrStat, ErrMsg) 
                               
END SUBROUTINE Cleanup
END SUBROUTINE IceD_RK4
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
SUBROUTINE IceD_AB4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           !< time step number
      TYPE(IceD_InputType),           INTENT(INOUT)  :: u(:)        !< Inputs at t
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   !< times of input
      TYPE(IceD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(IceD_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
      TYPE(IceD_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(IceD_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t
      TYPE(IceD_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states at t
      TYPE(IceD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


      ! local variables
      TYPE(IceD_ContinuousStateType)                 :: xdot       ! Continuous state derivs at t
      TYPE(IceD_InputType)                           :: u_interp


      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""

      
      CALL IceD_CopyInput( u(1), u_interp, MESH_NEWCOPY, ErrStat, ErrMsg) ! bjj: need to allocate space for u_interp so that IceD_Input_ExtrapInterp works
      
      ! need xdot at t
      CALL IceD_Input_ExtrapInterp(u, utimes, u_interp, t, ErrStat, ErrMsg)
      CALL IceD_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, m, xdot, ErrStat, ErrMsg )

      if (n .le. 2) then

         OtherState%n = n

         OtherState%xdot ( 3 - n ) = xdot

         CALL IceD_RK4(t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )

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

      CALL Cleanup()

CONTAINS
SUBROUTINE Cleanup()
   
   CALL IceD_DestroyInput( u_interp, ErrStat, ErrMsg) 
   CALL IceD_DestroyContState( xdot , ErrStat, ErrMsg) 
                               
END SUBROUTINE Cleanup
      
      
END SUBROUTINE IceD_AB4
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine implements the fourth-order Adams-Bashforth-Moulton Method (RK4) for numerically integrating ordinary
!! differential equations:
!!
!!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x).
!!
!!   Adams-Bashforth Predictor:
!!   x^p(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
!!
!!   Adams-Moulton Corrector:
!!   x(t+dt) = x(t)  + (dt / 24.) * ( 9.*f(t+dt,x^p) + 19.*f(t,x) - 5.*f(t-dt,x) + 1.*f(t-2.*dt,x) )
!!
!!  See, e.g.,
!!  http://en.wikipedia.org/wiki/Linear_multistep_method
!!
!!  or
!!
!!  K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
SUBROUTINE IceD_ABM4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           !< time step number
      TYPE(IceD_InputType),           INTENT(INOUT)  :: u(:)        !< Inputs at t
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   !< times of input
      TYPE(IceD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(IceD_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
      TYPE(IceD_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(IceD_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t
      TYPE(IceD_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states at t
      TYPE(IceD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! local variables

      TYPE(IceD_InputType)                           :: u_interp        ! Continuous states at t
      TYPE(IceD_ContinuousStateType)                 :: x_pred          ! Continuous states at t
      TYPE(IceD_ContinuousStateType)                 :: xdot_pred       ! Continuous states at t

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL IceD_CopyContState(x, x_pred, MESH_NEWCOPY, ErrStat, ErrMsg)

      CALL IceD_AB4( t, n, u, utimes, p, x_pred, xd, z, OtherState, m, ErrStat, ErrMsg )

      if (n .gt. 2) then

         CALL IceD_CopyInput( u(1), u_interp, MESH_NEWCOPY, ErrStat, ErrMsg) ! bjj: need to allocate space for u_interp so that IceD_Input_ExtrapInterp works
         CALL IceD_Input_ExtrapInterp(u, utimes, u_interp, t + p%dt, ErrStat, ErrMsg)

         CALL IceD_CalcContStateDeriv(t + p%dt, u_interp, p, x_pred, xd, z, OtherState, m, xdot_pred, ErrStat, ErrMsg )

         x%q    = x%q    + (p%dt / 24.) * ( 9. * xdot_pred%q +  19. * OtherState%xdot(1)%q - 5. * OtherState%xdot(2)%q &
                                          + 1. * OtherState%xdot(3)%q )

         x%dqdt = x%dqdt + (p%dt / 24.) * ( 9. * xdot_pred%dqdt + 19. * OtherState%xdot(1)%dqdt - 5. * OtherState%xdot(2)%dqdt &
                                          + 1. * OtherState%xdot(3)%dqdt )

      else

         x%q    = x_pred%q
         x%dqdt = x_pred%dqdt

      endif

      CALL Cleanup()

CONTAINS
SUBROUTINE Cleanup()
   
   CALL IceD_DestroyInput( u_interp, ErrStat, ErrMsg) 
   CALL IceD_DestroyContState( x_pred , ErrStat, ErrMsg) 
   CALL IceD_DestroyContState( xdot_pred, ErrStat, ErrMsg) 
                               
END SUBROUTINE Cleanup
      
END SUBROUTINE IceD_ABM4
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE IceDyn
!**********************************************************************************************************************************
