!>
!!
!! Abbreviations:
!!   - FVW: Free Vortex Wake
!!   - LL : Lifting Line
!!   - CP : Control point
!!   - NW : Near Wake
!!   - FW : Far Wake
!!
MODULE FVW
   USE NWTC_Library
   USE FVW_Types
   USE FVW_Subs
   use FVW_IO
   use FVW_Wings
   use FVW_BiotSavart
   use FVW_Tests

      ! NOTE:  this is a rough format that AD14 stores airfoil info.  This will need
      !        to be replaced by the AirFoilInfo module when we couple FVW into AD15
   USE AD14AeroConf_Types


   IMPLICIT NONE

   PRIVATE

   type(ProgDesc), PARAMETER  :: FVW_Ver = ProgDesc( 'FVW', '', '' )

   PUBLIC   :: FVW_Init             ! Initialization routine
   PUBLIC   :: FVW_End

   PUBLIC   :: FVW_CalcOutput
   PUBLIC   :: FVW_UpdateStates


CONTAINS

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the start of the simulation to perform initialization steps.
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
subroutine FVW_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )
!..................................................................................................................................
   type(FVW_InitInputType),         intent(inout)  :: InitInp        !< Input data for initialization routine  (inout so we can use MOVE_ALLOC)
   type(FVW_InputType),             intent(  out)  :: u              !< An initial guess for the input; input mesh must be defined
   type(FVW_ParameterType),         intent(  out)  :: p              !< Parameters
   type(FVW_ContinuousStateType),   intent(  out)  :: x              !< Initial continuous states
   type(FVW_DiscreteStateType),     intent(  out)  :: xd             !< Initial discrete states
   type(FVW_ConstraintStateType),   intent(  out)  :: z              !< Initial guess of the constraint states
   type(FVW_OtherStateType),        intent(  out)  :: OtherState     !< Initial other states
   type(FVW_OutputType),            intent(  out)  :: y              !< Initial system outputs (outputs are not calculated;
                                                                     !!   only the output mesh is initialized)
   type(FVW_MiscVarType),           intent(  out)  :: m              !< Initial misc/optimization variables
   real(DbKi),                      intent(inout)  :: interval       !< Coupling interval in seconds: the rate that
                                                                     !!   (1) FVW_UpdateStates() is called in loose coupling &
                                                                     !!   (2) FVW_UpdateDiscState() is called in tight coupling.
                                                                     !!   Input is the suggested time from the glue code;
                                                                     !!   Output is the actual coupling interval that will be used
                                                                     !!   by the glue code.
   type(FVW_InitOutputType),        intent(  out)  :: InitOut        !< Output for initialization routine
   integer(IntKi),                  intent(  out)  :: ErrStat        !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None
   ! Local variables
   integer(IntKi)          :: ErrStat2       ! temporary error status of the operation
   character(ErrMsgLen)    :: ErrMsg2        ! temporary error message
   integer(IntKi)          :: UnEcho         ! Unit number for the echo file
   character(*), parameter :: RoutineName = 'FVW_Init'
   type(FVW_InputFile)     :: InputFileData                                            !< Data stored in the module's input file

   ! Initialize variables for this routine
   ErrStat = ErrID_None
   ErrMsg  = ""
   UnEcho  = -1

   ! Initialize the NWTC Subroutine Library
   call NWTC_Init( EchoLibVer=.FALSE. )

   ! Display the module information
   call DispNVD( FVW_Ver )

   ! Set Parameters and *Misc* from inputs
   CALL FVW_SetParametersFromInputs(InitInp, p, m, ErrStat2, ErrMsg2); if(Failed()) return

   ! Read and parse the input file here to get other parameters and info
   CALL FVW_ReadInputFile(InitInp%FVWFileName, p, InputFileData, ErrStat2, ErrMsg2); if(Failed()) return

   ! Initialize Misc Vars (may depend on input file)
   p%nNWMax=100 ! TODO 
   p%nFWMax=100 ! TODO
   CALL FVW_InitMiscVars( p, m, ErrStat2, ErrMsg2 ); if(Failed()) return

   ! Move the InitInp%WingsMesh to u
   CALL MOVE_ALLOC( InitInp%WingsMesh, u%WingsMesh )     ! Move from InitInp to u

   ! This mesh is passed in as a cousin of the BladeMotion mesh.
   CALL Wings_Panelling_Init(u%WingsMesh, InitInp%zLocal, InitInp%chord, p, m, ErrStat2, ErrMsg2); if(Failed()) return

   ! Set parameters from InputFileData (need Misc allocated)
   CALL FVW_SetParametersFromInputFile(InputFileData, p, m, ErrStat2, ErrMsg2); if(Failed()) return
   CALL FVW_ToString(p, m) ! Print to screen

   ! Initialize States Vars
   CALL FVW_InitStates( x, p, m, ErrStat2, ErrMsg2 ); if(Failed()) return

   ! Initialize Constraints Vars
   CALL FVW_InitConstraint( z, p, m, ErrStat2, ErrMsg2 ); if(Failed()) return

   ! Panelling wings based on initial input mesh provided
   ! This mesh is now a cousin of the BladeMotion mesh from AD.
   CALL Wings_Panelling     (u%WingsMesh, p, m, ErrStat2, ErrMsg2); if(Failed()) return

   ! Returned guessed locations where wind will be required
   CALL SetRequestedWindPoints(y%r_wind, x, p, m, ErrStat2, ErrMsg2 ); if(Failed()) return
   ! Return anything in FVW_InitOutput that should be passed back to the calling code here

CONTAINS

   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FVW_Init') 
        Failed =  ErrStat >= AbortErrLev
        !if (Failed) call CleanUp()
   end function Failed

end subroutine FVW_Init

! ==============================================================================
subroutine FVW_InitMiscVars( p, m, ErrStat, ErrMsg )
   type(FVW_ParameterType),         intent(in   )  :: p              !< Parameters
   type(FVW_MiscVarType),           intent(inout)  :: m              !< Initial misc/optimization variables
   integer(IntKi),                  intent(  out)  :: ErrStat        !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None
   integer(IntKi)          :: ErrStat2       ! temporary error status of the operation
   character(ErrMsgLen)    :: ErrMsg2        ! temporary error message
   character(*), parameter :: RoutineName = 'FVW_InitMiscVars'
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

   m%FirstCall = .True.
   m%nNW       = 0      ! Number of active nearwake panels
   m%nFW       = 0      ! Number of active farwake  panels

   call AllocAry( m%LE      ,  3  ,  p%nSpan+1  , p%nWings, 'Leading Edge Points', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%LE = -999999_ReKi;
   call AllocAry( m%TE      ,  3  ,  p%nSpan+1  , p%nWings, 'TrailingEdge Points', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%TE = -999999_ReKi;
   call AllocAry( m%s_LL          ,  p%nSpan+1  , p%nWings, 'Spanwise coord LL  ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%s_LL= -999999_ReKi;
   call AllocAry( m%chord_LL      ,  p%nSpan+1  , p%nWings, 'Chord on LL        ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%chord_LL= -999999_ReKi;
   ! Variables at control points/elements
   call AllocAry( m%Gamma_LL,        p%nSpan    , p%nWings, 'Lifting line Circulation', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%Gamma_LL = -999999_ReKi;
   call AllocAry( m%chord_CP_LL   ,  p%nSpan    , p%nWings, 'Chord on CP LL     ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%chord_CP_LL= -999999_ReKi;
   call AllocAry( m%s_CP_LL       ,  p%nSpan    , p%nWings, 'Spanwise coord CPll', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%s_CP_LL= -999999_ReKi;
   call AllocAry( m%CP_LL   , 3   ,  p%nSpan    , p%nWings, 'Control points LL  ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%CP_LL= -999999_ReKi;
   call AllocAry( m%Tang    , 3   ,  p%nSpan    , p%nWings, 'Tangential vector  ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%Tang= -999999_ReKi;
   call AllocAry( m%Norm    , 3   ,  p%nSpan    , p%nWings, 'Normal     vector  ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%Norm= -999999_ReKi;
   call AllocAry( m%Orth    , 3   ,  p%nSpan    , p%nWings, 'Orthogonal vector  ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%Orth= -999999_ReKi;
   call AllocAry( m%Vind_LL , 3   ,  p%nSpan    , p%nWings, 'Vind on CP ll      ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%Vind_LL= -999999_ReKi;
   call AllocAry( m%Vtot_LL , 3   ,  p%nSpan    , p%nWings, 'Vtot on CP ll      ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%Vtot_LL= -999999_ReKi;
   call AllocAry( m%Vstr_LL , 3   ,  p%nSpan    , p%nWings, 'Vstr on CP ll      ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%Vstr_LL= -999999_ReKi;
   call AllocAry( m%Vwnd_LL , 3   ,  p%nSpan    , p%nWings, 'Wind on CP ll      ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%Vwnd_LL= -999999_ReKi;
   ! Variables at panels points
   call AllocAry( m%r_LL    , 3   ,  p%nSpan+1  , 2        ,  p%nWings, 'Lifting Line Panels', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%r_LL= -999999_ReKi;
   call AllocAry( m%Vwnd_NW , 3   ,  p%nSpan+1  ,p%nNWMax+1,  p%nWings, 'Wind on NW ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%Vwnd_NW= -999_ReKi;
   call AllocAry( m%Vwnd_FW , 3   ,  p%nSpan+1  ,p%nFWMax+1,  p%nWings, 'Wind on FW ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%Vwnd_FW= -999_ReKi;
   call AllocAry( m%Vind_NW , 3   ,  p%nSpan+1  ,p%nNWMax+1,  p%nWings, 'Vind on NW ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%Vind_NW= -999_ReKi;
   call AllocAry( m%Vind_FW , 3   ,  p%nSpan+1  ,p%nFWMax+1,  p%nWings, 'Vind on FW ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%Vind_FW= -999_ReKi;
   m%Vwnd_NW(1,:,:,:)= 8 
   m%Vwnd_NW(2,:,:,:)= 0
   m%Vwnd_NW(3,:,:,:)= 0
   m%Vwnd_FW(1,:,:,:)= 8 
   m%Vwnd_FW(2,:,:,:)= 0
   m%Vwnd_FW(3,:,:,:)= 0
end subroutine FVW_InitMiscVars
! ==============================================================================
subroutine FVW_InitStates( x, p, m, ErrStat, ErrMsg )
   type(FVW_ContinuousStateType),   intent(  out)  :: x              !< States
   type(FVW_ParameterType),         intent(in   )  :: p              !< Parameters
   type(FVW_MiscVarType),           intent(in   )  :: m              !< Initial misc/optimization variables
   integer(IntKi),                  intent(  out)  :: ErrStat        !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None
   integer(IntKi)          :: ErrStat2       ! temporary error status of the operation
   character(ErrMsgLen)    :: ErrMsg2        ! temporary error message
   character(*), parameter :: RoutineName = 'FVW_InitMiscVars'
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

   call AllocAry( x%Gamma_NW,    p%nSpan   , p%nNWMax  , p%nWings, 'NW Panels Circulation', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitStates' ); x%Gamma_NW = -999999_ReKi;
   call AllocAry( x%Gamma_FW,       1      , p%nFWMax  , p%nWings, 'FW Panels Circulation', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitStates' ); x%Gamma_FW = -999999_ReKi;
   call AllocAry( x%r_NW    , 3, p%nSpan+1 , p%nNWMax+1, p%nWings, 'NW Panels Points'     , ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitStates' ); x%r_NW     = -999999_ReKi;
   call AllocAry( x%r_FW    , 3,    2      , p%nFWMax+1, p%nWings, 'FW Panels Points'     , ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitStates' ); x%r_FW     = -999999_ReKi;


   if (ErrStat >= AbortErrLev) return
end subroutine FVW_InitStates
! ==============================================================================
subroutine FVW_InitConstraint( z, p, m, ErrStat, ErrMsg )
   type(FVW_ConstraintStateType),   intent(  out)  :: z              !< Constraints
   type(FVW_ParameterType),         intent(in   )  :: p              !< Parameters
   type(FVW_MiscVarType),           intent(in   )  :: m              !< Initial misc/optimization variables
   integer(IntKi),                  intent(  out)  :: ErrStat        !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None
   integer(IntKi)          :: ErrStat2       ! temporary error status of the operation
   character(ErrMsgLen)    :: ErrMsg2        ! temporary error message
   character(*), parameter :: RoutineName = 'FVW_InitMiscVars'
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""
   !
   call AllocAry( z%Gamma_LL,  p%nSpan, p%nWings, 'Lifting line Circulation', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitConstraint' ); z%Gamma_LL = -999999_ReKi;

   if (ErrStat >= AbortErrLev) return
end subroutine FVW_InitConstraint


! ==============================================================================
!> Setting parameters *and misc* from module inputs
SUBROUTINE FVW_SetParametersFromInputs( InitInp, p, m, ErrStat, ErrMsg )
   type(FVW_InitInputType),    intent(inout)  :: InitInp       !< Input data for initialization routine  (inout so we can use MOVE_ALLOC)
   type(FVW_ParameterType),    intent(inout) :: p             !< Parameters
   type(FVW_MiscVarType),      intent(inout) :: m             !< Misc
   integer(IntKi),             intent(  out) :: ErrStat       !< Error status of the operation
   character(*),               intent(  out) :: ErrMsg        !< Error message if ErrStat /= ErrID_None
   ! Local variables
   integer(IntKi)          :: iW             ! Index on wings
   integer(IntKi)          :: iSpan
   integer(IntKi)          :: ErrStat2
   character(ErrMsgLen)    :: ErrMsg2
   character(*), parameter :: RoutineName = 'FVW_SetParametersFromInputs'
   ErrStat = ErrID_None
   ErrMsg  = ""
   ! 
   p%nWings       =  InitInp%NumBlades
   
   ! NOTE: temporary limitation, all wings have the same nspan
   p%nSpan        =  InitInp%numBladeNodes-1

end subroutine FVW_SetParametersFromInputs
! ==============================================================================
!>
SUBROUTINE FVW_SetParametersFromInputFile( InputFileData, p, m, ErrStat, ErrMsg )
   type(FVW_InputFile),        intent(in   ) :: InputFileData !< Data stored in the module's input file
   type(FVW_ParameterType),    intent(inout) :: p             !< Parameters
   type(FVW_MiscVarType),      intent(inout) :: m             !< Misc
   integer(IntKi),             intent(  out) :: ErrStat       !< Error status of the operation
   character(*),               intent(  out) :: ErrMsg        !< Error message if ErrStat /= ErrID_None
   ! Local variables
   integer(IntKi)       :: ErrStat2
   character(ErrMsgLen) :: ErrMsg2
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Set parameters from input file
   p%IntMethod         = InputFileData%IntMethod
   p%CirculationMethod = InputFileData%CirculationMethod
   p%FreeWake          = InputFileData%FreeWake

   if (allocated(p%PrescribedCirculation)) deallocate(p%PrescribedCirculation)
   if (InputFileData%CirculationMethod==idCircPrescribed) then 
      call AllocAry( p%PrescribedCirculation,  p%nSpan, 'Prescribed Circulation', ErrStat2, ErrMsg2 );    call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_SetParameters' );    p%PrescribedCirculation = -999999_ReKi;
      if (.not. allocated(m%s_CP_LL)) then
         print*,'Spanwise coordinate not allocated, TODO'
         STOP
      endif
      call ReadAndInterpGamma(trim(InputFileData%CirculationFile), m%s_CP_LL(1:p%nSpan,1), m%s_LL(p%nSpan+1,1), p%PrescribedCirculation)
   endif

end subroutine FVW_SetParametersFromInputFile

subroutine FVW_ToString(p,m)
   type(FVW_ParameterType), intent(in)       :: p !< Parameters
   type(FVW_MiscVarType),      intent(inout) :: m !< Misc
   print*,'------------------------------------'
end subroutine FVW_ToString

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
subroutine FVW_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )

   type(FVW_InputType),             intent(inout)  :: u(2)        !< System inputs
   type(FVW_ParameterType),         intent(inout)  :: p           !< Parameters
   type(FVW_ContinuousStateType),   intent(inout)  :: x           !< Continuous states
   type(FVW_DiscreteStateType),     intent(inout)  :: xd          !< Discrete states
   type(FVW_ConstraintStateType),   intent(inout)  :: z           !< Constraint states
   type(FVW_OtherStateType),        intent(inout)  :: OtherState  !< Other states
   type(FVW_OutputType),            intent(inout)  :: y           !< System outputs
   type(FVW_MiscVarType),           intent(inout)  :: m           !< Misc/optimization variables
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""
   ! Place any last minute operations or calculations here:
   ! Close files here:
   ! Destroy the input data:
   call FVW_DestroyInput( u(1), ErrStat, ErrMsg )
   call FVW_DestroyInput( u(2), ErrStat, ErrMsg )

   ! Destroy the parameter data:
   call FVW_DestroyParam( p, ErrStat, ErrMsg )

   ! Destroy the state data:
   call FVW_DestroyContState(   x,           ErrStat, ErrMsg )
   call FVW_DestroyDiscState(   xd,          ErrStat, ErrMsg )
   call FVW_DestroyConstrState( z,           ErrStat, ErrMsg )
   call FVW_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )
   call FVW_DestroyMisc(        m,           ErrStat, ErrMsg )

   ! Destroy the output data:
   call FVW_DestroyOutput( y, ErrStat, ErrMsg )

end subroutine FVW_End



!----------------------------------------------------------------------------------------------------------------------------------
!> Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete and other states.
!! Continuous, constraint, discrete, and other states are updated for t + Interval
subroutine FVW_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, m, errStat, errMsg )
!..................................................................................................................................
   real(DbKi),                      intent(in   )  :: t           !< Current simulation time in seconds
   integer(IntKi),                  intent(in   )  :: n           !< Current simulation time step n = 0,1,...
   type(FVW_InputType),             intent(inout)  :: u(:)        !< Inputs at utimes (out only for mesh record-keeping in ExtrapInterp routine)
   real(DbKi),                      intent(in   )  :: utimes(:)   !< Times associated with u(:), in seconds
   type(FVW_ParameterType),         intent(in   )  :: p           !< Parameters
   type(FVW_ContinuousStateType),   intent(inout)  :: x           !< Input: Continuous states at t; Output: at t+dt
   type(FVW_DiscreteStateType),     intent(inout)  :: xd          !< Input: Discrete states at t;   Output: at t+dt
   type(FVW_ConstraintStateType),   intent(inout)  :: z           !< Input: Constraint states at t; Output: at t+dt
   type(FVW_OtherStateType),        intent(inout)  :: OtherState  !< Input: Other states at t;      Output: at t+dt
   type(FVW_MiscVarType),           intent(inout)  :: m           !< Misc/optimization variables
   integer(IntKi),                  intent(  out)  :: errStat     !< Error status of the operation
   character(*),                    intent(  out)  :: errMsg      !< Error message if ErrStat /= ErrID_None
   ! local variables
   type(FVW_InputType)           :: uInterp                                                            ! Interpolated/Extrapolated input
   integer(IntKi)                :: ErrStat2                                                           ! temporary Error status
   character(ErrMsgLen)          :: ErrMsg2                                                            ! temporary Error message
   type(FVW_ConstraintStateType) :: z_guess                                                                              ! <
   integer(IntKi) :: iW, iSpan, iAge

   ErrStat = ErrID_None
   ErrMsg  = ""

   print'(A,F10.3,A,F10.3,A,F10.3,A,I0,A,I0)','     Update states, t:',t,'  t_u:', utimes(1),' dt: ',utimes(1)-t,'   ',n,' nNW:',m%nNW

   ! Panelling wings based on input mesh provided
   CALL Wings_Panelling(u(1)%WingsMesh, p, m, ErrStat2, ErrMsg2); if(Failed()) return
   
   ! Distribute the Wind we requested to Inflow wind to storage Misc arrays
   ! TODO ANDY
   !CALL DistributeRequestedWind(u(1)%V_wind, x, p, m, ErrStat2, ErrMsg2);  if(Failed()) return

   ! Solve for circulation at t
   CALL FVW_CalcConstrStateResidual(t, u(1), p, x, xd, z_guess, OtherState, m, z, ErrStat2, ErrMsg2); if(Failed()) return

   ! Map circulation and positions between LL and NW 
   ! Changes: x only
   CALL Wings_Map_LL_NW(p, m, z, x, ErrStat, ErrMsg )

   if (p%IntMethod .eq. idEuler1) then 
     call FVW_Euler1( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat2, ErrMsg2); if(Failed()) return

   !elseif (p%IntMethod .eq. idRK4) then 
   !   call FVW_RK4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
   !elseif (p%IntMethod .eq. idAB4) then
   !   call FVW_AB4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
   !elseif (p%IntMethod .eq. idABM4) then
   !   call FVW_ABM4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
   else  
      call SetErrStat(ErrID_Fatal,'Invalid time integration method:'//Num2LStr(p%IntMethod),ErrStat,ErrMsg,'FVW_UpdateState') 
   end IF

   ! -- Propagate wake m%nNW+1+1
   do iW=1,p%nWings
      do iAge=p%nNWMax+1,2,-1
         do iSpan=1,p%nSpan+1
            x%r_NW(1:3,iSpan,iAge,iW) = x%r_NW(1:3,iSpan,iAge-1,iW)
         enddo
      enddo
      x%r_NW(1:3,:,1,iW) = -999.0_ReKi
   enddo
   do iW=1,p%nWings
      do iAge=p%nNWMax,2,-1
         do iSpan=1,p%nSpan
            x%Gamma_NW(iSpan,iAge,iW) = x%Gamma_NW(iSpan,iAge-1,iW)
         enddo
      enddo
      x%Gamma_NW(:,1,iW) = -999.0_ReKi
   enddo
!    do iAge=1,m%nNW+1+1
!       print*,'iAge',iAge
!       print*,'Prop', x%r_NW(1, 1:p%nSpan+1, iAge,1)
!       print*,'Prop', x%r_NW(2, 1:p%nSpan+1, iAge,1)
!       print*,'Prop', x%r_NW(3, 1:p%nSpan+1, iAge,1)
!    enddo
!    do iW=1,p%nWings
!       print*,'Wing',iW
!       do iAge=1,m%nNW+1+1
!          print*,'Prop',x%Gamma_NW(:,iAge,iW)
!       enddo
!    enddo

   ! Solve for circulation at t+dt
   CALL FVW_CalcConstrStateResidual(t, u(1), p, x, xd, z_guess, OtherState, m, z, ErrStat2, ErrMsg2); if(Failed()) return
   
   ! Map circulation and positions between LL and NW 
   ! Changes: x only
   CALL Wings_Map_LL_NW(p, m, z, x, ErrStat, ErrMsg )

!    do iW=1,p%nWings
!       print*,'Wing',iW
!       do iAge=1,m%nNW+1+1
!          print*,'Map',x%Gamma_NW(:,iAge,iW)
!       enddo
!    enddo
!    do iAge=1,m%nNW+1+1
!       print*,'iAge',iAge
!       print*,'Map', x%r_NW(1, 1:p%nSpan+1, iAge,1)
!       print*,'Map', x%r_NW(2, 1:p%nSpan+1, iAge,1)
!       print*,'Map', x%r_NW(3, 1:p%nSpan+1, iAge,1)
!    enddo

!    if (m%nNW>3) STOP

   if (m%FirstCall) then
      m%FirstCall=.False.
   endif


   call FVW_DestroyConstrState(z_guess, ErrStat2, ErrMsg2);

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FVW_UpdateStates') 
      Failed =  ErrStat >= AbortErrLev
      !if (Failed) call CleanUp()
   end function Failed

end subroutine FVW_UpdateStates


!----------------------------------------------------------------------------------------------------------------------------------
!> This is a tight coupling routine for computing derivatives of continuous states.
subroutine FVW_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )
   real(DbKi),                    intent(in   ) :: t          !< Current simulation time in seconds
   type(FVW_InputType),           intent(in   ) :: u          !< Inputs at t
   type(FVW_ParameterType),       intent(in   ) :: p          !< Parameters
   type(FVW_ContinuousStateType), intent(in   ) :: x          !< Continuous states at t
   type(FVW_DiscreteStateType),   intent(in   ) :: xd         !< Discrete states at t
   type(FVW_ConstraintStateType), intent(in   ) :: z          !< Constraint states at t
   type(FVW_OtherStateType),      intent(in   ) :: OtherState !< Other states at t
   type(FVW_MiscVarType),         intent(inout) :: m          !< Misc variables for optimization (not copied in glue code)
   type(FVW_ContinuousStateType), intent(  out) :: dxdt       !< Continuous state derivatives at t
   integer(IntKi),                intent(  out) :: ErrStat    !< Error status of the operation
   character(*),                  intent(  out) :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   ! Local variables
   integer(IntKi)           :: iSpan,iAge, iW, nSeg, nSegP
   real(ReKi), dimension(3) :: U_mean
   integer(IntKi),dimension(:,:), allocatable :: SegConnct !< Segment connectivity
   real(ReKi),    dimension(:,:), allocatable :: SegPoints !< Segment Points
   real(ReKi),    dimension(:)  , allocatable :: SegGamma  !< Segment Circulation
   real(ReKi),    dimension(:)  , allocatable :: SegSmooth  !< Segment smooth parameter
   real(ReKi),    dimension(:,:),   allocatable :: CPs   !< ControlPoints
   real(ReKi),    dimension(:,:),   allocatable :: Uind  !< Induced velocity
   integer(IntKi)           :: SmoothModel

   call AllocAry( dxdt%r_NW , 3   ,  p%nSpan+1  ,p%nNWMax+1,  p%nWings, 'Wind on NW ', ErrStat, ErrMsg ); dxdt%r_NW= -999999_ReKi;
   call AllocAry( dxdt%r_FW , 3   ,      2      ,p%nFWMax+1,  p%nWings, 'Wind on FW ', ErrStat, ErrMsg ); dxdt%r_FW= -999999_ReKi;

   if (p%FreeWake) then
      print*,'TODO free wake convection'
      STOP
   else

      ! --- Packing all vortex elements into a list of segments
      call PackAllPanelsToSegments(p, m, x, z, SegConnct, SegPoints, SegGamma, nSeg, nSegP)
      print*,'Number of segments',nSeg, 'Number of points',nSegP

      ! --- Computing induced velocity
      ! TODO for now point by point..
      allocate(SegSmooth(1:nSeg));
      SegSmooth=10
      allocate(Uind(1:3,1))
      allocate(CPs (1:3,1))
      CPs=0
      Uind=0
      SmoothModel=idSegSmoothLambOseen
      ! --- On NW
      do iW=1,p%nWings
         do iAge=1,m%nNW+1
            do iSpan=1,p%nSpan+1
               Uind(1:3,1)=0.0_ReKi
               CPs(1:3,1) = x%r_NW(1:3,iSpan,iAge,iW)

               CALL ui_seg(CPs, 1, 1, 1, &
                     SegPoints, SegConnct, SegGamma, 1, nSeg, nSeg, nSegP,   &
                     SmoothModel, SegSmooth, Uind)
!                print*,'Uind',Uind
               m%Vind_NW(1:3,iSpan,iAge,iW) = Uind(1:3,1)
            enddo
         enddo
      enddo
      deallocate(SegConnct)
      deallocate(SegGamma)
      deallocate(SegPoints)
      deallocate(SegSmooth)
      deallocate(Uind)
      deallocate(CPs)

      U_mean(1:3)=0
      do iW=1,p%nWings; do iAge=1,m%nNW+1; do iSpan=1,p%nSpan+1;  
         U_mean(1:3)= U_mean(1:3)+ m%Vind_NW(1:3, iSpan, iAge, iW)
      enddo; enddo; enddo
      U_mean(1:3) = U_mean(1:3)/ ((m%nNW+1)*(p%nSpan+1)*p%nWings)
      print*,'Mean convection velocity NW: ',U_mean(1:3)


      U_mean(1:3)=0
      do iW=1,p%nWings; do iAge=1,m%nNW+1; do iSpan=1,p%nSpan+1;  
         U_mean(1:3)= U_mean(1:3)+ m%Vwnd_NW(1:3, iSpan, iAge, iW)
      enddo; enddo; enddo
      U_mean(1:3) = U_mean(1:3)/ ((m%nNW+1)*(p%nSpan+1)*p%nWings)
      print*,'Mean convection velocity NW: ',U_mean(1:3)
      U_mean(1:3)=0
      do iW=1,p%nWings; do iAge=1,m%nFW+1; do iSpan=1,2;
         U_mean(1:3)= U_mean(1:3)+ m%Vwnd_FW(1:3, iSpan, iAge, iW)
      enddo; enddo; enddo
      U_mean(1:3) = U_mean(1:3)/ ((m%nFW+1)*(2)*p%nWings)
      print*,'Mean convection velocity FW: ',U_mean(1:3)


      ! --- Vortex points are convected with the free stream
      dxdt%r_NW(1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) = m%Vwnd_NW(1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) 
      dxdt%r_FW(1:3, 1:2        , 1:m%nFW+1, 1:p%nWings) = m%Vwnd_FW(1:3, 1:2        , 1:m%nFW+1, 1:p%nWings)

      dxdt%r_NW(1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) = m%Vwnd_NW(1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) +  m%Vind_NW(1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings)
   endif
end subroutine FVW_CalcContStateDeriv

!----------------------------------------------------------------------------------------------------------------------------------
subroutine FVW_Euler1( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
   real(DbKi),                    intent(in   ) :: t          !< Current simulation time in seconds
   integer(IntKi),                intent(in   ) :: n          !< time step number
   type(FVW_InputType),           intent(inout) :: u(:)       !< Inputs at t
   real(DbKi),                    intent(in   ) :: utimes(:)  !< times of input
   type(FVW_ParameterType),       intent(in   ) :: p          !< Parameters
   type(FVW_ContinuousStateType), intent(inout) :: x          !< Continuous states at t on input at t + dt on output
   type(FVW_DiscreteStateType),   intent(in   ) :: xd         !< Discrete states at t
   type(FVW_ConstraintStateType), intent(in   ) :: z          !< Constraint states at t (possibly a guess)
   type(FVW_OtherStateType),      intent(inout) :: OtherState !< Other states at t on input at t + dt on output
   type(FVW_MiscVarType),         intent(inout) :: m          !< Misc/optimization variables
   integer(IntKi),                intent(  out) :: ErrStat    !< Error status of the operation
   character(*),                  intent(  out) :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   ! local variables
   type(FVW_ContinuousStateType) :: dxdt        ! time derivatives of continuous states      
   real(ReKi)                    :: dt          ! Time step
   integer(IntKi) :: iAge
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = "" 

   ! Compute "right hand side"
   CALL FVW_CalcContStateDeriv( t, u(1), p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )

   dt = utimes(1)-t ! TODO TODO is this correct 

!    do iAge=1,m%nNW+1
!       print*,'iAge',iAge
!       print*,'Before', x%r_NW(1, 1:p%nSpan+1, iAge,1)
!       print*,'Before', x%r_NW(2, 1:p%nSpan+1, iAge,1)
!       print*,'Before', x%r_NW(3, 1:p%nSpan+1, iAge,1)
!    enddo

   ! Update of positions
   x%r_NW(1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) = x%r_NW(1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) +  dt * dxdt%r_NW(1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings)
   x%r_FW(1:3,    1:2     , 1:m%nFW+1, 1:p%nWings) = x%r_FW(1:3, 1:2        , 1:m%nFW+1, 1:p%nWings) +  dt * dxdt%r_FW(1:3, 1:2        , 1:m%nFW+1, 1:p%nWings)

!    do iAge=1,m%nNW+1
!       print*,'iAge',iAge
!       print*,'After', x%r_NW(1, 1:p%nSpan+1, iAge,1)
!       print*,'After', x%r_NW(2, 1:p%nSpan+1, iAge,1)
!       print*,'After', x%r_NW(3, 1:p%nSpan+1, iAge,1)
!    enddo
   ! Update of Gamma
   ! TODO, viscous diffusion, stretching

   call FVW_DestroyContState(dxdt, ErrStat, ErrMsg)

end subroutine FVW_Euler1
!----------------------------------------------------------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------------------------------------------------
!> This is a tight coupling routine for solving for the residual of the constraint state functions.
subroutine FVW_CalcConstrStateResidual( t, u, p, x, xd, z, OtherState, m, z_out, ErrStat, ErrMsg )
   real(DbKi),                    intent(in   )  :: t           !< Current simulation time in seconds
   type(FVW_InputType),           intent(in   )  :: u           !< Inputs at t
   type(FVW_ParameterType),       intent(in   )  :: p           !< Parameters
   type(FVW_ContinuousStateType), intent(in   )  :: x           !< Continuous states at t
   type(FVW_DiscreteStateType),   intent(in   )  :: xd          !< Discrete states at t
   type(FVW_ConstraintStateType), intent(in   )  :: z           !< Constraint states at t (possibly a guess)
   type(FVW_OtherStateType),      intent(in   )  :: OtherState  !< Other states at t
   type(FVW_MiscVarType),         intent(inout)  :: m           !< Misc variables for optimization (not copied in glue code)
   type(FVW_ConstraintStateType), intent(  out)  :: z_out            !< Residual of the constraint state functions using
                                                                    !!     the input values described above
   integer(IntKi),                    intent(  OUT)  :: ErrStat     !< Error status of the operation
   character(*),                      intent(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""
   ! Solve for the residual of the constraint state functions here:
   !z%residual = 0.0_ReKi
   !z%Gamma_LL = 0.0_ReKi
   call AllocAry( z_out%Gamma_LL,  p%nSpan, p%nWings, 'Lifting line Circulation', ErrStat, ErrMsg );
   z_out%Gamma_LL = -999999_ReKi;

   CALL Wings_ComputeCirculation(t, z_out%Gamma_LL, z%Gamma_LL, u, p, x, m, ErrStat, ErrMsg)

end subroutine FVW_CalcConstrStateResidual

!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
!! This subroutine is used to compute the output channels (motions and loads) and place them in the WriteOutput() array.
!! The descriptions of the output channels are not given here. Please see the included OutListParameters.xlsx sheet for
!! for a complete description of each output parameter.
subroutine FVW_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
! NOTE: no matter how many channels are selected for output, all of the outputs are calculated
! All of the calculated output channels are placed into the m%AllOuts(:), while the channels selected for outputs are
! placed in the y%WriteOutput(:) array.
!..................................................................................................................................
   real(DbKi),                      intent(in   )  :: t           !< Current simulation time in seconds
   type(FVW_InputType),             intent(in   )  :: u           !< Inputs at Time t
   type(FVW_ParameterType),         intent(in   )  :: p           !< Parameters
   type(FVW_ContinuousStateType),   intent(in   )  :: x           !< Continuous states at t
   type(FVW_DiscreteStateType),     intent(in   )  :: xd          !< Discrete states at t
!FIXME:TODO: AD15_CalcOutput has constraint states as intent(in) only. This is forcing me to store z in the AD15 miscvars for now.
   type(FVW_ConstraintStateType),   intent(inout)  :: z           !< Constraint states at t
   type(FVW_OtherStateType),        intent(in   )  :: OtherState  !< Other states at t
   type(FVW_OutputType),            intent(inout)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                                  !!   nectivity information does not have to be recalculated)
   type(FVW_MiscVarType),           intent(inout)  :: m           !< Misc/optimization variables
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   ! Local variables
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   character(*), parameter       :: RoutineName = 'FVW_CalcOutput'

   ErrStat = ErrID_None
   ErrMsg  = ""

   print'(A,F10.3,A,L1)','     CalcOutput     t:',t,'   ',m%FirstCall

   if (m%FirstCall) then
      print*,'>>> First Call of CalcOuput, calling panelling and constrstate'
      CALL Wings_Panelling(u%WingsMesh, p, m, ErrStat2, ErrMsg2); if(Failed()) return
      CALL Wings_ComputeCirculation(t, m%Gamma_LL, z%Gamma_LL, u, p, x, m, ErrStat, ErrMsg) ! For plotting only
   else
      m%Gamma_LL = z%Gamma_LL ! For plotting only
   endif

   if (.not. allocated(y%Vind)) then
      !call AllocAry( y%Vind ,  3, p%nSpan, p%nWings, 'Induced velocity vector',  ErrStat2, ErrMsg2 );
      call AllocAry( y%Vind ,  3, p%nSpan, 3, 'Induced velocity vector',  ErrStat2, ErrMsg2 ); ! TODO TODO TODO Hack nWings=3 for output
      if(Failed()) return
   endif
   ! Returned guessed locations where wind will be required
   CALL SetRequestedWindPoints(y%r_wind, x, p, m, ErrStat2, ErrMsg2 ); if(Failed()) return

   !if (m%FirstCall) then
   !endif

   ! For now returning 0
   y%Vind(1,:,:) = 0.0_ReKi
   y%Vind(2,:,:) = 0.0_ReKi
   y%Vind(3,:,:) = 0.0_ReKi

   ! We don't propagate the "Old"-> "New" if update states was not called once
   ! This is introduced since at init, CalcOutput is called before UpdateState
   if (.not. m%FirstCall) then 
       call PrepareNextTimeStep()
   endif

contains
   !==========================================================================
   !> Computes induction on the lifting line (3/4 chord point)
   ! Interpolate the values at the radial station of AeroDyn
   subroutine CalcInduction_LL()

   end subroutine CalcInduction_LL


   subroutine PrepareNextTimeStep()
      ! --- Propagate wake
      ! TODO
      
      ! 
      m%nNW=m%nNW+1
      if (m%nNW>p%nNWMax) m%nNW = p%nNWMax

   end subroutine PrepareNextTimeStep

   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FVW_CalcOutput') 
      Failed =  ErrStat >= AbortErrLev
      !if (Failed) call CleanUp()
   end function Failed
   !==========================================================================
end subroutine FVW_CalcOutput


END MODULE FVW
