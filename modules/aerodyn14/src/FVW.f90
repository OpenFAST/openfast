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

   ! Trigger required before allocations
   p%nNWMax  = max(InputFileData%nNWPanels,0)+1          ! +1 since LL panel included in NW
   p%nFWMax  = max(InputFileData%nFWPanels,0)
   p%nFWFree = max(InputFileData%nFWPanelsFree,0)

   if (InputFileData%HACK==1) then
      p%nWings=1 ! Elliptical wing temporary hack
   endif

   ! Initialize Misc Vars (may depend on input file)
   CALL FVW_InitMiscVars( p, m, ErrStat2, ErrMsg2 ); if(Failed()) return

   ! Wind Speed hack, TODO temporary
   m%Vwnd_LL(:,:,:)   = 0
   m%Vwnd_NW(:,:,:,:) = 0
   m%Vwnd_FW(:,:,:,:) = 0
   m%Vwnd_LL(1,:,:)   = InputFileData%Uinf
   m%Vwnd_NW(1,:,:,:) = InputFileData%Uinf
   m%Vwnd_FW(1,:,:,:) = InputFileData%Uinf
   if (InputFileData%HACK==1) then
      m%Vwnd_LL(3,:,:)   =0.1
      m%Vwnd_NW(3,:,:,:) =0.1
      m%Vwnd_FW(3,:,:,:) =0.1
   endif


   ! Preliminary meshing of the wings (may depend on input file)
   ! NOTE: the mesh is not located at the right position yet, the first call to calcoutput will redo some meshing
   CALL Wings_Panelling_Init(InitInp%WingsMesh, InitInp%RElm, InitInp%chord, p, m, ErrStat2, ErrMsg2); if(Failed()) return

   ! Set parameters from InputFileData (need Misc allocated)
   CALL FVW_SetParametersFromInputFile(InputFileData, p, m, ErrStat2, ErrMsg2); if(Failed()) return
   CALL FVW_ToString(p, m) ! Print to screen

   ! Initialize States Vars
   CALL FVW_InitStates( x, p, m, ErrStat2, ErrMsg2 ); if(Failed()) return

   ! Initialize Constraints Vars
   CALL FVW_InitConstraint( z, p, m, ErrStat2, ErrMsg2 ); if(Failed()) return

   ! Panelling wings based on initial input mesh provided
   ! NOTE: the mesh is not located at the right position yet, the first call to calcoutput will redo some meshing
   CALL Wings_Panelling     (InitInp%WingsMesh, p, m, ErrStat2, ErrMsg2); if(Failed()) return

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
   m%nNW       = iNWStart-1    ! Number of active nearwake panels
   m%nFW       = 0             ! Number of active farwake  panels
   m%iStep     = 0             ! Current step number

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
   call AllocAry( m%dl      , 3   ,  p%nSpan    , p%nWings, 'Orthogonal vector  ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%dl= -999999_ReKi;
   call AllocAry( m%Area    ,        p%nSpan    , p%nWings, 'LL Panel area      ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%Area = -999999_ReKi;
   call AllocAry( m%Vind_LL , 3   ,  p%nSpan    , p%nWings, 'Vind on CP ll      ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%Vind_LL= -999999_ReKi;
   call AllocAry( m%Vtot_LL , 3   ,  p%nSpan    , p%nWings, 'Vtot on CP ll      ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%Vtot_LL= -999999_ReKi;
   call AllocAry( m%Vstr_LL , 3   ,  p%nSpan    , p%nWings, 'Vstr on CP ll      ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%Vstr_LL= -999999_ReKi;
   call AllocAry( m%Vwnd_LL , 3   ,  p%nSpan    , p%nWings, 'Wind on CP ll      ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%Vwnd_LL= -999999_ReKi;
   ! Variables at panels points
   call AllocAry( m%r_LL    , 3   ,  p%nSpan+1  , 2        ,  p%nWings, 'Lifting Line Panels', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%r_LL= -999999_ReKi;
   call AllocAry( m%Vwnd_NW , 3   ,  p%nSpan+1  ,p%nNWMax+1,  p%nWings, 'Wind on NW ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%Vwnd_NW= -999_ReKi;
   call AllocAry( m%Vwnd_FW , 3   ,  p%nSpan+1  ,p%nFWMax+1,  p%nWings, 'Wind on FW ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%Vwnd_FW= -999_ReKi;
   call AllocAry( m%Vind_NW , 3   ,  p%nSpan+1  ,p%nNWMax+1,  p%nWings, 'Vind on NW ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%Vind_NW= -999_ReKi;
   call AllocAry( m%Vind_FW , 3   ,  FWnSpan+1  ,p%nFWMax+1,  p%nWings, 'Vind on FW ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitMisc' ); m%Vind_FW= -999_ReKi;
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
   call AllocAry( x%Gamma_FW,    FWnSpan   , p%nFWMax  , p%nWings, 'FW Panels Circulation', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitStates' ); x%Gamma_FW = -999999_ReKi;
   call AllocAry( x%r_NW    , 3, p%nSpan+1 , p%nNWMax+1, p%nWings, 'NW Panels Points'     , ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitStates' ); x%r_NW     = -99_ReKi;
   call AllocAry( x%r_FW    , 3, FWnSpan+1 , p%nFWMax+1, p%nWings, 'FW Panels Points'     , ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitStates' ); x%r_FW     = -99_ReKi;


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
   p%nWings       =  InitInp%NumBl
   
   ! NOTE: temporary limitation, all wings have the same nspan
   p%nSpan        =  size(InitInp%RElm)-1

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
   p%IntMethod            = InputFileData%IntMethod
   p%CirculationMethod    = InputFileData%CirculationMethod
   p%CircSolvConvCrit     = InputFileData%CircSolvConvCrit
   p%CircSolvRelaxation   = InputFileData%CircSolvRelaxation
   p%CircSolvMaxIter      = InputFileData%CircSolvMaxIter
   p%FreeWakeStart        = InputFileData%FreeWakeStart
   p%CircSolvPolar        = InputFileData%CircSolvPolar
   p%FullCirculationStart = InputFileData%FullCirculationStart
   p%RegFunction          = InputFileData%RegFunction
   p%WakeRegMethod        = InputFileData%WakeRegMethod
   p%WakeRegFactor        = InputFileData%WakeRegFactor
   p%WingRegFactor        = InputFileData%WingRegFactor
   p%WrVTK                = InputFileData%WrVTK
   p%VTKBlades            = min(max(InputFileData%VTKBlades,0),p%nWings)
   p%HACK                 = InputFileData%HACK

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

   type(FVW_InputType),             intent(inout)  :: u           !< System inputs
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
   call FVW_DestroyInput( u, ErrStat, ErrMsg )

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
   real(ReKi) :: dt

   ErrStat = ErrID_None
   ErrMsg  = ""

   dt=utimes(1)-t ! TODO TODO TODO
   m%iStep=n

   print'(A,F10.3,A,F10.3,A,F10.3,A,I0,A,I0,A,I0)','Update states, t:',t,'  t_u:', utimes(1),' dt: ',dt,'   ',n,' nNW:',m%nNW,' nFW:',m%nFW

   ! --- Evaluation at t
   ! Inputs at t
   ! TODO TODO TODO AD14 HACK!  inputs at other times are wrong with AD14
   !   call FVW_CopyInput( u(1), uInterp, MESH_NEWCOPY, ErrStat2, ErrMsg2); if(Failed()) return
   !    call FVW_Input_ExtrapInterp(u,utimes,uInterp,t, ErrStat2, ErrMsg2); if(Failed()) return
   if (m%FirstCall) then
      call FVW_CopyInput( u(2), uInterp, MESH_NEWCOPY, ErrStat2, ErrMsg2); if(Failed()) return
      ! NOTE AD14 HACK should be put outside!
      ! Panelling wings based on input mesh provided
      call Wings_Panelling(uInterp%WingsMesh, p, m, ErrStat2, ErrMsg2); if(Failed()) return
   else
      call FVW_CopyInput( u(1), uInterp, MESH_NEWCOPY, ErrStat2, ErrMsg2); if(Failed()) return
   endif
   
   ! Distribute the Wind we requested to Inflow wind to storage Misc arrays
   ! TODO ANDY
   !CALL DistributeRequestedWind(u(1)%V_wind, x, p, m, ErrStat2, ErrMsg2);  if(Failed()) return

   ! --- Solve for circulation at t
   ! Returns: z%Gamma_LL (at t)
   call AllocAry( z_guess%Gamma_LL,  p%nSpan, p%nWings, 'Lifting line Circulation', ErrStat, ErrMsg );
   z_guess%Gamma_LL = m%Gamma_LL
   call FVW_CalcConstrStateResidual(t, uInterp, p, x, xd, z_guess, OtherState, m, z, ErrStat2, ErrMsg2, 1); if(Failed()) return

   ! Map circulation and positions between LL and NW  and then NW and FW
   ! Changes: x only
   call Map_LL_NW(p, m, z, x, ErrStat2, ErrMsg2)
   call Map_NW_FW(p, m, z, x, ErrStat2, ErrMsg2)
   !call print_x_NW_FW(p, m, z, x,'Map_')

   ! --- Integration between t and t+dt
   if (p%IntMethod .eq. idEuler1) then 
     call FVW_Euler1( t, dt, uInterp, p, x, xd, z, OtherState, m, ErrStat2, ErrMsg2); if(Failed()) return
   !elseif (p%IntMethod .eq. idRK4) then 
   !   call FVW_RK4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
   !elseif (p%IntMethod .eq. idAB4) then
   !   call FVW_AB4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
   !elseif (p%IntMethod .eq. idABM4) then
   !   call FVW_ABM4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
   else  
      call SetErrStat(ErrID_Fatal,'Invalid time integration method:'//Num2LStr(p%IntMethod),ErrStat,ErrMsg,'FVW_UpdateState') 
   end IF
   !call print_x_NW_FW(p, m, z, x,'Conv')


   ! --- t+dt

   ! Propagation/creation of new layer of panels
   call PropagateWake(p, m, z, x, ErrStat2, ErrMsg2)
   !call print_x_NW_FW(p, m, z, x,'Prop_')

   ! Inputs at t+dt
   ! TODO TODO TODO  inputs at other times are wrong with AD14
   !call FVW_Input_ExtrapInterp(u,utimes,uInterp,t+dt, ErrStat2, ErrMsg2); if(Failed()) return
   call FVW_CopyInput( u(1), uInterp, MESH_NEWCOPY, ErrStat2, ErrMsg2); if(Failed()) return

   ! Panelling wings based on input mesh at t+dt
   call Wings_Panelling(uInterp%WingsMesh, p, m, ErrStat2, ErrMsg2); if(Failed()) return

   ! Updating positions of first NW and FW panels (Circulation also updated but irrelevant)
   ! Changes: x only
   call Map_LL_NW(p, m, z, x, ErrStat2, ErrMsg2)
   call Map_NW_FW(p, m, z, x, ErrStat2, ErrMsg2)

   ! --- Solve for circulation at t+dt
   ! Returns: z%Gamma_LL (at t+dt)
   z_guess%Gamma_LL = z%Gamma_LL ! We use as guess the circulation from the previous time step (see above)
   call FVW_CalcConstrStateResidual(t+dt, uInterp, p, x, xd, z_guess, OtherState, m, z, ErrStat2, ErrMsg2, 2); if(Failed()) return

   ! Updating circulation of near wake panel (and position but irrelevant)
   ! Changes: x only
   call Map_LL_NW(p, m, z, x, ErrStat2, ErrMsg2)
   call Map_NW_FW(p, m, z, x, ErrStat2, ErrMsg2)
   !call print_x_NW_FW(p, m, z, x,'Map_')

   !if (m%nFW>4) STOP
   !if (t>0.5) STOP

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
   integer(IntKi)          :: ErrStat2       ! temporary error status of the operation
   character(ErrMsgLen)    :: ErrMsg2        ! temporary error message
   integer(IntKi) :: nFWEff ! Number of farwake panels that are free at current tmie step

   ErrStat = ErrID_None
   ErrMsg  = ""

   call AllocAry( dxdt%r_NW , 3   ,  p%nSpan+1  ,p%nNWMax+1,  p%nWings, 'Wind on NW ', ErrStat, ErrMsg ); dxdt%r_NW= -999999_ReKi;
   call AllocAry( dxdt%r_FW , 3   ,  FWnSpan+1  ,p%nFWMax+1,  p%nWings, 'Wind on FW ', ErrStat, ErrMsg ); dxdt%r_FW= -999999_ReKi;

   if (t> p%FreeWakeStart) then
      nFWEff = min(m%nFW, p%nFWFree)
      
      ! --- Compute Induced velocities on the Near wake and far wake based on the marker postions:
      ! (expensive N^2 call)
      ! In  : x%r_NW,    r%r_FW 
      ! Out:  m%Vind_NW, m%Vind_FW
      call WakeInducedVelocities(p, x, m, ErrStat2, ErrMsg2)
      m%Vind_FW(1:3, 1:FWnSpan+1, p%nFWFree+1:p%nFWMax, 1:p%nWings) =0 ! Ensuring no convection velocity for panels that user doesnt want free

      call print_mean_4d( m%Vind_NW(:,:, 1:m%nNW+1,:), 'Mean induced vel. NW')
      if (nFWEff>0) then
         call print_mean_4d( m%Vind_FW(:,:, 1:nFWEff ,:), 'Mean induced vel. FW')
      endif
      call print_mean_4d( m%Vwnd_NW(:,:, 1:m%nNW+1,:), 'Mean wind vel.    NW')
      call print_mean_4d( m%Vwnd_FW(:,:, 1:m%nFW+1,:), 'Mean wind vel.    FW')

      ! --- Vortex points are convected with the free stream and induced velocity
      dxdt%r_NW(1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) = m%Vwnd_NW(1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) +  m%Vind_NW(1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings)
      dxdt%r_FW(1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings) = m%Vwnd_FW(1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings) +  m%Vind_FW(1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings)
   else

      call print_mean_4d( m%Vwnd_NW(:,:,1:m%nNW+1,:), 'Mean wind vel.    NW')
      !call print_mean_4d( m%Vwnd_FW(:,:,1:m%nFW+1,:), 'Mean wind vel.    FW')

      ! --- Vortex points are convected with the free stream
      dxdt%r_NW(1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) = m%Vwnd_NW(1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) 
      dxdt%r_FW(1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings) = m%Vwnd_FW(1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings)
   endif
   ! First NW point does not convect (bound to LL)
   dxdt%r_NW(1:3, :, 1:iNWStart-1, :)=0
   ! First FW point does not convect (bound to NW)
   !dxdt%r_FW(1:3, :, 1, :)=0
end subroutine FVW_CalcContStateDeriv

!----------------------------------------------------------------------------------------------------------------------------------
subroutine FVW_Euler1( t, dt, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
   real(DbKi),                    intent(in   ) :: t          !< Current simulation time in seconds
   real(ReKi),                    intent(in   ) :: dt         !< Time step
   type(FVW_InputType),           intent(in   ) :: u          !< Input at t
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
   integer(IntKi) :: iAge
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = "" 

   ! Compute "right hand side"
   CALL FVW_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )

   ! Update of positions
   x%r_NW(1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) = x%r_NW(1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) +  dt * dxdt%r_NW(1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings)
   x%r_FW(1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings) = x%r_FW(1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings) +  dt * dxdt%r_FW(1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings)

   ! Update of Gamma
   ! TODO, viscous diffusion, stretching

   call FVW_DestroyContState(dxdt, ErrStat, ErrMsg)

end subroutine FVW_Euler1
!----------------------------------------------------------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------------------------------------------------
!> This is a tight coupling routine for solving for the residual of the constraint state functions.
subroutine FVW_CalcConstrStateResidual( t, u, p, x, xd, z_guess, OtherState, m, z_out, ErrStat, ErrMsg, iLabel)
   real(DbKi),                    intent(in   )  :: t           !< Current simulation time in seconds
   type(FVW_InputType),           intent(in   )  :: u           !< Inputs at t
   type(FVW_ParameterType),       intent(in   )  :: p           !< Parameters
   type(FVW_ContinuousStateType), intent(in   )  :: x           !< Continuous states at t
   type(FVW_DiscreteStateType),   intent(in   )  :: xd          !< Discrete states at t
   type(FVW_ConstraintStateType), intent(in   )  :: z_guess     !< Constraint states at t (possibly a guess)
   type(FVW_OtherStateType),      intent(in   )  :: OtherState  !< Other states at t
   type(FVW_MiscVarType),         intent(inout)  :: m           !< Misc variables for optimization (not copied in glue code)
   type(FVW_ConstraintStateType), intent(  out)  :: z_out            !< Residual of the constraint state functions using
   integer(IntKi), intent(in) :: iLabel
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

   CALL Wings_ComputeCirculation(t, z_out%Gamma_LL, z_guess%Gamma_LL, u, p, x, m, ErrStat, ErrMsg, iLabel)

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
   type(FVW_ConstraintStateType),   intent(inout)  :: z           !< Constraint states at t
   type(FVW_OtherStateType),        intent(in   )  :: OtherState  !< Other states at t
   type(FVW_OutputType),            intent(inout)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                                  !!   nectivity information does not have to be recalculated)
   type(FVW_MiscVarType),           intent(inout)  :: m           !< Misc/optimization variables
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   ! Local variables
   integer(IntKi)                :: iSpan, iW
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   character(*), parameter       :: RoutineName = 'FVW_CalcOutput'

   ErrStat = ErrID_None
   ErrMsg  = ""

   print'(A,F10.3,A,L1,A,I0,A,I0)','CalcOutput     t:',t,'   ',m%FirstCall,'                                nNW:',m%nNW,' nFW:',m%nFW

   if (m%FirstCall) then
      print*,'>>> First Call of CalcOutput, calling panelling and constrstate'
      CALL Wings_Panelling(u%WingsMesh, p, m, ErrStat2, ErrMsg2); if(Failed()) return
      CALL Wings_ComputeCirculation(t, m%Gamma_LL, z%Gamma_LL, u, p, x, m, ErrStat, ErrMsg, 0) ! For plotting only
   else
      m%Gamma_LL = z%Gamma_LL ! For plotting only
   endif

   if (.not. allocated(y%Vind)) then
      !call AllocAry( y%Vind ,  3, p%nSpan+1, p%nWings, 'Induced velocity vector',  ErrStat2, ErrMsg2 ); ! TODO potentially nSpan+1 for AD15
      call AllocAry( y%Vind ,  3, p%nSpan+1, 3, 'Induced velocity vector',  ErrStat2, ErrMsg2 ); ! NOTE: temporary hack 3 blades
      if(Failed()) return
   endif
   ! Returned guessed locations where wind will be required
   CALL SetRequestedWindPoints(y%r_wind, x, p, m, ErrStat2, ErrMsg2 ); if(Failed()) return

   ! Induction on the lifting line control point
   ! Set m%Vind_LL
   m%Vind_LL=-9999.0_ReKi
   call LiftingLineInducedVelocities(p, x, 1, m, ErrStat2, ErrMsg2); if(Failed()) return

   ! Interpolation to AeroDyn radial station TODO TODO TODO
   y%Vind(1:3,:,:) = 0.0_ReKi
   do iW=1,p%nWings
      do iSpan=1,p%nSpan
         y%Vind(1:3,iSpan,iW) = m%Vind_LL(1:3,iSpan,iW)
      enddo
   enddo

   ! For plotting only
   m%Vtot_ll = m%Vind_LL + m%Vwnd_LL - m%Vstr_ll
   call print_mean_3d(m%Vind_LL,'Mean induced vel. LL')
   call print_mean_3d(m%Vtot_LL,'Mean relativevel. LL')

   ! We don't propagate the "Old"-> "New" if update states was not called once
   ! This is introduced since at init, CalcOutput is called before UpdateState
   if (.not. m%FirstCall) then 
       call PrepareNextTimeStep()
   endif


   ! --- Write to VTK 
   if (p%WrVTK==1) then
      if (m%FirstCall) then
         call MKDIR('vtk_out')
         call WrVTK_FVW(p, x, z, m, 'vtk_out/FVW', m%iStep, 9)
      else
         call WrVTK_FVW(p, x, z, m, 'vtk_out/FVW', m%iStep+1, 9)
      endif
   endif

contains

   subroutine PrepareNextTimeStep()
      ! --- Increase wake length if maximum not reached
      if (m%nNW==p%nNWMax) then ! a far wake exist
         m%nFW=min(m%nFW+1, p%nFWMax)
      endif
      m%nNW=min(m%nNW+1, p%nNWMax)

   end subroutine PrepareNextTimeStep

   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FVW_CalcOutput') 
      Failed =  ErrStat >= AbortErrLev
      !if (Failed) call CleanUp()
   end function Failed
   !==========================================================================
end subroutine FVW_CalcOutput


END MODULE FVW
