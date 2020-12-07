!>
!!
!! Abbreviations:
!!   - FVW: Free Vortex Wake
!!   - LL : Lifting Line
!!   - CP : Control point
!!   - NW : Near Wake
!!   - FW : Far Wake
!!
module FVW
   use NWTC_Library
   use FVW_Types
   use FVW_Subs
   use FVW_IO
   use FVW_Wings
   use FVW_BiotSavart
   use FVW_Tests
   use AirFoilInfo

   IMPLICIT NONE

   PRIVATE

   type(ProgDesc), parameter  :: FVW_Ver = ProgDesc( 'FVW', '', '' )

   public   :: FVW_Init             ! Initialization routine
   public   :: FVW_End

   public   :: FVW_CalcOutput
   public   :: FVW_UpdateStates

   ! parameter for deciding if enough time has elapsed (Wake calculation, and vtk output)
   real(DbKi), parameter      :: OneMinusEpsilon = 1 - 10000*EPSILON(1.0_DbKi)

contains

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the start of the simulation to perform initialization steps.
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
subroutine FVW_Init(AFInfo, InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )
   use OMP_LIB ! wrap with #ifdef _OPENMP if this causes an issue
   type(AFI_ParameterType),         intent(in   )  :: AFInfo(:)      !< The airfoil parameter data, temporary, for UA..
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
   character(len=1054) :: DirName

   ! Initialize variables for this routine
   ErrStat = ErrID_None
   ErrMsg  = ""
   UnEcho  = -1

   ! Initialize the NWTC Subroutine Library
   call NWTC_Init( EchoLibVer=.FALSE. )

   ! Display the module information
   call DispNVD( FVW_Ver )
   ! Display convenient info to screen, until this is one day displayed by OpenFAST
   call getcwd(DirName)
   call WrScr(' - Directory:         '//trim(DirName))
   call WrScr(' - RootName:          '//trim(InitInp%RootName))
#ifdef _OPENMP   
   call WrScr(' - Compiled with OpenMP')
   !$OMP PARALLEL default(shared)
   if (omp_get_thread_num()==0) then
        call WrScr('   Number of threads: '//trim(Num2LStr(omp_get_num_threads()))//'/'//trim(Num2LStr(omp_get_max_threads())))
   endif
   !$OMP END PARALLEL 
#else
   call WrScr(' - No OpenMP support')
#endif
   if (DEV_VERSION) then
      CALL FVW_RunTests(ErrStat2, ErrMsg2); if (Failed()) return
   endif

   ! Set Parameters and *Misc* from inputs
   CALL FVW_SetParametersFromInputs(InitInp, p, ErrStat2, ErrMsg2); if(Failed()) return

   ! Read and parse the input file here to get other parameters and info
   CALL FVW_ReadInputFile(InitInp%FVWFileName, p, InputFileData, ErrStat2, ErrMsg2); if(Failed()) return

   ! Trigger required before allocations
   p%nNWMax  = max(InputFileData%nNWPanels,0)+1          ! +1 since LL panel included in NW
   p%nFWMax  = max(InputFileData%nFWPanels,0)
   p%nFWFree = max(InputFileData%nFWPanelsFree,0)
   p%DTfvw   = InputFileData%DTfvw
   p%DTvtk   = InputFileData%DTvtk

   ! Initialize Misc Vars (may depend on input file)
   CALL FVW_InitMiscVars( p, m, ErrStat2, ErrMsg2 ); if(Failed()) return

   ! Move the InitInp%WingsMesh to u
   CALL MOVE_ALLOC( InitInp%WingsMesh, u%WingsMesh )     ! Move from InitInp to u

!NOTE: We do not have the windspeed until after the FVW initialization (IfW is not initialized until after AD15)
   ! Wind Speed hack, TODO temporary NOTE: it is still needed?
   m%Vwnd_LL(:,:,:)   = 0
   m%Vwnd_NW(:,:,:,:) = 0
   m%Vwnd_FW(:,:,:,:) = 0

   ! This mesh is passed in as a cousin of the BladeMotion mesh.
   CALL Wings_Panelling_Init(u%WingsMesh, p, m, ErrStat2, ErrMsg2); if(Failed()) return

   ! Set parameters from InputFileData (need Misc allocated)
   CALL FVW_SetParametersFromInputFile(InputFileData, p, m, ErrStat2, ErrMsg2); if(Failed()) return

   ! Initialize Misc Vars (after input file params)
   CALL FVW_InitMiscVarsPostParam( p, m, ErrStat2, ErrMsg2 ); if(Failed()) return

   ! Initialize States Vars
   CALL FVW_InitStates( x, p, ErrStat2, ErrMsg2 ); if(Failed()) return

   ! Initialize Constraints Vars
   CALL FVW_InitConstraint( z, p, m, ErrStat2, ErrMsg2 ); if(Failed()) return

   ! Panelling wings based on initial input mesh provided
   ! This mesh is now a cousin of the BladeMotion mesh from AD.
   CALL Wings_Panelling     (u%WingsMesh, p, m, ErrStat2, ErrMsg2); if(Failed()) return
   CALL FVW_InitRegularization(p, m, ErrStat2, ErrMsg2); if(Failed()) return
   CALL FVW_ToString(p, m) ! Print to screen

   ! Mapping NW and FW (purely for esthetics, and maybe wind) ! TODO, just points
   call Map_LL_NW(p, m, z, x, 1.0_ReKi, ErrStat2, ErrMsg2); if(Failed()) return
   call Map_NW_FW(p, m, z, x, ErrStat2, ErrMsg2); if(Failed()) return

   ! Initialize output
   CALL FVW_Init_Y( p, u, y, ErrStat2, ErrMsg2); if(Failed()) return

   ! Returned guessed locations where wind will be required
   CALL SetRequestedWindPoints(m%r_wind, x, p, m )
   ! Return anything in FVW_InitOutput that should be passed back to the calling code here


   ! --- UA 
   ! NOTE: quick and dirty since this should belong to AD
   interval = InitInp%DTAero ! important, UA, needs proper interval
   call UA_Init_Wrapper(AFInfo, InitInp, interval, p, x, xd, OtherState, m, ErrStat2, ErrMsg2); if (Failed()) return

   ! Framework types unused
   Interval = InitInp%DTAero
   OtherState%NULL = 0
   xd%NULL         = 0
   InitOut%NULL    = 0
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
   integer(IntKi)          :: nMax           ! Total number of wind points possible
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
   m%VTKStep   = -1            ! Counter of VTK outputs
   m%VTKlastTime = -HUGE(1.0_DbKi)
   m%tSpent    = 0

   call AllocAry( m%LE      ,  3  ,  p%nSpan+1  , p%nWings, 'Leading Edge Points', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%LE = -999999_ReKi;
   call AllocAry( m%TE      ,  3  ,  p%nSpan+1  , p%nWings, 'TrailingEdge Points', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%TE = -999999_ReKi;
   call AllocAry( m%s_LL          ,  p%nSpan+1  , p%nWings, 'Spanwise coord LL  ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%s_LL= -999999_ReKi;
   call AllocAry( m%chord_LL      ,  p%nSpan+1  , p%nWings, 'Chord on LL        ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%chord_LL= -999999_ReKi;
   call AllocAry( m%PitchAndTwist ,  p%nSpan+1  , p%nWings, 'Pitch and twist    ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%PitchAndTwist= -999999_ReKi;
   call AllocAry( m%alpha_LL,        p%nSpan    , p%nWings, 'Wind on CP ll      ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%alpha_LL= -999999_ReKi;
   call AllocAry( m%Vreln_LL,        p%nSpan    , p%nWings, 'Wind on CP ll      ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%Vreln_LL = -999999_ReKi;
   ! Variables at control points/elements
   call AllocAry( m%Gamma_LL,        p%nSpan    , p%nWings, 'Lifting line Circulation', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%Gamma_LL = -999999_ReKi;
   call AllocAry( m%chord_CP_LL   ,  p%nSpan    , p%nWings, 'Chord on CP LL     ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%chord_CP_LL= -999999_ReKi;
   call AllocAry( m%s_CP_LL       ,  p%nSpan    , p%nWings, 'Spanwise coord CPll', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%s_CP_LL= -999999_ReKi;
   call AllocAry( m%CP_LL   , 3   ,  p%nSpan    , p%nWings, 'Control points LL  ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%CP_LL= -999999_ReKi;
   call AllocAry( m%Tang    , 3   ,  p%nSpan    , p%nWings, 'Tangential vector  ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%Tang= -999999_ReKi;
   call AllocAry( m%Norm    , 3   ,  p%nSpan    , p%nWings, 'Normal     vector  ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%Norm= -999999_ReKi;
   call AllocAry( m%Orth    , 3   ,  p%nSpan    , p%nWings, 'Orthogonal vector  ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%Orth= -999999_ReKi;
   call AllocAry( m%dl      , 3   ,  p%nSpan    , p%nWings, 'Orthogonal vector  ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%dl= -999999_ReKi;
   call AllocAry( m%Area    ,        p%nSpan    , p%nWings, 'LL Panel area      ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%Area = -999999_ReKi;
   call AllocAry( m%diag_LL ,        p%nSpan    , p%nWings, 'LL Panel diagonals ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%diag_LL = -999999_ReKi;
   call AllocAry( m%Vind_LL , 3   ,  p%nSpan    , p%nWings, 'Vind on CP ll      ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%Vind_LL= -999999_ReKi;
   call AllocAry( m%Vtot_LL , 3   ,  p%nSpan    , p%nWings, 'Vtot on CP ll      ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%Vtot_LL= -999999_ReKi;
   call AllocAry( m%Vstr_LL , 3   ,  p%nSpan    , p%nWings, 'Vstr on CP ll      ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%Vstr_LL= -999999_ReKi;
   call AllocAry( m%Vwnd_LL , 3   ,  p%nSpan    , p%nWings, 'Wind on CP ll      ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%Vwnd_LL= -999999_ReKi;
   ! Variables at panels points
   call AllocAry( m%r_LL    , 3   ,  p%nSpan+1  , 2        ,  p%nWings, 'Lifting Line Panels', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%r_LL= -999999_ReKi;
   call AllocAry( m%Vwnd_NW , 3   ,  p%nSpan+1  ,p%nNWMax+1,  p%nWings, 'Wind on NW ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%Vwnd_NW= -999_ReKi;
   call AllocAry( m%Vwnd_FW , 3   ,  FWnSpan+1  ,p%nFWMax+1,  p%nWings, 'Wind on FW ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%Vwnd_FW= -999_ReKi;
   call AllocAry( m%Vind_NW , 3   ,  p%nSpan+1  ,p%nNWMax+1,  p%nWings, 'Vind on NW ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%Vind_NW= -999_ReKi;
   call AllocAry( m%Vind_FW , 3   ,  FWnSpan+1  ,p%nFWMax+1,  p%nWings, 'Vind on FW ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%Vind_FW= -999_ReKi;
   call AllocAry( m%dxdt_NW , 3   ,  p%nSpan+1 , p%nNWMax+1,  p%nWings, 'NW dxdt'    , ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%dxdt_NW = -999999_ReKi;
   call AllocAry( m%dxdt_FW , 3   ,  FWnSpan+1 , p%nFWMax+1,  p%nWings, 'FW dxdt'    , ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%dxdt_FW = -999999_ReKi;
   ! Variables for optimizing outputs at blade nodes
   call AllocAry( m%BN_UrelWind_s, 3, p%nSpan+1 , p%nWings, 'Relative wind in section coordinates',   ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%BN_UrelWind_s= -999999_ReKi;
   call AllocAry( m%BN_AxInd   ,      p%nSpan+1 , p%nWings, 'Axial induction',                        ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%BN_AxInd     = -999999_ReKi;
   call AllocAry( m%BN_TanInd  ,      p%nSpan+1 , p%nWings, 'Tangential induction',                   ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%BN_TanInd    = -999999_ReKi;
   call AllocAry( m%BN_Vrel    ,      p%nSpan+1 , p%nWings, 'Relative velocity',                      ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%BN_Vrel      = -999999_ReKi;
   call AllocAry( m%BN_alpha   ,      p%nSpan+1 , p%nWings, 'Angle of attack',                        ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%BN_alpha     = -999999_ReKi;
   call AllocAry( m%BN_phi     ,      p%nSpan+1 , p%nWings, 'angle between the plane local wind dir', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%BN_phi       = -999999_ReKi;
   call AllocAry( m%BN_Re      ,      p%nSpan+1 , p%nWings, 'Reynolds number',                        ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%BN_Re        = -999999_ReKi;
   call AllocAry( m%BN_Cl_Static ,    p%nSpan+1 , p%nWings, 'Coefficient lift - no UA',               ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%BN_Cl_Static = -999999_ReKi;
   call AllocAry( m%BN_Cd_Static ,    p%nSpan+1 , p%nWings, 'Coefficient drag - no UA',               ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%BN_Cd_Static = -999999_ReKi;
   call AllocAry( m%BN_Cm_Static ,    p%nSpan+1 , p%nWings, 'Coefficient moment - no UA',             ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%BN_Cm_Static = -999999_ReKi;
   call AllocAry( m%BN_Cl        ,    p%nSpan+1 , p%nWings, 'Coefficient lift - with UA',             ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%BN_Cl        = -999999_ReKi;
   call AllocAry( m%BN_Cd        ,    p%nSpan+1 , p%nWings, 'Coefficient drag - with UA',             ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%BN_Cd        = -999999_ReKi;
   call AllocAry( m%BN_Cm        ,    p%nSpan+1 , p%nWings, 'Coefficient moment - with UA',           ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%BN_Cm        = -999999_ReKi;
   call AllocAry( m%BN_Cx        ,    p%nSpan+1 , p%nWings, 'Coefficient normal (to plane)',          ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%BN_Cx        = -999999_ReKi;
   call AllocAry( m%BN_Cy        ,    p%nSpan+1 , p%nWings, 'Coefficient tangential (to plane)',      ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%BN_Cy        = -999999_ReKi;


   ! Wind request points
   nMax = 0
   nMax = nMax +  p%nSpan                   * p%nWings   ! Lifting line Control Points
   nMax = nMax + (p%nSpan+1) * (p%nNWMax+1) * p%nWings   ! Nearwake points
   nMax = nMax + (FWnSpan+1) * (p%nFWMax+1) * p%nWings   ! Far wake points
   call AllocAry( m%r_wind, 3, nMax, 'Requested wind points', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName )
   m%r_wind = 0.0_ReKi     ! set to zero so InflowWind can shortcut calculations
   m%OldWakeTime = -HUGE(1.0_DbKi)
   ! Temporary UA
   call AllocAry( m%Vwnd_ND, 3, p%nSpan+1, p%nWings, 'Vwnd_ND', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%TE = -999999_ReKi;

end subroutine FVW_InitMiscVars
! ==============================================================================
subroutine FVW_InitMiscVarsPostParam( p, m, ErrStat, ErrMsg )
   type(FVW_ParameterType),         intent(in   )  :: p              !< Parameters
   type(FVW_MiscVarType),           intent(inout)  :: m              !< Initial misc/optimization variables
   integer(IntKi),                  intent(  out)  :: ErrStat        !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None
   integer(IntKi)          :: ErrStat2       ! temporary error status of the operation
   character(ErrMsgLen)    :: ErrMsg2        ! temporary error message
   character(*), parameter :: RoutineName = 'FVW_InitMiscVarsPostParam'
   integer(IntKi) :: nSeg, nSegP, nSegNW  !< Total number of segments after packing
   integer(IntKi) :: nCPs                 !< Total number of control points
   logical :: bMirror
   ErrStat = ErrID_None
   ErrMsg  = ""
   ! --- Counting maximum number of segments and Control Points expected for the whole simulation
   call CountSegments(p, p%nNWMax, p%nFWMax, 1, nSeg, nSegP, nSegNW)
   nCPs = CountCPs(p, p%nNWMax, p%nFWFree)

   bMirror = p%ShearModel==idShearMirror ! Whether or not we mirror the vorticity wrt ground
   if (bMirror) then
      nSeg  = nSeg*2
      nSegP = nSegP*2
   endif
   call AllocAry( m%SegConnct, 4, nSeg , 'SegConnct' , ErrStat2, ErrMsg2 );call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%SegConnct = -999;
   call AllocAry( m%SegPoints, 3, nSegP, 'SegPoints' , ErrStat2, ErrMsg2 );call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%SegPoints = -999999_ReKi;
   call AllocAry( m%SegGamma ,    nSeg,  'SegGamma'  , ErrStat2, ErrMsg2 );call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%SegGamma  = -999999_ReKi;
   call AllocAry( m%SegEpsilon,   nSeg,  'SegEpsilon', ErrStat2, ErrMsg2 );call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%SegEpsilon= -999999_ReKi;

   call AllocAry( m%CPs      , 3,  nCPs, 'CPs'       , ErrStat2, ErrMsg2 );call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%CPs= -999999_ReKi;
   call AllocAry( m%Uind     , 3,  nCPs, 'Uind'      , ErrStat2, ErrMsg2 );call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%Uind= -999999_ReKi;

end subroutine FVW_InitMiscVarsPostParam
! ==============================================================================
subroutine FVW_InitStates( x, p, ErrStat, ErrMsg )
   type(FVW_ContinuousStateType),   intent(  out)  :: x              !< States
   type(FVW_ParameterType),         intent(in   )  :: p              !< Parameters
   integer(IntKi),                  intent(  out)  :: ErrStat        !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None
   integer(IntKi)          :: ErrStat2       ! temporary error status of the operation
   character(ErrMsgLen)    :: ErrMsg2        ! temporary error message
   character(*), parameter :: RoutineName = 'FVW_InitMiscVars'
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

   call AllocAry( x%Gamma_NW,    p%nSpan   , p%nNWMax  , p%nWings, 'NW Panels Circulation', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitStates' ); 
   call AllocAry( x%Gamma_FW,    FWnSpan   , p%nFWMax  , p%nWings, 'FW Panels Circulation', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitStates' ); 
   ! set x%r_NW and x%r_FW to (0,0,0) so that InflowWind can shortcut the calculations
   call AllocAry( x%r_NW    , 3, p%nSpan+1 , p%nNWMax+1, p%nWings, 'NW Panels Points'     , ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitStates' );
   call AllocAry( x%r_FW    , 3, FWnSpan+1 , p%nFWMax+1, p%nWings, 'FW Panels Points'     , ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitStates' );
   !if (DEV_VERSION) then
   !   x%r_NW     = -9999999_ReKi;
   !   x%r_FW     = -9999999_ReKi;
   !   x%Gamma_NW = -999999_ReKi;
   !   x%Gamma_FW = -999999_ReKi;
   !else
   x%r_NW     = 0.0_ReKi
   x%r_FW     = 0.0_ReKi
   x%Gamma_NW = 0.0_ReKi ! First call of calcoutput, states might not be set 
   x%Gamma_FW = 0.0_ReKi ! NOTE, these values might be mapped from z%Gamma_LL at init
   !endif
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
   call AllocAry( z%Gamma_LL,  p%nSpan, p%nWings, 'Lifting line Circulation', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitConstraint' );
   !z%Gamma_LL = -999999_ReKi
   z%Gamma_LL = 0.0_ReKi

   if (ErrStat >= AbortErrLev) return
   if(.false.) print*,m%nNW ! unused var for now
end subroutine FVW_InitConstraint
! ==============================================================================
subroutine FVW_Init_Y( p, u, y, ErrStat, ErrMsg )
   type(FVW_ParameterType),         intent(in   )  :: p              !< Parameters
   type(FVW_InputType),             intent(inout)  :: u              !< An initial guess for the input; input mesh must be defined
   type(FVW_OutputType),            intent(  out)  :: y              !< Constraints
   integer(IntKi),                  intent(  out)  :: ErrStat        !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None
   integer(IntKi)          :: nMax           ! Total number of wind points possible
   integer(IntKi)          :: ErrStat2       ! temporary error status of the operation
   character(ErrMsgLen)    :: ErrMsg2        ! temporary error message
   character(*), parameter :: RoutineName = 'FVW_Init_Y'
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""
   !
   nMax = 0
   nMax = nMax +  p%nSpan                   * p%nWings   ! Lifting line Control Points
   nMax = nMax + (p%nSpan+1) * (p%nNWMax+1) * p%nWings   ! Nearwake points
   nMax = nMax + (FWnSpan+1) * (p%nFWMax+1) * p%nWings   ! Far wake points

   call AllocAry( u%V_wind, 3, nMax, 'Wind Velocity at points', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName )
   call AllocAry( y%Vind ,  3, p%nSpan+1, p%nWings, 'Induced velocity vector',  ErrStat2, ErrMsg2 ); ! TODO potentially nSpan+1 for AD15
   !call AllocAry( y%Cl_KJ , 1, 1, 'Lift coefficient from circulation (Kutta-Joukowski)', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName )
   if (ErrStat >= AbortErrLev) return
   y%Vind   = 0.0_ReKi
   return
end subroutine FVW_Init_Y


! ==============================================================================
!> Setting parameters *and misc* from module inputs
SUBROUTINE FVW_SetParametersFromInputs( InitInp, p, ErrStat, ErrMsg )
   type(FVW_InitInputType),    intent(inout)  :: InitInp       !< Input data for initialization routine  (inout so we can use MOVE_ALLOC)
   type(FVW_ParameterType),    intent(inout) :: p             !< Parameters
   integer(IntKi),             intent(  out) :: ErrStat       !< Error status of the operation
   character(*),               intent(  out) :: ErrMsg        !< Error message if ErrStat /= ErrID_None
   ! Local variables
   character(1024)          :: rootDir, baseName  ! Simulation root dir and basename
   !integer(IntKi)          :: ErrStat2
   !character(ErrMsgLen)    :: ErrMsg2
   character(*), parameter :: RoutineName = 'FVW_SetParametersFromInputs'
   ErrStat = ErrID_None
   ErrMsg  = ""
   ! 
   p%nWings       = InitInp%NumBlades
   p%nSpan        = InitInp%numBladeNodes-1 ! NOTE: temporary limitation, all wings have the same nspan
   p%DTaero       = InitInp%DTaero          ! AeroDyn Time step
   p%KinVisc      = InitInp%KinVisc         ! Kinematic air viscosity
   p%RootName     = InitInp%RootName        ! Rootname for outputs
   call GetPath( p%RootName, rootDir, baseName ) 
   p%VTK_OutFileRoot = trim(rootDir) // 'vtk_fvw'  ! Directory for VTK outputs
   p%VTK_OutFileBase = trim(rootDir) // 'vtk_fvw' // PathSep // trim(baseName) ! Basename for VTK files
   ! Set indexing to AFI tables -- this is set from the AD15 calling code.
   call AllocAry(p%AFindx,size(InitInp%AFindx,1),size(InitInp%AFindx,2),'AFindx',ErrStat,ErrMsg)
   p%AFindx = InitInp%AFindx     ! Copying in case AD15 still needs these

   ! Set the Chord values
   call move_alloc(InitInp%Chord, p%Chord)

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
   p%FWShedVorticity      = InputFileData%FWShedVorticity
   p%DiffusionMethod      = InputFileData%DiffusionMethod
   p%RegFunction          = InputFileData%RegFunction
   p%RegDeterMethod       = InputFileData%RegDeterMethod
   p%WakeRegMethod        = InputFileData%WakeRegMethod
   p%WakeRegParam         = InputFileData%WakeRegParam
   p%WingRegParam         = InputFileData%WingRegParam
   p%CoreSpreadEddyVisc   = InputFileData%CoreSpreadEddyVisc
   p%ShearModel           = InputFileData%ShearModel
   p%TwrShadowOnWake      = InputFileData%TwrShadowOnWake
   p%VelocityMethod       = InputFileData%VelocityMethod
   p%TreeBranchFactor     = InputFileData%TreeBranchFactor
   p%PartPerSegment       = InputFileData%PartPerSegment
   p%WrVTK                = InputFileData%WrVTK
   p%VTKBlades            = min(InputFileData%VTKBlades,p%nWings) ! Note: allowing it to be negative for temporary hack
   p%VTKCoord             = InputFileData%VTKCoord

   if (allocated(p%PrescribedCirculation)) deallocate(p%PrescribedCirculation)
   if (InputFileData%CirculationMethod==idCircPrescribed) then 
      call AllocAry( p%PrescribedCirculation,  p%nSpan, 'Prescribed Circulation', ErrStat2, ErrMsg2 ); call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_SetParameters' );    p%PrescribedCirculation = -999999_ReKi;
      if (.not. allocated(m%s_CP_LL)) then
         ErrMsg  = 'Spanwise coordinate not allocated.'
         ErrStat = ErrID_Fatal
         return
      endif
      call ReadAndInterpGamma(trim(InputFileData%CirculationFile), m%s_CP_LL(1:p%nSpan,1), m%s_LL(p%nSpan+1,1), p%PrescribedCirculation, ErrStat2, ErrMsg2)
      call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_SetParameters' ); 
   endif

end subroutine FVW_SetParametersFromInputFile

subroutine FVW_ToString(p,m)
   type(FVW_ParameterType), intent(in)       :: p !< Parameters
   type(FVW_MiscVarType),      intent(inout) :: m !< Misc
   if (DEV_VERSION) then
      print*,'-----------------------------------------------------------------------------------------'
      if(.false.) print*,m%nNW ! unused var for now
   endif
end subroutine FVW_ToString

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
subroutine FVW_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )

   type(FVW_InputType),             intent(inout)  :: u(:)        !< System inputs
   type(FVW_ParameterType),         intent(inout)  :: p           !< Parameters
   type(FVW_ContinuousStateType),   intent(inout)  :: x           !< Continuous states
   type(FVW_DiscreteStateType),     intent(inout)  :: xd          !< Discrete states
   type(FVW_ConstraintStateType),   intent(inout)  :: z           !< Constraint states
   type(FVW_OtherStateType),        intent(inout)  :: OtherState  !< Other states
   type(FVW_OutputType),            intent(inout)  :: y           !< System outputs
   type(FVW_MiscVarType),           intent(inout)  :: m           !< Misc/optimization variables
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   integer(IntKi) :: i

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""
   ! Place any last minute operations or calculations here:
   ! Close files here:
   ! Destroy the input data:
   do i=1,size(u)
      call FVW_DestroyInput( u(i), ErrStat, ErrMsg )
   enddo

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

#ifdef UA_OUTS
   CLOSE(69)
#endif
end subroutine FVW_End



!----------------------------------------------------------------------------------------------------------------------------------
!> Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete and other states.
!! Continuous, constraint, discrete, and other states are updated for t + Interval
subroutine FVW_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, AFInfo, m, errStat, errMsg )
!..................................................................................................................................
   real(DbKi),                      intent(in   )  :: t           !< Current simulation time in seconds
   integer(IntKi),                  intent(in   )  :: n           !< Current simulation time step n = 0,1,...
   type(FVW_InputType),             intent(inout)  :: u(:)        !< Inputs at utimes (out only for mesh record-keeping in ExtrapInterp routine)
   real(DbKi),                      intent(in   )  :: utimes(:)   !< Times associated with u(:), in seconds
   type(FVW_ParameterType),         intent(in   )  :: p           !< Parameters
   type(FVW_ContinuousStateType),   intent(inout)  :: x           !< Input: Continuous states at t; Output: at t+DTaero
   type(FVW_DiscreteStateType),     intent(inout)  :: xd          !< Input: Discrete states at t;   Output: at t+DTaero
   type(FVW_ConstraintStateType),   intent(inout)  :: z           !< Input: Constraint states at t; Output: at t+DTaero
   type(FVW_OtherStateType),        intent(inout)  :: OtherState  !< Input: Other states at t;      Output: at t+DTaero
   type(AFI_ParameterType),         intent(in   )  :: AFInfo(:)   !< The airfoil parameter data
   type(FVW_MiscVarType),           intent(inout)  :: m           !< Misc/optimization variables
   integer(IntKi),                  intent(  out)  :: errStat     !< Error status of the operation
   character(*),                    intent(  out)  :: errMsg      !< Error message if ErrStat /= ErrID_None
   ! local variables
   type(FVW_InputType)           :: uInterp                                                            ! Interpolated/Extrapolated input
   integer(IntKi)                :: ErrStat2                                                           ! temporary Error status
   character(ErrMsgLen)          :: ErrMsg2                                                            ! temporary Error message
   type(FVW_ConstraintStateType) :: z_guess                                                                              ! <
   integer(IntKi) :: nP, nFWEff
   integer, dimension(8) :: time1, time2, time_diff
   real(ReKi) :: ShedScale !< Scaling factor for shed vorticity (for sub-cycling), 1 if no subcycling
   logical :: bReevaluation

   ErrStat = ErrID_None
   ErrMsg  = ""

      
   ! --- Handling of time step, and time compared to previous call
   m%iStep = n
   ! Reevaluation: two repetitive calls starting from the same time, we will roll back the wake emission
   bReevaluation=.False.
   if (abs(t-m%OldWakeTime)<0.25_ReKi* p%DTaero) then
      bReevaluation=.True.
   endif
   ! Compute Induced wake effects only if time since last compute is > DTfvw
   if ( (( t - m%OldWakeTime ) >= p%DTfvw*OneMinusEpsilon) )  then
      m%OldWakeTime = t
      m%ComputeWakeInduced = .TRUE.    ! It's time to update the induced velocities from wake
   else
      m%ComputeWakeInduced = .FALSE.
   endif
   if (bReevaluation) then
      print*,'[INFO] FVW: Update States: reevaluation at the same starting time'
      call RollBackPreviousTimeStep() ! Cancel wake emission done in previous call
      m%ComputeWakeInduced = .TRUE.
   endif
   if (m%ComputeWakeInduced) then
      call date_and_time(values=time1)
   endif


   nP = p%nWings * (  (p%nSpan+1)*(m%nNW-1+2) +(FWnSpan+1)*(m%nFW+1) )
   nFWEff = min(m%nFW, p%nFWFree)
   ! --- Display some status to screen
!FIXME: this conflicts with the SimStatus WrOver from the FAST_Subs.f90.  Leaving out for now.
!     Ideally we put this into a log file.
!   if (mod(n,10)==0) print'(A,F10.3,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,F7.2,A)','FVW status - t:',t,'  n:',n,'  nNW:',m%nNW-1,'/',p%nNWMax-1,'  nFW:',nFWEff, '+',m%nFW-nFWEff,'=',m%nFW,'/',p%nFWMax,'  nP:',nP,'  spent:', m%tSpent, 's'
   if (DEV_VERSION)  print'(A,F10.3,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,F7.2,A,L1)','FVW status - t:',t,'  n:',n,'  nNW:',m%nNW-1,'/',p%nNWMax-1,'  nFW:',nFWEff, '+',m%nFW-nFWEff,'=',m%nFW,'/',p%nFWMax,'  nP:',nP,'  spent:', m%tSpent, 's Comp:',m%ComputeWakeInduced

   ! --- Evaluation at t
   ! Inputs at t
   call FVW_CopyInput( u(2), uInterp, MESH_NEWCOPY, ErrStat2, ErrMsg2); if(Failed()) return
   call FVW_Input_ExtrapInterp(u(1:size(utimes)),utimes(:),uInterp,t, ErrStat2, ErrMsg2); if(Failed()) return
   call Wings_Panelling(uInterp%WingsMesh, p, m, ErrStat2, ErrMsg2); if(Failed()) return
   
   ! Distribute the Wind we requested to Inflow wind to storage Misc arrays
   CALL DistributeRequestedWind(u(1)%V_wind, p, m)

   ! --- Solve for circulation at t
   ! Returns: z%Gamma_LL (at t)
   call AllocAry( z_guess%Gamma_LL,  p%nSpan, p%nWings, 'Lifting line Circulation', ErrStat, ErrMsg );
   z_guess%Gamma_LL = m%Gamma_LL
   call FVW_CalcConstrStateResidual(t, uInterp, p, x, xd, z_guess, OtherState, m, z, AFInfo, ErrStat2, ErrMsg2, 1); if(Failed()) return

   call UA_UpdateState_Wrapper(AFInfo,t, n, uInterp, p, x, xd, OtherState, m, ErrStat2, ErrMsg2); if(Failed()) return

   ! Map circulation and positions between LL and NW  and then NW and FW
   ! Changes: x only
   ShedScale = 1.0_ReKi
   call Map_LL_NW(p, m, z, x, ShedScale, ErrStat2, ErrMsg2); if(Failed()) return
   call Map_NW_FW(p, m, z, x, ErrStat2, ErrMsg2); if(Failed()) return
   !call print_x_NW_FW(p, m, x,'Map_')

   ! --- Integration between t and t+DTaero
   ! NOTE: when sub-cycling, the previous convection velocity is used
   ! If dtfvw = n dtaero, we assume xdot_local dtaero = xdot_stored * dtfvw/n
   if (p%IntMethod .eq. idEuler1) then 
     call FVW_Euler1( t, uInterp, p, x, xd, z, OtherState, m, ErrStat2, ErrMsg2); if(Failed()) return
   !elseif (p%IntMethod .eq. idRK4) then 
   !   call FVW_RK4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
   !elseif (p%IntMethod .eq. idAB4) then
   !   call FVW_AB4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
   !elseif (p%IntMethod .eq. idABM4) then
   !   call FVW_ABM4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
   else  
      call SetErrStat(ErrID_Fatal,'Invalid time integration method:'//Num2LStr(p%IntMethod),ErrStat,ErrMsg,'FVW_UpdateState') 
   end IF
   !call print_x_NW_FW(p, m, x,'Conv')

   if (m%ComputeWakeInduced) then
      ! We extend the wake length, i.e. we emit a new panel of vorticity at the TE
      ! NOTE: this will be rolled back if UpdateState is called at the same starting time again
      call PrepareNextTimeStep()
      ! --- t+DTaero
      ! Propagation/creation of new layer of panels
      call PropagateWake(p, m, z, x, ErrStat2, ErrMsg2); if(Failed()) return
      !call print_x_NW_FW(p, m, x,'Prop_')
   endif

   ! Inputs at t+DTaero
   call FVW_Input_ExtrapInterp(u(1:size(utimes)),utimes,uInterp,t+p%DTaero, ErrStat2, ErrMsg2); if(Failed()) return

   ! Panelling wings based on input mesh at t+p%DTaero
   call Wings_Panelling(uInterp%WingsMesh, p, m, ErrStat2, ErrMsg2); if(Failed()) return

   ! Updating positions of first NW and FW panels (Circulation also updated but irrelevant)
   ! Changes: x only
   ShedScale = (t+p%DTaero - m%OldWakeTime)/p%DTfvw
   call Map_LL_NW(p, m, z, x, ShedScale, ErrStat2, ErrMsg2); if(Failed()) return
   call Map_NW_FW(p, m, z, x, ErrStat2, ErrMsg2); if(Failed()) return

   !call print_x_NW_FW(p, m, x,'Map2')

   ! --- Solve for circulation at t+p%DTaero
   ! Returns: z%Gamma_LL (at t+p%DTaero)
   z_guess%Gamma_LL = z%Gamma_LL ! We use as guess the circulation from the previous time step (see above)
   call FVW_CalcConstrStateResidual(t+p%DTaero, uInterp, p, x, xd, z_guess, OtherState, m, z, AFInfo, ErrStat2, ErrMsg2, 2); if(Failed()) return
!    print*,'US: z_Gamma',x%Gamma_NW(1,1,1)
!    print*,'US: x_Gamma',z%Gamma_LL(1,1)

   ! Updating circulation of near wake panel (and position but irrelevant)
   ! Changes: x only
   call Map_LL_NW(p, m, z, x, ShedScale, ErrStat2, ErrMsg2); if(Failed()) return
   call Map_NW_FW(p, m, z, x, ErrStat2, ErrMsg2); if(Failed()) return
   !call print_x_NW_FW(p, m, x,'Map3')

   ! --- Fake handling of ground effect
   call FakeGroundEffect(p, x, m, ErrStat, ErrMsg)

   ! set the wind points required for t+p%DTaero timestep
   CALL SetRequestedWindPoints(m%r_wind, x, p, m)

   if (m%FirstCall) then
      m%FirstCall=.False.
   endif
   if (m%ComputeWakeInduced) then
      ! Profiling of expensive time step
      call date_and_time(values=time2)
      time_diff=time2-time1
      m%tSpent = time_diff(5)*3600+time_diff(6)*60 +time_diff(7)+0.001*time_diff(8)
   endif
   call FVW_DestroyConstrState(z_guess, ErrStat2, ErrMsg2); if(Failed()) return

contains
   subroutine PrepareNextTimeStep()
      ! --- Increase wake length if maximum not reached
      if (m%nNW==p%nNWMax) then ! a far wake exist
         m%nFW=min(m%nFW+1, p%nFWMax)
      endif
      m%nNW=min(m%nNW+1, p%nNWMax)
   end subroutine PrepareNextTimeStep

   subroutine RollBackPreviousTimeStep()
      ! --- Decrease wake length if maximum not reached
      if (m%nNW==p%nNWMax) then ! a far wake exist
         m%nFW=max(m%nFW-1, 0)
      endif
      m%nNW=max(m%nNW-1, 0)
   end subroutine RollBackPreviousTimeStep

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
   integer(IntKi)       :: ErrStat2       ! temporary error status of the operation
   character(ErrMsgLen) :: ErrMsg2        ! temporary error message
   integer(IntKi)       :: nFWEff ! Number of farwake panels that are free at current time step
   integer(IntKi)       :: i,j,k
   real(ReKi), dimension(3) :: VmeanFW, VmeanNW ! Mean velocity of the near wake and far wake

   ErrStat = ErrID_None
   ErrMsg  = ""

   call AllocAry( dxdt%r_NW , 3   ,  p%nSpan+1  ,p%nNWMax+1,  p%nWings, 'Wind on NW ', ErrStat2, ErrMsg2); dxdt%r_NW= -999999_ReKi;
   call AllocAry( dxdt%r_FW , 3   ,  FWnSpan+1  ,p%nFWMax+1,  p%nWings, 'Wind on FW ', ErrStat2, ErrMsg2); dxdt%r_FW= -999999_ReKi;
   if(Failed()) return

   ! Only calculate freewake after start time and if on a timestep when it should be calculated.
   if ((t>= p%FreeWakeStart)) then
      nFWEff = min(m%nFW, p%nFWFree)

      ! --- Compute Induced velocities on the Near wake and far wake based on the marker postions:
      ! (expensive N^2 call)
      ! In  : x%r_NW,    r%r_FW 
      ! Out:  m%Vind_NW, m%Vind_FW 
      call WakeInducedVelocities(p, x, m, ErrStat2, ErrMsg2); if(Failed()) return

      ! --- Mean induced velocity over the near wake (NW)
      VmeanNW(1:3)=0
      if (m%nNW >1) then
         do i=1,size(m%Vind_NW,4); do j=2,m%nNW+1; do k=1,size(m%Vind_NW,2); 
            VmeanNW(1:3) = VmeanNW(1:3) + m%Vind_NW(1:3, k, j, i)
         enddo; enddo; enddo; 
         VmeanNW(1:3) = VmeanNW(1:3) / (size(m%Vind_NW,4)*m%nNW*size(m%Vind_NW,2))
      endif
      ! --- Induced velocity over the free far wake (FWEff)
      VmeanFW(1:3)=0
      if (nFWEff >0) then
         do i=1,size(m%Vind_FW,4); do j=1,nFWEff; do k=1,size(m%Vind_FW,2); 
            VmeanFW(1:3) = VmeanFW(1:3) + m%Vind_FW(1:3, k, j, i)
         enddo; enddo; enddo; 
         VmeanFW(1:3) = VmeanFW(1:3) / (size(m%Vind_FW,4)*nFWEff*size(m%Vind_FW,2))
      else
         VmeanFW=VmeanNW
         ! Since we convect the first FW point, we need a reasonable velocity there 
         ! NOTE: mostly needed for sub-cycling and when no NW
         m%Vind_FW(1, 1:FWnSpan+1, 1, 1:p%nWings) = VmeanFW(1)
         m%Vind_FW(2, 1:FWnSpan+1, 1, 1:p%nWings) = VmeanFW(2)
         m%Vind_FW(3, 1:FWnSpan+1, 1, 1:p%nWings) = VmeanFW(3)
      endif

      ! --- Convecting non-free FW with a constant induced velocity (and free stream)
      m%Vind_FW(1, 1:FWnSpan+1, p%nFWFree+1:p%nFWMax+1, 1:p%nWings) = VmeanFW(1) !
      m%Vind_FW(2, 1:FWnSpan+1, p%nFWFree+1:p%nFWMax+1, 1:p%nWings) = VmeanFW(2) !
      m%Vind_FW(3, 1:FWnSpan+1, p%nFWFree+1:p%nFWMax+1, 1:p%nWings) = VmeanFW(3) !

      if (DEV_VERSION) then
         call print_mean_4d( m%Vind_NW(:,:, 1:m%nNW+1,:), 'Mean induced vel. NW')
         if (nFWEff>0) then
            call print_mean_4d( m%Vind_FW(:,:, 1:nFWEff ,:), 'Mean induced vel. FW')
         endif
         print'(A25,3F12.4)','MeanFW (non free)',VmeanFW
         call print_mean_4d( m%Vwnd_NW(:,:, 1:m%nNW+1,:), 'Mean wind vel.    NW')
         call print_mean_4d( m%Vwnd_FW(:,:, 1:nFWEff+1,:), 'Mean wind vel. FWEff')
         call print_mean_4d( m%Vwnd_FW(:,:, (p%nFWFree+1):m%nFW+1,:), 'Mean wind vel.    FWNF')
         call print_mean_4d( m%Vwnd_FW(:,:, 1:m%nFW+1,:), 'Mean wind vel.    FW')
      endif

      ! --- Vortex points are convected with the free stream and induced velocity
      dxdt%r_NW(1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) = m%Vwnd_NW(1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) +  m%Vind_NW(1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings)
      dxdt%r_FW(1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings) = m%Vwnd_FW(1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings) +  m%Vind_FW(1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings)
   else
      if(DEV_VERSION) then
         call print_mean_4d( m%Vwnd_NW(:,:,1:m%nNW+1,:), 'Mean wind vel.    NW')
         !call print_mean_4d( m%Vwnd_FW(:,:,1:m%nFW+1,:), 'Mean wind vel.    FW')
      endif

      ! --- Vortex points are convected with the free stream
      dxdt%r_NW(1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) = m%Vwnd_NW(1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) 
      dxdt%r_FW(1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings) = m%Vwnd_FW(1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings)
   endif
   ! First NW point does not convect (bound to LL)
   dxdt%r_NW(1:3, :, 1:iNWStart-1, :)=0
   ! First FW point always convects (even if bound to NW)
   ! This is done for subcycling
   !dxdt%r_FW(1:3, :, 1, :)=0
contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FVW_CalcContStateDeriv') 
      Failed =  ErrStat >= AbortErrLev
   end function Failed
end subroutine FVW_CalcContStateDeriv

!----------------------------------------------------------------------------------------------------------------------------------
subroutine FVW_Euler1( t, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
   real(DbKi),                    intent(in   ) :: t          !< Current simulation time in seconds
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
   type(FVW_ContinuousStateType) :: dxdt ! time derivatives of continuous states
   real(ReKi)     :: dt
   integer(IntKi)       :: ErrStat2      ! temporary error status of the operation
   character(ErrMsgLen) :: ErrMsg2       ! temporary error message
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = "" 

   dt = real(p%DTaero,ReKi) ! NOTE: this is DTaero not DTfvw since we integrate at each sub time step
   ! Compute "right hand side"
   if (m%ComputeWakeInduced) then
      CALL FVW_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dxdt, ErrStat2, ErrMsg2); if (Failed()) return
      ! Storage of convection velocity, purely for sub-cycling for now
      ! Since Euler1 is linear we use partial increments of dtaero<dtfvw
      m%dxdt_NW = dxdt%r_NW
      m%dxdt_FW = dxdt%r_FW
   endif

   if (DEV_VERSION) then
      ! Additional checks
      if (any(m%dxdt_NW(1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings)<-999)) then
         print*,'FVW_Euler1: Attempting to convect NW with a wrong velocity'
         STOP
      endif
      if ( m%nFW>0) then
         if (any(m%dxdt_FW(1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings)<-999)) then
            call print_x_NW_FW(p, m, x, 'STP')
            print*,'FVW_Euler1: Attempting to convect FW with a wrong velocity'
            STOP
         endif
      endif
   endif

   ! Update of positions
   x%r_NW(1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) = x%r_NW(1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) +  dt * m%dxdt_NW(1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings)
   if ( m%nFW>0) then
      x%r_FW(1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings) = x%r_FW(1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings) +  dt * m%dxdt_FW(1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings)
   endif
   ! Update of Gamma
   ! TODO, viscous diffusion, stretching

   call FVW_DestroyContState(dxdt, ErrStat2, ErrMsg2); if(Failed()) return
contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FVW_Euler1') 
      Failed =  ErrStat >= AbortErrLev
   end function Failed
end subroutine FVW_Euler1
!----------------------------------------------------------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------------------------------------------------
!> This is a tight coupling routine for solving for the residual of the constraint state functions.
subroutine FVW_CalcConstrStateResidual( t, u, p, x, xd, z_guess, OtherState, m, z_out, AFInfo, ErrStat, ErrMsg, iLabel)
   real(DbKi),                    intent(in   )  :: t           !< Current simulation time in seconds
   type(FVW_InputType),           intent(in   )  :: u           !< Inputs at t
   type(FVW_ParameterType),       intent(in   )  :: p           !< Parameters
   type(FVW_ContinuousStateType), intent(in   )  :: x           !< Continuous states at t
   type(FVW_DiscreteStateType),   intent(in   )  :: xd          !< Discrete states at t
   type(FVW_ConstraintStateType), intent(in   )  :: z_guess     !< Constraint states at t (possibly a guess)
   type(FVW_OtherStateType),      intent(in   )  :: OtherState  !< Other states at t
   type(FVW_MiscVarType),         intent(inout)  :: m           !< Misc variables for optimization (not copied in glue code)
   type(FVW_ConstraintStateType), intent(  out)  :: z_out       !< Residual of the constraint state functions using
   type(AFI_ParameterType),       intent(in   )  :: AFInfo(:)   !< The airfoil parameter data
   integer(IntKi),                intent(in   )  :: iLabel
   integer(IntKi),                intent(  OUT)  :: ErrStat     !< Error status of the operation
   character(*),                  intent(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""
   ! Solve for the residual of the constraint state functions here:
   !z%residual = 0.0_ReKi
   !z%Gamma_LL = 0.0_ReKi
   call AllocAry( z_out%Gamma_LL,  p%nSpan, p%nWings, 'Lifting line Circulation', ErrStat, ErrMsg );
   z_out%Gamma_LL = -999999_ReKi;

   CALL Wings_ComputeCirculation(t, z_out%Gamma_LL, z_guess%Gamma_LL, u, p, x, m, AFInfo, ErrStat, ErrMsg, iLabel)

end subroutine FVW_CalcConstrStateResidual

!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
!! This subroutine is used to compute the output channels (motions and loads) and place them in the WriteOutput() array.
!! The descriptions of the output channels are not given here. Please see the included OutListParameters.xlsx sheet for
!! for a complete description of each output parameter.
! NOTE: no matter how many channels are selected for output, all of the outputs are calculated
! All of the calculated output channels are placed into the m%AllOuts(:), while the channels selected for outputs are
! placed in the y%WriteOutput(:) array.
subroutine FVW_CalcOutput( t, u, p, x, xd, z, OtherState, AFInfo, y, m, ErrStat, ErrMsg )
   use FVW_VTK, only: set_vtk_coordinate_transform
   use FVW_VortexTools, only: interpextrap_cp2node
   real(DbKi),                      intent(in   )  :: t           !< Current simulation time in seconds
   type(FVW_InputType),             intent(in   )  :: u           !< Inputs at Time t
   type(FVW_ParameterType),         intent(in   )  :: p           !< Parameters
   type(FVW_ContinuousStateType),   intent(in   )  :: x           !< Continuous states at t
   type(FVW_DiscreteStateType),     intent(in   )  :: xd          !< Discrete states at t
!FIXME:TODO: AD15_CalcOutput has constraint states as intent(in) only. This is forcing me to store z in the AD15 miscvars for now.
   type(FVW_ConstraintStateType),   intent(in   )  :: z           !< Constraint states at t
   type(FVW_OtherStateType),        intent(in   )  :: OtherState  !< Other states at t
   type(AFI_ParameterType),         intent(in   )  :: AFInfo(:)   !< The airfoil parameter data
   type(FVW_OutputType),            intent(inout)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                                  !!   nectivity information does not have to be recalculated)
   type(FVW_MiscVarType),           intent(inout)  :: m           !< Misc/optimization variables
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   ! Local variables
   integer(IntKi)                :: iW, n, i0, i1, i2
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   character(*), parameter       :: RoutineName = 'FVW_CalcOutput'

   ErrStat = ErrID_None
   ErrMsg  = ""

   if (DEV_VERSION) then
      print'(A,F10.3,A,L1,A,I0,A,I0)','CalcOutput     t:',t,'   ',m%FirstCall,'                                nNW:',m%nNW,' nFW:',m%nFW
   endif

   ! Set the wind velocity at vortex
   CALL DistributeRequestedWind(u%V_wind, p, m)

   ! if we are on a correction step, CalcOutput may be called again with different inputs
   ! Compute m%Gamma_LL
   CALL Wings_ComputeCirculation(t, m%Gamma_LL, z%Gamma_LL, u, p, x, m, AFInfo, ErrStat2, ErrMsg2, 0); if(Failed()) return ! For plotting only


   ! Induction on the lifting line control point
   ! Set m%Vind_LL
   m%Vind_LL=-9999.0_ReKi
   call LiftingLineInducedVelocities(p, x, 1, m, ErrStat2, ErrMsg2); if(Failed()) return

   !  Induction on the mesh points (AeroDyn nodes)
   n=p%nSpan
   y%Vind(1:3,:,:) = 0.0_ReKi
   do iW=1,p%nWings
      ! --- Linear interpolation for interior points and extrapolations at boundaries
      call interpextrap_cp2node(m%s_CP_LL(:,iW), m%Vind_LL(1,:,iW), m%s_LL(:,iW), y%Vind(1,:,iW))
      call interpextrap_cp2node(m%s_CP_LL(:,iW), m%Vind_LL(2,:,iW), m%s_LL(:,iW), y%Vind(2,:,iW))
      call interpextrap_cp2node(m%s_CP_LL(:,iW), m%Vind_LL(3,:,iW), m%s_LL(:,iW), y%Vind(3,:,iW))
   enddo

   ! For plotting only
   m%Vtot_ll = m%Vind_LL + m%Vwnd_LL - m%Vstr_ll
   if (DEV_VERSION) then
      call print_mean_3d(m%Vind_LL,'Mean induced vel. LL')
      call print_mean_3d(m%Vtot_LL,'Mean relativevel. LL')
   endif

   ! --- Write to local VTK at fps requested
   if (p%WrVTK==1) then
      if (m%VTKStep==-1) then 
         m%VTKStep = 0 ! Has never been called, special handling for init
      else
         m%VTKStep = m%iStep+1 ! We use glue code step number for outputs
      endif
      if (m%FirstCall) then
         call MKDIR(p%VTK_OutFileRoot)
      endif
      if ( ( t - m%VTKlastTime ) >= p%DTvtk*OneMinusEpsilon )  then
         m%VTKlastTime = t
         if ((p%VTKCoord==2).or.(p%VTKCoord==3)) then
            ! Hub reference coordinates, for export only, ALL VTK Will be exported in this coordinate system!
            ! Note: hubOrientation and HubPosition are optional, but required for bladeFrame==TRUE
            call WrVTK_FVW(p, x, z, m, trim(p%VTK_OutFileBase)//'FVW_Hub', m%VTKStep, 9, bladeFrame=.TRUE.,  &
                     HubOrientation=real(u%HubOrientation,ReKi),HubPosition=real(u%HubPosition,ReKi))
         endif
         if ((p%VTKCoord==1).or.(p%VTKCoord==3)) then
            ! Global coordinate system, ALL VTK will be exported in global
            call WrVTK_FVW(p, x, z, m, trim(p%VTK_OutFileBase)//'FVW_Glb', m%VTKStep, 9, bladeFrame=.FALSE.)
         endif
      endif
   endif


contains

   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FVW_CalcOutput') 
      Failed =  ErrStat >= AbortErrLev
   end function Failed

end subroutine FVW_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------   
! --- UA related, should be merged with AeroDyn 
!----------------------------------------------------------------------------------------------------------------------------------   
subroutine UA_Init_Wrapper(AFInfo, InitInp, interval, p, x, xd, OtherState, m, ErrStat, ErrMsg )
   use UnsteadyAero, only: UA_Init, UA_TurnOff_param
   type(AFI_ParameterType),         intent(in   )  :: AFInfo(:)      !< The airfoil parameter data, temporary, for UA..
   type(FVW_InitInputType),         intent(inout)  :: InitInp     !< Input data for initialization routine  (inout so we can use MOVE_ALLOC)
   real(DbKi),                      intent(inout)  :: interval    !< time interval  
   type(FVW_ParameterType),         intent(inout)  :: p           !< Parameters
   type(FVW_ContinuousStateType),   intent(inout)  :: x           !< Initial continuous states
   type(FVW_DiscreteStateType),     intent(inout)  :: xd          !< Initial discrete states
   type(FVW_OtherStateType),        intent(inout)  :: OtherState  !< Initial other states
   type(FVW_MiscVarType),           intent(inout)  :: m           !< Initial misc/optimization variables
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   !
   type(UA_InitInputType) :: Init_UA_Data
   type(UA_InputType)     :: u_UA
   type(UA_InitOutputType):: InitOutData_UA
   integer                :: i,j
   integer(intKi)         :: ErrStat2           ! temporary Error status
   character(ErrMsgLen)   :: ErrMsg2
   ErrStat = ErrID_None
   ErrMsg  = ""

   m%UA_Flag=InitInp%UA_Flag
   ! --- Condensed version of BEMT_Init_Otherstate
   allocate ( OtherState%UA_Flag( InitInp%numBladeNodes, InitInp%NumBlades ), STAT = ErrStat2 )
   OtherState%UA_Flag=m%UA_Flag
   if ( m%UA_Flag ) then
      ErrMsg2='Unsteady aerodynamic (`AFAeroMod>1`) cannot be used with the free wake code (`WakeMod=3`) for now.'; ErrStat2=ErrID_Fatal;
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'UA_Init_Wrapper'); return

      ! ---Condensed version of "BEMT_Set_UA_InitData"
      allocate(Init_UA_Data%c(InitInp%numBladeNodes,InitInp%numBlades), STAT = errStat2)
      do j = 1,InitInp%NumBlades; do i = 1,InitInp%numBladeNodes;
            Init_UA_Data%c(i,j)      = p%chord(i,j) ! NOTE: InitInp chord move-allocd to p
      end do; end do
      Init_UA_Data%dt              = interval          
      Init_UA_Data%OutRootName     = ''
      Init_UA_Data%numBlades       = InitInp%NumBlades 
      Init_UA_Data%nNodesPerBlade  = InitInp%numBladeNodes
      Init_UA_Data%NumOuts         = 0
      Init_UA_Data%UAMod           = InitInp%UAMod  
      Init_UA_Data%Flookup         = InitInp%Flookup
      Init_UA_Data%a_s             = InitInp%a_s ! m/s  
      ! --- UA init
      call UA_Init( Init_UA_Data, u_UA, m%p_UA, x%UA, xd%UA, OtherState%UA, m%y_UA, m%m_UA, interval, InitOutData_UA, ErrStat2, ErrMsg2); if(Failed())return
      m%p_UA%ShedEffect=.False. !< Important, when coupling UA wih vortex code, shed vorticity is inherently accounted for
      ! --- Condensed version of "BEMT_CheckInitUA"
      do j = 1,InitInp%numBlades; do i = 1,InitInp%numBladeNodes; ! Loop over blades and nodes
         call UA_TurnOff_param(m%p_UA, AFInfo(p%AFindx(i,j)), ErrStat2, ErrMsg2)
         if (ErrStat2 /= ErrID_None) then
            call WrScr( 'Warning: Turning off Unsteady Aerodynamics because '//trim(ErrMsg2)//' BladeNode = '//trim(num2lstr(i))//', Blade = '//trim(num2lstr(j)) )
            OtherState%UA_Flag(i,j) = .false.
         end if
      end do; end do;
#ifdef UA_OUTS   
      CALL OpenFOutFile ( 69, 'Debug.UA.out', errStat, errMsg ); IF (ErrStat >= AbortErrLev) RETURN
      WRITE (69,'(/,A)')  'This output information was generated by FVW'// ' on '//CurDate()//' at '//CurTime()//'.'
      WRITE (69,'(:,A20)', ADVANCE='no' ) 'Time'
      do i=1,size(InitOutData_UA%WriteOutputHdr)
         WRITE (69,'(:,A20)', ADVANCE='no' )  trim(InitOutData_UA%WriteOutputHdr(i))
      end do  
      write (69,'(A)')    ' '
      WRITE (69,'(:,A20)', ADVANCE='no' ) '(s)'
      do i=1,size(InitOutData_UA%WriteOutputUnt)
         WRITE (69,'(:,A20)', ADVANCE='no' )  trim(InitOutData_UA%WriteOutputUnt(i))
      end do  
      write (69,'(A)')    ' '   
#endif
      call UA_DestroyInput( u_UA, ErrStat2, ErrMsg2 ); if(Failed())return
      call UA_DestroyInitInput( Init_UA_Data, ErrStat2, ErrMsg2 ); if(Failed())return
      call UA_DestroyInitOutput( InitOutData_UA, ErrStat2, ErrMsg2 ); if(Failed())return

      ! --- FVW specific
      if (p%CirculationMethod/=idCircPolarData) then 
         ErrMsg2='Unsteady aerodynamic (`AFAeroMod>1`) is only available with a circulation solving using profile data (`CircSolvingMethod=1`)'; ErrStat2=ErrID_Fatal;
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'UA_Init_Wrapper'); return
      endif
   endif
contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'UA_Init_Wrapper') 
      Failed =  ErrStat >= AbortErrLev
   end function Failed
end subroutine  UA_Init_Wrapper

subroutine UA_UpdateState_Wrapper(AFInfo, t, n, u, p, x, xd, OtherState, m, ErrStat, ErrMsg )
   use FVW_VortexTools, only: interpextrap_cp2node
   use UnsteadyAero, only: UA_UpdateStates, UA_TurnOff_input
   type(AFI_ParameterType),         intent(in   )  :: AFInfo(:)      !< The airfoil parameter data, temporary, for UA..
   real(DbKi),                      intent(in   )  :: t           !< Current simulation time in seconds
   integer(IntKi),                  intent(in   )  :: n           !< time step  
   type(FVW_ParameterType),         intent(in   )  :: p           !< Parameters
   type(FVW_InputType),             intent(in   )  :: u           !< Inputs
   type(FVW_ContinuousStateType),   intent(inout)  :: x           !< Initial continuous states
   type(FVW_DiscreteStateType),     intent(inout)  :: xd          !< Initial discrete states
   type(FVW_OtherStateType),        intent(inout)  :: OtherState  !< Initial other states
   type(FVW_MiscVarType),           intent(inout)  :: m           !< Initial misc/optimization variables
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   ! Local
   type(UA_InputType)     :: u_UA(1)
   REAL(DbKi)             :: uTimes(1)
   integer                :: i,j
   integer, parameter     :: k=1 ! index for u_UA (in case it is ever dimensioned differently)
   integer(intKi)         :: ErrStat2           ! temporary Error status
   character(ErrMsgLen)   :: ErrMsg2
   real(ReKi), dimension(:,:), allocatable :: Vind_node
   ErrStat  = ErrID_None
   ErrStat2 = ErrID_None
   ErrMsg   = ""
   ErrMsg2  = ""
   
   uTimes = t
   
   if (m%UA_Flag) then

      ! --- Induction on the lifting line control point
      ! NOTE: this is expensive since it's an output for FVW but here we have to use it for UA
      ! Set m%Vind_LL
      m%Vind_LL=-9999.0_ReKi
      call LiftingLineInducedVelocities(p, x, 1, m, ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'UA_UpdateState_Wrapper'); if (ErrStat >= AbortErrLev) return

      allocate(Vind_node(3,1:p%nSpan+1))

      ! --- Condensed version of BEMT_Update States
      do j = 1,p%nWings  
         ! Induced velocity at Nodes (NOTE: we rely on storage done when computing Circulation)
         if (m%nNW>1) then
            call interpextrap_cp2node(m%s_CP_LL(:,j), m%Vind_LL(1,:,j), m%s_LL(:,j), Vind_node(1,:))
            call interpextrap_cp2node(m%s_CP_LL(:,j), m%Vind_LL(2,:,j), m%s_LL(:,j), Vind_node(2,:))
            call interpextrap_cp2node(m%s_CP_LL(:,j), m%Vind_LL(3,:,j), m%s_LL(:,j), Vind_node(3,:))
         else
            Vind_node=0.0_ReKi
         endif

         do i = 1,p%nSpan+1 
            ! We only update the UnsteadyAero states if we have unsteady aero turned on for this node      
            if (OtherState%UA_Flag(i,j) .and. n > 0) then
               !! ....... compute inputs to UA ...........
               ! NOTE: To be consistent with CalcOutput we take Vwind_ND that was set using m%DisturbedInflow from AeroDyn.. 
               ! This is not clean, but done to be consistent, waiting for AeroDyn to handle UA
               call AlphaVrel_Generic(u%WingsMesh(j)%Orientation(1:3,1:3,i), u%WingsMesh(j)%TranslationVel(1:3,i),  Vind_node(:,i), m%Vwnd_ND(:,i,j), &
                                       p%KinVisc, p%Chord(i,j), u_UA(k)%U, u_UA(k)%alpha, u_UA(k)%Re)
               ! FIX ME: this is copied 3 times!!!!
               u_UA%v_ac(1) = sin(u_UA%alpha)*u_UA%U
               u_UA%v_ac(2) = cos(u_UA%alpha)*u_UA%U
               u_UA%omega = 0.0_ReKi ! FIX ME!!!! dot_product( u%BladeMotion(j)%RotationVel(   :,i), m%WithoutSweepPitchTwist(3,:,i,j) ) ! rotation of no-sweep-pitch coordinate system around z of the jth node in the kth blade
               u_UA(k)%UserProp = 0 ! u1%UserProp(i,j) ! TODO
               
               !! ....... check inputs to UA ...........
               call UA_TurnOff_input(m%p_UA, AFInfo(p%AFIndx(i,j)), u_UA(k), ErrStat2, ErrMsg2)
               if (ErrStat2 /= ErrID_None) then
                  OtherState%UA_Flag(i,j) = .FALSE.
                  call WrScr( 'Warning: Turning off dynamic stall due to '//trim(ErrMsg2)//' '//trim(NodeText(i,j)))
               else
                 ! COMPUTE: xd%UA, OtherState%UA
                 call UA_UpdateStates( i, j, t, n, u_UA, uTimes, m%p_UA, x%UA, xd%UA, OtherState%UA, AFInfo(p%AFIndx(i,j)), m%m_UA, ErrStat2, ErrMsg2 )
                 if (ErrStat2 /= ErrID_None) then
                     call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'UA_UpdateState_Wrapper'//trim(NodeText(i,j)))
                     call WrScr(trim(ErrMsg))
                    if (ErrStat >= AbortErrLev) return
                 end if
               end if
            end if
         end do
      end do
      call UA_DestroyInput( u_UA(k), ErrStat2, ErrMsg2 ); 
      deallocate(Vind_node)
   endif
   call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'UA_UpdateState_Wrapper')
contains 
   function NodeText(i,j)
      integer(IntKi), intent(in) :: i ! node number
      integer(IntKi), intent(in) :: j ! blade number
      character(25)              :: NodeText
      NodeText = '(nd:'//trim(num2lstr(i))//' bld:'//trim(num2lstr(j))//')'
   end function NodeText
end subroutine UA_UpdateState_Wrapper

end module FVW
