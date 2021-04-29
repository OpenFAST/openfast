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
   CALL FVW_ReadInputFile(InitInp%FVWFileName, p, m, InputFileData, ErrStat2, ErrMsg2); if(Failed()) return

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
   CALL FVW_InitRegularization(x, p, m, ErrStat2, ErrMsg2); if(Failed()) return
   CALL FVW_ToString(p, m) ! Print to screen

   ! Mapping NW and FW (purely for esthetics, and maybe wind) ! TODO, just points
   call Map_LL_NW(p, m, z, x, 1.0_ReKi, ErrStat2, ErrMsg2); if(Failed()) return
   call Map_NW_FW(p, m, z, x, ErrStat2, ErrMsg2); if(Failed()) return

   ! Initialize input guess and output
   CALL FVW_Init_U_Y( p, u, y, m, ErrStat2, ErrMsg2); if(Failed()) return

   ! Returned guessed locations where wind will be required
   CALL SetRequestedWindPoints(m%r_wind, x, p, m )
   ! Return anything in FVW_InitOutput that should be passed back to the calling code here

   ! --- UA 
   ! NOTE: quick and dirty since this should belong to AD
   interval = InitInp%DTAero ! important, gluecode and UA, needs proper interval
   call UA_Init_Wrapper(AFInfo, InitInp, interval, p, x, xd, OtherState, m, ErrStat2, ErrMsg2); if (Failed()) return

   ! Framework types unused
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
   integer(IntKi)          :: iGrid          ! 
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
   call AllocAry(m%iTip,  p%nWings, 'iTip',  ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%iTip  = -1;! Important init
   call AllocAry(m%iRoot, p%nWings, 'iRoot', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%iRoot = -1;! Important init

   call AllocAry( m%LE      ,  3  ,  p%nSpan+1  , p%nWings, 'Leading Edge Points', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%LE = -999999_ReKi;
   call AllocAry( m%TE      ,  3  ,  p%nSpan+1  , p%nWings, 'TrailingEdge Points', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%TE = -999999_ReKi;
   call AllocAry( m%PitchAndTwist ,  p%nSpan+1  , p%nWings, 'Pitch and twist    ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%PitchAndTwist= -999999_ReKi;
   call AllocAry( m%alpha_LL,        p%nSpan    , p%nWings, 'Wind on CP ll      ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%alpha_LL= -999999_ReKi;
   call AllocAry( m%Vreln_LL,        p%nSpan    , p%nWings, 'Wind on CP ll      ', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%Vreln_LL = -999999_ReKi;
   ! Variables at control points/elements
   call AllocAry( m%Gamma_LL,        p%nSpan    , p%nWings, 'Lifting line Circulation', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%Gamma_LL = -999999_ReKi;
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


   ! dxdt, to avoid realloc all the time, and storage for subcycling 
   call AllocAry( m%dxdt%r_NW , 3   ,  p%nSpan+1 , p%nNWMax+1,  p%nWings, 'r NW dxdt'  , ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%dxdt%r_NW = -999999_ReKi;
   call AllocAry( m%dxdt%r_FW , 3   ,  FWnSpan+1 , p%nFWMax+1,  p%nWings, 'r FW dxdt'  , ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%dxdt%r_FW = -999999_ReKi;
   call AllocAry( m%dxdt%Eps_NW, 3  ,  p%nSpan    ,p%nNWMax  ,  p%nWings, 'Eps NW dxdt', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%dxdt%Eps_NW = -999999_ReKi;
   call AllocAry( m%dxdt%Eps_FW, 3  ,  FWnSpan    ,p%nFWMax  ,  p%nWings, 'Eps FW dxdt', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%dxdt%Eps_FW = -999999_ReKi;

   ! Wind request points
   nMax = 0
   nMax = nMax +  p%nSpan                   * p%nWings   ! Lifting line Control Points
   nMax = nMax + (p%nSpan+1) * (p%nNWMax+1) * p%nWings   ! Nearwake points
   nMax = nMax + (FWnSpan+1) * (p%nFWMax+1) * p%nWings   ! Far wake points
   do iGrid=1,p%nGridOut
      nMax = nMax + m%GridOutputs(iGrid)%nx * m%GridOutputs(iGrid)%ny * m%GridOutputs(iGrid)%nz
      call AllocAry(m%GridOutputs(iGrid)%uGrid, 3, m%GridOutputs(iGrid)%nx,  m%GridOutputs(iGrid)%ny, m%GridOutputs(iGrid)%nz, 'uGrid', ErrStat2, ErrMsg2);
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName)
      m%GridOutputs(iGrid)%tLastOutput = -HUGE(1.0_DbKi)
   enddo
   call AllocAry( m%r_wind, 3, nMax, 'Requested wind points', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName )
   m%r_wind = 0.0_ReKi     ! set to zero so InflowWind can shortcut calculations
   m%OldWakeTime = -HUGE(1.0_DbKi)
   ! Wind set to 0. TODO check if -99999 works now
   !NOTE: We do not have the windspeed until after the FVW initialization (IfW is not initialized until after AD15)
   m%Vwnd_LL(:,:,:)   = 0
   m%Vwnd_NW(:,:,:,:) = 0
   m%Vwnd_FW(:,:,:,:) = 0

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
   call AllocAry( m%Sgmt%Connct, 4, nSeg , 'SegConnct' , ErrStat2, ErrMsg2 );call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%Sgmt%Connct = -999;
   call AllocAry( m%Sgmt%Points, 3, nSegP, 'SegPoints' , ErrStat2, ErrMsg2 );call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%Sgmt%Points = -999999_ReKi;
   call AllocAry( m%Sgmt%Gamma ,    nSeg,  'SegGamma'  , ErrStat2, ErrMsg2 );call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%Sgmt%Gamma  = -999999_ReKi;
   call AllocAry( m%Sgmt%Epsilon,   nSeg,  'SegEpsilon', ErrStat2, ErrMsg2 );call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%Sgmt%Epsilon= -999999_ReKi;
   m%Sgmt%nAct        = -1  ! Active segments
   m%Sgmt%nActP       = -1
   m%Sgmt%RegFunction = p%RegFunction

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
   call AllocAry( x%Eps_NW  , 3, p%nSpan   , p%nNWMax  , p%nWings, 'NW Panels Reg Param'  , ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitStates' );
   call AllocAry( x%Eps_FW  , 3, FWnSpan   , p%nFWMax  , p%nWings, 'FW Panels Reg Param'  , ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitStates' );
   ! set x%r_NW and x%r_FW to (0,0,0) so that InflowWind can shortcut the calculations
   call AllocAry( x%r_NW    , 3, p%nSpan+1 , p%nNWMax+1, p%nWings, 'NW Panels Points'     , ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitStates' );
   call AllocAry( x%r_FW    , 3, FWnSpan+1 , p%nFWMax+1, p%nWings, 'FW Panels Points'     , ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitStates' );
   if (ErrStat >= AbortErrLev) return
   x%r_NW     = 0.0_ReKi
   x%r_FW     = 0.0_ReKi
   x%Gamma_NW = 0.0_ReKi ! First call of calcoutput, states might not be set 
   x%Gamma_FW = 0.0_ReKi ! NOTE, these values might be mapped from z%Gamma_LL at init
   x%Eps_NW   = 0.001_ReKi 
   x%Eps_FW   = 0.001_ReKi 
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
!> Init/allocate inputs and outputs
subroutine FVW_Init_U_Y( p, u, y, m, ErrStat, ErrMsg )
   type(FVW_ParameterType),         intent(in   )  :: p              !< Parameters
   type(FVW_InputType),             intent(inout)  :: u              !< An initial guess for the input; input mesh must be defined
   type(FVW_MiscVarType),           intent(inout)  :: m              !< Initial misc/optimization variables
   type(FVW_OutputType),            intent(  out)  :: y              !< Constraints
   integer(IntKi),                  intent(  out)  :: ErrStat        !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None
   integer(IntKi)          :: nMax           ! Total number of wind points possible
   integer(IntKi)          :: ErrStat2       ! temporary error status of the operation
   character(ErrMsgLen)    :: ErrMsg2        ! temporary error message
   character(*), parameter :: RoutineName = 'FVW_Init_U_Y'
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""
   !
   nMax = 0
   nMax = nMax +  p%nSpan                   * p%nWings   ! Lifting line Control Points
   nMax = nMax + (p%nSpan+1) * (p%nNWMax+1) * p%nWings   ! Nearwake points
   nMax = nMax + (FWnSpan+1) * (p%nFWMax+1) * p%nWings   ! Far wake points

   call AllocAry( u%V_wind,   3, nMax, 'Wind Velocity at points', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName)
   call AllocAry( y%Vind ,    3, p%nSpan+1, p%nWings, 'Induced velocity vector', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName)
   call AllocAry( u%omega_z,     p%nSpan+1, p%nWings, 'Section torsion rate'   , ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName)
   call AllocAry( u%Vwnd_LLMP,3, p%nSpan+1, p%nWings, 'Dist. wind at LL nodes',  ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName)
   y%Vind    = -9999.9_ReKi  
   u%V_wind  = -9999.9_ReKi
   u%Vwnd_LLMP = -9999.9_ReKi
   u%omega_z = -9999.9_ReKi
end subroutine FVW_Init_U_Y
! ==============================================================================
!> Setting parameters *and misc* from module inputs
SUBROUTINE FVW_SetParametersFromInputs( InitInp, p, ErrStat, ErrMsg )
   type(FVW_InitInputType),    intent(inout)  :: InitInp       !< Input data for initialization routine  (inout so we can use MOVE_ALLOC)
   type(FVW_ParameterType),    intent(inout) :: p             !< Parameters
   integer(IntKi),             intent(  out) :: ErrStat       !< Error status of the operation
   character(*),               intent(  out) :: ErrMsg        !< Error message if ErrStat /= ErrID_None
   ! Local variables
   character(1024)          :: rootDir, baseName  ! Simulation root dir and basename
   integer(IntKi)          :: ErrStat2
   character(ErrMsgLen)    :: ErrMsg2
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
   call AllocAry(p%s_LL          ,  p%nSpan+1  , p%nWings, 'Spanwise coord LL  ', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); p%s_LL= -999999_ReKi;
   call AllocAry(p%s_CP_LL       ,  p%nSpan    , p%nWings, 'Spanwise coord CPll', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); p%s_CP_LL= -999999_ReKi;
   call AllocAry(p%chord_LL      ,  p%nSpan+1  , p%nWings, 'Chord on LL        ', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); p%chord_LL= -999999_ReKi;
   call AllocAry(p%chord_CP_LL   ,  p%nSpan    , p%nWings, 'Chord on CP LL     ', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); p%chord_CP_LL= -999999_ReKi;

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
      call AllocAry(p%PrescribedCirculation,  p%nSpan, 'Prescribed Circulation', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_SetParameters'); p%PrescribedCirculation = -999999_ReKi;
      if (.not. allocated(p%s_CP_LL)) then
         ErrMsg  = 'Spanwise coordinate not allocated.'
         ErrStat = ErrID_Fatal
         return
      endif
      call ReadAndInterpGamma(trim(InputFileData%CirculationFile), p%s_CP_LL(1:p%nSpan,1), p%s_LL(p%nSpan+1,1), p%PrescribedCirculation, ErrStat2, ErrMsg2)
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

   type(FVW_InputType),allocatable, intent(inout)  :: u(:)        !< System inputs
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
   real(DbKi) :: t

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""
   ! Place any last minute operations or calculations here:
   if (p%WrVTK==2) then
      call WrScr('Outputs of VTK before FVW_END')
      t=-1.0_ReKi
      m%VTKStep=999999999 ! not pretty, but we know we have twidth=9
      call WriteVTKOutputs(t, .true., u(1), p, x, z, y, m, ErrStat, ErrMsg)
   endif

   ! Close files here:
   
   ! Destroy the input data:
   if (allocated(u)) then
      do i=1,size(u)
         call FVW_DestroyInput( u(i), ErrStat, ErrMsg )
      enddo
   endif

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
   logical :: bOverCycling
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! --- Handling of time step, and time compared to previous call
   m%iStep = n
   ! OverCycling DTfvw> DTaero
   bOverCycling = p%DTfvw > p%DTaero
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
   if (DEV_VERSION)  print'(A,F10.3,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,F7.2,A,L1)','FVW status - t:',t,'  n:',n,'  nNW:',m%nNW-1,'/',p%nNWMax-1,'  nFW:',nFWEff, '+',m%nFW-nFWEff,'=',m%nFW,'/',p%nFWMax,'  nP:',nP,'  spent:', m%tSpent, 's Comp:',m%ComputeWakeInduced

   ! --- Evaluation at t
   ! Inputs at t
   call FVW_CopyInput( u(2), uInterp, MESH_NEWCOPY, ErrStat2, ErrMsg2); if(Failed()) return
   call FVW_Input_ExtrapInterp(u(1:size(utimes)),utimes(:),uInterp,t, ErrStat2, ErrMsg2); if(Failed()) return
   call Wings_Panelling(uInterp%WingsMesh, p, m, ErrStat2, ErrMsg2); if(Failed()) return
   call Map_LL_NW(p, m, z, x, 1.0_ReKi, ErrStat2, ErrMsg2); if(Failed()) return ! needed at t=0 if wing moved after init
   call Map_NW_FW(p, m, z, x, ErrStat2, ErrMsg2); if(Failed()) return

   !  TODO convert quasi steady Gamma to unsteady gamma with UA states

   ! Compute UA inputs at t
   if (m%UA_Flag) then
      call CalculateInputsAndOtherStatesForUA(1, uInterp, p, x, xd, z, OtherState, AFInfo, m, ErrStat2, ErrMsg2); if(Failed()) return
   end if

   ! --- Integration between t and t+DTfvw
   if (m%ComputeWakeInduced) then

      ! TODO TODO: this should be in CCSD, but memory is changing between time steps, so for now we have to use u(1)..
      CALL DistributeRequestedWind_NWFW(u(1)%V_wind, p, m%Vwnd_NW, m%Vwnd_FW)

      if (bOverCycling) then
         ! Store states at t, and use this opportunity to store outputs at t
         call FVW_CopyContState(x, m%x1, 0, ErrStat2, ErrMsg2) ! Backup current state at t
         m%t1=t
      endif
      if (p%IntMethod .eq. idEuler1) then 
        call FVW_Euler1( t, uInterp, p, x, xd, z, OtherState, m, ErrStat2, ErrMsg2); if(Failed()) return
      elseif (p%IntMethod .eq. idRK4) then 
         call FVW_RK4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat2, ErrMsg2); if(Failed()) return
      !elseif (p%IntMethod .eq. idAB4) then
      !   call FVW_AB4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
      !elseif (p%IntMethod .eq. idABM4) then
      !   call FVW_ABM4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
      else  
         call SetErrStat(ErrID_Fatal,'Invalid time integration method:'//Num2LStr(p%IntMethod),ErrStat,ErrMsg,'FVW_UpdateState') 
      end if
      ! We extend the wake length, i.e. we emit a new panel of vorticity at the TE
      ! NOTE: this will be rolled back if UpdateState is called at the same starting time again
      call PrepareNextTimeStep()
      ! --- t+DTfvw
      ! Propagation/creation of new layer of panels
      call PropagateWake(p, m, z, x, ErrStat2, ErrMsg2); if(Failed()) return

      if (bOverCycling) then
         ! States x1 
         ! - we need to propagate the states at t to match the memory of state t+DTfvw
         ! - the positions and intensities for the LL and 1st NW panels are NaN for x1 and x2,
         !   so we need to remap them
         call PropagateWake(p, m, z, m%x1, ErrStat2, ErrMsg2); if(Failed()) return
         !call Map_LL_NW(p, m, z, m%x1, ShedScale, ErrStat2, ErrMsg2); if(Failed()) return
         !call Map_NW_FW(p, m, z, m%x1, ErrStat2, ErrMsg2); if(Failed()) return

         ! States x2
         call FVW_CopyContState(x, m%x2, 0, ErrStat2, ErrMsg2) ! Backup current state at t+DTfvw
         m%t2=t+p%DTfvw
         !! Inputs at t+DTfvw (Wings Panelling updates CP_LL, and Vstr_LL) 
         !call FVW_Input_ExtrapInterp(u(1:size(utimes)),utimes,uInterp,t+p%DTfvw, ErrStat2, ErrMsg2); if(Failed()) return
         !call Wings_Panelling(uInterp%WingsMesh, p, m, ErrStat2, ErrMsg2); if(Failed()) return
         !! Updating positions of first NW and FW panels (Circulation also updated but irrelevant)
         !call Map_LL_NW(p, m, z, m%x2, 1.0, ErrStat2, ErrMsg2); if(Failed()) return
         !call Map_NW_FW(p, m, z, m%x2, ErrStat2, ErrMsg2); if(Failed()) return
         !! --- Solve for quasi steady circulation at t+p%DTfvw
         !! Returns: z%Gamma_LL (at t+p%DTfvw)
         !z_guess%Gamma_LL = z%Gamma_LL ! We use as guess the circulation from the previous time step (see above)
         !call FVW_CalcConstrStateResidual(t+p%DTfvw, uInterp, p, m%x2, xd, z_guess, OtherState, m, z, AFInfo, ErrStat2, ErrMsg2, 2); if(Failed()) return
         !! Compute UA inputs at t+DTfvw and integrate UA states between t and t+dtAero
         !if (m%UA_Flag) then
         !   call CalculateInputsAndOtherStatesForUA(2, uInterp, p, m%x2, xd, z, OtherState, AFInfo, m, ErrStat2, ErrMsg2); if(Failed()) return
         !   call UA_UpdateState_Wrapper(AFInfo, t, n, (/t,t+p%DTfvw/), p, m%x2, xd, OtherState, m, ErrStat2, ErrMsg2); if(Failed()) return
         !end if
         !! Updating circulation of near wake panel (and position but irrelevant)
         !call Map_LL_NW(p, m, z, m%x2, ShedScale, ErrStat2, ErrMsg2); if(Failed()) return
         !call Map_NW_FW(p, m, z, m%x2, ErrStat2, ErrMsg2); if(Failed()) return
      endif
   endif

   ! --- Integration between t and t+DTaero if DTaero/=DTfvw
   if (bOverCycling) then
      ! Linear interpolation of states between t and dtaero
      call FVW_ContStates_Interp(t+p%DTaero, (/m%x1, m%x2/), (/m%t1, m%t2/), p, m, x, ErrStat2, ErrMsg2); if(Failed()) return
   endif

   ! Inputs at t+DTaero (Wings Panelling updates CP_LL, and Vstr_LL) 
   call FVW_Input_ExtrapInterp(u(1:size(utimes)),utimes,uInterp,t+p%DTaero, ErrStat2, ErrMsg2); if(Failed()) return
   call Wings_Panelling(uInterp%WingsMesh, p, m, ErrStat2, ErrMsg2); if(Failed()) return

   ! Updating positions of first NW and FW panels (Circulation also updated but irrelevant)
   ! Changes: x only
   ShedScale = (t+p%DTaero - m%OldWakeTime)/p%DTfvw
   call Map_LL_NW(p, m, z, x, ShedScale, ErrStat2, ErrMsg2); if(Failed()) return
   call Map_NW_FW(p, m, z, x, ErrStat2, ErrMsg2); if(Failed()) return
   !call print_x_NW_FW(p, m, x,'Map2')

   ! --- Solve for quasi steady circulation at t+p%DTaero
   ! Returns: z%Gamma_LL (at t+p%DTaero)
   z_guess%Gamma_LL = z%Gamma_LL ! We use as guess the circulation from the previous time step (see above)
   call FVW_CalcConstrStateResidual(t+p%DTaero, uInterp, p, x, xd, z_guess, OtherState, m, z, AFInfo, ErrStat2, ErrMsg2, 2); if(Failed()) return
   ! Compute UA inputs at t+DTaero and integrate UA states between t and t+dtAero
   if (m%UA_Flag) then
      call CalculateInputsAndOtherStatesForUA(2, uInterp, p, x, xd, z, OtherState, AFInfo, m, ErrStat2, ErrMsg2); if(Failed()) return
      call UA_UpdateState_Wrapper(AFInfo, t, n, (/t,t+p%DTaero/), p, x, xd, OtherState, m, ErrStat2, ErrMsg2); if(Failed()) return
   end if

   ! TODO compute unsteady Gamma here based on UA Cl

   ! Updating circulation of near wake panel (and position but irrelevant)
   ! Changes: x only
   call Map_LL_NW(p, m, z, x, ShedScale, ErrStat2, ErrMsg2); if(Failed()) return
   call Map_NW_FW(p, m, z, x, ErrStat2, ErrMsg2); if(Failed()) return
   !call print_x_NW_FW(p, m, x,'Map3')

   ! --- Fake handling of ground effect (ensure vorticies above ground)
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

   if (DEV_VERSION) then
      if(have_nan(p, m, x, u, 'End Update ')) then
         STOP
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
!> This is a tight coupling routine for computing derivatives of continuous states. (CCSD)
subroutine FVW_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )
!..................................................................................................................................
   real(DbKi),                    intent(in   ) :: t          !< Current simulation time in seconds
   type(FVW_InputType),           intent(in   ) :: u          !< Inputs at t
   type(FVW_ParameterType),       intent(in   ) :: p          !< Parameters
   type(FVW_ContinuousStateType), intent(in   ) :: x          !< Continuous states at t
   type(FVW_DiscreteStateType),   intent(in   ) :: xd         !< Discrete states at t
   type(FVW_ConstraintStateType), intent(in   ) :: z          !< Constraint states at t
   type(FVW_OtherStateType),      intent(in   ) :: OtherState !< Other states at t
   type(FVW_MiscVarType),         intent(inout) :: m          !< Misc variables for optimization (not copied in glue code)
   type(FVW_ContinuousStateType), intent(inout) :: dxdt       !< Continuous state derivatives at t
   integer(IntKi),                intent(  out) :: ErrStat    !< Error status of the operation
   character(*),                  intent(  out) :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   ! Local variables
   integer(IntKi)       :: ErrStat2       ! temporary error status of the operation
   character(ErrMsgLen) :: ErrMsg2        ! temporary error message
   integer(IntKi)       :: nFWEff ! Number of farwake panels that are free at current time step
   integer(IntKi)       :: i,j,k
   real(ReKi)           :: visc_fact, age  ! Viscosity factor for diffusion of reg param
   real(ReKi), dimension(3) :: VmeanFW, VmeanNW ! Mean velocity of the near wake and far wake

   ErrStat = ErrID_None
   ErrMsg  = ""

   if (.not.allocated(dxdt%r_NW)) then
      call AllocAry( dxdt%r_NW , 3   ,  p%nSpan+1  ,p%nNWMax+1,  p%nWings, 'Wind on NW ', ErrStat2, ErrMsg2); dxdt%r_NW= -999999_ReKi;
      call AllocAry( dxdt%r_FW , 3   ,  FWnSpan+1  ,p%nFWMax+1,  p%nWings, 'Wind on FW ', ErrStat2, ErrMsg2); dxdt%r_FW= -999999_ReKi;
      if(Failed()) return
   endif

   ! Distribute the Wind we requested to Inflow wind to storage Misc arrays
   ! TODO ANDY: replace with direct call to inflow wind at r_NW and r_FW locations
   ! NOTE: this has been commented out due to some information missing at some times (and memoery reindexing)
   !       Call to inflow wind sould be done here at actual positions.
   !CALL DistributeRequestedWind_NWFW(u%V_wind, p, m%Vwnd_NW, m%Vwnd_FW)

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
         ! NOTE: mostly needed for sub-cycling and when no FW
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
   dxdt%r_NW(1:3, :, 1:iNWStart-1, :)=0.0_ReKi
   ! First FW point always convects (even if bound to NW)
   ! This is done for overcycling
   !dxdt%r_FW(1:3, :, 1, :)=0

   ! --- Regularization
   if (.not.allocated(dxdt%Eps_NW)) then
      call AllocAry( dxdt%Eps_NW , 3   ,  p%nSpan    ,p%nNWMax  ,  p%nWings, 'Eps NW ', ErrStat2, ErrMsg2); 
      call AllocAry( dxdt%Eps_FW , 3   ,  FWnSpan    ,p%nFWMax  ,  p%nWings, 'Eps FW ', ErrStat2, ErrMsg2); 
      if(Failed()) return
   endif
   if (p%WakeRegMethod==idRegConstant) then
      dxdt%Eps_NW(1:3, :, :, :)=0.0_ReKi
      dxdt%Eps_FW(1:3, :, :, :)=0.0_ReKi

   else if (p%WakeRegMethod==idRegStretching) then
      ! TODO
   else if (p%WakeRegMethod==idRegAge) then
      visc_fact = 2.0_ReKi * CoreSpreadAlpha * p%CoreSpreadEddyVisc * p%KinVisc
      ! --- Method 1, use d(rc^2)/dt = 4 k 
      dxdt%Eps_NW(1:3, :, iNWStart:, :) = visc_fact/x%Eps_NW(1:3, :, iNWStart:, :)
      dxdt%Eps_FW(1:3, :,         :, :) = visc_fact/x%Eps_FW(1:3, :, :, :)
      ! --- Method 2, use rc(tau) = 2k/sqrt(r_c^2(tau=0) + 4 k tau)
      !dxdt%Eps_NW(1:3, :, :, :) = (visc_fact)/sqrt(x%Eps_NW(1:3, :, :, :)**2 + 2*visc_fact*p%DTaero)
      !dxdt%Eps_FW(1:3, :, :, :) = (visc_fact)/sqrt(x%Eps_FW(1:3, :, :, :)**2 + 4*visc_fact*p%DTaero)
   else
      ErrStat = ErrID_Fatal
      ErrMsg ='Regularization method not implemented'
   endif
   dxdt%Eps_NW(1:3,:,1:iNWStart,:) = 0.0_ReKi ! Important! LL and First NW panel epsilon does not change

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FVW_CalcContStateDeriv') 
      Failed =  ErrStat >= AbortErrLev
   end function Failed
end subroutine FVW_CalcContStateDeriv


!> Interpolate states to the current time
!! For now: linear interpolation, two states, with t1<t2
subroutine FVW_ContStates_Interp(t, states, times, p, m, x, ErrStat, ErrMsg )
   real(DbKi),                      intent(in   )  :: t         !< Current simulation time in seconds
   type(FVW_ContinuousStateType),   intent(in   )  :: states(:) !< States at times
   real(DbKi),                      intent(in   )  :: times(:)  !< Times associated with states(:), in seconds
   type(FVW_ParameterType),         intent(in   ) :: p          !< Parameters
   type(FVW_MiscVarType),           intent(inout) :: m          !< Misc/optimization variables
   type(FVW_ContinuousStateType),   intent(inout) :: x          !< Continuous states at t on input at t + dt on output
   integer(IntKi),                  intent(  out) :: ErrStat    !< Error status of the operation
   character(*),                    intent(  out) :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   real(ReKi) :: fact
   ErrStat = ErrID_None
   ErrMsg  = "" 
   if (size(times)/=2) then
      ErrStat = ErrID_Fatal
      ErrMsg  = "FVW_ContStates_Interp: Times must be of size 2 " 
   endif
   if (times(1)>=times(2)) then
      ErrStat = ErrID_Fatal
      ErrMsg  = "FVW_ContStates_Interp: t1 must be < t2" 
   endif

   fact = (t-times(1))/(times(2)-times(1))

   x%r_NW     = (1_ReKi-fact) * states(1)%r_NW     + fact * states(2)%r_NW
   x%r_FW     = (1_ReKi-fact) * states(1)%r_FW     + fact * states(2)%r_FW
   x%Eps_NW   = (1_ReKi-fact) * states(1)%Eps_NW   + fact * states(2)%Eps_NW
   x%Eps_FW   = (1_ReKi-fact) * states(1)%Eps_FW   + fact * states(2)%Eps_FW
   x%Gamma_NW = (1_ReKi-fact) * states(1)%Gamma_NW + fact * states(2)%Gamma_NW
   x%Gamma_FW = (1_ReKi-fact) * states(1)%Gamma_FW + fact * states(2)%Gamma_FW
   !print*,'fact',fact,states(1)%Gamma_NW(29,iNWStart+1,1),x%Gamma_NW(29,iNWStart+1,1),states(2)%Gamma_NW(29,iNWStart+1,1)
   !print*,'fact',fact,states(1)%r_NW(1,29,iNWStart+1,1),x%r_NW(1,29,iNWStart+1,1),states(2)%r_NW(1,29,iNWStart+1,1)

end subroutine FVW_ContStates_Interp

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
   real(ReKi)     :: dt
   integer(IntKi)       :: ErrStat2      ! temporary error status of the operation
   character(ErrMsgLen) :: ErrMsg2       ! temporary error message
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = "" 

   dt = real(p%DTfvw,ReKi)  ! NOTE: this is DTfvw
   ! Compute "right hand side"
   CALL FVW_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, m%dxdt, ErrStat2, ErrMsg2); if (Failed()) return

   ! Update of positions and reg param
   x%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) = x%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) + dt * m%dxdt%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings)
   x%Eps_NW(1:3, 1:p%nSpan  , 1:m%nNW  , 1:p%nWings) = x%Eps_NW(1:3, 1:p%nSpan  , 1:m%nNW  , 1:p%nWings) + dt * m%dxdt%Eps_NW(1:3, 1:p%nSpan  , 1:m%nNW  , 1:p%nWings)
   if ( m%nFW>0) then
      x%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings) = x%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings) + dt * m%dxdt%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings)
      x%Eps_FW(1:3, 1:FWnSpan  , 1:m%nFW  , 1:p%nWings) = x%Eps_FW(1:3, 1:FWnSpan  , 1:m%nFW  , 1:p%nWings) + dt * m%dxdt%Eps_FW(1:3, 1:FWnSpan  , 1:m%nFW  , 1:p%nWings)
   endif
   ! Update of Gamma TODO (viscous diffusion, stretching)


   if (DEV_VERSION) then
      ! Additional checks
      if (any(m%dxdt%r_NW(1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings)<-999)) then
         print*,'FVW_Euler1: Attempting to convect NW with a wrong velocity'
         STOP
      endif
      if ( m%nFW>0) then
         if (any(m%dxdt%r_FW(1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings)<-999)) then
            call print_x_NW_FW(p, m, x, 'STP')
            print*,'FVW_Euler1: Attempting to convect FW with a wrong velocity'
            STOP
         endif
      endif
      if (any(m%dxdt%Eps_NW(1:3, 1:p%nSpan, 1:m%nNW, 1:p%nWings)<-0)) then
         print*,'FVW_Euler1: Wrong Epsilon NW'
         STOP
      endif
      if ( m%nFW>0) then
         if (any(m%dxdt%Eps_FW(1:3, 1:FWnSpan, 1:m%nFW, 1:p%nWings)<-999)) then
            print*,'FVW_Euler1: Wrong Epsilon FW'
            STOP
         endif
      endif
   endif
contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FVW_Euler1') 
      Failed =  ErrStat >= AbortErrLev
   end function Failed
end subroutine FVW_Euler1

!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine implements the fourth-order Runge-Kutta Method (RK4) for
!numerically integrating ordinary differential equations:
!!
!!   Let f(t, x) = dxdt denote the time (t) derivative of the continuous states
!(x). 
!!   Define constants k1, k2, k3, and k4 as 
!!        k1 = dt * f(t        , x_t        )
!!        k2 = dt * f(t + dt/2 , x_t + k1/2 )
!!        k3 = dt * f(t + dt/2 , x_t + k2/2 ), and
!!        k4 = dt * f(t + dt   , x_t + k3   ).
!!   Then the continuous states at t = t + dt are
!!        x_(t+dt) = x_t + k1/6 + k2/3 + k3/3 + k4/6 + O(dt^5)
!!
!! For details, see:
!! Press, W. H.; Flannery, B. P.; Teukolsky, S. A.; and Vetterling, W. T.
!"Runge-Kutta Method" and "Adaptive Step Size Control for 
!!   Runge-Kutta." Sections 16.1 and 16.2 in Numerical Recipes in FORTRAN: The
!Art of Scientific Computing, 2nd ed. Cambridge, England: 
!!   Cambridge University Press, pp. 704-716, 1992.
SUBROUTINE FVW_RK4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg)
      REAL(DbKi),                   INTENT(IN   )  :: t           !< Current simulation time in seconds
      INTEGER(IntKi),               INTENT(IN   )  :: n           !< time step number
      TYPE(FVW_InputType),           INTENT(INOUT)  :: u(:)        !< Inputs at t (out only for mesh record-keeping in ExtrapInterp routine)
      REAL(DbKi),                   INTENT(IN   )  :: utimes(:)   !< times of input
      TYPE(FVW_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(FVW_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
      TYPE(FVW_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(FVW_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t (possibly a guess)
      TYPE(FVW_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states
      TYPE(FVW_MiscVarType),         INTENT(INOUT)  :: m           !< misc/optimization variables
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
      ! local variables
      real(ReKi)                                    :: dt   
      TYPE(FVW_ContinuousStateType)                 :: k1 ! RK4 constant; see above
      TYPE(FVW_ContinuousStateType)                 :: k2 ! RK4 constant; see above 
      TYPE(FVW_ContinuousStateType)                 :: k3 ! RK4 constant; see above 
      TYPE(FVW_ContinuousStateType)                 :: k4 ! RK4 constant; see above 
      TYPE(FVW_ContinuousStateType)                 :: x_tmp       ! Holds temporary modification to x
      TYPE(FVW_InputType)                           :: u_interp    ! interpolated value of inputs 
      INTEGER(IntKi)                               :: ErrStat2    ! local error status
      CHARACTER(ErrMsgLen)                         :: ErrMsg2     ! local error message (ErrMsg)
      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""

      dt = real(p%DTfvw,ReKi) ! NOTE: this is DTfvw

      CALL FVW_CopyContState( x, k1, MESH_NEWCOPY, ErrStat2, ErrMsg2 ); CALL CheckError(ErrStat2,ErrMsg2)
      CALL FVW_CopyContState( x, k2, MESH_NEWCOPY, ErrStat2, ErrMsg2 ); CALL CheckError(ErrStat2,ErrMsg2)
      CALL FVW_CopyContState( x, k3, MESH_NEWCOPY, ErrStat2, ErrMsg2 ); CALL CheckError(ErrStat2,ErrMsg2)
      CALL FVW_CopyContState( x, k4,    MESH_NEWCOPY, ErrStat2, ErrMsg2 ); CALL CheckError(ErrStat2,ErrMsg2)
      CALL FVW_CopyContState( x, x_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2 ); CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      CALL FVW_CopyInput( u(1), u_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2 ); CALL CheckError(ErrStat2,ErrMsg2); IF ( ErrStat >= AbortErrLev ) RETURN

      ! interpolate u to find u_interp = u(t)
      CALL FVW_Input_ExtrapInterp( u(1:size(utimes)),utimes(:),u_interp, t, ErrStat2, ErrMsg2 ); CALL CheckError(ErrStat2,ErrMsg2); IF ( ErrStat >= AbortErrLev ) RETURN        

      ! find dxdt at t
      CALL FVW_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, m, m%dxdt, ErrStat2, ErrMsg2 ); CALL CheckError(ErrStat2,ErrMsg2); IF ( ErrStat >= AbortErrLev ) RETURN

      if (DEV_VERSION) then
         ! Additional checks
         if (any(m%dxdt%r_NW(1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings)<-999)) then
            print*,'FVW_RK4: Attempting to convect NW with a wrong velocity'
            STOP
         endif
         if ( m%nFW>0) then
            if (any(m%dxdt%r_FW(1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings)<-999)) then
               call print_x_NW_FW(p, m, x, 'STP')
               print*,'FVW_RK4: Attempting to convect FW with a wrong velocity'
               STOP
            endif
         endif
      endif

      k1%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) = dt * m%dxdt%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) 
      k1%Eps_NW(1:3, 1:p%nSpan  , 1:m%nNW,   1:p%nWings) = dt * m%dxdt%Eps_NW(1:3, 1:p%nSpan,   1:m%nNW,   1:p%nWings)
      if ( m%nFW>0) then   
         k1%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings) = dt * m%dxdt%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings)
         k1%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW,   1:p%nWings) = dt * m%dxdt%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW,   1:p%nWings)
      endif

      x_tmp%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) = x%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) + 0.5 * k1%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings)
      x_tmp%Eps_NW(1:3, 1:p%nSpan  , 1:m%nNW,   1:p%nWings) = x%Eps_NW(1:3, 1:p%nSpan  , 1:m%nNW,   1:p%nWings) + 0.5 * k1%Eps_NW(1:3, 1:p%nSpan,   1:m%nNW,   1:p%nWings)
      if ( m%nFW>0) then
         x_tmp%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings) = x%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings)  + 0.5 * k1%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings)
         x_tmp%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW,   1:p%nWings) = x%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW,   1:p%nWings)  + 0.5 * k1%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW,   1:p%nWings)
      endif

      ! interpolate u to find u_interp = u(t + dt/2)
      CALL FVW_Input_ExtrapInterp(u(1:size(utimes)),utimes(:),u_interp, t+0.5*dt, ErrStat2, ErrMsg2); CALL CheckError(ErrStat2,ErrMsg2); IF ( ErrStat >= AbortErrLev ) RETURN        

      ! find dxdt at t + dt/2
      CALL FVW_CalcContStateDeriv( t + 0.5*dt, u_interp, p, x_tmp, xd, z, OtherState, m, m%dxdt, ErrStat2, ErrMsg2 ); CALL CheckError(ErrStat2,ErrMsg2); IF ( ErrStat >= AbortErrLev ) RETURN

      k2%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) = dt * m%dxdt%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings)
      k2%Eps_NW(1:3, 1:p%nSpan  , 1:m%nNW,   1:p%nWings) = dt * m%dxdt%Eps_NW(1:3, 1:p%nSpan,   1:m%nNW,   1:p%nWings)
      if ( m%nFW>0) then
         k2%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings) = dt * m%dxdt%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings)
         k2%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW,   1:p%nWings) = dt * m%dxdt%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW,   1:p%nWings)
      endif

      x_tmp%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) = x%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) + 0.5 * k2%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings)
      x_tmp%Eps_NW(1:3, 1:p%nSpan  , 1:m%nNW,   1:p%nWings) = x%Eps_NW(1:3, 1:p%nSpan  , 1:m%nNW,   1:p%nWings) + 0.5 * k2%Eps_NW(1:3, 1:p%nSpan,   1:m%nNW,   1:p%nWings)
      if ( m%nFW>0) then
         x_tmp%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings) = x%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings)  + 0.5 * k2%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings)
         x_tmp%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW,   1:p%nWings) = x%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW,   1:p%nWings)  + 0.5 * k2%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW,   1:p%nWings)
      endif

      ! find dxdt at t + dt/2       
      CALL FVW_CalcContStateDeriv( t + 0.5*dt, u_interp, p, x_tmp, xd, z, OtherState, m, m%dxdt, ErrStat2, ErrMsg2 ); CALL CheckError(ErrStat2,ErrMsg2); IF ( ErrStat >= AbortErrLev ) RETURN

      k3%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) = dt * m%dxdt%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings)
      k3%Eps_NW(1:3, 1:p%nSpan  , 1:m%nNW,   1:p%nWings) = dt * m%dxdt%Eps_NW(1:3, 1:p%nSpan,   1:m%nNW,   1:p%nWings)
      if ( m%nFW>0) then
         k3%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings) = dt * m%dxdt%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings)
         k3%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW,   1:p%nWings) = dt * m%dxdt%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW,   1:p%nWings)
      endif

      x_tmp%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) = x%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) + k3%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings)
      x_tmp%Eps_NW(1:3, 1:p%nSpan  , 1:m%nNW,   1:p%nWings) = x%Eps_NW(1:3, 1:p%nSpan  , 1:m%nNW,   1:p%nWings) + k3%Eps_NW(1:3, 1:p%nSpan,   1:m%nNW,   1:p%nWings)
      if ( m%nFW>0) then
         x_tmp%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings) = x%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings)  + k3%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings)
         x_tmp%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW,   1:p%nWings) = x%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW,   1:p%nWings)  + k3%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW,   1:p%nWings)
      endif

      ! interpolate u to find u_interp = u(t + dt)
      CALL FVW_Input_ExtrapInterp(u(1:size(utimes)),utimes(:),u_interp, t + dt, ErrStat2, ErrMsg2); CALL CheckError(ErrStat2,ErrMsg2); IF ( ErrStat >= AbortErrLev ) RETURN

      ! find dxdt at t + dt
      CALL FVW_CalcContStateDeriv( t + dt, u_interp, p, x_tmp, xd, z, OtherState, m, m%dxdt, ErrStat2, ErrMsg2 ); CALL CheckError(ErrStat2,ErrMsg2); IF ( ErrStat >= AbortErrLev ) RETURN

      k4%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) = dt * m%dxdt%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings)
      k4%Eps_NW(1:3, 1:p%nSpan  , 1:m%nNW,   1:p%nWings) = dt * m%dxdt%Eps_NW(1:3, 1:p%nSpan,   1:m%nNW,   1:p%nWings)
      if ( m%nFW>0) then
         k4%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings) = dt * m%dxdt%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings)
         k4%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW,   1:p%nWings) = dt * m%dxdt%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW,   1:p%nWings)
      endif

      ! Compute and store combined dx = (k1/6 + k2/3 + k3/3 + k4/6) ! NOTE: this has dt, it's not a true dxdt yet
      m%dxdt%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) = ( k1%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) + 2._ReKi * k2%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) + 2._ReKi * k3%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) + k4%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings)  ) / 6._ReKi
      m%dxdt%Eps_NW(1:3, 1:p%nSpan  , 1:m%nNW  , 1:p%nWings) = ( k1%Eps_NW(1:3, 1:p%nSpan  , 1:m%nNW  , 1:p%nWings) + 2._ReKi * k2%Eps_NW(1:3, 1:p%nSpan  , 1:m%nNW  , 1:p%nWings) + 2._ReKi * k3%Eps_NW(1:3, 1:p%nSpan  , 1:m%nNW  , 1:p%nWings) + k4%Eps_NW(1:3, 1:p%nSpan  , 1:m%nNW  , 1:p%nWings)  ) / 6._ReKi
      if ( m%nFW>0) then         
         m%dxdt%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings) = ( k1%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings) + 2._ReKi * k2%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings) + 2._ReKi * k3%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings) + k4%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings)  ) / 6._ReKi
         m%dxdt%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW,   1:p%nWings) = ( k1%Eps_FW(1:3, 1:FWnSpan  , 1:m%nFW  , 1:p%nWings) + 2._ReKi * k2%Eps_FW(1:3, 1:FWnSpan  , 1:m%nFW  , 1:p%nWings) + 2._ReKi * k3%Eps_FW(1:3, 1:FWnSpan  , 1:m%nFW  , 1:p%nWings) + k4%Eps_FW(1:3, 1:FWnSpan  , 1:m%nFW  , 1:p%nWings)  ) / 6._ReKi
      endif

      !update positions
      x%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) = x%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) + m%dxdt%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) 
      x%Eps_NW(1:3, 1:p%nSpan  , 1:m%nNW  , 1:p%nWings) = x%Eps_NW(1:3, 1:p%nSpan  , 1:m%nNW  , 1:p%nWings) + m%dxdt%Eps_NW(1:3, 1:p%nSpan  , 1:m%nNW  , 1:p%nWings) 
      if ( m%nFW>0) then         
         x%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings) = x%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings) + m%dxdt%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings)
         x%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW,   1:p%nWings) = x%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW,   1:p%nWings) + m%dxdt%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW,   1:p%nWings)
      endif

      ! Store true dxdt =  (k1/6 + k2/3 + k3/3 + k4/6)/dt 
      m%dxdt%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings) = m%dxdt%r_NW  (1:3, 1:p%nSpan+1, 1:m%nNW+1, 1:p%nWings)/dt
      m%dxdt%Eps_NW(1:3, 1:p%nSpan  , 1:m%nNW  , 1:p%nWings) = m%dxdt%Eps_NW(1:3, 1:p%nSpan  , 1:m%nNW  , 1:p%nWings)/dt
      if ( m%nFW>0) then         
         m%dxdt%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings) = m%dxdt%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1, 1:p%nWings)/dt
         m%dxdt%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW,   1:p%nWings) = m%dxdt%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW,   1:p%nWings)/dt
      endif

      ! clean up local variables:
      CALL ExitThisRoutine(  )

CONTAINS
   !> This subroutine destroys all the local variables
   SUBROUTINE ExitThisRoutine()
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)
      CALL FVW_DestroyContState( k1,       ErrStat3, ErrMsg3 )
      CALL FVW_DestroyContState( k2,       ErrStat3, ErrMsg3 )
      CALL FVW_DestroyContState( k3,       ErrStat3, ErrMsg3 )
      CALL FVW_DestroyContState( k4,       ErrStat3, ErrMsg3 )
      CALL FVW_DestroyContState( x_tmp,    ErrStat3, ErrMsg3 )
      CALL FVW_DestroyInput(     u_interp, ErrStat3, ErrMsg3 )
   END SUBROUTINE ExitThisRoutine
   !> This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   SUBROUTINE CheckError(ErrID,Msg)
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)
      IF ( ErrID /= ErrID_None ) THEN
         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'FVW_RK4:'//TRIM(Msg)
         ErrStat = MAX(ErrStat,ErrID)
         IF ( ErrStat >= AbortErrLev ) CALL ExitThisRoutine( )
      END IF
   END SUBROUTINE CheckError
END SUBROUTINE FVW_RK4


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

   ! Distribute the Wind we requested to Inflow wind to storage Misc arrays
   ! TODO ANDY: replace with direct call to inflow wind at m%CP_LL location
   CALL DistributeRequestedWind_LL(u%V_wind, p, m%Vwnd_LL)

   ! Solve for the residual of the constraint state functions here:
   !z%residual = 0.0_ReKi
   !z%Gamma_LL = 0.0_ReKi
   call AllocAry( z_out%Gamma_LL,  p%nSpan, p%nWings, 'Lifting line Circulation', ErrStat, ErrMsg );
   z_out%Gamma_LL = -999999_ReKi;

   CALL Wings_ComputeCirculation(t, z_out%Gamma_LL, z_guess%Gamma_LL, u, p, x, m, AFInfo, ErrStat, ErrMsg, iLabel)

end subroutine FVW_CalcConstrStateResidual


subroutine CalcOutputForAD(t, u, p, x, y, m, AFInfo, ErrStat, ErrMsg)
   real(DbKi),                      intent(in   )  :: t           !< Current simulation time in seconds
   type(FVW_InputType),             intent(in   )  :: u           !< Inputs at Time t
   type(FVW_ParameterType),         intent(in   )  :: p           !< Parameters
   type(FVW_ContinuousStateType),   intent(in   )  :: x           !< Continuous states at t
   type(FVW_OutputType),            intent(inout)  :: y           !< Outputs computed at t (Input only so that mesh con-
   type(FVW_MiscVarType),           intent(inout)  :: m           !< Misc/optimization variables
   type(AFI_ParameterType),         intent(in   )  :: AFInfo(:)   !< The airfoil parameter data
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   integer(IntKi) :: iW, n
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   character(*), parameter       :: RoutineName = 'FVW_CalcOutput'

   ErrStat = ErrID_None
   ErrMsg  = ""
!       ! --- NOTE: this below might not be needed
!       ! Distribute the Wind we requested to Inflow wind to storage Misc arrays
!       ! TODO ANDY: replace with direct call to inflow wind at m%CP_LL location
!       CALL DistributeRequestedWind_LL(u%V_wind, p, m%Vwnd_LL)
! 
!       ! Control points location and structrual velocity
   call Wings_Panelling(u%WingsMesh, p, m, ErrStat2, ErrMsg2);
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
! 
!    ! if we are on a correction step, CalcOutput may be called again with different inputs
!    ! Compute m%Gamma_LL
!    CALL Wings_ComputeCirculation(t, m%Gamma_LL, z%Gamma_LL, u, p, x, m, AFInfo, ErrStat2, ErrMsg2, 0); if(Failed()) return ! For plotting only
   !---

   ! Induction on the lifting line control point
   ! Compute m%Vind_LL
   m%Vind_LL=-9999.0_ReKi
   call LiftingLineInducedVelocities(m%CP_LL, p, x, 1, m, m%Vind_LL, ErrStat2, ErrMsg2); 
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ! Induction on the mesh points (AeroDyn nodes)
   n=p%nSpan
   y%Vind(1:3,:,:) = 0.0_ReKi
   do iW=1,p%nWings
      ! --- Linear interpolation for interior points and extrapolations at boundaries
      call interpextrap_cp2node(p%s_CP_LL(:,iW), m%Vind_LL(1,:,iW), p%s_LL(:,iW), y%Vind(1,:,iW))
      call interpextrap_cp2node(p%s_CP_LL(:,iW), m%Vind_LL(2,:,iW), p%s_LL(:,iW), y%Vind(2,:,iW))
      call interpextrap_cp2node(p%s_CP_LL(:,iW), m%Vind_LL(3,:,iW), p%s_LL(:,iW), y%Vind(3,:,iW))
   enddo
end subroutine CalcOutputForAD
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
subroutine FVW_CalcOutput(t, u, p, x, xd, z, OtherState, AFInfo, y, m, ErrStat, ErrMsg)
   use FVW_VTK, only: set_vtk_coordinate_transform
   use FVW_VortexTools, only: interpextrap_cp2node
   real(DbKi),                      intent(in   )  :: t           !< Current simulation time in seconds
   type(FVW_InputType),             intent(in   )  :: u           !< Inputs at Time t
   type(FVW_ParameterType),         intent(in   )  :: p           !< Parameters
   type(FVW_ContinuousStateType),   intent(in   )  :: x           !< Continuous states at t
   type(FVW_DiscreteStateType),     intent(in   )  :: xd          !< Discrete states at t
   type(FVW_ConstraintStateType),   intent(in   )  :: z           !< Constraint states at t
   type(FVW_OtherStateType),        intent(in   )  :: OtherState  !< Other states at t
   type(AFI_ParameterType),         intent(in   )  :: AFInfo(:)   !< The airfoil parameter data
   type(FVW_OutputType),            intent(inout)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                                  !!   nectivity information does not have to be recalculated)
   type(FVW_MiscVarType),           intent(inout)  :: m           !< Misc/optimization variables
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   ! Local variables
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   character(*), parameter       :: RoutineName = 'FVW_CalcOutput'
   logical :: bOverCycling
   real(ReKi) :: fact
   ErrStat = ErrID_None
   ErrMsg  = ""
   if (DEV_VERSION) then
      print'(A,F10.3,A,L1,A,I0,A,I0)','CalcOutput     t:',t,'   ',m%FirstCall,'                                nNW:',m%nNW,' nFW:',m%nFW
   endif

   ! OverCycling DTfvw> DTaero
   bOverCycling = p%DTfvw > p%DTaero

   ! Compute induced velocity at AD nodes
   call CalcOutputForAD(t,u,p,x,y,m,AFInfo, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   
   ! Export to VTK
   if (m%VTKStep==-1) then 
      m%VTKStep = 0 ! Has never been called, special handling for init
   else
      m%VTKStep = m%iStep+1 ! We use glue code step number for outputs
   endif
   call WriteVTKOutputs(t, .False., u, p, x, z, y, m, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

end subroutine FVW_CalcOutput

!> Write to  vtk_fvw folder at fps requested
subroutine WriteVTKOutputs(t, force, u, p, x, z, y, m, ErrStat, ErrMsg)
   real(DbKi),                      intent(in   )  :: t       !< Current simulation time in seconds
   logical,                         intent(in   )  :: force   !< force the writing
   type(FVW_InputType),             intent(in   )  :: u       !< Inputs at Time t
   type(FVW_ParameterType),         intent(in   )  :: p       !< Parameters
   type(FVW_ContinuousStateType),   intent(in   )  :: x       !< Continuous states at t
   type(FVW_ConstraintStateType),   intent(in   )  :: z       !< Constraint states at t
   type(FVW_OutputType),            intent(in   )  :: y       !< Outputs computed at t (Input only so that mesh con-
   type(FVW_MiscVarType),           intent(inout)  :: m       !< Misc/optimization variables
   integer(IntKi),                  intent(  out)  :: ErrStat !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg  !< Error message if ErrStat /= ErrID_None
   ! Local variables
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   character(*), parameter       :: RoutineName = 'FVW_CalcOutput'
   integer(IntKi) :: iW, iGrid
   integer(IntKi) :: nSeg, nSegP
   if (p%WrVTK>0) then
      if (m%FirstCall .or. force) then
         call MKDIR(p%VTK_OutFileRoot)
      endif
      ! For plotting only
      call PackPanelsToSegments(p, x, 1, (p%ShearModel==idShearMirror), m%nNW, m%nFW, m%Sgmt%Connct, m%Sgmt%Points, m%Sgmt%Gamma, m%Sgmt%Epsilon, nSeg, nSegP)
      m%Vtot_LL = m%Vind_LL + m%Vwnd_LL - m%Vstr_LL
      if (DEV_VERSION) then
         call print_mean_3d(m%Vind_LL,'Mean induced vel. LL')
         call print_mean_3d(m%Vtot_LL,'Mean relativevel. LL')
      endif
      if ( force .or. (( t - m%VTKlastTime ) >= p%DTvtk*OneMinusEpsilon ))  then
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
   ! --- Write VTK grids
   if (p%nGridOut>0) then
      if (m%FirstCall .or. force) then
         call MKDIR(p%VTK_OutFileRoot)
      endif
      ! Distribute the Wind we requested to Inflow wind to storage Misc arrays
      ! TODO ANDY: replace with direct call to inflow wind at Grid points
      CALL DistributeRequestedWind_Grid(u%V_wind, p, m)
      do iGrid=1,p%nGridOut
         if (force.or. (( t - m%GridOutputs(iGrid)%tLastOutput) >= m%GridOutputs(iGrid)%DTout * OneMinusEpsilon) )  then
            ! Compute induced velocity on grid, TODO use the same Tree for all CalcOutput
            call InducedVelocitiesAll_OnGrid(m%GridOutputs(iGrid), p, x, m, ErrStat2, ErrMsg2);
            m%GridOutputs(iGrid)%tLastOutput = t
            call WrVTK_FVW_Grid(p, x, z, m, iGrid, trim(p%VTK_OutFileBase)//'FVW_Grid', m%VTKStep, 9)
         endif
      enddo
   endif
end subroutine WriteVTKOutputs

!----------------------------------------------------------------------------------------------------------------------------------   
! --- UA related, should be merged with AeroDyn 
!----------------------------------------------------------------------------------------------------------------------------------   
!> Init UA
!! NOTE: UA is done at the "AeroDyn" nodes, not the control points
subroutine UA_Init_Wrapper(AFInfo, InitInp, interval, p, x, xd, OtherState, m, ErrStat, ErrMsg )
   use UnsteadyAero, only: UA_Init
   type(AFI_ParameterType),         intent(in   )  :: AFInfo(:)   !< The airfoil parameter data, temporary, for UA..
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
   type(UA_InitOutputType):: InitOutData_UA
   integer                :: i,j
   integer(intKi)         :: ErrStat2
   character(ErrMsgLen)   :: ErrMsg2
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   m%UA_Flag=InitInp%UA_Flag
   ! --- Condensed version of BEMT_Init_Otherstate
   if ( m%UA_Flag ) then
      ! ---Condensed version of "BEMT_Set_UA_InitData"
      allocate(Init_UA_Data%c(InitInp%numBladeNodes,InitInp%numBlades), STAT = errStat2)
      do j = 1,InitInp%numBlades
         do i = 1,InitInp%numBladeNodes
            Init_UA_Data%c(i,j)      = p%chord(i,j) ! NOTE: InitInp chord move-allocd to p
         end do
      end do
      Init_UA_Data%dt              = interval          

      Init_UA_Data%OutRootName     = 'Debug.UA'
      Init_UA_Data%numBlades       = InitInp%numBlades 
      Init_UA_Data%nNodesPerBlade  = InitInp%numBladeNodes ! At AeroDyn ndoes, not CP
      Init_UA_Data%UAMod           = InitInp%UAMod  
      Init_UA_Data%Flookup         = InitInp%Flookup
      Init_UA_Data%a_s             = InitInp%a_s ! Speed of sound, m/s  
      Init_UA_Data%ShedEffect      = .False. ! Important, when coupling UA wih vortex code, shed vorticity is inherently accounted for

      ! --- UA init
      allocate(m%u_UA( InitInp%numBladeNodes, InitInp%numBlades, 2), stat=errStat2) 
      call UA_Init( Init_UA_Data, m%u_UA(1,1,1), m%p_UA, x%UA, xd%UA, OtherState%UA, m%y_UA, m%m_UA, interval, AFInfo, p%AFIndx, InitOutData_UA, ErrStat2, ErrMsg2); if(Failed())return

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

!> Compute necessary inputs for UA at a given time step, stored in m%u_UA
!!  Inputs are AoA, U, Re, 
!!  See equivalent version in BEMT, and SetInputs_for_UA in BEMT
subroutine CalculateInputsAndOtherStatesForUA(InputIndex, u, p, x, xd, z, OtherState, AFInfo, m, ErrStat, ErrMsg)
   integer(IntKi),                     intent(in   ) :: InputIndex ! InputIndex= 1 or 2, depending on time step we are calculating inputs for
   type(FVW_InputType),                intent(in   ) :: u          ! Input
   type(FVW_ParameterType),            intent(in   ) :: p          ! Parameters   
   type(FVW_ContinuousStateType),      intent(in   ) :: x          ! Continuous states at given time step
   type(FVW_DiscreteStateType),        intent(in   ) :: xd         ! Discrete states at given time step
   type(FVW_ConstraintStateType),      intent(in   ) :: z          ! Constraint states at given time step
   type(FVW_OtherStateType),           intent(inout) :: OtherState ! Other states at given time step
   type(FVW_MiscVarType), target,      intent(inout) :: m          ! Misc/optimization variables
   type(AFI_ParameterType),            intent(in   ) :: AFInfo(:)  ! The airfoil parameter data
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   ! Local
   real(ReKi), dimension(:,:), allocatable :: Vind_node
   type(UA_InputType), pointer     :: u_UA ! Alias to shorten notations
   integer(IntKi)                                    :: i,j
   character(ErrMsgLen)                              :: errMsg2     ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)                                    :: errStat2    ! temporary Error status of the operation
   ErrStat  = ErrID_None
   ErrMsg   = ""

   ! --- Induction on the lifting line control points
   ! NOTE: this is expensive since it's an output for FVW but here we have to use it for UA
   ! Set m%Vind_LL
   m%Vind_LL=-9999.0_ReKi
   call LiftingLineInducedVelocities(m%CP_LL, p, x, 1, m, m%Vind_LL, ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'UA_UpdateState_Wrapper'); if (ErrStat >= AbortErrLev) return
   allocate(Vind_node(3,1:p%nSpan+1))

   do j = 1,p%nWings  
      ! Induced velocity at Nodes (NOTE: we rely on storage done when computing Circulation)
      if (m%nNW>1) then
         call interpextrap_cp2node(p%s_CP_LL(:,j), m%Vind_LL(1,:,j), p%s_LL(:,j), Vind_node(1,:))
         call interpextrap_cp2node(p%s_CP_LL(:,j), m%Vind_LL(2,:,j), p%s_LL(:,j), Vind_node(2,:))
         call interpextrap_cp2node(p%s_CP_LL(:,j), m%Vind_LL(3,:,j), p%s_LL(:,j), Vind_node(3,:))
      else
         Vind_node=0.0_ReKi
      endif
      do i = 1,p%nSpan+1 
         ! We only update the UnsteadyAero states if we have unsteady aero turned on for this node      
         u_UA => m%u_UA(i,j,InputIndex) ! Alias
         !! ....... compute inputs to UA ...........
         ! NOTE: To be consistent with CalcOutput we take Vwind_ND that was set using m%DisturbedInflow from AeroDyn.. 
         ! This is not clean, but done to be consistent, waiting for AeroDyn to handle UA
         call AlphaVrel_Generic(u%WingsMesh(j)%Orientation(1:3,1:3,i), u%WingsMesh(j)%TranslationVel(1:3,i),  Vind_node(1:3,i), u%Vwnd_LLMP(1:3,i,j), &
                                 p%KinVisc, p%Chord(i,j), u_UA%U, u_UA%alpha, u_UA%Re)
         u_UA%v_ac(1)  = sin(u_UA%alpha)*u_UA%U
         u_UA%v_ac(2)  = cos(u_UA%alpha)*u_UA%U
         u_UA%omega    = u%omega_z(i,j)
         u_UA%UserProp = 0 ! u1%UserProp(i,j) ! TODO
      end do ! i nSpan
   end do ! j nWings
   deallocate(Vind_node)

contains
   function NodeText(i,j)
      integer(IntKi), intent(in) :: i ! node number
      integer(IntKi), intent(in) :: j ! blade number
      character(25)              :: NodeText
      NodeText = '(nd:'//trim(num2lstr(i))//' bld:'//trim(num2lstr(j))//')'
   end function NodeText
end subroutine CalculateInputsAndOtherStatesForUA


subroutine UA_UpdateState_Wrapper(AFInfo, t, n, uTimes, p, x, xd, OtherState, m, ErrStat, ErrMsg )
   use FVW_VortexTools, only: interpextrap_cp2node
   use UnsteadyAero, only: UA_UpdateStates
   type(AFI_ParameterType),         intent(in   )  :: AFInfo(:)   !< The airfoil parameter data, temporary, for UA..
   real(DbKi),                      intent(in   )  :: t           !< Curent time
   real(DbKi),                      intent(in   )  :: uTimes(:)   !< Simulation times where 
   integer(IntKi),                  intent(in   )  :: n           !< time step  
   type(FVW_ParameterType),         intent(in   )  :: p           !< Parameters
   type(FVW_ContinuousStateType),   intent(inout)  :: x           !< Initial continuous states
   type(FVW_DiscreteStateType),     intent(inout)  :: xd          !< Initial discrete states
   type(FVW_OtherStateType),        intent(inout)  :: OtherState  !< Initial other states
   type(FVW_MiscVarType),           intent(inout)  :: m           !< Initial misc/optimization variables
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   ! Local
   integer                :: i,j
   integer(intKi)         :: ErrStat2           ! temporary Error status
   character(ErrMsgLen)   :: ErrMsg2
   ErrStat  = ErrID_None
   ErrStat2 = ErrID_None
   ErrMsg   = ""
   ErrMsg2  = ""

   ! --- Condensed version of BEMT_Update States
   do j = 1,p%nWings  
      do i = 1,p%nSpan+1 
         ! COMPUTE: x%UA, xd%UA
         call UA_UpdateStates( i, j, t, n, m%u_UA(i,j,:), uTimes, m%p_UA, x%UA, xd%UA, OtherState%UA, AFInfo(p%AFIndx(i,j)), m%m_UA, errStat2, errMsg2 )
         if (ErrStat2 /= ErrID_None) then
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'UA_UpdateState_Wrapper'//trim(NodeText(i,j)))
            call WrScr(trim(ErrMsg))
            if (ErrStat >= AbortErrLev) return
         end if
      end do ! i span
   end do !j wings
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
