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

   type(ProgDesc), parameter  :: FVW_Ver = ProgDesc( 'OLAF', '', '' )

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
   integer :: iW

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
   if (p%WakeAtTE) then
      p%iNWStart=2
   else
      p%iNWStart=1
   endif

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
   integer(IntKi) :: iW
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

   m%FirstCall = .True.
   m%nNW       = p%iNWStart-1  ! Number of active nearwake panels
   m%nFW       = 0             ! Number of active farwake  panels
   m%iStep     = 0             ! Current step number
   m%VTKStep   = -1            ! Counter of VTK outputs
   m%VTKlastTime = -HUGE(1.0_DbKi)
   m%OldWakeTime = -HUGE(1.0_DbKi)

   allocate(m%W(p%nWings))
   allocate(m%dxdt%W(p%nWings))
   do iW = 1,p%nWings
      m%W(iW)%iTip=-1 ! Imort init
      m%W(iW)%iRoot=-1 ! Imort init

      call AllocAry( m%W(iW)%LE           , 3,p%W(iW)%nSpan+1, 'Leading Edge Points', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%LE = -999999_ReKi;
      call AllocAry( m%W(iW)%TE           , 3,p%W(iW)%nSpan+1, 'TrailingEdge Points', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%TE = -999999_ReKi;
      call AllocAry( m%W(iW)%PitchAndTwist,   p%W(iW)%nSpan+1, 'Pitch and twist    ', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%PitchAndTwist= -999999_ReKi;
      call AllocAry( m%W(iW)%alpha_LL,        p%W(iW)%nSpan  , 'Wind on CP ll      ', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%alpha_LL= -999999_ReKi;
      call AllocAry( m%W(iW)%Vreln_LL,        p%W(iW)%nSpan  , 'Wind on CP ll      ', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%Vreln_LL = -999999_ReKi;
      ! Variables at control points/elements
      call AllocAry( m%W(iW)%CP      , 3   ,  p%W(iW)%nSpan, 'Control points LL  ', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%CP= -999999_ReKi;
      call AllocAry( m%W(iW)%Tang    , 3   ,  p%W(iW)%nSpan, 'Tangential vector  ', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%Tang= -999999_ReKi;
      call AllocAry( m%W(iW)%Norm    , 3   ,  p%W(iW)%nSpan, 'Normal     vector  ', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%Norm= -999999_ReKi;
      call AllocAry( m%W(iW)%Orth    , 3   ,  p%W(iW)%nSpan, 'Orthogonal vector  ', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%Orth= -999999_ReKi;
      call AllocAry( m%W(iW)%dl      , 3   ,  p%W(iW)%nSpan, 'Orthogonal vector  ', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%dl= -999999_ReKi;
      call AllocAry( m%W(iW)%Area    ,        p%W(iW)%nSpan, 'LL Panel area      ', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%Area = -999999_ReKi;
      call AllocAry( m%W(iW)%diag_LL ,        p%W(iW)%nSpan, 'LL Panel diagonals ', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%diag_LL = -999999_ReKi;
      call AllocAry( m%W(iW)%Vind_CP , 3   ,  p%W(iW)%nSpan, 'Vind on CP ll      ', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%Vind_CP= -999999_ReKi;
      call AllocAry( m%W(iW)%Vtot_CP , 3   ,  p%W(iW)%nSpan, 'Vtot on CP ll      ', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%Vtot_CP= -999999_ReKi;
      call AllocAry( m%W(iW)%Vstr_CP , 3   ,  p%W(iW)%nSpan, 'Vstr on CP ll      ', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%Vstr_CP= -999999_ReKi;
      call AllocAry( m%W(iW)%Vwnd_CP , 3   ,  p%W(iW)%nSpan, 'Wind on CP ll      ', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%Vwnd_CP= -999999_ReKi;
      ! Variables at panels points
      call AllocAry( m%W(iW)%r_LL    , 3   ,  p%W(iW)%nSpan+1  , 2        , 'Lifting Line Panels', ErrStat2, ErrMsg2 );call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%r_LL= -999999_ReKi;
      call AllocAry( m%W(iW)%Vind_LL , 3   ,  p%W(iW)%nSpan+1,              'Vind on CP ll      ', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%Vind_LL= -999999_ReKi;
      !call AllocAry( m%W(iW)%Vtot_LL , 3   ,  p%W(iW)%nSpan+1,              'Vtot on CP ll      ', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%Vtot_LL= -999999_ReKi;
      !call AllocAry( m%W(iW)%Vstr_LL , 3   ,  p%W(iW)%nSpan+1,              'Vstr on CP ll      ', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%Vstr_LL= -999999_ReKi;
      !call AllocAry( m%W(iW)%Vwnd_LL , 3   ,  p%W(iW)%nSpan+1,              'Wind on CP ll      ', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%Vwnd_LL= -999999_ReKi;
      call AllocAry( m%W(iW)%Vwnd_NW , 3   ,  p%W(iW)%nSpan+1  ,p%nNWMax+1, 'Wind on NW ', ErrStat2, ErrMsg2 );call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%Vwnd_NW= -999_ReKi;
      call AllocAry( m%W(iW)%Vwnd_FW , 3   ,  FWnSpan+1  ,p%nFWMax+1, 'Wind on FW ', ErrStat2, ErrMsg2 );call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%Vwnd_FW= -999_ReKi;
      call AllocAry( m%W(iW)%Vind_NW , 3   ,  p%W(iW)%nSpan+1  ,p%nNWMax+1, 'Vind on NW ', ErrStat2, ErrMsg2 );call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%Vind_NW= -999_ReKi;
      call AllocAry( m%W(iW)%Vind_FW , 3   ,  FWnSpan+1  ,p%nFWMax+1, 'Vind on FW ', ErrStat2, ErrMsg2 );call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%Vind_FW= -999_ReKi;
      ! Variables for optimizing outputs at blade nodes
      call AllocAry( m%W(iW)%BN_UrelWind_s, 3, p%W(iW)%nSpan+1 , 'Relative wind in section coordinates',   ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%BN_UrelWind_s= -999999_ReKi;
      call AllocAry( m%W(iW)%BN_AxInd   ,      p%W(iW)%nSpan+1 , 'Axial induction',                        ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%BN_AxInd     = -999999_ReKi;
      call AllocAry( m%W(iW)%BN_TanInd  ,      p%W(iW)%nSpan+1 , 'Tangential induction',                   ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%BN_TanInd    = -999999_ReKi;
      call AllocAry( m%W(iW)%BN_Vrel    ,      p%W(iW)%nSpan+1 , 'Relative velocity',                      ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%BN_Vrel      = -999999_ReKi;
      call AllocAry( m%W(iW)%BN_alpha   ,      p%W(iW)%nSpan+1 , 'Angle of attack',                        ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%BN_alpha     = -999999_ReKi;
      call AllocAry( m%W(iW)%BN_phi     ,      p%W(iW)%nSpan+1 , 'angle between the plane local wind dir', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%BN_phi       = -999999_ReKi;
      call AllocAry( m%W(iW)%BN_Re      ,      p%W(iW)%nSpan+1 , 'Reynolds number',                        ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%BN_Re        = -999999_ReKi;
      call AllocAry( m%W(iW)%BN_Cl_Static ,    p%W(iW)%nSpan+1 , 'Coefficient lift - no UA',               ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%BN_Cl_Static = -999999_ReKi;
      call AllocAry( m%W(iW)%BN_Cd_Static ,    p%W(iW)%nSpan+1 , 'Coefficient drag - no UA',               ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%BN_Cd_Static = -999999_ReKi;
      call AllocAry( m%W(iW)%BN_Cm_Static ,    p%W(iW)%nSpan+1 , 'Coefficient moment - no UA',             ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%BN_Cm_Static = -999999_ReKi;
      call AllocAry( m%W(iW)%BN_Cl        ,    p%W(iW)%nSpan+1 , 'Coefficient lift - with UA',             ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%BN_Cl        = -999999_ReKi;
      call AllocAry( m%W(iW)%BN_Cd        ,    p%W(iW)%nSpan+1 , 'Coefficient drag - with UA',             ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%BN_Cd        = -999999_ReKi;
      call AllocAry( m%W(iW)%BN_Cm        ,    p%W(iW)%nSpan+1 , 'Coefficient moment - with UA',           ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%BN_Cm        = -999999_ReKi;
      call AllocAry( m%W(iW)%BN_Cx        ,    p%W(iW)%nSpan+1 , 'Coefficient normal (to plane)',          ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%BN_Cx        = -999999_ReKi;
      call AllocAry( m%W(iW)%BN_Cy        ,    p%W(iW)%nSpan+1 , 'Coefficient tangential (to plane)',      ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName ); m%W(iW)%BN_Cy        = -999999_ReKi;
      ! dxdt, to avoid realloc all the time, and storage for subcycling 
      call AllocAry( m%dxdt%W(iW)%r_NW , 3   ,  p%W(iW)%nSpan+1 , p%nNWMax+1, 'r NW dxdt'  , ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%dxdt%W(iW)%r_NW = -999999_ReKi;
      call AllocAry( m%dxdt%W(iW)%r_FW , 3   ,  FWnSpan+1 , p%nFWMax+1, 'r FW dxdt'  , ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%dxdt%W(iW)%r_FW = -999999_ReKi;
      call AllocAry( m%dxdt%W(iW)%Eps_NW, 3  ,  p%W(iW)%nSpan    ,p%nNWMax  , 'Eps NW dxdt', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%dxdt%W(iW)%Eps_NW = -999999_ReKi;
      call AllocAry( m%dxdt%W(iW)%Eps_FW, 3  ,  FWnSpan    ,p%nFWMax  , 'Eps FW dxdt', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%dxdt%W(iW)%Eps_FW = -999999_ReKi;

      ! Wind set to 0. TODO check if -99999 works now
      !NOTE: We do not have the windspeed until after the FVW initialization (IfW is not initialized until after AD15)
      m%W(iW)%Vwnd_CP(:,:)   = 0
      m%W(iW)%Vwnd_NW(:,:,:) = 0
      m%W(iW)%Vwnd_FW(:,:,:) = 0
   enddo

   ! Wind request points
   nMax = 0
   do iW = 1,p%nWings
      nMax = nMax +  p%W(iW)%nSpan                    ! Lifting line Control Points
      nMax = nMax + (p%W(iW)%nSpan+1) * (p%nNWMax+1)  ! Nearwake points
      nMax = nMax + (FWnSpan+1) * (p%nFWMax+1)  ! Far wake points
   enddo
   ! Grid outputs
   do iGrid=1,p%nGridOut
      nMax = nMax + m%GridOutputs(iGrid)%nx * m%GridOutputs(iGrid)%ny * m%GridOutputs(iGrid)%nz
      call AllocAry(m%GridOutputs(iGrid)%uGrid, 3, m%GridOutputs(iGrid)%nx,  m%GridOutputs(iGrid)%ny, m%GridOutputs(iGrid)%nz, 'uGrid', ErrStat2, ErrMsg2);
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName)
      if (m%GridOutputs(iGrid)%type==idGridVelVorticity) then
         call AllocAry(m%GridOutputs(iGrid)%omGrid, 3, m%GridOutputs(iGrid)%nx,  m%GridOutputs(iGrid)%ny, m%GridOutputs(iGrid)%nz, 'omGrid', ErrStat2, ErrMsg2);
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName)
      endif
      m%GridOutputs(iGrid)%tLastOutput = -HUGE(1.0_DbKi)
   enddo
   call AllocAry( m%r_wind, 3, nMax, 'Requested wind points', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName )
   m%r_wind = 0.0_ReKi     ! set to zero so InflowWind can shortcut calculations

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
   integer :: iW
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

   allocate(x%W(p%nWings))
   do iW=1,p%nWings
      call AllocAry( x%W(iW)%Gamma_NW,    p%W(iW)%nSpan   , p%nNWMax  , 'NW Panels Circulation', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitStates' ); 
      call AllocAry( x%W(iW)%Gamma_FW,    FWnSpan   , p%nFWMax  , 'FW Panels Circulation', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitStates' ); 
      call AllocAry( x%W(iW)%Eps_NW  , 3, p%W(iW)%nSpan   , p%nNWMax  , 'NW Panels Reg Param'  , ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitStates' );
      call AllocAry( x%W(iW)%Eps_FW  , 3, FWnSpan   , p%nFWMax  , 'FW Panels Reg Param'  , ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitStates' );
      ! set x%W(iW)%r_NW and x%W(iW)%r_FW to (0,0,0) so that InflowWind can shortcut the calculations
      call AllocAry( x%W(iW)%r_NW    , 3, p%W(iW)%nSpan+1 , p%nNWMax+1, 'NW Panels Points'     , ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitStates' );
      call AllocAry( x%W(iW)%r_FW    , 3, FWnSpan+1 , p%nFWMax+1, 'FW Panels Points'     , ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitStates' );
      if (ErrStat >= AbortErrLev) return
      x%W(iW)%r_NW     = 0.0_ReKi
      x%W(iW)%r_FW     = 0.0_ReKi
      x%W(iW)%Gamma_NW = 0.0_ReKi ! First call of calcoutput, states might not be set 
      x%W(iW)%Gamma_FW = 0.0_ReKi ! NOTE, these values might be mapped from z%W(iW)%Gamma_LL at init
      x%W(iW)%Eps_NW   = 0.001_ReKi 
      x%W(iW)%Eps_FW   = 0.001_ReKi 
   enddo
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
   integer :: iW
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""
   !
   allocate(z%W(p%nWings))
   do iW=1,p%nWings 
      call AllocAry( z%W(iW)%Gamma_LL,  p%W(iW)%nSpan, 'Lifting line Circulation', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitConstraint' );
      !z%W(iW)%Gamma_LL = -999999_ReKi
      z%W(iW)%Gamma_LL = 0.0_ReKi
   enddo

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
   integer :: iW
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""
   ! Wind at requested points
   call AllocAry(u%V_wind, 3, size(m%r_wind,2), 'Wind Velocity at points', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName)
   u%V_wind  = -9999.9_ReKi

   allocate(y%W(p%nWings))
   allocate(u%W(p%nWings))
   do iW=1,p%nWings
      call AllocAry( y%W(iW)%Vind ,    3, p%W(iW)%nSpan+1, 'Induced velocity vector', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName)
      call AllocAry( u%W(iW)%omega_z,     p%W(iW)%nSpan+1, 'Section torsion rate'   , ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName)
      call AllocAry( u%W(iW)%Vwnd_LL,  3, p%W(iW)%nSpan+1, 'Dist. wind at LL nodes',  ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName)
      y%W(iW)%Vind    = -9999.9_ReKi  
      u%W(iW)%Vwnd_LL = -9999.9_ReKi
      u%W(iW)%omega_z = -9999.9_ReKi
   enddo
   ! Rotors, contain hub info
   allocate(u%rotors(p%nRotors))


end subroutine FVW_Init_U_Y
! ==============================================================================
!> Setting parameters *and misc* from module inputs
SUBROUTINE FVW_SetParametersFromInputs( InitInp, p, ErrStat, ErrMsg )
   type(FVW_InitInputType),    intent(inout)  :: InitInp       !< Input data for initialization routine  (inout so we can use MOVE_ALLOC)
   type(FVW_ParameterType),    intent(inout) :: p             !< Parameters
   integer(IntKi),             intent(  out) :: ErrStat       !< Error status of the operation
   character(*),               intent(  out) :: ErrMsg        !< Error message if ErrStat /= ErrID_None
   ! Local variables
   character(1024)         :: rootDir, baseName  ! Simulation root dir and basename
   integer(IntKi)          :: iW, nRotors, nBldMax
   integer(IntKi), allocatable :: nBldPerRot(:)
   integer(IntKi)          :: ErrStat2
   character(ErrMsgLen)    :: ErrMsg2
   character(*), parameter :: RoutineName = 'FVW_SetParametersFromInputs'
   ErrStat = ErrID_None
   ErrMsg  = ""
   ! 
   p%nWings       = size(InitInp%WingsMesh)
   p%DTaero       = InitInp%DTaero          ! AeroDyn Time step
   p%KinVisc      = InitInp%KinVisc         ! Kinematic air viscosity
   p%RootName     = InitInp%RootName        ! Rootname for outputs
   call GetPath( p%RootName, rootDir, baseName ) 
   p%VTK_OutFileRoot = trim(rootDir) // 'vtk_fvw'  ! Directory for VTK outputs
   p%VTK_OutFileBase = trim(rootDir) // 'vtk_fvw' // PathSep // trim(baseName) ! Basename for VTK files
   ! Set indexing to AFI tables -- this is set from the AD15 calling code.

   ! Set the Chord values
   allocate(p%W(p%nWings))
   do iW=1,p%nWings
      call AllocAry(p%W(iW)%AFindx, size(InitInp%W(iW)%AFindx,1), 1, 'AFindx',ErrStat,ErrMsg)
      p%W(iW)%AFindx = InitInp%W(iW)%AFindx     ! Copying in case AD15 still needs these
      p%W(iW)%iRotor = InitInp%W(iW)%iRotor    

      p%W(iW)%nSpan  = size(InitInp%W(iW)%chord)-1
      call move_alloc(InitInp%W(iW)%chord, p%W(iW)%chord_LL)
      call AllocAry(p%W(iW)%s_LL    , p%W(iW)%nSpan+1, 'Spanwise coord LL  ', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); p%W(iW)%s_LL= -999999_ReKi;
      call AllocAry(p%W(iW)%s_CP    , p%W(iW)%nSpan  , 'Spanwise coord CPll', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); p%W(iW)%s_CP= -999999_ReKi;
      !call AllocAry(p%W(iW)%chord_LL   , p%W(iW)%nSpan+1, 'Chord on LL        ', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); p%W(iW)%chord_LL= -999999_ReKi;
      call AllocAry(p%W(iW)%chord_CP, p%W(iW)%nSpan  , 'Chord on CP LL     ', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); p%W(iW)%chord_CP= -999999_ReKi;
   enddo

   ! --- Distributing wings to rotors
   p%nRotors = p%W(1)%iRotor
   do iW=2,p%nWings
      p%nRotors = max(p%nRotors,p%W(iW)%iRotor)
   end do
   ! Count number of blades per rotor
   call AllocAry(nBldPerRot, p%nRotors , 'nBldPerRot', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); 
   nBldPerRot=0
   do iW=1,p%nWings
      nBldPerRot(p%W(iW)%iRotor) = nBldPerRot(p%W(iW)%iRotor)+1
   enddo
   nBldMax   = maxval(nBldPerRot)
   ! Set mapping from (rotor,blades) to wings
   call AllocAry(p%Bld2Wings, p%nRotors , nBldMax, 'Bld2Wings', ErrStat2, ErrMsg2);call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); 
   p%Bld2Wings=-1 ! Import init to trigger some index errors
   nBldPerRot=0
   do iW=1,p%nWings
      nBldPerRot(p%W(iW)%iRotor) = nBldPerRot(p%W(iW)%iRotor)+1
      p%Bld2Wings(p%W(iW)%iRotor, nBldPerRot(p%W(iW)%iRotor)) = iW
   enddo

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
   integer(IntKi)       :: iW
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

   do iW=1,p%nWings
      if (allocated(p%W(iW)%PrescribedCirculation)) deallocate(p%W(iW)%PrescribedCirculation)
      if (InputFileData%CirculationMethod==idCircPrescribed) then 
         call AllocAry(p%W(iW)%PrescribedCirculation,  p%W(iW)%nSpan, 'Prescribed Circulation', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_SetParameters');
         p%W(iW)%PrescribedCirculation = -999999_ReKi;
         if (.not. allocated(p%W(iW)%s_CP)) then
            ErrMsg  = 'Spanwise coordinate not allocated.'
            ErrStat = ErrID_Fatal
            return
         endif
         call ReadAndInterpGamma(trim(InputFileData%CirculationFile), p%W(iW)%s_CP(1:p%W(iW)%nSpan), p%W(iW)%s_LL(p%W(iW)%nSpan+1), p%W(iW)%PrescribedCirculation, ErrStat2, ErrMsg2)
         call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_SetParameters' ); 
      endif
   enddo

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
   integer(IntKi) :: nP, nFWEff, iW
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

   nP=0
   do iW=1,p%nWings
      nP = np + (  (p%W(iW)%nSpan+1)*(m%nNW-1+2) +(FWnSpan+1)*(m%nFW+1) )
   enddo
   nFWEff = min(m%nFW, p%nFWFree)
   ! --- Display some status to screen
   if (DEV_VERSION)  print'(A,F10.3,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,L1)','FVW status - t:',t,'  n:',n,'  nNW:',m%nNW-1,'/',p%nNWMax-1,'  nFW:',nFWEff, '+',m%nFW-nFWEff,'=',m%nFW,'/',p%nFWMax,'  nP:',nP, 's Comp:',m%ComputeWakeInduced

   ! --- Evaluation at t
   ! Inputs at t
   call FVW_CopyInput( u(2), uInterp, MESH_NEWCOPY, ErrStat2, ErrMsg2); if(Failed()) return
   call FVW_Input_ExtrapInterp(u(1:size(utimes)),utimes(:),uInterp,t, ErrStat2, ErrMsg2); if(Failed()) return
   call Wings_Panelling(uInterp%WingsMesh, p, m, ErrStat2, ErrMsg2); if(Failed()) return
   call Map_LL_NW(p, m, z, x, 1.0_ReKi, ErrStat2, ErrMsg2); if(Failed()) return ! needed at t=0 if wing moved after init
   call Map_NW_FW(p, m, z, x, ErrStat2, ErrMsg2); if(Failed()) return

   ! Compute UA inputs at t
   if (m%UA_Flag) then
      call CalculateInputsAndOtherStatesForUA(1, uInterp, p, x, xd, z, OtherState, AFInfo, m, ErrStat2, ErrMsg2); if(Failed()) return
   end if

   ! --- Integration between t and t+DTfvw
   if (m%ComputeWakeInduced) then

      ! TODO TODO: this should be in CCSD, but memory is changing between time steps, so for now we have to use u(1)..
      ! inputs: V_wind, output: set m%W%Vwnd_NW, m%W%Vwnd_FW
      CALL DistributeRequestedWind_NWFW(u(1)%V_wind, p, m)

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
         !! Inputs at t+DTfvw (Wings Panelling updates CP, and VstW(iW)%r_LL) 
         !call FVW_Input_ExtrapInterp(u(1:size(utimes)),utimes,uInterp,t+p%DTfvw, ErrStat2, ErrMsg2); if(Failed()) return
         !call Wings_Panelling(uInterp%WingsMesh, p, m, ErrStat2, ErrMsg2); if(Failed()) return
         !! Updating positions of first NW and FW panels (Circulation also updated but irrelevant)
         !call Map_LL_NW(p, m, z, m%x2, 1.0, ErrStat2, ErrMsg2); if(Failed()) return
         !call Map_NW_FW(p, m, z, m%x2, ErrStat2, ErrMsg2); if(Failed()) return
         !! --- Solve for quasi steady circulation at t+p%DTfvw
         !! Returns: z%W(iW)%Gamma_LL (at t+p%DTfvw)
         !z_guess%W(iW)%Gamma_LL = z%W(iW)%Gamma_LL ! We use as guess the circulation from the previous time step (see above)
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

   ! Inputs at t+DTaero (Wings Panelling updates CP, and VstW(iW)%r_LL) 
   call FVW_Input_ExtrapInterp(u(1:size(utimes)),utimes,uInterp,t+p%DTaero, ErrStat2, ErrMsg2); if(Failed()) return
   call Wings_Panelling(uInterp%WingsMesh, p, m, ErrStat2, ErrMsg2); if(Failed()) return

   ! Updating positions of first NW and FW panels (Circulation also updated but irrelevant)
   ! Changes: x only
   ShedScale = (t+p%DTaero - m%OldWakeTime)/p%DTfvw
   call Map_LL_NW(p, m, z, x, ShedScale, ErrStat2, ErrMsg2); if(Failed()) return
   call Map_NW_FW(p, m, z, x, ErrStat2, ErrMsg2); if(Failed()) return
   !call print_x_NW_FW(p, m, x,'Map2')

   ! --- Solve for quasi steady circulation at t+p%DTaero
   ! Returns: z%W(iW)%Gamma_LL (at t+p%DTaero)
   allocate(z_guess%W(p%nWings))
   do iW=1,p%nWings 
      z_guess%W(iW)%Gamma_LL = z%W(iW)%Gamma_LL ! We use as guess the circulation from the previous time step (see above)
   enddo
   call FVW_CalcConstrStateResidual(t+p%DTaero, uInterp, p, x, xd, z_guess, OtherState, m, z, AFInfo, ErrStat2, ErrMsg2, 2); if(Failed()) return
   ! Updating circulation of near wake panel (need to be set for UA, Uind on LL) (and position but irrelevant)
   ! Changes: x only
   call Map_LL_NW(p, m, z, x, ShedScale, ErrStat2, ErrMsg2); if(Failed()) return
   call Map_NW_FW(p, m, z, x, ErrStat2, ErrMsg2); if(Failed()) return
   ! Compute UA inputs at t+DTaero and integrate UA states between t and t+dtAero
   if (m%UA_Flag) then
      call CalculateInputsAndOtherStatesForUA(2, uInterp, p, x, xd, z, OtherState, AFInfo, m, ErrStat2, ErrMsg2); if(Failed()) return
      call UA_UpdateState_Wrapper(AFInfo, t, n, (/t,t+p%DTaero/), p, x, xd, OtherState, m, ErrStat2, ErrMsg2); if(Failed()) return
      ! Compute unsteady Gamma based on UA Cl
      if (p%DStallOnWake .and. p%CirculationMethod/=idCircPrescribed) then 
         call UA_SetGammaDyn(t, uInterp, p, x, xd, OtherState, m, AFInfo, z, ErrStat, ErrMsg)
         ! Updating circulation of near wake panel again (and position but irrelevant)
         ! Changes: x only
         call Map_LL_NW(p, m, z, x, ShedScale, ErrStat2, ErrMsg2); if(Failed()) return
         call Map_NW_FW(p, m, z, x, ErrStat2, ErrMsg2); if(Failed()) return
      end if
   end if



   ! --- Fake handling of ground effect (ensure vorticies above ground)
   call FakeGroundEffect(p, x, m, ErrStat, ErrMsg)

   ! set the wind points required for t+p%DTaero timestep
   CALL SetRequestedWindPoints(m%r_wind, x, p, m)

   if (m%FirstCall) then
      m%FirstCall=.False.
   endif
   call CleanUp()

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

   subroutine CleanUp()
      call FVW_DestroyConstrState(z_guess, ErrStat2, ErrMsg2); if(Failed()) return
   end subroutine

   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FVW_UpdateStates') 
      Failed =  ErrStat >= AbortErrLev
      if (Failed) call CleanUp()
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
   integer(IntKi)       :: i,j,k,iW
   real(ReKi)           :: visc_fact, age  ! Viscosity factor for diffusion of reg param
   real(ReKi), dimension(3) :: VmeanFW, VmeanNW ! Mean velocity of the near wake and far wake

   ErrStat = ErrID_None
   ErrMsg  = ""


   if (.not.allocated(dxdt%W)) then
      allocate(dxdt%W(p%nWings))
      do iW=1,p%nWings
         call AllocAry( dxdt%W(iW)%r_NW , 3   ,  p%W(iW)%nSpan+1  ,p%nNWMax+1, 'Wind on NW ', ErrStat2, ErrMsg2); dxdt%W(iW)%r_NW= -999999_ReKi;
         call AllocAry( dxdt%W(iW)%r_FW , 3   ,  FWnSpan+1  ,p%nFWMax+1, 'Wind on FW ', ErrStat2, ErrMsg2); dxdt%W(iW)%r_FW= -999999_ReKi;
      enddo
      if(Failed()) return
   endif

   ! Distribute the Wind we requested to Inflow wind to storage Misc arrays
   ! TODO ANDY: replace with direct call to inflow wind at W(iW)%r_NW and W(iW)%r_FW locations
   ! NOTE: this has been commented out due to some information missing at some times (and memoery reindexing)
   !       Call to inflow wind sould be done here at actual positions.
   !CALL DistributeRequestedWind_NWFW(u%V_wind, p, m%Vwnd_NW, m%Vwnd_FW)

   ! Only calculate freewake after start time and if on a timestep when it should be calculated.
   if ((t>= p%FreeWakeStart)) then
      nFWEff = min(m%nFW, p%nFWFree)

      ! --- Compute Induced velocities on the Near wake and far wake based on the marker postions:
      ! (expensive N^2 call)
      ! In  : x%W(iW)%r_NW,    r%W(iW)%r_FW 
      ! Out:  m%W(iW)%Vind_NW, m%Vind_FW 
      call WakeInducedVelocities(p, x, m, ErrStat2, ErrMsg2); if(Failed()) return

      ! --- Mean induced velocity over the near wake (NW) TODO, store per wing
      VmeanNW(1:3)=0
      if (m%nNW >1) then
         do iW=1,size(m%W); do j=2,m%nNW+1; do k=1,size(m%W(iW)%Vind_NW,2); 
            VmeanNW(1:3) = VmeanNW(1:3) + m%W(iW)%Vind_NW(1:3, k, j)
         enddo; enddo; enddo; 
         VmeanNW(1:3) = VmeanNW(1:3) / (size(m%W)*m%nNW*size(m%W(1)%Vind_NW,2)) ! TODO TODO
      endif
      ! --- Induced velocity over the free far wake (FWEff)
      VmeanFW(1:3)=0
      if (nFWEff >0) then
         do iW=1,size(m%W); do j=1,nFWEff; do k=1,size(m%W(iW)%Vind_FW,2); 
            VmeanFW(1:3) = VmeanFW(1:3) + m%W(iW)%Vind_FW(1:3, k, j)
         enddo; enddo; enddo; 
         VmeanFW(1:3) = VmeanFW(1:3) / (size(m%W)*nFWEff*size(m%W(1)%Vind_FW,2)) ! TODO TODO
      else
         VmeanFW=VmeanNW
         ! Since we convect the first FW point, we need a reasonable velocity there 
         ! NOTE: mostly needed for sub-cycling and when no FW
         do iW=1,p%nWings
            m%W(iW)%Vind_FW(1, 1:FWnSpan+1, 1) = VmeanFW(1)
            m%W(iW)%Vind_FW(2, 1:FWnSpan+1, 1) = VmeanFW(2)
            m%W(iW)%Vind_FW(3, 1:FWnSpan+1, 1) = VmeanFW(3)
         enddo
      endif

      ! --- Convecting non-free FW with a constant induced velocity (and free stream)
      do iW=1,p%nWings
         m%W(iW)%Vind_FW(1, 1:FWnSpan+1, p%nFWFree+1:p%nFWMax+1) = VmeanFW(1) !
         m%W(iW)%Vind_FW(2, 1:FWnSpan+1, p%nFWFree+1:p%nFWMax+1) = VmeanFW(2) !
         m%W(iW)%Vind_FW(3, 1:FWnSpan+1, p%nFWFree+1:p%nFWMax+1) = VmeanFW(3) !
      enddo

      if (DEV_VERSION) then
         !call print_mean_4d( m%W(iW)%Vind_NW(:,:, 1:m%nNW+1,:), 'Mean induced vel. NW')
         !if (nFWEff>0) then
         !   call print_mean_4d( m%W(iW)%Vind_FW(:,:, 1:nFWEff ,:), 'Mean induced vel. FW')
         !endif
         !print'(A25,3F12.4)','MeanFW (non free)',VmeanFW
         !call print_mean_4d( m%Vwnd_NW(:,:, 1:m%nNW+1,:), 'Mean wind vel.    NW')
         !call print_mean_4d( m%Vwnd_FW(:,:, 1:nFWEff+1,:), 'Mean wind vel. FWEff')
         !call print_mean_4d( m%Vwnd_FW(:,:, (p%nFWFree+1):m%nFW+1,:), 'Mean wind vel.    FWNF')
         !call print_mean_4d( m%Vwnd_FW(:,:, 1:m%nFW+1,:), 'Mean wind vel.    FW')
      endif

      ! --- Vortex points are convected with the free stream and induced velocity
      do iW=1,p%nWings
         dxdt%W(iW)%r_NW(1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1) = m%W(iW)%Vwnd_NW(1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1) +  m%W(iW)%Vind_NW(1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1)
         dxdt%W(iW)%r_FW(1:3, 1:FWnSpan+1, 1:m%nFW+1) = m%W(iW)%Vwnd_FW(1:3, 1:FWnSpan+1, 1:m%nFW+1) +  m%W(iW)%Vind_FW(1:3, 1:FWnSpan+1, 1:m%nFW+1)
      enddo
   else
      if(DEV_VERSION) then
         !call print_mean_4d( m%Vwnd_NW(:,:,1:m%nNW+1,:), 'Mean wind vel.    NW')
         !call print_mean_4d( m%Vwnd_FW(:,:,1:m%nFW+1,:), 'Mean wind vel.    FW')
      endif

      ! --- Vortex points are convected with the free stream
      do iW=1,p%nWings
         dxdt%W(iW)%r_NW(1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1) = m%W(iW)%Vwnd_NW(1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1) 
         dxdt%W(iW)%r_FW(1:3, 1:FWnSpan+1, 1:m%nFW+1) = m%W(iW)%Vwnd_FW(1:3, 1:FWnSpan+1, 1:m%nFW+1)
      enddo
   endif
   ! First NW point does not convect (bound to LL)
   do iW=1,p%nWings
      dxdt%W(iW)%r_NW(1:3, :, 1:p%iNWStart-1)=0.0_ReKi
   enddo
   ! First FW point always convects (even if bound to NW)
   ! This is done for overcycling
   !dxdt%W(iW)%r_FW(1:3, :, 1, :)=0

   ! --- Regularization
   do iW=1,p%nWings
      if (.not.allocated(dxdt%W(iW)%Eps_NW)) then
         call AllocAry( dxdt%W(iW)%Eps_NW , 3   ,  p%W(iW)%nSpan    ,p%nNWMax, 'Eps NW ', ErrStat2, ErrMsg2); 
         call AllocAry( dxdt%W(iW)%Eps_FW , 3   ,  FWnSpan    ,p%nFWMax, 'Eps FW ', ErrStat2, ErrMsg2); 
         if(Failed()) return
      endif
   enddo
   if (p%WakeRegMethod==idRegConstant) then
      do iW=1,p%nWings
         dxdt%W(iW)%Eps_NW(1:3, :, :)=0.0_ReKi
         dxdt%W(iW)%Eps_FW(1:3, :, :)=0.0_ReKi
      enddo

   else if (p%WakeRegMethod==idRegStretching) then
      ! TODO
   else if (p%WakeRegMethod==idRegAge) then
      visc_fact = 2.0_ReKi * CoreSpreadAlpha * p%CoreSpreadEddyVisc * p%KinVisc
      ! --- Method 1, use d(rc^2)/dt = 4 k 
      do iW=1,p%nWings
         dxdt%W(iW)%Eps_NW(1:3, :,p%iNWStart:) = visc_fact/x%W(iW)%Eps_NW(1:3, :, p%iNWStart:)
         dxdt%W(iW)%Eps_FW(1:3, :,          :) = visc_fact/x%W(iW)%Eps_FW(1:3, :, :)
         ! --- Method 2, use rc(tau) = 2k/sqrt(r_c^2(tau=0) + 4 k tau)
         !dxdt%W(iW)%Eps_NW(1:3, :, :, :) = (visc_fact)/sqrt(x%W(iW)%Eps_NW(1:3, :, :, :)**2 + 2*visc_fact*p%DTaero)
         !dxdt%W(iW)%Eps_FW(1:3, :, :, :) = (visc_fact)/sqrt(x%W(iW)%Eps_FW(1:3, :, :, :)**2 + 4*visc_fact*p%DTaero)
      enddo
   else
      ErrStat = ErrID_Fatal
      ErrMsg ='Regularization method not implemented'
   endif
   do iW=1,p%nWings
      dxdt%W(iW)%Eps_NW(1:3,:,1:p%iNWStart) = 0.0_ReKi ! Important! LL and First NW panel epsilon does not change
   enddo

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
   integer :: iW
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

   do iW=1,p%nWings
      x%W(iW)%r_NW     = (1_ReKi-fact) * states(1)%W(iW)%r_NW     + fact * states(2)%W(iW)%r_NW
      x%W(iW)%r_FW     = (1_ReKi-fact) * states(1)%W(iW)%r_FW     + fact * states(2)%W(iW)%r_FW
      x%W(iW)%Eps_NW   = (1_ReKi-fact) * states(1)%W(iW)%Eps_NW   + fact * states(2)%W(iW)%Eps_NW
      x%W(iW)%Eps_FW   = (1_ReKi-fact) * states(1)%W(iW)%Eps_FW   + fact * states(2)%W(iW)%Eps_FW
      x%W(iW)%Gamma_NW = (1_ReKi-fact) * states(1)%W(iW)%Gamma_NW + fact * states(2)%W(iW)%Gamma_NW
      x%W(iW)%Gamma_FW = (1_ReKi-fact) * states(1)%W(iW)%Gamma_FW + fact * states(2)%W(iW)%Gamma_FW
   enddo
   !print*,'fact',fact,states(1)%W(iW)%Gamma_NW(29,iNWStart+1,1),x%W(iW)%Gamma_NW(29,iNWStart+1,1),states(2)%W(iW)%Gamma_NW(29,iNWStart+1,1)
   !print*,'fact',fact,states(1)%W(iW)%r_NW(1,29,iNWStart+1,1),x%W(iW)%r_NW(1,29,iNWStart+1,1),states(2)%W(iW)%r_NW(1,29,iNWStart+1,1)

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
   integer :: iW
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = "" 

   dt = real(p%DTfvw,ReKi)  ! NOTE: this is DTfvw
   ! Compute "right hand side"
   CALL FVW_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, m%dxdt, ErrStat2, ErrMsg2); if (Failed()) return

   ! Update of positions and reg param
   do iW = 1, p%nWings
      x%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1) = x%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1) + dt * m%dxdt%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1)
      x%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan  , 1:m%nNW  ) = x%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan  , 1:m%nNW  ) + dt * m%dxdt%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan  , 1:m%nNW  )
      if ( m%nFW>0) then
         x%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1) = x%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1) + dt * m%dxdt%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1)
         x%W(iW)%Eps_FW(1:3, 1:FWnSpan  , 1:m%nFW  ) = x%W(iW)%Eps_FW(1:3, 1:FWnSpan  , 1:m%nFW  ) + dt * m%dxdt%W(iW)%Eps_FW(1:3, 1:FWnSpan  , 1:m%nFW  )
      endif
   enddo
   ! Update of Gamma TODO (viscous diffusion, stretching)


   if (DEV_VERSION) then
      do iW = 1, p%nWings
         ! Additional checks
         if (any(m%dxdt%W(iW)%r_NW(1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1)<-999)) then
            print*,'FVW_Euler1: Attempting to convect NW with a wrong velocity'
            STOP
         endif
         if ( m%nFW>0) then
            if (any(m%dxdt%W(iW)%r_FW(1:3, 1:FWnSpan+1, 1:m%nFW+1)<-999)) then
               call print_x_NW_FW(p, m, x, 'STP')
               print*,'FVW_Euler1: Attempting to convect FW with a wrong velocity'
               STOP
            endif
         endif
         if (any(m%dxdt%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan, 1:m%nNW)<-0)) then
            print*,'FVW_Euler1: Wrong Epsilon NW'
            STOP
         endif
         if ( m%nFW>0) then
            if (any(m%dxdt%W(iW)%Eps_FW(1:3, 1:FWnSpan, 1:m%nFW)<-999)) then
               print*,'FVW_Euler1: Wrong Epsilon FW'
               STOP
            endif
         endif
      enddo
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
      integer :: iW
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
         do iW = 1,p%nWings
            if (any(m%dxdt%W(iW)%r_NW(1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1)<-999)) then
               print*,'FVW_RK4: Attempting to convect NW with a wrong velocity'
               STOP
            endif
            if ( m%nFW>0) then
               if (any(m%dxdt%W(iW)%r_FW(1:3, 1:FWnSpan+1, 1:m%nFW+1)<-999)) then
                  call print_x_NW_FW(p, m, x, 'STP')
                  print*,'FVW_RK4: Attempting to convect FW with a wrong velocity'
                  STOP
               endif
            endif
         enddo
      endif

      do iW = 1,p%nWings
         k1%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1) = dt * m%dxdt%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1) 
         k1%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan  , 1:m%nNW  ) = dt * m%dxdt%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan,   1:m%nNW  )
         if ( m%nFW>0) then   
            k1%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1) = dt * m%dxdt%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1)
            k1%W(iW)%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW  ) = dt * m%dxdt%W(iW)%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW  )
         endif

         x_tmp%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1) = x%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1) + 0.5 * k1%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1)
         x_tmp%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan  , 1:m%nNW  ) = x%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan  , 1:m%nNW  ) + 0.5 * k1%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan,   1:m%nNW  )
         if ( m%nFW>0) then
            x_tmp%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1) = x%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1)  + 0.5 * k1%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1)
            x_tmp%W(iW)%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW  ) = x%W(iW)%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW  )  + 0.5 * k1%W(iW)%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW  )
         endif
      enddo

      ! interpolate u to find u_interp = u(t + dt/2)
      CALL FVW_Input_ExtrapInterp(u(1:size(utimes)),utimes(:),u_interp, t+0.5*dt, ErrStat2, ErrMsg2); CALL CheckError(ErrStat2,ErrMsg2); IF ( ErrStat >= AbortErrLev ) RETURN        

      ! find dxdt at t + dt/2
      CALL FVW_CalcContStateDeriv( t + 0.5*dt, u_interp, p, x_tmp, xd, z, OtherState, m, m%dxdt, ErrStat2, ErrMsg2 ); CALL CheckError(ErrStat2,ErrMsg2); IF ( ErrStat >= AbortErrLev ) RETURN

      do iW = 1,p%nWings
         k2%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1) = dt * m%dxdt%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1)
         k2%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan  , 1:m%nNW  ) = dt * m%dxdt%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan,   1:m%nNW  )
         if ( m%nFW>0) then
            k2%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1) = dt * m%dxdt%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1)
            k2%W(iW)%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW  ) = dt * m%dxdt%W(iW)%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW  )
         endif

         x_tmp%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1) = x%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1) + 0.5 * k2%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1)
         x_tmp%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan  , 1:m%nNW  ) = x%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan  , 1:m%nNW  ) + 0.5 * k2%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan,   1:m%nNW  )
         if ( m%nFW>0) then
            x_tmp%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1) = x%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1)  + 0.5 * k2%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1)
            x_tmp%W(iW)%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW  ) = x%W(iW)%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW  )  + 0.5 * k2%W(iW)%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW  )
         endif
      enddo

      ! find dxdt at t + dt/2       
      CALL FVW_CalcContStateDeriv( t + 0.5*dt, u_interp, p, x_tmp, xd, z, OtherState, m, m%dxdt, ErrStat2, ErrMsg2 ); CALL CheckError(ErrStat2,ErrMsg2); IF ( ErrStat >= AbortErrLev ) RETURN

      do iW = 1,p%nWings
         k3%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1) = dt * m%dxdt%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1)
         k3%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan  , 1:m%nNW  ) = dt * m%dxdt%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan,   1:m%nNW)
         if ( m%nFW>0) then
            k3%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1) = dt * m%dxdt%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1)
            k3%W(iW)%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW  ) = dt * m%dxdt%W(iW)%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW)
         endif

         x_tmp%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1) = x%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1) + k3%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1)
         x_tmp%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan  , 1:m%nNW  ) = x%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan  , 1:m%nNW  ) + k3%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan,   1:m%nNW)
         if ( m%nFW>0) then
            x_tmp%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1) = x%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1)  + k3%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1)
            x_tmp%W(iW)%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW  ) = x%W(iW)%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW  )  + k3%W(iW)%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW)
         endif
      enddo

      ! interpolate u to find u_interp = u(t + dt)
      CALL FVW_Input_ExtrapInterp(u(1:size(utimes)),utimes(:),u_interp, t + dt, ErrStat2, ErrMsg2); CALL CheckError(ErrStat2,ErrMsg2); IF ( ErrStat >= AbortErrLev ) RETURN

      ! find dxdt at t + dt
      CALL FVW_CalcContStateDeriv( t + dt, u_interp, p, x_tmp, xd, z, OtherState, m, m%dxdt, ErrStat2, ErrMsg2 ); CALL CheckError(ErrStat2,ErrMsg2); IF ( ErrStat >= AbortErrLev ) RETURN

      do iW = 1,p%nWings
         k4%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1) = dt * m%dxdt%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1)
         k4%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan  , 1:m%nNW  ) = dt * m%dxdt%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan,   1:m%nNW)
         if ( m%nFW>0) then
            k4%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1) = dt * m%dxdt%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1)
            k4%W(iW)%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW  ) = dt * m%dxdt%W(iW)%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW)
         endif
      enddo

      ! Compute and store combined dx = (k1/6 + k2/3 + k3/3 + k4/6) ! NOTE: this has dt, it's not a true dxdt yet
      do iW = 1,p%nWings
         m%dxdt%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1) = ( k1%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1) + 2._ReKi * k2%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1) + 2._ReKi * k3%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1) + k4%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1)  ) / 6._ReKi
         m%dxdt%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan  , 1:m%nNW  ) = ( k1%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan  , 1:m%nNW  ) + 2._ReKi * k2%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan  , 1:m%nNW  ) + 2._ReKi * k3%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan  , 1:m%nNW  ) + k4%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan  , 1:m%nNW  )  ) / 6._ReKi
         if ( m%nFW>0) then         
            m%dxdt%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1) = ( k1%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1) + 2._ReKi * k2%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1) + 2._ReKi * k3%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1) + k4%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1)  ) / 6._ReKi
            m%dxdt%W(iW)%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW  ) = ( k1%W(iW)%Eps_FW(1:3, 1:FWnSpan  , 1:m%nFW  ) + 2._ReKi * k2%W(iW)%Eps_FW(1:3, 1:FWnSpan  , 1:m%nFW  ) + 2._ReKi * k3%W(iW)%Eps_FW(1:3, 1:FWnSpan  , 1:m%nFW  ) + k4%W(iW)%Eps_FW(1:3, 1:FWnSpan  , 1:m%nFW  )  ) / 6._ReKi
         endif
      enddo

      !update positions
      do iW = 1,p%nWings
         x%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1) = x%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1) + m%dxdt%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1) 
         x%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan  , 1:m%nNW  ) = x%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan  , 1:m%nNW  ) + m%dxdt%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan  , 1:m%nNW  ) 
         if ( m%nFW>0) then         
            x%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1) = x%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1) + m%dxdt%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1)
            x%W(iW)%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW  ) = x%W(iW)%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW  ) + m%dxdt%W(iW)%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW)
         endif
      enddo

      ! Store true dxdt =  (k1/6 + k2/3 + k3/3 + k4/6)/dt 
      do iW = 1,p%nWings
         m%dxdt%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1) = m%dxdt%W(iW)%r_NW  (1:3, 1:p%W(iW)%nSpan+1, 1:m%nNW+1)/dt
         m%dxdt%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan  , 1:m%nNW  ) = m%dxdt%W(iW)%Eps_NW(1:3, 1:p%W(iW)%nSpan  , 1:m%nNW  )/dt
         if ( m%nFW>0) then         
            m%dxdt%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1) = m%dxdt%W(iW)%r_FW  (1:3, 1:FWnSpan+1, 1:m%nFW+1)/dt
            m%dxdt%W(iW)%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW  ) = m%dxdt%W(iW)%Eps_FW(1:3, 1:FWnSpan,   1:m%nFW  )/dt
         endif
      enddo

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
   integer :: iW

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Distribute the Wind we requested to Inflow wind to storage Misc arrays
   ! TODO ANDY: replace with direct call to inflow wind at m%W(iW)%CP location
   ! input: V_wind, output: m%W%Vwnd_LL
   CALL DistributeRequestedWind_LL(u%V_wind, p, m)

   ! Solve for the residual of the constraint state functions here:
   !z%residual = 0.0_ReKi
   !z%W(iW)%Gamma_LL = 0.0_ReKi
   allocate(z_out%W(p%nWings))
   do iW= 1,p%nWings
      call AllocAry( z_out%W(iW)%Gamma_LL,  p%W(iW)%nSpan, 'Lifting line Circulation', ErrStat, ErrMsg );
      z_out%W(iW)%Gamma_LL = -999999_ReKi;
   enddo

   CALL Wings_ComputeCirculation(t, z_out, z_guess, u, p, x, m, AFInfo, ErrStat, ErrMsg, iLabel)

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
   integer(IntKi) :: iW
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   character(*), parameter       :: RoutineName = 'FVW_CalcOutput'

   ErrStat = ErrID_None
   ErrMsg  = ""
!       ! --- NOTE: this below might not be needed
!       ! Distribute the Wind we requested to Inflow wind to storage Misc arrays
!       ! TODO ANDY: replace with direct call to inflow wind at m%W(iW)%CP location
!       CALL DistributeRequestedWind_LL(u%V_wind, p, m%Vwnd_LL)
! 
!       ! Control points location and structural velocity
   call Wings_Panelling(u%WingsMesh, p, m, ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
! 
!    ! if we are on a correction step, CalcOutput may be called again with different inputs
!    ! Compute m%W(iW)%Gamma_LL
!    CALL Wings_ComputeCirculation(t, m%W(iW)%Gamma_LL, z%W(iW)%Gamma_LL, u, p, x, m, AFInfo, ErrStat2, ErrMsg2, 0); if(Failed()) return ! For plotting only

   !--- Induction on the lifting line control point
   ! if     InductionAtCP : In: m%W%CP,  Out:m%W%Vind_CP                 and m%W%Vind_LL (averaged)
   ! if not InductionAtCP : In: m%W%r_LL,   Out:m%W%Vind_CP (interp/extrap) and m%W%Vind_LL
   if (p%Induction) then
      call LiftingLineInducedVelocities(p, x, p%InductionAtCP, 1, m, ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      ! Transfer to output
      do iW=1,p%nWings
          y%W(iW)%Vind(1,:) = m%W(iW)%Vind_LL(1,:)
          y%W(iW)%Vind(2,:) = m%W(iW)%Vind_LL(2,:)
          y%W(iW)%Vind(3,:) = m%W(iW)%Vind_LL(3,:)
      enddo
   else
      do iW=1,p%nWings
          y%W(iW)%Vind(1:3,:) = 0.0_ReKi
      enddo
   endif
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
   logical                 :: bTimeToOutput
   logical                 :: bWithinTime
   integer(IntKi)          :: ErrStat2
   character(ErrMsgLen)    :: ErrMsg2
   character(*), parameter :: RoutineName = 'FVW_CalcOutput'
   integer(IntKi) :: iW, iGrid
   integer(IntKi) :: nSeg, nSegP
   ErrStat = ErrID_None
   ErrMsg  = ''

   ! --- Write VTK of wake/blade filaments/panels
   if (p%WrVTK>0) then
      if (m%FirstCall .or. force) then
         call MKDIR(p%VTK_OutFileRoot)
      endif
      ! For plotting only
      call PackPanelsToSegments(p, x, 1, (p%ShearModel==idShearMirror), m%nNW, m%nFW, m%Sgmt%Connct, m%Sgmt%Points, m%Sgmt%Gamma, m%Sgmt%Epsilon, nSeg, nSegP)
      do iW=1,p%nWings
         m%W(iW)%Vtot_CP = m%W(iW)%Vind_CP + m%W(iW)%Vwnd_CP - m%W(iW)%Vstr_CP
      enddo
      if ( force .or. (( t - m%VTKlastTime ) >= p%DTvtk*OneMinusEpsilon ))  then
         m%VTKlastTime = t
         if ((p%VTKCoord==2).or.(p%VTKCoord==3)) then
            ! Hub reference coordinates, for export only, ALL VTK Will be exported in this coordinate system!
            ! Note: hubOrientation and HubPosition are optional, but required for bladeFrame==TRUE
            if (size(u%rotors)>1) then
               ErrMsg2='VTK outputs in hub coordinates not implemented for multiple rotors'
               call SetErrStat(ErrID_Warn, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            endif
            ! We ouput in first rotor frame
            call WrVTK_FVW(p, x, z, m, trim(p%VTK_OutFileBase)//'FVW_Hub', m%VTKStep, 9, bladeFrame=.TRUE.,  &
                     HubOrientation=real(u%rotors(1)%HubOrientation,ReKi),HubPosition=real(u%rotors(1)%HubPosition,ReKi))
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
         bWithinTime   = t>=m%GridOutputs(iGrid)%tStart-p%DTaero/2. .and. t<= m%GridOutputs(iGrid)%tEnd+p%DTaero/2.
         bTimeToOutput = ( t - m%GridOutputs(iGrid)%tLastOutput) >= m%GridOutputs(iGrid)%DTout * OneMinusEpsilon
         if (force .or. (bWithinTime .and. bTimeToOutput) )  then
            ! Compute induced velocity on grid, TODO use the same Tree for all CalcOutput
            call InducedVelocitiesAll_OnGrid(m%GridOutputs(iGrid), p, x, m, ErrStat2, ErrMsg2);
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
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
   integer                :: i,iW
   integer(intKi)         :: ErrStat2
   character(ErrMsgLen)   :: ErrMsg2
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   m%UA_Flag=InitInp%UA_Flag
   ! --- Condensed version of BEMT_Init_Otherstate
   if ( m%UA_Flag ) then

      ! We store these per wings (they each contains element info for (nNodes x 1)
      allocate(x%UA(p%nWings), xd%UA(p%nWings), OtherState%UA(p%nWings))

      do iW=1,p%nWings
         ! ---Condensed version of "BEMT_Set_UA_InitData"
         allocate(Init_UA_Data%c(InitInp%numBladeNodes,1), STAT = errStat2)
         do i = 1,InitInp%numBladeNodes
            Init_UA_Data%c(i,1)      = p%W(iW)%chord_LL(i) ! NOTE: InitInp chord move-allocd to p
         end do
         Init_UA_Data%dt              = interval          
         Init_UA_Data%OutRootName     = trim(InitInp%RootName)//'W'//num2lstr(iW)
         Init_UA_Data%numBlades       = 1
         Init_UA_Data%nNodesPerBlade  = InitInp%numBladeNodes ! At AeroDyn ndoes, not CP

         Init_UA_Data%UAMod           = InitInp%UAMod  
         Init_UA_Data%Flookup         = InitInp%Flookup
         Init_UA_Data%a_s             = InitInp%a_s ! Speed of sound, m/s  
         Init_UA_Data%ShedEffect      = .False. ! Important, when coupling UA wih vortex code, shed vorticity is inherently accounted for
         Init_UA_Data%WrSum           = InitInp%SumPrint
         allocate(Init_UA_Data%UAOff_innerNode(1), stat=errStat2)
         allocate(Init_UA_Data%UAOff_outerNode(1), stat=errStat2)
         Init_UA_Data%UAOff_innerNode(1) = InitInp%W(iW)%UAOff_innerNode
         Init_UA_Data%UAOff_outerNode(1) = InitInp%W(iW)%UAOff_outerNode

         ! --- UA init
         allocate(m%W(iW)%u_UA(InitInp%numBladeNodes, 2), stat=errStat2) 
         call UA_Init( Init_UA_Data, m%W(iW)%u_UA(1,1), m%W(iW)%p_UA, x%UA(iW), xd%UA(iW), OtherState%UA(iW), m%W(iW)%y_UA, m%W(iW)%m_UA, interval, AFInfo, p%W(iW)%AFIndx, InitOutData_UA, ErrStat2, ErrMsg2); if(Failed())return

         call UA_DestroyInitInput( Init_UA_Data, ErrStat2, ErrMsg2 ); if(Failed())return
         call UA_DestroyInitOutput( InitOutData_UA, ErrStat2, ErrMsg2 ); if(Failed())return

      enddo

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
   type(UA_InputType), pointer     :: u_UA ! Alias to shorten notations
   integer(IntKi)                                    :: i,iW
   character(ErrMsgLen)                              :: errMsg2     ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)                                    :: errStat2    ! temporary Error status of the operation
   ErrStat  = ErrID_None
   ErrMsg   = ""
   !NOTE: UA happens at the LL nodes (different from the Control Points)

   ! --- Induction on the lifting line control points
   ! if     InductionAtCP : In: m%W%CP,     Out:m%W%Vind_CP                 and m%W%Vind_LL (averaged)
   ! if not InductionAtCP : In: m%W%r_LL,   Out:m%W%Vind_CP (interp/extrap) and m%W%Vind_LL
   if (p%Induction) then
      call LiftingLineInducedVelocities(p, x, p%InductionAtCP, 1, m, ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'CalculateInputsAndOtherStatesForUA'); if (ErrStat >= AbortErrLev) return
   else
      do iW = 1,p%nWings  
         m%W(iW)%Vind_LL(1:3,:)=0.0_ReKi
      enddo
   endif
   if (p%InductionAtCP) then
      if (m%nNW<=1) then
         do iW = 1,p%nWings  
            m%W(iW)%Vind_LL(1:3,:)=0.0_ReKi
         enddo
      endif
   endif
   ! --- UA inputs
   do iW = 1,p%nWings  
      do i = 1,p%W(iW)%nSpan+1 
         ! We only update the UnsteadyAero states if we have unsteady aero turned on for this node      
         u_UA => m%W(iW)%u_UA(i,InputIndex) ! Alias
         !! ....... compute inputs to UA ...........
         ! NOTE: To be consistent with CalcOutput we take Vwind_LL that was set using m%DisturbedInflow from AeroDyn.. 
         ! This is not clean, but done to be consistent, waiting for AeroDyn to handle UA
         call AlphaVrel_Generic(u%WingsMesh(iW)%Orientation(1:3,1:3,i), u%WingsMesh(iW)%TranslationVel(1:3,i),  m%W(iW)%Vind_LL(1:3,i), u%W(iW)%Vwnd_LL(1:3,i), &
                                 p%KinVisc, p%W(iW)%chord_LL(i), u_UA%U, u_UA%alpha, u_UA%Re)
         u_UA%v_ac(1)  = sin(u_UA%alpha)*u_UA%U
         u_UA%v_ac(2)  = cos(u_UA%alpha)*u_UA%U
         u_UA%omega    = u%W(iW)%omega_z(i)
         u_UA%UserProp = 0 ! u1%UserProp(i,j) ! TODO
      end do ! i nSpan
   end do ! iW nWings
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
      do i = 1,p%W(j)%nSpan+1 
         ! COMPUTE: x%UA, xd%UA
         call UA_UpdateStates( i, 1, t, n, m%W(j)%u_UA(i,:), uTimes, m%W(j)%p_UA, x%UA(j), xd%UA(j), OtherState%UA(j), AFInfo(p%W(j)%AFIndx(i,1)), m%W(j)%m_UA, errStat2, errMsg2 )
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

!> Set dynamic gamma based on dynamic stall states
!! NOTE: We use Vind_LL computed in CalculateInputsAndOtherStatesForUA
subroutine UA_SetGammaDyn(t, u, p, x, xd, OtherState, m, AFInfo, z, ErrStat, ErrMsg)
   use UnsteadyAero, only: UA_CalcOutput
   real(DbKi),                    intent(in   ) :: t           !< Curent time
   type(FVW_InputType),           intent(in   ) :: u          !< Inputs at Time t
   type(FVW_ParameterType),       intent(in   ) :: p          !< AD parameters
   type(FVW_ContinuousStateType), intent(in )   :: x          !< continuous states
   type(FVW_DiscreteStateType),   intent(in   ) :: xd         !< Discrete states
   type(FVW_OtherStateType),      intent(in   ) :: OtherState !< OtherState
   type(FVW_MiscVarType),target,  intent(inout) :: m          !< Misc/optimization variables
   type(AFI_ParameterType ),      intent(in   ) :: AFInfo(:)  !< The airfoil parameter data, temporary, for UA..
   type(FVW_ConstraintStateType), intent(inout) :: z          !< Constraint states
   integer(IntKi),                intent(  out) :: ErrStat    !< Error status of the operation
   character(*),                  intent(  out) :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   real(ReKi)                  :: Cl_dyn, Cl_dyn_prev, Cl_dyn_avg
   real(ReKi)                  :: Gamma_dyn, Gamma_dyn_prev, Gamma_dyn_avg
   type(UA_InputType), pointer :: u_UA ! Alias to shorten notations
   integer(IntKi), parameter   :: InputIndex=2 ! we will always use values at t+dt in this routine
   integer(intKi)              :: iW, j ! loop counter on wings and nodes
   integer(intKi)              :: errStat2
   character(ErrMsgLen)        :: errMsg2

   ErrStat = 0
   ErrMsg = ""

   do iW=1,p%nWings
      ! Gamma_LL is expressed at CP, so we average the dynamic gamma from both nodes
      ! NOTE: this is inconsistent with the Wings solving which occurs at the CPs
      j=1
      u_UA => m%W(iW)%u_UA(j,InputIndex) ! Alias
      call UA_CalcOutput(j, 1, t, u_UA, m%W(iW)%p_UA, x%UA(iW), xd%UA(iW), OtherState%UA(iW), AFInfo(p%W(iW)%AFindx(j,1)), m%W(iW)%y_UA, m%W(iW)%m_UA, errStat2, errMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'UA_SetGammaDyn')
      Gamma_dyn_prev = 0.5_ReKi * m%W(iW)%y_UA%Cl * u_UA%U * p%W(iW)%chord_LL(j)
      do j = 2,p%W(iW)%nSpan+1
         u_UA => m%W(iW)%u_UA(j,InputIndex) ! Alias
         call UA_CalcOutput(j, 1, t, u_UA, m%W(iW)%p_UA, x%UA(iW), xd%UA(iW), OtherState%UA(iW), AFInfo(p%W(iW)%AFindx(j,1)), m%W(iW)%y_UA, m%W(iW)%m_UA, errStat2, errMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'UA_SetGammaDyn')
         Gamma_dyn      = 0.5_ReKi * m%W(iW)%y_UA%Cl * u_UA%U * p%W(iW)%chord_LL(j)
         Gamma_dyn_avg  = (Gamma_dyn+Gamma_dyn_prev)*0.5_ReKi
         Gamma_dyn_prev = Gamma_dyn
         z%W(iW)%Gamma_LL(j-1) = Gamma_dyn_avg
      enddo
   enddo ! iW, Loop on wings
end subroutine UA_SetGammaDyn

end module FVW
