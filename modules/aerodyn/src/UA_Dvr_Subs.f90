module UA_Dvr_Subs
   
   use NWTC_Library
   use AirfoilInfo
   use AirfoilInfo_Types
   use UnsteadyAero_Types
   use UnsteadyAero
   use LinDyn
   use LinDyn_Types
   
   implicit none

   integer, parameter        :: NumAFfiles = 1
   integer(IntKi), parameter :: NumInp = 2           ! Number of inputs sent to UA_UpdateStates (must be at least 2)
   integer(IntKi), parameter :: InflowMod_Cst  = 1   ! Inflow is constant
   integer(IntKi), parameter :: InflowMod_File = 2   ! Inflow is read from file
   integer(IntKi), parameter, dimension(2) :: InflowMod_Valid  = (/InflowMod_Cst, InflowMod_File/)
   integer(IntKi), parameter :: MotionMod_Cst  = 1   ! Motion is constant
   integer(IntKi), parameter :: MotionMod_File = 2   ! Motion is read from file
   integer(IntKi), parameter, dimension(2) :: MotionMod_Valid  = (/MotionMod_Cst, MotionMod_File/)
   real(ReKi), parameter     :: myNaN = -9999.9_ReKi
   integer(IntKi), parameter :: idFmt_Ascii  = 1
   integer(IntKi), parameter :: idFmt_Binary = 2
   integer(IntKi), parameter :: idFmt_Both   = 3
   integer(IntKi), parameter, dimension(3) :: idFmt_Valid  = (/idFmt_Ascii, idFmt_Binary, idFmt_Both/)

   type Dvr_Parameters
      logical         :: Echo            
      ! Environment
      real(ReKi)      :: KinVisc
      real(ReKi)      :: FldDens
      real(ReKi)      :: SpdSound        
      ! 
      integer         :: UAMod           
      logical         :: Flookup        
      logical         :: UseCm         
      character(1024) :: AirFoil1 
      real(ReKi)      :: Chord
      ! 
      integer         :: SimMod          
      ! Reduced frequency - SimMod = 1
      real(ReKi)      :: InflowVel       
      real(ReKi)      :: NCycles         
      real(ReKi)      :: Frequency 
      real(ReKi)      :: Re
      integer         :: StepsPerCycle
      real(ReKi)      :: Amplitude       
      real(ReKi)      :: Mean            
      integer         :: Phase           
      ! Prescribed Aero - SimMod = 2
      real(ReKi)      :: TMax_PA
      real(DbKi)      :: dt_PA
      character(1024) :: AeroTSFile      
      ! AeroElastic Section - SimMod =3
      real(ReKi)      :: TMax
      real(DbKi)      :: dt
      real(ReKi)      :: MM(3,3)
      real(ReKi)      :: CC(3,3)
      real(ReKi)      :: KK(3,3)
      logical         :: activeDOFs(3)
      real(ReKi)      :: GFScaling(3,3)
      real(ReKi)      :: initPos(3)
      real(ReKi)      :: initVel(3)
      real(ReKi)      :: Vec_AQ(2)  ! Vector from A to quarter chord /aerodynamic center
      real(ReKi)      :: Vec_AT(2)  ! Vector from A to three quarter chord
      real(ReKi)      :: Twist      ! Twist of the airfoil section (input deg, but stored in rad afterwards)
      ! Inflow
      integer         :: InflowMod = InflowMod_Cst
      real(ReKi)      :: Inflow(2)
      character(1024) :: InflowTSFile
      ! Motion
      integer         :: MotionMod = MotionMod_Cst
      character(1024) :: MotionTSFile
      ! Outputs
      logical         :: SumPrint         
      logical         :: WrAFITables
      ! ---- Parameters
      real(ReKi)                             :: d_34_to_ac  
      !real(DbKi)                             :: dt
      real(DbKi)                             :: simTime  
      integer                                :: numSteps
      character(1024)                        :: OutRootName ! Automatically obtained from input file name
      ! Prescribed AoA simulations
      real(DbKi), allocatable                :: timeArr(:)
      real(ReKi), allocatable                :: vPrescrAero(:,:)   ! Aero as function of time, shape nt x 4:  Time, AOA, U, Omega
      ! Prescribed inflow simulations
      real(ReKi), allocatable                :: vU0(:,:)      ! Inflow as function of time, shape nt x 3 : Time, U0x, U0y
   end type Dvr_Parameters


   type :: Dvr_Outputs
      integer(intki)                                 :: unOutFile = -1        !< unit number for writing output file
      !integer(intki)                                :: actualchanlen         !< actual length of channels written to text file (less than or equal to chanlen) [-]
      integer(intki)                                :: ny                     !< total number of outputs for the driver 
      integer(intki)                                :: ny_dvr                 !< number of outputs for the driver (without UA and LD, and Time)
      integer(intki)                                :: ny_UA                  !< number of outputs for UA
      integer(intki)                                :: ny_LD                  !< number of outputs for LD
      !character(20)                                 :: fmt_t                 !< format specifier for time channel [-]
      !character(25)                                 :: fmt_a                 !< format specifier for each column (including delimiter) [-]
      !character(1)                                  :: delim                 !< column delimiter [-]
      !character(20)                                 :: outfmt                !< format specifier [-]
      integer(intki)                                 :: fileFmt = idFmt_Binary !< output format 1=text, 2=binary, 3=both [-]
      character(1024)                                :: root  = ''            !< output file rootname [-]
      character(ChanLen) , dimension(:), allocatable :: WriteOutputHdr        !< channel headers [-]
      character(ChanLen) , dimension(:), allocatable :: WriteOutputUnt        !< channel units [-]
      real(ReKi) , dimension(:,:), allocatable       :: storage               !< nchannel x ntime [-]
      real(ReKi) , dimension(:), allocatable         :: outline               !< output line to be written to disk [-]
      !real(dbki)  :: dt_outs      !< output time resolution [s]
      !integer(intki)  :: n_dt_out      !< number of time steps between writing a line in the time-marching output files [-]
   end type Dvr_Outputs

   type :: Dvr_Misc
      ! Reminder:
      ! Q: 1/4 chord / aerodynamic center
      ! T: 3/4 chord
      ! A: Airfoil origin
      real(ReKi) :: Vst_Q(2)        !< Structural velocity   at Q [m/s]
      real(ReKi) :: Vst_T(2)        !< Structural velocity   at T [m/s]
      real(ReKi) :: Vrel_Q(2)       !< Relative velocity     at Q [m/s]
      real(ReKi) :: Vrel_T(2)       !< Relative velocity     at T [m/s]
      real(ReKi) :: Vrel_norm2_T    !< Squared velocity norm at T [m^2/s^2]
      real(ReKi) :: Vrel_norm2_Q    !< Squared velocity norm at Q [m^2/s^2]
      real(ReKi) :: alpha_Q         !< Angle of attack       at Q [rad]
      real(ReKi) :: alpha_T         !< Angle of attack       at T [rad]
      real(ReKi) :: phi_Q           !< Flow angle            at Q [rad]
      real(ReKi) :: phi_T           !< Flow angle            at T [rad]
      real(ReKi) :: Re              !< Reynolds number (NOT in Million!)
      real(ReKi) :: L, D, tau_Q     !< Aerodynamic loads     at Q [N/m & Nm/m]
      real(ReKi) :: FxA, FyA, tau_A !< Aerodynamic loads     at A [N/m & Nm/m] 
      real(ReKi) :: GF(3)           !< Generalized force, Scaled aerodynamic loads to be representative of the blade
      real(ReKi) :: twist_full      !< Full twist (includes initial twist, potential pitch, and torsion)
      integer    :: iU0Last = 1     !< Index for faster interpolation of wind speed
      integer    :: iPALast = 1     !< Index for faster interpolation of prescribed aero
      real(ReKi) :: uPA(3)          !< Prescribed Aero inputs
   end type Dvr_Misc

   type Dvr_Data
      ! Time control
      real(DbKi)                             :: uTimes(NumInp)
      ! Parameters / initinp  set as the same...
      type(Dvr_Parameters)                   :: p ! Initialization/parameter data for the driver program
      type(Dvr_Misc)                         :: m ! Misc variables for aerodynamic calculations
      ! Outputs
      type(Dvr_Outputs)                      :: out
      ! Inflow
      real(ReKi)                             :: U0(NumInp, 2)           ! Inflow velocity vector at time t and t+dt
      ! AFI
      type(AFI_ParameterType)                :: AFI_Params(NumAFfiles)
      integer, allocatable                   :: AFIndx(:,:)
      ! UA
      type(UA_InitInputType)                 :: UA_InitInData           ! Input data for initialization
      type(UA_InitOutputType)                :: UA_InitOutData          ! Output data from initialization
      type(UA_ContinuousStateType)           :: UA_x                    ! Continuous states
      type(UA_DiscreteStateType)             :: UA_xd                   ! Discrete states
      type(UA_OtherStateType)                :: UA_OtherState           ! Other/optimization states
      type(UA_MiscVarType)                   :: UA_m                    ! Misc/optimization variables
      type(UA_ParameterType)                 :: UA_p                    ! Parameters
      type(UA_InputType)                     :: UA_u(NumInp)            ! System inputs
      type(UA_OutputType)                    :: UA_y                    ! System outputs
      ! Dynamics
      type(LD_InitInputType)                 :: LD_InitInData           ! Input data for initialization
      type(LD_InitOutputType)                :: LD_InitOutData          ! Output data from initialization
      type(LD_ContinuousStateType)           :: LD_x                    ! Continuous states
      type(LD_DiscreteStateType)             :: LD_xd                   ! Discrete states
      type(LD_OtherStateType)                :: LD_OtherState           ! Other/optimization states
      type(LD_ConstraintStateType)           :: LD_z                    ! Constraint states
      type(LD_MiscVarType)                   :: LD_m                    ! Misc/optimization variables
      type(LD_ParameterType)                 :: LD_p                    ! Parameters
      type(LD_InputType)                     :: LD_u(NumInp)            ! System inputs
      type(LD_OutputType)                    :: LD_y                    ! System outputs
      !
      type(LD_ContinuousStateType)           :: LD_x_swp                ! Continuous states
      type(LD_OtherStateType)                :: LD_OtherState_swp       ! Other/optimization states
      type(UA_ContinuousStateType)           :: UA_x_swp                ! Continuous states
      type(UA_DiscreteStateType)             :: UA_xd_swp               ! Discrete states
      type(UA_OtherStateType)                :: UA_OtherState_swp       ! Other/optimization states
   end type Dvr_Data

contains
   
!--------------------------------------------------------------------------------------------------------------
subroutine ReadDriverInputFile( FileName, InitInp, ErrStat, ErrMsg )
   character(1024),               intent( in    )   :: filename
   type(Dvr_Parameters),          intent(   out )   :: InitInp
   integer,                       intent(   out )   :: ErrStat              ! returns a non-zero value when an error occurs  
   character(*),                  intent(   out )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   ! Local variables  
   integer                 :: UnEcho   ! The local unit number for this module's echo file
   integer                 :: iLine
   character(1024)         :: EchoFile ! Name of HydroDyn echo file
   character(1024)         :: PriPath                         ! the path to the primary input file
   character(1024)         :: Line                         ! the path to the primary input file
   type(FileInfoType)      :: FI       !< The derived type for holding the file information.
   integer(IntKi)          :: errStat2 ! Status of error message
   character(1024)         :: errMsg2  ! Error message if ErrStat /= ErrID_None
   character(*), parameter :: RoutineName = 'ReadDriverInputFile'
   ! Initialize the echo file unit to -1 which is the default to prevent echoing, we will alter this based on user input
   UnEcho  = -1
   ErrStat = ErrID_None
   ErrMsg  = ''

   ! Read all input file lines into fileinfo
   call WrScr(' Opening UnsteadyAero Driver input file: '//trim(FileName) )
   call ProcessComFile(FileName, FI, errStat2, errMsg2); if (Failed()) return
   CALL GetPath( FileName, PriPath )    ! Input files will be relative to the path where the primary input file is located.
   !call GetRoot(FileName, dvr%root)      

   ! --- Header and echo
   iLine = 3 ! Skip the first two lines as they are known to be header lines and separators
   call ParseVar(FI, iLine, 'Echo', InitInp%Echo, errStat2, errMsg2); if (Failed()) return;
   if ( InitInp%Echo ) then
      EchoFile = trim(FileName)//'.ech'
      call OpenEcho (UnEcho, EchoFile, errStat2, errMsg2 ); if(Failed()) return
      do iLine = 1, 3
         write(UnEcho, '(A)') trim(FI%Lines(iLine))
      enddo
   end if

   iLine = 4
   ! --- Environmental conditions section
   call ParseCom(FI, iLine, Line                        , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'FldDens',  InitInp%FldDens , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'KinVisc',  InitInp%KinVisc , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'SpdSound', InitInp%SpdSound, errStat2, errMsg2, UnEcho); if(Failed()) return

   ! --- UNSTEADYAERO section
   call ParseCom(FI, iLine, Line                              , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'UAMod'      , InitInp%UAMod      , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'Flookup'    , InitInp%Flookup    , errStat2, errMsg2, UnEcho); if(Failed()) return
   
   ! --- AIRFOIL PROPERTIES section
   call ParseCom(FI, iLine, Line                        , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'AirFoil' , InitInp%AirFoil1, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'Chord'   , InitInp%Chord   , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'Vec_AQ'       , InitInp%Vec_AQ    , 2, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'Vec_AT'       , InitInp%Vec_AT    , 2, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'UseCm'   , InitInp%UseCm   , errStat2, errMsg2, UnEcho); if(Failed()) return
   
   ! --- SIMULATION CONTROL section
   call ParseCom(FI, iLine, Line                                  , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'SimMod'       , InitInp%SimMod       , errStat2, errMsg2, UnEcho); if(Failed()) return

   ! --- REDUCED FREQUENCY
   call ParseCom(FI, iLine, Line                                  , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'InflowVel'  , InitInp%InflowVel  , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'NCycles'      , InitInp%NCycles      , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'StepsPerCycle', InitInp%StepsPerCycle, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'Frequency'    , InitInp%Frequency    , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'Amplitude'    , InitInp%Amplitude    , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'Mean'         , InitInp%Mean         , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'Phase'        , InitInp%Phase        , errStat2, errMsg2, UnEcho); if(Failed()) return

   ! --- PRESCRIBED AERO section
   call ParseCom(FI, iLine, Line                                  , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'TMax_PA'      , InitInp%Tmax_PA      , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'DT_PA'        , InitInp%dt_PA        , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'AeroTSFile'   , InitInp%AeroTSFile   , errStat2, errMsg2, UnEcho); if(Failed()) return

   ! --- ELASTIC SECTION section
   call ParseCom(FI, iLine, Line                        , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'TMax'         , InitInp%Tmax         , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'DT'           , InitInp%dt           , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'activeDOFs'   , InitInp%activeDOFs, 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'initPos'      , InitInp%initPos   , 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'initVel'      , InitInp%initVel   , 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'GFScaling1'   , InitInp%GFScaling(1,:) , 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'GFScaling2'   , InitInp%GFScaling(2,:) , 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'GFScaling3'   , InitInp%GFScaling(3,:) , 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'MassMatrix1'  , InitInp%MM(1,:)   , 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'MassMatrix2'  , InitInp%MM(2,:)   , 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'MassMatrix3'  , InitInp%MM(3,:)   , 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'DampMatrix1'  , InitInp%CC(1,:)   , 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'DampMatrix2'  , InitInp%CC(2,:)   , 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'DampMatrix3'  , InitInp%CC(3,:)   , 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'StifMatrix1'  , InitInp%KK(1,:)   , 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'StifMatrix2'  , InitInp%KK(2,:)   , 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'StifMatrix3'  , InitInp%KK(3,:)   , 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'Twist'        , InitInp%Twist     ,    errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'InflowMod'    , InitInp%InflowMod ,    errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'Inflow'       , InitInp%Inflow    , 2, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'InflowTSFile' , InitInp%InflowTSFile,  errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'MotionMod'    , InitInp%MotionMod ,    errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'MotionTSFile' , InitInp%MotionTSFile,  errStat2, errMsg2, UnEcho); if(Failed()) return

   ! --- OUTPUT section
   call ParseCom(FI, iLine, Line                        , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'SumPrint'   , InitInp%SumPrint   , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'WrAFITables', InitInp%WrAFITables, errStat2, errMsg2, UnEcho); if(Failed()) return

   ! --- Triggers
   call GetRoot(FileName, InitInp%OutRootName) ! OutRootName is inferred from current filename.
   !InitInp%OutRootName=trim(InitInp%OutRootName)//'.UA' ! For backward compatibility
   !if (PathIsRelative(InitInp%OutRootName)) InitInp%OutRootName = TRIM(PriPath)//TRIM(InitInp%OutRootName)
   if (PathIsRelative(InitInp%Airfoil1))    InitInp%Airfoil1 = TRIM(PriPath)//TRIM(InitInp%Airfoil1)
   if (PathIsRelative(InitInp%AeroTSFile   )) InitInp%AeroTSFile   = TRIM(PriPath)//TRIM(InitInp%AeroTSFile  )
   if (PathIsRelative(InitInp%InflowTSFile )) InitInp%InflowTSFile = TRIM(PriPath)//TRIM(InitInp%InflowTSFile)
   if (PathIsRelative(InitInp%MotionTSFile )) InitInp%MotionTSFile = TRIM(PriPath)//TRIM(InitInp%MotionTSFile)

   ! --- Checks
   !if (Check(.not.(any(dvr%out%fileFmt==idFmt_Valid   )),   'FileFormat not implemented: '//trim(Num2LStr(InitInp%InflowMod)))) return
   if (Check(.not.(any(InitInp%InflowMod==InflowMod_Valid)), 'InflowMod not implemented: '//trim(Num2LStr(InitInp%MotionMod)))) return
   if (Check(.not.(any(InitInp%MotionMod==MotionMod_Valid)), 'MotionMod not implemented: '//trim(Num2LStr(InitInp%MotionMod)))) return
   if (InitInp%SimMod==3) then ! Temporary to avoid changing r-test for now
      !if (Check(.not.EqualRealNos(InitInp%MM(1,1), InitInp%MM(2,2), 'Mass matrix entries 11 and 22 should match.') return

      if (InitInp%Vec_AT(2)<0) call WrScr('[WARN] Vec_AT(2) is negative, but this value is usually positive (for A between T and Q)')
      if (InitInp%Vec_AQ(2)>0) call WrScr('[WARN] Vec_AQ(2) is positive, but this value is usually negative (for A between T and Q)')
   endif

   call Cleanup()
contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call Cleanup()
   end function Failed

   subroutine Cleanup()
      ! Close this module's echo file    
      if ( InitInp%Echo ) then
         close(UnEcho)
      end if
      Call NWTC_Library_Destroyfileinfotype(FI, errStat2, errMsg2)
   end subroutine Cleanup

   logical function Check(Condition, errMsg_in)
        logical, intent(in) :: Condition
        character(len=*), intent(in) :: errMsg_in
        Check=Condition
        if (Check) then
           call SetErrStat( ErrID_Fatal, errMsg_in, errStat, errMsg, RoutineName )
        endif
   end function Check

end subroutine ReadDriverInputFile
!--------------------------------------------------------------------------------------------------------------
subroutine Dvr_SetParameters(p, errStat, errMsg)
   type(Dvr_Parameters), intent(inout) :: p
   integer,              intent(out ) :: errStat              ! returns a non-zero value when an error occurs
   character(*),         intent(out ) :: errMsg               ! Error message if ErrStat /= ErrID_None
   integer(IntKi)  :: errStat2    ! Status of error message
   character(1024) :: errMsg2     ! Error message if ErrStat /= ErrID_None
   errStat2= ErrID_None
   errStat = ErrID_None
   errMsg  = ''
   ! Unit conversions
   p%Twist  = p%Twist * D2R
   p%d_34_to_ac = (-p%Vec_AQ(2) + p%Vec_AT(2)) !  d_34_to_ac = d_QT ~0.5 [-], Approximated using y coordinate
   p%Vec_AT = p%Vec_AT * p%chord
   p%Vec_AQ = p%Vec_AQ * p%chord

   if ( p%SimMod == 1 ) then
      call WrScr('[WARN] The behavior of SimMod=1 might change in the future.')

      ! We will use a constant Reynolds..
      p%Re = p%InflowVel * p%chord/ p%KinVisc  ! NOT IN MILLIONS
      print*,'    Re     ',p%Re
      ! Using the frequency and NCycles, determine how long the simulation needs to run
      p%simTime   = p%NCycles/p%Frequency
      p%numSteps  = p%StepsPerCycle*p%NCycles  ! we could add 1 here to make this a complete cycle
      p%dt        = p%simTime / p%numSteps
      
   else if ( p%SimMod == 2 ) then
      ! Read time-series data file with columns:( time,  Angle-of-attack, Vrel, omega )
      call WrScr( ' Opening prescribed-aero time-series input file: '//trim(p%AeroTSFile) )
      call ReadDelimFile(p%AeroTSFile, 4, p%vPrescrAero, errStat2, errMsg2); if(Failed()) return
      p%vPrescrAero(:,2) = p%vPrescrAero(:,2)*D2R ! Deg 2 rad
      p%dt       = p%dt_PA
      p%simTime  = p%TMax_PA
      p%numSteps = int(p%simTime/p%dt)

   elseif ( p%SimMod == 3 ) then
      p%simTime   = p%TMax
      p%numSteps = int(p%simTime/p%dt)

      if (p%InflowMod==InflowMod_File) then
         ! Read inflow file
         call ReadDelimFile(p%InflowTSFile, 3, p%vU0, errStat2, errMsg2); if(Failed()) return
      endif

   endif

   if(Failed()) return
contains
   logical function Failed()
      call setErrStat(errStat2, errMsg2, errStat, errMsg, 'Dvr_SetParameters')
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine Dvr_SetParameters
!--------------------------------------------------------------------------------------------------------------
subroutine driverInputsToUAInitData(p, InitInData, AFI_Params, AFIndx, errStat, errMsg)
   type(Dvr_Parameters)   , intent(in ) :: p           ! Initialization data for the driver program
   type(UA_InitInputType) , intent(out) :: InitInData           ! Input data for initialization
   type(AFI_ParameterType), intent(out) :: AFI_Params(NumAFfiles)
   integer, allocatable   , intent(out) :: AFIndx(:,:)
   integer(IntKi),          intent(out) :: errStat                       ! Error status.
   character(*),            intent(out) :: errMsg                        ! Error message.
   character(1024)         :: afNames(NumAFfiles)
   integer(IntKi)          :: errStat2    ! Status of error message
   character(1024)         :: errMsg2     ! Error message if ErrStat /= ErrID_None
   character(*), parameter  :: RoutineName = 'driverInputsToUAInitData'
   
   errStat     = ErrID_None
   errMsg      = ''
   
   InitInData%UA_OUTS = 1  ! 0=None, 1=Write Outputs, 2=Separate File
#ifdef ADD_UA_OUTS
   InitInData%UA_OUTS = 2 ! Compiler Flag Override,  2=Write a separate file
#endif
   

   ! -- UA Init Input Data
   InitInData%nNodesPerBlade  = 1 
   InitInData%numBlades       = 1
   call AllocAry(InitInData%c, InitInData%nNodesPerBlade, InitInData%numBlades, 'chord', errStat2, errMsg2); if(Failed()) return
   call AllocAry(InitInData%UAOff_innerNode             , InitInData%numBlades, 'UAO'  , errStat2, errMsg2); if(Failed()) return
   call AllocAry(InitInData%UAOff_outerNode             , InitInData%numBlades, 'UAO'  , errStat2, errMsg2); if(Failed()) return

   ! don't turn off UA based on span location:
   InitInData%UAOff_innerNode = 0
   InitInData%UAOff_outerNode = InitInData%nNodesPerBlade + 1
   InitInData%a_s          = p%SpdSound
   InitInData%c(1,1)       = p%Chord
   InitInData%UAMod        = p%UAMod 
   InitInData%IntegrationMethod = UA_Method_ABM4
   InitInData%Flookup      = p%Flookup
   InitInData%OutRootName  = trim(p%OutRootName)//'.UA'
   InitInData%WrSum        = p%SumPrint 
   InitInData%d_34_to_ac   = p%d_34_to_ac !  d_34_to_ac = d_QT ~0.5 [-], Approximated using y coordinate

   ! --- AFI
   allocate(AFIndx(InitInData%nNodesPerBlade,InitInData%numBlades), STAT = errStat2)
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error trying to allocate InitInData%AFIndx.', errStat, errMsg, RoutineName)
      return
   end if
   AFIndx(1,1) = 1

   afNames(1)  = p%AirFoil1 ! All nodes/blades are using the same 2D airfoil
   call Init_AFI( InitInData%UAMod, NumAFfiles, afNames, p%UseCm, AFI_Params, errStat2, errMsg2); if(Failed()) return

   if (p%WrAFITables) then
      call AFI_WrTables(AFI_Params(1), InitInData%UAMod, p%OutRootName)
   endif
contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
endsubroutine driverInputsToUAInitData
!--------------------------------------------------------------------------------------------------------------
subroutine Init_AFI(UAMod, NumAFfiles, afNames, UseCm, AFI_Params, ErrStat, ErrMsg)
   integer,             intent(in   )  :: UAMod
   integer,             intent(in   )  :: NumAFfiles
   CHARACTER(1024),     intent(in   )  :: afNames(NumAFfiles)
   logical,             intent(in   )  :: UseCm
   type(AFI_ParameterType), intent(  out)  :: AFI_Params(NumAFfiles)
   integer(IntKi),      intent(  out)  :: ErrStat                       ! Error status.
   character(*),        intent(  out)  :: ErrMsg                        ! Error message.

   type(AFI_InitInputType)  :: AFI_InitInputs
   integer                  :: UnEc
   integer                  :: i
   integer(IntKi)           :: errStat2    ! Status of error message
   character(1024)          :: errMsg2     ! Error message if ErrStat /= ErrID_None
   character(*), parameter  :: RoutineName = 'Init_AFI'

   UnEc = 0
   ErrStat = ErrID_None
   ErrMsg  = ""
      ! Set this to 1 to use the UA coefs
   !AFI_InitInputs%UA_Model    = 1
      ! This is the number of columns of coefs in the AOA table: Cl, Cd, Cm, for example, but doesn't include Alpha
   !AFI_InitInputs%NumCoefs    = 3
      !
   AFI_InitInputs%InCol_Alfa  = 1
   AFI_InitInputs%InCol_Cl    = 2
   AFI_InitInputs%InCol_Cd    = 3
   if (UseCm) then
      AFI_InitInputs%InCol_Cm    = 4
   else   
      AFI_InitInputs%InCol_Cm    = 0
   end if

   AFI_InitInputs%InCol_Cpmin = 0
   AFI_InitInputs%AFTabMod    = AFITable_1 ! 1D-interpolation (on AoA only)
   AFI_InitInputs%UAMod  = UAMod  ! We calculate some of the UA coefficients based on UA Model
   
   do i=1,NumAFfiles
      AFI_InitInputs%FileName = afNames(i) !InitInp%AF_File(i)
      
         ! Read in and process the airfoil files.
         ! This includes creating the spline coefficients to be used for interpolation.

      call AFI_Init ( AFI_InitInputs, AFI_Params(i), errStat2, errMsg2, UnEc )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if
   end do

   call Cleanup()

contains
   subroutine Cleanup()
      call AFI_DestroyInitInput(AFI_InitInputs, errStat2, errMsg2)
   end subroutine Cleanup
end subroutine Init_AFI  
!--------------------------------------------------------------------------------------------------------------
!> Set Inflow inputs
subroutine setInflow(t, p, U0, m)
   real(DbKi),               intent(in)    :: t
   type(Dvr_Parameters),     intent(in)    :: p
   type(Dvr_Misc      ),     intent(inout) :: m
   real(ReKi), dimension(:), intent(out)   :: U0
   if      (p%InflowMod == InflowMod_Cst) then
      U0(:) = p%Inflow
   else if (p%InflowMod == InflowMod_File) then
      call interpTimeValue(p%vU0, t, m%iU0Last, U0(:))
   else
      print*,'Should never happen'
      STOP
   endif
end subroutine setInflow
!--------------------------------------------------------------------------------------------------------------
!> Compute aerodynamic kinematics quantities (velocities and angles) at different points
subroutine AeroKinematics(U0, q, qd, p, m)
   real(ReKi),           intent(in)    :: U0(2)    !< Free stream
   real(ReKi),           intent(in)    :: q(3)     !< Elastic positions  x,y,th
   real(ReKi),           intent(in)    :: qd(3)    !< Elastic velocities
   type(Dvr_Parameters), intent(in   ) :: p        !< Parameters
   type(Dvr_Misc),       intent(inout) :: m        !< Misc aero var
   real(ReKi), parameter :: W(2) =0 ! Induced velocities
   real(ReKi) :: ST, CT

   ! Full twist
   m%twist_full = q(3) + p%Twist ! + Pitch if a controller is added
   ST = sin(m%twist_full)
   CT = cos(m%twist_full)

   ! Structual velocity including torsional velocity
   m%Vst_T(1) = qd(1) + qd(3) * (-p%Vec_AT(1)*ST + p%Vec_AT(2)*CT)
   m%Vst_T(2) = qd(2) + qd(3) * ( p%Vec_AT(1)*CT + p%Vec_AT(2)*ST)

   m%Vst_Q(1) = qd(1) + qd(3) * (-p%Vec_AQ(1)*ST + p%Vec_AQ(2)*CT)
   m%Vst_Q(2) = qd(2) + qd(3) * ( p%Vec_AQ(1)*CT + p%Vec_AQ(2)*ST)

   ! Relative velocity, Vrel = U0 - Vst + W
   m%Vrel_T = U0 - m%Vst_T + W
   m%Vrel_Q = U0 - m%Vst_Q + W

   ! Squared velocity norm
   m%Vrel_norm2_T = m%Vrel_T(1)**2 + m%Vrel_T(2)**2
   m%Vrel_norm2_Q = m%Vrel_Q(1)**2 + m%Vrel_Q(2)**2

   ! Flow angle
   m%phi_Q = atan2(m%Vrel_Q(1), m%Vrel_Q(2))
   m%phi_T = atan2(m%Vrel_T(1), m%Vrel_T(2))

   ! Angle of attack
   m%alpha_Q = m%phi_Q - m%twist_full
   m%alpha_T = m%phi_T - m%twist_full

   ! Reynolds at 1/4 chord
   m%Re = sqrt(m%Vrel_norm2_Q) * p%chord  / p%KinVisc
end subroutine AeroKinematics

!--------------------------------------------------------------------------------------------------------------
!> Compute aerodynamic kinetics quantities (loads)
subroutine AeroKinetics(U0, q, qd, C_dyn, p, m)
   real(ReKi),           intent(in)    :: U0(2)    !< Free stream
   real(ReKi),           intent(in)    :: q(3)     !< Elastic positions  x,y,th
   real(ReKi),           intent(in)    :: qd(3)    !< Elastic velocities
   real(ReKi),           intent(in)    :: C_dyn(3) !< Dynamic aerodynamic coefficients (Cl, Cd, Cm)
   type(Dvr_Parameters), intent(in   ) :: p        !< Parameters
   type(Dvr_Misc),       intent(inout) :: m        !< Misc aero var
   real(ReKi) :: ST, CT
   real(ReKi) :: SP, CP
   real(ReKi) :: q_dyn

   ! First get kinematics
   call AeroKinematics(U0, q, qd, p, m)

   ST = sin(m%twist_full)
   CT = cos(m%twist_full)

   ! Loads at Q
   q_dyn = 0.5_ReKi * p%FldDens * p%chord * m%Vrel_norm2_Q
   m%L     = q_dyn * C_dyn(1)
   m%D     = q_dyn * C_dyn(2)
   m%tau_Q = q_dyn * C_dyn(3) * p%chord

   ! Loads at A
   SP = sin(m%phi_Q)
   CP = cos(m%phi_Q)
   m%FxA   =  m%L * CP + m%D * SP
   m%FyA   = -m%L * SP + m%D * CP
   ! Tau A (Positive about "z") - version 1
   m%tau_A = m%tau_Q 
   m%tau_A = m%tau_A - m%FxA * (- p%Vec_AQ(1) * ST + p%Vec_AQ(2) * CT) 
   m%tau_A = m%tau_A + m%FyA * (  p%Vec_AQ(1) * CT + p%Vec_AQ(2) * ST) 
   ! Tau A (Positive about "z") - version 2
   !SA = sin(m%alpha_Q)
   !CA = cos(m%alpha_Q)
   !tau_A2 = m%tau_Q 
   !tau_A2 = tau_A2  - q_dyn *C_dyn(1)* ( p%Vec_AQ(1) * SA + p%Vec_AQ(2) * CA) 
   !tau_A2 = tau_A2  + q_dyn *C_dyn(2)* ( p%Vec_AQ(1) * CA - p%Vec_AQ(2) * SA) 

   ! Generalized loads
   m%GF(1) =  m%FxA    * p%GFScaling(1,1) +   m%FyA * p%GFScaling(1,2) - m%tau_A  * p%GFScaling(1,3) 
   m%GF(2) =  m%FxA    * p%GFScaling(2,1) +   m%FyA * p%GFScaling(2,2) - m%tau_A  * p%GFScaling(2,3) 
   m%GF(3) =  m%FxA    * p%GFScaling(3,1) +   m%FyA * p%GFScaling(3,2) - m%tau_A  * p%GFScaling(3,3) 

end subroutine AeroKinetics
!----------------------------------------------------------------------------------------------------  

!> Set LinDyn inputs (scaled aerodynamic forces at point A)
subroutine setLDinputs(U0, LD_x, UA_y, p, m, LD_u)
   real(ReKi)          ,         intent(in   ) :: U0(2) !< Parameters
   type(LD_ContinuousStateType), intent(in   ) :: LD_x  !< LinDyn states
   type(UA_OutputType),          intent(in   ) :: UA_y  !< UA outputs
   type(Dvr_Parameters),         intent(in   ) :: p     !< Parameters
   type(Dvr_Misc),               intent(inout) :: m     !< Misc aero var
   type(LD_InputType),           intent(inout) :: LD_u  !< LinDyn inputs

   call AeroKinetics  (U0, LD_x%q(1:3), LD_x%q(4:6), (/UA_y%Cl, UA_y%Cd, UA_y%Cm/), p, m)
   LD_u%Fext(1) = m%GF(1)
   LD_u%Fext(2) = m%GF(2)
   LD_u%Fext(3) = m%GF(3)

end subroutine setLDinputs

!----------------------------------------------------------------------------------------------------  
!> Set UA Inputs from Flow and LinDyn
subroutine setUAinputs(U0, LD_x, p, m, UA_u)
   real(ReKi)          ,         intent(in   ) :: U0(2) !< Parameters
   type(LD_ContinuousStateType), intent(in   ) :: LD_x  !< LinDyn states
   type(Dvr_Parameters),         intent(in   ) :: p     !< Parameters
   type(Dvr_Misc),               intent(inout) :: m     !< Misc aero var
   type(UA_InputType),           intent(inout) :: UA_u  !< UA inputs

   call AeroKinematics(U0, LD_x%q(1:3), LD_x%q(4:6), p, m)
   UA_u%UserProp = 0
   UA_u%Re       = m%Re
   UA_u%omega    =  -LD_x%q(6) ! NOTE: theta convention for the driver is negative along z, but UA expect an omega along z
   ! Angle of attack and relative velocity at 1/4 point/aerodynamic center point "Q"
   UA_u%alpha    = m%alpha_Q
   UA_u%U        = sqrt(m%Vrel_norm2_Q)
   UA_u%v_ac(1)  = UA_u%U * sin(UA_u%alpha) ! In airfoil coordinate system (a)
   UA_u%v_ac(2)  = UA_u%U * cos(UA_u%alpha) ! In airfoil coordinate system (a)
end subroutine setUAinputs

!----------------------------------------------------------------------------------------------------  
!> Set UA inptus for a simulation where the angle of attack is prescribed and the relative velocity is constant
subroutine setUAinputsAlphaSim(n, u, t, p, m, errStat, errMsg)
   integer,                intent(in)              :: n
   type(UA_InputType),     intent(inout)           :: u            ! System inputs
   real(DbKi),             intent(  out)           :: t
   type(Dvr_Parameters),   intent(in)              :: p           ! Initialization data for the driver program
   type(Dvr_Misc      ),   intent(inout)           :: m           ! Initialization data for the driver program
   integer,                intent(out)             :: errStat
   character(len=*),       intent(out)             :: errMsg
   real(ReKi) :: phase
   real(ReKi) :: d_ref2AC
   real(ReKi) :: alpha_ref
   real(ReKi) :: U_ref
   real(ReKi) :: v_ref(2)
   real(ReKi) :: v_34(2)
   logical, parameter :: OscillationAtMidChord=.true.  ! for legacy, use false
   logical, parameter :: VelocityAt34         =.true.  ! for legacy, use false

   ! Initialize error handling variables
   ErrMsg  = ''
   ErrStat = ErrID_None

   u%UserProp = 0
   t = (n-1)*p%dt

   if ( p%SimMod == 1 ) then
      if (OscillationAtMidChord) then
         d_ref2AC = -0.25_ReKi  ! -0.25: oscillations at mid_chord
         d_ref2AC = -p%d_34_to_ac/2. ! TODO
      else
         d_ref2AC = 0.0_ReKi   ! 0: oscillations at AC
      endif
      U_ref = p%InflowVel  ! m/s

      phase = (n+p%Phase-1)*2*pi/p%StepsPerCycle
      alpha_ref = (p%Amplitude * sin(phase) + p%Mean)*D2R   ! This needs to be in radians
      v_ref(1) = sin(alpha_ref)*U_ref
      v_ref(2) = cos(alpha_ref)*U_ref
      u%omega =  p%Amplitude * cos(phase) * 2*pi/p%StepsPerCycle / p%dt * D2R  ! This needs to be in radians derivative: d_alpha /d_t

      u%v_ac(1) = v_ref(1) + u%omega * d_ref2AC* p%Chord
      u%v_ac(2) = v_ref(2)

      u%alpha = atan2(u%v_ac(1), u%v_ac(2) )  ! 
      if (VelocityAt34) then
         v_34(1) = u%v_ac(1) + u%omega * 0.5* p%Chord
         v_34(2) = u%v_ac(2)

         u%U =  sqrt(v_34(1)**2 + v_34(2)**2) ! Using U at 3/4
      else
         u%U =  sqrt(u%v_ac(1)**2 + u%v_ac(2)**2) ! Using U at 1/4
      endif
      u%Re = p%Re ! Option for constant Reynolds or not?

   else
      ! Interpolate at current time
      call interpTimeValue(p%vPrescrAero, t, m%iPALast, m%uPA)
      u%alpha = m%uPA(1) ! rad
      u%U     = m%uPA(2)
      u%omega = m%uPA(3)
      u%v_ac(1) = sin(u%alpha)*u%U
      u%v_ac(2) = cos(u%alpha)*u%U
      u%Re = u%U * p%chord  / p%KinVisc
   end if

end subroutine setUAinputsAlphaSim

!----------------------------------------------------------------------------------------------------------------------------------
subroutine Dvr_EndSim(dvr, errStat, errMsg)
   type(Dvr_Data), target,  intent(inout) :: dvr       ! driver data
   integer(IntKi)         , intent(out)   :: errStat              ! Status of error message
   character(*)           , intent(out)   :: errMsg               ! Error message if errStat /= ErrID_None
   character(ErrMsgLen)       :: errMsg2                 ! temporary Error message if errStat /= ErrID_None
   integer(IntKi)             :: errStat2                ! temporary Error status of the operation
   character(*), parameter    :: RoutineName = 'Dvr_EndSim'
   type(Dvr_Outputs), pointer :: out       ! driver output, data
   out => dvr%out
   errStat = ErrID_None
   errMsg  = ''
   ! Close the output file
   if (out%fileFmt==idFmt_Both .or. out%fileFmt == idFmt_Ascii) then
      if (out%unOutFile > 0) close(out%unOutFile)
   endif
   if (out%fileFmt==idFmt_Both .or. out%fileFmt == idFmt_Binary) then
      call WrScr(' Writing output file: '//trim(out%Root)//'.outb')
      call WrBinFAST(trim(out%Root)//'.outb', FileFmtID_ChanLen_In, 'AeroDynDriver', out%WriteOutputHdr, out%WriteOutputUnt, (/0.0_DbKi, dvr%p%dt/), out%storage(:,:), errStat2, errMsg2)
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
   endif
end subroutine Dvr_EndSim



   
   
! --------------------------------------------------------------------------------
! --- IO 
! --------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
!> Concatenate new output channels info to the extisting ones in the driver
!! TODO COPY PASTED FROM AeroDyn_Inflow. Should be placed in NWTC_Lib NWTC_Str
subroutine concatOutputHeaders(WriteOutputHdr0, WriteOutputUnt0, WriteOutputHdr, WriteOutputUnt, errStat, errMsg)
   character(ChanLen), dimension(:), allocatable, intent(inout) ::  WriteOutputHdr0 !< Channel headers
   character(ChanLen), dimension(:), allocatable, intent(inout) ::  WriteOutputUnt0 !< Channel units
   character(ChanLen), dimension(:), allocatable, intent(inout) ::  WriteOutputHdr !< Channel headers
   character(ChanLen), dimension(:), allocatable, intent(inout) ::  WriteOutputUnt !< Channel units
   integer(IntKi)              , intent(  out) :: errStat       !< Status of error message
   character(*)                , intent(  out) :: errMsg        !< Error message if errStat /= ErrID_None
   ! Locals
   character(ChanLen), allocatable :: TmpHdr(:)
   character(ChanLen), allocatable :: TmpUnt(:)
   integer :: nOld, nAdd
   errStat = ErrID_None
   errMsg  = ''
   if (.not.allocated(WriteOutputHdr)) return
   if (.not.allocated(WriteOutputHdr0)) then
      call move_alloc(WriteOutputHdr, WriteOutputHdr0)
      call move_alloc(WriteOutputUnt, WriteOutputUnt0)   
   else
      nOld = size(WriteOutputHdr0)
      nAdd = size(WriteOutputHdr)

      call move_alloc(WriteOutputHdr0, TmpHdr)
      call move_alloc(WriteOutputUnt0, TmpUnt)   

      allocate(WriteOutputHdr0(nOld+nAdd))
      allocate(WriteOutputUnt0(nOld+nAdd))
      WriteOutputHdr0(1:nOld) = TmpHdr
      WriteOutputUnt0(1:nOld) = TmpUnt
      WriteOutputHdr0(nOld+1:nOld+nAdd) = WriteOutputHdr
      WriteOutputUnt0(nOld+1:nOld+nAdd) = WriteOutputUnt
      deallocate(TmpHdr)
      deallocate(TmpUnt)
   endif
end subroutine concatOutputHeaders
!----------------------------------------------------------------------------------------------------------------------------------
!> Initialize outputs to file for driver 
subroutine Dvr_InitializeOutputs(out, numSteps, errStat, errMsg)
   type(Dvr_Outputs),        intent(inout)   :: out 
   integer(IntKi)         ,  intent(in   )   :: numSteps             ! Number of time steps
   integer(IntKi)         ,  intent(  out)   :: errStat              ! Status of error message
   character(*)           ,  intent(  out)   :: errMsg               ! Error message if errStat /= ErrID_None
   ! locals
   integer(IntKi)     :: numOuts
!       integer(IntKi)     :: i
!       integer(IntKi)     :: numSpaces
!       integer(IntKi)     :: iWT
!       character(ChanLen) :: colTxt
!       character(ChanLen) :: caseTxt
! 
   numOuts = size(out%WriteOutputHdr)

   call AllocAry(out%outLine, numOuts-1, 'outLine', errStat, errMsg); ! NOTE: time not stored
   out%outLine=0.0_ReKi
! 
!       ! --- Ascii
!       if (out%fileFmt==idFmt_Both .or. out%fileFmt == idFmt_Ascii) then
! 
!          ! compute the width of the column output
!          numSpaces = out%ActualChanLen ! the size of column produced by OutFmt
!          out%ActualChanLen = max( out%ActualChanLen, MinChanLen ) ! set this to at least MinChanLen , or the size of the column produced by OutFmt
!          do i=1,numOuts
!             out%ActualChanLen = max(out%ActualChanLen, LEN_TRIM(out%WriteOutputHdr(i)))
!             out%ActualChanLen = max(out%ActualChanLen, LEN_TRIM(out%WriteOutputUnt(i)))
!          end do
! 
!          ! create format statements for time and the array outputs:
!          out%Fmt_t = '(F'//trim(num2lstr(out%ActualChanLen))//'.4)'
!          out%Fmt_a = '"'//out%delim//'"'//trim(out%outFmt)      ! format for array elements from individual modules
!          numSpaces = out%ActualChanLen - numSpaces  ! the difference between the size of the headers and what is produced by OutFmt
!          if (numSpaces > 0) then
!             out%Fmt_a = trim(out%Fmt_a)//','//trim(num2lstr(numSpaces))//'x'
!          end if
! 
!          ! --- Start writing to ascii input file 
!          do iWT=1,nWT
!             if (nWT>1) then
!                sWT = '.T'//trim(num2lstr(iWT))
!             else
!                sWT = ''
!             endif
!             call GetNewUnit(out%unOutFile(iWT), errStat, errMsg)
!             if ( errStat >= AbortErrLev ) then
!                out%unOutFile(iWT) = -1
!                return
!             end if
!             call OpenFOutFile ( out%unOutFile(iWT), trim(out%Root)//trim(sWT)//'.out', errStat, errMsg )
!             if ( errStat >= AbortErrLev ) return
!             write (out%unOutFile(iWT),'(/,A)')  'Predictions were generated on '//CurDate()//' at '//CurTime()//' using '//trim( version%Name )
!             write (out%unOutFile(iWT),'(1X,A)') trim(GetNVD(out%AD_ver))
!             write (out%unOutFile(iWT),'()' )    !print a blank line
!             write (out%unOutFile(iWT),'()' )    !print a blank line
!             write (out%unOutFile(iWT),'()' )    !print a blank line
! 
!             ! Write the names of the output parameters on one line:
!             do i=1,numOuts
!                call WrFileNR ( out%unOutFile(iWT), out%delim//out%WriteOutputHdr(i)(1:out%ActualChanLen) )
!             end do ! i
!             write (out%unOutFile(iWT),'()')
! 
!             ! Write the units of the output parameters on one line:
!             do i=1,numOuts
!                call WrFileNR ( out%unOutFile(iWT), out%delim//out%WriteOutputUnt(i)(1:out%ActualChanLen) )
!             end do ! i
!             write (out%unOutFile(iWT),'()')
!          enddo
!       endif
! 
      ! --- Binary
      if (out%fileFmt==idFmt_Both .or. out%fileFmt == idFmt_Binary) then
         call AllocAry(out%storage, numOuts-1, numSteps, 'storage', errStat, errMsg)
         out%storage= myNaN !0.0_ReKi ! Alternative: myNaN
      endif
end subroutine Dvr_InitializeOutputs
!----------------------------------------------------------------------------------------------------------------------------------
!> Initialize driver (not module-level) output channels 
subroutine Dvr_InitializeDriverOutputs(dvr, out, errStat, errMsg)
   type(Dvr_Data),        intent(inout) :: dvr              ! driver data
   type(Dvr_Outputs),     intent(inout) :: out              ! driver output data
   integer(IntKi)         ,  intent(  out) :: errStat              ! Status of error message
   character(*)           ,  intent(  out) :: errMsg               ! Error message if errStat /= ErrID_None
   integer(IntKi)       :: errStat2 ! Status of error message
   character(ErrMsgLen) :: errMsg2  ! Error message
   integer :: j
   errStat = ErrID_None
   errMsg  = ''

   out%ny_UA = size(dvr%UA_InitOutData%WriteOutputHdr)
   if (dvr%p%SimMod==3) then
      out%ny_dvr = 27 ! Driver only ! TODO
      out%ny_LD = size(dvr%LD_InitOutData%WriteOutputHdr)
   else
      out%ny_dvr = 0 
      out%ny_LD = 0
   endif


   ! --- Allocate driver-level outputs
   call AllocAry(out%WriteOutputHdr, 1+out%ny_dvr, 'WriteOutputHdr', errStat2, errMsg2); if(Failed()) return
   call AllocAry(out%WriteOutputUnt, 1+out%ny_dvr, 'WriteOutputUnt', errStat2, errMsg2); if(Failed()) return

   j=1
   out%WriteOutputHdr(j) = 'Time'        ; out%WriteOutputUnt(j) = '(s)'  ; j=j+1
   if (dvr%p%SimMod==3) then
   ! TODO SIMMOD HARMONIZATION
   ! Driver Variables
   out%WriteOutputHdr(j) = 'VUndx'           ; out%WriteOutputUnt(j) = '(m/s)'  ; j=j+1
   out%WriteOutputHdr(j) = 'VUndy'           ; out%WriteOutputUnt(j) = '(m/s)'  ; j=j+1
   out%WriteOutputHdr(j) = 'VSTx_Q'          ; out%WriteOutputUnt(j) = '(m/s)'  ; j=j+1
   out%WriteOutputHdr(j) = 'VSTy_Q'          ; out%WriteOutputUnt(j) = '(m/s)'  ; j=j+1
   out%WriteOutputHdr(j) = 'VSTx_T'          ; out%WriteOutputUnt(j) = '(m/s)'  ; j=j+1
   out%WriteOutputHdr(j) = 'VSTy_T'          ; out%WriteOutputUnt(j) = '(m/s)'  ; j=j+1
   out%WriteOutputHdr(j) = 'Vrelx_Q'         ; out%WriteOutputUnt(j) = '(m/s)'  ; j=j+1
   out%WriteOutputHdr(j) = 'Vrely_Q'         ; out%WriteOutputUnt(j) = '(m/s)'  ; j=j+1
   out%WriteOutputHdr(j) = 'Vrelx_T'         ; out%WriteOutputUnt(j) = '(m/s)'  ; j=j+1
   out%WriteOutputHdr(j) = 'Vrely_T'         ; out%WriteOutputUnt(j) = '(m/s)'  ; j=j+1
   out%WriteOutputHdr(j) = 'Vrel_Q'          ; out%WriteOutputUnt(j) = '(m/s)'  ; j=j+1
   out%WriteOutputHdr(j) = 'Vrel_T'          ; out%WriteOutputUnt(j) = '(m/s)'  ; j=j+1
   out%WriteOutputHdr(j) = 'alpha_Q'         ; out%WriteOutputUnt(j) = '(deg)'  ; j=j+1
   out%WriteOutputHdr(j) = 'alpha_T'         ; out%WriteOutputUnt(j) = '(deg)'  ; j=j+1
   out%WriteOutputHdr(j) = 'phi_Q'           ; out%WriteOutputUnt(j) = '(deg)'  ; j=j+1
   out%WriteOutputHdr(j) = 'phi_T'           ; out%WriteOutputUnt(j) = '(deg)'  ; j=j+1
   out%WriteOutputHdr(j) = 'twist_full'      ; out%WriteOutputUnt(j) = '(deg)'  ; j=j+1
   out%WriteOutputHdr(j) = 'Re_T'            ; out%WriteOutputUnt(j) = '(-)'  ; j=j+1
   out%WriteOutputHdr(j) = 'L'               ; out%WriteOutputUnt(j) = '(N/m)'  ; j=j+1
   out%WriteOutputHdr(j) = 'D'               ; out%WriteOutputUnt(j) = '(N/m)'  ; j=j+1
   out%WriteOutputHdr(j) = 'M'               ; out%WriteOutputUnt(j) = '(Nm/m)'  ; j=j+1
   out%WriteOutputHdr(j) = 'Fx_A'            ; out%WriteOutputUnt(j) = '(N/m)'  ; j=j+1
   out%WriteOutputHdr(j) = 'Fy_A'            ; out%WriteOutputUnt(j) = '(N/m)'  ; j=j+1
   out%WriteOutputHdr(j) = 'M_A'             ; out%WriteOutputUnt(j) = '(Nm/m)'  ; j=j+1
   out%WriteOutputHdr(j) = 'GFx'             ; out%WriteOutputUnt(j) = '(N/m)'  ; j=j+1
   out%WriteOutputHdr(j) = 'GFy'             ; out%WriteOutputUnt(j) = '(N/m)'  ; j=j+1
   out%WriteOutputHdr(j) = 'GFM'             ; out%WriteOutputUnt(j) = '(Nm/m)'  ; j=j+1
   ! Dynamics 
   call concatOutputHeaders(out%WriteOutputHdr, out%WriteOutputUnt, dvr%LD_InitOutData%WriteOutputHdr, dvr%LD_InitOutData%WriteOutputUnt, errStat2, errMsg2)
   endif
   ! UA
   call concatOutputHeaders(out%WriteOutputHdr, out%WriteOutputUnt, dvr%UA_InitOutData%WriteOutputHdr, dvr%UA_InitOutData%WriteOutputUnt, errStat2, errMsg2)

   out%ny = size(out%WriteOutputHdr)
   ! Debug Write
   !do j = 1, out%ny
   !   print*,'Write Out: ',j, trim(out%WriteOutputHdr(j)), ' ', trim(out%WriteOutputUnt(j))
   !enddo
contains
   logical function Failed()
      CALL SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Dvr_InitializeDriverOutputs' )
      Failed = errStat >= AbortErrLev
   end function Failed
end subroutine Dvr_InitializeDriverOutputs
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Dvr_WriteOutputs(nt, t, dvr, out, errStat, errMsg)
   integer(IntKi)         ,  intent(in   )   :: nt                   ! simulation time step
   real(DbKi)             ,  intent(in   )   :: t                    ! simulation time (s)
   type(Dvr_Data),           intent(inout)   :: dvr                  ! driver data
   type(Dvr_Outputs)      ,  intent(inout)   :: out                  ! driver uotput options
   integer(IntKi)         ,  intent(inout)   :: errStat              ! Status of error message
   character(*)           ,  intent(inout)   :: errMsg               ! Error message if errStat /= ErrID_None
!    ! Local variables.
!    character(ChanLen) :: tmpStr         ! temporary string to print the time output as text
   integer :: nDV , nUA, nLD,j
   errStat = ErrID_None
   errMsg  = ''
   out%outLine = myNaN ! Safety
! 
   ! --- Packing all outputs except time into one array
   nDV = out%ny_dvr
   nLD = out%ny_LD
   nUA = out%ny_UA
   ! Driver outputs
   j = 1 
   if (dvr%p%SimMod==3) then
      ! TODO harmonization
      out%outLine(j) = dvr%U0(1, 1)             ; j=j+1  ! Ux
      out%outLine(j) = dvr%U0(1, 2)             ; j=j+1  ! Uy
      out%outLine(j) = dvr%m%Vst_Q(1)           ; j=j+1  ! VSTx_Q
      out%outLine(j) = dvr%m%Vst_Q(2)           ; j=j+1  ! VSTy_Q
      out%outLine(j) = dvr%m%Vst_T(1)           ; j=j+1  ! VSTx_T
      out%outLine(j) = dvr%m%Vst_T(2)           ; j=j+1  ! VSTy_T
      out%outLine(j) = dvr%m%Vrel_Q(1)          ; j=j+1  ! Vrelx_Q
      out%outLine(j) = dvr%m%Vrel_Q(2)          ; j=j+1  ! Vrely_Q
      out%outLine(j) = dvr%m%Vrel_T(1)          ; j=j+1  ! Vrelx_T
      out%outLine(j) = dvr%m%Vrel_T(2)          ; j=j+1  ! Vrely_T
      out%outLine(j) = sqrt(dvr%m%Vrel_norm2_Q) ; j=j+1  ! Vrel_Q
      out%outLine(j) = sqrt(dvr%m%Vrel_norm2_T) ; j=j+1  ! Vrel_T
      out%outLine(j) = dvr%m%alpha_Q*R2D        ; j=j+1  ! alpha_Q
      out%outLine(j) = dvr%m%alpha_T*R2D        ; j=j+1  ! alpha_T
      out%outLine(j) = dvr%m%phi_Q  *R2D        ; j=j+1  ! phi_Q
      out%outLine(j) = dvr%m%phi_T  *R2D        ; j=j+1  ! phi_T
      out%outLine(j) = dvr%m%twist_full*R2D     ; j=j+1  ! twist_full
      out%outLine(j) = dvr%m%Re                 ; j=j+1  ! Re_T
      out%outLine(j) = dvr%m%L                  ; j=j+1  ! L
      out%outLine(j) = dvr%m%D                  ; j=j+1  ! D
      out%outLine(j) = dvr%m%tau_Q              ; j=j+1  ! M
      out%outLine(j) = dvr%m%FxA                ; j=j+1  ! Fx_A
      out%outLine(j) = dvr%m%FyA                ; j=j+1  ! Fy_A
      out%outLine(j) = dvr%m%tau_A              ; j=j+1  ! M_A
      out%outLine(j) = dvr%m%GF(1)              ; j=j+1  ! GFx
      out%outLine(j) = dvr%m%GF(2)              ; j=j+1  ! GFy
      out%outLine(j) = dvr%m%GF(3)              ; j=j+1  ! GFM
      ! LD Outputs
      out%outLine(nDV+1:nDV+nLD) = dvr%LD_y%WriteOutput(1:nLD)
   endif
   ! UA Outputs
   out%outLine(nDV+nLD+1:nDV+nLD+nUA) = dvr%UA_y%WriteOutput(1:nUA)

   !if (out%fileFmt==idFmt_Both .or. out%fileFmt == idFmt_Ascii) then
   !   ! ASCII
   !   ! time
   !   write( tmpStr, out%Fmt_t ) t  ! '(F15.4)'
   !   call WrFileNR( out%unOutFile, tmpStr(1:out%ActualChanLen) )
   !   call WrNumAryFileNR(out%unOutFile, out%outLine,  out%Fmt_a, errStat, errMsg)
   !   ! write a new line (advance to the next line)
   !   write(out%unOutFile,'()')
   !endif
   if (out%fileFmt==idFmt_Both .or. out%fileFmt == idFmt_Binary) then
      ! Store for binary
      out%storage(:, nt) = out%outLine(:)
      !out%storage(1:nDV+nAD+nIW, nt) = out%outLine(1:nDV+nAD+nIW)
   endif
end subroutine Dvr_WriteOutputs

end module UA_Dvr_Subs
   
