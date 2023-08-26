module UA_Dvr_Subs
   
   use NWTC_Library
   use AirfoilInfo
   use AirfoilInfo_Types
   use UnsteadyAero_Types
   use UnsteadyAero
   use LinDyn
   
   implicit none

   integer, parameter        :: NumAFfiles = 1
   integer(IntKi), parameter :: NumInp = 2           ! Number of inputs sent to UA_UpdateStates (must be at least 2)
   real(ReKi), parameter     :: myNaN = -99.9_ReKi
   integer(IntKi), parameter :: idFmtAscii  = 1
   integer(IntKi), parameter :: idFmtBinary = 2
   integer(IntKi), parameter :: idFmtBoth   = 3
   integer(IntKi), parameter, dimension(3) :: idFmtVALID  = (/idFmtAscii, idFmtBinary, idFmtBoth/)

   type UA_Dvr_InitInput
      logical         :: Echo            
      real(ReKi)      :: SpdSound        
      character(1024) :: OutRootName
      real(ReKi)      :: InflowVel       
      integer         :: UAMod           
      logical         :: Flookup        
      logical         :: UseCm         
      character(1024) :: AirFoil1 
      real(ReKi)      :: Chord
      integer         :: SimMod          
      real(ReKi)      :: NCycles         
      real(ReKi)      :: Frequency 
      real(ReKi)      :: Re
      integer         :: StepsPerCycle
      real(ReKi)      :: Amplitude       
      real(ReKi)      :: Mean            
      integer         :: Phase           
      character(1024) :: InputsFile      
      logical         :: SumPrint         
      logical         :: WrAFITables
      ! Section
      real(ReKi)      :: TMax
      real(DbKi)      :: dt
      real(ReKi)      :: MM(3,3)
      real(ReKi)      :: CC(3,3)
      real(ReKi)      :: KK(3,3)
      logical         :: activeDOFs(3)
      real(ReKi)      :: initPos(3)
      real(ReKi)      :: initVel(3)
   end type UA_Dvr_InitInput


   type :: Dvr_Outputs
      integer(intki)                                 :: unOutFile = -1        !< unit number for writing output file
      !integer(intki)                                :: actualchanlen         !< actual length of channels written to text file (less than or equal to chanlen) [-]
      integer(intki)                                :: nDvrOutputs=0           !< number of outputs for the driver (without ad and iw) [-]
      !character(20)                                 :: fmt_t                 !< format specifier for time channel [-]
      !character(25)                                 :: fmt_a                 !< format specifier for each column (including delimiter) [-]
      !character(1)                                  :: delim                 !< column delimiter [-]
      !character(20)                                 :: outfmt                !< format specifier [-]
      integer(intki)                                 :: fileFmt = idFmtBinary !< output format 1=text, 2=binary, 3=both [-]
      character(1024)                                :: root  = ''            !< output file rootname [-]
      character(chanlen) , dimension(:), allocatable :: writeoutputhdr        !< channel headers [-]
      character(chanlen) , dimension(:), allocatable :: writeoutputunt        !< channel units [-]
      real(ReKi) , dimension(:,:), allocatable       :: storage               !< nchannel x ntime [-]
      real(ReKi) , dimension(:), allocatable         :: outline               !< output line to be written to disk [-]
      !real(dbki)  :: dt_outs      !< output time resolution [s]
      !integer(intki)  :: n_dt_out      !< number of time steps between writing a line in the time-marching output files [-]
   end type Dvr_Outputs

   type Dvr_Data
      real(DbKi)                             :: dt
      type(Dvr_Outputs)                      :: out
      type(UA_InitInputType)      , pointer  :: UA_InitInData           ! Input data for initialization
      type(UA_InitOutputType)     , pointer  :: UA_InitOutData          ! Output data from initialization
      type(UA_ContinuousStateType), pointer  :: UA_x                    ! Continuous states
      type(UA_DiscreteStateType)  , pointer  :: UA_xd                   ! Discrete states
      type(UA_OtherStateType)     , pointer  :: UA_OtherState           ! Other/optimization states
      type(UA_MiscVarType)        , pointer  :: UA_m                    ! Misc/optimization variables
      type(UA_ParameterType)      , pointer  :: UA_p                    ! Parameters
      type(UA_InputType)          , pointer  :: UA_u(:)                 ! System inputs
      type(UA_OutputType)         , pointer  :: UA_y                    ! System outputs
      type(LD_InitInputType)      , pointer  :: LD_InitInData           ! Input data for initialization
      type(LD_InitOutputType)     , pointer  :: LD_InitOutData          ! Output data from initialization
      type(LD_ContinuousStateType), pointer  :: LD_x                    ! Continuous states
      type(LD_DiscreteStateType)  , pointer  :: LD_xd                   ! Discrete states
      type(LD_OtherStateType)     , pointer  :: LD_OtherState           ! Other/optimization states
      type(LD_ConstraintStateType), pointer  :: LD_z                    ! Constraint states
      type(LD_MiscVarType)        , pointer  :: LD_m                    ! Misc/optimization variables
      type(LD_ParameterType)      , pointer  :: LD_p                    ! Parameters
      type(LD_InputType)          , pointer  :: LD_u(:)                 ! System inputs
      type(LD_OutputType)         , pointer  :: LD_y                    ! System outputs
   end type Dvr_Data

contains
   
!--------------------------------------------------------------------------------------------------------------
subroutine ReadDriverInputFile( FileName, InitInp, ErrStat, ErrMsg )
   character(1024),               intent( in    )   :: filename
   type(UA_Dvr_InitInput),        intent(   out )   :: InitInp
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
   character(*), parameter :: RoutineName = 'ReadDriverfilename'
   ! Initialize the echo file unit to -1 which is the default to prevent echoing, we will alter this based on user input
   UnEcho  = -1
   ErrStat = ErrID_None
   ErrMsg  = ''

   ! Read all input file lines into fileinfo
   call WrScr( '  Opening UnsteadyAero Driver input file:  '//trim(FileName) )
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
   call ParseVar(FI, iLine, 'SpdSound', InitInp%SpdSound, errStat2, errMsg2, UnEcho); if(Failed()) return

   ! --- UNSTEADYAERO section
   call ParseCom(FI, iLine, Line                              , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'OutRootName', InitInp%OutRootName, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'InflowVel'  , InitInp%InflowVel  , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'Re'         , InitInp%Re         , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'UAMod'      , InitInp%UAMod      , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'Flookup'    , InitInp%Flookup    , errStat2, errMsg2, UnEcho); if(Failed()) return
   
   ! --- AIRFOIL PROPERTIES section
   call ParseCom(FI, iLine, Line                        , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'AirFoil' , InitInp%AirFoil1, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'Chord'   , InitInp%Chord   , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'UseCm'   , InitInp%UseCm   , errStat2, errMsg2, UnEcho); if(Failed()) return
   
   ! --- SIMULATION CONTROL section
   call ParseCom(FI, iLine, Line                                  , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'SimMod'       , InitInp%SimMod       , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'NCycles'      , InitInp%NCycles      , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'StepsPerCycle', InitInp%StepsPerCycle, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'Frequency'    , InitInp%Frequency    , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'Amplitude'    , InitInp%Amplitude    , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'Mean'         , InitInp%Mean         , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'Phase'        , InitInp%Phase        , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'InputsFile'   , InitInp%InputsFile   , errStat2, errMsg2, UnEcho); if(Failed()) return

   ! --- ELASTIC SECTION section
   if (InitInp%SimMod==3) then ! Temporary to avoid changing r-test for now
   call ParseCom(FI, iLine, Line                        , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'TMax'         , InitInp%Tmax         , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'DT'           , InitInp%dt           , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'activeDOFs'   , InitInp%activeDOFs, 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'initPos'      , InitInp%initPos   , 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'initVel'      , InitInp%initVel   , 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'MassMatrix1'  , InitInp%MM(1,:)   , 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'MassMatrix2'  , InitInp%MM(2,:)   , 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'MassMatrix3'  , InitInp%MM(3,:)   , 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'DampMatrix1'  , InitInp%CC(1,:)   , 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'DampMatrix2'  , InitInp%CC(2,:)   , 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'DampMatrix3'  , InitInp%CC(3,:)   , 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'StifMatrix1'  , InitInp%KK(1,:)   , 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'StifMatrix2'  , InitInp%KK(2,:)   , 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'StifMatrix3'  , InitInp%KK(3,:)   , 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   endif

   ! --- OUTPUT section
   call ParseCom(FI, iLine, Line                        , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'SumPrint'   , InitInp%SumPrint   , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'WrAFITables', InitInp%WrAFITables, errStat2, errMsg2, UnEcho); if(Failed()) return

   ! --- Triggers
   if (PathIsRelative(InitInp%OutRootName)) InitInp%OutRootName = TRIM(PriPath)//TRIM(InitInp%OutRootName)
   if (PathIsRelative(InitInp%Airfoil1)) InitInp%Airfoil1 = TRIM(PriPath)//TRIM(InitInp%Airfoil1)

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
      CALL NWTC_Library_Destroyfileinfotype(FI, errStat2, errMsg2)
   end subroutine Cleanup
end subroutine ReadDriverInputFile
!--------------------------------------------------------------------------------------------------------------
subroutine ReadTimeSeriesData( FileName, nSimSteps, timeArr, AOAarr, Uarr, OmegaArr, ErrStat, ErrMsg )
   character(1024),               intent( in    )   :: FileName
   integer,                       intent(   out )   :: nSimSteps
   real(DbKi),allocatable,        intent(   out )   :: timeArr(:)
   real(ReKi),allocatable,        intent(   out )   :: AOAarr(:)
   real(ReKi),allocatable,        intent(   out )   :: Uarr(:)
   real(ReKi),allocatable,        intent(   out )   :: OmegaArr(:)
   integer,                       intent(   out )   :: ErrStat              ! returns a non-zero value when an error occurs  
   character(*),                  intent(   out )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   
   real(SiKi)                                       :: dt
   real(DbKi)                                       :: tmpArr(4)
   integer(IntKi)                                   :: errStat2    ! Status of error message
   character(1024)                                  :: errMsg2     ! Error message if ErrStat /= ErrID_None
   character(*), parameter                          :: RoutineName = 'ReadTimeSeriesData'
   integer                                          :: UnIn
   integer                                          :: i
   integer, parameter                               :: hdrlines=8
   ErrStat     = ErrID_None
   ErrMsg      = ''
   nSimSteps   = 0 ! allocate here in case errors occur
   
   call WrScr( '  Opening UnsteadyAero time-series input file:  '//trim(FileName) )
   call GetNewUnit( UnIn )   
   call OpenFInpFile( UnIn, FileName, errStat2, errMsg2 ); if(Failed()) return

   ! --- Determine how many lines of data are in the file
   ! TODO use a more generic routine. For instane SubDyn has ReadDelimFile and line_count, which should be placed in NWTC_Lib
   do i=1,hdrlines
      call ReadCom( UnIn, FileName, '  UnsteadyAero time-series input file header line 1', errStat2, errMsg2 ); if(Failed()) return
   enddo
   do ! Loop on all lines..
      call  ReadAry( UnIn, FileName, tmpArr, 4, 'Data', 'Time-series data', errStat2, errMsg2  ); 
      ! The assumption is that the only parsing error occurs at the end of the file and therefore we stop reading data
      if (errStat2 > ErrID_None) then
         exit
      else
         nSimSteps = nSimSteps + 1
      end if        
   end do
   
   ! --- Allocate arrays to be read
   call AllocAry( timeArr , nSimSteps, 'timeArr' , errStat2, errMsg2); if(Failed()) return
   call AllocAry( AOAArr  , nSimSteps, 'AOAArr'  , errStat2, errMsg2); if(Failed()) return
   call AllocAry( Uarr    , nSimSteps, 'UArr'    , errStat2, errMsg2); if(Failed()) return
   call AllocAry( OmegaArr, nSimSteps, 'OmegaArr', errStat2, errMsg2); if(Failed()) return
      
   ! --- Read arrays from file
   rewind(UnIn)
   do i=1,hdrlines !RRD
       call ReadCom( UnIn, FileName, '  UnsteadyAero time-series input file header line 1', errStat2, errMsg2 ); if(Failed()) return
   enddo
   do i = 1,nSimSteps
      call  ReadAry( UnIn, FileName, tmpArr, 4, 'Data', 'Time-series data', errStat2, errMsg2  ); if(Failed()) return
      timeArr(i)  = tmpArr(1)
      AOAarr(i)   = real(tmpArr(2),ReKi)
      Uarr(i)     = real(tmpArr(3),ReKi)
      OmegaArr(i) = real(tmpArr(4),ReKi)
   end do

   ! --- Sanity checks 
   if (nSimSteps > 1) then
      ! TODO SubDyn allows for time interpolation of input array
      dt = timeArr(2) - timeArr(1)
      do i = 2,nSimSteps-1
         if (.not. EqualRealNos(dt, REAL(timeArr(i+1)-timeArr(i), SiKi) ) ) then
            call SetErrStat( ErrID_Fatal, 'Times in inputfile must be contain the same delta t.', ErrStat, ErrMsg, RoutineName)
            exit !exit the do loop
         end if
      end do
   end if
            
   call Cleanup()
      
contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call Cleanup()
   end function Failed
   subroutine Cleanup()
      close( UnIn )
   end subroutine Cleanup
end subroutine ReadTimeSeriesData
!--------------------------------------------------------------------------------------------------------------
subroutine driverInputsToUAInitData(dvrInitInp, InitInData, AFI_Params, AFIndx, errStat, errMsg)
   type(UA_Dvr_InitInput) , intent(in ) :: dvrInitInp           ! Initialization data for the driver program
   type(UA_InitInputType) , intent(out) :: InitInData           ! Input data for initialization
   type(AFI_ParameterType), intent(out) :: AFI_Params(NumAFfiles)
   integer, allocatable   , intent(out) :: AFIndx(:,:)
   integer(IntKi),          intent(out) :: errStat                       ! Error status.
   character(*),            intent(out) :: errMsg                        ! Error message.
   logical                 :: UA_f_cn ! Should the separation function be computed using Cn or Cl
   character(1024)         :: afNames(NumAFfiles)
   integer(IntKi)          :: errStat2    ! Status of error message
   character(1024)         :: errMsg2     ! Error message if ErrStat /= ErrID_None
   character(*), parameter  :: RoutineName = 'driverInputsToUAInitData'
   errStat     = ErrID_None
   errMsg      = ''

   ! -- UA Init Input Data
   InitInData%nNodesPerBlade  = 1 
   InitInData%numBlades       = 1
   call AllocAry(InitInData%c, InitInData%nNodesPerBlade, InitInData%numBlades, 'chord', errStat2, errMsg2); if(Failed()) return
   call AllocAry(InitInData%UAOff_innerNode             , InitInData%numBlades, 'UAO'  , errStat2, errMsg2); if(Failed()) return
   call AllocAry(InitInData%UAOff_outerNode             , InitInData%numBlades, 'UAO'  , errStat2, errMsg2); if(Failed()) return

   ! don't turn off UA based on span location:
   InitInData%UAOff_innerNode = 0
   InitInData%UAOff_outerNode = InitInData%nNodesPerBlade + 1
   InitInData%a_s          = dvrInitInp%SpdSound
   InitInData%c(1,1)       = dvrInitInp%Chord
   InitInData%UAMod        = dvrInitInp%UAMod 
   InitInData%Flookup      = dvrInitInp%Flookup
   InitInData%OutRootName  = dvrInitInp%OutRootName
   InitInData%WrSum        = dvrInitInp%SumPrint 

   ! --- AFI
   allocate(AFIndx(InitInData%nNodesPerBlade,InitInData%numBlades), STAT = errStat2)
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error trying to allocate InitInData%AFIndx.', errStat, errMsg, RoutineName)
      return
   end if
   AFIndx(1,1) = 1

   UA_f_cn  = (InitInData%UAMod /= UA_HGM).and.(InitInData%UAMod /= UA_OYE)  ! HGM and OYE use the separation function based on cl instead of cn

   afNames(1)  = dvrInitInp%AirFoil1 ! All nodes/blades are using the same 2D airfoil
   call Init_AFI( InitInData%UAMod, NumAFfiles, afNames, dvrInitInp%UseCm, UA_f_cn, AFI_Params, errStat2, errMsg2); if(Failed()) return

   if (dvrInitInp%WrAFITables) then
      call WriteAFITables(AFI_Params(1), dvrInitInp%OutRootName, dvrInitInp%UseCm, UA_f_cn)
   endif
contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
endsubroutine driverInputsToUAInitData
!--------------------------------------------------------------------------------------------------------------
subroutine Init_AFI(UAMod, NumAFfiles, afNames, UseCm, UA_f_cn, AFI_Params, ErrStat, ErrMsg)
   integer,             intent(in   )  :: UAMod
   integer,             intent(in   )  :: NumAFfiles
   CHARACTER(1024),     intent(in   )  :: afNames(NumAFfiles)
   logical,             intent(in   )  :: UseCm
   logical,             intent(in   )  :: UA_f_cn
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
   AFI_InitInputs%UA_f_cn     = UA_f_cn

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
   if (out%fileFmt==idFmtBoth .or. out%fileFmt == idFmtAscii) then
      if (out%unOutFile > 0) close(out%unOutFile)
   endif
   if (out%fileFmt==idFmtBoth .or. out%fileFmt == idFmtBinary) then
      print*,'>>>> OUTPUT',trim(out%Root)//'.outb'
      call WrBinFAST(trim(out%Root)//'.outb', FileFmtID_ChanLen_In, 'AeroDynDriver', out%WriteOutputHdr, out%WriteOutputUnt, (/0.0_DbKi, dvr%dt/), out%storage(:,:), errStat2, errMsg2)
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
   endif
end subroutine Dvr_EndSim

   
   
! --------------------------------------------------------------------------------
! --- IO 
! --------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
!> Concatenate new output channels info to the extisting ones in the driver
!! TODO COPY PASTED FROM AeroDyn_Inflow. Should be placed in NWTC_Lib
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
   !print*,'>>> Concat',allocated(WriteOutputHdr0), allocated(WriteOutputUnt0), allocated(WriteOutputHdr), allocated(WriteOutputUnt)
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
!       if (out%fileFmt==idFmtBoth .or. out%fileFmt == idFmtAscii) then
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
      if (out%fileFmt==idFmtBoth .or. out%fileFmt == idFmtBinary) then
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

   ! --- Allocate driver-level outputs
   out%nDvrOutputs =  6  ! temporary hack

   call AllocAry(out%WriteOutputHdr, 1+out%nDvrOutputs, 'WriteOutputHdr', errStat2, errMsg2); if(Failed()) return
   call AllocAry(out%WriteOutputUnt, 1+out%nDvrOutputs, 'WriteOutputUnt', errStat2, errMsg2); if(Failed()) return

   j=1
   out%WriteOutputHdr(j) = 'Time'        ; out%WriteOutputUnt(j) = '(s)'  ; j=j+1
   ! HACK
   out%WriteOutputHdr(j) = 'x'           ; out%WriteOutputUnt(j) = '(m)'  ; j=j+1
   out%WriteOutputHdr(j) = 'y'           ; out%WriteOutputUnt(j) = '(m)'  ; j=j+1
   out%WriteOutputHdr(j) = 'th'          ; out%WriteOutputUnt(j) = '(rad)'  ; j=j+1
   out%WriteOutputHdr(j) = 'dx'          ; out%WriteOutputUnt(j) = '(m/s)'  ; j=j+1
   out%WriteOutputHdr(j) = 'dy'          ; out%WriteOutputUnt(j) = '(m/s)'  ; j=j+1
   out%WriteOutputHdr(j) = 'dth'         ; out%WriteOutputUnt(j) = '(rad/s)'  ; j=j+1
contains
   logical function Failed()
      CALL SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Dvr_InitializeDriverOutputs' )
      Failed = errStat >= AbortErrLev
   end function Failed
end subroutine Dvr_InitializeDriverOutputs
!----------------------------------------------------------------------------------------------------------------------------------
! !> Store driver data
! subroutine Dvr_CalcOutputDriver(dvr, y_ADI, FED, errStat, errMsg)
!    type(Dvr_SimData), target,  intent(inout) :: dvr              ! driver data
!    type(FED_Data),    target,   intent(in   ) :: FED       !< Elastic wind turbine data (Fake ElastoDyn)
!    type(ADI_OutputType),        intent(in   ) :: y_ADI           ! ADI output data
!    integer(IntKi)           ,   intent(  out) :: errStat         ! Status of error message
!    character(*)             ,   intent(  out) :: errMsg          ! Error message if errStat /= ErrID_None
!    integer              :: maxNumBlades, k, j, iWT
!    real(ReKi)           :: rotations(3)
!    integer(IntKi)       :: errStat2        ! Status of error message
!    character(ErrMsgLen) :: errMsg2 ! Error message
!    real(ReKi), pointer  :: arr(:)
!    type(WTData), pointer :: wt ! Alias to shorten notation
!    type(RotFED), pointer :: y_ED ! Alias to shorten notation
! 
!    errStat = ErrID_None
!    errMsg  = ''
!    
!    maxNumBlades = 0
!    do iWT=1,size(dvr%WT)
!       maxNumBlades= max(maxNumBlades, dvr%WT(iWT)%numBlades)
!    end do
! 
!    ! Determine if a swap array is present
!    
!    do iWT = 1, dvr%numTurbines
!       wt => dvr%wt(iWT)
!       y_ED => FED%wt(iWT)
!       if (dvr%wt(iWT)%numBlades >0 ) then ! TODO, export for tower only
!          arr => dvr%wt(iWT)%WriteOutput
!          k=1
!          ! NOTE: to do this properly we would need to store at the previous time step and perform a rotation
!          arr(k) = dvr%iCase           ; k=k+1
!          ! Environment
!          arr(k) = y_ADI%HHVel(1, iWT) ; k=k+1  ! NOTE: stored at beginning of array
!          arr(k) = y_ADI%HHVel(2, iWT) ; k=k+1
!          arr(k) = y_ADI%HHVel(3, iWT) ; k=k+1 
!          arr(k) = y_ADI%PLExp         ; k=k+1 ! shear exp, not set if CompInflow=1
! 
!          ! 6 base DOF
!          rotations  = EulerExtract(y_ED%PlatformPtMesh%Orientation(:,:,1)); 
!          arr(k) = y_ED%PlatformPtMesh%TranslationDisp(1,1); k=k+1 ! surge
!          arr(k) = y_ED%PlatformPtMesh%TranslationDisp(2,1); k=k+1 ! sway
!          arr(k) = y_ED%PlatformPtMesh%TranslationDisp(3,1); k=k+1 ! heave
!          arr(k) = rotations(1) * R2D  ; k=k+1 ! roll
!          arr(k) = rotations(2) * R2D  ; k=k+1 ! pitch
!          arr(k) = rotations(3) * R2D  ; k=k+1 ! yaw
!          ! RNA motion
!          arr(k) = wt%nac%yaw*R2D         ; k=k+1 ! yaw [deg]
!          arr(k) = modulo(real(wt%hub%azimuth+(dvr%dt * wt%hub%rotSpeed)*R2D, ReKi), 360.0_ReKi); k=k+1 ! azimuth [deg], stored at nt-1
!          arr(k) = wt%hub%rotSpeed*RPS2RPM; k=k+1 ! rotspeed [rpm]
!          do j=1,maxNumBlades
!             if (j<= wt%numBlades) then
!                arr(k) = wt%bld(j)%pitch*R2D ! pitch [deg]
!             else
!                arr(k) = 0.0_ReKi ! myNaN
!             endif
!             k=k+1;
!          enddo
!          ! Swap array
!          if (wt%hub%motionType == idHubMotionUserFunction) then
!             do j=1,size(wt%userSwapArray)
!                arr(k) = wt%userSwapArray(j); k=k+1;
!             enddo
!          endif
! 
!       endif
!    enddo
! 
! end subroutine Dvr_CalcOutputDriver
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
   integer :: nDV , nUA, nLD
   errStat = ErrID_None
   errMsg  = ''
   out%outLine = myNaN ! Safety
! 
!    ! Packing all outputs excpet time into one array
   !nUA = size(yADI%AD%rotors(1)%WriteOutput)
   !nLD = size(yADI%IW_WriteOutput)
   nLD = 6 ! HACK
   nDV = out%nDvrOutputs
   !out%outLine(1:nDV)         = dvr%LD_x%q(1:nDV)  ! Driver Write Outputs
   out%outLine(1:nLD)         = dvr%LD_x%q(1:nDV)  ! Driver Write Outputs

 
   !out%outLine(nDV+1:nDV+nAD) = yADI%AD%rotors%WriteOutput     ! AeroDyn WriteOutputs
   !out%outLine(nDV+nAD+1:)    = yADI%IW_WriteOutput                 ! InflowWind WriteOutputs
   !if (out%fileFmt==idFmtBoth .or. out%fileFmt == idFmtAscii) then
   !   ! ASCII
   !   ! time
   !   write( tmpStr, out%Fmt_t ) t  ! '(F15.4)'
   !   call WrFileNR( out%unOutFile, tmpStr(1:out%ActualChanLen) )
   !   call WrNumAryFileNR(out%unOutFile, out%outLine,  out%Fmt_a, errStat, errMsg)
   !   ! write a new line (advance to the next line)
   !   write(out%unOutFile,'()')
   !endif
   if (out%fileFmt==idFmtBoth .or. out%fileFmt == idFmtBinary) then
      ! Store for binary
      out%storage(:, nt) = out%outLine(:)
      !out%storage(1:nDV+nAD+nIW, nt) = out%outLine(1:nDV+nAD+nIW)
   endif
end subroutine Dvr_WriteOutputs
! 


subroutine WriteAFITables(AFI_Params, OutRootName, UseCm, UA_f_cn)

   type(AFI_ParameterType), intent(in), target  :: AFI_Params
   character(len=*)       , intent(in)          :: OutRootName
   logical                , intent(in)          :: UseCm
   logical                , intent(in)          :: UA_f_cn
   
   integer(IntKi)                               :: unOutFile
   integer(IntKi)                               :: ErrStat
   character(ErrMsgLen)                         :: ErrMsg
   
   Real(ReKi), allocatable  :: cl_smooth(:)
   Real(ReKi), allocatable  :: cn_smooth(:)
   Real(ReKi), allocatable  :: cn(:)
   Real(ReKi), allocatable  :: cl_lin(:)
   Real(ReKi), allocatable  :: cn_lin(:)
   character(len=3) :: Prefix
   character(len=11) :: sFullyAtt
   character(len=8) :: sCm
   integer :: iTab, iRow, iStartUA
   type(AFI_Table_Type), pointer :: tab !< Alias

   if (UA_f_cn) then
      Prefix='Cn_'
      sFullyAtt='Cn_FullyAtt'
   else
      Prefix='Cl_'
      sFullyAtt='Dummy'
   endif
   if (UseCm) then
      sCm='Cm'
   else
      sCm='Cm_Dummy'
   endif


   ! Loop on tables, write a different file for each table.
   do iTab = 1, size(AFI_Params%Table)
      tab => AFI_Params%Table(iTab)

      ! Compute derived parameters from cl and cd, and UA_BL
      if(allocated(cl_smooth)) deallocate(cl_smooth)
      if(allocated(cn_smooth)) deallocate(cn_smooth)
      if(allocated(cn       )) deallocate(cn       )
      if(allocated(cl_lin   )) deallocate(cl_lin   )
      if(allocated(cn_lin   )) deallocate(cn_lin   )
      allocate(cl_smooth(tab%NumAlf))
      allocate(cn_smooth(tab%NumAlf))
      allocate(cn       (tab%NumAlf))
      allocate(cl_lin   (tab%NumAlf))
      allocate(cn_lin   (tab%NumAlf))

   
      cn     = tab%Coefs(:,AFI_Params%ColCl) * cos(tab%alpha) + (tab%Coefs(:,AFI_Params%ColCd) - tab%UA_BL%Cd0) * sin(tab%alpha);
      cn_lin = tab%UA_BL%C_nalpha * (tab%alpha - tab%UA_BL%alpha0)
      cl_lin = tab%UA_BL%C_lalpha * (tab%alpha - tab%UA_BL%alpha0)

      do iRow = 1, tab%NumAlf
         if ((tab%alpha(iRow)<tab%UA_BL%alphaLowerWrap).or. tab%alpha(iRow)>tab%UA_BL%alphaUpperWrap) then
            cl_lin(iRow) =0.0_ReKi
            cn_lin(iRow) =0.0_ReKi
         endif
      enddo

      ! Smoothing (used priot to compute slope in CalculateUACoeffs)
      call kernelSmoothing(tab%alpha, cn                           , kernelType_TRIWEIGHT, 2.0_ReKi*D2R, cn_smooth)
      call kernelSmoothing(tab%alpha, tab%Coefs(:,AFI_Params%ColCl), kernelType_TRIWEIGHT, 2.0_ReKi*D2R, cl_smooth)

      ! Write to file

      CALL GetNewUnit( unOutFile, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN

      CALL OpenFOutFile ( unOutFile, trim(OutRootName)//'.UA.Coefs.'//trim(num2lstr(iTab))//'.out', ErrStat, ErrMsg )
         if (ErrStat >= AbortErrLev) then
            call WrScr(Trim(ErrMsg))
            return
         end if
   
      WRITE (unOutFile,'(/,A/)') 'These predictions were generated by UnsteadyAero Driver on '//CurDate()//' at '//CurTime()//'.'
      WRITE (unOutFile,'(/,A/)')  ' '

      WRITE(unOutFile, '(20(A20,1x))') 'Alpha', 'Cl',  'Cd',  sCm,  'Cn', 'f_st', Prefix//'FullySep', sFullyAtt , 'Cl_lin','Cn_lin','Cl_smooth', 'Cn_smooth'
      WRITE(unOutFile, '(20(A20,1x))') '(deg)', '(-)', '(-)', '(-)', '(-)', '(-)', '(-)'             , '(-)'     ,  '(-)'  , '(-)'  , '(-)'    ,'(-)'

      ! TODO, we could do something with ColCpmim and ColUAf
      if (UseCm) then
         iStartUA = 4
         do iRow=1,size(tab%Alpha)
            WRITE(unOutFile, '(20(F20.6,1x))') tab%Alpha(iRow)*R2D, tab%Coefs(iRow,AFI_Params%ColCl), tab%Coefs(iRow,AFI_Params%ColCd), tab%Coefs(iRow,AFI_Params%ColCm), &
                                           cn(iRow),  tab%Coefs(iRow,iStartUA:), cl_lin(iRow), cn_lin(iRow), cl_smooth(iRow), cn_smooth(iRow)
         end do
      else
         iStartUA = 3
         do iRow=1,size(tab%Alpha)
            WRITE(unOutFile, '(20(F20.6,1x))') tab%Alpha(iRow)*R2D, tab%Coefs(iRow,AFI_Params%ColCl), tab%Coefs(iRow,AFI_Params%ColCd), 0.0_ReKi, &
                                           cn(iRow), tab%Coefs(iRow,iStartUA:), cl_lin(iRow), cn_lin(iRow), cl_smooth(iRow), cn_smooth(iRow)
         end do
      endif
      
      CLOSE(unOutFile)
   enddo
   
end subroutine WriteAFITables
   
end module UA_Dvr_Subs
   
