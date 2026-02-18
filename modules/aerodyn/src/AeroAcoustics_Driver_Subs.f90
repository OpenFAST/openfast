module AeroAcoustics_Driver_Subs
   
   use NWTC_Library
   use AirfoilInfo
   use AirfoilInfo_Types
   use AeroAcoustics
   use AeroAcoustics_Types
   
   implicit none

   integer, parameter        :: NumAFfiles = 1
   integer, parameter        :: NumBlades = 1
   integer, parameter        :: NumBlNds = 1
   logical, parameter        :: UseCm = .false.

   integer(IntKi), parameter :: idFmt_Ascii  = 1
   integer(IntKi), parameter :: idFmt_Binary = 2
   integer(IntKi), parameter :: idFmt_Both   = 3
   integer(IntKi), parameter :: idFmt_Valid(3)  = (/idFmt_Ascii, idFmt_Binary, idFmt_Both/)
   real(ReKi), parameter     :: myNaN = -9999.9
   character(1), parameter   :: delim = TAB

   real(DbKi), parameter :: RotGtoL(3,3) = reshape( [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], SHAPE=[3,3] )
   TYPE(ProgDesc), PARAMETER   :: version   = ProgDesc( 'AeroAcoustics_driver', '', '' )  ! The version number of this program.


   type Dvr_Data
      ! Environmental Conditions
      real(ReKi)                                     :: KinVisc                   !< Kinematic viscosity of working fluid (m^2/s)
      real(ReKi)                                     :: AirDens                   !< AirDens | Air density (kg/m^3)
      real(ReKi)                                     :: SpdSound                  !< Speed of sound in working fluid (m/s)

      ! Output data
      character(1024)                                :: OutRootName  = ''         !< output file rootname [-]
      integer(IntKi)                                 :: unOutFile = -1            !< unit number for writing text output file
      character(256)                                 :: OutFmt                    !< Format used for text tabular output, excluding the time channel.  Resulting field should be 10 characters. (quoted string)
      integer(IntKi)                                 :: OutFileFmt = idFmt_Binary !< Format for tabular (time-marching) output file (switch) {1: text file [<RootName>.out], 2: binary file [<RootName>.outb], 3: both}
      logical                                        :: WrBinaryOutput=.false.
      logical                                        :: WrTextOutput=.false.
      integer(IntKi)                                 :: NumOuts= 0                !< number of output channels, including time
      integer(IntKi)                                 :: NumSteps= 0               !< number of steps in output
      character(ChanLen) , dimension(:), allocatable :: WriteOutputHdr            !< channel headers [-]
      character(ChanLen) , dimension(:), allocatable :: WriteOutputUnt            !< channel units [-]
      real(ReKi) , dimension(:,:), allocatable       :: storage                   !< nchannel x ntime [-]
      real(ReKi) , dimension(:), allocatable         :: outline                   !< output line to be written to disk [-]
      integer(IntKi)                                 :: FmtWidth

      ! AeroAcoustics Input data
      REAL(DbKi)                                     :: AeroCent_G(3)             !< Global position of the blade node
      REAL(ReKi)                                     :: vRel                      !< Relative velocity (m/s)
      REAL(ReKi)                                     :: AoA                       !< Angle of attack (rad)
      REAL(ReKi)                                     :: WindSpeed                 !< Atmospheric undisturbed flow on blade [Inflow] (m/s)
      REAL(ReKi)                                     :: HubHeight                 !< hub height (m)
      REAL(ReKi)                                     :: BladeLength               !< effectively the element span (m) since we are running this with only one element

      ! Time control
      real(DbKi)                                     :: DT                        !< Simulation time step [used only when AnalysisType/=3] (s)
      real(DbKi)                                     :: TMax                      !< Total run time [used only when AnalysisType/=3] (s)

      ! AFI data
      character(1024)                                :: AirFoil_FileName
      real(ReKi)                                     :: Chord = 1.0
      type(AFI_ParameterType)                        :: AFInfo(NumAFfiles)
!      integer, allocatable                           :: AFIndx(:,:)

      ! AeroAcoustics data
      character(1024)                                :: AA_InputFileName    !< name of the AA input file
      type(AA_InitInputType)                         :: InitInp       !< Input data for initialization routine
      type(AA_InputType)                             :: u             !< An initial guess for the input; input mesh must be defined
      type(AA_ParameterType)                         :: p             !< Parameters
     !type(AA_ContinuousStateType)                   :: x             !< Initial continuous states
      type(AA_DiscreteStateType)                     :: xd            !< Initial discrete states
     !type(AA_ConstraintStateType)                   :: z             !< Initial guess of the constraint states
      type(AA_OtherStateType)                        :: OtherState    !< Initial other states
      type(AA_OutputType)                            :: y             !< Initial system outputs (outputs are not calculated)
      type(AA_MiscVarType)                           :: m             !< Initial misc/optimization variables
   end type Dvr_Data

contains
   
!--------------------------------------------------------------------------------------------------------------
subroutine ReadDriverInputFile( FileName, DriverData, ErrStat, ErrMsg )
   character(1024),               intent(in   )   :: FileName
   type(Dvr_Data),                intent(inout)   :: DriverData           ! driver data
   integer,                       intent(  out)   :: ErrStat              ! returns a non-zero value when an error occurs  
   character(*),                  intent(  out)   :: ErrMsg               ! Error message if ErrStat /= ErrID_None

   ! Local variables  
   integer                 :: UnEcho   ! The local unit number for this module's echo file
   integer                 :: iLine, i
   character(1024)         :: EchoFile ! Name of driver's echo file
   character(1024)         :: PriPath  ! the path to the primary input file
   character(1024)         :: Line     ! the path to the primary input file
   type(FileInfoType)      :: FI       !< The derived type for holding the file information.
   integer(IntKi)          :: errStat2 ! Status of error message
   character(1024)         :: errMsg2  ! Error message if ErrStat /= ErrID_None
   character(*), parameter :: RoutineName = 'ReadDriverInputFile'
   integer, parameter      :: NumHeaderLines = 3
   logical                 :: Echo
   
   ! Initialize the echo file unit to -1 which is the default to prevent echoing, we will alter this based on user input
   UnEcho  = -1
   Echo = .false.
   ErrStat = ErrID_None
   ErrMsg  = ''

   ! Read all input file lines into fileinfo
   call WrScr(' Opening AeroAcoustics Driver input file: '//trim(FileName) )
   call ProcessComFile(FileName, FI, errStat2, errMsg2); if (Failed()) return
   CALL GetPath( FileName, PriPath )    ! Input files will be relative to the path where the primary input file is located.
   !call GetRoot(FileName, dvr%root)      

   ! --- Header and echo
   iLine = NumHeaderLines ! Skip the first NumHeaderLines lines as they are known to be header lines and separators
   call ParseVar(FI, iLine, 'Echo', Echo, errStat2, errMsg2); if (Failed()) return;
   if ( Echo ) then
      EchoFile = trim(FileName)//'.ech'
      call OpenEcho (UnEcho, EchoFile, errStat2, errMsg2 ); if(Failed()) return
      do i = 1,iLine-1
         write(UnEcho, '(A)') trim(FI%Lines(i))
      enddo
   end if

   call ParseVar(FI, iLine, 'TMax',              DriverData%TMax ,            errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'DT',                DriverData%DT,               errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'AA_InputFile',      DriverData%AA_InputFileName ,    errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'AirFoil_FileName' , DriverData%AirFoil_FileName, errStat2, errMsg2, UnEcho); if(Failed()) return
   
   ! --- Environmental conditions section
   call ParseCom(FI, iLine, Line                           , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'AirDens',  DriverData%AirDens , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'KinVisc',  DriverData%KinVisc , errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'SpdSound', DriverData%SpdSound, errStat2, errMsg2, UnEcho); if(Failed()) return

   ! --- SIMULATION INPUTS section
   call ParseCom(FI, iLine, Line                                 ,    errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'WindSpeed' ,  DriverData%WindSpeed  ,    errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'AoA' ,        DriverData%AoA        ,    errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'vRel' ,       DriverData%vRel       ,    errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseAry(FI, iLine, 'AeroCent_G' , DriverData%AeroCent_G , 3, errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'HubHeight' ,  DriverData%HubHeight  ,    errStat2, errMsg2, UnEcho); if(Failed()) return
!  call ParseVar(FI, iLine, 'Chord'   ,    DriverData%Chord      ,    errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'Span' ,       DriverData%BladeLength,    errStat2, errMsg2, UnEcho); if(Failed()) return

   ! --- OUTPUT section
   call ParseCom(FI, iLine, Line                        ,        errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'OutFmt'   ,  DriverData%OutFmt   ,  errStat2, errMsg2, UnEcho); if(Failed()) return
   call ParseVar(FI, iLine, 'OutFileFmt', DriverData%OutFileFmt, errStat2, errMsg2, UnEcho); if(Failed()) return

   ! convert units:
   DriverData%AoA = DriverData%AoA * D2R
   
   
   ! --- Get relative path names
   call GetRoot(FileName, DriverData%OutRootName) ! OutRootName is inferred from current filename.
  !if (PathIsRelative(DriverData%OutRootName))      DriverData%OutRootName       = TRIM(PriPath)//TRIM(DriverData%OutRootName)
   if (PathIsRelative(DriverData%AA_InputFileName)) DriverData%AA_InputFileName  = TRIM(PriPath)//TRIM(DriverData%AA_InputFileName)
   if (PathIsRelative(DriverData%AirFoil_FileName)) DriverData%AirFoil_FileName  = TRIM(PriPath)//TRIM(DriverData%AirFoil_FileName  )

   ! --- Checks
   if (DriverData%OutFileFmt == idFmt_Both) then
      DriverData%WrBinaryOutput = .true.
      DriverData%WrTextOutput = .true.
   elseif (DriverData%OutFileFmt == idFmt_Ascii) then
      DriverData%WrBinaryOutput = .false.
      DriverData%WrTextOutput = .true.
   elseif (DriverData%OutFileFmt == idFmt_Binary) then
      DriverData%WrBinaryOutput = .true.
      DriverData%WrTextOutput = .false.
   else
      DriverData%WrBinaryOutput = .false.
      DriverData%WrTextOutput = .false.
   end if
   
   if (DriverData%WrTextOutput) then
      CALL ChkRealFmtStr( DriverData%OutFmt, 'OutFmt', DriverData%FmtWidth, ErrStat2, ErrMsg2 )
      IF ( DriverData%FmtWidth < 10 ) CALL SetErrStat( ErrID_Warn, 'OutFmt produces a column width of '// &
         TRIM(Num2LStr(DriverData%FmtWidth))//'), which may be too small.', ErrStat, ErrMsg, RoutineName )
   end if

   call Cleanup()
contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call Cleanup()
   end function Failed

   subroutine Cleanup()
      ! Close this module's echo file    
      if ( Echo ) then
         close(UnEcho)
      end if
      Call NWTC_Library_Destroyfileinfotype(FI, errStat2, errMsg2)
   end subroutine Cleanup


end subroutine ReadDriverInputFile
!--------------------------------------------------------------------------------------------------------------

subroutine Dvr_EndOutput(DriverData, nt, errStat, errMsg)
   type(Dvr_Data),          intent(inout) :: DriverData       ! driver data
   integer(IntKi),          intent(in   ) :: nt               ! number of time steps written
   integer(IntKi)         , intent(out)   :: errStat          ! Status of error message
   character(*)           , intent(out)   :: errMsg           ! Error message if errStat /= ErrID_None
   
   character(ErrMsgLen)       :: errMsg2                 ! temporary Error message if errStat /= ErrID_None
   integer(IntKi)             :: errStat2                ! temporary Error status of the operation
   character(*), parameter    :: RoutineName = 'Dvr_EndOutput'
   
   errStat = ErrID_None
   errMsg  = ''
   
   ! Close the output file
   if (DriverData%WrTextOutput) then
      if (DriverData%unOutFile > 0) close(DriverData%unOutFile)
      DriverData%unOutFile = -1
   endif
   if (DriverData%WrBinaryOutput .and. allocated(DriverData%storage)) then
      call WrScr(' Writing output file: '//trim(DriverData%OutRootName)//'.outb')
      call WrBinFAST(trim(DriverData%OutRootName)//'.outb', FileFmtID_ChanLen_In, version%Name, DriverData%WriteOutputHdr, DriverData%WriteOutputUnt, (/0.0_DbKi, DriverData%dt/), DriverData%storage(:,1:nt), errStat2, errMsg2)
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
   endif
end subroutine Dvr_EndOutput


!----------------------------------------------------------------------------------------------------------------------------------
!> Initialize outputs to file for driver 
subroutine Dvr_InitializeOutputs(DriverData, AA_InitOut, errStat, errMsg)
   TYPE(Dvr_Data) ,          intent(inout)   :: DriverData
   TYPE(AA_InitOutputTYpe),  intent(in   )   :: AA_InitOut
   integer(IntKi)         ,  intent(  out)   :: errStat              ! Status of error message
   character(*)           ,  intent(  out)   :: errMsg               ! Error message if errStat /= ErrID_None
   ! locals
   integer(IntKi)       :: errStat2 ! Status of error message
   character(ErrMsgLen) :: errMsg2  ! Error message
   integer              :: i, j
   integer              :: ActualChanLen

   errStat = ErrID_None
   errMsg  = ''

   DriverData%numSteps = ceiling(DriverData%TMax / DriverData%dt)
   DriverData%numOuts = sum(DriverData%p%numOutsAll) + 1 ! includes time channel
   if (DriverData%numOuts < 2) then
      ErrStat2=ErrID_Fatal
      ErrMsg2='AeroAcoustics module is not printing any outputs. Simulation will end.'
      if (Failed()) return
   end if

   ! --- Allocate driver-level outputs
   call AllocAry(DriverData%WriteOutputHdr, DriverData%numOuts, 'WriteOutputHdr', errStat2, errMsg2); if(Failed()) return
   call AllocAry(DriverData%WriteOutputUnt, DriverData%numOuts, 'WriteOutputUnt', errStat2, errMsg2); if(Failed()) return

   i=1
   DriverData%WriteOutputHdr(i) = 'Time'
   DriverData%WriteOutputUnt(i) = '(s)'
   
   if (DriverData%numOuts > 0) then
      do j=1,DriverData%p%numOutsAll(1)
         i = i + 1
         DriverData%WriteOutputHdr(i) = AA_InitOut%WriteOutputHdr(j)
         DriverData%WriteOutputUnt(i) = AA_InitOut%WriteOutputUnt(j)
      end do

      do j=1,DriverData%p%numOutsAll(2)
         i = i + 1
         DriverData%WriteOutputHdr(i) = AA_InitOut%WriteOutputHdrforPE(j)
         DriverData%WriteOutputUnt(i) = AA_InitOut%WriteOutputUntforPE(j)
      end do

      do j=1,DriverData%p%numOutsAll(3)
         i = i + 1
         DriverData%WriteOutputHdr(i) = AA_InitOut%WriteOutputHdrSep(j)
         DriverData%WriteOutputUnt(i) = AA_InitOut%WriteOutputUntSep(j)
      end do

      do j=1,DriverData%p%numOutsAll(4)
         i = i + 1
         DriverData%WriteOutputHdr(i) = AA_InitOut%WriteOutputHdrNodes(j)
         DriverData%WriteOutputUnt(i) = AA_InitOut%WriteOutputUntNodes(j)
      end do
   end if

   if (DriverData%WrTextOutput .or. DriverData%WrBinaryOutput) then
      call AllocAry(DriverData%outLine, DriverData%numOuts-1, 'outLine', errStat2, errMsg2); if(Failed()) return
      DriverData%outLine=0.0_ReKi
   end if

   if (DriverData%WrTextOutput) then
      ActualChanLen = min(ChanLen, max(10, DriverData%FmtWidth))
      do i=1,DriverData%NumOuts
         ActualChanLen = max( ActualChanLen, LEN_TRIM(DriverData%WriteOutputHdr(I)) )
         ActualChanLen = max( ActualChanLen, LEN_TRIM(DriverData%WriteOutputUnt(I)) )
      enddo ! I

      call GetNewUnit( DriverData%unOutFile, ErrStat2, ErrMsg2 )
      if (Failed()) return

      call OpenFOutFile ( DriverData%unOutFile, trim(DriverData%OutRootName)//'.out', ErrStat2, ErrMsg2 )
      if (Failed()) return

      write (DriverData%unOutFile,'(A)')  ''
      write (DriverData%unOutFile,'(A)')  'Predictions were generated on '//CurDate()//' at '//CurTime()//' using '//trim(GetNVD(version))
      write (DriverData%unOutFile,'(A)')  ''
      write (DriverData%unOutFile,'(A)')  ''
      write (DriverData%unOutFile,'(A)')  'Output from AeroAcoustics driver'
      write (DriverData%unOutFile,'(A)')  ''

      !......................................................
      ! Write the names of the output parameters on one line: line 7
      !......................................................
      call WrFileNR ( DriverData%unOutFile, DriverData%WriteOutputHdr(1)(1:min(15,ChanLen)) )
      do i=2,DriverData%NumOuts
         call WrFileNR ( DriverData%unOutFile, delim//DriverData%WriteOutputHdr(i)(1:ActualChanLen) )
      end do ! i
      write (DriverData%unOutFile,'()')
      
      !......................................................
      ! Write the units of the output parameters on one line: line 8
      !......................................................
      call WrFileNR ( DriverData%unOutFile, DriverData%WriteOutputUnt(1)(1:min(15,ChanLen)) )
      do i=2,DriverData%NumOuts
         call WrFileNR ( DriverData%unOutFile, delim//DriverData%WriteOutputUnt(i)(1:ActualChanLen) )
      end do ! i
      write (DriverData%unOutFile,'()')

   end if
   
   ! --- Binary
   if (DriverData%WrBinaryOutput) then
      ! we aren't storing time here
      call AllocAry(DriverData%storage, DriverData%numOuts-1, DriverData%numSteps, 'storage', errStat2, errMsg2); if(Failed()) return
      DriverData%storage= myNaN !0.0_ReKi ! Alternative: myNaN
   endif

contains
   logical function Failed()
      CALL SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Dvr_InitializeOutputs' )
      Failed = errStat >= AbortErrLev
   end function Failed
end subroutine Dvr_InitializeOutputs
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Dvr_WriteOutputs(t, nt, DriverData)
   Integer(IntKi)         ,  intent(in   )   :: nt                   ! time step number
   real(DbKi)             ,  intent(in   )   :: t                    ! simulation time (s)
   type(Dvr_Data),           intent(inout)   :: DriverData           ! driver data

   !    ! Local variables.
   integer :: i, j

   if (DriverData%WrTextOutput .or. DriverData%WrBinaryOutput) then
      i = 0

      ! Driver outputs
      if (DriverData%numOuts > 0) then
         do j=1,DriverData%p%numOutsAll(1)
            i = i + 1
            DriverData%outLine(i) = DriverData%y%WriteOutput(j)
         end do

         do j=1,DriverData%p%numOutsAll(2)
            i = i + 1
            DriverData%outLine(i) = DriverData%y%WriteOutputforPE(j)
         end do

         do j=1,DriverData%p%numOutsAll(3)
            i = i + 1
            DriverData%outLine(i) = DriverData%y%WriteOutputSep(j)
         end do

         do j=1,DriverData%p%numOutsAll(4)
            i = i + 1
            DriverData%outLine(i) = DriverData%y%WriteOutputNodes(j)
         end do
      end if
      
      if (DriverData%WrBinaryOutput) DriverData%storage(:,nt) = DriverData%outLine
      if (DriverData%WrTextOutput) then
         write(DriverData%unOutFile,'(F15.4,'//trim(num2lstr(DriverData%numOuts-1))//'("'//delim//'"'//trim(DriverData%outFmt)//'))') t, DriverData%outLine(:)
      end if

         
   end if
   
   
end subroutine Dvr_WriteOutputs
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Init_AFI(afName, AFInfo, ErrStat, ErrMsg)
   
   CHARACTER(1024),         intent(in   )  :: afName
   type(AFI_ParameterType), intent(  out)  :: AFInfo(NumAFfiles)
   integer(IntKi),          intent(  out)  :: ErrStat                       ! Error status.
   character(*),            intent(  out)  :: ErrMsg                        ! Error message.

   type(AFI_InitInputType)  :: AFI_InitInputs
   integer                  :: UnEc

   ErrStat = ErrID_None
   ErrMsg  = ""
   
   AFI_InitInputs%InCol_Alfa  = 1
   AFI_InitInputs%InCol_Cl    = 2
   AFI_InitInputs%InCol_Cd    = 3
   AFI_InitInputs%InCol_Cm    = 0
   AFI_InitInputs%InCol_Cpmin = 0
   AFI_InitInputs%AFTabMod    = AFITable_1 ! 1D-interpolation (on AoA only)
   AFI_InitInputs%UAMod       = 3  ! We calculate some of the UA coefficients based on UA Model, but AA doesn't care which
   AFI_InitInputs%FileName    = afName !InitInp%AF_File(i)
   
   UnEc = 0

   ! Read in and process the airfoil file.
   ! This includes creating the spline coefficients to be used for interpolation.

   call AFI_Init ( AFI_InitInputs, AFInfo(1), errStat, errMsg, UnEc )
   if (ErrStat >= AbortErrLev) return
   
end subroutine Init_AFI
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine initializes the Airfoil Noise module from within AeroDyn.
SUBROUTINE Init_AAmodule( DriverData, ErrStat, ErrMsg )
!..................................................................................................................................
   type(Dvr_Data),               intent(inout) :: DriverData    !< AeroDyn-level initialization inputs

   integer(IntKi),               intent(  out) :: errStat        !< Error status of the operation
   character(*),                 intent(  out) :: errMsg         !< Error message if ErrStat /= ErrID_None

   ! Local variables
   real(DbKi)                                  :: Interval       ! DT
   type(AA_InitInputType)                      :: InitInp        ! Input data for initialization routine
   type(AA_InitOutputType)                     :: InitOut        ! Output for initialization routine
   integer(intKi)                              :: j              ! node index
   integer(intKi)                              :: k              ! blade index
   integer(IntKi)                              :: ErrStat2
   character(ErrMsgLen)                        :: ErrMsg2
   character(*), parameter                     :: RoutineName = 'Init_AAmodule'
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   ! Transfer from parameters and input file to init input
   InitInp%InputFile        = DriverData%AA_InputFileName
   InitInp%NumBlades        = NumBlades
   InitInp%NumBlNds         = NumBlNds
   InitInp%RootName         = DriverData%OutRootName

! read from input file or set default value
   Interval                 = DriverData%DT
   InitInp%airDens          = DriverData%airDens  !(rho)
   InitInp%kinVisc          = DriverData%kinVisc  !(nu)
   InitInp%SpdSound         = DriverData%SpdSound !(co)
   InitInp%HubHeight        = DriverData%HubHeight

   ! --- Allocate and set AirfoilID, chord and Span for each blades
   ! note here that each blade is required to have the same number of nodes
   call AllocAry( InitInp%BlAFID,  NumBlNds, NumBlades,'InitInp%BlAFID', errStat2, ErrMsg2 ); call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   call AllocAry( InitInp%BlChord, NumBlNds, NumBlades, 'BlChord', errStat2, ErrMsg2 ); call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   call AllocAry( InitInp%BlSpn,   NumBlNds, NumBlades, 'BlSpn', errStat2, ErrMsg2 ); call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   if (ErrStat >= AbortErrLev) then
      call cleanup()
      return
   end if
   
   do k = 1, NumBlades
      do j=1, NumBlNds
         InitInp%BlChord(j,k)  = DriverData%Chord !RotInputFileData%BladeProps(k)%BlChord(j)
         InitInp%BlSpn  (j,k)  = real(j,ReKi)/real(NumBlNds,ReKi) * DriverData%BladeLength
         InitInp%BlAFID(j,k)   = NumAFfiles !RotInputFileData%BladeProps(k)%BlAFID(j)
      end do
   end do
   
   ! --- AeroAcoustics initialization call
   call AA_Init(InitInp, DriverData%u, DriverData%p, DriverData%xd, DriverData%OtherState,DriverData%y, DriverData%m, Interval, DriverData%AFInfo, InitOut, ErrStat2, ErrMsg2 )
   call SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)   

   if (ErrStat < AbortErrLev) then
      call Dvr_InitializeOutputs(DriverData, InitOut, errStat2, errMsg2)
         call SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
   end if
   
   call Cleanup()
   
contains   

   subroutine Cleanup()
      call AA_DestroyInitInput ( InitInp, ErrStat2, ErrMsg2 )   
      call AA_DestroyInitOutput ( InitOut, ErrStat2, ErrMsg2 )   
   end subroutine Cleanup
   
END SUBROUTINE Init_AAmodule
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine sets m%AA_u.
subroutine SetInputsForAA(DriverData)
   type(Dvr_Data),          intent(inout) :: DriverData    !< AeroDyn-level initialization inputs

   ! local variables
   integer(intKi)                         :: i        ! loop counter for nodes
   integer(intKi)                         :: j        ! loop counter for blades
   
   do j=1,NumBlades
      do i = 1,NumBlNds
         ! Get local orientation matrix to transform from blade element coordinates to global coordinates
         DriverData%u%RotGtoL(:,:,i,j) = RotGtoL ! default to identitiy orientation

         ! Get blade element aerodynamic center in global coordinates
         DriverData%u%AeroCent_G(:,i,j) = DriverData%AeroCent_G !BJJ: does this need to change with time? probably

         ! Set the blade element relative velocity (including induction)
         DriverData%u%Vrel(i,j) = DriverData%VRel
   
         ! Set the blade element angle of attack
         DriverData%u%AoANoise(i,j) = DriverData%AoA

         ! Set the blade element undisturbed flow
         DriverData%u%Inflow(:,i,j) = [DriverData%WindSpeed, 0.0_ReKi, 0.0_ReKi]
      end do
   end do
end subroutine SetInputsForAA
!----------------------------------------------------------------------------------------------------------------------------------

end module AeroAcoustics_Driver_Subs
   
