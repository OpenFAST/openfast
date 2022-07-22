module UA_Dvr_Subs
   
   use NWTC_Library
   use AirfoilInfo
   use AirfoilInfo_Types
   use UnsteadyAero_Types
   use UnsteadyAero
   
   implicit none

   
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
   end type UA_Dvr_InitInput
   
   contains
   
   subroutine ReadDriverInputFile( inputFile, InitInp, ErrStat, ErrMsg )

      character(1024),               intent( in    )   :: inputFile
      type(UA_Dvr_InitInput),       intent(   out )   :: InitInp
      integer,                       intent(   out )   :: ErrStat              ! returns a non-zero value when an error occurs  
      character(*),                  intent(   out )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   
         ! Local variables  
      integer                                          :: UnIn                 ! Unit number for the input file
      integer                                          :: UnEchoLocal          ! The local unit number for this module's echo file
      character(1024)                                  :: EchoFile             ! Name of HydroDyn echo file  
      character(1024)                                  :: FileName             ! Name of HydroDyn input file  

      integer(IntKi)                                   :: errStat2    ! Status of error message
      character(1024)                                  :: errMsg2     ! Error message if ErrStat /= ErrID_None
      character(*), parameter                          :: RoutineName = 'ReadDriverInputFile'

      character(1024) :: PriPath                         ! the path to the primary input file
      CALL GetPath( inputFile, PriPath )    ! Input files will be relative to the path where the primary input file is located.

   
         ! Initialize the echo file unit to -1 which is the default to prevent echoing, we will alter this based on user input
      UnEchoLocal = -1
      ErrStat     = ErrID_None
      ErrMsg      = ''
      FileName = trim(inputFile)
   
      call GetNewUnit( UnIn )   
      call OpenFInpFile( UnIn, FileName, errStat2, errMsg2 )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if

   
      call WrScr( '  Opening UnsteadyAero Driver input file:  '//FileName )
   
   
      !-------------------------------------------------------------------------------------------------
      ! File header
      !-------------------------------------------------------------------------------------------------
   
      call ReadCom( UnIn, FileName, '  UnsteadyAero Driver input file header line 1', errStat2, errMsg2 )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if


      call ReadCom( UnIn, FileName, 'UnsteadyAero Driver input file header line 2', errStat2, errMsg2 )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if

   
        ! Echo Input Files.
      
      call ReadVar ( UnIn, FileName, InitInp%Echo, 'Echo', 'Echo Input', errStat2, errMsg2 )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if
   
   
      ! If we are Echoing the input then we should re-read the first three lines so that we can echo them
      ! using the NWTC_Library routines.  The echoing is done inside those routines via a global variable
      ! which we must store, set, and then replace on error or completion.
      
      if ( InitInp%Echo ) then
      
         EchoFile = TRIM(FileName)//'.ech'
         call GetNewUnit( UnEchoLocal )   
         call OpenEcho ( UnEchoLocal, EchoFile, errStat2, errMsg2 )
            call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
            if (ErrStat >= AbortErrLev) then
               call Cleanup()
               return
            end if
      
         rewind(UnIn)
      
         call ReadCom( UnIn, FileName, 'UnsteadyAero Driver input file header line 1', errStat2, errMsg2, UnEchoLocal )
            call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
            if (ErrStat >= AbortErrLev) then
               call Cleanup()
               return
            end if


         call ReadCom( UnIn, FileName, 'UnsteadyAero Driver input file header line 2', errStat2, errMsg2, UnEchoLocal )
            call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
            if (ErrStat >= AbortErrLev) then
               call Cleanup()
               return
            end if

   
            ! Echo Input Files. Note this line is prevented from being echoed by the ReadVar routine.
      
         call ReadVar ( UnIn, FileName, InitInp%Echo, 'Echo', 'Echo the input file data', errStat2, errMsg2, UnEchoLocal )
            call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
            if (ErrStat >= AbortErrLev) then
               call Cleanup()
               return
            end if
      
      end if
      
   !-------------------------------------------------------------------------------------------------
   ! Environmental conditions section
   !-------------------------------------------------------------------------------------------------

      ! Header
      
      call ReadCom( UnIn, FileName, 'Environmental conditions header', errStat2, errMsg2, UnEchoLocal )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if


         ! SpdSound - Speed of Sound.
      
      call ReadVar ( UnIn, FileName, InitInp%SpdSound, 'SpdSound', 'SpdSound', errStat2, errMsg2, UnEchoLocal )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if

   
      !-------------------------------------------------------------------------------------------------
      ! UNSTEADYAERO section
      !-------------------------------------------------------------------------------------------------

         ! Header
      
      call ReadCom( UnIn, FileName, 'UNSTEADYAERO header', errStat2, errMsg2, UnEchoLocal )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if

      
         ! OutRootName  
      call ReadVar ( UnIn, FileName, InitInp%OutRootName, 'OutRootName', &
                                       'UnsteadyAero output root filename', errStat2, errMsg2, UnEchoLocal )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if
      if (PathIsRelative(InitInp%OutRootName)) InitInp%OutRootName = TRIM(PriPath)//TRIM(InitInp%OutRootName)
     
         ! InflowVel
   
      call ReadVar ( UnIn, FileName, InitInp%InflowVel, 'InflowVel', &
                                       'Inflow velocity', errStat2, errMsg2, UnEchoLocal )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if
   
         ! Re
   
      call ReadVar ( UnIn, FileName, InitInp%Re, 'Re', &
                                          'Reynolds number in millions', errStat2, errMsg2, UnEchoLocal )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if
         
         ! UAMod 
      call ReadVar ( UnIn, FileName, InitInp%UAMod, 'UAMod', &
                                       'Unsteady Aero Model Switch', errStat2, errMsg2, UnEchoLocal )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if
 
         ! Flookup 
      call ReadVar ( UnIn, FileName, InitInp%Flookup, 'Flookup', &
                                       "Lookup used to determine f'", errStat2, errMsg2, UnEchoLocal )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if
   
      
   
   
      !-------------------------------------------------------------------------------------------------
      ! AIRFOIL PROPERTIES section
      !-------------------------------------------------------------------------------------------------

         ! Header
      
      call ReadCom( UnIn, FileName, 'AIRFOIL PROPERTIES header', errStat2, errMsg2, UnEchoLocal )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if
 
   
   
         ! AirFoil1      
       
      call ReadVar ( UnIn, FileName, InitInp%AirFoil1, 'AirFoil1', &
                                       'Filename for the airfoil table and properties', errStat2, errMsg2, UnEchoLocal )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if
      if (PathIsRelative(InitInp%Airfoil1)) InitInp%Airfoil1 = TRIM(PriPath)//TRIM(InitInp%Airfoil1)
   
        ! Chord      
       
      call ReadVar ( UnIn, FileName, InitInp%Chord, 'Chord', &
                                       'Chord length', errStat2, errMsg2, UnEchoLocal )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if
      
         ! Using Cm column
      call ReadVar ( UnIn, FileName, InitInp%UseCm, 'UseCm', &
                                       "Using Cm Airfoil table data", errStat2, errMsg2, UnEchoLocal )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if
      
      !-------------------------------------------------------------------------------------------------
      ! SIMULATION CONTROL section
      !-------------------------------------------------------------------------------------------------

         ! Header
      
      call ReadCom( UnIn, FileName, 'SIMULATION CONTROL header', errStat2, errMsg2, UnEchoLocal )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if
   
      
         ! SimMod             
      call ReadVar ( UnIn, FileName, InitInp%SimMod, 'SimMod', &
                                       'Simulation model', errStat2, errMsg2, UnEchoLocal )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if
      
         ! NCycles             
      call ReadVar ( UnIn, FileName, InitInp%NCycles, 'NCycles', &
                                       'Number of cycles for angle-of-attack inputs', errStat2, errMsg2, UnEchoLocal )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if
      
         ! StepsPerCycle             
      call ReadVar ( UnIn, FileName, InitInp%StepsPerCycle, 'StepsPerCycle', &
                                       'Number of timesteps per cycle', errStat2, errMsg2, UnEchoLocal )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if
         
         ! Frequency             
      call ReadVar ( UnIn, FileName, InitInp%Frequency, 'Frequency', &
                                       'Frequency of angle-of-attack inputs', errStat2, errMsg2, UnEchoLocal )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if
      
         ! Amplitude             
      call ReadVar ( UnIn, FileName, InitInp%Amplitude, 'Amplitude', &
                                       'Amplitude for angle-of-attack inputs', errStat2, errMsg2, UnEchoLocal )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if
      
         ! Mean             
      call ReadVar ( UnIn, FileName, InitInp%Mean, 'Mean', &
                                       'Mean for angle-of-attack inputs', errStat2, errMsg2, UnEchoLocal )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if
      
         ! Phase             
      call ReadVar ( UnIn, FileName, InitInp%Phase, 'Phase', &
                                       'Initial phase for angle-of-attack inputs', errStat2, errMsg2, UnEchoLocal )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if
      
         ! InputsFile             
      call ReadVar ( UnIn, FileName, InitInp%InputsFile, 'InputsFile', &
                                       'Filename for Time series data in an ASCII input file', errStat2, errMsg2, UnEchoLocal )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if
         
      call Cleanup()
      
   
   contains
      !====================================================================================================
      subroutine Cleanup()
      !     The routine cleans up the module echo file and resets the NWTC_Library, reattaching it to 
      !     any existing echo information
      !----------------------------------------------------------------------------------------------------  
      !   logical,                       intent( in    )   :: EchoFlag             ! local version of echo flag
      !   integer,                       intent( in    )   :: UnEcho               !  echo unit number
   
            ! Close this module's echo file    
         if ( InitInp%Echo ) then
            close(UnEchoLocal)
         end if
         
         close( UnIn )
         
      end subroutine Cleanup
      

   end subroutine ReadDriverInputFile
   
   subroutine ReadTimeSeriesData( inputsFile, nSimSteps, timeArr, AOAarr, Uarr, OmegaArr, ErrStat, ErrMsg )
      character(1024),               intent( in    )   :: inputsFile
      integer,                       intent(   out )   :: nSimSteps
      real(DbKi),allocatable,        intent(   out )   :: timeArr(:)
      real(ReKi),allocatable,        intent(   out )   :: AOAarr(:)
      real(ReKi),allocatable,        intent(   out )   :: Uarr(:) !RRD
      real(ReKi),allocatable,        intent(   out )   :: OmegaArr(:)
      integer,                       intent(   out )   :: ErrStat              ! returns a non-zero value when an error occurs  
      character(*),                  intent(   out )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None
      
      real(SiKi)                                       :: dt
      real(DbKi)                                       :: tmpArr(4)
      integer(IntKi)                                   :: errStat2    ! Status of error message
      character(1024)                                  :: errMsg2     ! Error message if ErrStat /= ErrID_None
      character(*), parameter                          :: RoutineName = 'ReadTimeSeriesData'
      character(1024)                                  :: FileName
      integer                                          :: UnIn
      integer                                          :: i
      integer, PARAMETER                               ::hdrlines=8 ! RRD
      
      ErrStat     = ErrID_None
      ErrMsg      = ''
      nSimSteps   = 0 ! allocate here in case errors occur
      
      FileName = trim(inputsFile)
      
      call GetNewUnit( UnIn )   
      call OpenFInpFile( UnIn, FileName, errStat2, errMsg2 )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if

   
      call WrScr( '  Opening UnsteadyAero time-series input file:  '//FileName )
   
   
      !-------------------------------------------------------------------------------------------------
      ! Determine how many lines of data are in the file:
      !-------------------------------------------------------------------------------------------------
      do i=1,hdrlines !RRD
        call ReadCom( UnIn, FileName, '  UnsteadyAero time-series input file header line 1', errStat2, errMsg2 )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if
      enddo
   
      do
         call  ReadAry( UnIn, FileName, tmpArr, 4, 'Data', 'Time-series data', errStat2, errMsg2  )
            ! The assumption is that the only parsing error occurs at the end of the file and therefore we stop reading data
         if (errStat2 > ErrID_None) then
            exit
         else
            nSimSteps = nSimSteps + 1
         end if        
      end do
      
      !-------------------------------------------------------------------------------------------------
      ! Allocate arrays to be read
      !-------------------------------------------------------------------------------------------------
      allocate ( timeArr( nSimSteps ), STAT=ErrStat2 )
         if ( ErrStat2 /= 0 ) then
            call SetErrStat( ErrID_Fatal, 'Error trying to allocate timeArr.', ErrStat, ErrMsg, RoutineName)  
            call Cleanup()
            return
         end if
         
      allocate ( AOAarr( nSimSteps ), STAT=ErrStat2 )
         if ( ErrStat2 /= 0 ) then
            call SetErrStat( ErrID_Fatal, 'Error trying to allocate AOAarr.', ErrStat, ErrMsg, RoutineName)  
            call Cleanup()
            return
         end if

      allocate ( Uarr( nSimSteps ), OmegaArr( nSimSteps ), STAT=ErrStat2 )
         if ( ErrStat2 /= 0 ) then
            call SetErrStat( ErrID_Fatal, 'Error trying to allocate Uarr and OmegaArr.', ErrStat, ErrMsg, RoutineName)  
            call Cleanup()
            return
         end if
         
         
      !-------------------------------------------------------------------------------------------------
      ! Read arrays from file
      !-------------------------------------------------------------------------------------------------
      rewind(UnIn)
      do i=1,hdrlines !RRD
          call ReadCom( UnIn, FileName, '  UnsteadyAero time-series input file header line 1', errStat2, errMsg2 )
      enddo
      do i = 1,nSimSteps
         call  ReadAry( UnIn, FileName, tmpArr, 4, 'Data', 'Time-series data', errStat2, errMsg2  )
         timeArr(i)  = tmpArr(1)
         AOAarr(i)   = real(tmpArr(2),ReKi)
         Uarr(i)     = real(tmpArr(3),ReKi)
         OmegaArr(i) = real(tmpArr(4),ReKi)
      end do
      
      if (nSimSteps > 1) then
         dt = timeArr(2) - timeArr(1)
         
         do i = 2,nSimSteps-1
            if (.not. EqualRealNos(dt, REAL(timeArr(i+1)-timeArr(i), SiKi) ) ) then
               call SetErrStat( ErrID_Fatal, 'Times in InputsFile must be contain the same delta t.', ErrStat, ErrMsg, RoutineName)
               exit !exit the do loop
            end if
         end do
      end if
               
      call Cleanup()
         
      contains
      !====================================================================================================
      subroutine Cleanup()
      !     The routine cleans up the module echo file and resets the NWTC_Library, reattaching it to 
      !     any existing echo information
      !----------------------------------------------------------------------------------------------------  
      !   logical,                       intent( in    )   :: EchoFlag             ! local version of echo flag
      !   integer,                       intent( in    )   :: UnEcho               !  echo unit number
   
            ! Close this module's echo file    
         
         
         close( UnIn )
         
      end subroutine Cleanup
   end subroutine ReadTimeSeriesData
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
   
  ! Initialize the Airfoil Info module
      ! Setup Airfoil info
   
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
   AFI_InitInputs%AFTabMod = AFITable_1 ! 1D-interpolation (on AoA only)
   AFI_InitInputs%UA_f_cn  = UAMod /= UA_HGM ! HGM needs the separation function based on cl instead of cn
   
   do i=1,NumAFfiles
      AFI_InitInputs%FileName = afNames(i) !InitInp%AF_File(i)
      
         ! Call AFI_Init to read in and process the airfoil files.
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
   
      !====================================================================================================
      subroutine Cleanup()
      !     The routine cleans up data arrays framework structures 
      !      
      !----------------------------------------------------------------------------------------------------  
            !Clean up initialization inputs
         call AFI_DestroyInitInput(AFI_InitInputs, errStat2, errMsg2)
      
      
      
      end subroutine Cleanup
   end subroutine Init_AFI  
   
   
   subroutine WriteAFITables(AFI_Params,OutRootName)
   
      type(AFI_ParameterType), intent(in)          :: AFI_Params
      character(ErrMsgLen)   , intent(in)          :: OutRootName
      
      integer(IntKi)                               :: unOutFile
      integer(IntKi)                               :: row
      integer(IntKi)                               :: ErrStat
      character(ErrMsgLen)                         :: ErrMsg
      
      Real(ReKi)                                   :: cl_smooth(AFI_Params%Table(1)%NumAlf)
      Real(ReKi)                                   :: cn_smooth(AFI_Params%Table(1)%NumAlf)
      Real(ReKi)                                   :: cn(AFI_Params%Table(1)%NumAlf)
      
      cn = AFI_Params%Table(1)%Coefs(:,AFI_Params%ColCl) * cos(AFI_Params%Table(1)%alpha) + (AFI_Params%Table(1)%Coefs(:,AFI_Params%ColCd) - AFI_Params%Table(1)%UA_BL%Cd0) * sin(AFI_Params%Table(1)%alpha);

      call kernelSmoothing(AFI_Params%Table(1)%alpha, cn, kernelType_TRIWEIGHT, 2.0_ReKi*D2R, cn_smooth)
      call kernelSmoothing(AFI_Params%Table(1)%alpha, AFI_Params%Table(1)%Coefs(:,AFI_Params%ColCl), kernelType_TRIWEIGHT, 2.0_ReKi*D2R, cl_smooth)

      CALL GetNewUnit( unOutFile, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN

      CALL OpenFOutFile ( unOutFile, trim(OutRootName)//'.Coefs.out', ErrStat, ErrMsg )
         if (ErrStat >= AbortErrLev) then
            call WrScr(Trim(ErrMsg))
            return
         end if

   
      WRITE (unOutFile,'(/,A/)') 'These predictions were generated by UnsteadyAero Driver on '//CurDate()//' at '//CurTime()//'.'
      WRITE (unOutFile,'(/,A/)')  ' '
         ! note that this header assumes we have Cm and unsteady aero coefficients
      WRITE(unOutFile, '(20(A20,1x))') 'Alpha', 'Cl',  'Cd',  'Cm', 'f_st', 'FullySeparate', 'FullyAttached', 'smoothed_Cl', 'smoothed_Cn'
      WRITE(unOutFile, '(20(A20,1x))') '(deg)', '(-)', '(-)', '(-)', '(-)', '(-)',   '(-)', '(-)', '(-)'

      do row=1,size(AFI_Params%Table(1)%Alpha)
         WRITE(unOutFile, '(20(F20.6,1x))') AFI_Params%Table(1)%Alpha(row)*R2D, AFI_Params%Table(1)%Coefs(row,:), cl_smooth(Row), cn_smooth(Row)
      end do
      
      CLOSE(unOutFile)
      
   end subroutine WriteAFITables
   
end module UA_Dvr_Subs
   
