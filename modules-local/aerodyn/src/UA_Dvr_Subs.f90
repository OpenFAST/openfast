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
         
      integer                                          :: i                    ! generic integer for counting
      integer                                          :: j                    ! generic integer for counting
      character(   2)                                  :: strI                 ! string version of the loop counter

      integer                                          :: UnIn                 ! Unit number for the input file
      integer                                          :: UnEchoLocal          ! The local unit number for this module's echo file
      character(1024)                                  :: EchoFile             ! Name of HydroDyn echo file  
      character(1024)                                  :: Line                 ! String to temporarially hold value of read line   
      character(1024)                                  :: TmpPath              ! Temporary storage for relative path name
      character(1024)                                  :: TmpFmt               ! Temporary storage for format statement
      character(1024)                                  :: FileName             ! Name of HydroDyn input file  

      real(ReKi)                                       :: TmpRealVar2(2)       !< Temporary real    array size 2
      integer(IntKi)                                   :: TmpIntVar2(2)        !< Temporary integer array size 2
      integer(IntKi)                                   :: errStat2    ! Status of error message
      character(1024)                                  :: errMsg2     ! Error message if ErrStat /= ErrID_None
      character(1024)                                  :: RoutineName
   
         ! Initialize the echo file unit to -1 which is the default to prevent echoing, we will alter this based on user input
      UnEchoLocal = -1
      ErrStat     = ErrID_None
      ErrMsg      = ''
      RoutineName = 'ReadDriverInputFile'
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
   
   subroutine ReadTimeSeriesData( inputsFile, nSimSteps, timeArr, AOAarr, Uarr, ErrStat, ErrMsg )
      character(1024),               intent( in    )   :: inputsFile
      integer,                       intent(   out )   :: nSimSteps
      real(DbKi),allocatable,        intent(   out )   :: timeArr(:)
      real(ReKi),allocatable,        intent(   out )   :: AOAarr(:)
      real(ReKi),allocatable,        intent(   out )   :: Uarr(:) !RRD
      integer,                       intent(   out )   :: ErrStat              ! returns a non-zero value when an error occurs  
      character(*),                  intent(   out )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None
      
      real(DbKi)                                       :: tmpArr(3) !RRD changed from (2)
      integer(IntKi)                                   :: errStat2    ! Status of error message
      character(1024)                                  :: errMsg2     ! Error message if ErrStat /= ErrID_None
      character(1024)                                  :: RoutineName
      character(1024)                                  :: FileName
      integer                                          :: UnIn
      integer                                           :: i
      integer, PARAMETER                               ::hdrlines=8 ! RRD
      
      ErrStat     = ErrID_None
      ErrMsg      = ''
      RoutineName = 'ReadTimeSeriesData'
      FileName = trim(inputsFile)
      nSimSteps   = 0
      
      call GetNewUnit( UnIn )   
      call OpenFInpFile( UnIn, FileName, errStat2, errMsg2 )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if

   
      call WrScr( '  Opening UnsteadyAero time-series input file:  '//FileName )
   
   
      !-------------------------------------------------------------------------------------------------
      ! File header  7lines
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
         call  ReadAry( UnIn, FileName, tmpArr, 3, 'Data', 'Time-series data', errStat2, errMsg2  ) !RRD 2-->3
            ! The assumption is that the only parsing error occurs at the end of the file and therefor we stop reading data
         if (errStat2 > ErrID_None) then
            exit
         else
            nSimSteps = nSimSteps + 1
         end if        
      end do
      
      rewind(UnIn)
      do i=1,hdrlines !RRD
          call ReadCom( UnIn, FileName, '  UnsteadyAero time-series input file header line 1', errStat2, errMsg2 )
      enddo
      allocate ( timeArr( nSimSteps ), STAT=ErrStat )
         if ( ErrStat /= 0 ) then
            call SetErrStat( ErrID_Fatal, 'Error trying to allocate timeArr.', ErrStat, ErrMsg, RoutineName)  
            call Cleanup()
            stop       
         end if
         
      allocate ( AOAarr( nSimSteps ), STAT=ErrStat )
         if ( ErrStat /= 0 ) then
            call SetErrStat( ErrID_Fatal, 'Error trying to allocate AOAarr.', ErrStat, ErrMsg, RoutineName)  
            call Cleanup()
            stop       
         end if

      allocate ( Uarr( nSimSteps ), STAT=ErrStat ) !RRD
         if ( ErrStat /= 0 ) then
            call SetErrStat( ErrID_Fatal, 'Error trying to allocate Uarr.', ErrStat, ErrMsg, RoutineName)  
            call Cleanup()
            stop       
         end if
         
      do i = 1,nSimSteps
         call  ReadAry( UnIn, FileName, tmpArr, 3, 'Data', 'Time-series data', errStat2, errMsg2  ) !RRD 2-->3
         timeArr(i) = tmpArr(1)
         AOAarr(i)  = real(tmpArr(2))
         Uarr(i)  = real(tmpArr(3))
      end do
      
         
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
   
   subroutine Init_AFI(NumAFfiles, afNames, Flookup, UseCm, AFI_Params, ErrStat, ErrMsg)

   
   
   integer,             intent(in   )  :: NumAFfiles
   CHARACTER(1024),     intent(in   )  :: afNames(NumAFfiles)
   logical,             intent(in   )  :: Flookup
   logical,             intent(in   )  :: UseCm
   type(AFI_ParameterType), intent(  out)  :: AFI_Params
   integer(IntKi),      intent(  out)  :: ErrStat                       ! Error status.
   character(*),        intent(  out)  :: ErrMsg                        ! Error message.

   
   type(AFI_InitInputType)  :: AFI_InitInputs
   integer                  :: UnEc
   integer                  :: i
   integer(IntKi)           :: errStat2    ! Status of error message
   character(1024)          :: errMsg2     ! Error message if ErrStat /= ErrID_None
   character(1024)          :: RoutineName
   
  ! Initialize the Airfoil Info module
      ! Setup Airfoil info
   AFI_InitInputs%NumAFfiles = NumAFfiles
   UnEc = 0
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   RoutineName = 'Init_AFI'
   
   
   allocate ( AFI_InitInputs%FileNames( AFI_InitInputs%NumAFfiles ), STAT=ErrStat )
      if ( ErrStat /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'Error trying to allocate AFI_InitInputs%FileNames.', ErrStat, ErrMsg, RoutineName)  
         call Cleanup()
         stop       
      end if
   
   
   do i=1,AFI_InitInputs%NumAFfiles
      AFI_InitInputs%FileNames(i) = afNames(i) !InitInp%AF_File(i)
   end do
   
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
   
   !AFI_InitInputs%Flookup     = Flookup

   
   
      ! Call AFI_Init to read in and process the airfoil files.
      ! This includes creating the spline coefficients to be used for interpolation.

   call AFI_Init ( AFI_InitInputs, AFI_Params, errStat2, errMsg2, UnEc )
      call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if
 
   
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
   
   
   
   
   
   
   
   
   
   
   
   
   
   
end module UA_Dvr_Subs
   