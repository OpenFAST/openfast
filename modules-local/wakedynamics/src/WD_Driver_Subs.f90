module WD_Driver_Subs
   
   use WakeDynamics_Types
   
   implicit none
   

   
   contains
 
!!----------------------------------------------------------------------------------------------------------------------------------
!subroutine ReadWD_DvrInputFile( filename, OutFileRoot, InputFileData, , errStat, errMsg )
!! This subroutine reads the input file and stores all the data in the WD_InputFile structure.
!! It does not perform data validation.
!!..................................................................................................................................
!
!      ! Passed variables
!
!   character(*),            intent(in)    :: filename       ! Name of the input file
!   
!
!   type(WD_InputFile),      intent(out)   :: InputFileData   ! Data stored in the module's input file
!   real(DbKi),              intent(out)   :: Default_DT      ! The default DT (from glue code)
!
!   integer(IntKi),          intent(out)   :: errStat         ! The error status code
!   character(*),            intent(out)   :: errMsg          ! The error message, if an error occurred
!
!      ! local variables
!
!   integer(IntKi)                         :: I
!   integer(IntKi)                         :: UnEc
!   integer(IntKi)                         :: UnIn            ! Unit number for reading file
!   logical                                :: Echo            ! Determines if an echo file should be written
!   character(1024)                        :: PriPath         ! Path name of the primary file
!   character(1024)                        :: FTitle          ! "File Title": the 2nd line of the input file, which contains a description of its contents
!   character(200)                         :: Line            ! Temporary storage of a line from the input file (to compare with "default")
!   integer(IntKi)                         :: errStat2, IOS   ! The error status code
!   character(ErrMsgLen)                   :: errMsg2         ! The error message, if an error occurred
!   character(*)                           :: DvrFileRoot     ! The rootname of all the driver input file.
!   !character(*)                           :: OutFileRoot     ! The rootname of output files written by this routine.
!   character(*), parameter                :: RoutineName = 'ReadWD_DvrInputFile'
!   
!   
!      ! initialize values:
!
!   errStat = ErrID_None
!   errMsg  = ''
!   
!      ! get the driver input-file data
!   
!     ! Initialize some variables:
!   ErrStat = ErrID_None
!   ErrMsg  = ""
!      
!   UnEc = -1
!   Echo = .FALSE.   
!   call GetPath( filename, PriPath )     ! Input files will be relative to the path where the driver input file is located.
!   
!        
!      ! Get an available unit number for the file.
!
!   call GetNewUnit( UnIn, errStat2, errMsg2 )
!      call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!
!      ! Open the driver input file.
!
!   call OpenFInpFile ( UnIn, filename, errStat2, errMsg2 )
!      call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!      if ( ErrStat >= AbortErrLev ) then
!         call Cleanup()
!         return
!      end if
!      
!                  
!      
!   ! Read the lines up/including to the "Echo" simulation control variable
!   ! If echo is FALSE, don't write these lines to the echo file. 
!   ! If Echo is TRUE, rewind and write on the second try.
!   
!   I = 1 !set the number of times we've read the file
!   do 
!   !----------- HEADER -------------------------------------------------------------
!   
!      call ReadCom( UnIn, filename, 'File header: Module Version (line 1)', errStat2, errMsg2, UnEc )
!         call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!   
!      call ReadStr( UnIn, filename, FTitle, 'FTitle', 'File Header: File Description (line 2)', errStat2, errMsg2, UnEc )
!         call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!         if ( ErrStat >= AbortErrLev ) then
!            call Cleanup()
!            return
!         end if
!   
!   
!   !----------- GENERAL OPTIONS ----------------------------------------------------
!   
!      call ReadCom( UnIn, filename, 'Section Header: General Options', errStat2, errMsg2, UnEc )
!         call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!   
!         ! Echo - Echo input to "<RootName>.WD.ech".
!   
!      call ReadVar( UnIn, filename, Echo, 'Echo',   'Echo flag', errStat2, errMsg2, UnEc )
!         call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!   
!   
!      if (.NOT. Echo .OR. I > 1) exit !exit this loop
!   
!         ! Otherwise, open the echo file, then rewind the input file and echo everything we've read
!      call getroot(filename,DvrFileRoot)   
!      I = I + 1         ! make sure we do this only once (increment counter that says how many times we've read this file)
!   
!      call OpenEcho ( UnEc, trim(DvrFileRoot)//'.ech', errStat2, errMsg2, WD_Ver )
!         call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!         if ( ErrStat >= AbortErrLev ) then
!            call Cleanup()
!            return
!         end if
!   
!      if ( UnEc > 0 )  write (UnEc,'(/,A,/)')  'Data from '//trim(WD_Drv_Ver%Name)//' primary input file "'//trim( filename )//'":'
!         
!      rewind( UnIn, IOSTAT=errStat2 )  
!         if (errStat2 /= 0_IntKi ) then
!            call SetErrStat( ErrID_Fatal, 'Error rewinding file "'//trim(filename)//'".', ErrStat, ErrMsg, RoutineName )
!            call Cleanup()
!            return
!         end if         
!      
!   end do    
!
!   if (NWTC_VerboseLevel == NWTC_Verbose) then
!      call WrScr( ' Heading of the '//trim(WD_Drv_Ver%Name)//' input file: ' )      
!      call WrScr( '   '//trim( FTitle ) )
!   end if
!   
!   
!      ! DT - Time interval for Wake dynamics calculations {or default} (s):
!   Line = ""
!   call ReadVar( UnIn, filename, Line, "DT", "Time interval for Wake dynamics calculations {or default} (s)", errStat2, errMsg2, UnEc)
!      call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!         
!      call Conv2UC( Line )
!      if ( index(Line, "DEFAULT" ) /= 1 ) then ! If it's not "default", read this variable; otherwise use the value already stored in InputFileData%DT
!         READ( Line, *, IOSTAT=IOS) InputFileData%DT
!            call CheckIOS ( IOS, filename, 'DT', NumType, errStat2, errMsg2 )
!            call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!      end if   
!      
!      ! NumPlanes - Number of wake planes (-):
!   call ReadVar( UnIn, filename, InputFileData%NumPlanes, "NumPlanes", "Number of wake planes (-)", errStat2, errMsg2, UnEc)
!      call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!
!      ! NumRadii - Number of radii in the radial finite-difference grid (-):
!   call ReadVar( UnIn, filename, InputFileData%NumRadii, "NumRadii", "Number of radii in the radial finite-difference grid (-)", errStat2, errMsg2, UnEc)
!      call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!
!      ! dR - Radial increment of radial finite-difference grid (m) :
!   call ReadVar( UnIn, filename, InputFileData%dR, "dR", "Radial increment of radial finite-difference grid (m)", errStat2, errMsg2, UnEc)
!      call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!
!      ! filtParam - Low-pass time-filter parameter, with a value between 0 (minimum filtering) and 1 (maximum filtering) (exclusive) (-) :
!   call ReadVar( UnIn, filename, InputFileData%filtParam, "filtParam", "Low-pass time-filter parameter, with a value between 0 (minimum filtering) and 1 (maximum filtering) (exclusive) (-)", errStat2, errMsg2, UnEc)
!      call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!      
!      ! C_NearWake - Calibrated parameter for near-wake correction (-):
!   call ReadVar( UnIn, filename, InputFileData%C_NearWake, "C_NearWake", "Calibrated parameter for near-wake correction (-)", errStat2, errMsg2, UnEc)
!      call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!      
!      ! C_DMin_vAmb - Calibrated parameter in the eddy viscosity filter function for ambient turbulence defining the transitional diameter fraction between the minimum and exponential regions (-):
!   call ReadVar( UnIn, filename, InputFileData%C_DMin_vAmb, "C_DMin_vAmb", "Calibrated parameter in the eddy viscosity filter function for ambient turbulence defining the transitional diameter fraction between the minimum and exponential regions (-)", errStat2, errMsg2, UnEc)
!      call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!      
!      ! C_DMax_vAmb - Calibrated parameter in the eddy viscosity filter function for ambient turbulence defining the transitional diameter fraction between the exponential and maximum regions (-):
!   call ReadVar( UnIn, filename, InputFileData%C_DMax_vAmb, "C_DMax_vAmb", "Calibrated parameter in the eddy viscosity filter function for ambient turbulence defining the transitional diameter fraction between the exponential and maximum regions (-)", errStat2, errMsg2, UnEc)
!      call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!   
!      ! C_FMin_vAmb - Calibrated parameter in the eddy viscosity filter function for ambient turbulence defining the functional value in the minimum region (-):
!   call ReadVar( UnIn, filename, InputFileData%C_FMin_vAmb, "C_FMin_vAmb", "Calibrated parameter in the eddy viscosity filter function for ambient turbulence defining the functional value in the minimum region (-)", errStat2, errMsg2, UnEc)
!      call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!      
!      ! C_Exp_vAmb - Calibrated parameter in the eddy viscosity filter function for ambient turbulence defining the exponent in the exponential region (-):
!   call ReadVar( UnIn, filename, InputFileData%C_Exp_vAmb, "C_Exp_vAmb", "Calibrated parameter in the eddy viscosity filter function for ambient turbulence defining the exponent in the exponential region (-)", errStat2, errMsg2, UnEc)
!      call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!
!      ! C_DMin_vShr - Calibrated parameter in the eddy viscosity filter function for the shear layer defining the transitional diameter fraction between the minimum and exponential regions (-):
!   call ReadVar( UnIn, filename, InputFileData%C_DMin_vShr, "C_DMin_vShr", "Calibrated parameter in the eddy viscosity filter function for the shear layer defining the transitional diameter fraction between the minimum and exponential regions (-)", errStat2, errMsg2, UnEc)
!      call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!
!      ! C_DMax_vShr - Calibrated parameter in the eddy viscosity filter function for the shear layer defining the transitional diameter fraction between the exponential and maximum regions (-):
!   call ReadVar( UnIn, filename, InputFileData%C_DMax_vShr, "C_DMax_vShr", "Calibrated parameter in the eddy viscosity filter function for the shear layer defining the transitional diameter fraction between the exponential and maximum regions (-)", errStat2, errMsg2, UnEc)
!      call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!
!      ! C_FMin_vShr - Calibrated parameter in the eddy viscosity filter function for the shear layer defining the functional value in the minimum region (-):
!   call ReadVar( UnIn, filename, InputFileData%C_FMin_vShr, "C_FMin_vShr", "Calibrated parameter in the eddy viscosity filter function for the shear layer defining the functional value in the minimum region (-)", errStat2, errMsg2, UnEc)
!      call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!
!      ! C_Exp_vShr - Calibrated parameter in the eddy viscosity filter function for the shear layer defining the exponent in the exponential region (-):
!   call ReadVar( UnIn, filename, InputFileData%C_Exp_vShr, "C_Exp_vShr", "Calibrated parameter in the eddy viscosity filter function for the shear layer defining the exponent in the exponential region (-)", errStat2, errMsg2, UnEc)
!      call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!
!      ! k_vAmb - Calibrated parameter for the influence of ambient turbulence in the eddy viscosity (-):
!   call ReadVar( UnIn, filename, InputFileData%k_vAmb, "k_vAmb", "Calibrated parameter for the influence of ambient turbulence in the eddy viscosity (-)", errStat2, errMsg2, UnEc)
!      call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!
!      ! k_vShr - Calibrated parameter for the influence of the shear layer in the eddy viscosity (-):
!   call ReadVar( UnIn, filename, InputFileData%k_vShr, "k_vShr", "Calibrated parameter for the influence of the shear layer in the eddy viscosity (-)", errStat2, errMsg2, UnEc)
!      call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!
!      ! WakeDiam_Mod - Wake diameter calculation model (-):
!   call ReadVar( UnIn, filename, InputFileData%WakeDiam_Mod, "WakeDiam_Mod", "Wake diameter calculation model (-)", errStat2, errMsg2, UnEc)
!      call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!
!      ! C_WakeDiam - Calibrated parameter for wake diameter calculation (-):
!   call ReadVar( UnIn, filename, InputFileData%C_WakeDiam, "C_WakeDiam", "Calibrated parameter for wake diameter calculation (-)", errStat2, errMsg2, UnEc)
!      call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!
!
!
!      ! Return on error at end of section
!   if ( ErrStat >= AbortErrLev ) then
!      call Cleanup()
!      return
!   end if
!            
!      
!
!   call Cleanup ( )
!   return
!
!contains
!   !...............................................................................................................................
!   subroutine Cleanup()
!   ! This subroutine cleans up before exiting this subroutine
!   !...............................................................................................................................
!
!      ! if ( UnEcho > 0 ) close( UnEcho )
!      if (UnIn > 0) close ( UnIn )
!      
!   end subroutine Cleanup
!
!end subroutine ReadWD_DvrInputFile
!
!!----------------------------------------------------------------------------------------------------------------------------------
!subroutine ReadWDInputFile( InputFile, InputFileData, Default_DT, OutFileRoot, UnEcho, ErrStat, ErrMsg )
!! This subroutine reads the input file and stores all the data in the WD_InputFile structure.
!! It does not perform data validation.
!!..................................................................................................................................
!
!      ! Passed variables
!   real(DbKi),              intent(in)    :: Default_DT      ! The default DT (from glue code)
!
!   character(*),            intent(in)    :: InputFile       ! Name of the input file
!   character(*),            intent(in)    :: OutFileRoot     ! The rootname of all the output files written by this routine.
!
!   type(WD_InputFile),      intent(out)   :: InputFileData   ! Data stored in the module's input file
!   integer(IntKi),          intent(out)   :: UnEcho          ! Unit number for the echo file
!
!   integer(IntKi),          intent(out)   :: ErrStat         ! The error status code
!   character(*),            intent(out)   :: ErrMsg          ! The error message, if an error occurred
!
!      ! local variables
!
!   integer(IntKi)                         :: I
!   integer(IntKi)                         :: UnEc
!   integer(IntKi)                         :: UnIn            ! Unit number for reading file
!   logical                                :: Echo            ! Determines if an echo file should be written
!   character(1024)                        :: PriPath         ! Path name of the primary file
!   character(1024)                        :: FTitle          ! "File Title": the 2nd line of the input file, which contains a description of its contents
!   character(200)                         :: Line            ! Temporary storage of a line from the input file (to compare with "default")
!   integer(IntKi)                         :: errStat2, IOS   ! The error status code
!   character(ErrMsgLen)                   :: errMsg2         ! The error message, if an error occurred
!
!   character(*), parameter                :: RoutineName = 'ReadWDInputFile'
!   
!   
!      ! initialize values:
!
!   ErrStat = ErrID_None
!   ErrMsg  = ''
!   UnEcho  = -1
!   InputFileData%DT = Default_DT  ! the glue code's suggested DT for the module (may be overwritten in ReadPrimaryFile())
!
!      ! get the primary/platform input-file data
!      ! sets UnEcho, ADBlFile
!   
!     ! Initialize some variables:
!   ErrStat = ErrID_None
!   ErrMsg  = ""
!      
!   UnEc = -1
!   Echo = .FALSE.   
!   call GetPath( InputFile, PriPath )     ! Input files will be relative to the path where the primary input file is located.
!   
!
!   !call AllocAry( InputFileData%OutList, MaxOutPts, "Outlist", errStat2, errMsg2 )
!   !   call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!        
!      ! Get an available unit number for the file.
!
!   call GetNewUnit( UnIn, errStat2, errMsg2 )
!      call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!
!      ! Open the Primary input file.
!
!   call OpenFInpFile ( UnIn, InputFile, errStat2, errMsg2 )
!      call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!      if ( ErrStat >= AbortErrLev ) then
!         call Cleanup()
!         return
!      end if
!      
!                  
!      
!   ! Read the lines up/including to the "Echo" simulation control variable
!   ! If echo is FALSE, don't write these lines to the echo file. 
!   ! If Echo is TRUE, rewind and write on the second try.
!   
!   I = 1 !set the number of times we've read the file
!   do 
!   !----------- HEADER -------------------------------------------------------------
!   
!      call ReadCom( UnIn, InputFile, 'File header: Module Version (line 1)', errStat2, errMsg2, UnEc )
!         call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!   
!      call ReadStr( UnIn, InputFile, FTitle, 'FTitle', 'File Header: File Description (line 2)', errStat2, errMsg2, UnEc )
!         call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!         if ( ErrStat >= AbortErrLev ) then
!            call Cleanup()
!            return
!         end if
!   
!   
!   !----------- GENERAL OPTIONS ----------------------------------------------------
!   
!      call ReadCom( UnIn, InputFile, 'Section Header: General Options', errStat2, errMsg2, UnEc )
!         call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!   
!         ! Echo - Echo input to "<RootName>.WD.ech".
!   
!      call ReadVar( UnIn, InputFile, Echo, 'Echo',   'Echo flag', errStat2, errMsg2, UnEc )
!         call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!   
!   
!      if (.NOT. Echo .OR. I > 1) exit !exit this loop
!   
!         ! Otherwise, open the echo file, then rewind the input file and echo everything we've read
!      
!      I = I + 1         ! make sure we do this only once (increment counter that says how many times we've read this file)
!   
!      call OpenEcho ( UnEc, trim(OutFileRoot)//'.ech', errStat2, errMsg2, WD_Ver )
!         call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!         if ( ErrStat >= AbortErrLev ) then
!            call Cleanup()
!            return
!         end if
!   
!      if ( UnEc > 0 )  write (UnEc,'(/,A,/)')  'Data from '//trim(WD_Ver%Name)//' primary input file "'//trim( InputFile )//'":'
!         
!      rewind( UnIn, IOSTAT=errStat2 )  
!         if (errStat2 /= 0_IntKi ) then
!            call SetErrStat( ErrID_Fatal, 'Error rewinding file "'//trim(InputFile)//'".', ErrStat, ErrMsg, RoutineName )
!            call Cleanup()
!            return
!         end if         
!      
!   end do    
!
!   if (NWTC_VerboseLevel == NWTC_Verbose) then
!      call WrScr( ' Heading of the '//trim(WD_Ver%Name)//' input file: ' )      
!      call WrScr( '   '//trim( FTitle ) )
!   end if
!   
!   
!      ! DT - Time interval for Wake dynamics calculations {or default} (s):
!   Line = ""
!   call ReadVar( UnIn, InputFile, Line, "DT", "Time interval for Wake dynamics calculations {or default} (s)", errStat2, errMsg2, UnEc)
!      call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!         
!      call Conv2UC( Line )
!      if ( index(Line, "DEFAULT" ) /= 1 ) then ! If it's not "default", read this variable; otherwise use the value already stored in InputFileData%DT
!         READ( Line, *, IOSTAT=IOS) InputFileData%DT
!            call CheckIOS ( IOS, InputFile, 'DT', NumType, errStat2, errMsg2 )
!            call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
!      end if   
!      
!
!            
!      
!
!   call Cleanup ( )
!   return
!
!contains
!   !...............................................................................................................................
!   subroutine Cleanup()
!   ! This subroutine cleans up before exiting this subroutine
!   !...............................................................................................................................
!
!      ! if ( UnEcho > 0 ) close( UnEcho )
!      if (UnIn > 0) close ( UnIn )
!      
!   end subroutine Cleanup
!
!end subroutine ReadWDInputFile

subroutine WD_Dvr_Init( WD_InitInp,  WD_InitOut, WD_u,WD_p, WD_xd,  WD_y, errStat, errMsg)

   type(WD_InitInputType),   intent(  out)   :: WD_InitInp                         !< Input data for initialization routine
   type(WD_InitOutputType),     intent(  out)   :: WD_InitOut                         !< Input data for initialization routine
   type(WD_InputType),       intent(  out)   :: WD_u                            !< Input data for initialization routine
   type(WD_ParameterType),   intent(  out)   :: WD_p                         !< Input data for initialization routine
   type(WD_DiscreteStateType),   intent(  out)   :: WD_xd                         !< Input data for initialization routine
   type(WD_OutputType),      intent(  out)   :: WD_y                         !< Input data for initialization routine
                                                                                      
   integer(IntKi),           intent(out)     :: errStat                             !< Error status
   character(*),             intent(out)     :: errMsg                              !< Error message
   
   
      ! Local variables
   TYPE(ProgDesc), PARAMETER                    :: version = ProgDesc( 'WakeDynamics', 'v0.01', '4-Nov-2016')               ! The name, version, and date of AirfoilInfo.
   character(*), parameter                   :: RoutineName = 'WD_Dvr_Init'
   character(1024)                           :: OutFileRoot                         !< The rootname of the echo file, possibly opened in this routine
   integer(IntKi)                            :: errStat2                            ! local status of error message
   character(ErrMsgLen)                      :: errMsg2                             ! local error message if ErrStat /= ErrID_None
   integer(IntKi)                            :: i                                   ! loop counter
   character(1024)                           :: inputFile                           ! String to hold the driver input file name.
   character(1024)                           :: WD_InputFile                        ! String to hold the WakeDynamics input file name. 
   real(DbKi)                                :: interval                            ! Default timestep size
   
   character(1024)                           :: Rootname
   integer(IntKi) :: Unecho
   
   Unecho = 0
   Rootname = "Test"
   
      ! Initialize the library which handle file echos and WrScr, for example
   call NWTC_Init()
      
      ! Display the copyright notice
   CAlL DispCopyrightLicense( version )
   
      ! Tell our users what they're running
   call WrScr( ' Running '//GetNVD( version )//NewLine//' linked with '//trim( GetNVD( NWTC_Ver ))//NewLine )

   inputFile = ""  ! initialize to empty string to make sure it's input from the command line
   call CheckArgs( inputFile, ErrStat2 )
   if (len_trim(inputFile) == 0) then ! no input file was specified
      call SetErrStat(ErrID_Fatal, 'The required input file was not specified on the command line.', ErrStat, ErrMsg, RoutineName) 
                      
      if (BITS_IN_ADDR==32) then
         call NWTC_DisplaySyntax( inputFile, 'WakeDynamics_Driver_Win32.exe' )
      elseif( BITS_IN_ADDR == 64) then
         call NWTC_DisplaySyntax( inputFile, 'WakeDynamics_Driver_x64.exe' )
      else
         call NWTC_DisplaySyntax( inputFile, 'WakeDynamics_Driver.exe' )
      end if
         
      return
   end if       
   
      ! Read the driver input file, and populate some of the data for both the driver and the WD module
   !call ReadWD_DvrInputFile( inputFile, WD_InitInp, interval, RootName, UnEcho, errStat2, errMsg2 )
   
   
   !   ! Read the primary WakeDynamics input file
   !call ReadWDInputFile( WD_InputFile, WD_InitInp, interval, RootName, UnEcho, errStat2, errMsg2 )   
   !   call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName ) 
   !   if (ErrStat >= AbortErrLev) then
   !      call Cleanup()
   !      return
   !   end if
      
end subroutine WD_Dvr_Init




end module WD_Driver_Subs
   