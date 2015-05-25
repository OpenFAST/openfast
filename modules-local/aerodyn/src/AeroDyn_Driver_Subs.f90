module AeroDyn_Driver_Subs
   
   use AeroDyn_Driver_Types   
   use AeroDyn
    
   implicit none   
   
   TYPE(ProgDesc), PARAMETER   :: version   = ProgDesc( 'AeroDyn_driver', 'v1.00.00', '22-May-2015' )  ! The version number of this program.
                                                    
   contains

!----------------------------------------------------------------------------------------------------------------------------------
subroutine Dvr_Init(DvrData,errStat,errMsg )

   type(Dvr_SimData),            intent(  out) :: DvrData       ! driver data
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None

      ! local variables
   integer(IntKi)                              :: errStat2      ! local status of error message
   character(ErrMsgLen)                        :: errMsg2       ! local error message if ErrStat /= ErrID_None
   character(*), parameter                     :: RoutineName = 'Dvr_Init'

   CHARACTER(1000)                             :: inputFile     ! String to hold the file name.

   ErrStat = ErrID_None
   ErrMsg  = ""

   DvrData%OutFileData%unOutFile   = -1
   
      ! Initialize the library which handle file echos and WrScr, for example
   call NWTC_Init()
      
      ! Display the copyright notice
   CALL DispCopyrightLicense( version )
   
      ! Tell our users what they're running
   CALL WrScr( ' Running '//GetNVD( version )//NewLine//' linked with '//TRIM( GetNVD( NWTC_Ver ))//NewLine )

   InputFile = ""  ! initialize to empty string to make sure it's input from the command line
   CALL CheckArgs( InputFile, ErrStat2 )
   IF (LEN_TRIM(InputFile) == 0) THEN ! no input file was specified
      call SetErrStat(ErrID_Fatal, 'The required input file was not specified on the command line.', ErrStat, ErrMsg, RoutineName) 
      
         !bjj:  if people have compiled themselves, they should be able to figure out the file name, right?         
      IF (BITS_IN_ADDR==32) THEN
         CALL NWTC_DisplaySyntax( InputFile, 'AeroDyn_Driver_Win32.exe' )
      ELSEIF( BITS_IN_ADDR == 64) THEN
         CALL NWTC_DisplaySyntax( InputFile, 'AeroDyn_Driver_x64.exe' )
      ELSE
         CALL NWTC_DisplaySyntax( InputFile, 'AeroDyn_Driver.exe' )
      END IF
         
      return
   END IF        
         
      ! Read the AeroDyn driver input file
   call Dvr_ReadInputFile(inputFile, DvrData, errStat2, errMsg2 )
      call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName) 

end subroutine Dvr_Init 
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Init_AeroDyn(DvrData, AD, dt, errStat, errMsg)

   type(Dvr_SimData),            intent(inout) :: DvrData       ! Input data for initialization
   type(AeroDyn_Data),           intent(inout) :: AD            ! AeroDyn data 
   real(DbKi),                   intent(inout) :: dt            ! interval
      
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None

      ! locals
   integer(IntKi)                              :: k
   integer(IntKi)                              :: errStat2      ! local status of error message
   character(ErrMsgLen)                        :: errMsg2       ! local error message if ErrStat /= ErrID_None
   character(*), parameter                     :: RoutineName = 'Init_AeroDyn'
                                                  
   ! local data                                
   type(AD_InitInputType)                      :: InitInData     ! Input data for initialization
   type(AD_InitOutputType)                     :: InitOutData    ! Output data from initialization
      
      
      
   errStat = ErrID_None
   errMsg  = ''
   
   InitInData%InputFile      = DvrData%AD_InputFile
   InitInData%NumBlades      = DvrData%numBlade
   InitInData%RootName       = DvrData%outFileData%Root
                        
   
      ! set initialization data:
   call AllocAry( InitInData%BladeRootPosition, 3, InitInData%NumBlades, 'BladeRootPosition', errStat2, ErrMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   call AllocAry( InitInData%BladeRootOrientation, 3, 3, InitInData%NumBlades, 'BladeRootOrientation', errStat2, ErrMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         
   if (ErrStat >= AbortErrLev) then
      call Cleanup()
      return
   end if
      
     
      ! bjj: fix this!
   do k=1,InitInData%numBlades
         
      InitInData%BladeRootPosition(:,k)  = (/ 0.0_ReKi, 0.0_ReKi, DvrData%hubRad /)
         
      !AD_Data%BladeRootOrientation(:,:,k) = 
   end do
      
   !AD_Data%HubPosition = 
   !AD_Data%HubOrientation =
      
call AD_Init(InitInData, AD%u(1), AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%y, dt, InitOutData, ErrStat2, ErrMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      
      
   ! move AD initOut data to WTP
call move_alloc( InitOutData%twist,  DvrData%twist )
call move_alloc( InitOutData%rLocal, DvrData%rLoc )
DvrData%RotorRad = InitOutData%RotorRad
      
call move_alloc( InitOutData%WriteOutputHdr, DvrData%OutFileData%WriteOutputHdr )
call move_alloc( InitOutData%WriteOutputUnt, DvrData%OutFileData%WriteOutputUnt )   
     
DvrData%OutFileData%AD_ver = InitOutData%ver
   
contains
   subroutine cleanup()
      call AD_DestroyInitInput( InitInData, ErrStat2, ErrMsg2 )   
      call AD_DestroyInitOutput( InitOutData, ErrStat2, ErrMsg2 )      
   end subroutine cleanup
   
end subroutine Init_AeroDyn
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Set_AD_Inputs(iCase,nt,time,DvrData,AD,errStat,errMsg)

   integer(IntKi)              , intent(in   ) :: iCase         ! case number 
   integer(IntKi)              , intent(in   ) :: nt            ! time step number
   real(DbKi)                  , intent(  out) :: time          ! time in seconds
   
   type(Dvr_SimData),            intent(inout) :: DvrData       ! Driver data 
   type(AeroDyn_Data),           intent(inout) :: AD            ! AeroDyn data 
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None

      ! local variables
   integer(IntKi)                              :: errStat2      ! local status of error message
   character(ErrMsgLen)                        :: errMsg2       ! local error message if ErrStat /= ErrID_None
   character(*), parameter                     :: RoutineName = 'Set_AD_Inputs'

   integer(intKi)                              :: j             ! loop counter for nodes
   integer(intKi)                              :: k             ! loop counter for blades

   real(ReKi)                                  :: psiRotor, velocityHub, deltar
   real(ReKi)                                  :: z             ! heigh (m)
   
      time = (nt-1) * DvrData%Cases(iCase)%dT
      AD%inputTime(1) = time
      
      !bjj: TODO FIX ME: set TowerMotion, HubMotion, BladeRootMotion, BladeMotion, InflowOnBlade, InflowOnTower

      !bjj: FIX ME!!!!
      ! Tower motions:
      do j=1,AD%u(1)%TowerMotion%nnodes
         AD%u(1)%TowerMotion%Orientation(  :,:,j) = AD%u(1)%TowerMotion%RefOrientation(:,:,j)
         AD%u(1)%TowerMotion%TranslationDisp(:,j) = 0.0_ReKi
         AD%u(1)%TowerMotion%TranslationVel( :,j) = 0.0_ReKi
         AD%u(1)%TowerMotion%RotationVel(    :,j) = 0.0_ReKi
      end do !j=nnodes
      
      ! Hub motions:
      AD%u(1)%HubMotion%Orientation(  :,:,1) = AD%u(1)%HubMotion%RefOrientation(:,:,1)
      AD%u(1)%HubMotion%TranslationDisp(:,1) = 0.0_ReKi
      AD%u(1)%HubMotion%TranslationVel( :,1) = 0.0_ReKi
      AD%u(1)%HubMotion%RotationVel(    :,1) = 0.0_ReKi
      
      ! Blade root motions:
      do k=1,DvrData%numBlade
         AD%u(1)%BladeRootMotion(k)%Orientation(  :,:,1) = AD%u(1)%BladeRootMotion(k)%RefOrientation(:,:,1)
         AD%u(1)%BladeRootMotion(k)%TranslationDisp(:,1) = 0.0_ReKi
         AD%u(1)%BladeRootMotion(k)%TranslationVel( :,1) = 0.0_ReKi
         AD%u(1)%BladeRootMotion(k)%RotationVel(    :,1) = 0.0_ReKi
      end do !k=numBlade
      
      ! Blade motions:
      do k=1,DvrData%numBlade
         do j=1,AD%u(1)%BladeMotion(k)%nnodes
            AD%u(1)%BladeMotion(k)%Orientation(  :,:,j) = AD%u(1)%BladeMotion(k)%RefOrientation(:,:,j)
            AD%u(1)%BladeMotion(k)%TranslationDisp(:,j) = 0.0_ReKi
            AD%u(1)%BladeMotion(k)%TranslationVel( :,j) = 0.0_ReKi
            AD%u(1)%BladeMotion(k)%RotationVel(    :,j) = 0.0_ReKi
         end do !j=nnodes
      end do !k=numBlade
      
      ! Inflow wind velocities:
      ! InflowOnBlade
      do k=1,DvrData%numBlade
         do j=1,AD%u(1)%BladeMotion(k)%nnodes
            z = AD%u(1)%BladeMotion(k)%Position(3,j) + AD%u(1)%BladeMotion(k)%TranslationDisp(3,j)
            AD%u(1)%InflowOnBlade(1,j,k) = GetU(  DvrData%Cases(iCase)%WndSpeed, DvrData%HubHt, DvrData%Cases(iCase)%ShearExp, z )
            AD%u(1)%InflowOnBlade(2,j,k) = 0.0_ReKi !V
            AD%u(1)%InflowOnBlade(3,j,k) = 0.0_ReKi !W         
         end do !j=nnodes
      end do !k=numBlade
      
      !InflowOnTower
      do j=1,AD%u(1)%TowerMotion%nnodes
         z = AD%u(1)%TowerMotion%Position(3,j) + AD%u(1)%TowerMotion%TranslationDisp(3,j)
         AD%u(1)%InflowOnTower(1,j) = GetU(  DvrData%Cases(iCase)%WndSpeed, DvrData%HubHt, DvrData%Cases(iCase)%ShearExp, z )
         AD%u(1)%InflowOnTower(2,j) = 0.0_ReKi !V
         AD%u(1)%InflowOnTower(3,j) = 0.0_ReKi !W         
      end do !j=nnodes
      
      
      
      !! Angular position of the rotor, starts at 0.0 
      !psiRotor = time*DvrData%Cases(iCase)%RotSpeed ! in s * rad/s 
      !
      !   ! Velocity of inflow wind at the hub height
      !velocityHub = DvrData%Cases(iCase)%WndSpeed ! m/s
      !   ! Zero yaw for now
      !AD%u(1)%gamma = DvrData%Cases(iCase)%Yaw
      !   ! Rotor angular velocity
      !AD%u(1)%omega = DvrData%Cases(iCase)%RotSpeed  
      !   ! Average tip-speed ratio
      !AD%u(1)%lambda = DvrData%Cases(iCase)%TSR
      !  
      !AD%u(1)%Vinf   = velocityHub
      !
      !do k=1,DvrData%numBlade
      !   AD%u(1)%psi(k)   = psiRotor + (k-1)*2.0*pi/(DvrData%numBlade) ! find psi for each blade based on the rotor psi value
      !   AD%u(1)%rTip(k)  = DvrData%RotorRad  ! set the tip radius to the initialization distance along the blade (straight blade assumption)
      !
      !   do j=1,size(DvrData%Twist,1)
      !      AD%u(1)%theta (j,k)  = DvrData%Twist(j,k) + DvrData%Cases(iCase)%Pitch
      !         u(1)%rLocal(j,k)  = DvrData%rLoc(j,k)
      !         u(1)%Vx    (j,k)  =  velocityHub*cos(DvrData%Cases(iCase)%Precone)
      !         u(1)%Vy    (j,k)  =  u(1)%omega*u(1)%rLocal(j,k)*cos(DvrData%Cases(iCase)%Precone)                                        
      !   end do !i=1,size(DvrData%Twist,1) # blade nodes
      !end do !j=1,DvrData%numBlade
         
   
   
end subroutine Set_AD_Inputs
!----------------------------------------------------------------------------------------------------------------------------------
function GetU(  URef, ZRef, PLExp, z ) result (U)
   real(ReKi), intent(in) :: URef
   real(ReKi), intent(in) :: ZRef
   real(ReKi), intent(in) :: PLExp
   real(ReKi), intent(in) :: z
   real(ReKi)             :: U
   
   U = URef*(z/ZRef)**PLExp

end function GetU
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Dvr_ReadInputFile(fileName, DvrData, errStat, errMsg )
   ! This routine opens the gets the data from the input files.

   character(*),                  intent( in    )   :: fileName
   type(Dvr_SimData),             intent(   out )   :: DvrData
   integer,                       intent(   out )   :: errStat              ! returns a non-zero value when an error occurs  
   character(*),                  intent(   out )   :: errMsg               ! Error message if errStat /= ErrID_None
   

      ! Local variables
   character(1024)              :: PriPath
   character(1024)              :: inpVersion                               ! String containing the input-version information.
   character(1024)              :: line                                     ! String containing a line of input.
   integer                      :: unIn, unEc, NumElmPr
   integer                      :: ISeg, ICase
   integer                      :: IOS, Sttus
   character( 11)               :: DateNow                                  ! Date shortly after the start of execution.
   character(  8)               :: TimeNow                                  ! Time of day shortly after the start of execution.
   real(ReKi)                   :: InpCase(10)                              ! Temporary array to hold combined-case input parameters.
   logical                      :: TabDel      
   logical                      :: echo   

   INTEGER(IntKi)               :: ErrStat2                                        ! Temporary Error status
   CHARACTER(ErrMsgLen)         :: ErrMsg2                                         ! Temporary Err msg
   CHARACTER(*), PARAMETER      :: RoutineName = 'Dvr_ReadInputFile'
   
   
   ErrStat = ErrID_None
   ErrMsg  = ''
   UnIn = -1
   UnEc = -1
   
   ! Open the input file
   CALL GetPath( fileName, PriPath )     ! Input files will be relative to the path where the primary input file is located.

   call GetNewUnit( unIn )   
   call OpenFInpFile( unIn, fileName, errStat2, ErrMsg2 )
   call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   if ( errStat >= AbortErrLev ) then
      call cleanup()
      return
   end if

   
   call WrScr( 'Opening WT_Perf input file:  '//fileName )

      ! Skip a line, read the run title and the version information.

   CALL ReadStr( UnIn, fileName, inpVersion, 'inpVersion', 'File Header: (line 1)', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
   CALL ReadStr( UnIn, fileName, DvrData%OutFileData%runTitle, 'runTitle', 'File Header: File Description (line 2)', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
   call WrScr1 ( ' '//DvrData%OutFileData%runTitle )
   
      ! Read in the title line for the input-configuration subsection.
   CALL ReadStr( UnIn, fileName, line, 'line', 'File Header: (line 3)', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


      ! See if we should echo the output.     
   call ReadVar ( unIn, fileName, echo, 'Echo', 'Echo Input', errStat2, errMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   if ( echo )  then
         ! Get date and time.
      dateNow = CurDate()
      timeNow = CurTime()
      call GetNewUnit( unEc ) 
      call getroot(fileName,DvrData%OutFileData%Root)      
      call  OpenFOutFile ( unEc, trim( DvrData%OutFileData%Root )//'.ech', errStat2, errMsg2 )
         call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
         if ( errStat >= AbortErrLev ) then
            call cleanup()
            return
         end if
      
      write (unEc,'(A)')      'Echo of Input File:'
      write (unEc,'(A)')      ' "'//fileName//'"'
      write (unEc,'(A)')      'Generated on: '//trim( dateNow )//' at '//trim( timeNow )//'.'
      write (unEc,'(A)')      inpVersion
      write (unEc,'(A)')      DvrData%OutFileData%runTitle
      write (unEc,'(A)')      line
      write (unEc,Ec_LgFrmt)  echo, 'Echo', 'Echo input parameters to "rootname.ech"?'
   end if


      ! Read the rest of input-configuration section.
      
   call ReadVar ( unIn, fileName, DvrData%AD_InputFile,   'AD_InputFile',   'Name of the AeroDyn input file', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
      if ( errStat >= AbortErrLev ) then
         call cleanup()
         return
      end if
   IF ( PathIsRelative( DvrData%AD_InputFile ) ) DvrData%AD_InputFile = TRIM(PriPath)//TRIM(DvrData%AD_InputFile)



      ! Read the model-configuration section.

   call ReadCom ( unIn, fileName,                                 'the model-configuration subtitle', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )

   
      ! Read the turbine-data section.

   call ReadCom ( unIn, fileName, 'the turbine-data subtitle', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadVar ( unIn, fileName, DvrData%NumBlade, 'NumBlade', 'Number of blades.', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadVar ( unIn, fileName, DvrData%HubRad,   'HubRad',   'Hub radius.', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadVar ( unIn, fileName, DvrData%HubHt,    'HubHt',    'Hub height.', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )

      if ( errStat >= AbortErrLev ) then
         call cleanup()
         return
      end if           


      ! Read the I/O-configuration section.

   call ReadCom ( unIn, fileName, 'the I/O-configuration subtitle', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadVar ( unIn, fileName, DvrData%OutFileData%Root, 'OutFileRoot', 'Root name for any output files', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   if (len_trim(DvrData%OutFileData%Root) == 0) then
      call getroot(fileName,DvrData%OutFileData%Root)
   end if
   
   call ReadVar ( unIn, fileName, TabDel,   'TabDel',   'Make output tab-delimited (fixed-width otherwise)?', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
      if (TabDel) then
         DvrData%OutFileData%delim = TAB
      else
         DvrData%OutFileData%delim = " "
      end if
               
   call ReadVar ( unIn, fileName, Beep,               'Beep',     'Beep on exit?', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName ) !bjj: this is a global variable in NWTC_Library
   call ReadVar ( unIn, fileName, DvrData%InputTSR, 'InputTSR', 'Input speeds as TSRs?', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
      if ( errStat >= AbortErrLev ) then
         call cleanup()
         return
      end if


      ! Read the combined-case section.

   call ReadCom  ( unIn, fileName, 'the combined-case subtitle', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadVar  ( unIn, fileName, DvrData%NumCases, 'NumCases', 'Number of cases to run.', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadCom  ( unIn, fileName, 'the combined-case-block header', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )

      if ( errStat >= AbortErrLev ) then
         call cleanup()
         return
      end if
      
   IF ( DvrData%NumCases < 0 )  THEN

      call setErrStat( ErrID_Fatal,'Variable "NumCases" must be >= 0.  Instead, it is "'//TRIM( Int2LStr( DvrData%NumCases ) )//'".' ,errstat,errmsg,routinename)
      call cleanup()
      return

   ELSEIF ( DvrData%NumCases > 0 )  THEN

      ALLOCATE ( DvrData%Cases(DvrData%NumCases) , STAT=Sttus )
      IF ( Sttus /= 0 )  THEN
         call setErrStat( ErrID_Fatal,'Error allocating memory for the Cases array.',errstat,errmsg,routinename)
         call cleanup()
         return
      ENDIF

      DO ICase=1,DvrData%NumCases

         CALL ReadAry ( unIn, fileName, InpCase,  10, 'InpCase',  'Wind Speed or TSR, Rotor Speed, and Pitch for Case #' &
                       //TRIM( Int2LStr( ICase ) )//'.', errStat2, errMsg2, UnEc )
            call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
            
         DvrData%Cases(iCase)%TSR             = InpCase( 1) ! we'll calculate WS and TSR later, after we get RotorRad from AD 
         DvrData%Cases(iCase)%WndSpeed        = InpCase( 1) ! we'll calculate WS and TSR later, after we get RotorRad from AD          
         DvrData%Cases(ICase)%ShearExp        = InpCase( 2)
         DvrData%Cases(ICase)%RotSpeed        = InpCase( 3)*RPM2RPS
         DvrData%Cases(ICase)%Pitch           = InpCase( 4)*D2R
         DvrData%Cases(ICase)%Yaw             = InpCase( 5)*D2R
         DvrData%Cases(ICase)%Tilt            = InpCase( 6)*D2R
         DvrData%Cases(iCase)%PreCone         = InpCase( 7)*D2R
         DvrData%Cases(iCase)%AzAng0          = InpCase( 8)*D2R
         DvrData%Cases(iCase)%dT              = InpCase( 9)*D2R
         DvrData%Cases(iCase)%Tmax            = InpCase(10)*D2R
               
      ENDDO ! ICase
   END IF
   
   call cleanup ( )


   RETURN
contains
   subroutine cleanup()
      if (UnIn>0) close(UnIn)
      if (UnEc>0) close(UnEc)
   end subroutine cleanup
END SUBROUTINE Dvr_ReadInputFile
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Dvr_WriteOutputLine(OutFileData, t, output, errStat, errMsg)

   real(DbKi)             ,  intent(in   )   :: t                    ! simulation time (s)
   type(Dvr_OutputFile)   ,  intent(in   )   :: OutFileData
   real(ReKi)             ,  intent(in   )   :: output(:)            ! Rootname for the output file
   integer(IntKi)         ,  intent(inout)   :: errStat              ! Status of error message
   character(*)           ,  intent(inout)   :: errMsg               ! Error message if ErrStat /= ErrID_None
      
   ! Local variables.

   integer(IntKi)                   :: i                                         ! loop counter
   integer(IntKi)                   :: indxLast                                  ! The index of the last row value to be written to AllOutData for this time step (column).
   integer(IntKi)                   :: indxNext                                  ! The index of the next row value to be written to AllOutData for this time step (column).

   character(200)                   :: frmt                                      ! A string to hold a format specifier
   character(15)                    :: tmpStr                                    ! temporary string to print the time output as text
integer :: numOuts
   
   errStat = ErrID_None
   errMsg  = ''
   numOuts = size(output,1)
   frmt = '"'//OutFileData%delim//'"'//trim(OutFileData%outFmt)      ! format for array elements from individual modules
   
      ! time
   write( tmpStr, '(F15.4)' ) t
   call WrFileNR( OutFileData%unOutFile, tmpStr )
   call WrReAryFileNR ( OutFileData%unOutFile, output,  frmt, errStat, errMsg )
   if ( errStat >= AbortErrLev ) return
   
     ! write a new line (advance to the next line)
   write (OutFileData%unOutFile,'()')
      
end subroutine Dvr_WriteOutputLine
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Dvr_InitializeOutputFile( iCase, CaseData, OutFileData, errStat, errMsg)
      type(Dvr_OutputFile),     intent(inout)   :: OutFileData 
      
      integer(IntKi)         ,  intent(in   )   :: iCase                ! case number (to write in file description line and use for file name)
      type(Dvr_Case),           intent(in   )   :: CaseData
      
      integer(IntKi)         ,  intent(inout)   :: errStat              ! Status of error message
      character(*)           ,  intent(inout)   :: errMsg               ! Error message if ErrStat /= ErrID_None

         ! locals
      integer(IntKi)                            ::  j
      
      integer(IntKi)                            :: numOuts
      
      character(200)                            :: frmt ,frmt2         ! A string to hold a format specifier

      
      numOuts = size(OutFileData%WriteOutputHdr)
      
      call GetNewUnit( OutFileData%unOutFile, ErrStat, ErrMsg )
         if ( ErrStat >= AbortErrLev ) then
            OutFileData%unOutFile = -1
            return
         end if
         

      call OpenFOutFile ( OutFileData%unOutFile, trim(outFileData%Root)//'.ADdrv.'//trim(num2lstr(iCase))//'.out', ErrStat, ErrMsg )
         if ( ErrStat >= AbortErrLev ) return
         
      WRITE (OutFileData%unOutFile,'(/,A)')  'Predictions were generated on '//CurDate()//' at '//CurTime()//' using '//trim(GetNVD(version))
      WRITE (OutFileData%unOutFile,'(1X,A)') trim(GetNVD(OutFileData%AD_ver))
      WRITE (OutFileData%unOutFile,'()' )    !print a blank line
      WRITE (OutFileData%unOutFile,'(A,11(1x,A,"=",ES11.4e2,1x,A))'   ) 'Case '//trim(num2lstr(iCase))//':' &
         ,'WndSpeed', CaseData%WndSpeed, 'm/s' &
         ,'ShearExp', CaseData%ShearExp, '' &
         ,'TSR',      CaseData%TSR, '' &
         ,'RotSpeed', CaseData%RotSpeed*RPS2RPM,'rpm' &
         ,'Pitch',    CaseData%Pitch*R2D, 'deg' &
         ,'Yaw',      CaseData%Yaw*R2D, 'deg' &
         ,'Tilt',     CaseData%Tilt*R2D, 'deg' &
         ,'Precone',  CaseData%Precone*R2D, 'deg' &
         ,'AzAng0',   CaseData%AzAng0*R2D, 'deg' &
         ,'dT',       CaseData%dT, 's' &
         ,'Tmax',     CaseData%Tmax,'s'
      
      WRITE (OutFileData%unOutFile,'()' )    !print a blank line
         
      

         !......................................................
         ! Write the names of the output parameters on one line:
         !......................................................
      frmt = '"'//OutFileData%delim//'"A15'      ! format for array elements 
      frmt2 = '(A,'//trim(num2lstr(numOuts))//'('//trim(frmt)//'))'
               
      write (OutFileData%unOutFile,frmt2,IOSTAT=errStat)  '     Time           ', OutFileData%WriteOutputHdr
      if ( errStat /= 0 ) then
         errMsg = 'Error '//trim(Num2LStr(errStat))//' occurred while writing to file in Dvr_InitializeOutputFile() using this format: '&
                  //trim(frmt2)
         errStat = ErrID_Fatal
         return
      end if      

         !......................................................
         ! Write the units of the output parameters on one line:
         !......................................................

      write (OutFileData%unOutFile,frmt2,IOSTAT=errStat)  '      (s)           ', OutFileData%WriteOutputUnt
      if ( errStat /= 0 ) then
         errMsg = 'Error '//trim(Num2LStr(errStat))//' occurred while writing to file in Dvr_InitializeOutputFile() using this format: '&
                  //trim(frmt2)
         errStat = ErrID_Fatal
         return
      end if
      
end subroutine Dvr_InitializeOutputFile
!----------------------------------------------------------------------------------------------------------------------------------
end module AeroDyn_Driver_Subs