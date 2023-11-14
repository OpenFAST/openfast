module FAST_Farm_IO
   
   USE NWTC_Library
   USE FAST_Farm_Types
   USE FAST_Farm_IO_Params
   
   IMPLICIT NONE
   
   TYPE(ProgDesc), PARAMETER      :: Farm_Ver      = ProgDesc( 'FAST.Farm', '', '' ) !< module date/version information 

   

      contains

!----------------------------------------------------------------------------------------------------------------------------------
!> This function returns a string describing the glue code and some of the compilation options we're using.
FUNCTION GetVersion(ThisProgVer)

   ! Passed Variables:

   TYPE(ProgDesc), INTENT( IN    ) :: ThisProgVer     !< program name/date/version description
   CHARACTER(1024)                 :: GetVersion      !< String containing a description of the compiled precision.

   GetVersion = TRIM(GetNVD(ThisProgVer))//', compiled'
   
   
   GetVersion = TRIM(GetVersion)//' as a '//TRIM(Num2LStr(BITS_IN_ADDR))//'-bit application using'
   
   ! determine precision

      IF ( ReKi == SiKi )  THEN     ! Single precision
         GetVersion = TRIM(GetVersion)//' single'
      ELSEIF ( ReKi == R8Ki )  THEN ! Double precision
         GetVersion = TRIM(GetVersion)// ' double'
      ELSE                          ! Unknown precision
         GetVersion = TRIM(GetVersion)//' unknown'
      ENDIF

!   GetVersion = TRIM(GetVersion)//' precision with '//OS_Desc
   GetVersion = TRIM(GetVersion)//' precision'


   RETURN
END FUNCTION GetVersion

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine generates the summary file, which contains a regurgitation of  the input data and interpolated flexible body data.
SUBROUTINE Farm_PrintSum( farm, WD_InputFileData, ErrStat, ErrMsg )

      ! Passed variables
   type(All_FastFarm_Data),        INTENT(IN )          :: farm                                  !< FAST.Farm data
   type(WD_InputFileType),         INTENT(IN )          :: WD_InputFileData                      !< Wake Dynamics Input File data
   INTEGER(IntKi),                 INTENT(OUT)          :: ErrStat                               !< Error status
   CHARACTER(*),                   INTENT(OUT)          :: ErrMsg                                !< Error message corresponding to ErrStat


      ! Local variables.

   INTEGER(IntKi)               :: I,J                                             ! Index for the nodes.
   INTEGER(IntKi)               :: K                                               ! Generic index (also for the blade number).
   INTEGER(IntKi)               :: UnSum                                            ! I/O unit number for the summary output file

   CHARACTER(*), PARAMETER      :: Fmt1      = "(34X,3(6X,'Blade',I2,:))"          ! Format for outputting blade headings.
   CHARACTER(*), PARAMETER      :: Fmt2      = "(34X,3(6X,A,:))"                   ! Format for outputting blade headings.
   CHARACTER(*), PARAMETER      :: FmtDat    = '(A,T35,3(:,F13.3))'                ! Format for outputting mass and modal data.
   CHARACTER(*), PARAMETER      :: FmtDatT   = '(A,T35,1(:,F13.8))'                ! Format for outputting time steps.
   CHARACTER(100)               :: RotorType                                       ! Text description of rotor.
   CHARACTER(30)                :: Fmt
   CHARACTER(30)                :: OutPFmt                                         ! Format to print list of selected output channels to summary file
   CHARACTER(10)                :: DOFEnabled                                      ! String to say if a DOF is enabled or disabled
   CHARACTER(3)                 :: outStr
   CHARACTER(10)                :: CalWakeDiamStr
   CHARACTER(100)               :: strModDescr
   
   ! Open the summary file and give it a heading.
   
   CALL GetNewUnit( UnSum, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN

   CALL OpenFOutFile ( UnSum, TRIM( farm%p%OutFileRoot )//'.sum', ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN

   
   ! Heading:
   !.......................... Module Versions .....................................................


   WRITE (UnSum,'(A)') 'FAST.Farm Summary File'
   WRITE (UnSum,'(/A)')  TRIM( farm%p%FileDescLines(1) )

   WRITE (UnSum,'(2X,A)'   )  'compiled with'
   Fmt = '(4x,A)'
   WRITE (UnSum,Fmt)  TRIM( GetNVD( NWTC_Ver ) )
   WRITE (UnSum,Fmt)  TRIM( GetNVD( farm%p%Module_Ver( ModuleFF_SC    ) ) )
   WRITE (UnSum,Fmt)  TRIM( GetNVD( farm%p%Module_Ver( ModuleFF_FWrap ) ) )
   WRITE (UnSum,Fmt)  TRIM( GetNVD( farm%p%Module_Ver( ModuleFF_WD    ) ) )
   WRITE (UnSum,Fmt)  TRIM( GetNVD( farm%p%Module_Ver( ModuleFF_AWAE  ) ) )
   !WRITE (y_FAST%UnSum,Fmt)  TRIM( GetNVD(  ) )

   WRITE (UnSum,'(/,A)') 'Description from the FAST.Farm input file: '//trim(farm%p%FTitle)
   
   WRITE (UnSum,'(/,A)') 'Ambient Wind:'
   
   if (     farm%AWAE%p%mod_AmbWind == 1 ) then
      strModDescr = 'High-Fidelity Precursor'
   elseif ( farm%AWAE%p%mod_AmbWind == 2 ) then
      strModDescr = 'One InflowWind Module'
   else   ! farm%AWAE%p%mod_AmbWind == 3
      strModDescr = 'Multiple InflowWind Modules'
   end if
   
   WRITE (UnSum,'(2X,A)') 'Ambient wind model: '//trim(strModDescr)
   if ( farm%AWAE%p%mod_AmbWind == 1 ) then
      WRITE (UnSum,'(2X,A)') 'Ambient wind input filepath: '//trim(farm%p%WindFilePath)
   else
      WRITE (UnSum,'(2X,A)') 'InflowWind module input file: '//trim(farm%p%WindFilePath)
   end if
   
   !..................................
   ! Turbine information.
   !..................................
   
   WRITE (UnSum,'(/,A)'   )  'Wind Turbines: '//trim(Num2LStr(farm%p%NumTurbines))
   WRITE (UnSum,'(2X,A)')  'Turbine Number  Output Turbine Number      X          Y          Z   OpenFAST Time Step  OpenFAST SubCycles  OpenFAST Input File'
   WRITE (UnSum,'(2X,A)')  '      (-)                (-)              (m)        (m)        (m)          (S)             (-)                (-)'

   do I = 1,farm%p%NumTurbines
      if ( I < 10 ) then
         outStr = 'T'//(trim(Num2LStr(I)))
      else
         outStr = ' - '
      end if

      WRITE(UnSum,'(6X,I4,17X,A,8X,3(1X,F10.3),5X,F10.5,8X,I4,10X,A)')  I, outStr, farm%p%WT_Position(:,I), (farm%p%DT_low/real(farm%FWrap(I)%p%n_FAST_low)), farm%FWrap(I)%p%n_FAST_low, trim(farm%p%WT_FASTInFile(I))
                                                                      
   end do
   
   WRITE (UnSum,'(/,A)'   )  'Wake Dynamics Finite-Difference Grid: '//trim(Num2LStr(farm%WD(1)%p%NumRadii))//' Radii, '//trim(Num2LStr(farm%WD(1)%p%NumPlanes))//' Planes'
   WRITE (UnSum,'(2X,A)')      'Radial Node Number  Output Node Number    Radius'
   WRITE (UnSum,'(2X,A)')      '       (-)                (-)              (m) '  
   do I = 0, farm%WD(1)%p%NumRadii-1
      outStr = ' - '
      do J = 1, farm%p%NOutRadii
         if (farm%p%OutRadii(J) == I ) then
            outStr = 'N'//trim(Num2LStr(J))
            exit
         end if
      end do
      
      WRITE(UnSum,'(8X,I4,16X,A3,9X,F10.3)')  I, outStr, farm%WD(1)%p%r(I)
      
   end do
   
   WRITE (UnSum,'(/,A)'   )  'Wake Dynamics Parameters'
   WRITE (UnSum,'(2X,A)')      'Cut-off (corner) frequency of the low-pass time-filter for the wake advection, deflection, and meandering model (Hz): '//trim(Num2LStr(WD_InputFileData%f_c))
   WRITE (UnSum,'(4X,A)')        '( low-pass time-filter parameter (-): '//trim(Num2LStr(farm%WD(1)%p%filtParam))//' )'
   WRITE (UnSum,'(2X,A)')      'Calibrated parameter in the correction for wake deflection defining the horizontal offset at the rotor (m): '//trim(Num2LStr(farm%WD(1)%p%C_HWkDfl_O))
   WRITE (UnSum,'(2X,A)')      'Calibrated parameter in the correction for wake deflection defining the horizontal offset at the rotor scaled with yaw error (m/deg): '//trim(Num2LStr(farm%WD(1)%p%C_HWkDfl_OY/R2D)) 
   WRITE (UnSum,'(2X,A)')      'Calibrated parameter in the correction for wake deflection defining the horizontal offset scaled with downstream distance (-):  '//trim(Num2LStr(farm%WD(1)%p%C_HWkDfl_x))
   WRITE (UnSum,'(2X,A)')      'Calibrated parameter in the correction for wake deflection defining the horizontal offset scaled with downstream distance and yaw error (1/deg):  '//trim(Num2LStr(farm%WD(1)%p%C_HWkDfl_xY/R2D))
   WRITE (UnSum,'(2X,A)')      'Calibrated parameter for near-wake correction (-): '//trim(Num2LStr(farm%WD(1)%p%C_NearWake))
   WRITE (UnSum,'(2X,A)')      'Calibrated parameter for the influence of ambient turbulence in the eddy viscosity (-):  '//trim(Num2LStr(farm%WD(1)%p%k_vAmb))
   WRITE (UnSum,'(2X,A)')      'Calibrated parameter for the influence of the shear layer in the eddy viscosity (-):  '//trim(Num2LStr(farm%WD(1)%p%k_vShr))
   WRITE (UnSum,'(2X,A)')      'Calibrated parameter in the eddy viscosity filter function for ambient turbulence defining the transitional diameter fraction between the minimum and exponential regions (-):  '//trim(Num2LStr(farm%WD(1)%p%C_vAmb_DMin))
   WRITE (UnSum,'(2X,A)')      'Calibrated parameter in the eddy viscosity filter function for ambient turbulence defining the transitional diameter fraction between the exponential and maximum regions (-):  '//trim(Num2LStr(farm%WD(1)%p%C_vAmb_DMax))
   WRITE (UnSum,'(2X,A)')      'Calibrated parameter in the eddy viscosity filter function for ambient turbulence defining the value in the minimum region (-):  '//trim(Num2LStr(farm%WD(1)%p%C_vAmb_FMin))
   WRITE (UnSum,'(2X,A)')      'Calibrated parameter in the eddy viscosity filter function for ambient turbulence defining the exponent in the exponential region (-):  '//trim(Num2LStr(farm%WD(1)%p%C_vAmb_Exp))
   WRITE (UnSum,'(2X,A)')      'Calibrated parameter in the eddy viscosity filter function for the shear layer defining the transitional diameter fraction between the minimum and exponential regions (-):  '//trim(Num2LStr(farm%WD(1)%p%C_vShr_DMin))
   WRITE (UnSum,'(2X,A)')      'Calibrated parameter in the eddy viscosity filter function for the shear layer defining the transitional diameter fraction between the exponential and maximum regions (-): '//trim(Num2LStr(farm%WD(1)%p%C_vShr_DMax)) 
   WRITE (UnSum,'(2X,A)')      'Calibrated parameter in the eddy viscosity filter function for the shear layer defining the functional value in the minimum region (-):  '//trim(Num2LStr(farm%WD(1)%p%C_vShr_FMin))
   WRITE (UnSum,'(2X,A)')      'Calibrated parameter in the eddy viscosity filter function for the shear layer defining the exponent in the exponential region (-):  '//trim(Num2LStr(farm%WD(1)%p%C_vShr_Exp))
   WRITE (UnSum,'(2X,A)')      'Wake diameter calculation model (-): '//trim(Num2LStr(farm%WD(1)%p%Mod_WakeDiam))
   select case ( farm%WD(1)%p%Mod_WakeDiam )
   case (WakeDiamMod_RotDiam) 
      WRITE (UnSum,'(4X,A)')        '( rotor diameter )'
   case (WakeDiamMod_Velocity)  
      WRITE (UnSum,'(4X,A)')        '( velocity based )'
   case (WakeDiamMod_MassFlux) 
      WRITE (UnSum,'(4X,A)')        '( mass-flux based )'
   case (WakeDiamMod_MtmFlux)
      WRITE (UnSum,'(4X,A)')        '( momentum-flux based )'
   end select
   if ( farm%WD(1)%p%Mod_WakeDiam > 1 ) then
      CalWakeDiamStr = trim(Num2LStr(farm%WD(1)%p%C_WakeDiam)) 
   else
      CalWakeDiamStr = '-'
   end if
   WRITE (UnSum,'(2X,A)')      'Calibrated parameter for wake diameter calculation (-): '//CalWakeDiamStr
   WRITE (UnSum,'(2X,A)')      'Spatial filter model for wake meandering (-): '//trim(Num2LStr(farm%AWAE%p%Mod_Meander))
   select case ( farm%AWAE%p%Mod_Meander )
   case (MeanderMod_Uniform) 
      WRITE (UnSum,'(4X,A)')        '( uniform )'
   case (MeanderMod_TruncJinc)  
      WRITE (UnSum,'(4X,A)')        '( truncated jinc )'
   case (MeanderMod_WndwdJinc) 
      WRITE (UnSum,'(4X,A)')        '( windowed jinc )'
   end select
WRITE (UnSum,'(2X,A)')      'Calibrated parameter for wake meandering (-): '//trim(Num2LStr(farm%AWAE%p%C_Meander))
   
   WRITE (UnSum,'(/,A)'   )  'Time Steps'
   WRITE (UnSum,'(2X,A)')      'Component                        Time Step         Subcyles'
   WRITE (UnSum,'(2X,A)')      '  (-)                                (s)             (-)'
   WRITE (UnSum,'(2X,A,F10.4,13X,A)')      'FAST.Farm (glue code)          ',farm%p%dt_low, '1'
   WRITE (UnSum,'(2X,A,F10.4,13X,A)')      'Super Controller               ',farm%p%dt_low, '1'
   WRITE (UnSum,'(2X,A,F10.4,13X,A)')      'FAST Wrapper                   ',farm%p%dt_low, '1 (See table above for OpenFAST.)'
   WRITE (UnSum,'(2X,A,F10.4,13X,A)')      'Wake Dynamics                  ',farm%p%dt_low, '1'
   WRITE (UnSum,'(2X,A,F10.4,13X,A)')      'Ambient Wind and Array Effects ',farm%p%dt_low, '1'
   WRITE (UnSum,'(2X,A,F10.4,13X,A)')      'Low -resolution wind input     ',farm%p%dt_low, '1'
   WRITE (UnSum,'(2X,A,F10.4,12X,I2)')     'High-resolution wind input     ',farm%p%DT_high, farm%p%n_high_low
   WRITE (UnSum,'(2X,A,F10.4,12X,I2,A)')   'Wind visualization output      ',farm%AWAE%p%WrDisSkp1*farm%p%dt_low, farm%AWAE%p%WrDisSkp1, '^-1'
   WRITE (UnSum,'(2X,A,F10.4,13X,A)')      'FAST.Farm output files         ',farm%p%dt_low, '1'

   WRITE (UnSum,'(/,A)'   )  'Requested Channels in FAST.Farm Output Files: '//trim(Num2LStr(farm%p%NumOuts+1))
   WRITE (UnSum,'(2X,A)'  )    'Number     Name         Units'
   WRITE (UnSum,'(2X,A)'  )    '   0       Time          (s)'
   do I=1,farm%p%NumOuts
      WRITE (UnSum,'(2X,I4,7X,A14,A14)'  )  I, farm%p%OutParam(I)%Name,farm%p%OutParam(I)%Units
   end do

   CLOSE(UnSum)

RETURN
END SUBROUTINE Farm_PrintSum
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine initializes the output for the glue code, including writing the header for the primary output file.
SUBROUTINE Farm_InitOutput( farm, ErrStat, ErrMsg )

   IMPLICIT NONE

      ! Passed variables
   type(All_FastFarm_Data),        INTENT(INOUT)        :: farm                                  !< FAST.Farm data
   INTEGER(IntKi),                 INTENT(OUT)          :: ErrStat                               !< Error status
   CHARACTER(*),                   INTENT(OUT)          :: ErrMsg                                !< Error message corresponding to ErrStat


      ! Local variables.

   INTEGER(IntKi)                   :: I, J                                            ! Generic index for DO loops.
   INTEGER(IntKi)                   :: indxLast                                        ! The index of the last value to be written to an array
   INTEGER(IntKi)                   :: indxNext                                        ! The index of the next value to be written to an array
   INTEGER(IntKi)                   :: NumOuts                                         ! number of channels to be written to the output file(s)
   
   if ( (farm%p%NumOuts == 0) ) then ! .or. .not. ( (farm%p%WrTxtOutFile) .or. (farm%p%WrBinOutFile) ) ) then
      return
   end if
   
   ALLOCATE ( farm%m%AllOuts(0:Farm_MaxOutPts) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error allocating memory for the Fast.Farm AllOuts array.'
         RETURN
      ENDIF   
      
      
   farm%m%AllOuts = 0.0_ReKi
#ifdef _OPENMP  
   farm%p%FileDescLines(1)  = 'Predictions were generated on '//CurDate()//' at '//CurTime()//' using '//TRIM(GetVersion(Farm_Ver))//' and with OpenMP'
#else
   farm%p%FileDescLines(1)  = 'Predictions were generated on '//CurDate()//' at '//CurTime()//' using '//TRIM(GetVersion(Farm_Ver))
#endif 
   
   farm%p%FileDescLines(2)  = 'linked with ' //' '//TRIM(GetNVD(NWTC_Ver            ))  ! we'll get the rest of the linked modules in the section below
   farm%p%FileDescLines(3)  = 'Description from the FAST.Farm input file: '//TRIM(farm%p%FTitle)
   
 !......................................................
   ! Open the text output file and print the headers
   !......................................................

  ! IF (farm%p%WrTxtOutFile) THEN

      CALL GetNewUnit( farm%p%UnOu, ErrStat, ErrMsg )
         IF ( ErrStat >= AbortErrLev ) RETURN

      CALL OpenFOutFile ( farm%p%UnOu, TRIM(farm%p%OutFileRoot)//'.out', ErrStat, ErrMsg )
         IF ( ErrStat >= AbortErrLev ) RETURN

         ! Add some file information:

      WRITE (farm%p%UnOu,'(/,A)')  TRIM( farm%p%FileDescLines(1) )
      WRITE (farm%p%UnOu,'(1X,A)') TRIM( farm%p%FileDescLines(2) )
      WRITE (farm%p%UnOu,'()' )    !print a blank line
      WRITE (farm%p%UnOu,'(A)'   ) TRIM( farm%p%FileDescLines(3) )
      WRITE (farm%p%UnOu,'()' )    !print a blank line


         !......................................................
         ! Write the names of the output parameters on one line:
         !......................................................

      CALL WrFileNR ( farm%p%UnOu, farm%p%OutParam(0)%Name )

      DO I=1,farm%p%NumOuts
        CALL WrFileNR ( farm%p%UnOu, farm%p%Delim//farm%p%OutParam(I)%Name )
      ENDDO ! I
!============================================================
! DEBUG OUTPUTS HERE
!
!      DO I = 0,farm%WD(1)%p%NumPlanes-1  ! Loop through all selected output channels
!
!         WRITE( farm%p%UnOu,'(A14)',ADVANCE='NO')   'PPLANEX'//trim(num2lstr(I))
!         WRITE( farm%p%UnOu,'(A14)',ADVANCE='NO')   'PPLANEY'//trim(num2lstr(I))
!         WRITE( farm%p%UnOu,'(A14)',ADVANCE='NO')   'PPLANEZ'//trim(num2lstr(I))  
!
!         WRITE( farm%p%UnOu,'(A14)',ADVANCE='NO' )  'XPLANE'//trim(num2lstr(I))
!
!         WRITE( farm%p%UnOu,'(A14)',ADVANCE='NO' )  'VPLANEX'//trim(num2lstr(I))
!         WRITE( farm%p%UnOu,'(A14)',ADVANCE='NO' )  'VPLANEY'//trim(num2lstr(I))
!         WRITE( farm%p%UnOu,'(A14)',ADVANCE='NO' )  'VPLANEZ'//trim(num2lstr(I))
!
!      ENDDO             ! I - All selected output channels
!      
!      
! END DEBUG OUTPUTS 
!============================================================
      WRITE (farm%p%UnOu,'()')

         !......................................................
         ! Write the units of the output parameters on one line:
         !......................................................

      CALL WrFileNR ( farm%p%UnOu, farm%p%OutParam(0)%Units )

      DO I=1,farm%p%NumOuts
         WRITE (farm%p%UnOu,'(A14)',ADVANCE='NO')  farm%p%Delim//farm%p%OutParam(I)%Units
         !CALL WrFileNR ( farm%p%UnOu, farm%p%Delim//farm%p%OutParam(I)%Units )
      ENDDO ! I

!============================================================
! DEBUG OUTPUTS HERE
!
!      DO I = 0,farm%WD(1)%p%NumPlanes-1  ! Loop through all selected output channels
!
!         WRITE( farm%p%UnOu,'(A14)',ADVANCE='NO')   '      (m)     '
!         WRITE( farm%p%UnOu,'(A14)',ADVANCE='NO')   '      (m)     '
!         WRITE( farm%p%UnOu,'(A14)',ADVANCE='NO')   '      (m)     '
!
!         WRITE( farm%p%UnOu,'(A14)',ADVANCE='NO' )  '      (m)     '
!
!         WRITE( farm%p%UnOu,'(A14)',ADVANCE='NO' )  '    (m/s)     '
!         WRITE( farm%p%UnOu,'(A14)',ADVANCE='NO' )  '    (m/s)     '
!         WRITE( farm%p%UnOu,'(A14)',ADVANCE='NO' )  '    (m/s)     '
!
!         IF ( I < farm%WD(1)%p%NumPlanes-1 ) THEN
!            WRITE( farm%p%UnOu,'(A14)',ADVANCE='NO' )  '      (-)     '
!         END IF
!
!      ENDDO             ! I - All selected output channels
!      
!      
! END DEBUG OUTPUTS 
!============================================================
      
      WRITE (farm%p%UnOu,'()')

  ! END IF

   ! TODO: Add binary
   !......................................................
   ! Allocate data for binary output file
   !......................................................
   !IF (farm%p%WrBinOutFile) THEN
   !
   !      ! calculate the size of the array of outputs we need to store
   !   farm%p%NOutSteps = CEILING ( (farm%p%TMax - farm%p%TStart) / farm%p%DT_low ) + 1
   !
   !   CALL AllocAry( farm%m%AllOutData, farm%p%NumOuts-1, farm%p%NOutSteps, 'AllOutData', ErrStat, ErrMsg )
   !   IF ( ErrStat >= AbortErrLev ) RETURN
   !
   !  ! IF ( OutputFileFmtID == FileFmtID_WithoutTime ) THEN
   !
   !      CALL AllocAry( farm%m%TimeData, 2_IntKi, 'TimeData', ErrStat, ErrMsg )
   !      IF ( ErrStat >= AbortErrLev ) RETURN
   !
   !      farm%m%TimeData(1) = 0.0_DbKi           ! This is the first output time, which we will set later
   !      farm%m%TimeData(2) = farm%p%DT_low      ! This is the (constant) time between subsequent writes to the output file
   !
   !   !ELSE  ! we store the entire time array
   !   !
   !   !   CALL AllocAry( farm%m%TimeData, farm%p%NOutSteps, 'TimeData', ErrStat, ErrMsg )
   !   !   IF ( ErrStat >= AbortErrLev ) RETURN
   !   !
   !   !END IF
   !
   !   farm%m%n_Out = 0  !number of steps actually written to the file
   !
   !END IF



RETURN
END SUBROUTINE Farm_InitOutput   



!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine is called at program termination. It writes any additional output files,
!! deallocates variables for FAST file I/O and closes files.
SUBROUTINE Farm_EndOutput( farm, ErrStat, ErrMsg )
   type(All_FastFarm_Data),  INTENT(INOUT) :: farm                                  !< FAST.Farm data


   INTEGER(IntKi),           INTENT(OUT)   :: ErrStat                   !< Error status
   CHARACTER(*),             INTENT(OUT)   :: ErrMsg                    !< Message associated with errro status

      ! local variables
   CHARACTER(1024)  :: FileDesc                  ! The description of the run, to be written in the binary output file

   !CHARACTER(ChanLenFF):: ChannelNames(farm%p%NumOuts)
   !CHARACTER(ChanLenFF):: ChannelUnits(farm%p%NumOuts)
   !INTEGER(IntKi)  :: I
      ! Initialize some values

   ErrStat = ErrID_None
   ErrMsg  = ''

   !-------------------------------------------------------------------------------------------------
   ! Write the binary output file if requested
   !-------------------------------------------------------------------------------------------------
   ! TODO: The ChannelNames and ChannelUnits need to be length ChanLenFF for Fast.Farm, but the WrBinFAST subroutine needs these to be ChanLen long!
   !IF (farm%p%WrBinOutFile .AND. farm%m%n_Out > 0) THEN
   !
   !   FileDesc = TRIM(farm%p%FileDescLines(1))//' '//TRIM(farm%p%FileDescLines(2))//'; '//TRIM(farm%p%FileDescLines(3))
   !
   !   DO I = 1,farm%p%NumOuts
   !      ChannelNames(I) = farm%p%OutParam(I)%Name
   !      ChannelUnits(I) = farm%p%OutParam(I)%Units
   !   END DO
   !   
   !   CALL WrBinFAST(TRIM(farm%p%OutFileRoot)//'.outb', 2, TRIM(FileDesc), &
   !         ChannelNames, ChannelUnits, farm%m%TimeData(:),farm%m%AllOutData(:,1:farm%m%n_Out), ErrStat, ErrMsg)
   !
   !   IF ( ErrStat /= ErrID_None ) CALL WrScr( TRIM(GetErrStr(ErrStat))//' when writing binary output file: '//TRIM(ErrMsg) )
   !
   !END IF


   !-------------------------------------------------------------------------------------------------
   ! Close the text tabular output file and summary file (if opened)
   !-------------------------------------------------------------------------------------------------
   IF (farm%p%UnOu  > 0) THEN ! I/O unit number for the tabular output file
      CLOSE( farm%p%UnOu )         
      farm%p%UnOu = -1
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Deallocate arrays
   !-------------------------------------------------------------------------------------------------

      ! Output
   !IF ( ALLOCATED(y_FAST%AllOutData                  ) ) DEALLOCATE(y_FAST%AllOutData                  )
   !IF ( ALLOCATED(y_FAST%TimeData                    ) ) DEALLOCATE(y_FAST%TimeData                    )
   !IF ( ALLOCATED(y_FAST%ChannelNames                ) ) DEALLOCATE(y_FAST%ChannelNames                )
   !IF ( ALLOCATED(y_FAST%ChannelUnits                ) ) DEALLOCATE(y_FAST%ChannelUnits                )


END SUBROUTINE Farm_EndOutput

!----------------------------------------------------------------------------------------------------------------------------------
! ROUTINES TO OUTPUT WRITE DATA TO FILE AT EACH REQUSTED TIME STEP
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine determines if it's time to write to the output files, and calls the routine to write to the files
!! with the output data. It should be called after all the output solves for a given time have been completed.
SUBROUTINE WriteFarmOutputToFile( t_global, farm, ErrStat, ErrMsg )
!...............................................................................................................................
   REAL(DbKi),               INTENT(IN   ) :: t_global            !< Current global time
   type(All_FastFarm_Data),  INTENT(INOUT) :: farm                !< FAST.Farm data
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None


   REAL(DbKi)                              :: OutTime             ! Used to determine if output should be generated at this simulation time
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'WriteFarmOutputToFile'
   CHARACTER(200)                          :: Frmt                                      ! A string to hold a format specifier
   CHARACTER(farm%p%TChanLen)              :: TmpStr                                    ! temporary string to print the time output as text 
   CHARACTER(ChanLen)                      :: TmpStr2                                    ! temporary string to print the output as text 
   INTEGER(IntKi)                          :: I, J                                      ! loop counter
   REAL(ReKi)                              :: OutputAry(farm%p%NumOuts)
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   Frmt = '"'//farm%p%Delim//'"'//farm%p%OutFmt      ! format for array elements from individual modules
   
      ! Write time-series channel data
   IF ( farm%p%NumOuts == 0 ) return
   
   IF ( t_global >= farm%p%TStart )  THEN
      
      WRITE( TmpStr, '('//trim(farm%p%OutFmt_t)//')' ) t_global
      CALL WrFileNR( farm%p%UnOu, TmpStr )

            ! Generate fast.farm output file
      
      
      DO I = 1,farm%p%NumOuts  ! Loop through all selected output channels

         OutputAry(I) = farm%p%OutParam(I)%SignM * farm%m%AllOuts( farm%p%OutParam(I)%Indx )
        ! WRITE( TmpStr2, '('//trim(Frmt)//')' ) OutputAry(I)
        ! WRITE (farm%p%UnOu,'(A14)',ADVANCE='NO') TmpStr2
         
        ! CALL WrFileNR( farm%p%UnOu, TmpStr2 )
      
     
      ENDDO             ! I - All selected output channels
        ! write the individual module output (convert to SiKi if necessary, so that we don't need to print so many digits in the exponent)
      CALL WrNumAryFileNR ( farm%p%UnOu, REAL(OutputAry,SiKi), Frmt, ErrStat, ErrMsg ) 
!============================================================
! DEBUG OUTPUTS HERE
!
!      DO I = 0,farm%WD(1)%p%NumPlanes-1  ! Loop through all selected output channels
!
!         DO J = 1,3
!            WRITE( TmpStr2, '('//trim(farm%p%OutFmt)//')' )  farm%WD(1)%y%p_plane(J,I)
!            CALL WrFileNR( farm%p%UnOu, TmpStr2 )
!         ENDDO
!
!         WRITE( TmpStr2, '('//trim(farm%p%OutFmt)//')' )  farm%WD(1)%xd%x_plane(I)
!         CALL WrFileNR( farm%p%UnOu, TmpStr2 )
!
!         DO J = 1,3
!            WRITE( TmpStr2, '('//trim(farm%p%OutFmt)//')' )  farm%AWAE%y%V_plane(J,I,1)
!            CALL WrFileNR( farm%p%UnOu, TmpStr2 )
!         ENDDO
!
!      ENDDO             ! I - All selected output channels
!      
!      
! END DEBUG OUTPUTS 
!============================================================
         ! write a new line (advance to the next line)
      WRITE (farm%p%UnOu,'()')

      !IF (farm%p%WrBinOutFile) THEN
      !
      !      ! Write data to array for binary output file
      !
      !   IF ( farm%m%n_Out == farm%p%NOutSteps ) THEN
      !      CALL ProgWarn( 'Not all data could be written to the binary output file.' )
      !      !this really would only happen if we have an error somewhere else, right?
      !      !otherwise, we could allocate a new, larger array and move existing data
      !   ELSE
      !      farm%m%n_Out = farm%m%n_Out + 1
      !
      !         ! store time data
      !      IF ( farm%m%n_Out == 1_IntKi ) THEN !.OR. OutputFileFmtID == FileFmtID_WithTime ) THEN
      !         farm%m%TimeData(farm%m%n_Out) = t_global   ! Time associated with these outputs
      !      END IF
      !
      !         ! store individual module data
      !      farm%m%AllOutData(:, farm%m%n_Out) = OutputAry
      !   
      !   END IF      
      !
      !END IF  
   ENDIF
END SUBROUTINE WriteFarmOutputToFile  



end module FAST_Farm_IO
