module FAST_Farm_IO
   
   USE NWTC_Library
   USE VersionInfo
   USE FAST_Farm_Types
   USE FAST_Farm_IO_Params
   
   IMPLICIT NONE
   
   TYPE(ProgDesc), PARAMETER      :: Farm_Ver      = ProgDesc( 'FAST.Farm', '', '' ) !< module date/version information 

   integer, parameter :: maxOutputPoints = 9
   integer, parameter :: maxOutputPlanes = 999     ! Allow up to 99 outpt planes


      contains

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
   
!FIXME: add summary info about WAT

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

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine reads in the primary FAST.Farm input file, does some validation, and places the values it reads in the
!!   parameter structure (p). It prints to an echo file if requested.
SUBROUTINE Farm_ReadPrimaryFile( InputFile, p, WD_InitInp, AWAE_InitInp, OutList, ErrStat, ErrMsg )
   TYPE(Farm_ParameterType),       INTENT(INOUT) :: p                               !< The parameter data for the FAST (glue-code) simulation
   CHARACTER(*),                   INTENT(IN   ) :: InputFile                       !< Name of the file containing the primary input data
   TYPE(WD_InputFileType),         INTENT(  OUT) :: WD_InitInp                      !< input-file data for WakeDynamics module
   TYPE(AWAE_InputFileType),       INTENT(  OUT) :: AWAE_InitInp                    !< input-file data for AWAE module
   CHARACTER(ChanLen),             INTENT(  OUT) :: OutList(:)                      !< list of user-requested output channels
   INTEGER(IntKi),                 INTENT(  OUT) :: ErrStat                         !< Error status
   CHARACTER(*),                   INTENT(  OUT) :: ErrMsg                          !< Error message

      ! Local variables:
   REAL(DbKi)                    :: TmpTime                                   ! temporary variable to read SttsTime and ChkptTime before converting to #steps based on DT_low
   INTEGER(IntKi)                :: I                                         ! loop counter
   INTEGER(IntKi)                :: UnIn                                      ! Unit number for reading file
   INTEGER(IntKi)                :: UnEc                                      ! I/O unit for echo file. If > 0, file is open for writing.

   INTEGER(IntKi)                :: IOS                                       ! Temporary Error status
   INTEGER(IntKi)                :: OutFileFmt                                ! An integer that indicates what kind of tabular output should be generated (1=text, 2=binary, 3=both)
   INTEGER(IntKi)                :: NLinTimes                                 ! An integer that indicates how many times to linearize
   LOGICAL                       :: Echo                                      ! Determines if an echo file should be written
   LOGICAL                       :: TabDelim                                  ! Determines if text output should be delimited by tabs (true) or space (false)
   CHARACTER(1024)               :: PriPath                                   ! Path name of the primary file
   character(1024)               :: sDummy ! Dummy string

   CHARACTER(10)                 :: AbortLevel                                ! String that indicates which error level should be used to abort the program: WARNING, SEVERE, or FATAL
   CHARACTER(30)                 :: Line                                      ! string for default entry in input file

   INTEGER(IntKi)                :: ErrStat2                                  ! Temporary Error status
   CHARACTER(ErrMsgLen)          :: ErrMsg2                                   ! Temporary Error message
   CHARACTER(*),   PARAMETER     :: RoutineName = 'Farm_ReadPrimaryFile'
   Real(ReKi)                    :: DefaultReVal ! Default real value
   real(ReKi)                    :: TmpRAry5(5)    ! Temporary array for reading in array of 5

      ! Initialize some variables:
   UnEc = -1
   Echo = .FALSE.                        ! Don't echo until we've read the "Echo" flag
   CALL GetPath( InputFile, PriPath )    ! Input files will be relative to the path where the primary input file is located.

      ! Get an available unit number and open input file
   CALL GetNewUnit( UnIn, ErrStat, ErrMsg );  IF ( ErrStat >= AbortErrLev ) RETURN
   CALL OpenFInpFile ( UnIn, InputFile, ErrStat2, ErrMsg2 );   if (Failed()) return

   ! Read the lines up/including to the "Echo" simulation control variable
   ! If echo is FALSE, don't write these lines to the echo file.
   ! If Echo is TRUE, rewind and write on the second try.

   I = 1 !set the number of times we've read the file
   DO
   !-------------------------- HEADER ---------------------------------------------
      CALL ReadCom( UnIn, InputFile, 'File header: FAST.Farm Version (line 1)', ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
      CALL ReadStr( UnIn, InputFile, p%FTitle, 'FTitle', 'File Header: File Description (line 2)', ErrStat2, ErrMsg2, UnEc ); if (Failed()) return

   !---------------------- SIMULATION CONTROL --------------------------------------
      CALL ReadCom( UnIn, InputFile, 'Section Header: Simulation Control', ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
      CALL ReadVar( UnIn, InputFile, Echo, "Echo", "Echo input data to <RootName>.ech (flag)", ErrStat2, ErrMsg2, UnEc); if (Failed()) return

      IF (.NOT. Echo .OR. I > 1) EXIT !exit this loop

         ! Otherwise, open the echo file, then rewind the input file and echo everything we've read
      I = I + 1         ! make sure we do this only once (increment counter that says how many times we've read this file)
      CALL OpenEcho ( UnEc, TRIM(p%OutFileRoot)//'.ech', ErrStat2, ErrMsg2, Farm_Ver ); if (Failed()) return
      IF ( UnEc > 0 )  WRITE (UnEc,'(/,A,/)')  'Data from '//TRIM(Farm_Ver%Name)//' primary input file "'//TRIM( InputFile )//'":'
      REWIND( UnIn, IOSTAT=ErrStat2 )
         IF (ErrStat2 /= 0_IntKi ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error rewinding file "'//TRIM(InputFile)//'".',ErrStat,ErrMsg,RoutineName)
            call cleanup()
            RETURN
         END IF
   END DO

   CALL WrScr( ' Heading of the '//TRIM(Farm_Ver%Name)//' input file: ' )
   CALL WrScr( '   '//TRIM( p%FTitle ) )

   ! AbortLevel - Error level when simulation should abort:
   CALL ReadVar( UnIn, InputFile, AbortLevel, "AbortLevel", "Error level when simulation should abort (string)", &
                        ErrStat2, ErrMsg2, UnEc); if (Failed()) return

      ! Let's set the abort level here.... knowing that everything before this aborted only on FATAL errors!
      CALL Conv2UC( AbortLevel ) !convert to upper case
      SELECT CASE( TRIM(AbortLevel) )
         CASE ( "WARNING" )
            AbortErrLev = ErrID_Warn
         CASE ( "SEVERE" )
            AbortErrLev = ErrID_Severe
         CASE ( "FATAL" )
            AbortErrLev = ErrID_Fatal
         CASE DEFAULT
            CALL SetErrStat( ErrID_Fatal, 'Invalid AbortLevel specified in FAST.Farm input file. '// &
                             'Valid entries are "WARNING", "SEVERE", or "FATAL".',ErrStat,ErrMsg,RoutineName)
            call cleanup()
            RETURN
      END SELECT

   CALL ReadVar( UnIn, InputFile, p%TMax, "TMax", "Total run time (s)", ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   CALL ReadVar( UnIn, InputFile, AWAE_InitInp%Mod_AmbWind, "Mod_AmbWind", "Ambient wind model (-) (switch) {1: high-fidelity precursor in VTK format, 2: one InflowWind module, 3: multiple InflowWind modules}", ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   CALL ReadVar( UnIn, InputFile, p%WaveFieldMod, "Mod_WaveField", "Wave field handling (-) (switch) {1: use individual HydroDyn inputs without adjustment, 2: adjust wave phases based on turbine offsets from farm origin}", ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   CALL ReadVar( UnIn, InputFile, p%MooringMod, "Mod_SharedMooring", "Array-level mooring handling (-) (switch) {0: none; 3: array-level MoorDyn model}", ErrStat2, ErrMsg2, UnEc); if (Failed()) return

   !---------------------- SHARED MOORING SYSTEM ------------------------------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: SHARED MOORING SYSTEM', ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   CALL ReadVar( UnIn, InputFile, p%MD_FileName, "MD_FileName", "Name/location of the dynamic library {.dll [Windows] or .so [Linux]} containing the Super Controller algorithms (quoated string)", ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   IF ( PathIsRelative( p%MD_FileName ) ) p%MD_FileName = TRIM(PriPath)//TRIM(p%MD_FileName)
   CALL ReadVar( UnIn, InputFile, p%DT_mooring, "DT_Mooring", "Time step for farm-levem mooring coupling with each turbine [used only when Mod_SharedMooring > 0] (s) [>0.0]", ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   CALL ReadVar( UnIn, InputFile, p%WrMooringVis, "MooringVis","Write shared mooring visualization, at DT_Mooring timestep (-) [only used for Mod_SharedMooring=3]", ErrStat2, ErrMsg2, UnEc); if (Failed()) return

   !---------------------- AMBIENT WIND: PRECURSOR IN VTK FORMAT ---------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Ambient Wind: Precursor in VTK Format', ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   CALL ReadVar( UnIn, InputFile, p%DT_low, "DT_Low-VTK", "Time step for low-resolution wind data input files; will be used as the global FAST.Farm time step (s) [>0.0]", ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   CALL ReadVar( UnIn, InputFile, p%DT_high, "DT_High-VTK", "Time step for high-resolution wind data input files (s) [>0.0]", ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   CALL ReadVar( UnIn, InputFile, p%WindFilePath, "WindFilePath", "Path name of wind data files from ABLSolver precursor (string)", ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   IF ( PathIsRelative( p%WindFilePath ) ) p%WindFilePath = TRIM(PriPath)//TRIM(p%WindFilePath)
   CALL ReadVar( UnIn, InputFile, AWAE_InitInp%ChkWndFiles, "ChkWndFiles", "Check all the ambient wind files for data consistency? (flag)", ErrStat2, ErrMsg2, UnEc); if (Failed()) return

   !---------------------- AMBIENT WIND: INFLOWWIND MODULE ---------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Ambient Wind: InflowWind Module', ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   CALL ReadVar( UnIn, InputFile, AWAE_InitInp%DT_low, "DT_Low", "Time step for low-resolution wind data input files; will be used as the global FAST.Farm time step (s) [>0.0]", ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   CALL ReadVar( UnIn, InputFile, AWAE_InitInp%DT_high, "DT_High", "Time step for high-resolution wind data input files (s) [>0.0]", ErrStat2, ErrMsg2, UnEc); if (Failed()) return

   ! Ensure consistency between AWAE_Inputs and FAST.Farm time steps
   if ( AWAE_InitInp%Mod_AmbWind == 1) AWAE_InitInp%DT_high = p%DT_high
   if ( AWAE_InitInp%Mod_AmbWind == 1) AWAE_InitInp%DT_low  = p%DT_low
   if ( AWAE_InitInp%Mod_AmbWind > 1 ) p%DT_low = AWAE_InitInp%DT_low
   if ( AWAE_InitInp%Mod_AmbWind > 1 ) p%DT_high = AWAE_InitInp%DT_high

   ! low res
   CALL ReadVar( UnIn, InputFile, AWAE_InitInp%nX_Low,  "nX_Low",  "Number of low-resolution spatial nodes in X direction for wind data interpolation (-) [>=2]",   ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   CALL ReadVar( UnIn, InputFile, AWAE_InitInp%nY_Low,  "nY_Low",  "Number of low-resolution spatial nodes in Y direction for wind data interpolation (-) [>=2]",   ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   CALL ReadVar( UnIn, InputFile, AWAE_InitInp%nZ_Low,  "nZ_Low",  "Number of low-resolution spatial nodes in Z direction for wind data interpolation (-) [>=2]",   ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   CALL ReadVar( UnIn, InputFile, AWAE_InitInp%X0_Low,  "X0_Low",  "Origin of low-resolution spatial nodes in X direction for wind data interpolation (m)",         ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   CALL ReadVar( UnIn, InputFile, AWAE_InitInp%Y0_Low,  "Y0_Low",  "Origin of low-resolution spatial nodes in Y direction for wind data interpolation (m)",         ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   CALL ReadVar( UnIn, InputFile, AWAE_InitInp%Z0_Low,  "Z0_Low",  "Origin of low-resolution spatial nodes in Z direction for wind data interpolation (m)",         ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   CALL ReadVar( UnIn, InputFile, AWAE_InitInp%dX_Low,  "dX_Low",  "Spacing of low-resolution spatial nodes in X direction for wind data interpolation (m) [>0.0]", ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   CALL ReadVar( UnIn, InputFile, AWAE_InitInp%dY_Low,  "dY_Low",  "Spacing of low-resolution spatial nodes in Y direction for wind data interpolation (m) [>0.0]", ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   CALL ReadVar( UnIn, InputFile, AWAE_InitInp%dZ_Low,  "dZ_Low",  "Spacing of low-resolution spatial nodes in Z direction for wind data interpolation (m) [>0.0]", ErrStat2, ErrMsg2, UnEc); if (Failed()) return

   ! high res
   CALL ReadVar( UnIn, InputFile, AWAE_InitInp%nX_High, "nX_High", "Number of high-resolution spatial nodes in X direction for wind data interpolation (-) [>=2]",  ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   CALL ReadVar( UnIn, InputFile, AWAE_InitInp%nY_High, "nY_High", "Number of high-resolution spatial nodes in Y direction for wind data interpolation (-) [>=2]",  ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   CALL ReadVar( UnIn, InputFile, AWAE_InitInp%nZ_High, "nZ_High", "Number of high-resolution spatial nodes in Z direction for wind data interpolation (-) [>=2]",  ErrStat2, ErrMsg2, UnEc); if (Failed()) return

   ! inflow file
   CALL ReadVar( UnIn, InputFile, AWAE_InitInp%InflowFile, "InflowFile", "Name of file containing InflowWind module input parameters (quoted string)", ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   IF ( PathIsRelative( AWAE_InitInp%InflowFile ) ) AWAE_InitInp%InflowFile = TRIM(PriPath)//TRIM(AWAE_InitInp%InflowFile)
   if ( AWAE_InitInp%Mod_AmbWind > 1 ) p%WindFilePath = AWAE_InitInp%InflowFile  ! For the summary file

   !---------------------- WIND TURBINES ---------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Wind Turbines', ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   CALL ReadVar( UnIn, InputFile, p%NumTurbines, "NumTurbines", "Number of wind turbines (-) [>=1]", ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   CALL ReadCom( UnIn, InputFile, 'Section Header: WT column names', ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   CALL ReadCom( UnIn, InputFile, 'Section Header: WT column units', ErrStat2, ErrMsg2, UnEc ); if (Failed()) return

   call AllocAry( p%WT_Position, 3, p%NumTurbines, 'WT_Position',   ErrStat2, ErrMsg2);  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName); if (Failed()) return
   call AllocAry( p%WT_FASTInFile,  p%NumTurbines, 'WT_FASTInFile', ErrStat2, ErrMsg2);  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName); if (Failed()) return
   call AllocAry( AWAE_InitInp%WT_Position, 3, p%NumTurbines, 'AWAE_InitInp%WT_Position', ErrStat2, ErrMsg2);  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName); if (Failed()) return

   if ( AWAE_InitInp%Mod_AmbWind > 1 ) then     ! Using InflowWind
      call AllocAry(AWAE_InitInp%X0_high, p%NumTurbines, 'AWAE_InitInp%X0_high', ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(AWAE_InitInp%Y0_high, p%NumTurbines, 'AWAE_InitInp%Y0_high', ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(AWAE_InitInp%Z0_high, p%NumTurbines, 'AWAE_InitInp%Z0_high', ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(AWAE_InitInp%dX_high, p%NumTurbines, 'AWAE_InitInp%dX_high', ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(AWAE_InitInp%dY_high, p%NumTurbines, 'AWAE_InitInp%dY_high', ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(AWAE_InitInp%dZ_high, p%NumTurbines, 'AWAE_InitInp%dZ_high', ErrStat2, ErrMsg2); if (Failed()) return
   end if

      ! WT_Position (WT_X, WT_Y, WT_Z) and WT_FASTInFile
   do i=1,p%NumTurbines
      if ( AWAE_InitInp%Mod_AmbWind == 1 ) then
         READ (UnIn, *, IOSTAT=IOS) p%WT_Position(:,i), p%WT_FASTInFile(i)
      else
         READ (UnIn, *, IOSTAT=IOS) p%WT_Position(:,i), p%WT_FASTInFile(i), AWAE_InitInp%X0_high(i), AWAE_InitInp%Y0_high(i), AWAE_InitInp%Z0_high(i), AWAE_InitInp%dX_high(i), AWAE_InitInp%dY_high(i), AWAE_InitInp%dZ_high(i)
      end if
      AWAE_InitInp%WT_Position(:,i) = p%WT_Position(:,i)
      CALL CheckIOS ( IOS, InputFile, 'Wind Turbine Columns', NumType, ErrStat2, ErrMsg2 ); if (Failed()) return
      IF ( UnEc > 0 ) THEN
         if ( AWAE_InitInp%Mod_AmbWind == 1 ) then
            WRITE( UnEc, "(3(ES11.4e2,2X),'""',A,'""',T50,' - WT(',I5,')')" ) p%WT_Position(:,i), TRIM( p%WT_FASTInFile(i) ), I
         else
            WRITE( UnEc, "(3(ES11.4e2,2X),'""',A,'""',T50,6(ES11.4e2,2X),' - WT(',I5,')')" ) p%WT_Position(:,i), TRIM( p%WT_FASTInFile(i) ), AWAE_InitInp%X0_high(i), AWAE_InitInp%Y0_high(i), AWAE_InitInp%Z0_high(i), AWAE_InitInp%dX_high(i), AWAE_InitInp%dY_high(i), AWAE_InitInp%dZ_high(i), I
         end if

      END IF
      IF ( PathIsRelative( p%WT_FASTInFile(i) ) ) p%WT_FASTInFile(i) = TRIM(PriPath)//TRIM(p%WT_FASTInFile(i))
   end do


   !---------------------- WAKE DYNAMICS ---------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Wake Dynamics', ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   CALL ReadVar( UnIn, InputFile, WD_InitInp%Mod_Wake, "Mod_Wake",  "Wake model", ErrStat2, ErrMsg2, UnEc); if(failed()) return
   CALL ReadVar( UnIn, InputFile, p%RotorDiamRef     , "RotorDiamRef", "Reference turbine rotor diameter for wake calculations (m) [>0.0]", ErrStat2, ErrMsg2, UnEc); if(failed()) return
   CALL ReadVar( UnIn, InputFile, WD_InitInp%dr      , "dr"      ,  "Radial increment of radial finite-difference grid (m) [>0.0]", ErrStat2, ErrMsg2, UnEc); if(failed()) return
   CALL ReadVar( UnIn, InputFile, WD_InitInp%NumRadii, "NumRadii",  "Number of radii in the radial finite-difference grid (-) [>=2]", ErrStat2, ErrMsg2, UnEc); if(failed()) return
   CALL ReadVar( UnIn, InputFile, WD_InitInp%NumPlanes,"NumPlanes", "Number of wake planes (-) [>=2]", ErrStat2, ErrMsg2, UnEc); if(failed()) return

   ! f_c - Cut-off (corner) frequency of the low-pass time-filter for the wake advection, deflection, and meandering model (Hz) [>0.0] or DEFAULT [DEFAULT=0.0007]:
   DefaultReVal = 12.5_ReKi/(p%RotorDiamRef/2._ReKi) ! Eq. (32) of https://doi.org/10.1002/we.2785, with U=10, a=1/3
   CALL ReadVarWDefault( UnIn, InputFile, WD_InitInp%f_c, "f_c", &
      "Cut-off (corner) frequency of the low-pass time-filter for the wake advection, deflection, and meandering model (Hz) [>0.0] or DEFAULT [DEFAULT=0.0007]", &
      DefaultReVal, ErrStat2, ErrMsg2, UnEc); if (Failed()) return

   ! C_HWkDfl_O
   CALL ReadVarWDefault( UnIn, InputFile, WD_InitInp%C_HWkDfl_O, "C_HWkDfl_O", &
      "Calibrated parameter in the correction for wake deflection defining the horizontal offset at the rotor (m) or DEFAULT [DEFAULT=0.0]", &
      0.0_ReKi, ErrStat2, ErrMsg2, UnEc); if (Failed()) return

   ! C_HWkDfl_OY
   if (WD_InitInp%Mod_Wake == Mod_Wake_Curl) then
      DefaultReVal = 0.0_ReKi
   else
      DefaultReVal = 0.3_ReKi
   endif
   CALL ReadVarWDefault( UnIn, InputFile, WD_InitInp%C_HWkDfl_OY, "C_HWkDfl_OY", &
      "Calibrated parameter in the correction for wake deflection defining the horizontal offset at the rotor scaled with yaw error (m/deg) or DEFAULT [DEFAULT=0.3]", &
      DefaultReVal, ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   WD_InitInp%C_HWkDfl_OY = WD_InitInp%C_HWkDfl_OY/D2R !immediately convert to m/radians instead of m/degrees

   ! C_HWkDfl_x
   CALL ReadVarWDefault( UnIn, InputFile, WD_InitInp%C_HWkDfl_x, "C_HWkDfl_x", &
      "Calibrated parameter in the correction for wake deflection defining the horizontal offset scaled with downstream distance (-) or DEFAULT [DEFAULT=0.0]", &
      0.0_ReKi, ErrStat2, ErrMsg2, UnEc); if (Failed()) return

   ! C_HWkDfl_xY
   if (WD_InitInp%Mod_Wake == Mod_Wake_Curl) then
      DefaultReVal = 0.0_ReKi
   else
      DefaultReVal = -0.004_ReKi
   endif
   CALL ReadVarWDefault( UnIn, InputFile, WD_InitInp%C_HWkDfl_xY, "C_HWkDfl_xY", &
      "Calibrated parameter in the correction for wake deflection defining the horizontal offset scaled with downstream distance and yaw error (1/deg) or DEFAULT [DEFAULT=-0.004]", &
      DefaultReVal, ErrStat2, ErrMsg2, UnEc); if (Failed()) return
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if
   WD_InitInp%C_HWkDfl_xY = WD_InitInp%C_HWkDfl_xY/D2R !immediately convert to 1/radians instead of 1/degrees


   ! C_NearWake
   CALL ReadVarWDefault( UnIn, InputFile, WD_InitInp%C_NearWake, "C_NearWake", &
      "Calibrated parameter for the near-wake correction (-) [>1.0] or DEFAULT [DEFAULT=1.8]", &
      1.8_ReKi, ErrStat2, ErrMsg2, UnEc); if (Failed()) return

   ! k_vAmb - Calibrated parameters for the influence of the shear layer in the eddy viscosity (set of 5 parameters: k, FMin, DMin, DMax, Exp) (-) [>=0.0, >=0.0 and <=1.0, >=0.0, >DMin, >0.0] or DEFAULT [DEFAULT=0.05, 1.0, 0.0, 1.0, 0.01]
   call ReadAryWDefault( UnIn, InputFile, TmpRAry5,  5, "k_vAmb",  &
      "Calibrated parameters for the influence of the shear layer in the eddy viscosity (set of 5 parameters: k, FMin, DMin, DMax, Exp) (-) [>=0.0, >=0.0 and <=1.0, >=0.0, >DMin, >0.0] or DEFAULT [DEFAULT=0.05, 1.0, 0.0, 1.0, 0.01]", &
      (/ 0.05_ReKi, 1.0_ReKi, 0.0_ReKi, 1.0_ReKi, 0.01_ReKi /), ErrStat2, ErrMsg2, UnEc); if (Failed()) return
      WD_InitInp%k_vAmb      = TmpRAry5(1)   ! Calibrated parameter for the influence of ambient turbulence in the eddy viscosity (-) [>=0.0] or DEFAULT [DEFAULT=0.05]
      WD_InitInp%C_vAmb_FMin = TmpRAry5(2)   ! Calibrated parameter in the eddy viscosity filter function for ambient turbulence defining the value in the minimum region (-) [>=0.0 and <=1.0] or DEFAULT [DEFAULT=1.0]
      WD_InitInp%C_vAmb_DMin = TmpRAry5(3)   ! Calibrated parameter in the eddy viscosity filter function for ambient turbulence defining the transitional diameter fraction between the minimum and exponential regions (-) [>=0.0] or DEFAULT [DEFAULT=0.0]
      WD_InitInp%C_vAmb_DMax = TmpRAry5(4)   ! Calibrated parameter in the eddy viscosity filter function for ambient turbulence defining the transitional diameter fraction between the exponential and maximum regions (-) [> C_vAmb_DMin  ] or DEFAULT [DEFAULT=1.0]
      WD_InitInp%C_vAmb_Exp  = TmpRAry5(5)   ! Calibrated parameter in the eddy viscosity filter function for ambient turbulence defining the exponent in the exponential region (-) [> 0.0] or DEFAULT [DEFAULT=0.01]

   ! k_vShr - Calibrated parameters for the influence of the ambient turbulence in the eddy viscosity (set of 5 parameters: k, FMin, DMin, DMax, Exp) (-) [>=0.0, >=0.0 and <=1.0, >=0.0, >DMin, >0.0] or DEFAULT [DEFAULT=0.016, 0.2, 3.0, 25.0, 0.1]
   call ReadAryWDefault( UnIn, InputFile, TmpRAry5,  5, "k_vShr",  &
      "Calibrated parameters for the influence of the ambient turbulence in the eddy viscosity (set of 5 parameters: k, FMin, DMin, DMax, Exp) (-) [>=0.0, >=0.0 and <=1.0, >=0.0, >DMin, >0.0] or DEFAULT [DEFAULT=0.016, 0.2, 3.0, 25.0, 0.1]", &
      (/ 0.016_ReKi, 0.2_ReKi, 3.0_ReKi, 25.0_ReKi, 0.1_ReKi /), ErrStat2, ErrMsg2, UnEc); if (Failed()) return
      WD_InitInp%k_vShr      = TmpRAry5(1)   ! Calibrated parameter for the influence of the shear layer in the eddy viscosity (-) [>=0.0] or DEFAULT [DEFAULT=0.016]
      WD_InitInp%C_vShr_FMin = TmpRAry5(2)   ! Calibrated parameter in the eddy viscosity filter function for the shear layer defining the value in the minimum region (-) [>=0.0 and <=1.0] or DEFAULT [DEFAULT=0.2]
      WD_InitInp%C_vShr_DMin = TmpRAry5(3)   ! Calibrated parameter in the eddy viscosity filter function for the shear layer defining the transitional diameter fraction between the minimum and exponential regions (-) [>=0.0] or DEFAULT [DEFAULT=3.0]
      WD_InitInp%C_vShr_DMax = TmpRAry5(4)   ! Calibrated parameter in the eddy viscosity filter function for the shear layer defining the transitional diameter fraction between the exponential and maximum regions (-) [> C_vShr_DMin] or DEFAULT [DEFAULT=25.0]
      WD_InitInp%C_vShr_Exp  = TmpRAry5(5)   ! Calibrated parameter in the eddy viscosity filter function for the shear layer defining the exponent in the exponential region (-) [> 0.0] or DEFAULT [DEFAULT=0.1]

   ! Mod_WakeDiam - Wake diameter calculation model (-) (switch) {1: rotor diameter, 2: velocity-based, 3: mass-flux based, 4: momentum-flux based} or DEFAULT [DEFAULT=1]:
   CALL ReadVarWDefault( UnIn, InputFile, WD_InitInp%Mod_WakeDiam, "Mod_WakeDiam", &
      "Wake diameter calculation model (-) (switch) {1: rotor diameter, 2: velocity-based, 3: mass-flux based, 4: momentum-flux based} or DEFAULT [DEFAULT=1]", &
      WakeDiamMod_RotDiam, ErrStat2, ErrMsg2, UnEc); if (Failed()) return

   ! C_WakeDiam
   CALL ReadVarWDefault( UnIn, InputFile, WD_InitInp%C_WakeDiam, "C_WakeDiam", &
      "Calibrated parameter for wake diameter calculation (-) [>0.0 and <1.0] or DEFAULT [DEFAULT=0.95] [unused for Mod_WakeDiam=1]", &
      0.95_ReKi, ErrStat2, ErrMsg2, UnEc); if (Failed()) return

   ! Mod_Meander - Spatial filter model for wake meandering (-) (switch) {1: uniform, 2: truncated jinc, 3: windowed jinc} or DEFAULT [DEFAULT=3]:
   CALL ReadVarWDefault( UnIn, InputFile, AWAE_InitInp%Mod_Meander, "Mod_Meander", &
      "Spatial filter model for wake meandering (-) (switch) {1: uniform, 2: truncated jinc, 3: windowed jinc} or DEFAULT [DEFAULT=3]", &
      MeanderMod_WndwdJinc, ErrStat2, ErrMsg2, UnEc); if (Failed()) return

   ! C_Meander
   CALL ReadVarWDefault( UnIn, InputFile, AWAE_InitInp%C_Meander, "C_Meander", &
      "Calibrated parameter for wake meandering (-) [>=1.0] or DEFAULT [DEFAULT=1.9]", &
      1.9_ReKi, ErrStat2, ErrMsg2, UnEc); if (Failed()) return

   !----------------------- CURL WAKE PARAMETERS ------------------------------------------
   CALL ReadCom        ( UnIn, InputFile, "Section Header: Curl wake parameters", ErrStat2, ErrMsg2, UnEc ); if(failed()) return
   CALL ReadVarWDefault( UnIn, InputFile, WD_InitInp%Swirl        ,    "Swirl", "Swirl switch", .True., ErrStat2, ErrMsg2, UnEc); if(failed()) return
   CALL ReadVarWDefault( UnIn, InputFile, WD_InitInp%k_VortexDecay,    "k_VortexDecay", "Vortex decay constant", 0.0, ErrStat2, ErrMsg2, UnEc); if(failed()) return
   CALL ReadVarWDefault( UnIn, InputFile, WD_InitInp%NumVortices,      "NumVortices", "Number of vortices in the curled wake", 100, ErrStat2, ErrMsg2, UnEc); if(failed()) return
   CALL ReadVarWDefault( UnIn, InputFile, WD_InitInp%sigma_D,          "sigma_D", "Gaussian vortex width", 0.2, ErrStat2, ErrMsg2, UnEc); if(failed()) return
   CALL ReadVarWDefault( UnIn, InputFile, WD_InitInp%FilterInit,       "FilterInit", "Filter Init", 1 , ErrStat2, ErrMsg2, UnEc); if(failed()) return
   CALL ReadVarWDefault( UnIn, InputFile, WD_InitInp%k_vCurl,          "k_vCurl",    "Eddy viscosity for curl", 2.0 , ErrStat2, ErrMsg2, UnEc); if(failed()) return
   CALL ReadVarWDefault( UnIn, InputFile, AWAE_InitInp%Mod_Projection, "Mod_Projection", "Mod_Projection", -1 , ErrStat2, ErrMsg2, UnEc); if(failed()) return
   if (AWAE_InitInp%Mod_Projection==-1) then
      ! -1 means the user selected "default"
      if (WD_InitInp%Mod_Wake==Mod_Wake_Curl) then
           AWAE_InitInp%Mod_Projection=2
      else
           AWAE_InitInp%Mod_Projection=1
      endif
   endif
   !----------------------- WAKE-ADDED TURBULENCE ------------------------------------------
   ! Read WAT variables
   CALL ReadCom( UnIn, InputFile, 'Section Header: Wake-added turbulence', ErrStat2, ErrMsg2, UnEc ); if(failed()) return
   CALL ReadVar( UnIn, InputFile, p%WAT, "WAT", "Switch between wake-added turbulence box options {0: no wake added turbulence, 1: predefined turbulence box, 2: user defined turbulence box}", ErrStat2, ErrMsg2, UnEc); if(failed()) return
   CALL ReadVar( UnIn, InputFile, p%WAT_BoxFile, 'WAT_BoxFile', "Filepath to the file containing the u-component of the turbulence box (either predefined or user-defined) (quoted string)", ErrStat2, ErrMsg2, UnEc ); if(failed()) return
   call ReadAry( UnIn, InputFile, p%WAT_NxNyNz, 3, "WAT_NxNyNz", "Number of points in the x, y, and z directions of the WAT_BoxFile [used only if WAT=2] (m)", ErrStat2, ErrMsg2, UnEc ); if(failed()) return
   call ReadAry( UnIn, InputFile, p%WAT_DxDyDz, 3, "WAT_DxDyDz", "Distance (in meters) between points in the x, y, and z directions of the WAT_BoxFile [used only if WAT=2] (m)", ErrStat2, ErrMsg2, UnEc ); if(failed()) return
   call ReadVarWDefault( UnIn, InputFile, p%WAT_ScaleBox, "WAT_ScaleBox",   "Flag to scale the input turbulence box to zero mean and unit standard deviation at every node",  .False., ErrStat2, ErrMsg2, UnEc); if(failed()) return
   call ReadAryWDefault( UnIn, InputFile, TmpRAry5,  5, "WAT_k_Def",  &
         "Calibrated parameters for the influence of the maximum wake deficit on wake-added turbulence (set of 5 parameters: k_Def , DMin, DMax, FMin, Exp) (-) [>=0.0, >=0.0, >DMin, >=0.0 and <=1.0, >=0.0] or DEFAULT [DEFAULT=[0.6, 0.0, 0.0,  2.0, 1.0 ]]", &
         (/0.6_ReKi, 0.0_ReKi, 0.0_ReKi, 2.0_ReKi, 1.00_ReKi/), ErrStat2, ErrMsg2, UnEc); if(failed()) return
      WD_InitInp%WAT_k_Def_k_c  = TmpRAry5(1)
      WD_InitInp%WAT_k_Def_FMin = TmpRAry5(2)
      WD_InitInp%WAT_k_Def_DMin = TmpRAry5(3)
      WD_InitInp%WAT_k_Def_DMax = TmpRAry5(4)
      WD_InitInp%WAT_k_Def_Exp  = TmpRAry5(5)
   call ReadAryWDefault( UnIn, InputFile, TmpRAry5, 5, "WAT_k_Grad", &
         "Calibrated parameters for the influence of the radial velocity gradient of the wake deficit on wake-added turbulence (set of 5 parameters: k_Grad, DMin, DMax, FMin, Exp) (-) [>=0.0, >=0.0, >DMin, >=0.0 and <=1.0, >=0.0] or DEFAULT [DEFAULT=[3.0, 0.0, 0.0, 12.0, 0.65]",  &
         (/3.0_ReKi, 0.0_ReKi, 0.0_ReKi,12.0_ReKi, 0.65_ReKi/), ErrStat2, ErrMsg2, UnEc); if(failed()) return
      WD_InitInp%WAT_k_Grad_k_c  = TmpRAry5(1)
      WD_InitInp%WAT_k_Grad_FMin = TmpRAry5(2)
      WD_InitInp%WAT_k_Grad_DMin = TmpRAry5(3)
      WD_InitInp%WAT_k_Grad_DMax = TmpRAry5(4)
      WD_InitInp%WAT_k_Grad_Exp  = TmpRAry5(5)
   if ( PathIsRelative( p%WAT_BoxFile ) ) p%WAT_BoxFile = TRIM(PriPath)//TRIM(p%WAT_BoxFile)
   if (p%WAT > 0_IntKi) WD_InitInp%WAT = .true.


   !---------------------- VISUALIZATION --------------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Visualization', ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   CALL ReadVar( UnIn, InputFile, AWAE_InitInp%WrDisWind, "WrDisWind", "Write disturbed wind data to <OutFileRoot>.Low.Dis.t<n/n_low-out>.vtk etc.? (flag)", ErrStat2, ErrMsg2, UnEc); if (Failed()) return

   ! XY planes
   CALL ReadVar( UnIn, InputFile, AWAE_InitInp%NOutDisWindXY, "NOutDisWindXY", "Number of XY planes for output of disturbed wind data across the low-resolution domain to <OutFileRoot>.Low.DisXY.<n_out>.t<n/n_low-out>.vtk (-) [0 to 999]", ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   call allocAry( AWAE_InitInp%OutDisWindZ, AWAE_InitInp%NOutDisWindXY, "OutDisWindZ", ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL ReadAry( UnIn, InputFile, AWAE_InitInp%OutDisWindZ, AWAE_InitInp%NOutDisWindXY, "OutDisWindZ", "Z coordinates of XY planes for output of disturbed wind data across the low-resolution domain (m) [1 to NOutDisWindXY] [unused for NOutDisWindXY=0]", ErrStat2, ErrMsg2, UnEc); if (Failed()) return

   ! YZ planes
   CALL ReadVar( UnIn, InputFile, AWAE_InitInp%NOutDisWindYZ, "NOutDisWindYZ", "Number of YZ planes for output of disturbed wind data across the low-resolution domain to <OutFileRoot>.Low.DisYZ.<n_out>.t<n/n_low-out>.vtk (-) [0 to 999]", ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   call allocAry( AWAE_InitInp%OutDisWindX, AWAE_InitInp%NOutDisWindYZ, "OutDisWindX", ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL ReadAry( UnIn, InputFile, AWAE_InitInp%OutDisWindX, AWAE_InitInp%NOutDisWindYZ, "OutDisWindX", "X coordinates of YZ planes for output of disturbed wind data across the low-resolution domain (m) [1 to NOutDisWindYZ] [unused for NOutDisWindYZ=0]", ErrStat2, ErrMsg2, UnEc); if (Failed()) return

   ! XZ planes
   CALL ReadVar( UnIn, InputFile, AWAE_InitInp%NOutDisWindXZ, "NOutDisWindXZ", "Number of XZ planes for output of disturbed wind data across the low-resolution domain to <OutFileRoot>.Low/DisXZ.<n_out>.t<n/n_low-out>.vtk (-) [0 to 999]", ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   call allocAry( AWAE_InitInp%OutDisWindY, AWAE_InitInp%NOutDisWindXZ, "OutDisWindY", ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL ReadAry( UnIn, InputFile, AWAE_InitInp%OutDisWindY, AWAE_InitInp%NOutDisWindXZ, "OutDisWindY", "Y coordinates of XZ planes for output of disturbed wind data across the low-resolution domain (m) [1 to NOutDisWindXZ] [unused for NOutDisWindXZ=0]", ErrStat2, ErrMsg2, UnEc); if (Failed()) return

   CALL ReadVarWDefault( UnIn, InputFile, AWAE_InitInp%WrDisDT, "WrDisDT", "The time between vtk outputs [must be a multiple of the low resolution time step]", p%DT_low, ErrStat2, ErrMsg2, UnEc); if (Failed()) return

   !---------------------- OUTPUT --------------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Output', ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   CALL ReadVar( UnIn, InputFile, p%SumPrint, "SumPrint", "Print summary data to <RootName>.sum? (flag)", ErrStat2, ErrMsg2, UnEc); if (Failed()) return

   CALL ReadVar( UnIn, InputFile, TmpTime, "ChkptTime", "Amount of time between creating checkpoint files for potential restart (s) [>0.0]", ErrStat2, ErrMsg2, UnEc); if (Failed()) return
      IF (TmpTime > p%TMax) THEN
         p%n_ChkptTime = HUGE(p%n_ChkptTime)
      ELSE
         p%n_ChkptTime = NINT( TmpTime / p%DT_low )
      END IF

   CALL ReadVar( UnIn, InputFile, p%TStart, "TStart", "Time to begin tabular output (s) [>=0.0]", ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   CALL ReadVar( UnIn, InputFile, OutFileFmt, "OutFileFmt", "Format for tabular (time-marching) output file (switch) {1: text file [<RootName>.out], 2: binary file [<RootName>.outb], 3: both}", ErrStat2, ErrMsg2, UnEc); if (Failed()) return
      SELECT CASE (OutFileFmt)
         CASE (1_IntKi)
            p%WrBinOutFile = .FALSE.
            p%WrTxtOutFile = .TRUE.
         CASE (2_IntKi)
            p%WrBinOutFile = .TRUE.
            p%WrTxtOutFile = .FALSE.
         CASE (3_IntKi)
            p%WrBinOutFile = .TRUE.
            p%WrTxtOutFile = .TRUE.
         CASE DEFAULT
            ! we'll check this later....
            !CALL SetErrStat( ErrID_Fatal, "FAST.Farm's OutFileFmt must be 1, 2, or 3.",ErrStat,ErrMsg,RoutineName)
            !if ( ErrStat >= AbortErrLev ) then
            !   call cleanup()
            !   RETURN
            !end if
      END SELECT

      if ( OutFileFmt /= 1_IntKi ) then ! TODO: Only allow text format for now; add binary format later.
         CALL SetErrStat( ErrID_Fatal, "FAST.Farm's OutFileFmt must be 1.",ErrStat,ErrMsg,RoutineName)
         call cleanup()
         RETURN
      end if

   CALL ReadVar( UnIn, InputFile, TabDelim, "TabDelim", "Use tab delimiters in text tabular output file? (flag) {uses spaces if False}", ErrStat2, ErrMsg2, UnEc); if (Failed()) return
      IF ( TabDelim ) THEN
         p%Delim = TAB
      ELSE
         p%Delim = ' '
      END IF

   CALL ReadVar( UnIn, InputFile, p%OutFmt, "OutFmt", "Format used for text tabular output, excluding the time channel. Resulting field should be 10 characters. (quoted string)", ErrStat2, ErrMsg2, UnEc); if (Failed()) return

   CALL ReadVarWDefault( UnIn, InputFile, WD_InitInp%OutAllPlanes,    "OutAllPlanes", "Output all planes", .False., ErrStat2, ErrMsg2, UnEc); if(failed()) return

   ! OutRadii
   CALL ReadVar( UnIn, InputFile, p%NOutRadii, "NOutRadii", "Number of radial nodes for wake output for an individual rotor (-) [0 to 20]", ErrStat2, ErrMsg2, UnEc); if (Failed()) return
      call allocary( p%OutRadii,  p%NOutRadii, "OutRadii", ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL ReadAry( UnIn, InputFile, p%OutRadii, p%NOutRadii, "OutRadii", "List of radial nodes for wake output for an individual rotor (-) [1 to NOutRadii]", ErrStat2, ErrMsg2, UnEc); if (Failed()) return

   ! OutDist
   CALL ReadVar( UnIn, InputFile, p%NOutDist, "NOutDist", "Number of downstream distances for wake output for an individual rotor (-) [0 to 9]", ErrStat2, ErrMsg2, UnEc); if (Failed()) return
      call allocary( p%OutDist,   p%NOutDist, "OutDist", ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL ReadAry( UnIn, InputFile, p%OutDist, p%NOutDist, "OutDist", "List of downstream distances for wake output for an individual rotor (m) [1 to NOutDist] [unused for NOutDist=0]", ErrStat2, ErrMsg2, UnEc); if (Failed()) return

   ! WindVel
   CALL ReadVar( UnIn, InputFile, p%NWindVel, "NWindVel", "Number of points for wind output (-) [0 to 9]", ErrStat2, ErrMsg2, UnEc); if (Failed()) return
      call allocAry( p%WindVelX,  p%NWindVel, "WindVelX", ErrStat2, ErrMsg2 );  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName); if (Failed()) return
      call allocAry( p%WindVelY,  p%NWindVel, "WindVelY", ErrStat2, ErrMsg2 );  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName); if (Failed()) return
      call allocAry( p%WindVelZ,  p%NWindVel, "WindVelZ", ErrStat2, ErrMsg2 );  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName); if (Failed()) return
   CALL ReadAry( UnIn, InputFile, p%WindVelX, p%NWindVel, "WindVelX", "List of coordinates in the X direction for wind output (m) [1 to NWindVel] [unused for NWindVel=0]", ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   CALL ReadAry( UnIn, InputFile, p%WindVelY, p%NWindVel, "WindVelY", "List of coordinates in the Y direction for wind output (m) [1 to NWindVel] [unused for NWindVel=0]", ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   CALL ReadAry( UnIn, InputFile, p%WindVelZ, p%NWindVel, "WindVelZ", "List of coordinates in the Z direction for wind output (m) [1 to NWindVel] [unused for NWindVel=0]", ErrStat2, ErrMsg2, UnEc); if (Failed()) return


      !---------------------- OUTLIST  --------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: OutList', ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   CALL ReadOutputList ( UnIn, InputFile, OutList, p%NumOuts, 'OutList', "List of user-requested output channels", ErrStat2, ErrMsg2, UnEc  ); if (Failed()) return     ! Routine in NWTC Subroutine Library


   !---------------------- END OF FILE -----------------------------------------

   call cleanup()
   RETURN

CONTAINS
   !...............................................................................................................................
   subroutine cleanup()
      CLOSE( UnIn )
      IF ( UnEc > 0 ) CLOSE ( UnEc )
   end subroutine cleanup

   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed =  ErrStat >= AbortErrLev
      if (Failed) call cleanup()
   end function Failed
   !...............................................................................................................................
END SUBROUTINE Farm_ReadPrimaryFile
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Farm_ValidateInput( p, WD_InitInp, AWAE_InitInp, ErrStat, ErrMsg )
      ! Passed variables
   TYPE(Farm_ParameterType), INTENT(INOUT) :: p                               !< The parameter data for the FAST (glue-code) simulation
   TYPE(WD_InputFileType),   INTENT(IN   ) :: WD_InitInp                      !< input-file data for WakeDynamics module
   TYPE(AWAE_InputFileType), INTENT(INOUT) :: AWAE_InitInp                    !< input-file data for AWAE module
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat                         !< Error status
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg                          !< Error message

      ! Local variables:
   INTEGER(IntKi)                :: i
   INTEGER(IntKi)                :: ErrStat2                                  ! Temporary Error status
   CHARACTER(ErrMsgLen)          :: ErrMsg2                                   ! Temporary Error message
   CHARACTER(*),   PARAMETER     :: RoutineName = 'Farm_ValidateInput'
   INTEGER(IntKi)                :: n_disDT_dt
   character(60)                 :: tmpStr

   ErrStat = ErrID_None
   ErrMsg  = ""


   ! --- SIMULATION CONTROL ---
   IF ((p%WaveFieldMod .ne. 1) .and. (p%WaveFieldMod .ne. 2)) CALL SetErrStat(ErrID_Fatal,'WaveFieldMod must be 1 or 2.',ErrStat,ErrMsg,RoutineName)
   IF ((p%MooringMod .ne. 0) .and. (p%MooringMod .ne. 3)) CALL SetErrStat(ErrID_Fatal,'MooringMod must be 0 or 3.',ErrStat,ErrMsg,RoutineName)


   IF (p%DT_low <= 0.0_ReKi) CALL SetErrStat(ErrID_Fatal,'DT_low must be positive.',ErrStat,ErrMsg,RoutineName)
   IF (p%DT_high <= 0.0_ReKi) CALL SetErrStat(ErrID_Fatal,'DT_high must be positive.',ErrStat,ErrMsg,RoutineName)
   IF (p%TMax < 0.0_ReKi) CALL SetErrStat(ErrID_Fatal,'TMax must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF (p%NumTurbines < 1) CALL SetErrStat(ErrID_Fatal,'FAST.Farm requires at least 1 turbine. Set NumTurbines > 0.',ErrStat,ErrMsg,RoutineName)

   ! --- SUPER CONTROLLER ---
   ! TODO : Verify that the DLL file exists

   ! --- SHARED MOORING SYSTEM ---
   ! TODO : Verify that p%MD_FileName file exists
   if ((p%DT_mooring <= 0.0_ReKi) .or. (p%DT_mooring > p%DT_high)) CALL SetErrStat(ErrID_Fatal,'DT_mooring must be greater than zero and no greater than dt_high.',ErrStat,ErrMsg,RoutineName)

   ! --- AMBIENT WIND: INFLOWWIND MODULE --- [used only for Mod_AmbWind=2 or 3] ---
   ! FIXME: this really should be checked with the turbine specific size diameter -- maybe relocate this check to AWAE or in FF after initializing all turbines?
   if (AWAE_InitInp%Mod_AmbWind > 1) then
      ! check that the grid is large enough to contain the turbine (only check Y and Z)
      do i=1,p%NumTurbines
         if (AWAE_InitInp%nY_High*AWAE_InitInp%dY_high(i) < p%RotorDiamRef) call SetErrStat(ErrID_Warn,'High res domain for turbine '//trim(Num2LStr(i))//' may be too small in Y (nY_High*dY_High < RotorDiamRef)',ErrStat,ErrMsg,RoutineName)
         if (AWAE_InitInp%nZ_High*AWAE_InitInp%dZ_high(i) < p%RotorDiamRef) call SetErrStat(ErrID_Warn,'High res domain for turbine '//trim(Num2LStr(i))//' may be too small in Z (nZ_High*dZ_High < RotorDiamRef)',ErrStat,ErrMsg,RoutineName)
      enddo
   endif

   ! --- WAKE DYNAMICS ---
   IF (WD_InitInp%Mod_Wake < 1 .or. WD_InitInp%Mod_Wake >3 ) CALL SetErrStat(ErrID_Fatal,'Mod_Wake needs to be 1,2 or 3',ErrStat,ErrMsg,RoutineName)
   IF (WD_InitInp%dr <= 0.0_ReKi) CALL SetErrStat(ErrID_Fatal,'dr (radial increment) must be larger than 0.',ErrStat,ErrMsg,RoutineName)
   IF (WD_InitInp%NumRadii < 2) CALL SetErrStat(ErrID_Fatal,'NumRadii (number of radii) must be at least 2.',ErrStat,ErrMsg,RoutineName)
   IF (WD_InitInp%NumPlanes < 2) CALL SetErrStat(ErrID_Fatal,'NumPlanes (number of wake planes) must be at least 2.',ErrStat,ErrMsg,RoutineName)

   IF (WD_InitInp%k_VortexDecay < 0.0_ReKi) CALL SetErrStat(ErrID_Fatal,'k_VortexDecay must be >= 0',ErrStat,ErrMsg,RoutineName)
   IF (WD_InitInp%NumVortices < 2) CALL SetErrStat(ErrID_Fatal,'NumVorticies must be greater than 1',ErrStat,ErrMsg,RoutineName)
   IF (WD_InitInp%sigma_D < 0.0_ReKi) CALL SetErrStat(ErrID_Fatal,'sigma_D must be postive',ErrStat,ErrMsg,RoutineName)
   IF (WD_InitInp%f_c <= 0.0_ReKi) CALL SetErrStat(ErrID_Fatal,'f_c (cut-off [corner] frequency) must be more than 0 Hz.',ErrStat,ErrMsg,RoutineName)
   IF (WD_InitInp%C_NearWake <= 1.0_Reki) CALL SetErrStat(ErrID_Fatal,'C_NearWake parameter must be greater than 1.',ErrStat,ErrMsg,RoutineName)
   IF (WD_InitInp%k_vCurl < 0.0_Reki) CALL SetErrStat(ErrID_Fatal,'k_vCurl parameter must not be negative.',ErrStat,ErrMsg,RoutineName)

   IF (WD_InitInp%k_vAmb      <  0.0_Reki)                                        CALL SetErrStat(ErrID_Fatal,'k_vAmb(1) (k_vAmb) must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF (WD_InitInp%C_vAmb_FMin <  0.0_Reki .or. WD_InitInp%C_vAmb_FMin > 1.0_Reki) CALL SetErrStat(ErrID_Fatal,'k_vAmb(2) (FMin) must be between 0 and 1 (inclusive).',ErrStat,ErrMsg,RoutineName)
   IF (WD_InitInp%C_vAmb_DMin <  0.0_Reki)                                        CALL SetErrStat(ErrID_Fatal,'k_vAmb(3) (DMin) must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF (WD_InitInp%C_vAmb_DMax <= WD_InitInp%C_vAmb_DMin)                          CALL SetErrStat(ErrID_Fatal,'k_vAmb(4) (DMax) must be larger than DMin.',ErrStat,ErrMsg,RoutineName)
   IF (WD_InitInp%C_vAmb_Exp  <  0.0_Reki)                                        CALL SetErrStat(ErrID_Fatal,'k_vAmb(5) (e) must be >=0.',ErrStat,ErrMsg,RoutineName)

   IF (WD_InitInp%k_vShr      <  0.0_Reki)                                        CALL SetErrStat(ErrID_Fatal,'k_vShr(1) (k_vShr) must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF (WD_InitInp%C_vShr_FMin <  0.0_Reki .or. WD_InitInp%C_vShr_FMin > 1.0_ReKi) CALL SetErrStat(ErrID_Fatal,'k_vShr(2) (FMin) must be between 0 and 1 (inclusive).',ErrStat,ErrMsg,RoutineName)
   IF (WD_InitInp%C_vShr_DMin <  0.0_Reki)                                        CALL SetErrStat(ErrID_Fatal,'k_vShr(3) (DMin) must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF (WD_InitInp%C_vShr_DMax <= WD_InitInp%C_vShr_DMin)                          CALL SetErrStat(ErrID_Fatal,'k_vShr(4) (DMax) must be larger than DMin.',ErrStat,ErrMsg,RoutineName)
   IF (WD_InitInp%C_vShr_Exp  <  0.0_Reki)                                        CALL SetErrStat(ErrID_Fatal,'k_vShr(5) (e) must be >=0.',ErrStat,ErrMsg,RoutineName)

   IF (WD_InitInp%Mod_WakeDiam < WakeDiamMod_RotDiam .or. WD_InitInp%Mod_WakeDiam > WakeDiamMod_MtmFlux) THEN
      call SetErrStat(ErrID_Fatal,'Wake diameter calculation model, Mod_WakeDiam, must be 1 (rotor diameter), 2 (velocity-based), 3 (mass-flux based), 4 (momentum-flux based) or DEFAULT.',ErrStat,ErrMsg,RoutineName)
   END IF

   IF (WD_InitInp%Mod_WakeDiam /= WakeDiamMod_RotDiam) THEN
      IF (WD_InitInp%C_WakeDiam <= 0.0_Reki .or. WD_InitInp%C_WakeDiam >= 1.0_ReKi) THEN
         CALL SetErrStat(ErrID_Fatal,'C_WakeDiam parameter must be between 0 and 1 (exclusive).',ErrStat,ErrMsg,RoutineName)
      END IF
   END IF



   IF (AWAE_InitInp%C_Meander < 1.0_Reki) THEN
      CALL SetErrStat(ErrID_Fatal,'C_Meander parameter must not be less than 1.',ErrStat,ErrMsg,RoutineName)
   END IF

   ! --- CURL
   IF (WD_InitInp%FilterInit < 0  ) CALL SetErrStat(ErrID_Fatal,'FilterInit needs to >= 0',ErrStat,ErrMsg,RoutineName)
   IF (AWAE_InitInp%Mod_Meander < MeanderMod_Uniform .or. AWAE_InitInp%Mod_Meander > MeanderMod_WndwdJinc) THEN
      call SetErrStat(ErrID_Fatal,'Spatial filter model for wake meandering, Mod_Meander, must be 1 (uniform), 2 (truncated jinc), 3 (windowed jinc) or DEFAULT.',ErrStat,ErrMsg,RoutineName)
   END IF
   IF (.not.(ANY((/1,2,3/)==AWAE_InitInp%Mod_Projection))) CALL SetErrStat(ErrID_Fatal,'Mod_Projection needs to be 1, 2 or 3',ErrStat,ErrMsg,RoutineName)

   ! --- WAT
   if (p%WAT < 0_IntKi .or. p%WAT > 2_IntKi) CALL SetErrStat(ErrID_Fatal,'WAT option must be 0: no wake added turbulence, 1: predefined turbulence box, or 2: user defined turbulence box.',ErrStat,ErrMsg,RoutineName)
   if (p%WAT>0) then
      ! Checks on k_Def
      if (WD_InitInp%WAT_k_Def_k_c  <= 0.0_ReKi)                                           call SetErrStat(ErrID_Fatal,'WAT_k_Def(1) (k_def) must be >0.',ErrStat,ErrMsg,RoutineName)
      if (WD_InitInp%WAT_k_Def_FMin <  0.0_ReKi .or. WD_InitInp%WAT_k_Def_FMin > 1.0_ReKi) call SetErrStat(ErrID_Fatal,'WAT_k_Def(2) (f_min) must be >=0 and <=1.',ErrStat,ErrMsg,RoutineName)
      if (WD_InitInp%WAT_k_Def_DMin <  0.0_ReKi)                                           call SetErrStat(ErrID_Fatal,'WAT_k_Def(3) (D_min) must be >=0.',ErrStat,ErrMsg,RoutineName)
      if (WD_InitInp%WAT_k_Def_DMax <= WD_InitInp%WAT_k_Def_DMin)                          call SetErrStat(ErrID_Fatal,'WAT_k_Def(4) (D_max) must be greater than D_min.',ErrStat,ErrMsg,RoutineName)
      if (WD_InitInp%WAT_k_Def_Exp  <  0.0_ReKi)                                           call SetErrStat(ErrID_Fatal,'WAT_k_Def(5) (e) must be >=0.',ErrStat,ErrMsg,RoutineName)
      ! Tests on k_Grad
      if (WD_InitInp%WAT_k_Grad_k_c  <= 0.0_ReKi)                                            call SetErrStat(ErrID_Fatal,'WAT_k_Grad(1) (k_def) must be >0.',ErrStat,ErrMsg,RoutineName)
      if (WD_InitInp%WAT_k_Grad_FMin <  0.0_ReKi .or. WD_InitInp%WAT_k_Grad_FMin > 1.0_ReKi) call SetErrStat(ErrID_Fatal,'WAT_k_Grad(2) (f_min) must be >=0 and <=1.',ErrStat,ErrMsg,RoutineName)
      if (WD_InitInp%WAT_k_Grad_DMin <  0.0_ReKi)                                            call SetErrStat(ErrID_Fatal,'WAT_k_Grad(3) (D_min) must be >=0.',ErrStat,ErrMsg,RoutineName)
      if (WD_InitInp%WAT_k_Grad_DMax <= WD_InitInp%WAT_k_Grad_DMin)                          call SetErrStat(ErrID_Fatal,'WAT_k_Grad(4) (D_max) must be greater than D_min.',ErrStat,ErrMsg,RoutineName)
      if (WD_InitInp%WAT_k_Grad_Exp  <  0.0_ReKi)                                            call SetErrStat(ErrID_Fatal,'WAT_k_Grad(5) (e) must be >=0.',ErrStat,ErrMsg,RoutineName)
      ! summary table
      call WrScr('  Wake-Added Turbulence (WAT): coefficients:')
      call WrScr('                k_c      f_min    D_min    D_max    e')
      write(tmpStr,'(A6,A6,6(f9.3))') '','k_Def', WD_InitInp%WAT_k_Def_k_c, WD_InitInp%WAT_k_Def_FMin, WD_InitInp%WAT_k_Def_DMin, WD_InitInp%WAT_k_Def_DMax, WD_InitInp%WAT_k_Def_Exp
      call WrScr(tmpStr)
      write(tmpStr,'(A6,A6,6(f9.3))') '','k_Grad',WD_InitInp%WAT_k_Grad_k_c,WD_InitInp%WAT_k_Grad_FMin,WD_InitInp%WAT_k_Grad_DMin,WD_InitInp%WAT_k_Grad_DMax,WD_InitInp%WAT_k_Grad_Exp
      call WrScr(tmpStr)
   endif

   !--- OUTPUT ---
   IF ( p%n_ChkptTime < 1_IntKi   ) CALL SetErrStat( ErrID_Fatal, 'ChkptTime must be greater than 0 seconds.', ErrStat, ErrMsg, RoutineName )
   IF (p%TStart < 0.0_ReKi) CALL SetErrStat(ErrID_Fatal,'TStart must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF (.not. p%WrBinOutFile .and. .not. p%WrTxtOutFile) CALL SetErrStat( ErrID_Fatal, "FAST.Farm's OutFileFmt must be 1, 2, or 3.",ErrStat,ErrMsg,RoutineName)

   if (AWAE_InitInp%WrDisDT < p%DT_low) CALL SetErrStat(ErrID_Fatal,'WrDisDT must greater than or equal to dt_low.',ErrStat,ErrMsg,RoutineName)

      ! let's make sure the FAST.Farm DT_low is an exact integer divisor of AWAE_InitInp%WrDisDT
   n_disDT_dt = nint( AWAE_InitInp%WrDisDT / p%DT_low )
      ! (i'm doing this outside of Farm_ValidateInput so we know that dt_low/=0 before computing n_high_low):
   IF ( .NOT. EqualRealNos( real(p%DT_low,SiKi)* n_disDT_dt, real(AWAE_InitInp%WrDisDT,SiKi)  )  ) THEN
      CALL SetErrStat(ErrID_Fatal, "WrDisDT ("//TRIM(Num2LStr(AWAE_InitInp%WrDisDT))//" s) must be an integer multiple of dt_low ("//TRIM(Num2LStr(p%DT_low))//" s).", ErrStat, ErrMsg, RoutineName )
   END IF
   AWAE_InitInp%WrDisDT =  p%DT_low * n_disDT_dt


   if (AWAE_InitInp%NOutDisWindXY < 0 .or. AWAE_InitInp%NOutDisWindXY > maxOutputPlanes ) CALL SetErrStat( ErrID_Fatal, 'NOutDisWindXY must be in the range [0, 999].', ErrStat, ErrMsg, RoutineName )
   if (AWAE_InitInp%NOutDisWindYZ < 0 .or. AWAE_InitInp%NOutDisWindYZ > maxOutputPlanes ) CALL SetErrStat( ErrID_Fatal, 'NOutDisWindYZ must be in the range [0, 999].', ErrStat, ErrMsg, RoutineName )
   if (AWAE_InitInp%NOutDisWindXZ < 0 .or. AWAE_InitInp%NOutDisWindXZ > maxOutputPlanes ) CALL SetErrStat( ErrID_Fatal, 'NOutDisWindXZ must be in the range [0, 999].', ErrStat, ErrMsg, RoutineName )
   if (p%NOutDist < 0 .or. p%NOutDist > maxOutputPoints ) then
      CALL SetErrStat( ErrID_Fatal, 'NOutDist must be in the range [0, 9].', ErrStat, ErrMsg, RoutineName )
   else
      do i=1,p%NOutDist
         if (p%OutDist(i) <  0.0_ReKi) then
            CALL SetErrStat( ErrID_Fatal, 'OutDist values must be greater than or equal to zero.', ErrStat, ErrMsg, RoutineName )
            exit
         end if
      end do
   end if

   if (p%NWindVel < 0 .or. p%NWindVel > maxOutputPoints ) CALL SetErrStat( ErrID_Fatal, 'NWindVel must be in the range [0, 9].', ErrStat, ErrMsg, RoutineName )
   if (p%NOutRadii < 0 .or. p%NOutRadii > 20 ) then
      CALL SetErrStat( ErrID_Fatal, 'NOutRadii must be in the range [0, 20].', ErrStat, ErrMsg, RoutineName )
   else
      do i=1,p%NOutRadii
         if (p%OutRadii(i) > WD_InitInp%NumRadii - 1 .or. p%OutRadii(i) < 0) then
            CALL SetErrStat( ErrID_Fatal, 'OutRadii must be in the range [0, NumRadii - 1].', ErrStat, ErrMsg, RoutineName )
            exit
         end if
      end do
   end if


      ! Check that OutFmt is a valid format specifier and will fit over the column headings
  CALL ChkRealFmtStr( p%OutFmt, 'OutFmt', p%FmtWidth, ErrStat2, ErrMsg2 ) !this sets p%FmtWidth!
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   IF ( p%FmtWidth > ChanLen ) CALL SetErrStat( ErrID_Warn, 'OutFmt produces a column width of '// &
         TRIM(Num2LStr(p%FmtWidth))//' instead of '//TRIM(Num2LStr(ChanLen))//' characters.', ErrStat, ErrMsg, RoutineName )


END SUBROUTINE Farm_ValidateInput


end module FAST_Farm_IO
