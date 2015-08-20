!**********************************************************************************************************************************
! The FAST_Prog.f90, FAST_IO.f90, and FAST_Mods.f90 make up the FAST glue code in the FAST Modularization Framework.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2015  National Renewable Energy Laboratory
!
!    This file is part of FAST.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
!**********************************************************************************************************************************
! File last committed: $Date$
! (File) Revision #: $Rev$
! URL: $HeadURL$
!**********************************************************************************************************************************
MODULE FAST_Subs

   USE NWTC_Library
   USE NWTC_LAPACK

   USE FAST_ModTypes
      
   USE AeroDyn
   USE AeroDyn14
   USE InflowWind
   USE ElastoDyn
   USE BeamDyn
   USE FEAMooring
   USE MoorDyn
   USE HydroDyn
   USE IceDyn
   USE IceFloe
   USE MAP
   USE ServoDyn
   USE SubDyn
   USE OpenFOAM
   

   IMPLICIT NONE

CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
FUNCTION GetVersion()
! This function returns a string describing the glue code and some of the compilation options we're using.
!..................................................................................................................................

   IMPLICIT                        NONE


   ! Passed Variables:

   CHARACTER(1024)  :: GetVersion                      ! String containing a description of the compiled precision.



   GetVersion = TRIM(GetNVD(FAST_Ver))//', compiled'

   IF ( Cmpl4SFun )  THEN     ! FAST has been compiled as an S-Function for Simulink
      GetVersion = TRIM(GetVersion)//' as a DLL S-Function for Simulink'
   ELSEIF ( Cmpl4LV )  THEN     ! FAST has been compiled as a DLL for Labview
      GetVersion = TRIM(GetVersion)//' as a DLL for LabVIEW'
   ENDIF   
   
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
SUBROUTINE FAST_EndOutput( p_FAST, y_FAST, ErrStat, ErrMsg )
! This subroutine is called at program termination. It writes any additional output files,
! deallocates variables and closes files.
!----------------------------------------------------------------------------------------------------

   TYPE(FAST_ParameterType), INTENT(INOUT) :: p_FAST                    ! FAST Parameters
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST                    ! FAST Output

   INTEGER(IntKi),           INTENT(OUT)   :: ErrStat                   ! Error status
   CHARACTER(*),             INTENT(OUT)   :: ErrMsg                    ! Message associated with errro status

      ! local variables
   CHARACTER(LEN(y_FAST%FileDescLines)*3)  :: FileDesc                  ! The description of the run, to be written in the binary output file


      ! Initialize some values

   ErrStat = ErrID_None
   ErrMsg  = ''

   !-------------------------------------------------------------------------------------------------
   ! Write the binary output file if requested
   !-------------------------------------------------------------------------------------------------

   IF (p_FAST%WrBinOutFile .AND. y_FAST%n_Out > 0) THEN

      FileDesc = TRIM(y_FAST%FileDescLines(1))//' '//TRIM(y_FAST%FileDescLines(2))//'; '//TRIM(y_FAST%FileDescLines(3))

      CALL WrBinFAST(TRIM(p_FAST%OutFileRoot)//'.outb', OutputFileFmtID, TRIM(FileDesc), &
            y_FAST%ChannelNames, y_FAST%ChannelUnits, y_FAST%TimeData, y_FAST%AllOutData(:,1:y_FAST%n_Out), ErrStat, ErrMsg)

      IF ( ErrStat /= ErrID_None ) CALL WrScr( TRIM(GetErrStr(ErrStat))//' when writing binary output file: '//TRIM(ErrMsg) )

   END IF


   !-------------------------------------------------------------------------------------------------
   ! Close the text tabular output file and summary file (if opened)
   !-------------------------------------------------------------------------------------------------
   IF (y_FAST%UnOu  > 0) THEN ! I/O unit number for the tabular output file
      CLOSE( y_FAST%UnOu )         
      y_FAST%UnOu = -1
   END IF
   
   IF (y_FAST%UnSum > 0) THEN ! I/O unit number for the tabular output file
      CLOSE( y_FAST%UnSum )        
      y_FAST%UnSum = -1
   END IF

   IF (y_FAST%UnGra > 0) THEN ! I/O unit number for the graphics output file
      CLOSE( y_FAST%UnGra )        
      y_FAST%UnGra = -1
   END IF
   
   !-------------------------------------------------------------------------------------------------
   ! Deallocate arrays
   !-------------------------------------------------------------------------------------------------

      ! Output
   IF ( ALLOCATED(y_FAST%AllOutData                  ) ) DEALLOCATE(y_FAST%AllOutData                  )
   IF ( ALLOCATED(y_FAST%TimeData                    ) ) DEALLOCATE(y_FAST%TimeData                    )
   IF ( ALLOCATED(y_FAST%ChannelNames                ) ) DEALLOCATE(y_FAST%ChannelNames                )
   IF ( ALLOCATED(y_FAST%ChannelUnits                ) ) DEALLOCATE(y_FAST%ChannelUnits                )


END SUBROUTINE FAST_EndOutput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FAST_Init( p, y_FAST, t_initial, ErrStat, ErrMsg, InFile, TMax  )
! This subroutine checks for command-line arguments, gets the root name of the input files
! (including full path name), and creates the names of the output files.
!..................................................................................................................................

      IMPLICIT                        NONE

   ! Passed variables

   TYPE(FAST_ParameterType), INTENT(INOUT)         :: p                 ! The parameter data for the FAST (glue-code) simulation
   TYPE(FAST_OutputFileType),INTENT(INOUT)         :: y_FAST            ! The output data for the FAST (glue-code) simulation
   REAL(DbKi),               INTENT(IN)            :: t_initial         ! the beginning time of the simulation
   INTEGER(IntKi),           INTENT(OUT)           :: ErrStat           ! Error status
   CHARACTER(*),             INTENT(OUT)           :: ErrMsg            ! Error message
   CHARACTER(*),             INTENT(IN), OPTIONAL  :: InFile            ! A CHARACTER string containing the name of the primary FAST input file (if not present, we'll get it from the command line)
   REAL(DbKi),               INTENT(IN), OPTIONAL  :: TMax              ! the length of the simulation (from Simulink)
      ! Local variables

   REAL(DbKi)                   :: TmpTime                              ! A temporary variable for error checking
   INTEGER                      :: i                                    ! loop counter
   !CHARACTER(1024)              :: DirName                              ! A CHARACTER string containing the path of the current working directory
   CHARACTER(1024)              :: LastArg                              ! A second command-line argument that will allow DWM module to be used in AeroDyn
   CHARACTER(1024)              :: InputFile                            ! A CHARACTER string containing the name of the primary FAST input file


   CHARACTER(*), PARAMETER      :: RoutineName = "FAST_Init"
   
   INTEGER(IntKi)               :: ErrStat2
   CHARACTER(1024)              :: ErrMsg2
   
      ! Initialize some variables
   ErrStat = ErrID_None
   ErrMsg = ''

      ! Display the copyright notice
   CALL DispCopyrightLicense( FAST_Ver )


      ! Tell our users what they're running
   CALL WrScr( ' Running '//GetVersion()//NewLine//' linked with '//TRIM( GetNVD( NWTC_Ver ))//NewLine )

   !...............................................................................................................................
   ! Get the name of the input file from the command line if it isn't an input to this routine
   ! and set the root name of the output files based on the input file name
   !...............................................................................................................................

   p%UseDWM = .FALSE.  ! by default, we're not going to use the DWM module
   
   IF ( PRESENT(InFile) ) THEN
      InputFile = InFile
      p%UseDWM  = .FALSE.
   ELSE ! get it from the command line
      InputFile = ""  ! initialize to empty string to make sure it's input from the command line
      CALL CheckArgs( InputFile, ErrStat2, LastArg )  ! if ErrStat2 /= ErrID_None, we'll ignore and deal with the problem when we try to read the input file
      
      IF (LEN_TRIM(InputFile) == 0) THEN ! no input file was specified
         CALL SetErrStat( ErrID_Fatal, 'The required input file was not specified on the command line.', ErrStat, ErrMsg, RoutineName )

            !bjj:  if people have compiled themselves, they should be able to figure out the file name, right?         
         IF (BITS_IN_ADDR==32) THEN
            CALL NWTC_DisplaySyntax( InputFile, 'FAST_Win32.exe' )
         ELSEIF( BITS_IN_ADDR == 64) THEN
            CALL NWTC_DisplaySyntax( InputFile, 'FAST_x64.exe' )
         ELSE
            CALL NWTC_DisplaySyntax( InputFile, 'FAST.exe' )
         END IF
         
         RETURN
      END IF            
      
      IF (LEN_TRIM(LastArg) > 0) THEN ! see if DWM was specified as the second option
         CALL Conv2UC( LastArg )
         IF ( TRIM(LastArg) == "DWM" ) THEN
            p%UseDWM    = .TRUE.
         END IF
      END IF
            
   END IF
      

      ! Determine the root name of the primary file (will be used for output files)
   CALL GetRoot( InputFile, p%OutFileRoot )
   IF ( Cmpl4SFun )  p%OutFileRoot = TRIM( p%OutFileRoot )//'.SFunc'
   
   !...............................................................................................................................
   ! Initialize the module name/date/version info:
   !...............................................................................................................................

   DO i=1,NumModules
      y_FAST%Module_Ver(i)%Date = 'unknown date'
      y_FAST%Module_Ver(i)%Ver  = 'unknown version'
   END DO       
   y_FAST%Module_Ver( Module_IfW  )%Name = 'InflowWind'
   y_FAST%Module_Ver( Module_OpFM )%Name = 'OpenFOAM integration'
   y_FAST%Module_Ver( Module_ED   )%Name = 'ElastoDyn'
   y_FAST%Module_Ver( Module_BD   )%Name = 'BeamDyn'
   y_FAST%Module_Ver( Module_AD14 )%Name = 'AeroDyn14'
   y_FAST%Module_Ver( Module_AD   )%Name = 'AeroDyn'
   y_FAST%Module_Ver( Module_SrvD )%Name = 'ServoDyn'
   y_FAST%Module_Ver( Module_HD   )%Name = 'HydroDyn'
   y_FAST%Module_Ver( Module_SD   )%Name = 'SubDyn'
   y_FAST%Module_Ver( Module_MAP  )%Name = 'MAP'
   y_FAST%Module_Ver( Module_FEAM )%Name = 'FEAMooring'
   y_FAST%Module_Ver( Module_MD   )%Name = 'MoorDyn'
   y_FAST%Module_Ver( Module_IceF )%Name = 'IceFloe'
   y_FAST%Module_Ver( Module_IceD )%Name = 'IceDyn'
         
   p%n_substeps = 1                                                ! number of substeps for between modules and global/FAST time
         
   !...............................................................................................................................
   ! Read the primary file for the glue code:
   !...............................................................................................................................
   CALL FAST_ReadPrimaryFile( InputFile, p, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
      
      ! overwrite TMax if necessary)
   IF (PRESENT(TMax)) THEN
      p%TMax = MAX( TMax, p%TMax )
   END IF
   
   IF ( ErrStat >= AbortErrLev ) RETURN


   p%KMax = 1                 ! after more checking, we may put this in the input file...
   !IF (p%CompIce == Module_IceF) p%KMax = 2
   p%SizeJac_ED_SD_HD_BD = 0  ! initialize this vector to zero; after we figure out what size the ED/SD/HD/BD meshes are, we'll fill this
   
   p%numIceLegs = 0           ! initialize number of support-structure legs in contact with ice (IceDyn will set this later)
   
   p%nBeams = 0               ! initialize number of BeamDyn instances (will be set later)
   
      ! determine what kind of turbine we're modeling:
   IF ( p%CompHydro == Module_HD ) THEN
      IF ( p%CompSub == Module_SD ) THEN
         p%TurbineType = Type_Offshore_Fixed
      ELSE
         p%TurbineType = Type_Offshore_Floating
      END IF
   ELSE
      p%TurbineType = Type_LandBased
   END IF   
         
    
    p%WrGraphics = .FALSE. !.TRUE.
#ifdef DEBUG_BEAMDYN
    p%WrGraphics = .TRUE.
#endif

    p%n_TMax_m1  = CEILING( ( (p%TMax - t_initial) / p%DT ) ) - 1 ! We're going to go from step 0 to n_TMax (thus the -1 here)

   !...............................................................................................................................
   ! Do some error checking on the inputs (validation):
   !...............................................................................................................................    
   IF ( p%TMax < 0.0_DbKi  )  THEN
      CALL SetErrStat( ErrID_Fatal, 'TMax must not be a negative number.', ErrStat, ErrMsg, RoutineName )
   ELSE IF ( p%TMax < p%TStart )  THEN
      CALL SetErrStat( ErrID_Fatal, 'TMax must not be less than TStart.', ErrStat, ErrMsg, RoutineName )
   END IF

   IF ( p%n_ChkptTime < p%n_TMax_m1 .AND. .NOT. p%WrBinOutFile) THEN
      CALL SetErrStat( ErrID_Severe, 'It is highly recommended that time-marching output files be generated in binary format when generating checkpoint files.', ErrStat, ErrMsg, RoutineName )
   END IF
      
   IF ( p%DT <= 0.0_DbKi )  THEN
      CALL SetErrStat( ErrID_Fatal, 'DT must be greater than 0.', ErrStat, ErrMsg, RoutineName )
   ELSE ! Test DT and TMax to ensure numerical stability -- HINT: see the use of OnePlusEps
      TmpTime = p%TMax*EPSILON(p%DT)
      IF ( p%DT <= TmpTime ) THEN
         CALL SetErrStat( ErrID_Fatal, 'DT must be greater than '//TRIM ( Num2LStr( TmpTime ) )//' seconds.', ErrStat, ErrMsg, RoutineName )
      END IF
   END IF

   IF ( p%WrTxtOutFile .AND. ( p%TMax > 9999.999_DbKi ) )  THEN
      CALL SetErrStat( ErrID_Fatal, 'TMax must not exceed 9999.999 seconds with text tabular (time-marching) output files.', ErrStat, ErrMsg, RoutineName )
   END IF

   IF ( p%TStart      <  0.0_DbKi ) CALL SetErrStat( ErrID_Fatal, 'TStart must not be less than 0 seconds.', ErrStat, ErrMsg, RoutineName )
!  IF ( p%SttsTime    <= 0.0_DbKi ) CALL SetErrStat( ErrID_Fatal, 'SttsTime must be greater than 0 seconds.', ErrStat, ErrMsg, RoutineName )
   IF ( p%n_SttsTime  < 1_IntKi   ) CALL SetErrStat( ErrID_Fatal, 'SttsTime must be greater than 0 seconds.', ErrStat, ErrMsg, RoutineName )
   IF ( p%n_ChkptTime < 1_IntKi   ) CALL SetErrStat( ErrID_Fatal, 'ChkptTime must be greater than 0 seconds.', ErrStat, ErrMsg, RoutineName )
   IF ( p%KMax        < 1_IntKi   ) CALL SetErrStat( ErrID_Fatal, 'KMax must be greater than 0.', ErrStat, ErrMsg, RoutineName )
   
   IF (p%CompElast   == Module_Unknown) CALL SetErrStat( ErrID_Fatal, 'CompElast must be 1 (ElastoDyn) or 2 (BeamDyn).', ErrStat, ErrMsg, RoutineName )   
   IF (p%CompAero    == Module_Unknown) CALL SetErrStat( ErrID_Fatal, 'CompAero must be 0 (None), 1 (AeroDyn14), or 2 (AeroDyn).', ErrStat, ErrMsg, RoutineName )
   IF (p%CompServo   == Module_Unknown) CALL SetErrStat( ErrID_Fatal, 'CompServo must be 0 (None) or 1 (ServoDyn).', ErrStat, ErrMsg, RoutineName )
   IF (p%CompHydro   == Module_Unknown) CALL SetErrStat( ErrID_Fatal, 'CompHydro must be 0 (None) or 1 (HydroDyn).', ErrStat, ErrMsg, RoutineName )
   IF (p%CompSub     == Module_Unknown) CALL SetErrStat( ErrID_Fatal, 'CompSub must be 0 (None) or 1 (SubDyn).', ErrStat, ErrMsg, RoutineName )
   IF (p%CompMooring == Module_Unknown) CALL SetErrStat( ErrID_Fatal, 'CompMooring must be 0 (None), 1 (MAP), 2 (FEAMooring), or 3 (MoorDyn).', ErrStat, ErrMsg, RoutineName )
   IF (p%CompIce     == Module_Unknown) CALL SetErrStat( ErrID_Fatal, 'CompIce must be 0 (None) or 1 (IceFloe).', ErrStat, ErrMsg, RoutineName )
   IF (p%CompHydro /= Module_HD) THEN
      IF (p%CompMooring == Module_MAP) THEN
         CALL SetErrStat( ErrID_Fatal, 'HydroDyn must be used when MAP is used. Set CompHydro > 0 or CompMooring = 0 in the FAST input file.', ErrStat, ErrMsg, RoutineName )
      ELSEIF (p%CompMooring == Module_FEAM) THEN
         CALL SetErrStat( ErrID_Fatal, 'HydroDyn must be used when FEAMooring is used. Set CompHydro > 0 or CompMooring = 0 in the FAST input file.', ErrStat, ErrMsg, RoutineName )
      END IF
   END IF
   
   IF (p%CompIce == Module_IceF) THEN
      IF (p%CompSub   /= Module_SD) CALL SetErrStat( ErrID_Fatal, 'SubDyn must be used when IceFloe is used. Set CompSub > 0 or CompIce = 0 in the FAST input file.', ErrStat, ErrMsg, RoutineName )
      IF (p%CompHydro /= Module_HD) CALL SetErrStat( ErrID_Fatal, 'HydroDyn must be used when IceFloe is used. Set CompHydro > 0 or CompIce = 0 in the FAST input file.', ErrStat, ErrMsg, RoutineName )
   ELSEIF (p%CompIce == Module_IceD) THEN
      IF (p%CompSub   /= Module_SD) CALL SetErrStat( ErrID_Fatal, 'SubDyn must be used when IceDyn is used. Set CompSub > 0 or CompIce = 0 in the FAST input file.', ErrStat, ErrMsg, RoutineName )
      IF (p%CompHydro /= Module_HD) CALL SetErrStat( ErrID_Fatal, 'HydroDyn must be used when IceDyn is used. Set CompHydro > 0 or CompIce = 0 in the FAST input file.', ErrStat, ErrMsg, RoutineName )
   END IF
   
   IF (p%CompElast == Module_BD .and. p%CompAero == Module_AD14 ) CALL SetErrStat( ErrID_Fatal, 'AeroDyn14 cannot be used when BeamDyn is used. Change CompAero or CompElast in the FAST input file.', ErrStat, ErrMsg, RoutineName )
   
!   IF ( p%InterpOrder < 0 .OR. p%InterpOrder > 2 ) THEN
   IF ( p%InterpOrder < 1 .OR. p%InterpOrder > 2 ) THEN
      CALL SetErrStat( ErrID_Fatal, 'InterpOrder must be 1 or 2.', ErrStat, ErrMsg, RoutineName ) ! 5/13/14 bjj: MAS and JMJ compromise for certain integrators is that InterpOrder cannot be 0
      p%InterpOrder = 1    ! Avoid problems in error handling by setting this to 0
   END IF

   IF ( p%NumCrctn < 0_IntKi ) THEN
      CALL SetErrStat( ErrID_Fatal, 'NumCrctn must be 0 or greater.', ErrStat, ErrMsg, RoutineName )
   END IF   
   
   
   IF ( ErrStat >= AbortErrLev ) RETURN

   
   !...............................................................................................................................

      ! temporary check on p_FAST%DT_out (bjj: fix this later)

   IF ( .NOT. EqualRealNos( p%DT_out, p%DT ) ) THEN
      IF ( p%DT_out < p%DT ) THEN
         CALL SetErrStat( ErrID_Fatal, 'DT_out must be at least DT ('//TRIM(Num2LStr(p%DT))//' s).', ErrStat, ErrMsg, RoutineName )
      ELSEIF ( .NOT. EqualRealNos( p%DT_out, p%DT * NINT(p%DT_out / p%DT ) )  ) THEN
         CALL SetErrStat( ErrID_Fatal, 'DT_out must currently be an integer multiple of DT.', ErrStat, ErrMsg, RoutineName )
      END IF
   END IF

   IF ( p%CompUserTwrLd ) THEN
      CALL SetErrStat( ErrID_Info, 'CompUserTwrLd will be ignored in this version of FAST.', ErrStat, ErrMsg, RoutineName )
      p%CompUserTwrLd = .FALSE.
   END IF

   IF ( p%CompUserPtfmLd ) THEN
      CALL SetErrStat( ErrID_Info, 'CompUserPtfmLd will be ignored in this version of FAST.', ErrStat, ErrMsg, RoutineName )
      p%CompUserPtfmLd = .FALSE.
   END IF
   
   RETURN
END SUBROUTINE FAST_Init
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FAST_InitOutput( p_FAST, y_FAST, InitOutData_ED, InitOutData_BD, InitOutData_SrvD, InitOutData_AD14, InitOutData_AD, &
                            InitOutData_IfW, InitOutData_OpFM, InitOutData_HD, InitOutData_SD, &
                            InitOutData_MAP, InitOutData_FEAM, InitOutData_MD, InitOutData_IceF, InitOutData_IceD, ErrStat, ErrMsg )
! This routine initializes the output for the glue code, including writing the header for the primary output file.
! was previously called WrOutHdr()
!..................................................................................................................................

   IMPLICIT NONE

      ! Passed variables
   TYPE(FAST_ParameterType),       INTENT(IN)           :: p_FAST                                ! Glue-code simulation parameters
   TYPE(FAST_OutputFileType),      INTENT(INOUT)        :: y_FAST                                ! Glue-code simulation outputs

   TYPE(ED_InitOutputType),        INTENT(IN)           :: InitOutData_ED                        ! Initialization output for ElastoDyn
   TYPE(BD_InitOutputType),        INTENT(IN)           :: InitOutData_BD(:)                     ! Initialization output for BeamDyn (each instance)
   TYPE(SrvD_InitOutputType),      INTENT(IN)           :: InitOutData_SrvD                      ! Initialization output for ServoDyn
   TYPE(AD14_InitOutputType),      INTENT(IN)           :: InitOutData_AD14                      ! Initialization output for AeroDyn14
   TYPE(AD_InitOutputType),        INTENT(IN)           :: InitOutData_AD                        ! Initialization output for AeroDyn
   TYPE(InflowWind_InitOutputType),INTENT(IN)           :: InitOutData_IfW                       ! Initialization output for InflowWind
   TYPE(OpFM_InitOutputType),      INTENT(IN)           :: InitOutData_OpFM                      ! Initialization output for OpenFOAM
   TYPE(HydroDyn_InitOutputType),  INTENT(IN)           :: InitOutData_HD                        ! Initialization output for HydroDyn
   TYPE(SD_InitOutputType),        INTENT(IN)           :: InitOutData_SD                        ! Initialization output for SubDyn
   TYPE(MAP_InitOutputType),       INTENT(IN)           :: InitOutData_MAP                       ! Initialization output for MAP
   TYPE(FEAM_InitOutputType),      INTENT(IN)           :: InitOutData_FEAM                      ! Initialization output for FEAMooring
   TYPE(MD_InitOutputType),        INTENT(IN)           :: InitOutData_MD                        ! Initialization output for MoorDyn
   TYPE(IceFloe_InitOutputType),   INTENT(IN)           :: InitOutData_IceF                      ! Initialization output for IceFloe
   TYPE(IceD_InitOutputType),      INTENT(IN)           :: InitOutData_IceD                      ! Initialization output for IceDyn

   INTEGER(IntKi),                 INTENT(OUT)          :: ErrStat                               ! Error status
   CHARACTER(*),                   INTENT(OUT)          :: ErrMsg                                ! Error message corresponding to ErrStat


      ! Local variables.

   INTEGER(IntKi)                   :: I, J                                            ! Generic index for DO loops.
   INTEGER(IntKi)                   :: indxLast                                        ! The index of the last value to be written to an array
   INTEGER(IntKi)                   :: indxNext                                        ! The index of the next value to be written to an array
   INTEGER(IntKi)                   :: NumOuts                                         ! number of channels to be written to the output file(s)



   !......................................................
   ! Set the description lines to be printed in the output file
   !......................................................
   y_FAST%FileDescLines(1)  = 'Predictions were generated on '//CurDate()//' at '//CurTime()//' using '//TRIM(GetVersion())
   y_FAST%FileDescLines(2)  = 'linked with ' //' '//TRIM(GetNVD(NWTC_Ver            ))  ! we'll get the rest of the linked modules in the section below
   y_FAST%FileDescLines(3)  = 'Description from the FAST input file: '//TRIM(p_FAST%FTitle)
   
   !......................................................
   ! We'll fill out the rest of FileDescLines(2), 
   ! and save the module version info for later use, too:
   !......................................................

   y_FAST%Module_Ver( Module_ED )   = InitOutData_ED%Ver
   y_FAST%FileDescLines(2)          = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_ED )  ))

   IF ( p_FAST%CompElast == Module_BD )  THEN
      y_FAST%Module_Ver( Module_BD ) = InitOutData_BD(1)%Ver ! call copy routine for this type if it every uses dynamic memory     
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_BD ))) 
   END IF   
   
   
   IF ( p_FAST%CompInflow == Module_IfW )  THEN
      y_FAST%Module_Ver( Module_IfW ) = InitOutData_IfW%Ver ! call copy routine for this type if it every uses dynamic memory     
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_IfW ))) 
   ELSEIF ( p_FAST%CompInflow == Module_OpFM )  THEN
      y_FAST%Module_Ver( Module_OpFM ) = InitOutData_OpFM%Ver ! call copy routine for this type if it every uses dynamic memory     
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_OpFM ))) 
   END IF   
   
   IF ( p_FAST%CompAero == Module_AD14 )  THEN
      y_FAST%Module_Ver( Module_AD14  ) = InitOutData_AD14%Ver     
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_AD14  ) ))                  
   ELSEIF ( p_FAST%CompAero == Module_AD )  THEN
      y_FAST%Module_Ver( Module_AD  ) = InitOutData_AD%Ver     
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_AD  ) ))                  
   END IF

   IF ( p_FAST%CompServo == Module_SrvD ) THEN
      y_FAST%Module_Ver( Module_SrvD ) = InitOutData_SrvD%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_SrvD )))
   END IF
         
   IF ( p_FAST%CompHydro == Module_HD ) THEN
      y_FAST%Module_Ver( Module_HD )   = InitOutData_HD%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_HD )))
   END IF

   IF ( p_FAST%CompSub == Module_SD ) THEN
      y_FAST%Module_Ver( Module_SD )   = InitOutData_SD%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_SD )))
   END IF

   IF ( p_FAST%CompMooring == Module_MAP ) THEN
      y_FAST%Module_Ver( Module_MAP )   = InitOutData_MAP%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_MAP )))
   ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
      y_FAST%Module_Ver( Module_MD )   = InitOutData_MD%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_MD )))
   ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
      y_FAST%Module_Ver( Module_FEAM )   = InitOutData_FEAM%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_FEAM )))
   END IF   
   
   IF ( p_FAST%CompIce == Module_IceF ) THEN
      y_FAST%Module_Ver( Module_IceF )   = InitOutData_IceF%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_IceF )))
   ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
      y_FAST%Module_Ver( Module_IceD )   = InitOutData_IceD%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_IceD )))   
   END IF      
   
   !......................................................
   ! Set the number of output columns from each module
   !......................................................
   y_FAST%numOuts = 0    ! Inintialize entire array
   
   
   
   !y_FAST%numOuts(Module_InfW)  = 3  !hack for now: always output 3 wind speeds at hub-height
   IF ( ALLOCATED( InitOutData_IfW%WriteOutputHdr  ) ) y_FAST%numOuts(Module_IfW)  = SIZE(InitOutData_IfW%WriteOutputHdr)
   IF ( ALLOCATED( InitOutData_OpFM%WriteOutputHdr ) ) y_FAST%numOuts(Module_OpFM) = SIZE(InitOutData_OpFM%WriteOutputHdr)
   IF ( ALLOCATED( InitOutData_ED%WriteOutputHdr   ) ) y_FAST%numOuts(Module_ED)   = SIZE(InitOutData_ED%WriteOutputHdr)
do i=1,p_FAST%nBeams
   IF ( ALLOCATED( InitOutData_BD(i)%WriteOutputHdr) ) y_FAST%numOuts(Module_BD)   = y_FAST%numOuts(Module_BD) + SIZE(InitOutData_BD(i)%WriteOutputHdr)
end do   
                                                       y_FAST%numOuts(Module_AD14) = 0
                                                       
   IF ( ALLOCATED( InitOutData_AD%WriteOutputHdr   ) ) y_FAST%numOuts(Module_AD)   = SIZE(InitOutData_AD%WriteOutputHdr)
   IF ( ALLOCATED( InitOutData_SrvD%WriteOutputHdr ) ) y_FAST%numOuts(Module_SrvD) = SIZE(InitOutData_SrvD%WriteOutputHdr)
   IF ( ALLOCATED( InitOutData_HD%WriteOutputHdr   ) ) y_FAST%numOuts(Module_HD)   = SIZE(InitOutData_HD%WriteOutputHdr)
   IF ( ALLOCATED( InitOutData_SD%WriteOutputHdr   ) ) y_FAST%numOuts(Module_SD)   = SIZE(InitOutData_SD%WriteOutputHdr)
   IF ( ALLOCATED( InitOutData_MAP%WriteOutputHdr  ) ) y_FAST%numOuts(Module_MAP)  = SIZE(InitOutData_MAP%WriteOutputHdr)
   IF ( ALLOCATED( InitOutData_FEAM%WriteOutputHdr ) ) y_FAST%numOuts(Module_FEAM) = SIZE(InitOutData_FEAM%WriteOutputHdr)
   IF ( ALLOCATED( InitOutData_MD%WriteOutputHdr   ) ) y_FAST%numOuts(Module_MD)   = SIZE(InitOutData_MD%WriteOutputHdr)
   IF ( ALLOCATED( InitOutData_IceF%WriteOutputHdr ) ) y_FAST%numOuts(Module_IceF) = SIZE(InitOutData_IceF%WriteOutputHdr)
   IF ( ALLOCATED( InitOutData_IceD%WriteOutputHdr ) ) y_FAST%numOuts(Module_IceD) = SIZE(InitOutData_IceD%WriteOutputHdr)*p_FAST%numIceLegs         
   
   !......................................................
   ! Initialize the output channel names and units
   !......................................................
   NumOuts   = 1 + SUM( y_FAST%numOuts )

   CALL AllocAry( y_FAST%ChannelNames,NumOuts, 'ChannelNames', ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( y_FAST%ChannelUnits,NumOuts, 'ChannelUnits', ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN

   y_FAST%ChannelNames(1) = 'Time'
   y_FAST%ChannelUnits(1) = '(s)'

   indxLast = 1
   indxNext = 2

   IF ( y_FAST%numOuts(Module_IfW) > 0_IntKi ) THEN  
      indxLast = indxNext + y_FAST%numOuts(Module_IfW) - 1
      y_FAST%ChannelNames(indxNext:indxLast) = InitOutData_IfW%WriteOutputHdr
      y_FAST%ChannelUnits(indxNext:indxLast) = InitOutData_IfW%WriteOutputUnt      
      indxNext = indxLast + 1
   ELSEIF ( y_FAST%numOuts(Module_OpFM) > 0_IntKi ) THEN  
      indxLast = indxNext + y_FAST%numOuts(Module_OpFM) - 1
      y_FAST%ChannelNames(indxNext:indxLast) = InitOutData_OpFM%WriteOutputHdr
      y_FAST%ChannelUnits(indxNext:indxLast) = InitOutData_OpFM%WriteOutputUnt      
      indxNext = indxLast + 1
   END IF


   IF ( y_FAST%numOuts(Module_ED) > 0_IntKi ) THEN !ElastoDyn
      indxLast = indxNext + y_FAST%numOuts(Module_ED) - 1
      y_FAST%ChannelNames(indxNext:indxLast) = InitOutData_ED%WriteOutputHdr
      y_FAST%ChannelUnits(indxNext:indxLast) = InitOutData_ED%WriteOutputUnt
      indxNext = indxLast + 1
   END IF

   IF ( y_FAST%numOuts(Module_BD) > 0_IntKi ) THEN !BeamDyn
      do i=1,p_FAST%nBeams
         if ( allocated(InitOutData_BD(i)%WriteOutputHdr) ) then            
            do j=1,size(InitOutData_BD(i)%WriteOutputHdr) 
               y_FAST%ChannelNames(indxNext) = 'B'//TRIM(Num2Lstr(i))//trim(InitOutData_BD(i)%WriteOutputHdr(j))
               y_FAST%ChannelUnits(indxNext) = InitOutData_BD(i)%WriteOutputUnt(j)
               indxNext = indxNext + 1
            end do ! j            
         end if         
      end do                 
   END IF
   
   
      ! none for AeroDyn14 
   
   IF ( y_FAST%numOuts(Module_AD) > 0_IntKi ) THEN !AeroDyn
      indxLast = indxNext + y_FAST%numOuts(Module_AD) - 1
      y_FAST%ChannelNames(indxNext:indxLast) = InitOutData_AD%WriteOutputHdr
      y_FAST%ChannelUnits(indxNext:indxLast) = InitOutData_AD%WriteOutputUnt
      indxNext = indxLast + 1
   END IF

   
   IF ( y_FAST%numOuts(Module_SrvD) > 0_IntKi ) THEN !ServoDyn
      indxLast = indxNext + y_FAST%numOuts(Module_SrvD) - 1
      y_FAST%ChannelNames(indxNext:indxLast) = InitOutData_SrvD%WriteOutputHdr
      y_FAST%ChannelUnits(indxNext:indxLast) = InitOutData_SrvD%WriteOutputUnt
      indxNext = indxLast + 1
   END IF

   IF ( y_FAST%numOuts(Module_HD) > 0_IntKi ) THEN !HydroDyn
      indxLast = indxNext + y_FAST%numOuts(Module_HD) - 1
      y_FAST%ChannelNames(indxNext:indxLast) = InitOutData_HD%WriteOutputHdr
      y_FAST%ChannelUnits(indxNext:indxLast) = InitOutData_HD%WriteOutputUnt
      indxNext = indxLast + 1
   END IF

   
   IF ( y_FAST%numOuts(Module_SD) > 0_IntKi ) THEN !SubDyn
      indxLast = indxNext + y_FAST%numOuts(Module_SD) - 1
      y_FAST%ChannelNames(indxNext:indxLast) = InitOutData_SD%WriteOutputHdr
      y_FAST%ChannelUnits(indxNext:indxLast) = InitOutData_SD%WriteOutputUnt
      indxNext = indxLast + 1
   END IF

   
   IF ( y_FAST%numOuts(Module_MAP) > 0_IntKi ) THEN !MAP
      indxLast = indxNext + y_FAST%numOuts(Module_MAP) - 1
      y_FAST%ChannelNames(indxNext:indxLast) = InitOutData_MAP%WriteOutputHdr
      y_FAST%ChannelUnits(indxNext:indxLast) = InitOutData_MAP%WriteOutputUnt
      indxNext = indxLast + 1
   ELSEIF ( y_FAST%numOuts(Module_MD) > 0_IntKi ) THEN !MoorDyn
      indxLast = indxNext + y_FAST%numOuts(Module_MD) - 1
      y_FAST%ChannelNames(indxNext:indxLast) = InitOutData_MD%WriteOutputHdr
      y_FAST%ChannelUnits(indxNext:indxLast) = InitOutData_MD%WriteOutputUnt
      indxNext = indxLast + 1
   ELSEIF ( y_FAST%numOuts(Module_FEAM) > 0_IntKi ) THEN !FEAMooring
      indxLast = indxNext + y_FAST%numOuts(Module_FEAM) - 1
      y_FAST%ChannelNames(indxNext:indxLast) = InitOutData_FEAM%WriteOutputHdr
      y_FAST%ChannelUnits(indxNext:indxLast) = InitOutData_FEAM%WriteOutputUnt
      indxNext = indxLast + 1
   END IF
   
   
   IF ( y_FAST%numOuts(Module_IceF) > 0_IntKi ) THEN !IceFloe
      indxLast = indxNext + y_FAST%numOuts(Module_IceF) - 1
      y_FAST%ChannelNames(indxNext:indxLast) = InitOutData_IceF%WriteOutputHdr
      y_FAST%ChannelUnits(indxNext:indxLast) = InitOutData_IceF%WriteOutputUnt
      indxNext = indxLast + 1
   ELSEIF ( y_FAST%numOuts(Module_IceD) > 0_IntKi ) THEN !IceDyn
      DO I=1,p_FAST%numIceLegs         
         DO J=1,SIZE(InitOutData_IceD%WriteOutputHdr) 
            y_FAST%ChannelNames(indxNext) =TRIM(InitOutData_IceD%WriteOutputHdr(J))//'L'//TRIM(Num2Lstr(I))  !bjj: do we want this "Lx" at the end?
            y_FAST%ChannelUnits(indxNext) = InitOutData_IceD%WriteOutputUnt(J)
            indxNext = indxNext + 1
         END DO ! J
      END DO ! I
   END IF   
      
   
   !......................................................
   ! Open the text output file and print the headers
   !......................................................

   IF (p_FAST%WrTxtOutFile) THEN

      CALL GetNewUnit( y_FAST%UnOu, ErrStat, ErrMsg )
         IF ( ErrStat >= AbortErrLev ) RETURN

      CALL OpenFOutFile ( y_FAST%UnOu, TRIM(p_FAST%OutFileRoot)//'.out', ErrStat, ErrMsg )
         IF ( ErrStat >= AbortErrLev ) RETURN

         ! Add some file information:

      WRITE (y_FAST%UnOu,'(/,A)')  TRIM( y_FAST%FileDescLines(1) )
      WRITE (y_FAST%UnOu,'(1X,A)') TRIM( y_FAST%FileDescLines(2) )
      WRITE (y_FAST%UnOu,'()' )    !print a blank line
      WRITE (y_FAST%UnOu,'(A)'   ) TRIM( y_FAST%FileDescLines(3) )
      WRITE (y_FAST%UnOu,'()' )    !print a blank line


         !......................................................
         ! Write the names of the output parameters on one line:
         !......................................................

      CALL WrFileNR ( y_FAST%UnOu, y_FAST%ChannelNames(1) )

      DO I=2,NumOuts
         CALL WrFileNR ( y_FAST%UnOu, p_FAST%Delim//y_FAST%ChannelNames(I) )
      ENDDO ! I

      WRITE (y_FAST%UnOu,'()')

         !......................................................
         ! Write the units of the output parameters on one line:
         !......................................................

      CALL WrFileNR ( y_FAST%UnOu, y_FAST%ChannelUnits(1) )

      DO I=2,NumOuts
         CALL WrFileNR ( y_FAST%UnOu, p_FAST%Delim//y_FAST%ChannelUnits(I) )
      ENDDO ! I

      WRITE (y_FAST%UnOu,'()')

   END IF

   !......................................................
   ! Allocate data for binary output file
   !......................................................
   IF (p_FAST%WrBinOutFile) THEN

         ! calculate the size of the array of outputs we need to store
      y_FAST%NOutSteps = NINT ( (p_FAST%TMax - p_FAST%TStart) / p_FAST%DT_OUT ) + 1

      CALL AllocAry( y_FAST%AllOutData, NumOuts-1, y_FAST%NOutSteps, 'AllOutData', ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) RETURN

      IF ( OutputFileFmtID == FileFmtID_WithoutTime ) THEN

         CALL AllocAry( y_FAST%TimeData, 2_IntKi, 'TimeData', ErrStat, ErrMsg )
         IF ( ErrStat >= AbortErrLev ) RETURN

         y_FAST%TimeData(1) = 0.0_DbKi           ! This is the first output time, which we will set later
         y_FAST%TimeData(2) = p_FAST%DT_out      ! This is the (constant) time between subsequent writes to the output file

      ELSE  ! we store the entire time array

         CALL AllocAry( y_FAST%TimeData, y_FAST%NOutSteps, 'TimeData', ErrStat, ErrMsg )
         IF ( ErrStat >= AbortErrLev ) RETURN

      END IF

      y_FAST%n_Out = 0  !number of steps actually written to the file

   END IF


RETURN
END SUBROUTINE FAST_InitOutput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FAST_WrSum( p_FAST, y_FAST, MeshMapData, ErrStat, ErrMsg )
! This subroutine opens and writes data to the FAST summary file. The file gets closed at the end of program (not in this 
! subroutine).
!..................................................................................................................................

   TYPE(FAST_ParameterType), INTENT(IN)    :: p_FAST                             ! Glue-code simulation parameters
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST                             ! Glue-code simulation outputs (changes value of UnSum)
   TYPE(FAST_ModuleMapType), INTENT(IN)    :: MeshMapData                        ! Data for mapping between modules
   INTEGER(IntKi),           INTENT(OUT)   :: ErrStat                            ! Error status (level)
   CHARACTER(*),             INTENT(OUT)   :: ErrMsg                             ! Message describing error reported in ErrStat

      ! local variables
   INTEGER(IntKi)                          :: I                                  ! temporary counter
   INTEGER(IntKi)                          :: J                                  ! temporary counter
   INTEGER(IntKi)                          :: Module_Number                      ! loop counter through the modules
   CHARACTER(200)                          :: Fmt                                ! temporary format string
   CHARACTER(200)                          :: DescStr                            ! temporary string to write text
   CHARACTER(*), PARAMETER                 :: NotUsedTxt = " [not called]"       ! text written if a module is not called

      ! Get a unit number and open the file:

   CALL GetNewUnit( y_FAST%UnSum, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) RETURN

   CALL OpenFOutFile ( y_FAST%UnSum, TRIM(p_FAST%OutFileRoot)//'.sum', ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) RETURN

         ! Add some file information:

   !.......................... Module Versions .....................................................
   !bjj: modules in this list are ordered by the order they are specified in the FAST input file

   WRITE (y_FAST%UnSum,'(/A)') 'FAST Summary File'
   WRITE (y_FAST%UnSum,'(/A)')  TRIM( y_FAST%FileDescLines(1) )

   WRITE (y_FAST%UnSum,'(2X,A)'   )  'compiled with'
   Fmt = '(4x,A)'
   WRITE (y_FAST%UnSum,Fmt)  TRIM( GetNVD(        NWTC_Ver ) )
   WRITE (y_FAST%UnSum,Fmt)  TRIM( GetNVD( y_FAST%Module_Ver( Module_ED )   ) )

   DescStr = GetNVD( y_FAST%Module_Ver( Module_BD ) )
   IF ( p_FAST%CompElast /= Module_BD ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )
      
   DescStr = GetNVD( y_FAST%Module_Ver( Module_IfW ) )
   IF ( p_FAST%CompInflow /= Module_IfW ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )
   
   ! I'm not going to write the openfoam module info to the summary file
   !DescStr = GetNVD( y_FAST%Module_Ver( Module_OpFM ) )
   !IF ( p_FAST%CompInflow /= Module_OpFM ) DescStr = TRIM(DescStr)//NotUsedTxt
   !WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )
      
   DescStr = GetNVD( y_FAST%Module_Ver( Module_AD14 ) )
   IF ( p_FAST%CompAero /= Module_AD14 ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )
      
   DescStr = GetNVD( y_FAST%Module_Ver( Module_AD ) )
   IF ( p_FAST%CompAero /= Module_AD ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )
     
   DescStr = GetNVD( y_FAST%Module_Ver( Module_SrvD ) )
   IF ( p_FAST%CompServo /= Module_SrvD ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )  
   
   DescStr = GetNVD( y_FAST%Module_Ver( Module_HD ) )
   IF ( p_FAST%CompHydro /= Module_HD  ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )
   
   DescStr = GetNVD( y_FAST%Module_Ver( Module_SD ) )
   IF ( p_FAST%CompSub /= Module_SD ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )
   
   DescStr = GetNVD( y_FAST%Module_Ver( Module_MAP ) )
   IF ( p_FAST%CompMooring /= Module_MAP ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )

   DescStr = GetNVD( y_FAST%Module_Ver( Module_FEAM ) )
   IF ( p_FAST%CompMooring /= Module_FEAM ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )
   
   DescStr = GetNVD( y_FAST%Module_Ver( Module_MD ) )
   IF ( p_FAST%CompMooring /= Module_MD ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )
   
   DescStr = GetNVD( y_FAST%Module_Ver( Module_IceF ) )
   IF ( p_FAST%CompIce /= Module_IceF ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )
   
   DescStr = GetNVD( y_FAST%Module_Ver( Module_IceD ) )
   IF ( p_FAST%CompIce /= Module_IceD ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )
   
   
   !.......................... Information from FAST input File ......................................
! OTHER information we could print here:   
! current working directory
! output file root name
! output file time step
! output file format (text/binary)
! coupling method

   SELECT CASE ( p_FAST%TurbineType )
   CASE ( Type_LandBased )
      DescStr = 'Modeling a land-based turbine'
   CASE ( Type_Offshore_Fixed )
      DescStr = 'Modeling a fixed-bottom offshore turbine'
   CASE ( Type_Offshore_Floating )
      DescStr = 'Modeling a floating offshore turbine'
   CASE DEFAULT ! This should never happen
      DescStr=""
   END SELECT                  
   WRITE(y_FAST%UnSum,'(//A)') TRIM(DescStr)

   WRITE (y_FAST%UnSum,'(A)' )   'Description from the FAST input file: '
   WRITE (y_FAST%UnSum,'(2X,A)')  TRIM(p_FAST%FTitle)

   !.......................... Requested Features ...................................................
   
   SELECT CASE ( p_FAST%InterpOrder )
   CASE (0)
      DescStr = ' (nearest neighbor)'
   CASE (1)
      DescStr = ' (linear)'
   CASE (2)
      DescStr = ' (quadratic)'
   CASE DEFAULT 
      DescStr = ' ( )'
   END SELECT               
   
   WRITE(y_FAST%UnSum,'(/A,I1,A)'  ) 'Interpolation order for input/output time histories: ', p_FAST%InterpOrder, TRIM(DescStr)
   WRITE(y_FAST%UnSum,'( A,I2)'    ) 'Number of correction iterations: ', p_FAST%NumCrctn
   
      
   !.......................... Information About Coupling ...................................................
      
   IF ( ALLOCATED( MeshMapData%Jacobian_ED_SD_HD_BD ) ) then !ED-SD-HD-BD
      
      IF ( p_FAST%CompSub == Module_SD .OR. p_FAST%CompElast == Module_BD ) THEN  ! SubDyn-BeamDyn-HydroDyn-ElastoDyn
         DescStr = 'ElastoDyn, SubDyn, HydroDyn, and/or BeamDyn'                  
      ELSE ! IF ( p_FAST%CompHydro == Module_HD ) THEN
         DescStr = "ElastoDyn to HydroDyn"
      END IF
                  
      WRITE(y_FAST%UnSum,'( A,I6)'  ) 'Number of rows in Jacobian matrix used for coupling '//TRIM(DescStr)//': ', &
                                       SIZE(MeshMapData%Jacobian_ED_SD_HD_BD, 1)
   END IF

   !.......................... Time step information: ...................................................
   
   WRITE (y_FAST%UnSum,'(//,2X,A)') " Requested Time Steps  "
   WRITE (y_FAST%UnSum,   '(2X,A)') "-------------------------------------------------"
   Fmt = '(2X,A17,2X,A15,2X,A13)'
   WRITE (y_FAST%UnSum, Fmt ) "Component        ", "Time Step (s)  ", "Subcycles (-)"
   WRITE (y_FAST%UnSum, Fmt ) "-----------------", "---------------", "-------------"
   Fmt = '(2X,A17,2X,'//TRIM(p_FAST%OutFmt)//',:,T37,2X,I8,:,A)'
   WRITE (y_FAST%UnSum, Fmt ) "FAST (glue code) ", p_FAST%DT
   DO Module_Number=1,NumModules
      IF (p_FAST%ModuleInitialized(Module_Number)) THEN
         WRITE (y_FAST%UnSum, Fmt ) y_FAST%Module_Ver(Module_Number)%Name, p_FAST%DT_module(Module_Number), p_FAST%n_substeps(Module_Number)
      END IF
   END DO
   IF ( NINT( p_FAST%DT_out / p_FAST%DT )  == 1_IntKi ) THEN
      WRITE (y_FAST%UnSum, Fmt ) "FAST output files", p_FAST%DT_out, 1_IntKi   ! we'll write "1" instead of "1^-1"
   ELSE
      WRITE (y_FAST%UnSum, Fmt ) "FAST output files", p_FAST%DT_out, NINT( p_FAST%DT_out / p_FAST%DT ),"^-1"
   END IF
   
   
   !.......................... Requested Output Channels ............................................

   WRITE (y_FAST%UnSum,'(//,2X,A)') " Requested Channels in FAST Output File(s)  "
   WRITE (y_FAST%UnSum,   '(2X,A)') "--------------------------------------------"
   Fmt = '(2X,A6,2(2X,A'//TRIM(num2lstr(ChanLen))//'),2X,A)'
   WRITE (y_FAST%UnSum, Fmt ) "Number", "Name      ", "Units     ", "Generated by"
   WRITE (y_FAST%UnSum, Fmt ) "------", "----------", "----------", "------------"

   Fmt = '(4X,I4,2(2X,A'//TRIM(num2lstr(ChanLen))//'),2X,A)'
   I = 1
   WRITE (y_FAST%UnSum, Fmt ) I, y_FAST%ChannelNames(I), y_FAST%ChannelUnits(I), TRIM(FAST_Ver%Name)

   
   DO Module_Number = 1,NumModules
      DO J = 1,y_FAST%numOuts( Module_Number )
         I = I + 1
         WRITE (y_FAST%UnSum, Fmt ) I, y_FAST%ChannelNames(I), y_FAST%ChannelUnits(I), TRIM(y_FAST%Module_Ver( Module_Number )%Name)
      END DO
   END DO
      
   
   !.......................... End of Summary File ............................................
   
   ! bjj: note that I'm not closing the summary file here, though at the present time we don't write to this file again.
   ! In the future, we may want to write additional information to this file during the simulation.
   ! bjj 4/21/2015: closing the file now because of restart. If it needs to be open later, we can change it again.
   
   CLOSE( y_FAST%UnSum )        
   y_FAST%UnSum = -1

END SUBROUTINE FAST_WrSum
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FAST_ReadPrimaryFile( InputFile, p, ErrStat, ErrMsg )
! This routine reads in the primary FAST input file, does some validation, and places the values it reads in the
!   parameter structure (p). It prints to an echo file if requested.
!..................................................................................................................................


   IMPLICIT                        NONE

      ! Passed variables
   TYPE(FAST_ParameterType), INTENT(INOUT) :: p                               ! The parameter data for the FAST (glue-code) simulation
   INTEGER(IntKi),           INTENT(OUT)   :: ErrStat                         ! Error status

   CHARACTER(*),             INTENT(IN)    :: InputFile                       ! Name of the file containing the primary input data
   CHARACTER(*),             INTENT(OUT)   :: ErrMsg                          ! Error message

      ! Local variables:
   REAL(DbKi)                    :: TmpTime                                   ! temporary variable to read SttsTime and ChkptTime before converting to #steps based on DT
   INTEGER(IntKi)                :: I                                         ! loop counter
   INTEGER(IntKi)                :: UnIn                                      ! Unit number for reading file
   INTEGER(IntKi)                :: UnEc                                      ! I/O unit for echo file. If > 0, file is open for writing.

   INTEGER(IntKi)                :: IOS                                       ! Temporary Error status
   INTEGER(IntKi)                :: ErrStat2                                  ! Temporary Error status
   INTEGER(IntKi)                :: FmtWidth                                  ! width of the field returned by the specified OutFmt
   INTEGER(IntKi)                :: OutFileFmt                                ! An integer that indicates what kind of tabular output should be generated (1=text, 2=binary, 3=both)
   LOGICAL                       :: Echo                                      ! Determines if an echo file should be written
   LOGICAL                       :: TabDelim                                  ! Determines if text output should be delimited by tabs (true) or space (false)
   CHARACTER(ErrMsgLen)          :: ErrMsg2                                   ! Temporary Error message
   CHARACTER(1024)               :: PriPath                                   ! Path name of the primary file

   CHARACTER(10)                 :: AbortLevel                                ! String that indicates which error level should be used to abort the program: WARNING, SEVERE, or FATAL
   CHARACTER(30)                 :: Line                                      ! string for default entry in input file


      ! Initialize some variables:
   UnEc = -1
   Echo = .FALSE.                        ! Don't echo until we've read the "Echo" flag
   CALL GetPath( InputFile, PriPath )    ! Input files will be relative to the path where the primary input file is located.


      ! Get an available unit number for the file.

   CALL GetNewUnit( UnIn, ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN


      ! Open the Primary input file.

   CALL OpenFInpFile ( UnIn, InputFile, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


   ! Read the lines up/including to the "Echo" simulation control variable
   ! If echo is FALSE, don't write these lines to the echo file.
   ! If Echo is TRUE, rewind and write on the second try.

   I = 1 !set the number of times we've read the file
   DO
   !-------------------------- HEADER ---------------------------------------------

      CALL ReadCom( UnIn, InputFile, 'File header: Module Version (line 1)', ErrStat2, ErrMsg2, UnEc )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF ( ErrStat >= AbortErrLev ) RETURN

      CALL ReadStr( UnIn, InputFile, p%FTitle, 'FTitle', 'File Header: File Description (line 2)', ErrStat2, ErrMsg2, UnEc )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF ( ErrStat >= AbortErrLev ) RETURN


   !---------------------- SIMULATION CONTROL --------------------------------------
      CALL ReadCom( UnIn, InputFile, 'Section Header: Simulation Control', ErrStat2, ErrMsg2, UnEc )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF ( ErrStat >= AbortErrLev ) RETURN


         ! Echo - Echo input data to <RootName>.ech (flag):
      CALL ReadVar( UnIn, InputFile, Echo, "Echo", "Echo input data to <RootName>.ech (flag)", ErrStat2, ErrMsg2, UnEc)
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF ( ErrStat >= AbortErrLev ) RETURN


      IF (.NOT. Echo .OR. I > 1) EXIT !exit this loop

         ! Otherwise, open the echo file, then rewind the input file and echo everything we've read

      I = I + 1         ! make sure we do this only once (increment counter that says how many times we've read this file)

      CALL OpenEcho ( UnEc, TRIM(p%OutFileRoot)//'.ech', ErrStat2, ErrMsg2, FAST_Ver )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF ( ErrStat >= AbortErrLev ) RETURN

      IF ( UnEc > 0 )  WRITE (UnEc,'(/,A,/)')  'Data from '//TRIM(FAST_Ver%Name)//' primary input file "'//TRIM( InputFile )//'":'

      REWIND( UnIn, IOSTAT=ErrStat2 )
         IF (ErrStat2 /= 0_IntKi ) THEN
            CALL CheckError( ErrID_Fatal, 'Error rewinding file "'//TRIM(InputFile)//'".' )
            IF ( ErrStat >= AbortErrLev ) RETURN
         END IF

   END DO

   CALL WrScr( ' Heading of the '//TRIM(FAST_Ver%Name)//' input file: ' )
   CALL WrScr( '   '//TRIM( p%FTitle ) )


      ! AbortLevel - Error level when simulation should abort:
   CALL ReadVar( UnIn, InputFile, AbortLevel, "AbortLevel", "Error level when simulation should abort (string)", &
                        ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

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
            CALL CheckError( ErrID_Fatal, 'Invalid AbortLevel specified in FAST input file. '// &
                                 'Valid entries are "WARNING", "SEVERE", or "FATAL".' )
            RETURN
      END SELECT


      ! TMax - Total run time (s):
   CALL ReadVar( UnIn, InputFile, p%TMax, "TMax", "Total run time (s)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! DT - Recommended module time step (s):
   CALL ReadVar( UnIn, InputFile, p%DT, "DT", "Recommended module time step (s)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! InterpOrder - Interpolation order for inputs and outputs {0=nearest neighbor ,1=linear, 2=quadratic}
   CALL ReadVar( UnIn, InputFile, p%InterpOrder, "InterpOrder", "Interpolation order "//&
                   "for inputs and outputs {0=nearest neighbor ,1=linear, 2=quadratic} (-)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! NumCrctn - Number of predictor-corrector iterations {1=explicit calculation, i.e., no corrections}
   CALL ReadVar( UnIn, InputFile, p%NumCrctn, "NumCrctn", "Number of corrections"//&
                   "{0=explicit calculation, i.e., no corrections} (-)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! DT_UJac - Time between calls to get Jacobians (s)
   CALL ReadVar( UnIn, InputFile, p%DT_UJac, "DT_UJac", "Time between calls to get Jacobians (s)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! UJacSclFact - Scaling factor used in Jacobians (-)
   CALL ReadVar( UnIn, InputFile, p%UJacSclFact, "UJacSclFact", "Scaling factor used in Jacobians (-)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN              
                  
   !---------------------- FEATURE SWITCHES AND FLAGS --------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Feature Switches and Flags', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! CompElast - Compute structural dynamics (switch) {1=ElastoDyn; 2=ElastoDyn + BeamDyn for blades}:
   CALL ReadVar( UnIn, InputFile, p%CompElast, "CompElast", "Compute structural dynamics (switch) {1=ElastoDyn; 2=ElastoDyn + BeamDyn for blades}", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
          ! immediately convert to values used inside the code:
         IF ( p%CompElast == 1 ) THEN 
            p%CompElast = Module_ED
         ELSEIF ( p%CompElast == 2 ) THEN
            p%CompElast = Module_BD
         ELSE
            p%CompElast = Module_Unknown
         END IF
                                  
      ! CompInflow - inflow wind velocities (switch) {0=still air; 1=InflowWind}:
   CALL ReadVar( UnIn, InputFile, p%CompInflow, "CompInflow", "inflow wind velocities (switch) {0=still air; 1=InflowWind}", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
          ! immediately convert to values used inside the code:
         IF ( p%CompInflow == 0 ) THEN 
            p%CompInflow = Module_NONE
         ELSEIF ( p%CompInflow == 1 ) THEN
            p%CompInflow = Module_IfW
         ELSEIF ( p%CompInflow == 2 ) THEN
            p%CompInflow = Module_OpFM
         ELSE
            p%CompInflow = Module_Unknown
         END IF

      ! CompAero - Compute aerodynamic loads (switch) {0=None; 1=AeroDyn}:
   CALL ReadVar( UnIn, InputFile, p%CompAero, "CompAero", "Compute aerodynamic loads (switch) {0=None; 1=AeroDyn}", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
          ! immediately convert to values used inside the code:
         IF ( p%CompAero == 0 ) THEN 
            p%CompAero = Module_NONE
         ELSEIF ( p%CompAero == 1 ) THEN
            p%CompAero = Module_AD14
         ELSEIF ( p%CompAero == 2 ) THEN
            p%CompAero = Module_AD
         ELSE
            p%CompAero = Module_Unknown
         END IF

      ! CompServo - Compute control and electrical-drive dynamics (switch) {0=None; 1=ServoDyn}:
   CALL ReadVar( UnIn, InputFile, p%CompServo, "CompServo", "Compute control and electrical-drive dynamics (switch) {0=None; 1=ServoDyn}", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
          ! immediately convert to values used inside the code:
         IF ( p%CompServo == 0 ) THEN 
            p%CompServo = Module_NONE
         ELSEIF ( p%CompServo == 1 ) THEN
            p%CompServo = Module_SrvD
         ELSE
            p%CompServo = Module_Unknown
         END IF
      
      
      ! CompHydro - Compute hydrodynamic loads (switch) {0=None; 1=HydroDyn}:
   CALL ReadVar( UnIn, InputFile, p%CompHydro, "CompHydro", "Compute hydrodynamic loads (switch) {0=None; 1=HydroDyn}", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
          ! immediately convert to values used inside the code:
         IF ( p%CompHydro == 0 ) THEN 
            p%CompHydro = Module_NONE
         ELSEIF ( p%CompHydro == 1 ) THEN
            p%CompHydro = Module_HD
         ELSE
            p%CompHydro = Module_Unknown
         END IF
         
      ! CompSub - Compute sub-structural dynamics (switch) {0=None; 1=SubDyn}:
   CALL ReadVar( UnIn, InputFile, p%CompSub, "CompSub", "Compute sub-structural dynamics (switch) {0=None; 1=SubDyn}", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
          ! immediately convert to values used inside the code:
         IF ( p%CompSub == 0 ) THEN 
            p%CompSub = Module_NONE
         ELSEIF ( p%CompSub == 1 ) THEN
            p%CompSub = Module_SD
         ELSE
            p%CompSub = Module_Unknown
         END IF
         
      ! CompMooring - Compute mooring line dynamics (flag):
   CALL ReadVar( UnIn, InputFile, p%CompMooring, "CompMooring", "Compute mooring system (switch) {0=None; 1=MAP, 2=FEAMooring}", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN         
          ! immediately convert to values used inside the code:
         IF ( p%CompMooring == 0 ) THEN 
            p%CompMooring = Module_NONE
         ELSEIF ( p%CompMooring == 1 ) THEN
            p%CompMooring = Module_MAP
         ELSEIF ( p%CompMooring == 2 ) THEN
            p%CompMooring = Module_FEAM
         ELSEIF ( p%CompMooring == 3 ) THEN
            p%CompMooring = Module_MD
         ELSE
            p%CompMooring = Module_Unknown
         END IF      
      
      ! CompIce - Compute ice loads (switch) {0=None; 1=IceFloe}:
   CALL ReadVar( UnIn, InputFile, p%CompIce, "CompIce", "Compute ice loads (switch) {0=None; 1=IceFloe}", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
          ! immediately convert to values used inside the code:
         IF ( p%CompIce == 0 ) THEN 
            p%CompIce = Module_NONE
         ELSEIF ( p%CompIce == 1 ) THEN
            p%CompIce = Module_IceF
         ELSEIF ( p%CompIce == 2 ) THEN
            p%CompIce = Module_IceD
         ELSE
            p%CompIce = Module_Unknown
         END IF
         
         
      ! CompUserPtfmLd - Compute additional platform loading {false: none, true: user-defined from routine UserPtfmLd} (flag):
   CALL ReadVar( UnIn, InputFile, p%CompUserPtfmLd, "CompUserPtfmLd", "Compute additional platform loading (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! CompUserTwrLd - Compute additional tower loading {false: none, true: user-defined from routine UserTwrLd} (flag):
   CALL ReadVar( UnIn, InputFile, p%CompUserTwrLd, "CompUserTwrLd", "Compute additional tower loading (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   !---------------------- INPUT FILES ---------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Input Files', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! EDFile - Name of file containing ElastoDyn input parameters (-):
   CALL ReadVar( UnIn, InputFile, p%EDFile, "EDFile", "Name of file containing ElastoDyn input parameters (-)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   IF ( PathIsRelative( p%EDFile ) ) p%EDFile = TRIM(PriPath)//TRIM(p%EDFile)

DO i=1,MaxNBlades
      ! BDBldFile - Name of file containing BeamDyn blade input parameters (-):
   CALL ReadVar( UnIn, InputFile, p%BDBldFile(i), "BDBldFile("//TRIM(num2LStr(i))//")", "Name of file containing BeamDyn blade "//trim(num2lstr(i))//"input parameters (-)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   IF ( PathIsRelative( p%BDBldFile(i) ) ) p%BDBldFile(i) = TRIM(PriPath)//TRIM(p%BDBldFile(i))
END DO
   
      ! InflowFile - Name of file containing inflow wind input parameters (-):
   CALL ReadVar( UnIn, InputFile, p%InflowFile, "InflowFile", "Name of file containing inflow wind input parameters (-)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   IF ( PathIsRelative( p%InflowFile ) ) p%InflowFile = TRIM(PriPath)//TRIM(p%InflowFile)

      ! AeroFile - Name of file containing aerodynamic input parameters (-):
   CALL ReadVar( UnIn, InputFile, p%AeroFile, "AeroFile", "Name of file containing aerodynamic input parameters (-)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   IF ( PathIsRelative( p%AeroFile ) ) p%AeroFile = TRIM(PriPath)//TRIM(p%AeroFile)

      ! ServoFile - Name of file containing control and electrical-drive input parameters (-):
   CALL ReadVar( UnIn, InputFile, p%ServoFile, "ServoFile", "Name of file containing control and electrical-drive input parameters (-)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   IF ( PathIsRelative( p%ServoFile ) ) p%ServoFile = TRIM(PriPath)//TRIM(p%ServoFile)

      ! HydroFile - Name of file containing hydrodynamic input parameters (-):
   CALL ReadVar( UnIn, InputFile, p%HydroFile, "HydroFile", "Name of file containing hydrodynamic input parameters (-)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   IF ( PathIsRelative( p%HydroFile ) ) p%HydroFile = TRIM(PriPath)//TRIM(p%HydroFile)

      ! SubFile - Name of file containing sub-structural input parameters (-):
   CALL ReadVar( UnIn, InputFile, p%SubFile, "SubFile", "Name of file containing sub-structural input parameters (-)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   IF ( PathIsRelative( p%SubFile ) ) p%SubFile = TRIM(PriPath)//TRIM(p%SubFile)

      ! MooringFile - Name of file containing mooring system input parameters (-):
   CALL ReadVar( UnIn, InputFile, p%MooringFile, "MooringFile", "Name of file containing mooring system input parameters (-)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   IF ( PathIsRelative( p%MooringFile ) ) p%MooringFile = TRIM(PriPath)//TRIM(p%MooringFile)
  
      ! IceFile - Name of file containing ice input parameters (-):
   CALL ReadVar( UnIn, InputFile, p%IceFile, "IceFile", "Name of file containing ice input parameters (-)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   IF ( PathIsRelative( p%IceFile ) ) p%IceFile = TRIM(PriPath)//TRIM(p%IceFile)
   
   
   !---------------------- OUTPUT --------------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Output', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! SumPrint - Print summary data to <RootName>.sum (flag):
   CALL ReadVar( UnIn, InputFile, p%SumPrint, "SumPrint", "Print summary data to <RootName>.sum (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! SttsTime - Amount of time between screen status messages (s):
   CALL ReadVar( UnIn, InputFile, TmpTime, "SttsTime", "Amount of time between screen status messages (s)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN   
      IF (TmpTime > p%TMax) THEN
         p%n_SttsTime = HUGE(p%n_SttsTime)
      ELSE         
         p%n_SttsTime = NINT( TmpTime / p%DT )
      END IF
      
      p%n_SttsTime = NINT( TmpTime / p%DT )

      ! ChkptTime - Amount of time between creating checkpoint files for potential restart (s):
   CALL ReadVar( UnIn, InputFile, TmpTime, "ChkptTime", "Amount of time between creating checkpoint files for potential restart (s)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN        
      IF (TmpTime > p%TMax) THEN
         p%n_ChkptTime = HUGE(p%n_ChkptTime)
      ELSE         
         p%n_ChkptTime = NINT( TmpTime / p%DT )
      END IF

      ! DT_Out - Time step for tabular output (s):
   CALL ReadVar( UnIn, InputFile, Line, "DT_Out", "Time step for tabular output (s)", ErrStat2, ErrMsg2, UnEc)
   !CALL ReadVar( UnIn, InputFile, p%DT_Out, "DT_Out", "Time step for tabular output (s)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
      CALL Conv2UC( Line )
      IF ( INDEX(Line, "DEFAULT" ) == 1 ) THEN 
         p%DT_Out = p%DT
      ELSE
         ! If it's not "default", read this variable; otherwise use the value already stored in InputFileData%DTBeam
         READ( Line, *, IOSTAT=IOS) p%DT_Out
            CALL CheckIOS ( IOS, InputFile, 'DT_Out', NumType, ErrStat2, ErrMsg2 )
            CALL CheckError( ErrStat2, ErrMsg2 )
      END IF

      
      
      ! TStart - Time to begin tabular output (s):
   CALL ReadVar( UnIn, InputFile, p%TStart, "TStart", "Time to begin tabular output (s)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! OutFileFmt - Format for tabular (time-marching) output file(s) (1: text file [<RootName>.out], 2: binary file [<RootName>.outb], 3: both) (-):
   CALL ReadVar( UnIn, InputFile, OutFileFmt, "OutFileFmt", "Format for tabular (time-marching) output file(s) (1: text file [<RootName>.out], 2: binary file [<RootName>.outb], 3: both) (-)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

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
           CALL CheckError( ErrID_Fatal, " FAST's OutFileFmt must be 1, 2, or 3." )
           RETURN
      END SELECT

      ! TabDelim - Use tab delimiters in text tabular output file? (flag):
   CALL ReadVar( UnIn, InputFile, TabDelim, "TabDelim", "Use tab delimiters in text tabular output file? (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      IF ( TabDelim ) THEN
         p%Delim = TAB
      ELSE
         p%Delim = ' '
      END IF


      ! OutFmt - Format used for text tabular output (except time).  Resulting field should be 10 characters. (-):
   CALL ReadVar( UnIn, InputFile, p%OutFmt, "OutFmt", "Format used for text tabular output (except time).  Resulting field should be 10 characters. (-)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! Check that InputFileData%OutFmt is a valid format specifier and will fit over the column headings
   CALL ChkRealFmtStr( p%OutFmt, 'OutFmt', FmtWidth, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      IF ( FmtWidth /= ChanLen ) CALL CheckError( ErrID_Warn, 'OutFmt produces a column width of '// &
            TRIM(Num2LStr(FmtWidth))//' instead of '//TRIM(Num2LStr(ChanLen))//' characters.' )
      IF ( ErrStat >= AbortErrLev ) RETURN

   !---------------------- END OF FILE -----------------------------------------

   CLOSE ( UnIn )
   IF ( UnEc > 0 ) CLOSE ( UnEc )
   RETURN


CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)

      CHARACTER(*),   PARAMETER  :: RoutineName = 'FAST_ReadPrimaryFile'

      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      CALL SetErrStat( ErrID, Msg, ErrStat, ErrMsg, RoutineName )      

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close file, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
            CLOSE( UnIn )
            IF ( UnEc > 0 ) CLOSE ( UnEc )
         END IF



   END SUBROUTINE CheckError
   !...............................................................................................................................
END SUBROUTINE FAST_ReadPrimaryFile
!----------------------------------------------------------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE WrOutputLine( t, p_FAST, y_FAST, IfWOutput, OpFMOutput, EDOutput, ADOutput, SrvDOutput, HDOutput, SDOutput, MAPOutput, FEAMOutput, &
                        MDOutput, IceFOutput, y_IceD, y_BD, ErrStat, ErrMsg)
! This routine writes the module output to the primary output file(s).
!..................................................................................................................................

   IMPLICIT                        NONE

!bjj: add Module_AD here
   
      ! Passed variables
   REAL(DbKi), INTENT(IN)                  :: t                                  ! Current simulation time, in seconds
   TYPE(FAST_ParameterType), INTENT(IN)    :: p_FAST                             ! Glue-code simulation parameters
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST                             ! Glue-code simulation outputs


   REAL(ReKi),               INTENT(IN)    :: IfWOutput (:)                      ! InflowWind WriteOutput values
   REAL(ReKi),               INTENT(IN)    :: OpFMOutput (:)                     ! OpenFOAM WriteOutput values
   REAL(ReKi),               INTENT(IN)    :: EDOutput (:)                       ! ElastoDyn WriteOutput values
   REAL(ReKi),               INTENT(IN)    :: ADOutput (:)                       ! AeroDyn WriteOutput values
   REAL(ReKi),               INTENT(IN)    :: SrvDOutput (:)                     ! ServoDyn WriteOutput values
   REAL(ReKi),               INTENT(IN)    :: HDOutput (:)                       ! HydroDyn WriteOutput values
   REAL(ReKi),               INTENT(IN)    :: SDOutput (:)                       ! SubDyn WriteOutput values
   REAL(ReKi),               INTENT(IN)    :: MAPOutput (:)                      ! MAP WriteOutput values
   REAL(ReKi),               INTENT(IN)    :: FEAMOutput (:)                     ! FEAMooring WriteOutput values
   REAL(ReKi),               INTENT(IN)    :: MDOutput (:)                       ! MoorDyn WriteOutput values
   REAL(ReKi),               INTENT(IN)    :: IceFOutput (:)                     ! IceFloe WriteOutput values
   TYPE(IceD_OutputType),    INTENT(IN)    :: y_IceD (:)                         ! IceDyn outputs (WriteOutput values are subset)
   TYPE(BD_OutputType),      INTENT(IN)    :: y_BD (:)                           ! BeamDyn outputs (WriteOutput values are subset)

   INTEGER(IntKi),           INTENT(OUT)   :: ErrStat
   CHARACTER(*),             INTENT(OUT)   :: ErrMsg

      ! Local variables.

   CHARACTER(200)                   :: Frmt                                      ! A string to hold a format specifier
   CHARACTER(ChanLen)               :: TmpStr                                    ! temporary string to print the time output as text

   REAL(ReKi)                       :: OutputAry(SIZE(y_FAST%ChannelNames)-1)

   ErrStat = ErrID_None
   ErrMsg  = ''
   
   CALL FillOutputAry(p_FAST, y_FAST, IfWOutput, OpFMOutput, EDOutput, ADOutput, SrvDOutput, HDOutput, SDOutput, MAPOutput, FEAMOutput, &
                           MDOutput, IceFOutput, y_IceD, y_BD, OutputAry)   

   IF (p_FAST%WrTxtOutFile) THEN

         ! Write one line of tabular output:
   !   Frmt = '(F8.3,'//TRIM(Num2LStr(p%NumOuts))//'(:,A,'//TRIM( p%OutFmt )//'))'
      Frmt = '"'//p_FAST%Delim//'"'//p_FAST%OutFmt      ! format for array elements from individual modules

            ! time
      WRITE( TmpStr, '(F10.4)' ) t
      CALL WrFileNR( y_FAST%UnOu, TmpStr )

         ! write the individual module output
      CALL WrReAryFileNR ( y_FAST%UnOu, OutputAry,   Frmt, ErrStat, ErrMsg )
         !IF ( ErrStat >= AbortErrLev ) RETURN
      
         ! write a new line (advance to the next line)
      WRITE (y_FAST%UnOu,'()')

   END IF


   IF (p_FAST%WrBinOutFile) THEN

         ! Write data to array for binary output file

      IF ( y_FAST%n_Out == y_FAST%NOutSteps ) THEN
         CALL ProgWarn( 'Not all data could be written to the binary output file.' )
         !this really would only happen if we have an error somewhere else, right?
         !otherwise, we could allocate a new, larger array and move existing data
      ELSE
         y_FAST%n_Out = y_FAST%n_Out + 1

            ! store time data
         IF ( y_FAST%n_Out == 1_IntKi .OR. OutputFileFmtID == FileFmtID_WithTime ) THEN
            y_FAST%TimeData(y_FAST%n_Out) = t   ! Time associated with these outputs
         END IF

            ! store individual module data
         y_FAST%AllOutData(:, y_FAST%n_Out) = OutputAry
         
      END IF      

   END IF

   RETURN
END SUBROUTINE WrOutputLine
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FillOutputAry_T(Turbine, Outputs)
                   
   TYPE(FAST_TurbineType),   INTENT(IN   ) :: Turbine
   REAL(ReKi),               INTENT(  OUT) :: Outputs(:)                       ! single array of output 
   

      CALL FillOutputAry(Turbine%p_FAST, Turbine%y_FAST, Turbine%IfW%y%WriteOutput, Turbine%OpFM%y%WriteOutput, &
                Turbine%ED%Output(1)%WriteOutput, Turbine%AD%y%WriteOutput, Turbine%SrvD%y%WriteOutput, &
                Turbine%HD%y%WriteOutput, Turbine%SD%y%WriteOutput, Turbine%MAP%y%WriteOutput, Turbine%FEAM%y%WriteOutput, Turbine%MD%y%WriteOutput, &
                Turbine%IceF%y%WriteOutput, Turbine%IceD%y, Turbine%BD%y, Outputs)   
                        
END SUBROUTINE FillOutputAry_T                        
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FillOutputAry(p_FAST, y_FAST, IfWOutput, OpFMOutput, EDOutput, ADOutput, SrvDOutput, HDOutput, SDOutput, MAPOutput, FEAMOutput, &
                        MDOutput, IceFOutput, y_IceD, y_BD, OutputAry)

   TYPE(FAST_ParameterType), INTENT(IN)    :: p_FAST                             ! Glue-code simulation parameters
   TYPE(FAST_OutputFileType),INTENT(IN)    :: y_FAST                             ! Glue-code simulation outputs


   REAL(ReKi),               INTENT(IN)    :: IfWOutput (:)                      ! InflowWind WriteOutput values
   REAL(ReKi),               INTENT(IN)    :: OpFMOutput (:)                     ! OpenFOAM WriteOutput values
   REAL(ReKi),               INTENT(IN)    :: EDOutput (:)                       ! ElastoDyn WriteOutput values
   REAL(ReKi),               INTENT(IN)    :: ADOutput (:)                       ! AeroDyn WriteOutput values
   REAL(ReKi),               INTENT(IN)    :: SrvDOutput (:)                     ! ServoDyn WriteOutput values
   REAL(ReKi),               INTENT(IN)    :: HDOutput (:)                       ! HydroDyn WriteOutput values
   REAL(ReKi),               INTENT(IN)    :: SDOutput (:)                       ! SubDyn WriteOutput values
   REAL(ReKi),               INTENT(IN)    :: MAPOutput (:)                      ! MAP WriteOutput values
   REAL(ReKi),               INTENT(IN)    :: FEAMOutput (:)                     ! FEAMooring WriteOutput values
   REAL(ReKi),               INTENT(IN)    :: MDOutput (:)                       ! MoorDyn WriteOutput values
   REAL(ReKi),               INTENT(IN)    :: IceFOutput (:)                     ! IceFloe WriteOutput values
   TYPE(IceD_OutputType),    INTENT(IN)    :: y_IceD (:)                         ! IceDyn outputs (WriteOutput values are subset)
   TYPE(BD_OutputType),      INTENT(IN)    :: y_BD (:)                           ! BeamDyn outputs (WriteOutput values are subset)

   REAL(ReKi),               INTENT(OUT)   :: OutputAry(:)                       ! single array of output 
   
   INTEGER(IntKi)                          :: i                                  ! loop counter
   INTEGER(IntKi)                          :: indxLast                           ! The index of the last row value to be written to AllOutData for this time step (column).
   INTEGER(IntKi)                          :: indxNext                           ! The index of the next row value to be written to AllOutData for this time step (column).
   
   
            ! store individual module data into one array for output

      indxLast = 0
      indxNext = 1

      IF ( y_FAST%numOuts(Module_IfW) > 0 ) THEN
         indxLast = indxNext + SIZE(IfWOutput) - 1
         OutputAry(indxNext:indxLast) = IfWOutput
         indxNext = IndxLast + 1
      ELSEIF ( y_FAST%numOuts(Module_OpFM) > 0 ) THEN
         indxLast = indxNext + SIZE(OpFMOutput) - 1
         OutputAry(indxNext:indxLast) = OpFMOutput
         indxNext = IndxLast + 1
      END IF

      IF ( y_FAST%numOuts(Module_ED) > 0 ) THEN
         indxLast = indxNext + SIZE(EDOutput) - 1
         OutputAry(indxNext:indxLast) = EDOutput
         indxNext = IndxLast + 1
      END IF
         
      IF ( y_FAST%numOuts(Module_BD) > 0 ) THEN
         do i=1,SIZE(y_BD)
            indxLast = indxNext + SIZE(y_BD(i)%WriteOutput) - 1
            OutputAry(indxNext:indxLast) = y_BD(i)%WriteOutput
            indxNext = IndxLast + 1
         end do         
      END IF            
      
      IF ( y_FAST%numOuts(Module_AD) > 0 ) THEN
         indxLast = indxNext + SIZE(ADOutput) - 1
         OutputAry(indxNext:indxLast) = ADOutput
         indxNext = IndxLast + 1
      END IF
         
      IF ( y_FAST%numOuts(Module_SrvD) > 0 ) THEN
         indxLast = indxNext + SIZE(SrvDOutput) - 1
         OutputAry(indxNext:indxLast) = SrvDOutput
         indxNext = IndxLast + 1
      END IF

      IF ( y_FAST%numOuts(Module_HD) > 0 ) THEN
         indxLast = indxNext + SIZE(HDOutput) - 1
         OutputAry(indxNext:indxLast) = HDOutput
         indxNext = IndxLast + 1
      END IF

      IF ( y_FAST%numOuts(Module_SD) > 0 ) THEN
         indxLast = indxNext + SIZE(SDOutput) - 1
         OutputAry(indxNext:indxLast) = SDOutput
         indxNext = IndxLast + 1
      END IF
                  
      IF ( y_FAST%numOuts(Module_MAP) > 0 ) THEN
         indxLast = indxNext + SIZE(MAPOutput) - 1
         OutputAry(indxNext:indxLast) = MAPOutput
         indxNext = IndxLast + 1
      ELSEIF ( y_FAST%numOuts(Module_MD) > 0 ) THEN
         indxLast = indxNext + SIZE(MDOutput) - 1
         OutputAry(indxNext:indxLast) = MDOutput
         indxNext = IndxLast + 1
      ELSEIF ( y_FAST%numOuts(Module_FEAM) > 0 ) THEN
         indxLast = indxNext + SIZE(FEAMOutput) - 1
         OutputAry(indxNext:indxLast) = FEAMOutput
         indxNext = IndxLast + 1
      END IF
         
      IF ( y_FAST%numOuts(Module_IceF) > 0 ) THEN
         indxLast = indxNext + SIZE(IceFOutput) - 1
         OutputAry(indxNext:indxLast) = IceFOutput
         indxNext = IndxLast + 1
      ELSEIF ( y_FAST%numOuts(Module_IceD) > 0 ) THEN
         DO i=1,p_FAST%numIceLegs
            indxLast = indxNext + SIZE(y_IceD(i)%WriteOutput) - 1
            OutputAry(indxNext:indxLast) = y_IceD(i)%WriteOutput
            indxNext = IndxLast + 1
         END DO            
      END IF     
         
END SUBROUTINE FillOutputAry
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_InputSolve( p_FAST, BD, y_AD, u_AD, MeshMapData, ErrStat, ErrMsg )
! This routine sets the inputs required for BD--using the Option 2 solve method; currently the only inputs solved in this routine
! are the blade distributed loads from AD15; other inputs are solved in option 1.
!..................................................................................................................................

   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST                   ! Glue-code simulation parameters
   TYPE(BeamDyn_Data),             INTENT(INOUT)  :: BD                       ! BD Inputs at t
   TYPE(AD_OutputType),            INTENT(IN   )  :: y_AD                     ! AeroDyn outputs
   TYPE(AD_InputType),             INTENT(INOUT)  :: u_AD                     ! AD inputs (for AD-BD load transfer)
   
   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData              ! Data for mapping between modules
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                  ! Error status
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                   ! Error message
   
      ! local variables
   INTEGER(IntKi)                                 :: J                        ! Loops through nodes / elements
   INTEGER(IntKi)                                 :: K                        ! Loops through blades
   INTEGER(IntKi)                                 :: ErrStat2                 ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                           :: ErrMsg2                  ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*), PARAMETER                        :: RoutineName = 'BD_InputSolve' 

      ! Initialize error status
      
   ErrStat = ErrID_None
   ErrMsg = ""

            
      ! BD inputs on blade from AeroDyn
   IF (p_FAST%CompElast == Module_BD) THEN 
      
      IF ( p_FAST%CompAero == Module_AD ) THEN
         
         DO K = 1,p_FAST%nBeams ! Loop through all blades
            
            ! need to transfer the BD output blade motions to nodes on a sibling of the BD blade motion mesh:
            CALL Transfer_Line2_to_Line2( BD%y(k)%BldMotion, MeshMapData%y_BD_BldMotion_4Loads(k), MeshMapData%BD_L_2_BD_L(k), ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                        
            CALL Transfer_Line2_to_Line2( y_AD%BladeLoad(k), BD%Input(1,k)%DistrLoad, MeshMapData%AD_L_2_BDED_L_B(k), ErrStat2, ErrMsg2, u_AD%BladeMotion(k), MeshMapData%y_BD_BldMotion_4Loads(k) )
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               
         END DO
                  
      ELSE

         DO K = 1,p_FAST%nBeams ! Loop through all blades
            BD%Input(1,k)%DistrLoad%Force  = 0.0_ReKi
            BD%Input(1,k)%DistrLoad%Moment = 0.0_ReKi
         END DO         
         
      END IF
      
   END IF
      
                  
     
END SUBROUTINE BD_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ED_InputSolve( p_FAST, u_ED, y_ED, y_AD14, y_AD, y_SrvD, u_AD, u_SrvD, MeshMapData, ErrStat, ErrMsg )
! This routine sets the inputs required for ED--using the Option 2 solve method; currently the only input not solved in this routine
! are the fields on PlatformPtMesh and HubPtLoad,  which are solved in option 1.
!..................................................................................................................................

   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST                   ! Glue-code simulation parameters
   TYPE(ED_InputType),             INTENT(INOUT)  :: u_ED                     ! ED Inputs at t
   TYPE(ED_OutputType),            INTENT(IN   )  :: y_ED                     ! ElastoDyn outputs (need translation displacement on meshes for loads mapping)
   TYPE(AD14_OutputType),          INTENT(IN   )  :: y_AD14                   ! AeroDyn14 outputs
   TYPE(AD_OutputType),            INTENT(IN   )  :: y_AD                     ! AeroDyn outputs
   TYPE(AD_InputType),             INTENT(INOUT)  :: u_AD                     ! AD inputs (for AD-ED load transfer)
   TYPE(SrvD_OutputType),          INTENT(IN   )  :: y_SrvD                   ! ServoDyn outputs
   TYPE(SrvD_InputType),           INTENT(IN   )  :: u_SrvD                   ! ServoDyn inputs
   
   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData              ! Data for mapping between modules
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                  ! Error status
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                   ! Error message
   
      ! local variables
   INTEGER(IntKi)                                 :: J                        ! Loops through nodes / elements
   INTEGER(IntKi)                                 :: K                        ! Loops through blades
   INTEGER(IntKi)                                 :: ErrStat2                 ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                           :: ErrMsg2                  ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*), PARAMETER                        :: RoutineName = 'ED_InputSolve' 

      ! Initialize error status
      
   ErrStat = ErrID_None
   ErrMsg = ""

   
      ! ED inputs from ServoDyn
   IF ( p_FAST%CompServo == Module_SrvD ) THEN

      u_ED%GenTrq     = y_SrvD%GenTrq
      u_ED%HSSBrTrqC  = y_SrvD%HSSBrTrqC
      u_ED%BlPitchCom = y_SrvD%BlPitchCom
      u_ED%YawMom     = y_SrvD%YawMom
   !   u_ED%TBDrCon    = y_SrvD%TBDrCon !array
   
      IF (y_SrvD%NTMD%Mesh%Committed) THEN
      
         CALL Transfer_Point_to_Point( y_SrvD%NTMD%Mesh, u_ED%NacelleLoads, MeshMapData%SrvD_P_2_ED_P_N, ErrStat2, ErrMsg2, u_SrvD%NTMD%Mesh, y_ED%NacelleMotion )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//':u_ED%NacelleLoads' )      
            
      END IF
   
   ELSE !we'll just take the initial guesses..
   END IF

            
      ! ED inputs on blade from AeroDyn
   IF (p_FAST%CompElast == Module_ED) THEN 
      
      IF ( p_FAST%CompAero == Module_AD14 ) THEN   
      
         DO K = 1,SIZE(u_ED%BladeLn2Mesh,1) ! Loop through all blades (p_ED%NumBl)
            DO J = 1,y_AD14%OutputLoads(K)%Nnodes ! Loop through the blade nodes / elements (p_ED%BldNodes)

               u_ED%BladeLn2Mesh(K)%Force(:,J)  = y_AD14%OutputLoads(K)%Force(:,J)
               u_ED%BladeLn2Mesh(K)%Moment(:,J) = y_AD14%OutputLoads(K)%Moment(:,J)
            
            END DO !J
         END DO   !K 
      ELSEIF ( p_FAST%CompAero == Module_AD ) THEN
         
         DO K = 1,SIZE(u_ED%BladeLn2Mesh,1) ! Loop through all blades (p_ED%NumBl)
            CALL Transfer_Line2_to_Line2( y_AD%BladeLoad(k), u_ED%BladeLn2Mesh(k), MeshMapData%AD_L_2_BDED_L_B(k), ErrStat2, ErrMsg2, u_AD%BladeMotion(k), y_ED%BladeLn2Mesh(k) )
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         END DO
                  
      ELSE
         !p_FAST%CompAero = Module_None
         DO K = 1,SIZE(u_ED%BladeLn2Mesh,1) ! Loop through all blades (p_ED%NumBl)
            u_ED%BladeLn2Mesh(K)%Force  = 0.0_ReKi
            u_ED%BladeLn2Mesh(K)%Moment = 0.0_ReKi
         END DO         
         
      END IF
      
   END IF
      
                  
   IF ( p_FAST%CompAero == Module_AD14 ) THEN   
      u_ED%TowerLn2Mesh%Force  = 0.0_ReKi
      u_ED%TowerLn2Mesh%Moment = 0.0_ReKi
            
         ! add aero force to the tower, if it's provided:
      IF ( y_AD14%Twr_OutputLoads%Committed ) THEN
      
         ! we're mapping loads, so we also need the sibling meshes' displacements:
         
   !      CALL Transfer_Line2_to_Line2( )
      
         J = y_AD14%Twr_OutputLoads%NNodes
         
         IF ( y_AD14%Twr_OutputLoads%FIELDMASK(MASKID_FORCE) ) &
            u_ED%TowerLn2Mesh%Force(:,1:J)  = u_ED%TowerLn2Mesh%Force( :,1:J) + y_AD14%Twr_OutputLoads%Force
         
         IF ( y_AD14%Twr_OutputLoads%FIELDMASK(MASKID_MOMENT) ) &
            u_ED%TowerLn2Mesh%Moment(:,1:J) = u_ED%TowerLn2Mesh%Moment(:,1:J) + y_AD14%Twr_OutputLoads%Moment 
      
      END IF   
      
   ELSEIF ( p_FAST%CompAero == Module_AD ) THEN
      
      IF ( y_AD%TowerLoad%Committed ) THEN
         CALL Transfer_Line2_to_Line2( y_AD%TowerLoad, u_ED%TowerLn2Mesh, MeshMapData%AD_L_2_ED_L_T, ErrStat2, ErrMsg2, u_AD%TowerMotion, y_ED%TowerLn2Mesh )
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)         
      END IF
            
   ELSE
      u_ED%TowerLn2Mesh%Force  = 0.0_ReKi
      u_ED%TowerLn2Mesh%Moment = 0.0_ReKi      
   END IF
   
   u_ED%TwrAddedMass  = 0.0_ReKi
   u_ED%PtfmAddedMass = 0.0_ReKi
               
END SUBROUTINE ED_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SrvD_InputSolve( p_FAST, m_FAST, u_SrvD, y_ED, y_IfW, y_OpFM, MeshMapData, ErrStat, ErrMsg, y_SrvD_prev )
! This routine sets the inputs required for ServoDyn
!..................................................................................................................................

   TYPE(FAST_ParameterType),         INTENT(IN)     :: p_FAST       ! Glue-code simulation parameters
   TYPE(FAST_MiscVarType),           INTENT(IN)     :: m_FAST       ! Glue-code misc variables (including inputs from external sources like Simulink)
   TYPE(SrvD_InputType),             INTENT(INOUT)  :: u_SrvD       ! ServoDyn Inputs at t
   TYPE(ED_OutputType),              INTENT(IN)     :: y_ED         ! ElastoDyn outputs
   TYPE(InflowWind_OutputType),      INTENT(IN)     :: y_IfW        ! InflowWind outputs
   TYPE(OpFM_OutputType),            INTENT(IN)     :: y_OpFM       ! InflowWind outputs
   TYPE(SrvD_OutputType), OPTIONAL,  INTENT(IN)     :: y_SrvD_prev  ! ServoDyn outputs from t - dt
   TYPE(FAST_ModuleMapType),         INTENT(INOUT)  :: MeshMapData  ! Data for mapping between modules
   INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat      ! Error status
   CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg       ! Error message
!  TYPE(AD_OutputType),              INTENT(IN)     :: y_AD         ! AeroDyn outputs



   ErrStat = ErrID_None
   ErrMsg  = ""

      ! ServoDyn inputs from combination of InflowWind and ElastoDyn

   u_SrvD%YawAngle  = y_ED%YawAngle !nacelle yaw plus platform yaw

      ! Calculate horizontal hub-height wind direction and the nacelle yaw error estimate (both positive about zi-axis); these are
      !   zero if there is no wind input when InflowWind is not used:

      !bjj: rename pass YawAngle (not YawErr from ED)
   IF ( p_FAST%CompInflow == Module_IfW )  THEN 

      u_SrvD%WindDir  = ATAN2( y_IfW%VelocityUVW(2,1), y_IfW%VelocityUVW(1,1) )
      u_SrvD%YawErr   = u_SrvD%WindDir - y_ED%YawAngle
      u_SrvD%HorWindV = SQRT( y_IfW%VelocityUVW(1,1)**2 + y_IfW%VelocityUVW(2,1)**2 )

   ELSEIF ( p_FAST%CompInflow == Module_OpFM )  THEN 
      
      u_SrvD%WindDir  = ATAN2( y_OpFM%v(1), y_OpFM%u(1) )
      u_SrvD%YawErr   = u_SrvD%WindDir - y_ED%YawAngle
      u_SrvD%HorWindV = SQRT( y_OpFM%u(1)**2 + y_OpFM%v(1)**2 )
      
   ELSE  ! No wind inflow

      u_SrvD%WindDir  = 0.0
      u_SrvD%YawErr   = 0.0
      u_SrvD%HorWindV = 0.0

   ENDIF

      ! ServoDyn inputs from ServoDyn outputs at previous step
      ! Jason says this violates the framework, but it's only for the Bladed DLL, which itself violates the framework, so I don't care.
   IF (PRESENT(y_SrvD_prev)) THEN
      u_SrvD%ElecPwr_prev = y_SrvD_prev%ElecPwr  ! we want to know the electrical power from the previous time step  (for the Bladed DLL)
      u_SrvD%GenTrq_prev  = y_SrvD_prev%GenTrq   ! we want to know the electrical generator torque from the previous time step  (for the Bladed DLL)
   ! Otherwise, we'll use the guess provided by the module (this only happens at Step=0)
   END IF

      ! ServoDyn inputs from ElastoDyn
   u_SrvD%Yaw       = y_ED%Yaw  !nacelle yaw
   u_SrvD%YawRate   = y_ED%YawRate
   u_SrvD%BlPitch   = y_ED%BlPitch
   u_SrvD%LSS_Spd   = y_ED%LSS_Spd
   u_SrvD%HSS_Spd   = y_ED%HSS_Spd
   u_SrvD%RotSpeed  = y_ED%RotSpeed
   
   IF ( p_FAST%CompElast == Module_BD )  THEN    
   !bjj: FIX ME! TODO: need to get this from BeamDyn!!!!
      u_SrvD%RootMxc(:) = 0.0_ReKi
      u_SrvD%RootMyc(:) = 0.0_ReKi
   ELSE
      u_SrvD%RootMxc = y_ED%RootMxc
      u_SrvD%RootMyc = y_ED%RootMyc
   END IF

   
   u_SrvD%YawBrTAxp = y_ED%YawBrTAxp
   u_SrvD%YawBrTAyp = y_ED%YawBrTAyp
   u_SrvD%LSSTipPxa = y_ED%LSSTipPxa

   u_SrvD%LSSTipMxa = y_ED%LSSTipMxa
   u_SrvD%LSSTipMya = y_ED%LSSTipMya
   u_SrvD%LSSTipMza = y_ED%LSSTipMza
   u_SrvD%LSSTipMys = y_ED%LSSTipMys
   u_SrvD%LSSTipMzs = y_ED%LSSTipMzs
   
   u_SrvD%YawBrMyn  = y_ED%YawBrMyn
   u_SrvD%YawBrMzn  = y_ED%YawBrMzn
   u_SrvD%NcIMURAxs = y_ED%NcIMURAxs
   u_SrvD%NcIMURAys = y_ED%NcIMURAys
   u_SrvD%NcIMURAzs = y_ED%NcIMURAzs

   u_SrvD%RotPwr    = y_ED%RotPwr

   !   ! ServoDyn inputs from AeroDyn
   !IF ( p_FAST%CompAero == Module_AD ) THEN
   !ELSE
   !END IF
   !
   
   IF (u_SrvD%NTMD%Mesh%Committed) THEN
      
         !bjj: watch error handling if SrvD_InputSolve ever gets more complicated
      CALL Transfer_Point_to_Point( y_ED%NacelleMotion, u_SrvD%NTMD%Mesh, MeshMapData%ED_P_2_SrvD_P_N, ErrStat, ErrMsg )
            
   END IF
   
#ifdef SIMULINK_TIMESHIFT   
      ! we're going to use the extrapolated values instead of the old values (Simulink inputs are from t, not t+dt)
   CALL SrvD_SetExternalInputs( p_FAST, m_FAST, u_SrvD )
#endif
      
                        
END SUBROUTINE SrvD_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SrvD_SetExternalInputs( p_FAST, m_FAST, u_SrvD )
! This routine sets the inputs required for ServoDyn from an external source (Simulink)
!..................................................................................................................................

   TYPE(FAST_ParameterType),         INTENT(IN)     :: p_FAST       ! Glue-code simulation parameters
   TYPE(FAST_MiscVarType),           INTENT(IN)     :: m_FAST       ! Glue-code misc variables (including inputs from external sources like Simulink)
   TYPE(SrvD_InputType),             INTENT(INOUT)  :: u_SrvD       ! ServoDyn Inputs at t

      ! local variables
   INTEGER(IntKi)                                   :: i            ! loop counter
   
      ! we are going to use extrapolated values because these external values from Simulink are at n instead of n+1
   u_SrvD%ExternalGenTrq       =  m_FAST%ExternInput%GenTrq     
   u_SrvD%ExternalElecPwr      =  m_FAST%ExternInput%ElecPwr    
   u_SrvD%ExternalYawPosCom    =  m_FAST%ExternInput%YawPosCom  
   u_SrvD%ExternalYawRateCom   =  m_FAST%ExternInput%YawRateCom 
   u_SrvD%ExternalHSSBrFrac    =  m_FAST%ExternInput%HSSBrFrac 

   do i=1,SIZE(u_SrvD%ExternalBlPitchCom)
      u_SrvD%ExternalBlPitchCom(i)   = m_FAST%ExternInput%BlPitchCom(i)
   end do

END SUBROUTINE SrvD_SetExternalInputs
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Transfer_SD_to_HD( y_SD, u_HD_M_LumpedMesh, u_HD_M_DistribMesh, MeshMapData, ErrStat, ErrMsg )
! This routine transfers the SD outputs into inputs required for HD
!..................................................................................................................................
   TYPE(SD_OutputType),         INTENT(IN   ) :: y_SD                         ! The outputs of the structural dynamics module
   TYPE(MeshType),              INTENT(INOUT) :: u_HD_M_LumpedMesh            ! HydroDyn input mesh (separated here so that we can use temp meshes in ED_SD_HD_InputSolve)
   TYPE(MeshType),              INTENT(INOUT) :: u_HD_M_DistribMesh           ! HydroDyn input mesh (separated here so that we can use temp meshes in ED_SD_HD_InputSolve)
   TYPE(FAST_ModuleMapType),    INTENT(INOUT) :: MeshMapData

   INTEGER(IntKi),              INTENT(  OUT) :: ErrStat                      ! Error status of the operation
   CHARACTER(*),                INTENT(  OUT) :: ErrMsg                       ! Error message if ErrStat /= ErrID_None
   
      ! local variables
   INTEGER(IntKi)                             :: ErrStat2                     ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                       :: ErrMsg2                      ! temporary Error message if ErrStat /= ErrID_None
      
      
   ErrStat = ErrID_None
   ErrMsg = ""
   
       
   IF ( u_HD_M_LumpedMesh%Committed ) THEN 

      ! These are the motions for the lumped point loads associated viscous drag on the WAMIT body and/or filled/flooded lumped forces of the WAMIT body
      CALL Transfer_Point_to_Point( y_SD%y2Mesh, u_HD_M_LumpedMesh, MeshMapData%SD_P_2_HD_M_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,'Transfer_SD_to_HD (u_HD%Morison%LumpedMesh)' )      
         
   END IF
   
   IF ( u_HD_M_DistribMesh%Committed ) THEN 
         
      ! These are the motions for the HD line2 (distributed) loads associated viscous drag on the WAMIT body and/or filled/flooded distributed forces of the WAMIT body
      CALL Transfer_Point_to_Line2( y_SD%y2Mesh, u_HD_M_DistribMesh, MeshMapData%SD_P_2_HD_M_L, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,'Transfer_SD_to_HD (u_HD%Morison%DistribMesh)' )      
   END IF
   
END SUBROUTINE Transfer_SD_to_HD
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Transfer_ED_to_HD( y_ED, u_HD, MeshMapData, ErrStat, ErrMsg )
! This routine transfers the ED outputs into inputs required for HD
!..................................................................................................................................
   TYPE(ED_OutputType),         INTENT(IN   ) :: y_ED                         ! The outputs of the structural dynamics module
   TYPE(HydroDyn_InputType),    INTENT(INOUT) :: u_HD                         ! HydroDyn input
   TYPE(FAST_ModuleMapType),    INTENT(INOUT) :: MeshMapData

   INTEGER(IntKi),              INTENT(OUT)   :: ErrStat                      ! Error status of the operation
   CHARACTER(*),                INTENT(OUT)   :: ErrMsg                       ! Error message if ErrStat /= ErrID_None
   
      ! local variables
   INTEGER(IntKi)                             :: ErrStat2                     ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                       :: ErrMsg2                      ! temporary Error message if ErrStat /= ErrID_None
      
      
   ErrStat = ErrID_None
   ErrMsg = ""
   
   !bjj: We do this without all the extra meshcopy/destroy calls with u_mapped because these inputs are only from one mesh
   
   IF ( u_HD%Mesh%Committed ) THEN

      ! These are the motions for the lumped point loads associated the WAMIT body and include: hydrostatics, radiation memory effect,
      !    wave kinematics, additional preload, additional stiffness, additional linear damping, additional quadratic damping,
      !    hydrodynamic added mass

      CALL Transfer_Point_to_Point( y_ED%PlatformPtMesh, u_HD%Mesh, MeshMapData%ED_P_2_HD_W_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,'Transfer_ED_to_HD (u_HD%Mesh)' )

   END IF !WAMIT
   
   
   IF ( u_HD%Morison%LumpedMesh%Committed ) THEN 

      ! These are the motions for the lumped point loads associated viscous drag on the WAMIT body and/or filled/flooded lumped forces of the WAMIT body
      CALL Transfer_Point_to_Point( y_ED%PlatformPtMesh, u_HD%Morison%LumpedMesh, MeshMapData%ED_P_2_HD_M_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,'Transfer_ED_to_HD (u_HD%Morison%LumpedMesh)' )
         
   END IF
   
   IF ( u_HD%Morison%DistribMesh%Committed ) THEN 
         
      ! These are the motions for the line2 (distributed) loads associated viscous drag on the WAMIT body and/or filled/flooded distributed forces of the WAMIT body
      CALL Transfer_Point_to_Line2( y_ED%PlatformPtMesh, u_HD%Morison%DistribMesh, MeshMapData%ED_P_2_HD_M_L, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,'Transfer_ED_to_HD (u_HD%Morison%DistribMesh)' )

   END IF
   
END SUBROUTINE Transfer_ED_to_HD
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Transfer_ED_to_HD_SD_BD_Mooring( p_FAST, y_ED, u_HD, u_SD, u_MAP, u_FEAM, u_MD, u_BD, MeshMapData, ErrStat, ErrMsg )
! This routine transfers the ED outputs into inputs required for HD, SD, BD, MAP, and/or FEAM
!..................................................................................................................................
   TYPE(FAST_ParameterType),    INTENT(IN)    :: p_FAST                       ! Glue-code simulation parameters
   TYPE(ED_OutputType),         INTENT(IN   ) :: y_ED                         ! The outputs of the structural dynamics module
   TYPE(HydroDyn_InputType),    INTENT(INOUT) :: u_HD                         ! HydroDyn input
   TYPE(SD_InputType),          INTENT(INOUT) :: u_SD                         ! SubDyn input
   TYPE(MAP_InputType),         INTENT(INOUT) :: u_MAP                        ! MAP input
   TYPE(FEAM_InputType),        INTENT(INOUT) :: u_FEAM                       ! FEAM input
   TYPE(MD_InputType),          INTENT(INOUT) :: u_MD                         ! MoorDyn input
   TYPE(BD_InputType),          INTENT(INOUT) :: u_BD(:)                      ! BeamDyn inputs
   TYPE(FAST_ModuleMapType),    INTENT(INOUT) :: MeshMapData

   INTEGER(IntKi),              INTENT(OUT)   :: ErrStat                      ! Error status of the operation
   CHARACTER(*),                INTENT(OUT)   :: ErrMsg                       ! Error message if ErrStat /= ErrID_None
   
      ! local variables
   INTEGER(IntKi)                             :: ErrStat2                     ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                       :: ErrMsg2                      ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*), PARAMETER                    :: RoutineName = 'Transfer_ED_to_HD_SD_BD_Mooring'   
      
   ErrStat = ErrID_None
   ErrMsg = ""
     
      ! transfer ED outputs to other modules used in option 1:
            
   IF ( p_FAST%CompSub == Module_SD  ) THEN
      
         ! Map ED (motion) outputs to SD inputs:                     
      CALL Transfer_Point_to_Point( y_ED%PlatformPtMesh, u_SD%TPMesh, MeshMapData%ED_P_2_SD_TP, ErrStat2, ErrMsg2 ) 
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//':u_SD%TPMesh' )
                                      
                  
   ELSEIF ( p_FAST%CompHydro == Module_HD ) THEN
         ! Map ED outputs to HD inputs:
      CALL Transfer_ED_to_HD( y_ED, u_HD, MeshMapData, ErrStat2, ErrMsg2 )                        
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName )
            
   END IF      
   
   IF ( p_FAST%CompElast == Module_BD ) THEN
      ! map ED root motion outputs to BeamDyn:
      CALL Transfer_ED_to_BD(y_ED, u_BD, MeshMapData, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName )
   END IF
      
   
   IF ( p_FAST%CompMooring == Module_MAP ) THEN
      
         ! motions:
      CALL Transfer_Point_to_Point( y_ED%PlatformPtMesh, u_MAP%PtFairDisplacement, MeshMapData%ED_P_2_Mooring_P, ErrStat, ErrMsg )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//'u_MAP%PtFairDisplacement' )
                                 
   ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
         ! motions:
      CALL Transfer_Point_to_Point( y_ED%PlatformPtMesh, u_MD%PtFairleadDisplacement, MeshMapData%ED_P_2_Mooring_P, ErrStat, ErrMsg )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//'u_MD%PtFairleadDisplacement' )
                        
   ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
         ! motions:
      CALL Transfer_Point_to_Point( y_ED%PlatformPtMesh, u_FEAM%PtFairleadDisplacement, MeshMapData%ED_P_2_Mooring_P, ErrStat, ErrMsg )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//'u_FEAM%PtFairleadDisplacement' )
                        
   END IF
   
            
END SUBROUTINE Transfer_ED_to_HD_SD_BD_Mooring
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE MAP_InputSolve(  u_MAP, y_ED, MeshMapData, ErrStat, ErrMsg )
! This routine sets the inputs required for MAP.
!..................................................................................................................................

      ! Passed variables
   TYPE(MAP_InputType),         INTENT(INOUT) :: u_MAP                        ! MAP input
   TYPE(ED_OutputType),         INTENT(IN   ) :: y_ED                         ! The outputs of the structural dynamics module
   TYPE(FAST_ModuleMapType),    INTENT(INOUT) :: MeshMapData

   INTEGER(IntKi),              INTENT(  OUT) :: ErrStat                      ! Error status of the operation
   CHARACTER(*)  ,              INTENT(  OUT) :: ErrMsg                       ! Error message if ErrStat /= ErrID_None


      !----------------------------------------------------------------------------------------------------
      ! Map ED outputs to MAP inputs
      !----------------------------------------------------------------------------------------------------
      ! motions:
   CALL Transfer_Point_to_Point( y_ED%PlatformPtMesh, u_MAP%PtFairDisplacement, MeshMapData%ED_P_2_Mooring_P, ErrStat, ErrMsg )


END SUBROUTINE MAP_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FEAM_InputSolve(  u_FEAM, y_ED, MeshMapData, ErrStat, ErrMsg )
! This routine sets the inputs required for FEAM.
!..................................................................................................................................

      ! Passed variables
   TYPE(FEAM_InputType),        INTENT(INOUT) :: u_FEAM                       ! FEAM input
   TYPE(ED_OutputType),         INTENT(IN   ) :: y_ED                         ! The outputs of the structural dynamics module
   TYPE(FAST_ModuleMapType),    INTENT(INOUT) :: MeshMapData

   INTEGER(IntKi),              INTENT(  OUT) :: ErrStat                      ! Error status of the operation
   CHARACTER(*)  ,              INTENT(  OUT) :: ErrMsg                       ! Error message if ErrStat /= ErrID_None


      !----------------------------------------------------------------------------------------------------
      ! Map ED outputs to FEAM inputs
      !----------------------------------------------------------------------------------------------------
      ! motions:
   CALL Transfer_Point_to_Point( y_ED%PlatformPtMesh, u_FEAM%PtFairleadDisplacement, MeshMapData%ED_P_2_Mooring_P, ErrStat, ErrMsg )


END SUBROUTINE FEAM_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE MD_InputSolve(  u_MD, y_ED, MeshMapData, ErrStat, ErrMsg )
! This routine sets the inputs required for MoorDyn.
!..................................................................................................................................

      ! Passed variables
   TYPE(MD_InputType),          INTENT(INOUT) :: u_MD                         ! MoorDyn input
   TYPE(ED_OutputType),         INTENT(IN   ) :: y_ED                         ! The outputs of the structural dynamics module
   TYPE(FAST_ModuleMapType),    INTENT(INOUT) :: MeshMapData

   INTEGER(IntKi),              INTENT(  OUT) :: ErrStat                      ! Error status of the operation
   CHARACTER(*)  ,              INTENT(  OUT) :: ErrMsg                       ! Error message if ErrStat /= ErrID_None


      !----------------------------------------------------------------------------------------------------
      ! Map ED outputs to MoorDyn inputs
      !----------------------------------------------------------------------------------------------------
      ! motions:
   CALL Transfer_Point_to_Point( y_ED%PlatformPtMesh, u_MD%PtFairleadDisplacement, MeshMapData%ED_P_2_Mooring_P, ErrStat, ErrMsg )


END SUBROUTINE MD_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE IceFloe_InputSolve(  u_IceF, y_SD, MeshMapData, ErrStat, ErrMsg )
! This routine sets the inputs required for IceFloe.
!..................................................................................................................................

      ! Passed variables
   TYPE(IceFloe_InputType),     INTENT(INOUT) :: u_IceF                       ! IceFloe input
   TYPE(SD_OutputType),         INTENT(IN   ) :: y_SD                         ! SubDyn outputs
   TYPE(FAST_ModuleMapType),    INTENT(INOUT) :: MeshMapData

   INTEGER(IntKi),              INTENT(  OUT) :: ErrStat                      ! Error status of the operation
   CHARACTER(*)  ,              INTENT(  OUT) :: ErrMsg                       ! Error message if ErrStat /= ErrID_None


      !----------------------------------------------------------------------------------------------------
      ! Map SD outputs to IceFloe inputs
      !----------------------------------------------------------------------------------------------------
      ! motions:
   CALL Transfer_Point_to_Point( y_SD%y2Mesh, u_IceF%IceMesh, MeshMapData%SD_P_2_IceF_P, ErrStat, ErrMsg )

END SUBROUTINE IceFloe_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE IceD_InputSolve(  u_IceD, y_SD, MeshMapData, legNum, ErrStat, ErrMsg )
! This routine sets the inputs required for IceFloe.
!..................................................................................................................................

      ! Passed variables
   TYPE(IceD_InputType),        INTENT(INOUT) :: u_IceD                       ! IceDyn input
   TYPE(SD_OutputType),         INTENT(IN   ) :: y_SD                         ! SubDyn outputs
   TYPE(FAST_ModuleMapType),    INTENT(INOUT) :: MeshMapData
   INTEGER(IntKi),              INTENT(IN   ) :: legNum

   INTEGER(IntKi),              INTENT(  OUT) :: ErrStat                      ! Error status of the operation
   CHARACTER(*)  ,              INTENT(  OUT) :: ErrMsg                       ! Error message if ErrStat /= ErrID_None


      !----------------------------------------------------------------------------------------------------
      ! Map SD outputs to IceFloe inputs
      !----------------------------------------------------------------------------------------------------
      ! motions:
   CALL Transfer_Point_to_Point( y_SD%y2Mesh, u_IceD%PointMesh, MeshMapData%SD_P_2_IceD_P(legNum), ErrStat, ErrMsg )

END SUBROUTINE IceD_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Transfer_ED_to_BD(  y_ED, u_BD, MeshMapData, ErrStat, ErrMsg )
! This routine sets the inputs required for IceFloe.
!..................................................................................................................................

      ! Passed variables
   TYPE(BD_InputType),          INTENT(INOUT) :: u_BD(:)                      ! BeamDyn inputs
   TYPE(ED_OutputType),         INTENT(IN   ) :: y_ED                         ! ElastoDyn outputs
   TYPE(FAST_ModuleMapType),    INTENT(INOUT) :: MeshMapData
   
   INTEGER(IntKi),              INTENT(  OUT) :: ErrStat                      ! Error status of the operation
   CHARACTER(*)  ,              INTENT(  OUT) :: ErrMsg                       ! Error message if ErrStat /= ErrID_None

      

   
      ! local variables
   INTEGER(IntKi)                             :: k
   INTEGER(IntKi)                             :: ErrStat2                  ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                       :: ErrMsg2                   ! temporary Error message if ErrStat /= ErrID_None

   CHARACTER(*), PARAMETER                    :: RoutineName = 'Transfer_ED_to_BD'   
   
   ErrStat = ErrID_None
   ErrMsg = ""   
   
      !----------------------------------------------------------------------------------------------------
      ! Map ED outputs to BeamDyn inputs
      !----------------------------------------------------------------------------------------------------
      ! motions:
   do k = 1,size(y_ED%BladeRootMotion)
      CALL Transfer_Point_to_Point( y_ED%BladeRootMotion(k), u_BD(k)%RootMotion, MeshMapData%ED_P_2_BD_P(k), ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
#ifdef DEBUG_BEAMDYN2
   call meshprintinfo(80+k,y_ED%BladeRootMotion(k))
   call meshprintinfo(60+k,u_BD(k)%RootMotion)   
#endif
      
      
   end do
   

END SUBROUTINE Transfer_ED_to_BD
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Transfer_ED_to_BD_tmp(  y_ED, MeshMapData, ErrStat, ErrMsg )
! This routine sets the inputs required for IceFloe.
!..................................................................................................................................

      ! Passed variables
   TYPE(ED_OutputType),         INTENT(IN   ) :: y_ED                         ! ElastoDyn outputs
   TYPE(FAST_ModuleMapType),    INTENT(INOUT) :: MeshMapData
   
   INTEGER(IntKi),              INTENT(  OUT) :: ErrStat                      ! Error status of the operation
   CHARACTER(*)  ,              INTENT(  OUT) :: ErrMsg                       ! Error message if ErrStat /= ErrID_None

      

   
      ! local variables
   INTEGER(IntKi)                             :: k
   INTEGER(IntKi)                             :: ErrStat2                  ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                       :: ErrMsg2                   ! temporary Error message if ErrStat /= ErrID_None

   CHARACTER(*), PARAMETER                    :: RoutineName = 'Transfer_ED_to_BD_tmp'   
   
   ErrStat = ErrID_None
   ErrMsg = ""   
   
      !----------------------------------------------------------------------------------------------------
      ! Map ED outputs to BeamDyn inputs
      !----------------------------------------------------------------------------------------------------
      ! motions:
   do k = 1,size(y_ED%BladeRootMotion)
      CALL Transfer_Point_to_Point( y_ED%BladeRootMotion(k), MeshMapData%u_BD_RootMotion(k), MeshMapData%ED_P_2_BD_P(k), ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   end do
   

END SUBROUTINE Transfer_ED_to_BD_tmp
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Transfer_HD_to_SD( u_mapped, u_SD_LMesh, u_mapped_positions, y_HD, u_HD_M_LumpedMesh, u_HD_M_DistribMesh, MeshMapData, ErrStat, ErrMsg )
! This routine transfers the HD outputs into inputs required for ED. Note that this *adds* to the values already in 
! u_SD_LMesh (so initialize it before calling this routine).
!..................................................................................................................................
   TYPE(MeshType),                 INTENT(INOUT)  :: u_mapped                 ! temporary copy of SD mesh (an argument to avoid another temporary mesh copy)
   TYPE(MeshType),                 INTENT(INOUT)  :: u_SD_LMesh               ! SD Inputs on LMesh at t (separate so we can call from ED_SD_HD_InputOutput solve with temp meshes)
   TYPE(MeshType),                 INTENT(IN   )  :: u_mapped_positions       ! Mesh sibling of u_mapped, with displaced positions
   TYPE(HydroDyn_OutputType),      INTENT(IN   )  :: y_HD                     ! HydroDyn outputs
   TYPE(MeshType),                 INTENT(IN   )  :: u_HD_M_LumpedMesh        ! HydroDyn input mesh (separate so we can call from ED_SD_HD_InputOutput solve with temp meshes)
   TYPE(MeshType),                 INTENT(IN   )  :: u_HD_M_DistribMesh       ! HydroDyn input mesh (separate so we can call from ED_SD_HD_InputOutput solve with temp meshes)
   
   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData              ! Data for mapping between modules
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                  ! Error status
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                   ! Error message
   
      ! local variables
   INTEGER(IntKi)                                 :: ErrStat2                  ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                           :: ErrMsg2                   ! temporary Error message if ErrStat /= ErrID_None

   CHARACTER(*), PARAMETER                        :: RoutineName = 'Transfer_HD_to_SD'   
   
   ErrStat = ErrID_None
   ErrMsg = ""
         
   !assumes u_SD%LMesh%Committed   (i.e., u_SD_LMesh%Committed)
         
      IF ( y_HD%Morison%LumpedMesh%Committed ) THEN      
         ! we're mapping loads, so we also need the sibling meshes' displacements:
         CALL Transfer_Point_to_Point( y_HD%Morison%LumpedMesh, u_mapped, MeshMapData%HD_M_P_2_SD_P, ErrStat2, ErrMsg2, u_HD_M_LumpedMesh, u_mapped_positions )   
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
         
         u_SD_LMesh%Force  = u_SD_LMesh%Force  + u_mapped%Force
         u_SD_LMesh%Moment = u_SD_LMesh%Moment + u_mapped%Moment     
         
#ifdef DEBUG_MESH_TRANSFER              
         CALL WrScr('********************************************************')
         CALL WrScr('****   SD to HD point-to-point (morison lumped)    *****')
         CALL WrScr('********************************************************')
         CALL WriteMappingTransferToFile(u_mapped, u_mapped_positions, u_HD_M_LumpedMesh, y_HD%Morison%LumpedMesh,&
               MeshMapData%SD_P_2_HD_M_P, MeshMapData%HD_M_P_2_SD_P, &
               'SD_y2_HD_ML_Meshes_t'//TRIM(Num2LStr(0))//'.bin' )
         !print *
         !pause
         
#endif         
                  
      END IF
      
      IF ( y_HD%Morison%DistribMesh%Committed ) THEN      
         ! we're mapping loads, so we also need the sibling meshes' displacements:
         CALL Transfer_Line2_to_Point( y_HD%Morison%DistribMesh, u_mapped, MeshMapData%HD_M_L_2_SD_P, ErrStat2, ErrMsg2, u_HD_M_DistribMesh, u_mapped_positions )   
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
         
         u_SD_LMesh%Force  = u_SD_LMesh%Force  + u_mapped%Force
         u_SD_LMesh%Moment = u_SD_LMesh%Moment + u_mapped%Moment

#ifdef DEBUG_MESH_TRANSFER        
         CALL WrScr('********************************************************')
         CALL WrScr('**** SD to HD point-to-line2 (morison distributed) *****')
         CALL WrScr('********************************************************')
         CALL WriteMappingTransferToFile(u_mapped, u_mapped_positions, u_HD_M_DistribMesh,y_HD%Morison%DistribMesh,&
               MeshMapData%SD_P_2_HD_M_L, MeshMapData%HD_M_L_2_SD_P, &
               'SD_y2_HD_MD_Meshes_t'//TRIM(Num2LStr(0))//'.bin' )         
         !print *
        ! pause
#endif         
                  
      END IF

END SUBROUTINE Transfer_HD_to_SD
!----------------------------------------------------------------------------------------------------------------------------------
REAL(ReKi) FUNCTION GetPerturb(x)
   REAL(ReKi), INTENT(IN) :: x
      
   !GetPerturb = sqrt( EPSILON(x)) * max( abs(x), 1._ReKi)  
!      GetPerturb = 1.0e6
   GetPerturb = 1.0
      
END FUNCTION GetPerturb
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ED_HD_InputOutputSolve(  this_time, p_FAST, calcJacobian &
                                  , u_ED, p_ED, x_ED, xd_ED, z_ED, OtherSt_ED, y_ED &
                                  , u_HD, p_HD, x_HD, xd_HD, z_HD, OtherSt_HD, y_HD & 
                                  , u_MAP, y_MAP, u_FEAM, y_FEAM, u_MD, y_MD & 
                                  , MeshMapData , ErrStat, ErrMsg )
! This routine performs the Input-Output solve for ED and HD.
! Note that this has been customized for the physics in the problems and is not a general solution.
!..................................................................................................................................

   USE ElastoDyn
   USE HydroDyn

      ! Passed variables

   REAL(DbKi)                        , INTENT(IN   ) :: this_time                 ! The current simulation time (actual or time of prediction)
   TYPE(FAST_ParameterType)          , INTENT(IN   ) :: p_FAST                    ! Glue-code simulation parameters
   LOGICAL                           , INTENT(IN   ) :: calcJacobian              ! Should we calculate Jacobians this time? (should be TRUE on initialization, then can be false [significantly reducing computational time])
   
      !ElastoDyn:                    
   TYPE(ED_ContinuousStateType)      , INTENT(IN   ) :: x_ED                      ! Continuous states
   TYPE(ED_DiscreteStateType)        , INTENT(IN   ) :: xd_ED                     ! Discrete states
   TYPE(ED_ConstraintStateType)      , INTENT(IN   ) :: z_ED                      ! Constraint states
   TYPE(ED_OtherStateType)           , INTENT(INOUT) :: OtherSt_ED                ! Other/optimization states
   TYPE(ED_ParameterType)            , INTENT(IN   ) :: p_ED                      ! Parameters
   TYPE(ED_InputType)                , INTENT(INOUT) :: u_ED                      ! System inputs
   TYPE(ED_OutputType)               , INTENT(INOUT) :: y_ED                      ! System outputs
   
      !HydroDyn: 
   TYPE(HydroDyn_ContinuousStateType), INTENT(IN   ) :: x_HD                      ! Continuous states
   TYPE(HydroDyn_DiscreteStateType)  , INTENT(IN   ) :: xd_HD                     ! Discrete states
   TYPE(HydroDyn_ConstraintStateType), INTENT(IN   ) :: z_HD                      ! Constraint states
   TYPE(HydroDyn_OtherStateType)     , INTENT(INOUT) :: OtherSt_HD                ! Other/optimization states
   TYPE(HydroDyn_ParameterType)      , INTENT(IN   ) :: p_HD                      ! Parameters
   TYPE(HydroDyn_InputType)          , INTENT(INOUT) :: u_HD                      ! System inputs
   TYPE(HydroDyn_OutputType)         , INTENT(INOUT) :: y_HD                      ! System outputs

      ! MAP/FEAM/MoorDyn:
   TYPE(MAP_OutputType),              INTENT(IN   )  :: y_MAP                     ! MAP outputs
   TYPE(MAP_InputType),               INTENT(INOUT)  :: u_MAP                     ! MAP inputs (INOUT just because I don't want to use another tempoarary mesh and we'll overwrite this later)
   TYPE(FEAM_OutputType),             INTENT(IN   )  :: y_FEAM                    ! FEAM outputs
   TYPE(FEAM_InputType),              INTENT(INOUT)  :: u_FEAM                    ! FEAM inputs (INOUT just because I don't want to use another tempoarary mesh and we'll overwrite this later)
   TYPE(MD_OutputType),               INTENT(IN   )  :: y_MD                      ! MoorDyn outputs
   TYPE(MD_InputType),                INTENT(INOUT)  :: u_MD                      ! MoorDyn inputs (INOUT just because I don't want to use another tempoarary mesh and we'll overwrite this later)
      
   TYPE(FAST_ModuleMapType)          , INTENT(INOUT) :: MeshMapData
   INTEGER(IntKi)                    , INTENT(  OUT) :: ErrStat                   ! Error status of the operation
   CHARACTER(*)                      , INTENT(  OUT) :: ErrMsg                    ! Error message if ErrStat /= ErrID_None

   ! Local variables:
   INTEGER,                                PARAMETER :: NumInputs = SizeJac_ED_HD !12
   REAL(ReKi),                             PARAMETER :: TOL_Squared = (1.0E-4)**2 !not currently used because KMax = 1
   REAL(ReKi)                                        :: ThisPerturb               ! an arbitrary perturbation (these are linear, so it shouldn't matter)
   
   REAL(ReKi)                                        :: u(           NumInputs)   ! 6 loads, 6 accelerations
   REAL(ReKi)                                        :: u_perturb(   NumInputs)   ! 6 loads, 6 accelerations
   REAL(ReKi)                                        :: u_delta(     NumInputs)   !
   REAL(ReKi)                                        :: Fn_U_perturb(NumInputs)   ! value of U with perturbations
   REAL(ReKi)                                        :: Fn_U_Resid(  NumInputs)   ! Residual of U
   
                                                                                  
   TYPE(ED_OutputType)                               :: y_ED_input                ! Copy of system outputs sent to this routine (routine input value)
   TYPE(ED_InputType)                                :: u_ED_perturb              ! Perturbed system inputs
   TYPE(ED_OutputType)                               :: y_ED_perturb              ! Perturbed system outputs
   TYPE(HydroDyn_InputType)                          :: u_HD_perturb              ! Perturbed system inputs
   TYPE(HydroDyn_OutputType)                         :: y_HD_perturb              ! Perturbed system outputs
                                                                                  
                                                                                  
   INTEGER(IntKi)                                    :: i                         ! loop counter (jacobian column number)
   INTEGER(IntKi)                                    :: K                         ! Input-output-solve iteration counter
   INTEGER(IntKi)                                    :: ErrStat2                  ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                              :: ErrMsg2                   ! temporary Error message if ErrStat /= ErrID_None
   
   CHARACTER(*), PARAMETER                           :: RoutineName = 'ED_HD_InputOutputSolve'
   
#ifdef OUTPUT_ADDEDMASS   
   REAL(ReKi)                                        :: AddedMassMatrix(6,6)
   INTEGER                                           :: UnAM
#endif
#ifdef OUTPUT_JACOBIAN
   INTEGER                                           :: UnJac
#endif

   ! Note: p_FAST%UJacSclFact is a scaling factor that gets us similar magnitudes between loads and accelerations...
 
!bjj: note, that this routine may have a problem if there is remapping done
    
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! note this routine should be called only
   ! IF ( p_FAST%CompHydro == Module_HD .AND. p_FAST%CompSub /= Module_SD ) 
                           
      !----------------------------------------------------------------------------------------------------
      ! Some more record keeping stuff:
      !---------------------------------------------------------------------------------------------------- 
         
         ! We need to know the outputs that were sent to this routine:
      CALL ED_CopyOutput( y_ED, y_ED_input, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
         ! Local copies for perturbing inputs and outputs (computing Jacobian):
      IF ( calcJacobian ) THEN         
         CALL ED_CopyInput(  u_ED, u_ED_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                 
         CALL ED_CopyOutput( y_ED, y_ED_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )         
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL HydroDyn_CopyInput(  u_HD, u_HD_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )           
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL HydroDyn_CopyOutput( y_HD, y_HD_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )  
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END IF
         
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF
      
      !----------------------------------------------------------------------------------------------------
      ! set up u vector, using local initial guesses:
      !----------------------------------------------------------------------------------------------------                      
      
         ! make hydrodyn inputs consistant with elastodyn outputs 
         ! (do this because we're using outputs in the u vector):
      CALL Transfer_ED_to_HD(y_ED_input,  u_HD, MeshMapData, ErrStat2, ErrMsg2 ) ! get u_HD from y_ED_input
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
      
      u( 1: 3) = u_ED%PlatformPtMesh%Force(:,1) / p_FAST%UJacSclFact
      u( 4: 6) = u_ED%PlatformPtMesh%Moment(:,1) / p_FAST%UJacSclFact  
      u( 7: 9) = y_ED_input%PlatformPtMesh%TranslationAcc(:,1)
      u(10:12) = y_ED_input%PlatformPtMesh%RotationAcc(:,1)
            
      K = 0
      
      DO
         
         !-------------------------------------------------------------------------------------------------
         ! Calculate outputs at this_time, based on inputs at this_time
         !-------------------------------------------------------------------------------------------------
         
         CALL ED_CalcOutput( this_time, u_ED, p_ED, x_ED, xd_ED, z_ED, OtherSt_ED, y_ED, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                                 
         CALL HydroDyn_CalcOutput( this_time, u_HD, p_HD, x_HD, xd_HD, z_HD, OtherSt_HD, y_HD, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )      
            
         IF (ErrStat >= AbortErrLev) THEN
            CALL CleanUp()
            RETURN
         END IF
      
         IF ( K >= p_FAST%KMax ) EXIT
         
                                                            
         !-------------------------------------------------------------------------------------------------
         ! Calculate Jacobian: partial U/partial u:
         ! (note that we don't want to change u_ED or u_HD here)
         !-------------------------------------------------------------------------------------------------
         
         CALL U_ED_HD_Residual(y_ED, y_HD, u, Fn_U_Resid)   ! U_ED_HD_Residual checks for error
            IF (ErrStat >= AbortErrLev) THEN
               CALL CleanUp()
               RETURN
            END IF         
         
         IF ( calcJacobian ) THEN
            
            !...............................
            ! Get ElastoDyn's contribution:
            !...............................
            DO i=1,6 !call ED_CalcOutput
                  
               CALL ED_CopyInput(  u_ED, u_ED_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )            
               u_perturb = u            
               CALL Perturb_u( i, u_perturb, u_ED_perturb=u_ED_perturb, perturb=ThisPerturb ) ! perturb u and u_ED by ThisPerturb [routine sets ThisPerturb]
                  
               ! calculate outputs with perturbed inputs:
               CALL ED_CalcOutput( this_time, u_ED_perturb, p_ED, x_ED, xd_ED, z_ED, OtherSt_ED, y_ED_perturb, ErrStat2, ErrMsg2 ) !calculate y_ED_perturb
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )            
                  
                  
               CALL U_ED_HD_Residual(y_ED_perturb, y_HD, u_perturb, Fn_U_perturb) ! get this perturbation, U_perturb
                  IF ( ErrStat >= AbortErrLev ) RETURN ! U_ED_HD_Residual checks for error
                  
               IF (ErrStat >= AbortErrLev) THEN
                  CALL CleanUp()
                  RETURN
               END IF         
                  
                  
               MeshMapData%Jacobian_ED_SD_HD_BD(:,i) = (Fn_U_perturb - Fn_U_Resid) / ThisPerturb
                  
            END DO ! ElastoDyn contribution ( columns 1-6 )
               
            !...............................
            ! Get HydroDyn's contribution:
            !...............................               
            DO i=7,12 !call HD_CalcOutput
                  
               ! we want to perturb u_HD, but we're going to perturb the input y_ED and transfer that to HD to get u_HD
               CALL ED_CopyOutput( y_ED_input, y_ED_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )         
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                                   
               u_perturb = u            
               CALL Perturb_u( i, u_perturb, y_ED_perturb=y_ED_perturb, perturb=ThisPerturb ) ! perturb u and y_ED by ThisPerturb [routine sets ThisPerturb]
               CALL Transfer_ED_to_HD( y_ED_perturb, u_HD_perturb, MeshMapData, ErrStat2, ErrMsg2 ) ! get u_HD_perturb
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                                   
                  
               ! calculate outputs with perturbed inputs:
               CALL HydroDyn_CalcOutput( this_time, u_HD_perturb, p_HD, x_HD, xd_HD, z_HD, OtherSt_HD, y_HD_perturb, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  
                  
               CALL U_ED_HD_Residual(y_ED, y_HD_perturb, u_perturb, Fn_U_perturb) ! get this perturbation  ! U_ED_HD_Residual checks for error                      
                  IF ( ErrStat >= AbortErrLev ) THEN
                     CALL CleanUp()
                     RETURN 
                  END IF
                  
               MeshMapData%Jacobian_ED_SD_HD_BD(:,i) = (Fn_U_perturb - Fn_U_Resid) / ThisPerturb
                                                                  
            END DO ! HydroDyn contribution ( columns 7-12 )
               
#ifdef OUTPUT_ADDEDMASS  
   UnAM = -1
   CALL GetNewUnit( UnAM, ErrStat, ErrMsg )
   CALL OpenFOutFile( UnAM, TRIM(p_FAST%OutFileRoot)//'.AddedMassMatrix', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )               
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN 
      END IF

   AddedMassMatrix = MeshMapData%Jacobian_ED_SD_HD_BD(1:6,7:12) * p_FAST%UJacSclFact   
   CALL WrMatrix(AddedMassMatrix,UnAM, p_FAST%OutFmt)
   CLOSE( UnAM )
#endif   
#ifdef OUTPUT_JACOBIAN
   UnJac = -1
   CALL GetNewUnit( UnJac, ErrStat2, ErrMsg2 )
   CALL OpenFOutFile( UnJac, TRIM(p_FAST%OutFileRoot)//'.'//TRIM(num2lstr(this_time))//'.Jacobian2', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )               
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN 
      END IF
      
   CALL WrFileNR(UnJac, '  ')
      CALL WrFileNR(UnJac, ' ElastoDyn_Force_X') 
      CALL WrFileNR(UnJac, ' ElastoDyn_Force_Y') 
      CALL WrFileNR(UnJac, ' ElastoDyn_Force_Z') 
      CALL WrFileNR(UnJac, ' ElastoDyn_Moment_X') 
      CALL WrFileNR(UnJac, ' ElastoDyn_Moment_Y') 
      CALL WrFileNR(UnJac, ' ElastoDyn_Moment_Z') 
   
      CALL WrFileNR(UnJac, ' y_ED_TranslationAcc_X') 
      CALL WrFileNR(UnJac, ' y_ED_TranslationAcc_Y') 
      CALL WrFileNR(UnJac, ' y_ED_TranslationAcc_Z') 
      CALL WrFileNR(UnJac, ' y_ED_RotationAcc_X') 
      CALL WrFileNR(UnJac, ' y_ED_RotationAcc_Y') 
      CALL WrFileNR(UnJac, ' y_ED_RotationAcc_Z') 
   WRITE(UnJac,'()')    
      
   CALL WrMatrix(MeshMapData%Jacobian_ED_SD_HD_BD,UnJac, p_FAST%OutFmt)
   CLOSE( UnJac )      
#endif   
            
            
               ! Get the LU decomposition of this matrix using a LAPACK routine: 
               ! The result is of the form MeshMapDat%Jacobian_ED_SD_HD_BD = P * L * U 

            CALL LAPACK_getrf( M=NumInputs, N=NumInputs, A=MeshMapData%Jacobian_ED_SD_HD_BD, IPIV=MeshMapData%Jacobian_pivot, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN 
               END IF               
         END IF         
            
         !-------------------------------------------------------------------------------------------------
         ! Solve for delta u: Jac*u_delta = - Fn_U_Resid
         !  using the LAPACK routine 
         !-------------------------------------------------------------------------------------------------
         
         u_delta = -Fn_U_Resid
         CALL LAPACK_getrs( TRANS='N', N=NumInputs, A=MeshMapData%Jacobian_ED_SD_HD_BD, IPIV=MeshMapData%Jacobian_pivot, B=u_delta, &
                            ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN 
               END IF
      
         !-------------------------------------------------------------------------------------------------
         ! check for error, update inputs (u_ED and u_HD), and iterate again
         !-------------------------------------------------------------------------------------------------
                  
!         IF ( DOT_PRODUCT(u_delta, u_delta) <= TOL_Squared ) EXIT
         
         u = u + u_delta
                  
         u_ED%PlatformPtMesh%Force( :,1)               = u_ED%PlatformPtMesh%Force( :,1)               + u_delta( 1: 3) * p_FAST%UJacSclFact 
         u_ED%PlatformPtMesh%Moment(:,1)               = u_ED%PlatformPtMesh%Moment(:,1)               + u_delta( 4: 6) * p_FAST%UJacSclFact
         y_ED_input%PlatformPtMesh%TranslationAcc(:,1) = y_ED_input%PlatformPtMesh%TranslationAcc(:,1) + u_delta( 7: 9)
         y_ED_input%PlatformPtMesh%RotationAcc(   :,1) = y_ED_input%PlatformPtMesh%RotationAcc(   :,1) + u_delta(10:12)
                  
         CALL Transfer_ED_to_HD( y_ED_input, u_HD, MeshMapData, ErrStat2, ErrMsg2 ) ! get u_HD with u_delta changes
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
         
         K = K + 1
         
      END DO ! K
            
   
      CALL CleanUp()
      
CONTAINS                                 
   !...............................................................................................................................
   SUBROUTINE Perturb_u( n, u_perturb, u_ED_perturb, y_ED_perturb, perturb )
   ! This routine perturbs the nth element of the u array (and ED input/output it corresponds to)
   !...............................................................................................................................
!   REAL( ReKi ),                       INTENT(IN)    :: this_U(NumInputs)
   INTEGER( IntKi )                  , INTENT(IN   ) :: n
   REAL( ReKi )                      , INTENT(INOUT) :: u_perturb(numInputs)
   TYPE(ED_InputType) , OPTIONAL     , INTENT(INOUT) :: u_ED_perturb           ! System inputs   (needed only when 1 <= n <=  6)
   TYPE(ED_OutputType), OPTIONAL     , INTENT(INOUT) :: y_ED_perturb           ! System outputs  (needed only when 7 <= n <= 12)
   REAL( ReKi )                      , INTENT(  OUT) :: perturb
   
   if ( n <= 6 ) then ! ED u
   
      if ( n <= 3 ) then         
         perturb = GetPerturb( u_ED_perturb%PlatformPtMesh%Force(n   ,1) )         
         u_ED_perturb%PlatformPtMesh%Force(n   ,1) = u_ED_perturb%PlatformPtMesh%Force(n   ,1) + perturb * p_FAST%UJacSclFact 
      else
         perturb = GetPerturb( u_ED_perturb%PlatformPtMesh%Moment(n-3,1) )         
         u_ED_perturb%PlatformPtMesh%Moment(n-3,1) = u_ED_perturb%PlatformPtMesh%Moment(n-3,1) + perturb * p_FAST%UJacSclFact 
      end if
                  
   else ! ED y = HD u
      
      if ( n <= 9 ) then         
         perturb = GetPerturb( y_ED_perturb%PlatformPtMesh%TranslationAcc(n-6,1) )         
         y_ED_perturb%PlatformPtMesh%TranslationAcc(n-6,1) = y_ED_perturb%PlatformPtMesh%TranslationAcc(n-6,1) + perturb
      else
         perturb = GetPerturb( y_ED_perturb%PlatformPtMesh%RotationAcc(n-9,1) )         
         y_ED_perturb%PlatformPtMesh%RotationAcc(   n-9,1) = y_ED_perturb%PlatformPtMesh%RotationAcc(   n-9,1) + perturb
      end if
                  
   end if
           
   u_perturb(n) = u_perturb(n) + perturb
   
        
   END SUBROUTINE Perturb_u
   !...............................................................................................................................
   SUBROUTINE U_ED_HD_Residual( y_ED2, y_HD2, u_IN, U_Resid)
   !...............................................................................................................................
                                  
   TYPE(ED_OutputType)               , INTENT(IN   ) :: y_ED2                  ! System outputs
   TYPE(HydroDyn_OutputType)         , INTENT(IN   ) :: y_HD2                  ! System outputs
   REAL(ReKi)                        , INTENT(IN   ) :: u_in(NumInputs)
   REAL(ReKi)                        , INTENT(  OUT) :: U_Resid(NumInputs)


   
   
      !   ! Transfer motions:
      
   !..................
   ! Set mooring line inputs (which don't have acceleration fields)
   !..................
   
      IF ( p_FAST%CompMooring == Module_MAP ) THEN
         
         ! note: MAP_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)      
         CALL MAP_InputSolve( u_map, y_ED2, MeshMapData, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )       
                                 
         CALL Transfer_Point_to_Point( y_MAP%PtFairleadLoad, MeshMapData%u_ED_PlatformPtMesh_2, MeshMapData%Mooring_P_2_ED_P, ErrStat2, ErrMsg2, u_MAP%PtFairDisplacement, y_ED2%PlatformPtMesh ) !u_MAP and y_ED contain the displacements needed for moment calculations
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )               
                 
      ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
         
         ! note: MD_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)      
         CALL MD_InputSolve( u_MD, y_ED2, MeshMapData, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )       
                 
         CALL Transfer_Point_to_Point( y_MD%PtFairleadLoad, MeshMapData%u_ED_PlatformPtMesh_2, MeshMapData%Mooring_P_2_ED_P, ErrStat2, ErrMsg2, u_MD%PtFairleadDisplacement, y_ED2%PlatformPtMesh ) !u_MD and y_ED contain the displacements needed for moment calculations
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )              
            
      ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
         
         ! note: FEAM_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)      
         CALL FEAM_InputSolve( u_FEAM, y_ED2, MeshMapData, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )       
                 
         CALL Transfer_Point_to_Point( y_FEAM%PtFairleadLoad, MeshMapData%u_ED_PlatformPtMesh_2, MeshMapData%Mooring_P_2_ED_P, ErrStat2, ErrMsg2, u_FEAM%PtFairleadDisplacement, y_ED2%PlatformPtMesh ) !u_FEAM and y_ED contain the displacements needed for moment calculations
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )              
            
      ELSE
         
         MeshMapData%u_ED_PlatformPtMesh_2%Force  = 0.0_ReKi
         MeshMapData%u_ED_PlatformPtMesh_2%Moment = 0.0_ReKi
         
      END IF        

   ! we use copies of the input meshes (we don't need to update values in the original data structures):            
      
!bjj: why don't we update u_HD2 here? shouldn't we update before using it to transfer the loads?
   
      CALL Transfer_Point_to_Point( y_ED2%PlatformPtMesh, MeshMapData%u_HD_Mesh, MeshMapData%ED_P_2_HD_W_P, ErrStat2, ErrMsg2) 
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         
  
         
         ! we're mapping loads, so we also need the sibling meshes' displacements:
      CALL Transfer_Point_to_Point( y_HD2%AllHdroOrigin, MeshMapData%u_ED_PlatformPtMesh, MeshMapData%HD_W_P_2_ED_P, ErrStat2, ErrMsg2, MeshMapData%u_HD_Mesh, y_ED2%PlatformPtMesh) !u_HD and u_mapped_positions contain the displaced positions for load calculations
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      MeshMapData%u_ED_PlatformPtMesh%Force  = MeshMapData%u_ED_PlatformPtMesh%Force  + MeshMapData%u_ED_PlatformPtMesh_2%Force
      MeshMapData%u_ED_PlatformPtMesh%Moment = MeshMapData%u_ED_PlatformPtMesh%Moment + MeshMapData%u_ED_PlatformPtMesh_2%Moment
                                       
            
      U_Resid( 1: 3) = u_in( 1: 3) - MeshMapData%u_ED_PlatformPtMesh%Force(:,1) / p_FAST%UJacSclFact
      U_Resid( 4: 6) = u_in( 4: 6) - MeshMapData%u_ED_PlatformPtMesh%Moment(:,1) / p_FAST%UJacSclFact      
      U_Resid( 7: 9) = u_in( 7: 9) - y_ED2%PlatformPtMesh%TranslationAcc(:,1)
      U_Resid(10:12) = u_in(10:12) - y_ED2%PlatformPtMesh%RotationAcc(:,1)
   
            
   END SUBROUTINE U_ED_HD_Residual   
   !...............................................................................................................................
   SUBROUTINE CleanUp()
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(1024)            :: ErrMsg3     ! The error message (ErrMsg)
                  
      CALL ED_DestroyOutput(y_ED_input, ErrStat3, ErrMsg3 )
         IF (ErrStat3 /= ErrID_None) CALL WrScr(RoutineName//'/ED_DestroyOutput: '//TRIM(ErrMsg3) )
         
      IF ( calcJacobian ) THEN
         CALL ED_DestroyInput( u_ED_perturb, ErrStat3, ErrMsg3 )
            IF (ErrStat3 /= ErrID_None) CALL WrScr(RoutineName//'/ED_DestroyInput: '//TRIM(ErrMsg3) )
         CALL ED_DestroyOutput(y_ED_perturb, ErrStat3, ErrMsg3 )
            IF (ErrStat3 /= ErrID_None) CALL WrScr(RoutineName//'/ED_DestroyOutput: '//TRIM(ErrMsg3) )
         
         CALL HydroDyn_DestroyInput( u_HD_perturb, ErrStat3, ErrMsg3 )
            IF (ErrStat3 /= ErrID_None) CALL WrScr(RoutineName//'/HydroDyn_DestroyInput: '//TRIM(ErrMsg3) )
         CALL HydroDyn_DestroyOutput(y_HD_perturb, ErrStat3, ErrMsg3 )
            IF (ErrStat3 /= ErrID_None) CALL WrScr(RoutineName//'/HydroDyn_DestroyOutput: '//TRIM(ErrMsg3) )                        
      END IF
      
   
      
   END SUBROUTINE CleanUp
   !...............................................................................................................................
END SUBROUTINE ED_HD_InputOutputSolve
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ED_SD_HD_BD_InputOutputSolve(  this_time, p_FAST, calcJacobian &
                                     , u_ED, p_ED, x_ED, xd_ED, z_ED, OtherSt_ED, y_ED &
                                     , u_SD, p_SD, x_SD, xd_SD, z_SD, OtherSt_SD, y_SD & 
                                     , u_HD, p_HD, x_HD, xd_HD, z_HD, OtherSt_HD, y_HD & 
                                     , u_BD, p_BD, x_BD, xd_BD, z_BD, OtherSt_BD, y_BD & 
                                     , u_MAP,  y_MAP  &
                                     , u_FEAM, y_FEAM & 
                                     , u_MD,   y_MD   & 
                                     , u_IceF, y_IceF & 
                                     , u_IceD, y_IceD & 
                                     , MeshMapData , ErrStat, ErrMsg )
! This routine performs the Input-Output solve for ED, SD, and HD.
! Note that this has been customized for the physics in the problems and is not a general solution.
!..................................................................................................................................

   USE ElastoDyn
   USE SubDyn
   USE HydroDyn
   USE BeamDyn

      ! Passed variables

   REAL(DbKi)                        , INTENT(IN   ) :: this_time                 ! The current simulation time (actual or time of prediction)
   TYPE(FAST_ParameterType)          , INTENT(IN   ) :: p_FAST                    ! Glue-code simulation parameters
   LOGICAL                           , INTENT(IN   ) :: calcJacobian              ! Should we calculate Jacobians this time?
                                                                                  
      !ElastoDyn:                                                                 
   TYPE(ED_ContinuousStateType)      , INTENT(IN   ) :: x_ED                      ! Continuous states
   TYPE(ED_DiscreteStateType)        , INTENT(IN   ) :: xd_ED                     ! Discrete states
   TYPE(ED_ConstraintStateType)      , INTENT(IN   ) :: z_ED                      ! Constraint states
   TYPE(ED_OtherStateType)           , INTENT(INOUT) :: OtherSt_ED                ! Other/optimization states
   TYPE(ED_ParameterType)            , INTENT(IN   ) :: p_ED                      ! Parameters
   TYPE(ED_InputType)                , INTENT(INOUT) :: u_ED                      ! System inputs
   TYPE(ED_OutputType)               , INTENT(INOUT) :: y_ED                      ! System outputs
         
      !BeamDyn (one instance per blade):                                                                 
   TYPE(BD_ContinuousStateType)      , INTENT(IN   ) :: x_BD(:)                   ! Continuous states
   TYPE(BD_DiscreteStateType)        , INTENT(IN   ) :: xd_BD(:)                  ! Discrete states
   TYPE(BD_ConstraintStateType)      , INTENT(IN   ) :: z_BD(:)                   ! Constraint states
   TYPE(BD_OtherStateType)           , INTENT(INOUT) :: OtherSt_BD(:)             ! Other/optimization states
   TYPE(BD_ParameterType)            , INTENT(IN   ) :: p_BD(:)                   ! Parameters
   TYPE(BD_InputType)                , INTENT(INOUT) :: u_BD(:)                   ! System inputs
   TYPE(BD_OutputType)               , INTENT(INOUT) :: y_BD(:)                   ! System outputs
   
      !SubDyn:                                                                    
   TYPE(SD_ContinuousStateType)      , INTENT(IN   ) :: x_SD                      ! Continuous states
   TYPE(SD_DiscreteStateType)        , INTENT(IN   ) :: xd_SD                     ! Discrete states
   TYPE(SD_ConstraintStateType)      , INTENT(IN   ) :: z_SD                      ! Constraint states
   TYPE(SD_OtherStateType)           , INTENT(INOUT) :: OtherSt_SD                ! Other/optimization states
   TYPE(SD_ParameterType)            , INTENT(IN   ) :: p_SD                      ! Parameters
   TYPE(SD_InputType)                , INTENT(INOUT) :: u_SD                      ! System inputs
   TYPE(SD_OutputType)               , INTENT(INOUT) :: y_SD                      ! System outputs
          
      !HydroDyn: 
   TYPE(HydroDyn_ContinuousStateType), INTENT(IN   ) :: x_HD                      ! Continuous states
   TYPE(HydroDyn_DiscreteStateType)  , INTENT(IN   ) :: xd_HD                     ! Discrete states
   TYPE(HydroDyn_ConstraintStateType), INTENT(IN   ) :: z_HD                      ! Constraint states
   TYPE(HydroDyn_OtherStateType)     , INTENT(INOUT) :: OtherSt_HD                ! Other/optimization states
   TYPE(HydroDyn_ParameterType)      , INTENT(IN   ) :: p_HD                      ! Parameters
   TYPE(HydroDyn_InputType)          , INTENT(INOUT) :: u_HD                      ! System inputs
   TYPE(HydroDyn_OutputType)         , INTENT(INOUT) :: y_HD                      ! System outputs
   
      ! MAP/FEAM/MoorDyn/IceFloe/IceDyn:
   TYPE(MAP_OutputType),              INTENT(IN   )  :: y_MAP                     ! MAP outputs 
   TYPE(MAP_InputType),               INTENT(INOUT)  :: u_MAP                     ! MAP inputs (INOUT just because I don't want to use another tempoarary mesh and we'll overwrite this later)
   TYPE(FEAM_OutputType),             INTENT(IN   )  :: y_FEAM                    ! FEAM outputs  
   TYPE(FEAM_InputType),              INTENT(INOUT)  :: u_FEAM                    ! FEAM inputs (INOUT just because I don't want to use another tempoarary mesh and we'll overwrite this later)
   TYPE(MD_OutputType),               INTENT(IN   )  :: y_MD                      ! MoorDyn outputs  
   TYPE(MD_InputType),                INTENT(INOUT)  :: u_MD                      ! MoorDyn inputs (INOUT just because I don't want to use another tempoarary mesh and we'll overwrite this later)
   TYPE(IceFloe_OutputType),          INTENT(IN   )  :: y_IceF                    ! IceFloe outputs  
   TYPE(IceFloe_InputType),           INTENT(INOUT)  :: u_IceF                    ! IceFloe inputs (INOUT just because I don't want to use another tempoarary mesh and we'll overwrite this later)
   TYPE(IceD_OutputType),             INTENT(IN   )  :: y_IceD(:)                 ! IceDyn outputs  
   TYPE(IceD_InputType),              INTENT(INOUT)  :: u_IceD(:)                 ! IceDyn inputs (INOUT just because I don't want to use another tempoarary mesh and we'll overwrite this later)
      
   TYPE(FAST_ModuleMapType)          , INTENT(INOUT) :: MeshMapData
   INTEGER(IntKi)                    , INTENT(  OUT) :: ErrStat                   ! Error status of the operation
   CHARACTER(*)                      , INTENT(  OUT) :: ErrMsg                    ! Error message if ErrStat /= ErrID_None

   ! Local variables:
   REAL(ReKi),                             PARAMETER :: TOL_Squared = (1.0E-4)**2 !not currently used because KMax = 1
   REAL(ReKi)                                        :: ThisPerturb               ! an arbitrary perturbation (these are linear, so it shouldn't matter)
   
   CHARACTER(*),                           PARAMETER :: RoutineName = 'ED_SD_HD_BD_InputOutputSolve'
   
!bjj: store these so that we don't reallocate every time?   
   REAL(ReKi)                                        :: u(           p_FAST%SizeJac_ED_SD_HD_BD(1))   ! size of loads/accelerations passed between the 4 modules
   REAL(ReKi)                                        :: u_perturb(   p_FAST%SizeJac_ED_SD_HD_BD(1))   ! size of loads/accelerations passed between the 4 modules
   REAL(ReKi)                                        :: u_delta(     p_FAST%SizeJac_ED_SD_HD_BD(1))   ! size of loads/accelerations passed between the 4 modules
   REAL(ReKi)                                        :: Fn_U_perturb(p_FAST%SizeJac_ED_SD_HD_BD(1))   ! value of U with perturbations
   REAL(ReKi)                                        :: Fn_U_Resid(  p_FAST%SizeJac_ED_SD_HD_BD(1))   ! Residual of U
                                                                                           
   TYPE(ED_InputType)                                :: u_ED_perturb              ! Perturbed system inputs
   TYPE(ED_OutputType)                               :: y_ED_perturb              ! Perturbed system outputs
   TYPE(SD_InputType)                                :: u_SD_perturb              ! Perturbed system inputs
   TYPE(SD_OutputType)                               :: y_SD_perturb              ! Perturbed system outputs
   TYPE(HydroDyn_InputType)                          :: u_HD_perturb              ! Perturbed system inputs
   TYPE(HydroDyn_OutputType)                         :: y_HD_perturb              ! Perturbed system outputs
   TYPE(BD_InputType)                                :: u_BD_perturb              ! Perturbed system inputs
   TYPE(BD_OutputType), ALLOCATABLE                  :: y_BD_perturb(:)           ! Perturbed system outputs
                                                                                  
                                                                                  
   INTEGER(IntKi)                                    :: i,j                       ! loop counters (jacobian column number)
   INTEGER(IntKi)                                    :: nb                        ! loop counter (blade number)
   INTEGER(IntKi)                                    :: K                         ! Input-output-solve iteration counter
   INTEGER(IntKi)                                    :: ErrStat2                  ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                              :: ErrMsg2                   ! temporary Error message if ErrStat /= ErrID_None
   
#ifdef OUTPUT_ADDEDMASS   
   REAL(ReKi)                                        :: AddedMassMatrix(6,6)
   INTEGER                                           :: UnAM
   INTEGER                                           :: AMIndx   
#endif
#ifdef OUTPUT_JACOBIAN
   INTEGER                                           :: UnJac
   INTEGER                                           :: TmpIndx   
#endif
   
      
   ! Note: p_FAST%UJacSclFact is a scaling factor that gets us similar magnitudes between loads and accelerations...
 
!bjj: note, that this routine may have a problem if there is remapping done
    
   ErrStat = ErrID_None
   ErrMsg  = ""
   
      !----------------------------------------------------------------------------------------------------
      ! Some record keeping stuff:
      !----------------------------------------------------------------------------------------------------      
                  
         ! Local copies for perturbing inputs and outputs (computing Jacobian):
      IF ( calcJacobian ) THEN         
         
         CALL ED_CopyInput(  u_ED, u_ED_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, 'u_ED_perturb:'//ErrMsg2, ErrStat, ErrMsg, RoutineName  )
         CALL ED_CopyOutput( y_ED, y_ED_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )         
            CALL SetErrStat( ErrStat2, 'y_ED_perturb:'//ErrMsg2, ErrStat, ErrMsg, RoutineName  )
            
         IF ( p_FAST%CompSub == Module_SD ) THEN   
            CALL SD_CopyInput(  u_SD, u_SD_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )           
               CALL SetErrStat( ErrStat2, 'u_SD_perturb:'//ErrMsg2, ErrStat, ErrMsg, RoutineName  )
            CALL SD_CopyOutput( y_SD, y_SD_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )  
               CALL SetErrStat( ErrStat2, 'y_SD_perturb:'//ErrMsg2, ErrStat, ErrMsg, RoutineName  )
         END IF
            
         IF ( p_FAST%CompHydro == Module_HD ) THEN            
            CALL HydroDyn_CopyInput(  u_HD, u_HD_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )           
               CALL SetErrStat( ErrStat2, 'u_HD_perturb:'//ErrMsg2, ErrStat, ErrMsg, RoutineName  )
            CALL HydroDyn_CopyOutput( y_HD, y_HD_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )  
               CALL SetErrStat( ErrStat2, 'y_HD_perturb:'//ErrMsg2, ErrStat, ErrMsg, RoutineName  )
         END IF
         
         IF ( p_FAST%CompElast == Module_BD ) THEN    
            ALLOCATE( y_BD_perturb(p_FAST%nBeams) , STAT = ErrStat2 )
            if (ErrStat2 /= 0) then
               call SetErrStat( ErrID_Fatal, 'Error allocating y_BD_perturb.', ErrStat, ErrMsg, RoutineName )
               call CleanUp()
               return
            end if
            
            !bjj: if this mesh could have different number of nodes per instance, we'd need an array. (but, it's a single point) 
            CALL BD_CopyInput(  u_BD(1), u_BD_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )           
               CALL SetErrStat( ErrStat2, 'u_BD_perturb:'//ErrMsg2, ErrStat, ErrMsg, RoutineName  )
            do nb=1,p_FAST%nBeams
               CALL BD_CopyOutput( y_BD(nb), y_BD_perturb(nb), MESH_NEWCOPY, ErrStat2, ErrMsg2 )  
                  CALL SetErrStat( ErrStat2, 'y_BD_perturb:'//ErrMsg2, ErrStat, ErrMsg, RoutineName  )
            end do
            
         END IF
         
         IF (ErrStat >= AbortErrLev) THEN
            CALL CleanUp()
            RETURN
         END IF
         
      END IF
         
      !----------------------------------------------------------------------------------------------------
      ! set up u vector, using local initial guesses:
      !---------------------------------------------------------------------------------------------------- 
      
      ! we need BeamDyn input mesh to be an array of meshes, so first we'll copy into temporary storage:
      DO nb=1,p_FAST%nBeams
         call MeshCopy( u_BD(nb)%RootMotion, MeshMapData%u_BD_RootMotion(nb), MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
      END DO
      
      CALL Create_ED_SD_HD_BD_UVector(u, u_ED%PlatformPtMesh, u_SD%TPMesh, u_SD%LMesh, u_HD%Morison%LumpedMesh, &
               u_HD%Morison%DistribMesh, u_HD%Mesh, u_ED%HubPtLoad, MeshMapData%u_BD_RootMotion, p_FAST )
                  
      K = 0
      
      DO
         
         !-------------------------------------------------------------------------------------------------
         ! Calculate outputs at this_time, based on inputs at this_time
         !-------------------------------------------------------------------------------------------------
         
         CALL ED_CalcOutput( this_time, u_ED, p_ED, x_ED, xd_ED, z_ED, OtherSt_ED, y_ED, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
                                 
         IF ( p_FAST%CompSub == Module_SD ) THEN            
            CALL SD_CalcOutput( this_time, u_SD, p_SD, x_SD, xd_SD, z_SD, OtherSt_SD, y_SD, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
         END IF
            
         IF ( p_FAST%CompHydro == Module_HD ) THEN 
            CALL HydroDyn_CalcOutput( this_time, u_HD, p_HD, x_HD, xd_HD, z_HD, OtherSt_HD, y_HD, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
         END IF
         
         IF ( p_FAST%CompElast == Module_BD ) THEN 
            do nb=1,p_FAST%nBeams
               CALL BD_CalcOutput( this_time, u_BD(nb), p_BD(nb), x_BD(nb), xd_BD(nb), z_BD(nb), OtherSt_BD(nb), y_BD(nb), ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
                  
            end do
            
#ifdef DEBUG_BEAMDYN            
            if (k>=0) then
               do nb=1,p_FAST%nBeams
                  write(74,'(ES10.3E2,2(1x,i3),100(1x,ES10.3E2))') this_time, k, nb, &
                     u_BD(nb)%RootMotion%Orientation,    u_BD(nb)%RootMotion%TranslationDisp, &
                     u_BD(nb)%RootMotion%TranslationVel, u_BD(nb)%RootMotion%RotationVel, &
                     u_BD(nb)%RootMotion%TranslationAcc, u_BD(nb)%RootMotion%RotationAcc, &
                     y_BD(nb)%ReactionForce%Force,       y_BD(nb)%ReactionForce%Moment, &
                     y_BD(nb)%BldMotion%TranslationDisp(:,y_BD(nb)%BldMotion%Nnodes)
               end do
            end if
#endif            
            
         END IF
         
         
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN      
         END IF
                  
         
         IF ( K >= p_FAST%KMax ) EXIT
         
                                                            
         !-------------------------------------------------------------------------------------------------
         ! Calculate Jacobian: partial U/partial u:
         ! (note that we don't want to change u_ED, u_SD, u_HD, or u_BD here)
         !-------------------------------------------------------------------------------------------------
         
         CALL U_ED_SD_HD_BD_Residual(y_ED, y_SD, y_HD, y_BD, u, Fn_U_Resid)    !May set errors here...              
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN      
            END IF
         
         IF ( calcJacobian ) THEN
            
            !...............................
            ! Get ElastoDyn's contribution:
            !...............................
            DO i=1,p_FAST%SizeJac_ED_SD_HD_BD(2) !call ED_CalcOutput
                  
               ! perturb u_ED:
               CALL ED_CopyInput(  u_ED, u_ED_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
               u_perturb = u            
               CALL Perturb_u_ED_SD_HD_BD( p_FAST, MeshMapData%Jac_u_indx, i, u_perturb, u_ED_perturb=u_ED_perturb, perturb=ThisPerturb ) ! perturb u and u_ED by ThisPerturb [routine sets ThisPerturb]
                  
               ! calculate outputs with perturbed inputs:
               CALL ED_CalcOutput( this_time, u_ED_perturb, p_ED, x_ED, xd_ED, z_ED, OtherSt_ED, y_ED_perturb, ErrStat2, ErrMsg2 ) !calculate y_ED_perturb
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
                  
                  
               CALL U_ED_SD_HD_BD_Residual(y_ED_perturb, y_SD, y_HD, y_BD, u_perturb, Fn_U_perturb) ! get this perturbation, U_perturb
               
               
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN      
               END IF
            
                MeshMapData%Jacobian_ED_SD_HD_BD(:,i) = (Fn_U_perturb - Fn_U_Resid) / ThisPerturb
                  
            END DO ! ElastoDyn contribution ( columns 1-p_FAST%SizeJac_ED_SD_HD_BD(2) )
               
            i = p_FAST%SizeJac_ED_SD_HD_BD(2)
            !...............................
            ! Get SubDyn's contribution:  (note if p_FAST%CompSub /= Module_SD, SizeJac_ED_SD_HD_BD(3) = 0)
            !...............................               
            DO j=1,p_FAST%SizeJac_ED_SD_HD_BD(3) !call SD_CalcOutput
               i = i + 1 ! i = j + p_FAST%SizeJac_ED_SD_HD_BD(2)
                              
               ! perturb u_SD:
               CALL SD_CopyInput(  u_SD, u_SD_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
               u_perturb = u            
               CALL Perturb_u_ED_SD_HD_BD( p_FAST, MeshMapData%Jac_u_indx, i, u_perturb, u_SD_perturb=u_SD_perturb, perturb=ThisPerturb ) ! perturb u and u_SD by ThisPerturb [routine sets ThisPerturb]
                  
               ! calculate outputs with perturbed inputs:
               CALL SD_CalcOutput( this_time, u_SD_perturb, p_SD, x_SD, xd_SD, z_SD, OtherSt_SD, y_SD_perturb, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
                  
                  
               CALL U_ED_SD_HD_BD_Residual(y_ED, y_SD_perturb, y_HD, y_BD, u_perturb, Fn_U_perturb) ! get this perturbation    
               
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN      
               END IF                  
               
               MeshMapData%Jacobian_ED_SD_HD_BD(:,i) = (Fn_U_perturb - Fn_U_Resid) / ThisPerturb
                                    
            END DO ! SubDyn contribution
            
                  
            !...............................
            ! Get HydroDyn's contribution: (note if p_FAST%CompHydro /= Module_HD, SizeJac_ED_SD_HD_BD(4) = 0)
            !...............................             
            DO j=1,p_FAST%SizeJac_ED_SD_HD_BD(4) !call HydroDyn_CalcOutput            
               i = i + 1 ! i = j + p_FAST%SizeJac_ED_SD_HD_BD(2) + p_FAST%SizeJac_ED_SD_HD_BD(3) 

               ! perturb u_HD:
               CALL HydroDyn_CopyInput(  u_HD, u_HD_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
               u_perturb = u            
               CALL Perturb_u_ED_SD_HD_BD( p_FAST, MeshMapData%Jac_u_indx, i, u_perturb, u_HD_perturb=u_HD_perturb, perturb=ThisPerturb ) ! perturb u and u_HD by ThisPerturb [routine sets ThisPerturb]
                  
               ! calculate outputs with perturbed inputs:
               CALL HydroDyn_CalcOutput( this_time, u_HD_perturb, p_HD, x_HD, xd_HD, z_HD, OtherSt_HD, y_HD_perturb, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
                  
               CALL U_ED_SD_HD_BD_Residual(y_ED, y_SD, y_HD_perturb, y_BD, u_perturb, Fn_U_perturb) ! get this perturbation  
               
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN      
               END IF                  
               
               MeshMapData%Jacobian_ED_SD_HD_BD(:,i) = (Fn_U_perturb - Fn_U_Resid) / ThisPerturb
                                             
            END DO !HydroDyn contribution
            
            !...............................
            ! Get BeamDyn's contribution: (note if p_FAST%CompElast /= Module_BD, SizeJac_ED_SD_HD_BD(5) = 0)
            !............................... 
            DO nb=1,p_FAST%nBeams
               CALL BD_CopyOutput(y_BD(nb),y_BD_perturb(nb),MESH_UPDATECOPY,ErrStat2,ErrMsg2)
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
            END DO
                        
            DO nb=1,p_FAST%nBeams
               
               ! make sure we perturb the outputs from only the current blade (overwrite the previous perturbation with y_BD values)
               if (nb > 1) then               
                  CALL BD_CopyOutput(  y_BD(nb-1), y_BD_perturb(nb-1), MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
                     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
               end if
               
               
               DO j=1,p_FAST%SizeJac_ED_SD_HD_BD(4+nb) !call BD_CalcOutput            
                  i = i + 1 ! i = j + sum(p_FAST%SizeJac_ED_SD_HD_BD(2:3+nb))
                  ! perturb u_BD:
                  CALL BD_CopyInput(  u_BD(nb), u_BD_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
                     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
                  u_perturb = u            
                  CALL Perturb_u_ED_SD_HD_BD( p_FAST, MeshMapData%Jac_u_indx, i, u_perturb, u_BD_perturb=u_BD_perturb, perturb=ThisPerturb ) ! perturb u and u_HD by ThisPerturb [routine sets ThisPerturb]
                  
                  ! calculate outputs with perturbed inputs:
                  CALL BD_CalcOutput( this_time, u_BD_perturb, p_BD(nb), x_BD(nb), xd_BD(nb), z_BD(nb), OtherSt_BD(nb), y_BD_perturb(nb), ErrStat2, ErrMsg2 )
                     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
                  CALL U_ED_SD_HD_BD_Residual(y_ED, y_SD, y_HD, y_BD_perturb, u_perturb, Fn_U_perturb) ! get this perturbation  
               
                  IF ( ErrStat >= AbortErrLev ) THEN
                     CALL CleanUp()
                     RETURN      
                  END IF                  
               
                  MeshMapData%Jacobian_ED_SD_HD_BD(:,i) = (Fn_U_perturb - Fn_U_Resid) / ThisPerturb
               END DO
               
               
            END DO !BeamDyn contribution          
            
            
#ifdef OUTPUT_ADDEDMASS  
IF (p_FAST%CompHydro == Module_HD ) THEN
   UnAM = -1
   CALL GetNewUnit( UnAM, ErrStat2, ErrMsg2 )
   CALL OpenFOutFile( UnAM, TRIM(p_FAST%OutFileRoot)//'.AddedMassMatrix', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
      IF ( ErrStat >= AbortErrLev ) RETURN               
   
   AMIndx = p_FAST%SizeJac_ED_SD_HD_BD(1) - 5 - p_FAST%SizeJac_ED_SD_HD_BD(5) &
          - p_FAST%SizeJac_ED_SD_HD_BD(6) - p_FAST%SizeJac_ED_SD_HD_BD(7) !the start of the HydroDyn Mesh inputs in the Jacobian
   AddedMassMatrix = MeshMapData%Jacobian_ED_SD_HD_BD(1:6,AMIndx:(AMIndx+5)) * p_FAST%UJacSclFact   
   CALL WrMatrix(AddedMassMatrix,UnAM, p_FAST%OutFmt)
   CLOSE( UnAM )
END IF
#endif
#ifdef OUTPUT_JACOBIAN
   UnJac = -1
   CALL GetNewUnit( UnJac, ErrStat2, ErrMsg2 )
   CALL OpenFOutFile( UnJac, TRIM(p_FAST%OutFileRoot)//'.'//TRIM(num2lstr(this_time))//'.Jacobian', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
      IF ( ErrStat >= AbortErrLev ) RETURN               
      
   CALL WrFileNR(UnJac, '  ')
   IF (p_FAST%CompHydro == Module_HD .or. p_FAST%CompSub == Module_SD) then   
      CALL WrFileNR(UnJac, ' ElastoDyn_Ptfm_Force_X') 
      CALL WrFileNR(UnJac, ' ElastoDyn_Ptfm_Force_Y') 
      CALL WrFileNR(UnJac, ' ElastoDyn_Ptfm_Force_Z') 
      CALL WrFileNR(UnJac, ' ElastoDyn_Ptfm_Moment_X') 
      CALL WrFileNR(UnJac, ' ElastoDyn_Ptfm_Moment_Y') 
      CALL WrFileNR(UnJac, ' ElastoDyn_Ptfm_Moment_Z') 
   end if
   
   IF (p_FAST%CompElast == Module_BD) then
      CALL WrFileNR(UnJac, ' ElastoDyn_hub_Force_X') 
      CALL WrFileNR(UnJac, ' ElastoDyn_hub_Force_Y') 
      CALL WrFileNR(UnJac, ' ElastoDyn_hub_Force_Z') 
      CALL WrFileNR(UnJac, ' ElastoDyn_hub_Moment_X') 
      CALL WrFileNR(UnJac, ' ElastoDyn_hub_Moment_Y') 
      CALL WrFileNR(UnJac, ' ElastoDyn_hub_Moment_Z')       
   END IF
   
   
   DO TmpIndx=1,u_SD%TPMesh%NNodes
      CALL WrFileNR(UnJac, ' SD_TPMesh_TranslationAcc_X_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' SD_TPMesh_TranslationAcc_Y_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' SD_TPMesh_TranslationAcc_Z_'//TRIM(Num2LStr(TmpIndx))) 
   END DO

   DO TmpIndx=1,u_SD%TPMesh%NNodes
      CALL WrFileNR(UnJac, ' SD_TPMesh_RotationAcc_X_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' SD_TPMesh_RotationAcc_Y_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' SD_TPMesh_RotationAcc_Z_'//TRIM(Num2LStr(TmpIndx))) 
   END DO
            
   IF ( p_FAST%CompHydro == Module_HD ) THEN   ! this SD mesh linked only when HD is enabled
      DO TmpIndx=1,u_SD%LMesh%NNodes
         CALL WrFileNR(UnJac, ' SD_LMesh_Force_X_'//TRIM(Num2LStr(TmpIndx))) 
         CALL WrFileNR(UnJac, ' SD_LMesh_Force_Y_'//TRIM(Num2LStr(TmpIndx))) 
         CALL WrFileNR(UnJac, ' SD_LMesh_Force_Z_'//TRIM(Num2LStr(TmpIndx))) 
      END DO      
      DO TmpIndx=1,u_SD%LMesh%NNodes
         CALL WrFileNR(UnJac, ' SD_LMesh_Moment_X_'//TRIM(Num2LStr(TmpIndx))) 
         CALL WrFileNR(UnJac, ' SD_LMesh_Moment_Y_'//TRIM(Num2LStr(TmpIndx))) 
         CALL WrFileNR(UnJac, ' SD_LMesh_Moment_Z_'//TRIM(Num2LStr(TmpIndx))) 
      END DO                  
   END IF
   
   DO TmpIndx=1,u_HD%Morison%LumpedMesh%NNodes
      CALL WrFileNR(UnJac, ' HD_M_Lumped_TranslationAcc_X_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_M_Lumped_TranslationAcc_Y_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_M_Lumped_TranslationAcc_Z_'//TRIM(Num2LStr(TmpIndx))) 
   END DO   
   DO TmpIndx=1,u_HD%Morison%LumpedMesh%NNodes
      CALL WrFileNR(UnJac, ' HD_M_Lumped_RotationAcc_X_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_M_Lumped_RotationAcc_Y_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_M_Lumped_RotationAcc_Z_'//TRIM(Num2LStr(TmpIndx))) 
   END DO   

   DO TmpIndx=1,u_HD%Morison%DistribMesh%NNodes
      CALL WrFileNR(UnJac, ' HD_M_Distrib_TranslationAcc_X_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_M_Distrib_TranslationAcc_Y_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_M_Distrib_TranslationAcc_Z_'//TRIM(Num2LStr(TmpIndx))) 
   END DO   
   DO TmpIndx=1,u_HD%Morison%DistribMesh%NNodes
      CALL WrFileNR(UnJac, ' HD_M_Distrib_RotationAcc_X_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_M_Distrib_RotationAcc_Y_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_M_Distrib_RotationAcc_Z_'//TRIM(Num2LStr(TmpIndx))) 
   END DO     
       
   DO TmpIndx=1,u_HD%Mesh%NNodes
      CALL WrFileNR(UnJac, ' HD_Mesh_TranslationAcc_X_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_Mesh_TranslationAcc_Y_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_Mesh_TranslationAcc_Z_'//TRIM(Num2LStr(TmpIndx))) 
   END DO   
   DO TmpIndx=1,u_HD%Mesh%NNodes
      CALL WrFileNR(UnJac, ' HD_Mesh_RotationAcc_X_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_Mesh_RotationAcc_Y_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_Mesh_RotationAcc_Z_'//TRIM(Num2LStr(TmpIndx))) 
   END DO 
   
   DO nb=1,p_FAST%nBeams
      CALL WrFileNR(UnJac, ' BD_Root'//trim(num2lstr(nb))//'Motion_TranslationAcc_X') 
      CALL WrFileNR(UnJac, ' BD_Root'//trim(num2lstr(nb))//'Motion_TranslationAcc_Y') 
      CALL WrFileNR(UnJac, ' BD_Root'//trim(num2lstr(nb))//'Motion_TranslationAcc_Z') 
      CALL WrFileNR(UnJac, ' BD_Root'//trim(num2lstr(nb))//'Motion_RotationAcc_X') 
      CALL WrFileNR(UnJac, ' BD_Root'//trim(num2lstr(nb))//'Motion_RotationAcc_Y') 
      CALL WrFileNR(UnJac, ' BD_Root'//trim(num2lstr(nb))//'Motion_RotationAcc_Z') 
   END DO
         
   WRITE(UnJac,'()')    
      
   CALL WrMatrix(MeshMapData%Jacobian_ED_SD_HD_BD,UnJac, p_FAST%OutFmt)
   CLOSE( UnJac )

#endif               
            
            
               ! Get the LU decomposition of this matrix using a LAPACK routine: 
               ! The result is of the form MeshMapDat%Jacobian_ED_SD_HD_BD = P * L * U 

            CALL LAPACK_getrf( M=p_FAST%SizeJac_ED_SD_HD_BD(1), N=p_FAST%SizeJac_ED_SD_HD_BD(1), &
                              A=MeshMapData%Jacobian_ED_SD_HD_BD, IPIV=MeshMapData%Jacobian_pivot, &
                              ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN 
               END IF
            
         END IF         
            
         !-------------------------------------------------------------------------------------------------
         ! Solve for delta u: Jac*u_delta = - Fn_U_Resid
         !  using the LAPACK routine 
         !-------------------------------------------------------------------------------------------------
         
         u_delta = -Fn_U_Resid
         CALL LAPACK_getrs( TRANS="N", N=p_FAST%SizeJac_ED_SD_HD_BD(1), A=MeshMapData%Jacobian_ED_SD_HD_BD, &
                            IPIV=MeshMapData%Jacobian_pivot, B=u_delta, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN 
               END IF

         !-------------------------------------------------------------------------------------------------
         ! check for error, update inputs (u_ED and u_HD), and iterate again
         !-------------------------------------------------------------------------------------------------
                  
!         IF ( DOT_PRODUCT(u_delta, u_delta) <= TOL_Squared ) EXIT
         
         u = u + u_delta                  
         CALL Add_ED_SD_HD_BD_u_delta( p_FAST, MeshMapData%Jac_u_indx, u_delta, u_ED, u_SD, u_HD, u_BD )
                           
         K = K + 1
         
      END DO ! K
               
      !...............................................
      ! This is effectively doing option 2, where we set the input velocities and displacements based on the outputs we just calculated
      !...............................................
      
      ! BD motion inputs: (from ED)
      IF (p_FAST%CompElast == Module_BD ) THEN
         
            ! Make copies of the accelerations we just solved for (so we don't overwrite them)         
         do nb = 1,p_FAST%nBeams            
            MeshMapData%u_BD_RootMotion(nb)%RotationAcc     = u_BD(nb)%RootMotion%RotationAcc
            MeshMapData%u_BD_RootMotion(nb)%TranslationAcc  = u_BD(nb)%RootMotion%TranslationAcc                  
         end do
         
         call Transfer_ED_to_BD(y_ED, u_BD, MeshMapData, ErrStat2, ErrMsg2 )
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
                  
            ! put the acceleration data (calucluted in this routine) back         
         do nb = 1,p_FAST%nBeams            
            u_BD(nb)%RootMotion%RotationAcc    = MeshMapData%u_BD_RootMotion(nb)%RotationAcc   
            u_BD(nb)%RootMotion%TranslationAcc = MeshMapData%u_BD_RootMotion(nb)%TranslationAcc                  
         end do
         
      END IF
            
      
      !...............
      ! HD motion inputs: (from SD and ED)
      IF (p_FAST%CompHydro == Module_HD ) THEN
      
            ! Make copies of the accelerations we just solved for (so we don't overwrite them)         
         MeshMapData%u_HD_M_LumpedMesh%RotationAcc     = u_HD%Morison%LumpedMesh%RotationAcc
         MeshMapData%u_HD_M_LumpedMesh%TranslationAcc  = u_HD%Morison%LumpedMesh%TranslationAcc
         MeshMapData%u_HD_M_DistribMesh%RotationAcc    = u_HD%Morison%DistribMesh%RotationAcc
         MeshMapData%u_HD_M_DistribMesh%TranslationAcc = u_HD%Morison%DistribMesh%TranslationAcc
         MeshMapData%u_HD_Mesh%RotationAcc             = u_HD%Mesh%RotationAcc   
         MeshMapData%u_HD_Mesh%TranslationAcc          = u_HD%Mesh%TranslationAcc

            ! transfer the output data to inputs
            
         IF ( p_FAST%CompSub == Module_SD ) THEN
               ! Map SD outputs to HD inputs (keeping the accelerations we just calculated)
         
            CALL Transfer_SD_to_HD( y_SD, u_HD%Morison%LumpedMesh, u_HD%Morison%DistribMesh, MeshMapData, ErrStat2, ErrMsg2 )      
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
               
               ! Map ED outputs to HD inputs (keeping the accelerations we just calculated):
               
            CALL Transfer_Point_to_Point( y_ED%PlatformPtMesh, u_HD%Mesh, MeshMapData%ED_P_2_HD_W_P, ErrStat2, ErrMsg2 ) 
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
                                             
         ELSE
            
            CALL Transfer_ED_to_HD( y_ED, u_HD, MeshMapData, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
         
         END IF

         
            ! put the acceleration data (calucluted in this routine) back         
         u_HD%Morison%LumpedMesh%RotationAcc     = MeshMapData%u_HD_M_LumpedMesh%RotationAcc
         u_HD%Morison%LumpedMesh%TranslationAcc  = MeshMapData%u_HD_M_LumpedMesh%TranslationAcc  
         u_HD%Morison%DistribMesh%RotationAcc    = MeshMapData%u_HD_M_DistribMesh%RotationAcc    
         u_HD%Morison%DistribMesh%TranslationAcc = MeshMapData%u_HD_M_DistribMesh%TranslationAcc 
         u_HD%Mesh%RotationAcc                   = MeshMapData%u_HD_Mesh%RotationAcc    
         u_HD%Mesh%TranslationAcc                = MeshMapData%u_HD_Mesh%TranslationAcc 
         
         !......
                          
      END IF
      
      IF ( p_FAST%CompSub == Module_SD ) THEN       
         !...............
         ! SD motion inputs: (from ED)
                
            ! Map ED outputs to SD inputs (keeping the accelerations we just calculated):
      
         MeshMapData%u_SD_TPMesh%RotationAcc    = u_SD%TPMesh%RotationAcc   
         MeshMapData%u_SD_TPMesh%TranslationAcc = u_SD%TPMesh%TranslationAcc
         
         CALL Transfer_Point_to_Point( y_ED%PlatformPtMesh, u_SD%TPMesh, MeshMapData%ED_P_2_SD_TP, ErrStat2, ErrMsg2 ) 
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
                  
               
         u_SD%TPMesh%RotationAcc    = MeshMapData%u_SD_TPMesh%RotationAcc    
         u_SD%TPMesh%TranslationAcc = MeshMapData%u_SD_TPMesh%TranslationAcc    
      END IF
         
      !...............................................
      ! We're finished
      !...............................................
      CALL CleanUp()
      
CONTAINS                                 
   !...............................................................................................................................
   SUBROUTINE U_ED_SD_HD_BD_Residual( y_ED2, y_SD2, y_HD2, y_BD2, u_IN, U_Resid)
   ! transfer outputs of ED, HD, and SD (and any additional loads that get summed with them) into inputs for ED, HD, and SD
   !...............................................................................................................................
                                  
   TYPE(ED_OutputType)               , INTENT(IN   ) :: y_ED2                  ! System outputs
   TYPE(SD_OutputType)               , INTENT(IN   ) :: y_SD2                  ! System outputs
   TYPE(HydroDyn_OutputType)         , INTENT(IN   ) :: y_HD2                  ! System outputs
   TYPE(BD_OutputType)               , INTENT(IN   ) :: y_BD2(:)               ! System outputs
   REAL(ReKi)                        , INTENT(IN   ) :: u_in(:)
   REAL(ReKi)                        , INTENT(  OUT) :: U_Resid(:)

   INTEGER(IntKi)                                    :: i                      ! counter for ice leg and beamdyn loops
   
   !..................
   ! Set mooring line and ice inputs (which don't have acceleration fields and aren't used elsewhere in this routine, thus we're using the actual inputs (not a copy) 
   ! Note that these values get overwritten at the completion of this routine.)
   !..................
   
      IF ( p_FAST%CompMooring == Module_MAP ) THEN
         
         ! note: MAP_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)      
         CALL MAP_InputSolve( u_map, y_ED2, MeshMapData, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
                                 
      ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
         
         ! note: MD_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)      
         CALL MD_InputSolve( u_MD, y_ED2, MeshMapData, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)       
                        
      ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
         
         ! note: FEAM_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)      
         CALL FEAM_InputSolve( u_FEAM, y_ED2, MeshMapData, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)       
                        
      END IF     
   
      
      IF ( p_FAST%CompIce == Module_IceF ) THEN
         
         CALL IceFloe_InputSolve(  u_IceF, y_SD2, MeshMapData, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)       
                                 
      ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
         
         DO i=1,p_FAST%numIceLegs
            
            CALL IceD_InputSolve(  u_IceD(i), y_SD2, MeshMapData, i, ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)  
               
         END DO
         
      END IF  
      
      MeshMapData%u_ED_HubPtLoad%Force  = 0.0_ReKi
      MeshMapData%u_ED_HubPtLoad%Moment = 0.0_ReKi
      IF ( p_FAST%CompElast == Module_BD ) THEN
         
         ! Transfer ED motions to BD inputs:
         call Transfer_ED_to_BD_tmp( y_ED2, MeshMapData, ErrStat2, ErrMsg2 )   ! sets MeshMapData%u_BD_RootMotion(:)
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
         ! Transfer BD loads to ED hub input:
            ! we're mapping loads, so we also need the sibling meshes' displacements:
         do i=1,p_FAST%nBeams            
            CALL Transfer_Point_to_Point( y_BD2(i)%ReactionForce, MeshMapData%u_ED_HubPtLoad_2, MeshMapData%BD_P_2_ED_P(i), ErrStat2, ErrMsg2, MeshMapData%u_BD_RootMotion(i), y_ED2%HubPtMotion) !u_BD_RootMotion and y_ED2%HubPtMotion contain the displaced positions for load calculations
               CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
               
            MeshMapData%u_ED_HubPtLoad%Force  = MeshMapData%u_ED_HubPtLoad%Force  + MeshMapData%u_ED_HubPtLoad_2%Force  
            MeshMapData%u_ED_HubPtLoad%Moment = MeshMapData%u_ED_HubPtLoad%Moment + MeshMapData%u_ED_HubPtLoad_2%Moment                
         end do
      END IF
      
      
      
      IF ( p_FAST%CompSub == Module_SD ) THEN
      
         IF ( p_FAST%CompHydro == Module_HD ) THEN
            
         ! initialize these SD loads inputs here in case HD is used  (note from initialiazation that these meshes don't exist if HD isn't used)       
         MeshMapData%u_SD_LMesh%Force  = 0.0_ReKi
         MeshMapData%u_SD_LMesh%Moment = 0.0_ReKi
      
            
      !..................
      ! Get HD inputs on Morison%LumpedMesh and Morison%DistribMesh
      !..................
   
               ! SD motions to HD:
            CALL Transfer_SD_to_HD( y_SD2, MeshMapData%u_HD_M_LumpedMesh, MeshMapData%u_HD_M_DistribMesh, MeshMapData, ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)       
      
   
               ! Map ED motion output to HD inputs:
            CALL Transfer_Point_to_Point( y_ED2%PlatformPtMesh, MeshMapData%u_HD_Mesh, MeshMapData%ED_P_2_HD_W_P, ErrStat2, ErrMsg2 ) 
               CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)      
               
      !..................
      ! Get SD loads inputs (MeshMapData%u_HD_M_LumpedMesh and MeshMapData%u_HD_M_DistribMesh meshes must be set first)
      !..................
         
            ! Loads (outputs) from HD meshes transfered to SD LMesh (zero them out first because they get summed in Transfer_HD_to_SD)
         
            CALL Transfer_HD_to_SD( MeshMapData%u_SD_LMesh_2, MeshMapData%u_SD_LMesh, y_SD2%Y2Mesh, y_HD2, MeshMapData%u_HD_M_LumpedMesh, MeshMapData%u_HD_M_DistribMesh, MeshMapData, ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)               
         

            IF ( p_FAST%CompIce == Module_IceF ) THEN
               
               ! SD loads from IceFloe:
               IF ( y_IceF%iceMesh%Committed ) THEN      
                  ! we're mapping loads, so we also need the sibling meshes' displacements:
                  CALL Transfer_Point_to_Point( y_IceF%iceMesh, MeshMapData%u_SD_LMesh_2, MeshMapData%IceF_P_2_SD_P, ErrStat2, ErrMsg2, u_IceF%iceMesh, y_SD2%Y2Mesh )   
                     CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)

                  MeshMapData%u_SD_LMesh%Force  = MeshMapData%u_SD_LMesh%Force  + MeshMapData%u_SD_LMesh_2%Force
                  MeshMapData%u_SD_LMesh%Moment = MeshMapData%u_SD_LMesh%Moment + MeshMapData%u_SD_LMesh_2%Moment    
                  
!...          
#ifdef DEBUG_MESH_TRANSFER_ICE
   if (.not. calcJacobian) then
         CALL WrScr('********************************************************')
         CALL WrScr('****   IceF to SD point-to-point                   *****')
         CALL WrScr('********************************************************')
         CALL WriteMappingTransferToFile(MeshMapData%u_SD_LMesh_2, y_SD2%Y2Mesh, u_IceF%iceMesh, y_IceF%iceMesh,&
               MeshMapData%SD_P_2_IceF_P, MeshMapData%IceF_P_2_SD_P, &
               'SD_y2_IceF_Meshes_t'//TRIM(Num2LStr(this_time))//'.I.bin' )
         !print *
         !pause         
   end IF         
#endif                
                                    
               END IF !Module_IceF
               
            ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
               
                ! SD loads from IceDyn:
               DO i=1,p_FAST%numIceLegs
                  
                  IF ( y_IceD(i)%PointMesh%Committed ) THEN      
                     ! we're mapping loads, so we also need the sibling meshes' displacements:
                     CALL Transfer_Point_to_Point( y_IceD(i)%PointMesh, MeshMapData%u_SD_LMesh_2, MeshMapData%IceD_P_2_SD_P(i), ErrStat2, ErrMsg2, u_IceD(i)%PointMesh, y_SD2%Y2Mesh )   
                        CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)

                     MeshMapData%u_SD_LMesh%Force  = MeshMapData%u_SD_LMesh%Force  + MeshMapData%u_SD_LMesh_2%Force
                     MeshMapData%u_SD_LMesh%Moment = MeshMapData%u_SD_LMesh%Moment + MeshMapData%u_SD_LMesh_2%Moment    
                     
                  END IF
                  
               END DO
               
            END IF   ! Ice loading
               
         END IF ! HD is used (IceFloe/IceDyn can't be used unless HydroDyn is used)


         
      !..................
      ! Get SD motions input
      !..................
         
         ! Motions (outputs) at ED platform ref point transfered to SD transition piece (input):
         CALL Transfer_Point_to_Point( y_ED2%PlatformPtMesh, MeshMapData%u_SD_TPMesh, MeshMapData%ED_P_2_SD_TP, ErrStat2, ErrMsg2 )   
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
      !..................
      ! Get ED loads input (from SD and possibly HD)
      !..................
            
         ! Loads (outputs) on the SD transition piece transfered to ED input location/mesh:
            ! we're mapping loads, so we also need the sibling meshes' displacements:
         CALL Transfer_Point_to_Point( y_SD2%Y1Mesh, MeshMapData%u_ED_PlatformPtMesh, MeshMapData%SD_TP_2_ED_P, ErrStat2, ErrMsg2, MeshMapData%u_SD_TPMesh, y_ED2%PlatformPtMesh ) !MeshMapData%u_SD_TPMesh contains the orientations needed for moment calculations
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
                        
               ! WAMIT loads from HD get added to this load:
         IF ( y_HD2%Mesh%Committed  ) THEN
   
            ! we're mapping loads, so we also need the sibling meshes' displacements:
            CALL Transfer_Point_to_Point( y_HD2%Mesh, MeshMapData%u_ED_PlatformPtMesh_2, MeshMapData%HD_W_P_2_ED_P, ErrStat2, ErrMsg2, MeshMapData%u_HD_Mesh, y_ED2%PlatformPtMesh ) !u_SD contains the orientations needed for moment calculations
               CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
         
            MeshMapData%u_ED_PlatformPtMesh%Force  = MeshMapData%u_ED_PlatformPtMesh%Force  + MeshMapData%u_ED_PlatformPtMesh_2%Force
            MeshMapData%u_ED_PlatformPtMesh%Moment = MeshMapData%u_ED_PlatformPtMesh%Moment + MeshMapData%u_ED_PlatformPtMesh_2%Moment
   
         END IF              
            
      ELSEIF ( p_FAST%CompHydro == Module_HD ) THEN 

      !..................
      ! Get HD inputs on 3 meshes
      !..................
         
         ! Map ED motion outputs to HD inputs:
         ! basically, we want to call Transfer_ED_to_HD, except we have the meshes in a different data structure (not a copy of u_HD)
         ! CALL Transfer_ED_to_HD( y_ED2, u_HD, MeshMapData, ErrStat2, ErrMsg2 ) 
         ! so, here are the transfers, again.
         
         
         ! These are the motions for the lumped point loads associated the WAMIT body:
         CALL Transfer_Point_to_Point( y_ED2%PlatformPtMesh, MeshMapData%u_HD_Mesh, MeshMapData%ED_P_2_HD_W_P, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
   
         ! These are the motions for the lumped point loads associated viscous drag on the WAMIT body and/or filled/flooded lumped forces of the WAMIT body
         CALL Transfer_Point_to_Point( y_ED2%PlatformPtMesh, MeshMapData%u_HD_M_LumpedMesh, MeshMapData%ED_P_2_HD_M_P, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
                  
         ! These are the motions for the line2 (distributed) loads associated viscous drag on the WAMIT body and/or filled/flooded distributed forces of the WAMIT body
         CALL Transfer_Point_to_Line2( y_ED2%PlatformPtMesh, MeshMapData%u_HD_M_DistribMesh, MeshMapData%ED_P_2_HD_M_L, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
      !..................
      ! Get ED loads input (from HD only)
      !..................
            
            ! we're mapping loads, so we also need the sibling meshes' displacements:
         CALL Transfer_Point_to_Point( y_HD2%AllHdroOrigin, MeshMapData%u_ED_PlatformPtMesh, MeshMapData%HD_W_P_2_ED_P, ErrStat2, ErrMsg2, MeshMapData%u_HD_Mesh, y_ED2%PlatformPtMesh) !u_HD and u_mapped_positions contain the displaced positions for load calculations
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
                                                                           
      END IF
      
   !..................
   ! Get remaining portion of ED loads input on MeshMapData%u_ED_PlatformPtMesh (must do this after MeshMapData%u_SD_TPMesh and MeshMapData%u_HD_Mesh are set)
   !   at this point, MeshMapData%u_ED_PlatformPtMesh contains the portion of loads from SD and/or HD
   !..................
                     
         
         ! Get the loads for ED from a mooring module and add them:
      IF ( p_FAST%CompMooring == Module_MAP ) THEN
         CALL Transfer_Point_to_Point( y_MAP%PtFairleadLoad, MeshMapData%u_ED_PlatformPtMesh_2, MeshMapData%Mooring_P_2_ED_P, ErrStat2, ErrMsg2, u_MAP%PtFairDisplacement, y_ED2%PlatformPtMesh ) !u_MAP and y_ED contain the displacements needed for moment calculations
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)               
            
         MeshMapData%u_ED_PlatformPtMesh%Force  = MeshMapData%u_ED_PlatformPtMesh%Force  + MeshMapData%u_ED_PlatformPtMesh_2%Force
         MeshMapData%u_ED_PlatformPtMesh%Moment = MeshMapData%u_ED_PlatformPtMesh%Moment + MeshMapData%u_ED_PlatformPtMesh_2%Moment
            
      ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
         CALL Transfer_Point_to_Point( y_MD%PtFairleadLoad, MeshMapData%u_ED_PlatformPtMesh_2, MeshMapData%Mooring_P_2_ED_P, ErrStat2, ErrMsg2, u_MD%PtFairleadDisplacement, y_ED2%PlatformPtMesh ) !u_MD and y_ED contain the displacements needed for moment calculations
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)  
            
         MeshMapData%u_ED_PlatformPtMesh%Force  = MeshMapData%u_ED_PlatformPtMesh%Force  + MeshMapData%u_ED_PlatformPtMesh_2%Force
         MeshMapData%u_ED_PlatformPtMesh%Moment = MeshMapData%u_ED_PlatformPtMesh%Moment + MeshMapData%u_ED_PlatformPtMesh_2%Moment            

      ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
         CALL Transfer_Point_to_Point( y_FEAM%PtFairleadLoad, MeshMapData%u_ED_PlatformPtMesh_2, MeshMapData%Mooring_P_2_ED_P, ErrStat2, ErrMsg2, u_FEAM%PtFairleadDisplacement, y_ED2%PlatformPtMesh ) !u_FEAM and y_ED contain the displacements needed for moment calculations
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)  
            
         MeshMapData%u_ED_PlatformPtMesh%Force  = MeshMapData%u_ED_PlatformPtMesh%Force  + MeshMapData%u_ED_PlatformPtMesh_2%Force
         MeshMapData%u_ED_PlatformPtMesh%Moment = MeshMapData%u_ED_PlatformPtMesh%Moment + MeshMapData%u_ED_PlatformPtMesh_2%Moment            
      END IF
                                                   
   !..................
   ! Calculate the residual with these new inputs:
   !..................                  
         
      CALL Create_ED_SD_HD_BD_UVector(U_Resid, MeshMapData%u_ED_PlatformPtMesh, MeshMapData%u_SD_TPMesh, MeshMapData%u_SD_LMesh, &
                                      MeshMapData%u_HD_M_LumpedMesh,   MeshMapData%u_HD_M_DistribMesh, MeshMapData%u_HD_Mesh, &
                                      MeshMapData%u_ED_HubPtLoad, MeshMapData%u_BD_RootMotion, p_FAST ) 
         
      U_Resid = u_in - U_Resid
   
            
   END SUBROUTINE U_ED_SD_HD_BD_Residual   
   !...............................................................................................................................
   SUBROUTINE CleanUp()
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(1024)            :: ErrMsg3     ! The error message (ErrMsg)
         
      IF ( calcJacobian ) THEN
         CALL ED_DestroyInput( u_ED_perturb, ErrStat3, ErrMsg3 )
            IF (ErrStat3 /= ErrID_None) CALL WrScr(' '//RoutineName//TRIM(ErrMsg3) )
         CALL ED_DestroyOutput(y_ED_perturb, ErrStat3, ErrMsg3 )
            IF (ErrStat3 /= ErrID_None) CALL WrScr(' '//RoutineName//TRIM(ErrMsg3) )
         
         CALL SD_DestroyInput( u_SD_perturb, ErrStat3, ErrMsg3 )
            IF (ErrStat3 /= ErrID_None) CALL WrScr(' '//RoutineName//TRIM(ErrMsg3) )
         CALL SD_DestroyOutput(y_SD_perturb, ErrStat3, ErrMsg3 )
            IF (ErrStat3 /= ErrID_None) CALL WrScr(' '//RoutineName//TRIM(ErrMsg3) )

         CALL HydroDyn_DestroyInput( u_HD_perturb, ErrStat3, ErrMsg3 )
            IF (ErrStat3 /= ErrID_None) CALL WrScr(' '//RoutineName//TRIM(ErrMsg3) )
         CALL HydroDyn_DestroyOutput(y_HD_perturb, ErrStat3, ErrMsg3 )
            IF (ErrStat3 /= ErrID_None) CALL WrScr(' '//RoutineName//TRIM(ErrMsg3) )
            
         CALL BD_DestroyInput( u_BD_perturb, ErrStat3, ErrMsg3 )
            IF (ErrStat3 /= ErrID_None) CALL WrScr(' '//RoutineName//TRIM(ErrMsg3) )
         if (allocated(y_BD_perturb)) then
            do nb=1,size(y_BD_perturb) 
               CALL BD_DestroyOutput(y_BD_perturb(nb), ErrStat3, ErrMsg3 )
                  IF (ErrStat3 /= ErrID_None) CALL WrScr(' '//RoutineName//TRIM(ErrMsg3) )
            end do
            deallocate(y_BD_perturb)
         end if
                        
      END IF
      
   
   END SUBROUTINE CleanUp
   !...............................................................................................................................
END SUBROUTINE ED_SD_HD_BD_InputOutputSolve
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Init_ED_SD_HD_BD_Jacobian( p_FAST, MeshMapData, ED_PlatformPtMesh, SD_TPMesh, SD_LMesh, HD_M_LumpedMesh, HD_M_DistribMesh, &
                                   HD_WAMIT_Mesh, ED_HubPtLoad, u_BD, ErrStat, ErrMsg)

   TYPE(FAST_ParameterType)          , INTENT(INOUT) :: p_FAST                     
   TYPE(FAST_ModuleMapType)          , INTENT(INOUT) :: MeshMapData
   
      ! input meshes for each of the 4 modules:
   TYPE(MeshType)                    , INTENT(IN   ) :: ED_PlatformPtMesh
   TYPE(MeshType)                    , INTENT(IN   ) :: ED_HubPtLoad
   TYPE(MeshType)                    , INTENT(IN   ) :: SD_TPMesh
   TYPE(MeshType)                    , INTENT(IN   ) :: SD_LMesh
   TYPE(MeshType)                    , INTENT(IN   ) :: HD_M_LumpedMesh
   TYPE(MeshType)                    , INTENT(IN   ) :: HD_M_DistribMesh
   TYPE(MeshType)                    , INTENT(IN   ) :: HD_WAMIT_Mesh
   TYPE(BD_InputType)                , INTENT(IN   ) :: u_BD(:)
   
   INTEGER(IntKi)                    , INTENT(  OUT) :: ErrStat                   ! Error status of the operation
   CHARACTER(*)                      , INTENT(  OUT) :: ErrMsg                    ! Error message if ErrStat /= ErrID_None
   
   CHARACTER(*), PARAMETER                           :: RoutineName = 'Init_ED_SD_HD_BD_Jacobian'
   
      ! local variables:
   INTEGER(IntKi)                :: i, j, k, index
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
      ! determine how many inputs there are between the 4 modules (ED, SD, HD, BD)
   
   if (p_FAST%CompHydro == Module_HD .or. p_FAST%CompSub == Module_SD) then
      p_FAST%SizeJac_ED_SD_HD_BD(2) = ED_PlatformPtMesh%NNodes*6        ! ED inputs: 3 forces and 3 moments per node (only 1 node)
   else
      p_FAST%SizeJac_ED_SD_HD_BD(2) = 0
   end if
   
                  
   p_FAST%SizeJac_ED_SD_HD_BD(3) = SD_TPMesh%NNodes*6                ! SD inputs: 6 accelerations per node (size of SD input from ED) 
   IF ( p_FAST%CompHydro == Module_HD ) THEN   
      p_FAST%SizeJac_ED_SD_HD_BD(3) = p_FAST%SizeJac_ED_SD_HD_BD(3) &   
                                    + SD_LMesh%NNodes *6             ! SD inputs: 6 loads per node (size of SD input from HD)       
   END IF
               
   p_FAST%SizeJac_ED_SD_HD_BD(4) = HD_M_LumpedMesh%NNodes *6 &    ! HD inputs: 6 accelerations per node (on each Morison mesh) 
                                 + HD_M_DistribMesh%NNodes*6 &    ! HD inputs: 6 accelerations per node (on each Morison mesh)
                                 + HD_WAMIT_Mesh%NNodes*6         ! HD inputs: 6 accelerations per node (on the WAMIT mesh)      
   
   IF ( p_FAST%CompElast == Module_BD ) THEN   
      p_FAST%SizeJac_ED_SD_HD_BD(2) = p_FAST%SizeJac_ED_SD_HD_BD(2) &   
                                     + ED_HubPtLoad%NNodes *6     ! ED inputs: 6 loads per node (size of ED input from BD)
      
      p_FAST%SizeJac_ED_SD_HD_BD(5:7) = 0 ! assumes a max of 3 blades
      do k=1,size(u_BD)
         p_FAST%SizeJac_ED_SD_HD_BD(4+k) = u_BD(k)%RootMotion%NNodes *6     ! BD inputs: 6 accelerations per node (size of BD input from ED)         
      end do
            
   END IF
                              
                              
                              
   p_FAST%SizeJac_ED_SD_HD_BD(1) = sum( p_FAST%SizeJac_ED_SD_HD_BD )   ! all the inputs from these modules
                  

      ! allocate matrix to store jacobian 
   CALL AllocAry( MeshMapData%Jacobian_ED_SD_HD_BD, p_FAST%SizeJac_ED_SD_HD_BD(1), p_FAST%SizeJac_ED_SD_HD_BD(1), "Jacobian for ED-SD-HD-BD", ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
         
      ! allocate matrix to store index to help us figure out what the ith value of the u vector really means
   ALLOCATE ( MeshMapData%Jac_u_indx( p_FAST%SizeJac_ED_SD_HD_BD(1), 3 ), STAT = ErrStat )
      IF ( ErrStat /= 0 ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg = 'Cannot allocate Jac_u_indx.'
         RETURN
      END IF
         
   ! fill matrix to store index to help us figure out what the ith value of the u vector really means
   ! ( see Create_ED_SD_HD_BD_UVector() ... these MUST match )
   ! column 1 indicates module's mesh and field
   ! column 2 indicates the first index of the acceleration/load field
   ! column 3 is the node
      
   !...............
   ! ED inputs:   
   !...............
   
   index = 1
   if (p_FAST%CompHydro == Module_HD .or. p_FAST%CompSub == Module_SD) then
   
      do i=1,ED_PlatformPtMesh%NNodes
         do j=1,3
            MeshMapData%Jac_u_indx(index,1) =  1 !Module/Mesh/Field: u_ED%PlatformPtMesh%Force = 1
            MeshMapData%Jac_u_indx(index,2) =  j !index:  j
            MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
            index = index + 1
         end do !j      
      end do !i
   
      do i=1,ED_PlatformPtMesh%NNodes
         do j=1,3
            MeshMapData%Jac_u_indx(index,1) =  2 !Module/Mesh/Field: u_ED%PlatformPtMesh%Moment = 2
            MeshMapData%Jac_u_indx(index,2) =  j !index:  j
            MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
            index = index + 1
         end do !j      
      end do !i
      
   end if
   
   
   if (p_FAST%CompElast == Module_BD) then
      
      do i=1,ED_HubPtLoad%NNodes
         do j=1,3
            MeshMapData%Jac_u_indx(index,1) =  3 !Module/Mesh/Field: u_ED%HubPtMesh%Force = 3
            MeshMapData%Jac_u_indx(index,2) =  j !index:  j
            MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
            index = index + 1
         end do !j      
      end do !i
      
      
      do i=1,ED_HubPtLoad%NNodes
         do j=1,3
            MeshMapData%Jac_u_indx(index,1) =  4 !Module/Mesh/Field: u_ED%HubPtMesh%Moment = 4
            MeshMapData%Jac_u_indx(index,2) =  j !index:  j
            MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
            index = index + 1
         end do !j      
      end do !i

   end if
   
      
   !...............
   ! SD inputs:   
   !...............
      
   ! SD_TPMesh                        
   do i=1,SD_TPMesh%NNodes
      do j=1,3
         MeshMapData%Jac_u_indx(index,1) =  5 !Module/Mesh/Field: u_SD%TPMesh%TranslationAcc = 5
         MeshMapData%Jac_u_indx(index,2) =  j !index:  j
         MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
         index = index + 1                  
      end do !j                             
   end do !i                                
                                            
   do i=1,SD_TPMesh%NNodes                  
      do j=1,3                              
         MeshMapData%Jac_u_indx(index,1) =  6 !Module/Mesh/Field:  u_SD%TPMesh%RotationAcc = 6
         MeshMapData%Jac_u_indx(index,2) =  j !index:  j
         MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
         index = index + 1
      end do !j      
   end do !i   
   
   IF ( p_FAST%CompHydro == Module_HD ) THEN   ! this SD mesh linked only when HD is enabled
   
      ! SD_LMesh
      do i=1,SD_LMesh%NNodes
         do j=1,3
            MeshMapData%Jac_u_indx(index,1) =  7 !Module/Mesh/Field: u_SD%LMesh%Force = 7
            MeshMapData%Jac_u_indx(index,2) =  j !index:  j
            MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
            index = index + 1                  
         end do !j                             
      end do !i                                
                                            
      do i=1,SD_LMesh%NNodes                   
         do j=1,3                              
            MeshMapData%Jac_u_indx(index,1) =  8 !Module/Mesh/Field: u_SD%LMesh%Moment = 8
            MeshMapData%Jac_u_indx(index,2) =  j !index:  j
            MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
            index = index + 1
         end do !j      
      end do !i 
      
   END IF
   
   !...............
   ! HD inputs:
   !...............
         
   !(Morison%LumpedMesh)
   do i=1,HD_M_LumpedMesh%NNodes
      do j=1,3
         MeshMapData%Jac_u_indx(index,1) =  9 !Module/Mesh/Field: u_HD%Morison%LumpedMesh%TranslationAcc = 9
         MeshMapData%Jac_u_indx(index,2) =  j !index:  j
         MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
         index = index + 1
      end do !j      
   end do !i
   
   do i=1,HD_M_LumpedMesh%NNodes
      do j=1,3
         MeshMapData%Jac_u_indx(index,1) = 10 !Module/Mesh/Field:  u_HD%Morison%LumpedMesh%RotationAcc = 10
         MeshMapData%Jac_u_indx(index,2) =  j !index:  j
         MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
         index = index + 1
      end do !j      
   end do !i     
   
   
   !(Morison%DistribMesh)
   do i=1,HD_M_DistribMesh%NNodes
      do j=1,3
         MeshMapData%Jac_u_indx(index,1) = 11 !Module/Mesh/Field: u_HD%Morison%DistribMesh%TranslationAcc = 11
         MeshMapData%Jac_u_indx(index,2) =  j !index:  j
         MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
         index = index + 1
      end do !j      
   end do !i
   
   do i=1,HD_M_DistribMesh%NNodes
      do j=1,3
         MeshMapData%Jac_u_indx(index,1) = 12 !Module/Mesh/Field:  u_HD%Morison%DistribMesh%RotationAcc = 12
         MeshMapData%Jac_u_indx(index,2) =  j !index:  j
         MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
         index = index + 1
      end do !j      
   end do !i     
   
   !(Mesh)
   do i=1,HD_WAMIT_Mesh%NNodes
      do j=1,3
         MeshMapData%Jac_u_indx(index,1) = 13 !Module/Mesh/Field: u_HD%Mesh%TranslationAcc = 13
         MeshMapData%Jac_u_indx(index,2) =  j !index:  j
         MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
         index = index + 1
      end do !j      
   end do !i
   
   do i=1,HD_WAMIT_Mesh%NNodes
      do j=1,3
         MeshMapData%Jac_u_indx(index,1) = 14 !Module/Mesh/Field:  u_HD%Mesh%RotationAcc = 14
         MeshMapData%Jac_u_indx(index,2) =  j !index:  j
         MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
         index = index + 1
      end do !j      
   end do !i        
   
   !...............
   ! BD inputs:
   !...............
   
   if (p_FAST%CompElast == Module_BD) then
                 
      do k=1,size(u_BD)
         
         do i=1,u_BD(k)%RootMotion%NNodes
            do j=1,3
               MeshMapData%Jac_u_indx(index,1) =  13 + 2*k !Module/Mesh/Field: u_BD(k)%RootMotion%TranslationAcc = 15 (k=1), 17 (k=2), 19 (k=2)
               MeshMapData%Jac_u_indx(index,2) =  j !index:  j
               MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
               index = index + 1
            end do !j      
         end do !i
      
         do i=1,u_BD(k)%RootMotion%NNodes
            do j=1,3
               MeshMapData%Jac_u_indx(index,1) =  14 + 2*k !Module/Mesh/Field: u_BD(k)%RootMotion%RotationAcc = 16 (k=1), 18 (k=2), 20 (k=2)
               MeshMapData%Jac_u_indx(index,2) =  j !index:  j
               MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
               index = index + 1
            end do !j      
         end do !i
                  
      end do !k
                  
   end if
   
   
END SUBROUTINE Init_ED_SD_HD_BD_Jacobian
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Create_ED_SD_HD_BD_UVector(u, ED_PlatformPtMesh, SD_TPMesh, SD_LMesh, HD_M_LumpedMesh, HD_M_DistribMesh, HD_WAMIT_Mesh, &
                                         ED_HubPtLoad,  BD_RootMotion, p_FAST )
! This routine basically packs the relevant parts of the modules' input meshes for use in this InputOutput solve.
! Do not change the order of this packing without changing subroutine Init_ED_SD_HD_Jacobian()!
!..................................................................................................................................
   
   REAL(ReKi)                        , INTENT(INOUT) :: u(:)                     ! output u vector
   
      ! input meshes for each of the 3 modules:
   TYPE(MeshType)                    , INTENT(IN   ) :: ED_PlatformPtMesh
   TYPE(MeshType)                    , INTENT(IN   ) :: SD_TPMesh
   TYPE(MeshType)                    , INTENT(IN   ) :: SD_LMesh
   TYPE(MeshType)                    , INTENT(IN   ) :: HD_M_LumpedMesh
   TYPE(MeshType)                    , INTENT(IN   ) :: HD_M_DistribMesh
   TYPE(MeshType)                    , INTENT(IN   ) :: HD_WAMIT_Mesh
   TYPE(MeshType)                    , INTENT(IN   ) :: ED_HubPtLoad
   TYPE(MeshType)                    , INTENT(IN   ) :: BD_RootMotion(:)
   
   TYPE(FAST_ParameterType)          , INTENT(IN   ) :: p_FAST
   
   
      ! local variables:
   INTEGER(IntKi)                :: i, k, indx_first, indx_last
   
   !...............
   ! ED inputs:   
   !...............
   if (p_FAST%CompHydro == Module_HD .or. p_FAST%CompSub == Module_SD) then
      u( 1: 3) = ED_PlatformPtMesh%Force(:,1) / p_FAST%UJacSclFact
      u( 4: 6) = ED_PlatformPtMesh%Moment(:,1) / p_FAST%UJacSclFact    
      indx_first = 7   
   else
      indx_first = 1   
   end if
   
   
   if (p_FAST%CompElast == Module_BD) then
      
      do i=1,ED_HubPtLoad%NNodes
         indx_last  = indx_first + 2 
         u(indx_first:indx_last) = ED_HubPtLoad%Force(:,i) / p_FAST%UJacSclFact
         indx_first = indx_last + 1
      end do
     
      do i=1,ED_HubPtLoad%NNodes
         indx_last  = indx_first + 2 
         u(indx_first:indx_last) = ED_HubPtLoad%Moment(:,i) / p_FAST%UJacSclFact
         indx_first = indx_last + 1
      end do
      
   end if
   

   !...............
   ! SD inputs (SD_TPMesh):      
   !...............
   do i=1,SD_TPMesh%NNodes 
      indx_last  = indx_first + 2
      u(indx_first:indx_last) = SD_TPMesh%TranslationAcc(:,i) 
      indx_first = indx_last + 1
   end do

   do i=1,SD_TPMesh%NNodes 
      indx_last  = indx_first + 2
      u(indx_first:indx_last) = SD_TPMesh%RotationAcc(:,i) 
      indx_first = indx_last + 1
   end do
         
   if ( p_FAST%CompHydro == Module_HD ) then   ! this SD mesh linked only when HD is enabled
      ! SD inputs (SD_LMesh):        
      do i=1,SD_LMesh%NNodes
         indx_last  = indx_first + 2 
         u(indx_first:indx_last) = SD_LMesh%Force(:,i) / p_FAST%UJacSclFact
         indx_first = indx_last + 1
      end do
     
      do i=1,SD_LMesh%NNodes
         indx_last  = indx_first + 2 
         u(indx_first:indx_last) = SD_LMesh%Moment(:,i) / p_FAST%UJacSclFact
         indx_first = indx_last + 1
      end do
   end if
   
   !...............
   ! HD inputs (Morison%LumpedMesh):
   !...............
   do i=1,HD_M_LumpedMesh%NNodes
      indx_last  = indx_first + 2 
      u(indx_first:indx_last) = HD_M_LumpedMesh%TranslationAcc(:,i)
      indx_first = indx_last + 1
   end do
      
   do i=1,HD_M_LumpedMesh%NNodes
      indx_last  = indx_first + 2 
      u(indx_first:indx_last) = HD_M_LumpedMesh%RotationAcc(:,i)
      indx_first = indx_last + 1
   end do
      
   ! HD inputs (Morison%DistribMesh):
   do i=1,HD_M_DistribMesh%NNodes
      indx_last  = indx_first + 2 
      u(indx_first:indx_last) = HD_M_DistribMesh%TranslationAcc(:,i)
      indx_first = indx_last + 1
   end do
      
   do i=1,HD_M_DistribMesh%NNodes
      indx_last  = indx_first + 2 
      u(indx_first:indx_last) = HD_M_DistribMesh%RotationAcc(:,i)
      indx_first = indx_last + 1
   end do
                  
   ! HD inputs (Mesh):
   do i=1,HD_WAMIT_Mesh%NNodes
      indx_last  = indx_first + 2 
      u(indx_first:indx_last) = HD_WAMIT_Mesh%TranslationAcc(:,i)
      indx_first = indx_last + 1
   end do
      
   do i=1,HD_WAMIT_Mesh%NNodes
      indx_last  = indx_first + 2 
      u(indx_first:indx_last) = HD_WAMIT_Mesh%RotationAcc(:,i)
      indx_first = indx_last + 1
   end do   
   
   !...............
   ! BD inputs:
   !...............   
   if (p_FAST%CompElast == Module_BD) then
      
      do k=1,p_FAST%nBeams
         
         do i=1,BD_RootMotion(k)%NNodes
            indx_last  = indx_first + 2 
            u(indx_first:indx_last) = BD_RootMotion(k)%TranslationAcc(:,i)
            indx_first = indx_last + 1
         end do !i
      
         do i=1,BD_RootMotion(k)%NNodes
            indx_last  = indx_first + 2 
            u(indx_first:indx_last) = BD_RootMotion(k)%RotationAcc(:,i)
            indx_first = indx_last + 1
         end do !i
                  
      end do !k      
      
   end if
   
END SUBROUTINE Create_ED_SD_HD_BD_UVector
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Add_ED_SD_HD_BD_u_delta( p_FAST, Jac_u_indx, u_delta, u_ED, u_SD, u_HD, u_BD )
! This routine adds u_delta to the corresponding mesh field and scales it as appropriate
!..................................................................................................................................
   TYPE(FAST_ParameterType)            , INTENT(IN   ) :: p_FAST         ! Glue-code simulation parameters
   INTEGER( IntKi )                    , INTENT(IN   ) :: Jac_u_indx(:,:)
   REAL( ReKi )                        , INTENT(IN   ) :: u_delta(:)
   TYPE(ED_InputType)                  , INTENT(INOUT) :: u_ED           ! ED System inputs 
   TYPE(SD_InputType)                  , INTENT(INOUT) :: u_SD           ! SD System inputs 
   TYPE(HydroDyn_InputType)            , INTENT(INOUT) :: u_HD           ! SD System inputs 
   TYPE(BD_InputType)                  , INTENT(INOUT) :: u_BD(:)        ! BD System inputs 
   
   ! local variables
   INTEGER                                             :: n
   INTEGER                                             :: fieldIndx
   INTEGER                                             :: node
   
   
   DO n = 1,SIZE(u_delta)
      
      fieldIndx = Jac_u_indx(n,2) 
      node      = Jac_u_indx(n,3) 
   
         ! determine which mesh we're trying to perturb and perturb the input:
      SELECT CASE( Jac_u_indx(n,1) )
      
      CASE ( 1) !Module/Mesh/Field: u_ED%PlatformPtMesh%Force = 1
         u_ED%PlatformPtMesh%Force( fieldIndx,node) = u_ED%PlatformPtMesh%Force( fieldIndx,node) + u_delta(n) * p_FAST%UJacSclFact       
      CASE ( 2) !Module/Mesh/Field: u_ED%PlatformPtMesh%Moment = 2
         u_ED%PlatformPtMesh%Moment(fieldIndx,node) = u_ED%PlatformPtMesh%Moment(fieldIndx,node) + u_delta(n) * p_FAST%UJacSclFact       
      CASE ( 3) !Module/Mesh/Field: u_ED%HubPtMesh%Force = 3
         u_ED%HubPtLoad%Force( fieldIndx,node) = u_ED%HubPtLoad%Force( fieldIndx,node) + u_delta(n) * p_FAST%UJacSclFact     
      CASE ( 4) !Module/Mesh/Field: u_ED%HubPtMesh%Moment = 4
         u_ED%HubPtLoad%Moment(fieldIndx,node) = u_ED%HubPtLoad%Moment(fieldIndx,node) + u_delta(n) * p_FAST%UJacSclFact     
               
      CASE ( 5) !Module/Mesh/Field: u_SD%TPMesh%TranslationAcc = 5
         u_SD%TPMesh%TranslationAcc(fieldIndx,node) = u_SD%TPMesh%TranslationAcc(fieldIndx,node) + u_delta(n)        
      CASE ( 6) !Module/Mesh/Field: u_SD%TPMesh%RotationAcc = 6
         u_SD%TPMesh%RotationAcc(   fieldIndx,node) = u_SD%TPMesh%RotationAcc(   fieldIndx,node) + u_delta(n)        
      CASE ( 7) !Module/Mesh/Field: u_SD%LMesh%Force = 7
         u_SD%LMesh%Force( fieldIndx,node) = u_SD%LMesh%Force( fieldIndx,node) + u_delta(n) * p_FAST%UJacSclFact     
      CASE ( 8) !Module/Mesh/Field: u_SD%LMesh%Moment = 8
         u_SD%LMesh%Moment(fieldIndx,node) = u_SD%LMesh%Moment(fieldIndx,node) + u_delta(n) * p_FAST%UJacSclFact     
      
      CASE ( 9) !Module/Mesh/Field: u_HD%Morison%LumpedMesh%TranslationAcc = 9
         u_HD%Morison%LumpedMesh%TranslationAcc(fieldIndx,node) = u_HD%Morison%LumpedMesh%TranslationAcc(fieldIndx,node) + u_delta(n)        
      CASE (10) !Module/Mesh/Field: u_HD%Morison%LumpedMesh%RotationAcc = 10
         u_HD%Morison%LumpedMesh%RotationAcc(   fieldIndx,node) = u_HD%Morison%LumpedMesh%RotationAcc(   fieldIndx,node) + u_delta(n)        
      CASE (11) !Module/Mesh/Field: u_HD%Morison%DistribMesh%TranslationAcc = 11
         u_HD%Morison%DistribMesh%TranslationAcc(fieldIndx,node) = u_HD%Morison%DistribMesh%TranslationAcc(fieldIndx,node) + u_delta(n)        
      CASE (12) !Module/Mesh/Field: u_HD%Morison%DistribMesh%RotationAcc = 12
         u_HD%Morison%DistribMesh%RotationAcc(   fieldIndx,node) = u_HD%Morison%DistribMesh%RotationAcc(   fieldIndx,node) + u_delta(n)      
      CASE (13) !Module/Mesh/Field: u_HD%Mesh%TranslationAcc = 13
         u_HD%Mesh%TranslationAcc(   fieldIndx,node) = u_HD%Mesh%TranslationAcc(   fieldIndx,node) + u_delta(n)      
      CASE (14) !Module/Mesh/Field: u_HD%Mesh%RotationAcc = 14
         u_HD%Mesh%RotationAcc(   fieldIndx,node) = u_HD%Mesh%RotationAcc(   fieldIndx,node) + u_delta(n)     
         
      CASE (15) !Module/Mesh/Field: u_BD(1)%RootMotion%TranslationAcc = 15 (k=1)
         u_BD(1)%RootMotion%TranslationAcc(fieldIndx,node) = u_BD(1)%RootMotion%TranslationAcc(fieldIndx,node) + u_delta(n)      
      CASE (16) !Module/Mesh/Field: u_BD(1)%RootMotion%RotationAcc = 16 (k=1)
         u_BD(1)%RootMotion%RotationAcc(   fieldIndx,node) = u_BD(1)%RootMotion%RotationAcc(   fieldIndx,node) + u_delta(n)     
      CASE (17) !Module/Mesh/Field: u_BD(2)%RootMotion%TranslationAcc = 17 (k=2)
         u_BD(2)%RootMotion%TranslationAcc(fieldIndx,node) = u_BD(2)%RootMotion%TranslationAcc(fieldIndx,node) + u_delta(n)      
      CASE (18) !Module/Mesh/Field: u_BD(2)%RootMotion%RotationAcc = 18 (k=2)
         u_BD(2)%RootMotion%RotationAcc(   fieldIndx,node) = u_BD(2)%RootMotion%RotationAcc(   fieldIndx,node) + u_delta(n)     
      CASE (19) !Module/Mesh/Field: u_BD(3)%RootMotion%TranslationAcc = 19 (k=3)
         u_BD(3)%RootMotion%TranslationAcc(fieldIndx,node) = u_BD(3)%RootMotion%TranslationAcc(fieldIndx,node) + u_delta(n)      
      CASE (20) !Module/Mesh/Field: u_BD(3)%RootMotion%RotationAcc = 20 (k=3)
         u_BD(3)%RootMotion%RotationAcc(   fieldIndx,node) = u_BD(3)%RootMotion%RotationAcc(   fieldIndx,node) + u_delta(n)     
         
      
      END SELECT
                                   
   END DO
   
        
END SUBROUTINE Add_ED_SD_HD_BD_u_delta
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Perturb_u_ED_SD_HD_BD( p_FAST, Jac_u_indx, n, u_perturb, u_ED_perturb, u_SD_perturb, u_HD_perturb, u_BD_perturb, perturb )
! This routine perturbs the nth element of the u array (and ED input/output it corresponds to)
!...............................................................................................................................
   TYPE(FAST_ParameterType)            , INTENT(IN   ) :: p_FAST                ! Glue-code simulation parameters
   INTEGER( IntKi )                    , INTENT(IN   ) :: Jac_u_indx(:,:)
   INTEGER( IntKi )                    , INTENT(IN   ) :: n
   REAL( ReKi )                        , INTENT(INOUT) :: u_perturb(:)
   TYPE(ED_InputType),        OPTIONAL , INTENT(INOUT) :: u_ED_perturb           ! ED System inputs  (needed only when                                  1 <= n <= NumEDNodes)
   TYPE(SD_InputType),        OPTIONAL , INTENT(INOUT) :: u_SD_perturb           ! SD System inputs  (needed only when NumEDNodes                      +1 <= n <= NumEDNodes+NumSDNodes) [if SD is used] 
   TYPE(HydroDyn_InputType),  OPTIONAL , INTENT(INOUT) :: u_HD_perturb           ! HD System inputs  (needed only when NumEDNodes+NumSDNodes           +1 <= n <= NumEDNodes+NumSDNodes+NumHDNodes) [if HD is used and SD is used. if SD not used, 
   TYPE(BD_InputType),        OPTIONAL , INTENT(INOUT) :: u_BD_perturb           ! BD System inputs  (needed only when NumEDNodes+NumSDNodes+NumHDNodes+1 <= n <= inf) [if BD is used]
   REAL( ReKi )                        , INTENT(  OUT) :: perturb
   
   ! local variables
   INTEGER                                             :: fieldIndx
   INTEGER                                             :: node
   
      
   fieldIndx = Jac_u_indx(n,2) 
   node      = Jac_u_indx(n,3) 
   
      ! determine which mesh we're trying to perturb and perturb the input:
   SELECT CASE( Jac_u_indx(n,1) )
      
   CASE ( 1) !Module/Mesh/Field: u_ED%PlatformPtMesh%Force = 1
      perturb = GetPerturb( u_ED_perturb%PlatformPtMesh%Force(fieldIndx , node) )
      u_ED_perturb%PlatformPtMesh%Force( fieldIndx,node) = u_ED_perturb%PlatformPtMesh%Force( fieldIndx,node) + perturb * p_FAST%UJacSclFact       
   CASE ( 2) !Module/Mesh/Field: u_ED%PlatformPtMesh%Moment = 2
      perturb = GetPerturb( u_ED_perturb%PlatformPtMesh%Moment(fieldIndx , node) )
      u_ED_perturb%PlatformPtMesh%Moment(fieldIndx,node) = u_ED_perturb%PlatformPtMesh%Moment(fieldIndx,node) + perturb * p_FAST%UJacSclFact             
   CASE ( 3) !Module/Mesh/Field: u_ED%HubPtMesh%Force = 3
      perturb = GetPerturb( u_ED_perturb%HubPtLoad%Force(fieldIndx , node) )
      u_ED_perturb%HubPtLoad%Force(   fieldIndx,node) = u_ED_perturb%HubPtLoad%Force(   fieldIndx,node) + perturb * p_FAST%UJacSclFact     
   CASE ( 4) !Module/Mesh/Field: u_ED%HubPtMesh%Moment = 4
      perturb = GetPerturb( u_ED_perturb%HubPtLoad%Moment(fieldIndx , node) )
      u_ED_perturb%HubPtLoad%Moment(  fieldIndx,node) = u_ED_perturb%HubPtLoad%Moment(  fieldIndx,node) + perturb * p_FAST%UJacSclFact     
               
   CASE ( 5) !Module/Mesh/Field: u_SD%TPMesh%TranslationAcc = 5
      perturb = GetPerturb( u_SD_perturb%TPMesh%TranslationAcc(fieldIndx , node) )
      u_SD_perturb%TPMesh%TranslationAcc(fieldIndx,node) = u_SD_perturb%TPMesh%TranslationAcc(fieldIndx,node) + perturb        
   CASE ( 6) !Module/Mesh/Field: u_SD%TPMesh%RotationAcc = 6
      perturb = GetPerturb( u_SD_perturb%TPMesh%RotationAcc(fieldIndx , node) )
      u_SD_perturb%TPMesh%RotationAcc(   fieldIndx,node) = u_SD_perturb%TPMesh%RotationAcc(   fieldIndx,node) + perturb        
   CASE ( 7) !Module/Mesh/Field: u_SD%LMesh%Force = 7
      perturb = GetPerturb( u_SD_perturb%LMesh%Force(fieldIndx , node) )
      u_SD_perturb%LMesh%Force( fieldIndx,node) = u_SD_perturb%LMesh%Force( fieldIndx,node) + perturb * p_FAST%UJacSclFact     
   CASE ( 8) !Module/Mesh/Field: u_SD%LMesh%Moment = 8
      perturb = GetPerturb( u_SD_perturb%LMesh%Moment(fieldIndx , node) )
      u_SD_perturb%LMesh%Moment(fieldIndx,node) = u_SD_perturb%LMesh%Moment(fieldIndx,node) + perturb * p_FAST%UJacSclFact     
      
   CASE ( 9) !Module/Mesh/Field: u_HD%Morison%LumpedMesh%TranslationAcc = 9
      perturb = GetPerturb( u_HD_perturb%Morison%LumpedMesh%TranslationAcc(fieldIndx , node) )
      u_HD_perturb%Morison%LumpedMesh%TranslationAcc(fieldIndx,node) = u_HD_perturb%Morison%LumpedMesh%TranslationAcc(fieldIndx,node) + perturb        
   CASE ( 10) !Module/Mesh/Field: u_HD%Morison%LumpedMesh%RotationAcc = 10
      perturb = GetPerturb( u_HD_perturb%Morison%LumpedMesh%RotationAcc(fieldIndx , node) )
      u_HD_perturb%Morison%LumpedMesh%RotationAcc(   fieldIndx,node) = u_HD_perturb%Morison%LumpedMesh%RotationAcc(   fieldIndx,node) + perturb        
   CASE (11) !Module/Mesh/Field: u_HD%Morison%DistribMesh%TranslationAcc = 11
      perturb = GetPerturb( u_HD_perturb%Morison%DistribMesh%TranslationAcc(fieldIndx , node) )
      u_HD_perturb%Morison%DistribMesh%TranslationAcc(fieldIndx,node) = u_HD_perturb%Morison%DistribMesh%TranslationAcc(fieldIndx,node) + perturb        
   CASE (12) !Module/Mesh/Field: u_HD%Morison%DistribMesh%RotationAcc = 12
      perturb = GetPerturb( u_HD_perturb%Morison%DistribMesh%RotationAcc(fieldIndx , node) )
      u_HD_perturb%Morison%DistribMesh%RotationAcc(   fieldIndx,node) = u_HD_perturb%Morison%DistribMesh%RotationAcc(   fieldIndx,node) + perturb      
   CASE (13) !Module/Mesh/Field: u_HD%Mesh%TranslationAcc = 13
      perturb = GetPerturb( u_HD_perturb%Mesh%TranslationAcc(fieldIndx , node) )
      u_HD_perturb%Mesh%TranslationAcc(fieldIndx,node) = u_HD_perturb%Mesh%TranslationAcc(fieldIndx,node) + perturb      
   CASE (14) !Module/Mesh/Field: u_HD%Mesh%RotationAcc = 14
      perturb = GetPerturb( u_HD_perturb%Mesh%RotationAcc(fieldIndx , node) )
      u_HD_perturb%Mesh%RotationAcc(   fieldIndx,node) = u_HD_perturb%Mesh%RotationAcc(   fieldIndx,node) + perturb      
            
   CASE (15) !Module/Mesh/Field: u_BD(1)%RootMotion%TranslationAcc = 15 (k=1)
      perturb = GetPerturb( u_BD_perturb%RootMotion%TranslationAcc(fieldIndx , node) )
      u_BD_perturb%RootMotion%TranslationAcc(fieldIndx,node) = u_BD_perturb%RootMotion%TranslationAcc(fieldIndx,node) + perturb      
   CASE (16) !Module/Mesh/Field: u_BD(1)%RootMotion%RotationAcc = 16 (k=1)
      perturb = GetPerturb( u_BD_perturb%RootMotion%RotationAcc(fieldIndx , node) )
      u_BD_perturb%RootMotion%RotationAcc(   fieldIndx,node) = u_BD_perturb%RootMotion%RotationAcc(   fieldIndx,node) + perturb      
   CASE (17) !Module/Mesh/Field: u_BD(2)%RootMotion%TranslationAcc = 17 (k=2)
      perturb = GetPerturb( u_BD_perturb%RootMotion%TranslationAcc(fieldIndx , node) )
      u_BD_perturb%RootMotion%TranslationAcc(fieldIndx,node) = u_BD_perturb%RootMotion%TranslationAcc(fieldIndx,node) + perturb      
   CASE (18) !Module/Mesh/Field: u_BD(2)%RootMotion%RotationAcc = 18 (k=2)
      perturb = GetPerturb( u_BD_perturb%RootMotion%RotationAcc(fieldIndx , node) )
      u_BD_perturb%RootMotion%RotationAcc(   fieldIndx,node) = u_BD_perturb%RootMotion%RotationAcc(   fieldIndx,node) + perturb      
   CASE (19) !Module/Mesh/Field: u_BD(3)%RootMotion%TranslationAcc = 19 (k=3)
      perturb = GetPerturb( u_BD_perturb%RootMotion%TranslationAcc(fieldIndx , node) )
      u_BD_perturb%RootMotion%TranslationAcc(fieldIndx,node) = u_BD_perturb%RootMotion%TranslationAcc(fieldIndx,node) + perturb      
   CASE (20) !Module/Mesh/Field: u_BD(3)%RootMotion%RotationAcc = 20 (k=3)
      perturb = GetPerturb( u_BD_perturb%RootMotion%RotationAcc(fieldIndx , node) )
      u_BD_perturb%RootMotion%RotationAcc(   fieldIndx,node) = u_BD_perturb%RootMotion%RotationAcc(   fieldIndx,node) + perturb      
            
   END SELECT
                                   
   u_perturb(n) = u_perturb(n) + perturb
   
        
END SUBROUTINE Perturb_u_ED_SD_HD_BD
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE IfW_InputSolve( p_FAST, m_FAST, u_IfW, p_IfW, u_AD14, u_AD, y_ED, ErrStat, ErrMsg )

   TYPE(InflowWind_InputType),     INTENT(INOUT)   :: u_IfW       ! The inputs to InflowWind
   TYPE(InflowWind_ParameterType), INTENT(IN   )   :: p_IfW       ! The parameters to InflowWind   
   TYPE(AD14_InputType),           INTENT(IN)      :: u_AD14      ! The input meshes (already calculated) from AeroDyn14
   TYPE(AD_InputType),             INTENT(IN)      :: u_AD        ! The input meshes (already calculated) from AeroDyn
   TYPE(ED_OutputType),            INTENT(IN)      :: y_ED        ! The outputs of the structural dynamics module
   TYPE(FAST_ParameterType),       INTENT(IN   )   :: p_FAST      ! FAST parameter data 
   TYPE(FAST_MiscVarType),         INTENT(IN   )   :: m_FAST      ! misc FAST data, including inputs from external codes like Simulink      
   
   INTEGER(IntKi)                                  :: ErrStat     ! Error status of the operation
   CHARACTER(*)                                    :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Local variables:

   INTEGER(IntKi)                                  :: J           ! Loops through nodes / elements.
   INTEGER(IntKi)                                  :: K           ! Loops through blades.
   INTEGER(IntKi)                                  :: Node        ! Node number for blade/node on mesh


   ErrStat = ErrID_None
   ErrMsg  = ""
      
   
      ! Fill input array for InflowWind
   
   !-------------------------------------------------------------------------------------------------
   ! bjj: this is a hack for now, just to get IfW separate from AD
   ! bjj: it assumes AD 14 is used!!!
   !-------------------------------------------------------------------------------------------------
   Node = 0      
   IF (p_FAST%CompServo == MODULE_SrvD) THEN
      Node = Node + 1
      u_IfW%PositionXYZ(:,Node) = y_ED%HubPtMotion%Position(:,1) ! undisplaced position. Maybe we want to use the displaced position (y_ED%HubPtMotion%TranslationDisp) at some point in time.
   END IF       
            
   IF (p_FAST%CompAero == MODULE_AD14) THEN   
         
      DO K = 1,SIZE(u_AD14%InputMarkers)
         DO J = 1,u_AD14%InputMarkers(K)%nnodes  !this mesh isn't properly set up (it's got the global [absolute] position and no reference position)
            Node = Node + 1
            u_IfW%PositionXYZ(:,Node) = u_AD14%InputMarkers(K)%Position(:,J)
         END DO !J = 1,p%BldNodes ! Loop through the blade nodes / elements
      END DO !K = 1,p%NumBl         
                  
      DO J=1,u_AD14%Twr_InputMarkers%nnodes
         Node = Node + 1      
         u_IfW%PositionXYZ(:,Node) = u_AD14%Twr_InputMarkers%TranslationDisp(:,J) + u_AD14%Twr_InputMarkers%Position(:,J)
      END DO      
         
   ELSEIF (p_FAST%CompAero == MODULE_AD) THEN               
      
      DO K = 1,SIZE(u_AD%BladeMotion)
         DO J = 1,u_AD%BladeMotion(k)%Nnodes
            
            Node = Node + 1
            u_IfW%PositionXYZ(:,Node) = u_AD%BladeMotion(k)%TranslationDisp(:,j) + u_AD%BladeMotion(k)%Position(:,j)
            
         END DO !J = 1,p%BldNodes ! Loop through the blade nodes / elements
      END DO !K = 1,p%NumBl         

      if (allocated(u_AD%InflowOnTower)) then

         DO J=1,u_AD%TowerMotion%nnodes
            Node = Node + 1      
            u_IfW%PositionXYZ(:,Node) = u_AD%TowerMotion%TranslationDisp(:,J) + u_AD%TowerMotion%Position(:,J)
         END DO      
                  
      end if
                        
   END IF
   
               
   CALL IfW_SetExternalInputs( p_IfW, m_FAST, y_ED, u_IfW )


END SUBROUTINE IfW_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE IfW_SetExternalInputs( p_IfW, m_FAST, y_ED, u_IfW )
! This routine sets the inputs required for InflowWind from an external source (Simulink)
!..................................................................................................................................

   TYPE(InflowWind_ParameterType),   INTENT(IN)     :: p_IfW        ! InflowWind parameters
   TYPE(FAST_MiscVarType),           INTENT(IN)     :: m_FAST       ! Glue-code misc variables (including inputs from external sources like Simulink)
   TYPE(ED_OutputType),              INTENT(IN)     :: y_ED         ! The outputs of the structural dynamics module
   TYPE(InflowWind_InputType),       INTENT(INOUT)  :: u_IfW        ! InflowWind Inputs at t

      ! local variables
   INTEGER(IntKi)                                   :: i            ! loop counter
   
   ! bjj: this is a total hack to get the lidar inputs into InflowWind. We should use a mesh to take care of this messiness (and, really this Lidar Focus should come
   ! from Fortran (a scanning pattern or file-lookup inside InflowWind), not MATLAB.
            
   u_IfW%lidar%LidPosition = y_ED%HubPtMotion%Position(:,1) + y_ED%HubPtMotion%TranslationDisp(:,1) & ! rotor apex position (absolute)
                                                            + p_IfW%lidar%RotorApexOffsetPos            ! lidar offset-from-rotor-apex position
      
   u_IfW%lidar%MsrPosition = m_FAST%ExternInput%LidarFocus + u_IfW%lidar%LidPosition


END SUBROUTINE IfW_SetExternalInputs
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AD_InputSolve_IfW( p_FAST, u_AD, y_IfW, y_OpFM, ErrStat, ErrMsg )

      ! Passed variables
   TYPE(FAST_ParameterType),    INTENT(IN   )   :: p_FAST      ! parameter FAST data    
   TYPE(AD_InputType),          INTENT(INOUT)   :: u_AD        ! The inputs to AeroDyn14
   TYPE(InflowWind_OutputType), INTENT(IN)      :: y_IfW       ! The outputs from InflowWind
   TYPE(OpFM_OutputType),       INTENT(IN)      :: y_OpFM      ! outputs from the OpenFOAM integration module
   
   INTEGER(IntKi)                               :: ErrStat     ! Error status of the operation
   CHARACTER(*)                                 :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Local variables:

   INTEGER(IntKi)                               :: J           ! Loops through nodes / elements.
   INTEGER(IntKi)                               :: K           ! Loops through blades.
   INTEGER(IntKi)                               :: NodeNum     ! Node number for blade/node on mesh
   INTEGER(IntKi)                               :: NumBl
   INTEGER(IntKi)                               :: NNodes
   INTEGER(IntKi)                               :: node

   
   ErrStat = ErrID_None
   ErrMsg  = ""
               
   !-------------------------------------------------------------------------------------------------
   ! Set the inputs from inflow wind:
   !-------------------------------------------------------------------------------------------------
   IF (p_FAST%CompInflow == MODULE_IfW) THEN

      if (p_FAST%CompServo == MODULE_SrvD) then
         node = 2
      else
         node = 1
      end if
      
      
      NumBl  = size(u_AD%InflowOnBlade,3)
      Nnodes = size(u_AD%InflowOnBlade,2)
      
      do k=1,NumBl
         do j=1,Nnodes
            u_AD%InflowOnBlade(:,j,k) = y_IfW%VelocityUVW(:,node)
            node = node + 1
         end do
      end do
                  
      if ( allocated(u_AD%InflowOnTower) ) then
         Nnodes = size(u_AD%InflowOnTower,2)
         do j=1,Nnodes
            u_AD%InflowOnTower(:,j) = y_IfW%VelocityUVW(:,node)
            node = node + 1
         end do      
      end if
         
   ELSEIF ( p_FAST%CompInflow == MODULE_OpFM ) THEN
      node = 2 !start of inputs to AD15

      NumBl  = size(u_AD%InflowOnBlade,3)
      Nnodes = size(u_AD%InflowOnBlade,2)
      
      do k=1,NumBl
         do j=1,Nnodes
            u_AD%InflowOnBlade(1,j,k) = y_OpFM%u(node)
            u_AD%InflowOnBlade(2,j,k) = y_OpFM%v(node)
            u_AD%InflowOnBlade(3,j,k) = y_OpFM%w(node)
            node = node + 1
         end do
      end do
                  
      if ( allocated(u_AD%InflowOnTower) ) then
         Nnodes = size(u_AD%InflowOnTower,2)
         do j=1,Nnodes
            u_AD%InflowOnTower(1,j) = y_OpFM%u(node)
            u_AD%InflowOnTower(2,j) = y_OpFM%v(node)
            u_AD%InflowOnTower(3,j) = y_OpFM%w(node)
            node = node + 1
         end do      
      end if      
      
   ELSE
      
      u_AD%InflowOnBlade = 0.0_ReKi ! whole array
      
   END IF
   
   
END SUBROUTINE AD_InputSolve_IfW
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AD_InputSolve_NoIfW( p_FAST, u_AD, y_ED, BD, MeshMapData, ErrStat, ErrMsg )

      ! Passed variables
   TYPE(FAST_ParameterType),    INTENT(IN   )   :: p_FAST      ! parameter FAST data    
   TYPE(AD_InputType),          INTENT(INOUT)   :: u_AD        ! The inputs to AeroDyn14
   TYPE(ED_OutputType),         INTENT(IN)      :: y_ED        ! The outputs from the structural dynamics module
   TYPE(BeamDyn_Data),          INTENT(IN)      :: BD          ! The data from BeamDyn (want the outputs only, but it's in an array)
   TYPE(FAST_ModuleMapType),    INTENT(INOUT)   :: MeshMapData ! Data for mapping between modules
   
   INTEGER(IntKi)                               :: ErrStat     ! Error status of the operation
   CHARACTER(*)                                 :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Local variables:

   INTEGER(IntKi)                               :: K           ! Loops through blades
   INTEGER(IntKi)                               :: NodeNum     ! Node number for blade/node on mesh
   INTEGER(IntKi)                               :: NumBl
   INTEGER(IntKi)                               :: node
   INTEGER(IntKi)                               :: ErrStat2
   CHARACTER(ErrMsgLen)                         :: ErrMsg2 
   CHARACTER(*), PARAMETER                      :: RoutineName = 'AD_InputSolve_NoIfW'

   
   ErrStat = ErrID_None
   ErrMsg  = ""
               
   
   !-------------------------------------------------------------------------------------------------
   ! Set the inputs from ElastoDyn and/or BeamDyn:
   !-------------------------------------------------------------------------------------------------
   
      ! tower
   IF (u_AD%TowerMotion%Committed) THEN
      
      CALL Transfer_Line2_to_Line2( y_ED%TowerLn2Mesh, u_AD%TowerMotion, MeshMapData%ED_L_2_AD_L_T, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%TowerMotion' )      
            
   END IF
   
      
      ! hub
   CALL Transfer_Point_to_Point( y_ED%HubPtMotion, u_AD%HubMotion, MeshMapData%ED_P_2_AD_P_H, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%HubMotion' )      
   
   
      ! blade root   
   DO k=1,size(y_ED%BladeRootMotion)
      CALL Transfer_Point_to_Point( y_ED%BladeRootMotion(k), u_AD%BladeRootMotion(k), MeshMapData%ED_P_2_AD_P_R(k), ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%BladeRootMotion('//trim(num2lstr(k))//')' )      
   END DO
      
   
      ! blades
   IF (p_FAST%CompElast == Module_ED ) THEN
      
!debug_print2 = .true.
      DO k=1,size(y_ED%BladeLn2Mesh)
         CALL Transfer_Line2_to_Line2( y_ED%BladeLn2Mesh(k), u_AD%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%BladeMotion('//trim(num2lstr(k))//')' )   
!debug_print2 = .false.            
      END DO
      
   ELSEIF (p_FAST%CompElast == Module_BD ) THEN
      
         ! get them from BeamDyn
      DO k=1,size(u_AD%BladeMotion)
         CALL Transfer_Line2_to_Line2( BD%y(k)%BldMotion, u_AD%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%BladeMotion('//trim(num2lstr(k))//')' )   
      END DO
      
            
   END IF
   
   
END SUBROUTINE AD_InputSolve_NoIfW
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AD14_InputSolve_IfW( p_FAST, u_AD14, y_IfW, y_OpFM, ErrStat, ErrMsg )
! THIS ROUTINE IS A HACK TO GET THE OUTPUTS FROM ELASTODYN INTO AERODYN14. IT WILL BE REPLACED WHEN THIS CODE LINKS WITH
! AERODYN IN THE NEW FRAMEWORK
!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

      ! Passed variables
   TYPE(FAST_ParameterType),    INTENT(IN   )   :: p_FAST      ! parameter FAST data    
   TYPE(AD14_InputType),        INTENT(INOUT)   :: u_AD14      ! The inputs to AeroDyn14
   TYPE(InflowWind_OutputType), INTENT(IN)      :: y_IfW       ! The outputs from InflowWind
   TYPE(OpFM_OutputType),       INTENT(IN)      :: y_OpFM      ! outputs from the OpenFOAM integration module
   
   INTEGER(IntKi)                               :: ErrStat     ! Error status of the operation
   CHARACTER(*)                                 :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Local variables:

   INTEGER(IntKi)                               :: J           ! Loops through nodes / elements.
   INTEGER(IntKi)                               :: K           ! Loops through blades.
   INTEGER(IntKi)                               :: NodeNum     ! Node number for blade/node on mesh
   INTEGER(IntKi)                               :: NumBl
   INTEGER(IntKi)                               :: BldNodes

   
   ErrStat = ErrID_None
   ErrMsg  = ""

   NumBl    = SIZE(u_AD14%InputMarkers,1)
   BldNodes = u_AD14%InputMarkers(1)%Nnodes
               
   !-------------------------------------------------------------------------------------------------
   ! Set the inputs from inflow wind:
   !-------------------------------------------------------------------------------------------------
   IF (p_FAST%CompInflow == MODULE_IfW) THEN
      IF (p_FAST%CompServo == MODULE_SrvD) THEN
         u_AD14%InflowVelocity = y_IfW%VelocityUVW(:,2:)  ! first point is used for ServoDyn input
      ELSE
         u_AD14%InflowVelocity = y_IfW%VelocityUVW(:,:)  
      END IF               
   ELSEIF ( p_FAST%CompInflow == MODULE_OpFM ) THEN
      u_AD14%InflowVelocity(1,:) = y_OpFM%u(2:)      
      u_AD14%InflowVelocity(2,:) = y_OpFM%v(2:)
      u_AD14%InflowVelocity(3,:) = y_OpFM%w(2:)  
   ELSE
      u_AD14%InflowVelocity = 0.0_ReKi           ! whole array
   END IF
      
   u_AD14%AvgInfVel = y_IfW%DiskVel
   
   
END SUBROUTINE AD14_InputSolve_IfW
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AD14_InputSolve_NoIfW( p_FAST, u_AD14, y_ED, MeshMapData, ErrStat, ErrMsg )
! THIS ROUTINE IS A HACK TO GET THE OUTPUTS FROM ELASTODYN INTO AERODYN14. IT WILL BE REPLACED WHEN THIS CODE LINKS WITH
! AERODYN IN THE NEW FRAMEWORK
!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

      ! Passed variables
   TYPE(FAST_ParameterType),    INTENT(IN   )   :: p_FAST      ! parameter FAST data    
   TYPE(AD14_InputType),        INTENT(INOUT)   :: u_AD14      ! The inputs to AeroDyn14
   TYPE(ED_OutputType),         INTENT(IN)      :: y_ED        ! The outputs from the structural dynamics module
   TYPE(FAST_ModuleMapType),    INTENT(INOUT)   :: MeshMapData ! Data for mapping between modules
   
   INTEGER(IntKi)                               :: ErrStat     ! Error status of the operation
   CHARACTER(*)                                 :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Local variables:

   INTEGER(IntKi)                               :: J           ! Loops through nodes / elements.
   INTEGER(IntKi)                               :: K           ! Loops through blades.
   INTEGER(IntKi)                               :: NodeNum     ! Node number for blade/node on mesh
   INTEGER(IntKi)                               :: NumBl
   INTEGER(IntKi)                               :: BldNodes

   
   ErrStat = ErrID_None
   ErrMsg  = ""

   NumBl    = SIZE(u_AD14%InputMarkers,1)
   BldNodes = u_AD14%InputMarkers(1)%Nnodes
   
               
   !-------------------------------------------------------------------------------------------------
   ! Blade positions, orientations, and velocities:
   !-------------------------------------------------------------------------------------------------
   IF (p_FAST%CompElast == Module_ED) THEN
      DO K = 1,NumBl !p%NumBl ! Loop through all blades
      
         !CALL Transfer_Line2_to_Line2( y_ED%BladeLn2Mesh(K), u_AD%InputMarkers(K), MeshMapData%BDED_L_2_AD_L_B(K), ErrStat, ErrMsg )
         !   IF (ErrStat >= AbortErrLev ) RETURN
         
         u_AD14%InputMarkers(K)%RotationVel = 0.0_ReKi ! bjj: we don't need this field
      
         DO J = 1,BldNodes !p%BldNodes ! Loop through the blade nodes / elements

            NodeNum = J         ! note that this assumes ED has same discretization as AD
         
            u_AD14%InputMarkers(K)%Position(:,J)       = y_ED%BladeLn2Mesh(K)%TranslationDisp(:,NodeNum) + y_ED%BladeLn2Mesh(K)%Position(:,NodeNum) 
            u_AD14%InputMarkers(K)%Orientation(:,:,J)  = y_ED%BladeLn2Mesh(K)%Orientation(:,:,NodeNum)
            u_AD14%InputMarkers(K)%TranslationVel(:,J) = y_ED%BladeLn2Mesh(K)%TranslationVel(:,NodeNum)
            u_AD14%InputMarkers(K)%TranslationAcc(:,J) = y_ED%BladeLn2Mesh(K)%TranslationAcc(:,NodeNum)
                  
         END DO !J = 1,p%BldNodes ! Loop through the blade nodes / elements
      END DO !K = 1,p%NumBl
   ELSE
      ! just leave them as the initial guesses?
      DO K = 1,NumBl
         u_AD14%InputMarkers(K)%RotationVel    = 0.0_ReKi
         u_AD14%InputMarkers(K)%TranslationVel = 0.0_ReKi
         u_AD14%InputMarkers(K)%TranslationAcc = 0.0_ReKi
      END DO
      
   END IF
   
   !-------------------------------------------------------------------------------------------------
   ! Hub positions, orientations, and velocities:
   !  (note that these may have to be adjusted in ElastoDyn as AeroDyn gets rewritten)
   !-------------------------------------------------------------------------------------------------
   u_AD14%TurbineComponents%Hub%Position    = y_ED%HubPtMotion%TranslationDisp(:,1) +  y_ED%HubPtMotion%Position(:,1)
   u_AD14%TurbineComponents%Hub%Orientation = y_ED%HubPtMotion%Orientation(:,:,1)   
   u_AD14%TurbineComponents%Hub%RotationVel = y_ED%HubPtMotion%RotationVel(:,1)
   
   u_AD14%TurbineComponents%Hub%TranslationVel = 0.0_ReKi !bjj we don't need this field
   !-------------------------------------------------------------------------------------------------
   ! Blade root orientations:
   !-------------------------------------------------------------------------------------------------
   
   DO K=1,NumBl
      u_AD14%TurbineComponents%Blade(K)%Orientation = y_ED%BladeRootMotion14%Orientation(:,:,K)
      
      u_AD14%TurbineComponents%Blade(K)%Position       = 0.0_ReKi !bjj we don't need this field
      u_AD14%TurbineComponents%Blade(K)%RotationVel    = 0.0_ReKi !bjj we don't need this field
      u_AD14%TurbineComponents%Blade(K)%TranslationVel = 0.0_ReKi !bjj we don't need this field
   END DO
            

   !-------------------------------------------------------------------------------------------------
   ! RotorFurl position, orientation, rotational velocity:
   !-------------------------------------------------------------------------------------------------

   u_AD14%TurbineComponents%RotorFurl%Position    = y_ED%RotorFurlMotion14%TranslationDisp(:,1) + y_ED%RotorFurlMotion14%Position(:,1)  
   u_AD14%TurbineComponents%RotorFurl%Orientation = y_ED%RotorFurlMotion14%Orientation(:,:,1)         
   u_AD14%TurbineComponents%RotorFurl%RotationVel = y_ED%RotorFurlMotion14%RotationVel(:,1)
   u_AD14%TurbineComponents%RotorFurl%TranslationVel = 0.0_ReKi !bjj we don't need this field
   
   !-------------------------------------------------------------------------------------------------
   ! Nacelle position, orientation, rotational velocity:
   !-------------------------------------------------------------------------------------------------      

   u_AD14%TurbineComponents%Nacelle%Position    = y_ED%NacelleMotion%TranslationDisp(:,1) + y_ED%NacelleMotion%Position(:,1)
   u_AD14%TurbineComponents%Nacelle%Orientation = y_ED%NacelleMotion%Orientation(:,:,1)      
   u_AD14%TurbineComponents%Nacelle%RotationVel = y_ED%NacelleMotion%RotationVel(:,1)  
   u_AD14%TurbineComponents%Nacelle%TranslationVel = 0.0_ReKi !bjj we don't need this field
   
   !-------------------------------------------------------------------------------------------------
   ! Tower base position, rotational velocity:
   !-------------------------------------------------------------------------------------------------      
   
   
      ! Tower base position should be rT(0) instead of rZ, but AeroDyn needs this for
      ! the HubVDue2Yaw calculation:
   u_AD14%TurbineComponents%Tower%Position     = y_ED%TowerBaseMotion14%TranslationDisp(:,1) + y_ED%TowerBaseMotion14%Position(:,1)
   u_AD14%TurbineComponents%Tower%RotationVel  = y_ED%TowerBaseMotion14%RotationVel(:,1)
   u_AD14%TurbineComponents%Tower%Orientation    = 0.0_ReKi !bjj we don't need this field
   u_AD14%TurbineComponents%Tower%TranslationVel = 0.0_ReKi !bjj we don't need this field
  

   !-------------------------------------------------------------------------------------------------
   ! Tower mesh info: Twr_InputMarkers
   !-------------------------------------------------------------------------------------------------      
   
   IF ( u_AD14%Twr_InputMarkers%Committed ) THEN
      
      !CALL Transfer_Line2_to_Line2( y_ED%TowerLn2Mesh, u_AD%Twr_InputMarkers, MeshMapData%ED_L_2_AD_L_T, ErrStat, ErrMsg )
      !   IF (ErrStat >= AbortErrLev ) RETURN   
      
      J = u_AD14%Twr_InputMarkers%NNodes
      u_AD14%Twr_InputMarkers%TranslationDisp = y_ED%TowerLn2Mesh%TranslationDisp(:,1:J)
      u_AD14%Twr_InputMarkers%Orientation     = y_ED%TowerLn2Mesh%Orientation    (:,:,1:J)
      
   END IF
      
   !-------------------------------------------------------------------------------------------------
   ! If using MulTabLoc feature, set it here:
   !-------------------------------------------------------------------------------------------------      
   
   !  u_AD14%MulTabLoc(IElements,IBlades) = ???
   
END SUBROUTINE AD14_InputSolve_NoIfW
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AD_SetInitInput(InitInData_AD14, InitOutData_ED, y_ED, p_FAST, ErrStat, ErrMsg)
! This subroutine sets up the information needed to initialize AeroDyn, then initializes AeroDyn
!----------------------------------------------------------------------------------------------------

   ! Passed variables:
   TYPE(AD14_InitInputType),INTENT(INOUT) :: InitInData_AD14  ! The initialization input to AeroDyn14
   TYPE(ED_InitOutputType), INTENT(IN)    :: InitOutData_ED   ! The initialization output from structural dynamics module
   TYPE(ED_OutputType),     INTENT(IN)    :: y_ED             ! The outputs of the structural dynamics module (meshes with position/RefOrientation set)
   TYPE(FAST_ParameterType),INTENT(IN)    :: p_FAST           ! The parameters of the glue code
   INTEGER(IntKi)                         :: ErrStat          ! Error status of the operation
   CHARACTER(*)                           :: ErrMsg           ! Error message if ErrStat /= ErrID_None

      ! Local variables

   !TYPE(AD_InitOptions)       :: ADOptions                  ! Options for AeroDyn

   INTEGER                    :: K


   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
      ! Set up the AeroDyn parameters
   InitInData_AD14%ADFileName   = p_FAST%AeroFile
   InitInData_AD14%OutRootName  = p_FAST%OutFileRoot
   InitInData_AD14%WrSumFile    = p_FAST%SumPrint      
   InitInData_AD14%NumBl        = InitOutData_ED%NumBl
   InitInData_AD14%UseDWM       = p_FAST%UseDWM
   
   InitInData_AD14%DWM_InitInputs%IfW_InitInputs%InputFileName   = p_FAST%InflowFile
   
      ! Hub position and orientation (relative here, but does not need to be)

   InitInData_AD14%TurbineComponents%Hub%Position(:)      = y_ED%HubPtMotion%Position(:,1) - y_ED%HubPtMotion%Position(:,1)  ! bjj: was 0; mesh was changed by adding p_ED%HubHt to 3rd component
   InitInData_AD14%TurbineComponents%Hub%Orientation(:,:) = y_ED%HubPtMotion%RefOrientation(:,:,1)
   InitInData_AD14%TurbineComponents%Hub%TranslationVel   = 0.0_ReKi ! bjj: we don't need this field
   InitInData_AD14%TurbineComponents%Hub%RotationVel      = 0.0_ReKi ! bjj: we don't need this field

      ! Blade root position and orientation (relative here, but does not need to be)

   IF (.NOT. ALLOCATED( InitInData_AD14%TurbineComponents%Blade ) ) THEN
      ALLOCATE( InitInData_AD14%TurbineComponents%Blade( InitInData_AD14%NumBl ), STAT = ErrStat )
      IF ( ErrStat /= 0 ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg = ' Error allocating space for InitInData_AD%TurbineComponents%Blade.'
         RETURN
      ELSE
         ErrStat = ErrID_None !reset to ErrID_None, just in case ErrID_None /= 0
      END IF
   END IF

   DO K=1, InitInData_AD14%NumBl
      InitInData_AD14%TurbineComponents%Blade(K)%Position        = y_ED%BladeRootMotion14%Position(:,K)
      InitInData_AD14%TurbineComponents%Blade(K)%Orientation     = y_ED%BladeRootMotion14%RefOrientation(:,:,K)
      InitInData_AD14%TurbineComponents%Blade(K)%TranslationVel  = 0.0_ReKi ! bjj: we don't need this field
      InitInData_AD14%TurbineComponents%Blade(K)%RotationVel     = 0.0_ReKi ! bjj: we don't need this field      
   END DO
  

      ! Blade length
   IF (p_FAST%CompElast == Module_ED) THEN
      InitInData_AD14%TurbineComponents%BladeLength = InitOutData_ED%BladeLength
   ELSE
      InitInData_AD14%TurbineComponents%BladeLength = 12.573 !13.75700  ! FIX THIS!!!!! GET THIS FROM BEAMDYN (right now is value from Test01)
   END IF
   
   
      ! Tower mesh ( here only because we currently need line2 meshes to contain the same nodes/elements )
      
   InitInData_AD14%NumTwrNodes = y_ED%TowerLn2Mesh%NNodes - 2
   IF (.NOT. ALLOCATED( InitInData_AD14%TwrNodeLocs ) ) THEN
      ALLOCATE( InitInData_AD14%TwrNodeLocs( 3, InitInData_AD14%NumTwrNodes ), STAT = ErrStat )
      IF ( ErrStat /= 0 ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg = ' Error allocating space for InitInData_AD%TwrNodeLocs.'
         RETURN
      ELSE
         ErrStat = ErrID_None
      END IF
   END IF   
   
   IF ( InitInData_AD14%NumTwrNodes > 0 ) THEN
      InitInData_AD14%TwrNodeLocs = y_ED%TowerLn2Mesh%Position(:,1:InitInData_AD14%NumTwrNodes)  ! ED has extra nodes at beginning and top and bottom of tower
   END IF
   
      ! hub height         
   InitInData_AD14%HubHt = InitOutData_ED%HubHt
             

   RETURN
END SUBROUTINE AD_SetInitInput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE WriteInputMeshesToFile(u_ED, u_AD, u_SD, u_HD, u_MAP, u_BD, FileName, ErrStat, ErrMsg) 
   TYPE(ED_InputType),        INTENT(IN)  :: u_ED         
   TYPE(AD_InputType),        INTENT(IN)  :: u_AD  
   TYPE(SD_InputType),        INTENT(IN)  :: u_SD         
   TYPE(HydroDyn_InputType),  INTENT(IN)  :: u_HD         
   TYPE(MAP_InputType),       INTENT(IN)  :: u_MAP         
   TYPE(BD_InputType),        INTENT(IN)  :: u_BD(:)         
   CHARACTER(*),              INTENT(IN)  :: FileName
   
   INTEGER(IntKi)                         :: ErrStat          ! Error status of the operation
   CHARACTER(*)                           :: ErrMsg           ! Error message if ErrStat /= ErrID_None
   
   
      !FileName = TRIM(p_FAST%OutFileRoot)//'.InputMeshes.bin'
   
      
   INTEGER(IntKi)           :: unOut
   !!!INTEGER(IntKi)           :: unOut2
   INTEGER(IntKi)           :: K_local
   INTEGER(B4Ki), PARAMETER :: File_ID = 3
   INTEGER(B4Ki)            :: NumBl
      

      ! Open the binary output file:
   unOut=-1      
   CALL GetNewUnit( unOut, ErrStat, ErrMsg )
   CALL OpenBOutFile ( unOut, TRIM(FileName), ErrStat, ErrMsg )
      IF (ErrStat /= ErrID_None) RETURN
  
   !!!   ! get unit for text output
   !!!CALL GetNewUnit( unOut2, ErrStat, ErrMsg )
            
   ! note that I'm not doing anything with the errors here, so it won't tell
   ! you there was a problem writing the data unless it was the last call.
          
      ! Add a file identification number (in case we ever have to change this):
   WRITE( unOut, IOSTAT=ErrStat )   File_ID
   

      ! Add how many blade meshes there are:
   NumBl =  SIZE(u_ED%BladeLn2Mesh,1)   ! Note that NumBl is B4Ki 
   WRITE( unOut, IOSTAT=ErrStat )   NumBl
      
      ! Add all of the input meshes:
   DO K_local = 1,NumBl
      CALL MeshWrBin( unOut, u_ED%BladeLn2Mesh(K_local), ErrStat, ErrMsg )
      
      !!!write(unOut2,'(A)') '_____________________________________________________________________________'
      !!!write(unOut2,*)     'ElastoDyn blade ', k_local, ' mesh'
      !!!write(unOut2,'(A)') '_____________________________________________________________________________'      
      !!!CALL MeshPrintInfo(unOut2, u_ED%BladeLn2Mesh(K_local))
      !!!write(unOut2,'(A)') '_____________________________________________________________________________'
   END DO            
   CALL MeshWrBin( unOut, u_ED%TowerLn2Mesh,            ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_ED%PlatformPtMesh,          ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_SD%TPMesh,                  ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_SD%LMesh,                   ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_HD%Morison%distribMesh,     ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_HD%Morison%lumpedMesh,      ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_HD%Mesh,                    ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_MAP%PtFairDisplacement,     ErrStat, ErrMsg )
      ! Add how many BD blade meshes there are:
   NumBl =  SIZE(u_BD,1)   ! Note that NumBl is B4Ki 
   WRITE( unOut, IOSTAT=ErrStat )   NumBl
   
   DO K_local = 1,NumBl
      CALL MeshWrBin( unOut, u_BD(K_local)%RootMotion, ErrStat, ErrMsg )
      CALL MeshWrBin( unOut, u_BD(K_local)%DistrLoad, ErrStat, ErrMsg )
   END DO            
      
      ! Add how many AD blade meshes there are:
   NumBl =  SIZE(u_AD%BladeMotion,1)   ! Note that NumBl is B4Ki 
   WRITE( unOut, IOSTAT=ErrStat )   NumBl
   
   DO K_local = 1,NumBl
      CALL MeshWrBin( unOut, u_AD%BladeMotion(k_local), ErrStat, ErrMsg )
   END DO    
      
      ! Close the file
   CLOSE(unOut)
   !!!CLOSE(unOut2)
         
END SUBROUTINE WriteInputMeshesToFile   
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE WriteMotionMeshesToFile(time, y_ED, u_SD, y_SD, u_HD, u_MAP, y_BD, u_BD, UnOut, ErrStat, ErrMsg, FileName) 
   REAL(DbKi),                 INTENT(IN)    :: time
   TYPE(ED_OutputType),        INTENT(IN)    :: y_ED         
   TYPE(SD_InputType),         INTENT(IN)    :: u_SD         
   TYPE(SD_OutputType),        INTENT(IN)    :: y_SD         
   TYPE(HydroDyn_InputType),   INTENT(IN)    :: u_HD         
   TYPE(MAP_InputType),        INTENT(IN)    :: u_MAP         
   TYPE(BD_OutputType),        INTENT(IN)    :: y_BD(:)
   TYPE(BD_InputType),         INTENT(IN)    :: u_BD(:)
   INTEGER(IntKi) ,            INTENT(INOUT) :: unOut
   CHARACTER(*),   OPTIONAL,   INTENT(IN)    :: FileName
   
   INTEGER(IntKi), INTENT(OUT)               :: ErrStat          ! Error status of the operation
   CHARACTER(*)  , INTENT(OUT)               :: ErrMsg           ! Error message if ErrStat /= ErrID_None
   
   
   REAL(R8Ki)               :: t
      
   INTEGER(IntKi)           :: K_local
   INTEGER(B4Ki), PARAMETER :: File_ID = 101
   INTEGER(B4Ki)            :: NumBl
      
   t = time  ! convert to 8-bytes if necessary (DbKi might not be R8Ki)
   
   ! note that I'm not doing anything with the errors here, so it won't tell
   ! you there was a problem writing the data unless it was the last call.
   
   
      ! Open the binary output file and write a header:
   if (unOut<0) then
      CALL GetNewUnit( unOut, ErrStat, ErrMsg )
      
      CALL OpenBOutFile ( unOut, TRIM(FileName), ErrStat, ErrMsg )
         IF (ErrStat /= ErrID_None) RETURN
               
         ! Add a file identification number (in case we ever have to change this):
      WRITE( unOut, IOSTAT=ErrStat )   File_ID
      
         ! Add how many blade meshes there are:
      NumBl =  SIZE(y_ED%BladeLn2Mesh,1)   ! Note that NumBl is B4Ki 
      WRITE( unOut, IOSTAT=ErrStat )   NumBl
      NumBl =  SIZE(y_BD,1)   ! Note that NumBl is B4Ki 
      WRITE( unOut, IOSTAT=ErrStat )   NumBl
   end if
   
   WRITE( unOut, IOSTAT=ErrStat ) t          
   
      ! Add all of the meshes with motions:
   DO K_local = 1,SIZE(y_ED%BladeLn2Mesh,1)
      CALL MeshWrBin( unOut, y_ED%BladeLn2Mesh(K_local), ErrStat, ErrMsg )
   END DO            
   CALL MeshWrBin( unOut, y_ED%TowerLn2Mesh,            ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, y_ED%PlatformPtMesh,          ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_SD%TPMesh,                  ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, y_SD%y2Mesh,                  ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_HD%Morison%distribMesh,     ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_HD%Morison%lumpedMesh,      ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_HD%Mesh,                    ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_MAP%PtFairDisplacement,     ErrStat, ErrMsg )
   DO K_local = 1,SIZE(y_BD,1)
      CALL MeshWrBin( unOut, u_BD(K_local)%RootMotion, ErrStat, ErrMsg )
      CALL MeshWrBin( unOut, y_BD(K_local)%BldMotion,  ErrStat, ErrMsg )
   END DO            
      
   !   
   !   ! Close the file
   !CLOSE(unOut)
   !      
END SUBROUTINE WriteMotionMeshesToFile   
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE WriteMappingTransferToFile(Mesh1_I,Mesh1_O,Mesh2_I,Mesh2_O,Map_Mod1_Mod2,Map_Mod2_Mod1,BinOutputName)
! this routine is used for debugging mesh mapping
!...............................................................................................................................

   TYPE(meshtype),    intent(in) :: mesh1_I
   TYPE(meshtype),    intent(in) :: mesh1_O
   TYPE(meshtype),    intent(in) :: mesh2_I
   TYPE(meshtype),    intent(in) :: mesh2_O 
   
   TYPE(MeshMapType), intent(in) :: Map_Mod1_Mod2        ! Data for mapping meshes from mod1 to mod2
   TYPE(MeshMapType), intent(in) :: Map_Mod2_Mod1        ! Data for mapping meshes from mod1 to mod2

   CHARACTER(*),      INTENT(IN) :: BinOutputName
   
   
   ! local variables:
   TYPE(meshtype)                         :: mesh_Motion_1PT, mesh1_I_1PT, mesh2_O_1PT
   TYPE(MeshMapType)                      :: Map_Mod2_O_1PT, Map_Mod1_I_1PT
      
   INTEGER(IntKi)                         :: i
   INTEGER(IntKi)                         :: un_out
   INTEGER(IntKi)                         :: ErrStat          ! Error status of the operation
   CHARACTER(1024)                        :: ErrMsg           ! Error message if ErrStat /= ErrID_None
   CHARACTER(256)                         :: PrintWarnF, PrintWarnM, TmpValues

   !------------------------------------------------------------------------
   ! Make sure the meshes are committed before checking them:
   !------------------------------------------------------------------------
   
   IF (.NOT. mesh1_I%Committed .OR. .NOT. mesh1_O%Committed ) RETURN
   IF (.NOT. mesh2_I%Committed .OR. .NOT. mesh2_O%Committed ) RETURN
      
   !------------------------------------------------------------------------
   !lump the loads to one point and compare:
   !------------------------------------------------------------------------       
   
   ! create one loads mesh with one point:
   CALL MeshCreate( BlankMesh       = mesh1_I_1PT        &
                  , IOS              = COMPONENT_INPUT   &
                  , NNodes           = 1                 &
                  , Force            = .TRUE.            &
                  , Moment           = .TRUE.            &
                  , ErrStat          = ErrStat           &
                  , ErrMess          = ErrMsg            )
      
   IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
         
         
   CALL MeshPositionNode ( mesh1_I_1PT, 1, (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi/), ErrStat, ErrMsg ) ; IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))             
   CALL MeshConstructElement ( mesh1_I_1PT, ELEMENT_POINT, ErrStat, ErrMsg, 1 );                 IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                                       
   CALL MeshCommit ( mesh1_I_1PT, ErrStat, ErrMsg )                                       ;      IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg)) 
      
   !.....         
   ! create a corresponding motion mesh with one point:
   
   CALL MeshCopy( mesh1_I_1PT, mesh_Motion_1PT, MESH_SIBLING, ErrStat, ErrMsg &
                  , IOS              = COMPONENT_OUTPUT  &
                  , TranslationDisp  = .TRUE.            ) ;                                     IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                     
   !.....         
   ! create a second loads mesh with one point:
   CALL MeshCopy( mesh1_I_1PT, mesh2_O_1PT, MESH_NEWCOPY, ErrStat, ErrMsg )  ! This thinks it's for input, but really it's for output. I don't think it matters...       
       
   !.....         
   ! create the mapping data structures:       
   CALL MeshMapCreate( Mesh1_I, Mesh1_I_1PT, Map_Mod1_I_1PT,  ErrStat, ErrMsg );                 IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   CALL MeshMapCreate( Mesh2_O, Mesh2_O_1PT, Map_Mod2_O_1PT,  ErrStat, ErrMsg );                 IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))

   !.....         
   ! transfer MESH1_I (loads) to single point:
   
   IF ( mesh1_I%ElemTable(ELEMENT_POINT)%nelem > 0 ) THEN
      CALL Transfer_Point_to_Point( Mesh1_I, Mesh1_I_1PT, Map_Mod1_I_1PT,ErrStat,ErrMsg,mesh1_O,mesh_Motion_1PT);  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))         
   ELSEIF ( mesh1_I%ElemTable(ELEMENT_LINE2)%nelem > 0 ) THEN
      CALL Transfer_Line2_to_Point( Mesh1_I, Mesh1_I_1PT, Map_Mod1_I_1PT,ErrStat,ErrMsg,mesh1_O,mesh_Motion_1PT);  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))          
   END IF
   
   !.....         
   ! transfer Mesh2_O (loads) to single point:      
   IF ( Mesh2_O%ElemTable(ELEMENT_POINT)%nelem > 0 ) THEN
      CALL Transfer_Point_to_Point( Mesh2_O, Mesh2_O_1PT, Map_Mod2_O_1PT,ErrStat,ErrMsg,mesh2_I,mesh_Motion_1PT);  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg)) 
   ELSEIF ( Mesh2_O%ElemTable(ELEMENT_LINE2)%nelem > 0 ) THEN 
      CALL Transfer_Line2_to_Point( Mesh2_O, Mesh2_O_1PT, Map_Mod2_O_1PT,ErrStat,ErrMsg,mesh2_I,mesh_Motion_1PT);  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))          
   END IF
   !............
            
   ! dsisplay a warning if the point loads are not equal:         
   PrintWarnF=""
   PrintWarnM=""
   do i=1,3
      if (.NOT. equalrealnos(mesh1_I_1PT%Force( i,1),mesh2_O_1PT%Force( i,1)) ) PrintWarnF=NewLine//"  <----------- WARNING: Forces are not equal ----------->  "//NewLine//NewLine
      if (.NOT. equalrealnos(mesh1_I_1PT%Moment(i,1),mesh2_O_1PT%Moment(i,1)) ) PrintWarnM=NewLine//"  <----------- WARNING: Moments are not equal ----------->  "//NewLine//NewLine
   end do

   
   call wrscr(TRIM(PrintWarnF)//'Total Force:' )
   write(TmpValues,*) mesh1_I_1PT%Force;   call wrscr('     Mesh 1: '//TRIM(TmpValues))
   write(TmpValues,*) mesh2_O_1PT%Force;   call wrscr('     Mesh 2: '//TRIM(TmpValues))
   call wrscr(TRIM(PrintWarnM)//'Total Moment:' )
   write(TmpValues,*) mesh1_I_1PT%Moment;  call wrscr('     Mesh 1: '//TRIM(TmpValues))
   write(TmpValues,*) mesh2_O_1PT%Moment;  call wrscr('     Mesh 2: '//TRIM(TmpValues))
   !............
   
   !------------------------------------------------------------------------
   ! now we'll write all the mesh info to a file for debugging:   
   !------------------------------------------------------------------------

   un_out = -1
   CALL MeshWrBin ( un_out, Mesh1_I,         ErrStat, ErrMsg, BinOutputName);  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   CALL MeshWrBin ( un_out, Mesh1_O,         ErrStat, ErrMsg, BinOutputName);  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   CALL MeshWrBin ( un_out, Mesh2_I,         ErrStat, ErrMsg, BinOutputName);  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   CALL MeshWrBin ( un_out, Mesh2_O,         ErrStat, ErrMsg, BinOutputName);  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      
   CALL MeshWrBin ( un_out, mesh1_I_1PT,     ErrStat, ErrMsg, BinOutputName);  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   CALL MeshWrBin ( un_out, mesh2_O_1PT,     ErrStat, ErrMsg, BinOutputName);  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))         
   CALL MeshWrBin ( un_out, mesh_Motion_1PT, ErrStat, ErrMsg, BinOutputName);  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))   
      
   CALL MeshMapWrBin( un_out, Mesh1_O, Mesh2_I, Map_Mod1_Mod2, ErrStat, ErrMsg, BinOutputName );  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg)) 
   CALL MeshMapWrBin( un_out, Mesh2_O, Mesh1_I, Map_Mod2_Mod1, ErrStat, ErrMsg, BinOutputName );  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg)) 
      
   CALL MeshMapWrBin( un_out, Mesh1_I, Mesh1_I_1PT, Map_Mod1_I_1PT, ErrStat, ErrMsg, BinOutputName );  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg)) 
   CALL MeshMapWrBin( un_out, Mesh2_O, Mesh2_O_1PT, Map_Mod2_O_1PT, ErrStat, ErrMsg, BinOutputName );  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg)) 

   close( un_out )
   
   !------------------------------------------------------------------------
   ! destroy local copies:
   !------------------------------------------------------------------------
   
   CALL MeshDestroy( mesh_Motion_1PT, ErrStat, ErrMsg ); IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   CALL MeshDestroy( mesh1_I_1PT, ErrStat, ErrMsg );     IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   CALL MeshDestroy( mesh2_O_1PT, ErrStat, ErrMsg );     IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   
   call MeshMapDestroy(Map_Mod1_I_1PT, ErrStat, ErrMsg);IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   call MeshMapDestroy(Map_Mod2_O_1PT, ErrStat, ErrMsg);IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   
   

END SUBROUTINE WriteMappingTransferToFile 
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ResetRemapFlags(p_FAST, ED, BD, AD14, AD, HD, SD, SrvD, MAPp, FEAM, MD, IceF, IceD )
! This routine resets the remap flags on all of the meshes
!...............................................................................................................................

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              ! Parameters for the glue code

   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  ! ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  ! BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                ! ServoDyn data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  ! AeroDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                ! AeroDyn14 data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  ! HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  ! SubDyn data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                ! MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                ! FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  ! MoorDyn data
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                ! IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                ! All the IceDyn data used in time-step loop

   !local variable(s)

   INTEGER(IntKi) :: i  ! counter for ice legs
   INTEGER(IntKi) :: k  ! counter for blades
         
   !.....................................................................
   ! Reset each mesh's RemapFlag (after calling all InputSolve routines):
   !.....................................................................     
   
   ! ElastoDyn meshes
   ED%Input( 1)%PlatformPtMesh%RemapFlag        = .FALSE.
   ED%Output(1)%PlatformPtMesh%RemapFlag        = .FALSE.
   ED%Input( 1)%TowerLn2Mesh%RemapFlag          = .FALSE.
   ED%Output(1)%TowerLn2Mesh%RemapFlag          = .FALSE.
   DO K=1,SIZE(ED%Output(1)%BladeRootMotion)
      ED%Output(1)%BladeRootMotion(K)%RemapFlag = .FALSE.      
   END DO
   if (allocated(ED%Input(1)%BladeLn2Mesh)) then   
      DO K=1,SIZE(ED%Input(1)%BladeLn2Mesh)
         ED%Input( 1)%BladeLn2Mesh(K)%RemapFlag    = .FALSE.
         ED%Output(1)%BladeLn2Mesh(K)%RemapFlag    = .FALSE.      
      END DO
   end if
   
   ED%Input( 1)%NacelleLoads%RemapFlag          = .FALSE.
   ED%Output(1)%NacelleMotion%RemapFlag         = .FALSE.
   ED%Input( 1)%HubPtLoad%RemapFlag             = .FALSE.
   ED%Output(1)%HubPtMotion%RemapFlag           = .FALSE.
            
   ! BeamDyn meshes
   IF ( p_FAST%CompElast == Module_BD ) THEN
      DO i=1,p_FAST%nBeams            
         BD%Input(1,i)%RootMotion%RemapFlag = .FALSE.
         BD%Input(1,i)%PointLoad%RemapFlag  = .FALSE.
         BD%Input(1,i)%DistrLoad%RemapFlag  = .FALSE.
             
         BD%y(i)%ReactionForce%RemapFlag    = .FALSE.
         BD%y(i)%BldForce%RemapFlag         = .FALSE.
         BD%y(i)%BldMotion%RemapFlag        = .FALSE.
      END DO                  
   END IF
      
   ! AeroDyn meshes
   IF ( p_FAST%CompAero == Module_AD14 ) THEN
         
      DO k=1,SIZE(AD14%Input(1)%InputMarkers)
         AD14%Input(1)%InputMarkers(k)%RemapFlag = .FALSE.
               AD14%y%OutputLoads(  k)%RemapFlag = .FALSE.
      END DO
                  
      IF (AD14%Input(1)%Twr_InputMarkers%Committed) THEN
         AD14%Input(1)%Twr_InputMarkers%RemapFlag = .FALSE.
                AD14%y%Twr_OutputLoads%RemapFlag  = .FALSE.
      END IF
   ELSEIF ( p_FAST%CompAero == Module_AD ) THEN
               
      AD%Input(1)%HubMotion%RemapFlag = .FALSE.

      IF (AD%Input(1)%TowerMotion%Committed) THEN
          AD%Input(1)%TowerMotion%RemapFlag = .FALSE.
          
         IF (AD%y%TowerLoad%Committed) THEN
                  AD%y%TowerLoad%RemapFlag = .FALSE.
         END IF      
      END IF      
      
      DO k=1,SIZE(AD%Input(1)%BladeMotion)
         AD%Input(1)%BladeRootMotion(k)%RemapFlag = .FALSE.
         AD%Input(1)%BladeMotion(    k)%RemapFlag = .FALSE.
                AD%y%BladeLoad(      k)%RemapFlag = .FALSE.
      END DO
                                    
   END IF
             
   ! HydroDyn
   IF ( p_FAST%CompHydro == Module_HD ) THEN
      IF (HD%Input(1)%Mesh%Committed) THEN
         HD%Input(1)%Mesh%RemapFlag               = .FALSE.
                HD%y%Mesh%RemapFlag               = .FALSE.  
                HD%y%AllHdroOrigin%RemapFlag      = .FALSE.
      END IF
      IF (HD%Input(1)%Morison%LumpedMesh%Committed) THEN
         HD%Input(1)%Morison%LumpedMesh%RemapFlag  = .FALSE.
                HD%y%Morison%LumpedMesh%RemapFlag  = .FALSE.
      END IF
      IF (HD%Input(1)%Morison%DistribMesh%Committed) THEN
         HD%Input(1)%Morison%DistribMesh%RemapFlag = .FALSE.
                HD%y%Morison%DistribMesh%RemapFlag = .FALSE.
      END IF
   END IF

   ! SubDyn
   IF ( p_FAST%CompSub == Module_SD ) THEN
      IF (SD%Input(1)%TPMesh%Committed) THEN
         SD%Input(1)%TPMesh%RemapFlag = .FALSE.
                SD%y%Y1Mesh%RemapFlag = .FALSE.
      END IF    
         
      IF (SD%Input(1)%LMesh%Committed) THEN
         SD%Input(1)%LMesh%RemapFlag  = .FALSE.
                SD%y%Y2Mesh%RemapFlag = .FALSE.
      END IF    
   END IF
      
      
   ! MAP , FEAM , MoorDyn
   IF ( p_FAST%CompMooring == Module_MAP ) THEN
      MAPp%Input(1)%PtFairDisplacement%RemapFlag      = .FALSE.
             MAPp%y%PtFairleadLoad%RemapFlag          = .FALSE.
   ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
      MD%Input(1)%PtFairleadDisplacement%RemapFlag    = .FALSE.
           MD%y%PtFairleadLoad%RemapFlag              = .FALSE.         
   ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
      FEAM%Input(1)%PtFairleadDisplacement%RemapFlag  = .FALSE.
             FEAM%y%PtFairleadLoad%RemapFlag          = .FALSE.         
   END IF
         
   ! IceFloe, IceDyn
   IF ( p_FAST%CompIce == Module_IceF ) THEN
      IF (IceF%Input(1)%iceMesh%Committed) THEN
         IceF%Input(1)%iceMesh%RemapFlag = .FALSE.
                IceF%y%iceMesh%RemapFlag = .FALSE.
      END IF    
   ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
      DO i=1,p_FAST%numIceLegs
         IF (IceD%Input(1,i)%PointMesh%Committed) THEN
            IceD%Input(1,i)%PointMesh%RemapFlag = .FALSE.
                  IceD%y(i)%PointMesh%RemapFlag = .FALSE.
         END IF    
      END DO         
   END IF
      
END SUBROUTINE ResetRemapFlags  
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SetModuleSubstepTime(ModuleID, p_FAST, y_FAST, ErrStat, ErrMsg)
! This module sets the number of subcycles (substeps) for modules, checking to make sure that their requested time step is valid 
!...............................................................................................................................
   INTEGER(IntKi),           INTENT(IN   ) :: ModuleID            ! ID of the module to check time step and set
   TYPE(FAST_ParameterType), INTENT(INOUT) :: p_FAST              ! Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(IN   ) :: y_FAST              ! Output variables for the glue code
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              ! Error message if ErrStat /= ErrID_None

      
   ErrStat = ErrID_None
   ErrMsg  = "" 
   
   IF ( EqualRealNos( p_FAST%dt_module( ModuleID ), p_FAST%dt ) ) THEN
      p_FAST%n_substeps(ModuleID) = 1
   ELSE
      IF ( p_FAST%dt_module( ModuleID ) > p_FAST%dt ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg = "The "//TRIM(y_FAST%Module_Ver(ModuleID)%Name)//" module time step ("//&
                          TRIM(Num2LStr(p_FAST%dt_module( ModuleID )))// &
                    " s) cannot be larger than FAST time step ("//TRIM(Num2LStr(p_FAST%dt))//" s)."
      ELSE
            ! calculate the number of subcycles:
         p_FAST%n_substeps(ModuleID) = NINT( p_FAST%dt / p_FAST%dt_module( ModuleID ) )
            
            ! let's make sure THE module DT is an exact integer divisor of the global (FAST) time step:
         IF ( .NOT. EqualRealNos( p_FAST%dt, p_FAST%dt_module( ModuleID ) * p_FAST%n_substeps(ModuleID) )  ) THEN
            ErrStat = ErrID_Fatal
            ErrMsg  = "The "//TRIM(y_FAST%Module_Ver(ModuleID)%Name)//" module time step ("//&
                              TRIM(Num2LStr(p_FAST%dt_module( ModuleID )))// &
                              " s) must be an integer divisor of the FAST time step ("//TRIM(Num2LStr(p_FAST%dt))//" s)."
         END IF
            
      END IF
   END IF      
                 
   RETURN
      
END SUBROUTINE SetModuleSubstepTime   
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE InitModuleMappings(p_FAST, ED, BD, AD14, AD, HD, SD, SrvD, MAPp, FEAM, MD, IceF, IceD, MeshMapData, ErrStat, ErrMsg)
! This routine initializes all of the mapping data structures needed between the various modules.
!...............................................................................................................................
   
   TYPE(FAST_ParameterType), INTENT(INOUT) :: p_FAST              ! Parameters for the glue code

   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  ! ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  ! BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                ! ServoDyn data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  ! AeroDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                ! AeroDyn14 data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  ! HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  ! SubDyn data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                ! MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                ! FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  ! MoorDyn data
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                ! IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                ! All the IceDyn data used in time-step loop

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         ! Data for mapping between modules
   
   
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              ! Error message if ErrStat /= ErrID_None
   

   INTEGER                                 :: K, i    ! loop counters
   INTEGER                                 :: NumBl   ! number of blades
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'InitModuleMappings'
   !............................................................................................................................
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   NumBl   = SIZE(ED%Output(1)%BladeRootMotion,1)            

   !............................................................................................................................
   ! Create the data structures and mappings in MeshMapType 
   !............................................................................................................................
   
!-------------------------
!  ElastoDyn <-> BeamDyn
!-------------------------
   IF ( p_FAST%CompElast == Module_BD ) THEN
      
      ! Blade meshes: (allocate two mapping data structures to number of blades, then allocate data inside the structures)
      ALLOCATE( MeshMapData%ED_P_2_BD_P(NumBl), MeshMapData%BD_P_2_ED_P(NumBl), STAT=ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error allocating MeshMapData%ED_P_2_BD_P and MeshMapData%BD_P_2_ED_P.', &
                            ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
         
      DO K=1,NumBl         
         CALL MeshMapCreate( ED%Output(1)%BladeRootMotion(K), BD%Input(1,k)%RootMotion, MeshMapData%ED_P_2_BD_P(K), ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_P_2_BD_P('//TRIM(Num2LStr(K))//')' )
         CALL MeshMapCreate( BD%y(k)%ReactionForce, ED%Input(1)%HubPtLoad,  MeshMapData%BD_P_2_ED_P(K), ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':BD_P_2_ED_P('//TRIM(Num2LStr(K))//')' )
      END DO      
      
   END IF
   
   
   
!-------------------------
!  ElastoDyn <-> ServoDyn
!-------------------------
   IF ( SrvD%Input(1)%NTMD%Mesh%Committed ) THEN ! ED-SrvD
         
      CALL MeshMapCreate( ED%Output(1)%NacelleMotion, SrvD%Input(1)%NTMD%Mesh, MeshMapData%ED_P_2_SrvD_P_N, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_P_2_SrvD_P_N' )
      CALL MeshMapCreate( SrvD%y%NTMD%Mesh, ED%Input(1)%NacelleLoads,  MeshMapData%SrvD_P_2_ED_P_N, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SrvD_P_2_ED_P_N' )
   
   END IF   
   
!-------------------------
!  ElastoDyn <-> AeroDyn
!-------------------------
   
   IF ( p_FAST%CompAero == Module_AD14 ) THEN ! ED-AD14
         
      ! Blade meshes: (allocate two mapping data structures to number of blades, then allocate data inside the structures)
      ALLOCATE( MeshMapData%BDED_L_2_AD_L_B(NumBl), MeshMapData%AD_L_2_BDED_L_B(NumBl), STAT=ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error allocating MeshMapData%BDED_L_2_AD_L_B and MeshMapData%AD_L_2_BDED_L_B.', &
                            ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
         
      DO K=1,NumBl         
         CALL MeshMapCreate( ED%Output(1)%BladeLn2Mesh(K), AD14%Input(1)%InputMarkers(K), MeshMapData%BDED_L_2_AD_L_B(K), ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':BDED_L_2_AD_L_B('//TRIM(Num2LStr(K))//')' )
         CALL MeshMapCreate( AD14%y%OutputLoads(K), ED%Input(1)%BladeLn2Mesh(K),  MeshMapData%AD_L_2_BDED_L_B(K), ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':AD_L_2_BDED_L_B('//TRIM(Num2LStr(K))//')' )
      END DO
         
      ! Tower mesh:
      IF ( AD14%Input(1)%Twr_InputMarkers%Committed ) THEN
         CALL MeshMapCreate( ED%Output(1)%TowerLn2Mesh, AD14%Input(1)%Twr_InputMarkers, MeshMapData%ED_L_2_AD_L_T, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_L_2_AD_L_T' )
         CALL MeshMapCreate( AD14%y%Twr_OutputLoads, ED%Input(1)%TowerLn2Mesh,  MeshMapData%AD_L_2_ED_L_T, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':AD_L_2_ED_L_T' )
      END IF
               
      IF (ErrStat >= AbortErrLev ) RETURN
      
   ELSEIF ( p_FAST%CompAero == Module_AD ) THEN ! ED-AD and/or BD-AD

      NumBl = SIZE(AD%Input(1)%BladeRootMotion)            
      
         ! Blade root meshes
      ALLOCATE( MeshMapData%ED_P_2_AD_P_R(NumBl), STAT=ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error allocating MeshMapData%ED_P_2_AD_P_R.', ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF      
         
      DO K=1,NumBl         
         CALL MeshMapCreate( ED%Output(1)%BladeRootMotion(K), AD%Input(1)%BladeRootMotion(K), MeshMapData%ED_P_2_AD_P_R(K), ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_P_2_AD_P_R('//TRIM(Num2LStr(K))//')' )
      END DO
      
      
         ! Hub point mesh
      CALL MeshMapCreate( ED%Output(1)%HubPtMotion, AD%Input(1)%HubMotion, MeshMapData%ED_P_2_AD_P_H, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_P_2_AD_P_H' )
      
      
      ! Blade meshes: (allocate two mapping data structures to number of blades, then allocate data inside the structures)                  
      ALLOCATE( MeshMapData%BDED_L_2_AD_L_B(NumBl), MeshMapData%AD_L_2_BDED_L_B(NumBl), STAT=ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error allocating MeshMapData%BDED_L_2_AD_L_B and MeshMapData%AD_L_2_BDED_L_B.', &
                            ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
         
      IF ( p_FAST%CompElast == Module_ED ) then
         
         DO K=1,NumBl         
            CALL MeshMapCreate( ED%Output(1)%BladeLn2Mesh(K), AD%Input(1)%BladeMotion(K), MeshMapData%BDED_L_2_AD_L_B(K), ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':BDED_L_2_AD_L_B('//TRIM(Num2LStr(K))//')' )
            CALL MeshMapCreate( AD%y%BladeLoad(K), ED%Input(1)%BladeLn2Mesh(K),  MeshMapData%AD_L_2_BDED_L_B(K), ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':AD_L_2_BDED_L_B('//TRIM(Num2LStr(K))//')' )
         END DO
         
      ELSEIF ( p_FAST%CompElast == Module_BD ) then
         
         ! connect AD mesh with BeamDyn
         DO K=1,NumBl         
            CALL MeshMapCreate( BD%y(k)%BldMotion, AD%Input(1)%BladeMotion(K), MeshMapData%BDED_L_2_AD_L_B(K), ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':BDED_L_2_AD_L_B('//TRIM(Num2LStr(K))//')' )
            CALL MeshMapCreate( AD%y%BladeLoad(K), BD%Input(1,k)%DistrLoad,  MeshMapData%AD_L_2_BDED_L_B(K), ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':AD_L_2_BDED_L_B('//TRIM(Num2LStr(K))//')' )
         END DO
         
         ! Blade meshes for load transfer: (allocate meshes at BD input locations for motions transferred from BD output locations)                  
         ALLOCATE( MeshMapData%BD_L_2_BD_L(NumBl), MeshMapData%y_BD_BldMotion_4Loads(NumBl), STAT=ErrStat2 )
            IF ( ErrStat2 /= 0 ) THEN
               CALL SetErrStat( ErrID_Fatal, 'Error allocating MeshMapData%BD_L_2_BD_L and MeshMapData%y_BD_BldMotion_4Loads.', &
                               ErrStat, ErrMsg, RoutineName )
               RETURN
            END IF
         
         DO K=1,NumBl         
               ! create the new mesh:
            CALL MeshCopy ( SrcMesh  = BD%Input(1,k)%DistrLoad &
                          , DestMesh = MeshMapData%y_BD_BldMotion_4Loads(k) &
                          , CtrlCode = MESH_SIBLING     &
                          , IOS      = COMPONENT_OUTPUT &
                          , TranslationDisp = .TRUE.    &
                          , Orientation     = .TRUE.    &
                          , RotationVel     = .TRUE.    &
                          , TranslationVel  = .TRUE.    &
                          , RotationAcc     = .TRUE.    &
                          , TranslationAcc  = .TRUE.    &
                          , ErrStat  = ErrStat2         &
                          , ErrMess  = ErrMsg2          ) 
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         
               IF (ErrStat >= AbortErrLev) RETURN
                                    
               ! create the mapping:
            CALL MeshMapCreate( BD%y(k)%BldMotion, MeshMapData%y_BD_BldMotion_4Loads(k), MeshMapData%BD_L_2_BD_L(K), ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':BD_L_2_BD_L('//TRIM(Num2LStr(K))//')' )         
         END DO
         
      END IF
      
         
      ! Tower mesh:
      IF ( AD%Input(1)%TowerMotion%Committed ) THEN
         CALL MeshMapCreate( ED%Output(1)%TowerLn2Mesh, AD%Input(1)%TowerMotion, MeshMapData%ED_L_2_AD_L_T, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_L_2_AD_L_T' )
            
         IF ( AD%y%TowerLoad%Committed ) THEN            
            CALL MeshMapCreate( AD%y%TowerLoad, ED%Input(1)%TowerLn2Mesh,  MeshMapData%AD_L_2_ED_L_T, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':AD_L_2_ED_L_T' )
         END IF         
      END IF
      
      
   END IF
   
      
      
   IF ( p_FAST%CompHydro == Module_HD ) THEN ! HydroDyn-{ElastoDyn or SubDyn}
         
         
!-------------------------
!  HydroDyn <-> ElastoDyn
!-------------------------            
      IF ( p_FAST%CompSub /= Module_SD ) THEN ! all of these get mapped to ElastoDyn
            
            ! we're just going to assume ED%Input(1)%PlatformPtMesh is committed
               
         IF ( HD%y%AllHdroOrigin%Committed  ) THEN ! meshes for floating
               ! HydroDyn WAMIT point mesh to/from ElastoDyn point mesh
            CALL MeshMapCreate( HD%y%AllHdroOrigin, ED%Input(1)%PlatformPtMesh, MeshMapData%HD_W_P_2_ED_P, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':HD_W_P_2_ED_P' )
            CALL MeshMapCreate( ED%Output(1)%PlatformPtMesh, HD%Input(1)%Mesh, MeshMapData%ED_P_2_HD_W_P, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_P_2_HD_W_P' )
         END IF            
            
            ! ElastoDyn point mesh HydroDyn Morison point mesh (ED sets inputs, but gets outputs from HD%y%AllHdroOrigin in floating case)
         IF ( HD%Input(1)%Morison%LumpedMesh%Committed  ) THEN            
            CALL MeshMapCreate( ED%Output(1)%PlatformPtMesh,  HD%Input(1)%Morison%LumpedMesh, MeshMapData%ED_P_2_HD_M_P, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_P_2_HD_M_P' )                  
         END IF
            
            ! ElastoDyn point mesh to HydroDyn Morison line mesh (ED sets inputs, but gets outputs from  HD%y%AllHdroOriginin floating case)
         IF ( HD%Input(1)%Morison%DistribMesh%Committed ) THEN
            CALL MeshMapCreate( ED%Output(1)%PlatformPtMesh,  HD%Input(1)%Morison%DistribMesh, MeshMapData%ED_P_2_HD_M_L, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_P_2_HD_M_L' )                  
         END IF

                        
      ELSE ! these get mapped to ElastoDyn AND SubDyn (in ED_SD_HD coupling)  ! offshore fixed
            
            ! HydroDyn WAMIT mesh to ElastoDyn point mesh               
         IF ( HD%y%Mesh%Committed  ) THEN

               ! HydroDyn WAMIT point mesh to ElastoDyn point mesh ! meshes for fixed-bottom
            CALL MeshMapCreate( HD%y%Mesh, ED%Input(1)%PlatformPtMesh, MeshMapData%HD_W_P_2_ED_P, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':HD_W_P_2_ED_P' )                  
            CALL MeshMapCreate( ED%Output(1)%PlatformPtMesh, HD%Input(1)%Mesh, MeshMapData%ED_P_2_HD_W_P, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_P_2_HD_W_P' )                  
         END IF             
            
!-------------------------
!  HydroDyn <-> SubDyn
!-------------------------                     
                     
            ! HydroDyn Morison point mesh to SubDyn point mesh
         IF ( HD%y%Morison%LumpedMesh%Committed ) THEN
            
            CALL MeshMapCreate( HD%y%Morison%LumpedMesh, SD%Input(1)%LMesh,  MeshMapData%HD_M_P_2_SD_P, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':HD_M_P_2_SD_P' )                  
            CALL MeshMapCreate( SD%y%y2Mesh,  HD%Input(1)%Morison%LumpedMesh, MeshMapData%SD_P_2_HD_M_P, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SD_P_2_HD_M_P' )                  
                              
         END IF
            
            ! HydroDyn Morison line mesh to SubDyn point mesh
         IF ( HD%y%Morison%DistribMesh%Committed ) THEN
               
            CALL MeshMapCreate( HD%y%Morison%DistribMesh, SD%Input(1)%LMesh,  MeshMapData%HD_M_L_2_SD_P, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':HD_M_L_2_SD_P' )                  
            CALL MeshMapCreate( SD%y%y2Mesh,  HD%Input(1)%Morison%DistribMesh, MeshMapData%SD_P_2_HD_M_L, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SD_P_2_HD_M_L' )                  
                  
         END IF

         
      END IF ! HydroDyn-SubDyn
      
      IF (ErrStat >= AbortErrLev ) RETURN
    
   END IF !HydroDyn-{ElastoDyn or SubDyn}

      
!-------------------------
!  ElastoDyn <-> SubDyn
!-------------------------
   IF ( p_FAST%CompSub == Module_SD ) THEN
                           
      ! NOTE: the MeshMapCreate routine returns fatal errors if either mesh is not committed
      
         ! SubDyn transition piece point mesh to/from ElastoDyn point mesh
      CALL MeshMapCreate( SD%y%Y1mesh, ED%Input(1)%PlatformPtMesh,  MeshMapData%SD_TP_2_ED_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SD_TP_2_ED_P' )                  
      CALL MeshMapCreate( ED%Output(1)%PlatformPtMesh, SD%Input(1)%TPMesh,  MeshMapData%ED_P_2_SD_TP, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_P_2_SD_TP' )                  
   
   END IF ! SubDyn-ElastoDyn      
      
      
   IF ( p_FAST%CompMooring == Module_MAP ) THEN
!-------------------------
!  ElastoDyn <-> MAP
!-------------------------      
      
         ! MAP point mesh to/from ElastoDyn point mesh
      CALL MeshMapCreate( MAPp%y%PtFairleadLoad, ED%Input(1)%PlatformPtMesh,  MeshMapData%Mooring_P_2_ED_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':Mooring_P_2_ED_P' )                  
      CALL MeshMapCreate( ED%Output(1)%PlatformPtMesh, MAPp%Input(1)%PtFairDisplacement,  MeshMapData%ED_P_2_Mooring_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_P_2_Mooring_P' )                  

   ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
!-------------------------
!  ElastoDyn <-> MoorDyn
!-------------------------      
      
         ! MAP point mesh to/from ElastoDyn point mesh
      CALL MeshMapCreate( MD%y%PtFairleadLoad, ED%Input(1)%PlatformPtMesh,  MeshMapData%Mooring_P_2_ED_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':Mooring_P_2_ED_P' )                  
      CALL MeshMapCreate( ED%Output(1)%PlatformPtMesh, MD%Input(1)%PtFairleadDisplacement,  MeshMapData%ED_P_2_Mooring_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_P_2_Mooring_P' )                  
                   
   ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
!-------------------------
!  ElastoDyn <-> FEAMooring
!-------------------------      
      
         ! MAP point mesh to/from ElastoDyn point mesh
      CALL MeshMapCreate( FEAM%y%PtFairleadLoad, ED%Input(1)%PlatformPtMesh,  MeshMapData%Mooring_P_2_ED_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':Mooring_P_2_ED_P' )                  
      CALL MeshMapCreate( ED%Output(1)%PlatformPtMesh, FEAM%Input(1)%PtFairleadDisplacement,  MeshMapData%ED_P_2_Mooring_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_P_2_Mooring_P' )                           
         
   END IF   ! MAP-ElastoDyn ; FEAM-ElastoDyn
            
         
!-------------------------
!  SubDyn <-> IceFloe
!-------------------------      
      
   IF ( p_FAST%CompIce == Module_IceF ) THEN
   
         ! IceFloe iceMesh point mesh to SubDyn LMesh point mesh              
      CALL MeshMapCreate( IceF%y%iceMesh, SD%Input(1)%LMesh,  MeshMapData%IceF_P_2_SD_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':IceF_P_2_SD_P' )                  
         ! SubDyn y2Mesh point mesh to IceFloe iceMesh point mesh 
      CALL MeshMapCreate( SD%y%y2Mesh, IceF%Input(1)%iceMesh,  MeshMapData%SD_P_2_IceF_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SD_P_2_IceF_P' )                  
                              
!-------------------------
!  SubDyn <-> IceDyn
!-------------------------      
      
   ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
   
      ALLOCATE( MeshMapData%IceD_P_2_SD_P( p_FAST%numIceLegs )  , & 
                MeshMapData%SD_P_2_IceD_P( p_FAST%numIceLegs )  , Stat=ErrStat2 )
      IF (ErrStat2 /= 0 ) THEN
         CALL SetErrStat( ErrID_Fatal, 'Unable to allocate IceD_P_2_SD_P and SD_P_2_IceD_P', ErrStat, ErrMsg, RoutineName )                  
         RETURN
      END IF
         
      DO i = 1,p_FAST%numIceLegs
            
            ! IceDyn PointMesh point mesh to SubDyn LMesh point mesh              
         CALL MeshMapCreate( IceD%y(i)%PointMesh, SD%Input(1)%LMesh,  MeshMapData%IceD_P_2_SD_P(i), ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':IceD_P_2_SD_P('//TRIM(num2LStr(i))//')' )                  
            ! SubDyn y2Mesh point mesh to IceDyn PointMesh point mesh 
         CALL MeshMapCreate( SD%y%y2Mesh, IceD%Input(1,i)%PointMesh,  MeshMapData%SD_P_2_IceD_P(i), ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SD_P_2_IceD_P('//TRIM(num2LStr(i))//')' )                  
               
      END DO
                        
   END IF   ! SubDyn-IceFloe
      
   IF (ErrStat >= AbortErrLev ) RETURN   
      
   !............................................................................................................................
   ! Initialize the Jacobian structures:
   !............................................................................................................................
   !IF ( p_FAST%TurbineType == Type_Offshore_Fixed ) THEN ! p_FAST%CompSub == Module_SD .AND. p_FAST%CompHydro == Module_HD 
   IF ( p_FAST%CompSub == Module_SD .OR. p_FAST%CompElast == Module_BD ) THEN  !.OR. p_FAST%CompHydro == Module_HD ) THEN         
      CALL Init_ED_SD_HD_BD_Jacobian( p_FAST, MeshMapData, ED%Input(1)%PlatformPtMesh, SD%Input(1)%TPMesh, SD%Input(1)%LMesh, &
                                    HD%Input(1)%Morison%LumpedMesh, HD%Input(1)%Morison%DistribMesh, HD%Input(1)%Mesh, &
                                    ED%Input(1)%HubPtLoad, BD%Input(1,:),  ErrStat2, ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                 
   ELSEIF ( p_FAST%CompHydro == Module_HD ) THEN
         CALL AllocAry( MeshMapData%Jacobian_ED_SD_HD_BD, SizeJac_ED_HD, SizeJac_ED_HD, 'Jacobian for ED-HD coupling', ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                 
   END IF
   
   IF ( ALLOCATED( MeshMapData%Jacobian_ED_SD_HD_BD ) ) THEN   
      CALL AllocAry( MeshMapData%Jacobian_pivot, SIZE(MeshMapData%Jacobian_ED_SD_HD_BD), 'Pivot array for Jacobian LU decomposition', ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                 
   END IF
   
   IF (ErrStat >= AbortErrLev ) RETURN   
   
   !............................................................................................................................
   ! reset the remap flags (do this before making the copies else the copies will always have remap = true)
   !............................................................................................................................
   CALL ResetRemapFlags(p_FAST, ED, BD, AD14, AD, HD, SD, SrvD, MAPp, FEAM, MD, IceF, IceD )      
            
   !............................................................................................................................
   ! initialize the temporary input meshes (for input-output solves):
   ! (note that we do this after ResetRemapFlags() so that the copies have remap=false)
   !............................................................................................................................
   IF ( p_FAST%CompHydro == Module_HD .OR. p_FAST%CompSub == Module_SD .OR. p_FAST%CompElast == Module_BD ) THEN
                  
         ! Temporary meshes for transfering inputs to ED, HD, BD, and SD
      CALL MeshCopy ( ED%Input(1)%HubPtLoad, MeshMapData%u_ED_HubPtLoad, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_ED_HubPtLoad' )                 
            
      CALL MeshCopy ( ED%Input(1)%PlatformPtMesh, MeshMapData%u_ED_PlatformPtMesh, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_ED_PlatformPtMesh' )                 

      CALL MeshCopy ( ED%Input(1)%PlatformPtMesh, MeshMapData%u_ED_PlatformPtMesh_2, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_ED_PlatformPtMesh_2' )                 
                        
             
      IF ( p_FAST%CompElast == Module_BD ) THEN
      
            ! Temporary meshes for transfering inputs to ED and BD
         CALL MeshCopy ( ED%Input(1)%HubPtLoad, MeshMapData%u_ED_HubPtLoad_2, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_ED_PlatformPtMesh' )     
            
         allocate( MeshMapData%u_BD_RootMotion( p_FAST%nBeams ), STAT = ErrStat2 )
         if (ErrStat2 /= 0) then
            CALL SetErrStat( ErrID_Fatal, "Error allocating u_BD_RootMotion", ErrStat, ErrMsg, RoutineName )     
            return
         end if
         
         do k=1,p_FAST%nBeams
            CALL MeshCopy ( BD%Input(1,k)%RootMotion, MeshMapData%u_BD_RootMotion(k), MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_BD_RootMotion('//trim(num2lstr(k))//')' )                 
         end do
         
                              
      END IF         
         
      IF ( p_FAST%CompSub == Module_SD ) THEN
         
         CALL MeshCopy ( SD%Input(1)%TPMesh, MeshMapData%u_SD_TPMesh, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_SD_TPMesh' )                 
               
         IF ( p_FAST%CompHydro == Module_HD ) THEN
               
            CALL MeshCopy ( SD%Input(1)%LMesh, MeshMapData%u_SD_LMesh, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_SD_LMesh' )                 
                  
            CALL MeshCopy ( SD%Input(1)%LMesh, MeshMapData%u_SD_LMesh_2, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_SD_LMesh_2' )                 
                              
         END IF
               
      END IF
         
      IF ( p_FAST%CompHydro == Module_HD ) THEN
            
         CALL MeshCopy ( HD%Input(1)%Mesh, MeshMapData%u_HD_Mesh, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_HD_Mesh' )                 
                  
         CALL MeshCopy ( HD%Input(1)%Morison%LumpedMesh, MeshMapData%u_HD_M_LumpedMesh, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_HD_M_LumpedMesh' )                 

         CALL MeshCopy ( HD%Input(1)%Morison%DistribMesh, MeshMapData%u_HD_M_DistribMesh, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_HD_M_DistribMesh' )                 
                                    
      END IF
                              
   END IF
   
   

   !............................................................................................................................

      
END SUBROUTINE InitModuleMappings
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE WriteOutputToFile(t_global, p_FAST, y_FAST, ED, BD, AD14, AD, IfW, OpFM, HD, SD, SrvD, MAPp, FEAM, MD, IceF, IceD, ErrStat, ErrMsg)
! This routine determines if it's time to write to the output files, and calls the routine to write to the files
! with the output data. It should be called after all the output solves for a given time have been completed.
!...............................................................................................................................
   REAL(DbKi),               INTENT(IN   ) :: t_global            ! Current global (glue) time step
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              ! Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              ! Output variables for the glue code

   TYPE(ElastoDyn_Data),     INTENT(IN   ) :: ED                  ! ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(IN   ) :: BD                  ! BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(IN   ) :: SrvD                ! ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(IN   ) :: AD14                ! AeroDyn14 data
   TYPE(AeroDyn_Data),       INTENT(IN   ) :: AD                  ! AeroDyn data
   TYPE(InflowWind_Data),    INTENT(IN   ) :: IfW                 ! InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                ! OpenFOAM data
   TYPE(HydroDyn_Data),      INTENT(IN   ) :: HD                  ! HydroDyn data
   TYPE(SubDyn_Data),        INTENT(IN   ) :: SD                  ! SubDyn data
   TYPE(MAP_Data),           INTENT(IN   ) :: MAPp                ! MAP data
   TYPE(FEAMooring_Data),    INTENT(IN   ) :: FEAM                ! FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(IN   ) :: MD                  ! MoorDyn data
   TYPE(IceFloe_Data),       INTENT(IN   ) :: IceF                ! IceFloe data
   TYPE(IceDyn_Data),        INTENT(IN   ) :: IceD                ! All the IceDyn data used in time-step loop

   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              ! Error message if ErrStat /= ErrID_None


   REAL(DbKi)                              :: OutTime             ! Used to determine if output should be generated at this simulation time
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
      
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   IF ( t_global >= p_FAST%TStart )  THEN

         !bjj FIX THIS algorithm!!! this assumes dt_out is an integer multiple of dt; we will probably have to do some interpolation to get these outputs at the times we want them....
         !bjj: perhaps we should do this with integer math on n_t_global now...
      OutTime = NINT( t_global / p_FAST%DT_out ) * p_FAST%DT_out
      IF ( EqualRealNos( t_global, OutTime ) )  THEN

            ! Generate glue-code output file

            CALL WrOutputLine( t_global, p_FAST, y_FAST, IfW%y%WriteOutput, OpFM%y%WriteOutput, ED%Output(1)%WriteOutput, AD%y%WriteOutput, SrvD%y%WriteOutput, &
                  HD%y%WriteOutput, SD%y%WriteOutput, MAPp%y%WriteOutput, FEAM%y%WriteOutput, MD%y%WriteOutput, IceF%y%WriteOutput, IceD%y, BD%y, ErrStat, ErrMsg )
              
            
            IF (p_FAST%WrGraphics) THEN
               CALL WriteMotionMeshesToFile(t_global, ED%Output(1), SD%Input(1), SD%y, HD%Input(1), MAPp%Input(1), BD%y, BD%Input(1,:), y_FAST%UnGra, ErrStat2, ErrMsg2, TRIM(p_FAST%OutFileRoot)//'.gra') 
                  CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, "WriteOutputToFile" )
            END IF
                        
      END IF

   ENDIF
      
            
END SUBROUTINE WriteOutputToFile     
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE CalcOutputs_And_SolveForInputs( n_t_global, this_time, this_state, calcJacobian, NextJacCalcTime, &
               p_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, OpFM, HD, SD, MAPp, FEAM, MD, IceF, IceD, MeshMapData, ErrStat, ErrMsg )
! This subroutine solves the input-output relations for all of the modules. It is a subroutine because it gets done twice--
! once at the start of the n_t_global loop and once in the j_pc loop, using different states.
! *** Note that modules that do not have direct feedthrough should be called first. ***
! also note that this routine uses variables from the main routine (not declared as arguments)
!...............................................................................................................................
   REAL(DbKi)              , intent(in   ) :: this_time           ! The current simulation time (actual or time of prediction)
   INTEGER(IntKi)          , intent(in   ) :: this_state          ! Index into the state array (current or predicted states)
   INTEGER(IntKi)          , intent(in   ) :: n_t_global          ! current time step (used only for SrvD hack)
   LOGICAL                 , intent(inout) :: calcJacobian        ! Should we calculate Jacobians in Option 1?
   REAL(DbKi)              , intent(in   ) :: NextJacCalcTime     ! Time between calculating Jacobians in the HD-ED and SD-ED simulations
      
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              ! Parameters for the glue code
   TYPE(FAST_MiscVarType),   INTENT(IN   ) :: m_FAST              ! Misc variables (including external inputs) for the glue code

   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  ! ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  ! BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                ! ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                ! AeroDyn14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  ! AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 ! InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                ! OpenFOAM data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  ! HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  ! SubDyn data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                ! MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                ! FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  ! Data for the MoorDyn module
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                ! IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                ! All the IceDyn data used in time-step loop

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         ! Data for mapping between modules
   
   
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              ! Error message if ErrStat /= ErrID_None
   
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'CalcOutputs_And_SolveForInputs'
   
   integer :: k
   
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Option 1: solve for consistent inputs and outputs, which is required when Y has direct feedthrough in 
   !           modules coupled together
   ! If you are doing this option at the beginning as well as the end (after option 2), you must initialize the values of
   ! MAPp%y,
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ErrStat = ErrID_None
   ErrMsg  = ""
   
   IF ( EqualRealNos( this_time, NextJacCalcTime ) .OR. NextJacCalcTime < this_time )  THEN
      calcJacobian = .TRUE.
   ELSE         
      calcJacobian = .FALSE.
   END IF
      
#ifdef SOLVE_OPTION_1_BEFORE_2      

   ! This is OPTION 1 before OPTION 2
      
   ! For cases with HydroDyn and/or SubDyn, it calls ED_CalcOuts (a time-sink) 2 times per step/correction (plus the 6 calls when calculating the Jacobian).
   ! For cases without HydroDyn or SubDyn, it calls ED_CalcOuts 1 time per step/correction.
      
   CALL SolveOption1(this_time, this_state, calcJacobian, p_FAST, ED, BD, HD, SD, MAPp, FEAM, MD, IceF, IceD, MeshMapData, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
   CALL SolveOption2(this_time, this_state, p_FAST, m_FAST, ED, BD, AD14, AD, SrvD, IfW, OpFM, MeshMapData, ErrStat2, ErrMsg2, n_t_global < 0)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
                  
#else

   ! This is OPTION 2 before OPTION 1
      
   ! For cases with HydroDyn and/or SubDyn, it calls ED_CalcOuts (a time-sink) 3 times per step/correction (plus the 6 calls when calculating the Jacobian).
   ! In cases without HydroDyn or SubDyn, it is the same as Option 1 before 2 (with 1 call to ED_CalcOuts either way).
      
   ! Option 1 before 2 usually requires a correction step, whereas Option 2 before Option 1 often does not. Thus we are using this option, calling ED_CalcOuts 
   ! 3 times (option 2 before 1 with no correction step) instead of 4 times (option1 before 2 with one correction step). 
   ! Note that this analyisis may change if/when AeroDyn (and ServoDyn?) generate different outputs on correction steps. (Currently, AeroDyn returns old
   ! values until time advances.)

   CALL ED_CalcOutput( this_time, ED%Input(1), ED%p, ED%x(this_state), ED%xd(this_state), ED%z(this_state), ED%OtherSt, ED%Output(1), ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
         
   CALL SolveOption2(this_time, this_state, p_FAST, m_FAST, ED, BD, AD14, AD, SrvD, IfW, OpFM, MeshMapData, ErrStat2, ErrMsg2, n_t_global < 0)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
               
      ! transfer ED outputs to other modules used in option 1:
   CALL Transfer_ED_to_HD_SD_BD_Mooring( p_FAST, ED%Output(1), HD%Input(1), SD%Input(1), MAPp%Input(1), FEAM%Input(1), MD%Input(1), BD%Input(1,:), MeshMapData, ErrStat2, ErrMsg2 )         
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                                     
      
   CALL SolveOption1(this_time, this_state, calcJacobian, p_FAST, ED, BD, HD, SD, MAPp, FEAM, MD, IceF, IceD, MeshMapData, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  

      
      ! use the ElastoDyn outputs from option1 to update the inputs for InflowWind, AeroDyn, and ServoDyn (necessary only if they have states)
                     
   IF ( p_FAST%CompAero == Module_AD14 ) THEN
      
      CALL AD14_InputSolve_NoIfW( p_FAST, AD14%Input(1), ED%Output(1), MeshMapData, ErrStat2, ErrMsg2 )   
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )        
               
         ! because we're not calling InflowWind_CalcOutput or getting new values from OpenFOAM, 
         ! this probably can be skipped
      CALL AD14_InputSolve_IfW( p_FAST, AD14%Input(1), IfW%y, OpFM%y, ErrStat2, ErrMsg2 )   
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                       
         
   ELSEIF ( p_FAST%CompAero == Module_AD ) THEN
      
!debug_print2_unit = 92   
      CALL AD_InputSolve_NoIfW( p_FAST, AD%Input(1), ED%Output(1), BD, MeshMapData, ErrStat2, ErrMsg2 )   
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )        
!debug_print2_unit = 90      

         ! because we're not calling InflowWind_CalcOutput or getting new values from OpenFOAM, 
         ! this probably can be skipped; 
         ! TODO: alternatively, we could call InflowWind_CalcOutput, too.
      CALL AD_InputSolve_IfW( p_FAST, AD%Input(1), IfW%y, OpFM%y, ErrStat2, ErrMsg2 )   
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                       

   END IF

   IF ( p_FAST%CompInflow == Module_IfW ) THEN
      CALL IfW_InputSolve( p_FAST, m_FAST, IfW%Input(1), IfW%p, AD14%Input(1), AD%Input(1), ED%Output(1), ErrStat2, ErrMsg2 )       
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
   ELSE IF ( p_FAST%CompInflow == Module_OpFM ) THEN
   ! OpenFOAM is the driver and it sets these inputs outside of this solve; the OpenFOAM inputs and outputs thus don't change 
   !   in this scenario until OpenFOAM takes another step  **this is a source of error, but it is the way the OpenFOAM-FAST7 coupling
   !   works, so I'm not going to spend time that I don't have now to fix it**
      CALL OpFM_SetInputs( p_FAST, AD14%p, AD14%Input(1), AD14%y, AD%Input(1), AD%y, ED%Output(1), OpFM, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )        
   END IF
   
   
   IF ( p_FAST%CompServo == Module_SrvD  ) THEN         
      CALL SrvD_InputSolve( p_FAST, m_FAST, SrvD%Input(1), ED%Output(1), IfW%y, OpFM%y, MeshmapData, ErrStat2, ErrMsg2 )    ! At initialization, we don't have a previous value, so we'll use the guess inputs instead. note that this violates the framework.... (done for the Bladed DLL)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
   END IF         
                     
#endif
                                                                                                      
   !.....................................................................
   ! Reset each mesh's RemapFlag (after calling all InputSolve routines):
   !.....................................................................              
         
   CALL ResetRemapFlags(p_FAST, ED, BD, AD14, AD, HD, SD, SrvD, MAPp, FEAM, MD, IceF, IceD)         
         
                        
END SUBROUTINE CalcOutputs_And_SolveForInputs  
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SolveOption1(this_time, this_state, calcJacobian, p_FAST, ED, BD, HD, SD, MAPp, FEAM, MD, IceF, IceD, MeshMapData, ErrStat, ErrMsg )
! This routine implements the "option 1" solve for all inputs with direct links to HD, SD, MAP, and the ED platform reference 
! point
!...............................................................................................................................
   REAL(DbKi)              , intent(in   ) :: this_time           ! The current simulation time (actual or time of prediction)
   INTEGER(IntKi)          , intent(in   ) :: this_state          ! Index into the state array (current or predicted states)
   LOGICAL                 , intent(in   ) :: calcJacobian        ! Should we calculate Jacobians in Option 1?

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              ! Parameters for the glue code

   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  ! ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  ! BeamDyn data
   !TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                ! ServoDyn data
   !TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                ! AeroDyn14 data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  ! HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  ! SubDyn data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                ! MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                ! FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  ! MoorDyn data
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                ! IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                ! All the IceDyn data used in time-step loop

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         ! Data for mapping between modules
   
   
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              ! Error message if ErrStat /= ErrID_None
   

   INTEGER                                 :: i                   ! loop counter
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   
   CHARACTER(*), PARAMETER                 :: RoutineName = 'SolveOption1'       
   
   !............................................................................................................................   
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Option 1: solve for consistent inputs and outputs, which is required when Y has direct feedthrough in 
   !           modules coupled together
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   ! Because MAP, FEAM, MoorDyn, IceDyn, and IceFloe do not contain acceleration inputs, we do this outside the DO loop in the ED{_SD}_HD_InputOutput solves.       
   IF ( p_FAST%CompMooring == Module_MAP ) THEN
                  
      CALL MAP_CalcOutput( this_time, MAPp%Input(1), MAPp%p, MAPp%x(this_state), MAPp%xd(this_state), MAPp%z(this_state), &
                            MAPp%OtherSt, MAPp%y, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
         
      CALL MD_CalcOutput( this_time, MD%Input(1), MD%p, MD%x(this_state), MD%xd(this_state), MD%z(this_state), &
                            MD%OtherSt, MD%y, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
   ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
         
      CALL FEAM_CalcOutput( this_time, FEAM%Input(1), FEAM%p, FEAM%x(this_state), FEAM%xd(this_state), FEAM%z(this_state), &
                            FEAM%OtherSt, FEAM%y, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                        
   END IF
      
   IF ( p_FAST%CompIce == Module_IceF ) THEN
                  
      CALL IceFloe_CalcOutput( this_time, IceF%Input(1), IceF%p, IceF%x(this_state), IceF%xd(this_state), IceF%z(this_state), &
                                 IceF%OtherSt, IceF%y, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
         
      DO i=1,p_FAST%numIceLegs                  
         CALL IceD_CalcOutput( this_time, IceD%Input(1,i), IceD%p(i), IceD%x(i,this_state), IceD%xd(i,this_state), &
                                 IceD%z(i,this_state), IceD%OtherSt(i), IceD%y(i), ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
         
   END IF
      
   !
   !   ! User Platform Loading
   !IF ( p_FAST%CompUserPtfmLd ) THEN !bjj: array below won't work... routine needs to be converted to UsrPtfm_CalcOutput()
   !!
   !!   CALL UserPtfmLd ( ED%x(STATE_CURR)%QT(1:6), ED%x(STATE_CURR)%QDT(1:6), t, p_FAST%DirRoot, y_UsrPtfm%AddedMass, (/ y_UsrPtfm%Force,y_UsrPtfm%Moment /) )
   !!   CALL UserPtfmLd ( ED%Output(1)%PlatformPtMesh, t, p_FAST%DirRoot, y_UsrPtfm%AddedMass, ED%u%PlatformPtMesh )
   !!
   !!      ! Ensure that the platform added mass matrix returned by UserPtfmLd, PtfmAM, is symmetric; Abort if necessary:
   !!   IF ( .NOT. IsSymmetric( y_UsrPtfm%AddedMass ) ) THEN
   !!      CALL SetErrStat ( ErrID_Fatal, ' The user-defined platform added mass matrix is unsymmetric.'// &
   !!                        '  Make sure AddedMass returned by UserPtfmLd() is symmetric.', ErrStat, ErrMsg, RoutineName )
   !!   END IF
   !!
   !END IF
      
   IF (ErrStat >= AbortErrLev) RETURN      
   
   IF ( p_FAST%CompSub == Module_SD .OR. p_FAST%CompElast == Module_BD ) THEN !.OR. p_FAST%CompHydro == Module_HD ) THEN
                                 
      CALL ED_SD_HD_BD_InputOutputSolve(  this_time, p_FAST, calcJacobian &
                                    , ED%Input(1),   ED%p, ED%x(  this_state), ED%xd(  this_state), ED%z(  this_state), ED%OtherSt,               ED%Output(1) &
                                    , SD%Input(1),   SD%p, SD%x(  this_state), SD%xd(  this_state), SD%z(  this_state), SD%OtherSt,               SD%y & 
                                    , HD%Input(1),   HD%p, HD%x(  this_state), HD%xd(  this_state), HD%z(  this_state), HD%OtherSt,               HD%y & 
                                    , BD%Input(1,:), BD%p, BD%x(:,this_state), BD%xd(:,this_state), BD%z(:,this_state), BD%OtherSt(:,this_state), BD%y & 
                                    , MAPp%Input(1),   MAPp%y &
                                    , FEAM%Input(1),   FEAM%y &   
                                    , MD%Input(1),     MD%y   &   
                                    , IceF%Input(1),   IceF%y &
                                    , IceD%Input(1,:), IceD%y &    ! bjj: I don't really want to make temp copies of input types. perhaps we should pass the whole Input() structure? (likewise for BD)...
                                    , MeshMapData , ErrStat2, ErrMsg2 )         
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                        
               
   ELSEIF ( p_FAST%CompHydro == Module_HD ) THEN
                                                    
      CALL ED_HD_InputOutputSolve(  this_time, p_FAST, calcJacobian &
                                    , ED%Input(1), ED%p, ED%x(this_state), ED%xd(this_state), ED%z(this_state), ED%OtherSt, ED%Output(1) &
                                    , HD%Input(1), HD%p, HD%x(this_state), HD%xd(this_state), HD%z(this_state), HD%OtherSt, HD%y & 
                                    , MAPp%Input(1), MAPp%y, FEAM%Input(1), FEAM%y, MD%Input(1), MD%y &          
                                    , MeshMapData , ErrStat2, ErrMsg2 )         
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                                                                  
#ifdef SOLVE_OPTION_1_BEFORE_2      
   ELSE 
         
      CALL ED_CalcOutput( this_time, ED%Input(1), ED%p, ED%x(this_state), ED%xd(this_state), ED%z(this_state), &
                           ED%OtherSt, ED%Output(1), ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
#endif         
   END IF ! HD, BD, and/or SD coupled to ElastoDyn
                         
!..................
! Set mooring line and ice inputs (which don't have acceleration fields)
!..................
   
   IF ( p_FAST%CompMooring == Module_MAP ) THEN
         
      ! note: MAP_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)      
      CALL MAP_InputSolve( MAPp%Input(1), ED%Output(1), MeshMapData, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                                 
   ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
         
      ! note: MD_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)      
      CALL MD_InputSolve( MD%Input(1), ED%Output(1), MeshMapData, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                        
   ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
         
      ! note: FEAM_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)      
      CALL FEAM_InputSolve( FEAM%Input(1), ED%Output(1), MeshMapData, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                        
   END IF        
      
   IF ( p_FAST%CompIce == Module_IceF ) THEN
         
      CALL IceFloe_InputSolve(  IceF%Input(1), SD%y, MeshMapData, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                                 
   ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
         
      DO i=1,p_FAST%numIceLegs
            
         CALL IceD_InputSolve(  IceD%Input(1,i), SD%y, MeshMapData, i, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':IceD_InputSolve' )
               
      END DO
         
   END IF        
      
#ifdef DEBUG_MESH_TRANSFER_ICE
      CALL WrScr('********************************************************')
      CALL WrScr('****   IceF to SD point-to-point                   *****')
      CALL WrScr('********************************************************')
      CALL WriteMappingTransferToFile(SD%Input(1)%LMesh, SD%y%Y2Mesh, IceF%Input(1)%iceMesh, IceF%y%iceMesh,&
            MeshMapData%SD_P_2_IceF_P, MeshMapData%IceF_P_2_SD_P, &
            'SD_y2_IceF_Meshes_t'//TRIM(Num2LStr(0))//'.PI.bin' )

         
      CALL WriteMappingTransferToFile(SD%Input(1)%LMesh, SD%y%Y2Mesh, HD%Input(1)%Morison%LumpedMesh, HD%y%Morison%LumpedMesh,&
            MeshMapData%SD_P_2_HD_M_P, MeshMapData%HD_M_P_2_SD_P, &
            'SD_y2_HD_M_L_Meshes_t'//TRIM(Num2LStr(0))//'.PHL.bin' )
         
      CALL WriteMappingTransferToFile(SD%Input(1)%LMesh, SD%y%Y2Mesh, HD%Input(1)%Morison%DistribMesh, HD%y%Morison%DistribMesh,&
            MeshMapData%SD_P_2_HD_M_L, MeshMapData%HD_M_L_2_SD_P, &
            'SD_y2_HD_M_D_Meshes_t'//TRIM(Num2LStr(0))//'.PHD.bin' )
         
         
      !print *
      !pause         
#endif         
                  
END SUBROUTINE SolveOption1
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SolveOption2(this_time, this_state, p_FAST, m_FAST, ED, BD, AD14, AD, SrvD, IfW, OpFM, MeshMapData, ErrStat, ErrMsg, firstCall)
! This routine implements the "option 2" solve for all inputs without direct links to HD, SD, MAP, or the ED platform reference 
! point
!...............................................................................................................................
   LOGICAL                 , intent(in   ) :: firstCall           ! flag to determine how to call ServoDyn (a hack)
   REAL(DbKi)              , intent(in   ) :: this_time           ! The current simulation time (actual or time of prediction)
   INTEGER(IntKi)          , intent(in   ) :: this_state          ! Index into the state array (current or predicted states)

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              ! Parameters for the glue code
   TYPE(FAST_MiscVarType),   INTENT(IN   ) :: m_FAST              ! Misc variables for the glue code (including external inputs)

   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  ! ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  ! BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                ! ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                ! AeroDyn14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  ! AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 ! InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                ! OpenFOAM data

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         ! Data for mapping between modules
   
   
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              ! Error message if ErrStat /= ErrID_None
   

   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   
   CHARACTER(*), PARAMETER                 :: RoutineName = 'SolveOption2'       
   
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Option 2: Solve for inputs based only on the current outputs. This is much faster than option 1 when the coupled modules
   !           do not have direct feedthrough.
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
   ErrStat = ErrID_None
   ErrMsg  = ""
   

      ! find the positions where we want inflow wind in AeroDyn (i.e., set all the motion inputs to AeroDyn)
   IF ( p_FAST%CompAero == Module_AD14 ) THEN 
      
      CALL AD14_InputSolve_NoIfW( p_FAST, AD14%Input(1), ED%Output(1), MeshMapData, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )      
      
   ELSE IF ( p_FAST%CompAero == Module_AD ) THEN 
                        
!debug_print2_unit = 91      
      CALL AD_InputSolve_NoIfW( p_FAST, AD%Input(1), ED%Output(1), BD, MeshMapData, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
!debug_print2_unit = 90      
         
   END IF
   
         
   IF (p_FAST%CompInflow == Module_IfW) THEN
      ! must be done after ED_CalcOutput and before AD_CalcOutput and SrvD
      CALL IfW_InputSolve( p_FAST, m_FAST, IfW%Input(1), IfW%p, AD14%Input(1), AD%Input(1), ED%Output(1), ErrStat2, ErrMsg2 )       
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
      CALL InflowWind_CalcOutput( this_time, IfW%Input(1), IfW%p, IfW%x(this_state), IfW%xd(this_state), IfW%z(this_state), &
                                  IfW%OtherSt, IfW%y, ErrStat2, ErrMsg2 )         
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         
   ELSE IF ( p_FAST%CompInflow == Module_OpFM ) THEN
   ! OpenFOAM is the driver and it computes outputs outside of this solve; the OpenFOAM inputs and outputs thus don't change 
   !   in this scenario until OpenFOAM takes another step  **this is a source of error, but it is the way the OpenFOAM-FAST7 coupling
   !   works, so I'm not going to spend time that I don't have now to fix it**
      CALL OpFM_SetInputs( p_FAST, AD14%p, AD14%Input(1), AD14%y, AD%Input(1), AD%y, ED%Output(1), OpFM, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
      CALL OpFM_SetWriteOutput(OpFM)
      
   END IF
   
   
   IF ( p_FAST%CompAero == Module_AD14 ) THEN 
                        
      CALL AD14_InputSolve_IfW( p_FAST, AD14%Input(1), IfW%y, OpFM%y, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
      CALL AD14_CalcOutput( this_time, AD14%Input(1), AD14%p, AD14%x(this_state), AD14%xd(this_state), AD14%z(this_state), &
                       AD14%OtherSt, AD14%y, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
        
   ELSE IF ( p_FAST%CompAero == Module_AD ) THEN 
                        
      CALL AD_InputSolve_IfW( p_FAST, AD%Input(1), IfW%y, OpFM%y, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
      CALL AD_CalcOutput( this_time, AD%Input(1), AD%p, AD%x(this_state), AD%xd(this_state), AD%z(this_state), &
                       AD%OtherSt, AD%y, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF
      
                       
   IF ( p_FAST%CompServo == Module_SrvD  ) THEN
         
         ! note that the inputs at step(n) for ServoDyn include the outputs from step(n-1)
      IF ( firstCall ) THEN
         CALL SrvD_InputSolve( p_FAST, m_FAST, SrvD%Input(1), ED%Output(1), IfW%y, OpFM%y, MeshMapData, ErrStat2, ErrMsg2 )    ! At initialization, we don't have a previous value, so we'll use the guess inputs instead. note that this violates the framework.... (done for the Bladed DLL)
      ELSE
         CALL SrvD_InputSolve( p_FAST, m_FAST, SrvD%Input(1), ED%Output(1), IfW%y, OpFM%y, MeshMapData, ErrStat2, ErrMsg2, SrvD%y_prev   ) 
      END IF
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      CALL SrvD_CalcOutput( this_time, SrvD%Input(1), SrvD%p, SrvD%x(this_state), SrvD%xd(this_state), SrvD%z(this_state), SrvD%OtherSt, SrvD%y, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   END IF
      
      
      ! User Tower Loading
   IF ( p_FAST%CompUserTwrLd ) THEN !bjj: array below won't work... routine needs to be converted to UsrTwr_CalcOutput()
   !   CALL UserTwrLd ( JNode, X, XD, t, p_FAST%DirRoot, y_UsrTwr%AddedMass(1:6,1:6,J), (/ y_UsrTwr%Force(:,J),y_UsrTwr%Moment(:,J) /) )
   END IF

        
      
   !bjj: note ED%Input(1) may be a sibling mesh of output, but ED%u is not (routine may update something that needs to be shared between siblings)      
   CALL ED_InputSolve( p_FAST, ED%Input(1), ED%Output(1), AD14%y, AD%y, SrvD%y, AD%Input(1), SrvD%Input(1), MeshMapData, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
   CALL BD_InputSolve( p_FAST, BD, AD%y, AD%Input(1), MeshMapData, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
            
END SUBROUTINE SolveOption2
!----------------------------------------------------------------------------------------------------------------------------------  
SUBROUTINE FAST_InitializeAll_T( t_initial, TurbID, Turbine, ErrStat, ErrMsg, InFile, ExternInitData )
! a wrapper routine to call FAST_Initialize a the full-turbine simulation level (makes easier to write top-level driver)

   REAL(DbKi),                        INTENT(IN   ) :: t_initial      ! initial time
   INTEGER(IntKi),                    INTENT(IN   ) :: TurbID         ! turbine Identifier (1-NumTurbines)
   TYPE(FAST_TurbineType),            INTENT(INOUT) :: Turbine        ! all data for one instance of a turbine
   INTEGER(IntKi),                    INTENT(  OUT) :: ErrStat        ! Error status of the operation
   CHARACTER(*),                      INTENT(  OUT) :: ErrMsg         ! Error message if ErrStat /= ErrID_None
   CHARACTER(*),             OPTIONAL,INTENT(IN   ) :: InFile         ! A CHARACTER string containing the name of the primary FAST input file (if not present, we'll get it from the command line)   
   TYPE(FAST_ExternInitType),OPTIONAL,INTENT(IN   ) :: ExternInitData ! Initialization input data from an external source (Simulink)
   
   Turbine%TurbID = TurbID  
   
   
   IF (PRESENT(InFile)) THEN
      IF (PRESENT(ExternInitData)) THEN
         CALL FAST_InitializeAll( t_initial, Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
                     Turbine%ED, Turbine%BD, Turbine%SrvD, Turbine%AD14, Turbine%AD, Turbine%IfW, Turbine%OpFM, &
                     Turbine%HD, Turbine%SD, Turbine%MAP, Turbine%FEAM, Turbine%MD, &
                     Turbine%IceF, Turbine%IceD, Turbine%MeshMapData, ErrStat, ErrMsg, InFile, ExternInitData )
      ELSE         
         CALL FAST_InitializeAll( t_initial, Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
                     Turbine%ED, Turbine%BD, Turbine%SrvD, Turbine%AD14, Turbine%AD, Turbine%IfW, Turbine%OpFM, &
                     Turbine%HD, Turbine%SD, Turbine%MAP, Turbine%FEAM, Turbine%MD, &
                     Turbine%IceF, Turbine%IceD, Turbine%MeshMapData, ErrStat, ErrMsg, InFile  )
      END IF
   ELSE
      CALL FAST_InitializeAll( t_initial, Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
                     Turbine%ED, Turbine%BD, Turbine%SrvD, Turbine%AD14, Turbine%AD, Turbine%IfW, Turbine%OpFM, &
                     Turbine%HD, Turbine%SD, Turbine%MAP, Turbine%FEAM, Turbine%MD, &
                     Turbine%IceF, Turbine%IceD, Turbine%MeshMapData, ErrStat, ErrMsg )
   END IF
   
         
END SUBROUTINE FAST_InitializeAll_T
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FAST_InitializeAll( t_initial, p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, OpFM, HD, SD, MAPp, FEAM, MD, IceF, IceD, &
                               MeshMapData, ErrStat, ErrMsg, InFile, ExternInitData )

   REAL(DbKi),               INTENT(IN   ) :: t_initial           ! initial time
   TYPE(FAST_ParameterType), INTENT(INOUT) :: p_FAST              ! Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              ! Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              ! Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  ! ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  ! BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                ! ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                ! AeroDyn14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  ! AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 ! InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                ! OpenFOAM data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  ! HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  ! SubDyn data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                ! MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                ! FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  ! Data for the MoorDyn module
   
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                ! IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                ! All the IceDyn data used in time-step loop

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         ! Data for mapping between modules
      
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              ! Error message if ErrStat /= ErrID_None
   CHARACTER(*), OPTIONAL,   INTENT(IN   ) :: InFile              ! A CHARACTER string containing the name of the primary FAST input file (if not present, we'll get it from the command line)
   
   TYPE(FAST_ExternInitType), OPTIONAL, INTENT(IN) :: ExternInitData ! Initialization input data from an external source (Simulink)
   
   ! local variables      
   TYPE(ED_InitInputType)                  :: InitInData_ED       ! Initialization input data
   TYPE(ED_InitOutputType)                 :: InitOutData_ED      ! Initialization output data
   
   TYPE(BD_InitInputType)                  :: InitInData_BD       ! Initialization input data
   TYPE(BD_InitOutputType), ALLOCATABLE    :: InitOutData_BD(:)   ! Initialization output data
                                           
   TYPE(SrvD_InitInputType)                :: InitInData_SrvD     ! Initialization input data
   TYPE(SrvD_InitOutputType)               :: InitOutData_SrvD    ! Initialization output data
                                           
   TYPE(AD14_InitInputType)                :: InitInData_AD14     ! Initialization input data
   TYPE(AD14_InitOutputType)               :: InitOutData_AD14    ! Initialization output data
                                           
   TYPE(AD_InitInputType)                  :: InitInData_AD       ! Initialization input data
   TYPE(AD_InitOutputType)                 :: InitOutData_AD      ! Initialization output data
      
   TYPE(InflowWind_InitInputType)          :: InitInData_IfW      ! Initialization input data
   TYPE(InflowWind_InitOutputType)         :: InitOutData_IfW     ! Initialization output data
   
   TYPE(OpFM_InitOutputType)               :: InitOutData_OpFM    ! Initialization output data
      
   TYPE(HydroDyn_InitInputType)            :: InitInData_HD       ! Initialization input data
   TYPE(HydroDyn_InitOutputType)           :: InitOutData_HD      ! Initialization output data
                                           
   TYPE(SD_InitInputType)                  :: InitInData_SD       ! Initialization input data
   TYPE(SD_InitOutputType)                 :: InitOutData_SD      ! Initialization output data
                                           
   TYPE(MAP_InitInputType)                 :: InitInData_MAP      ! Initialization input data
   TYPE(MAP_InitOutputType)                :: InitOutData_MAP     ! Initialization output data
                                           
   TYPE(FEAM_InitInputType)                :: InitInData_FEAM     ! Initialization input data
   TYPE(FEAM_InitOutputType)               :: InitOutData_FEAM    ! Initialization output data
                                           
   TYPE(MD_InitInputType)                  :: InitInData_MD       ! Initialization input data
   TYPE(MD_InitOutputType)                 :: InitOutData_MD      ! Initialization output data
                                           
   TYPE(IceFloe_InitInputType)             :: InitInData_IceF     ! Initialization input data
   TYPE(IceFloe_InitOutputType)            :: InitOutData_IceF    ! Initialization output data
                                           
   TYPE(IceD_InitInputType)                :: InitInData_IceD     ! Initialization input data
   TYPE(IceD_InitOutputType)               :: InitOutData_IceD    ! Initialization output data (each instance will have the same output channels)
       
   REAL(ReKi)                              :: AirDens             ! air density for initialization/normalization of OpenFOAM data
   REAL(DbKi)                              :: dt_IceD             ! tmp dt variable to ensure IceDyn doesn't specify different dt values for different legs (IceDyn instances)
   REAL(DbKi)                              :: dt_BD               ! tmp dt variable to ensure BeamDyn doesn't specify different dt values for different instances
   INTEGER(IntKi)                          :: ErrStat2
   INTEGER(IntKi)                          :: IceDim              ! dimension we're pre-allocating for number of IceDyn legs/instances
   INTEGER(IntKi)                          :: I                   ! generic loop counter
   INTEGER(IntKi)                          :: k                   ! blade loop counter
                                           
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
                                           
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_InitializeAll'       
   
   
   !..........
   ErrStat = ErrID_None
   ErrMsg  = ""
      
   y_FAST%UnSum = -1                                                    ! set the summary file unit to -1 to indicate it's not open
   y_FAST%UnOu  = -1                                                    ! set the text output file unit to -1 to indicate it's not open
   y_FAST%UnGra = -1                                                    ! set the binary graphics output file unit to -1 to indicate it's not open
      
   y_FAST%n_Out = 0                                                     ! set the number of ouptut channels to 0 to indicate there's nothing to write to the binary file
   p_FAST%ModuleInitialized = .FALSE.                                   ! (array initialization) no modules are initialized 
   
      ! Get the current time
   CALL DATE_AND_TIME ( Values=m_FAST%StrtTime )                        ! Let's time the whole simulation
   CALL CPU_TIME ( m_FAST%UsrTime1 )                                    ! Initial time (this zeros the start time when used as a MATLAB function)
   m_FAST%UsrTime1 = MAX( 0.0_ReKi, m_FAST%UsrTime1 )                   ! CPU_TIME: If a meaningful time cannot be returned, a processor-dependent negative value is returned
   

   AbortErrLev            = ErrID_Fatal                                 ! Until we read otherwise from the FAST input file, we abort only on FATAL errors
   m_FAST%t_global        = t_initial - 20.                             ! initialize this to a number < t_initial for error message in ProgAbort
   m_FAST%calcJacobian    = .TRUE.                                      ! we need to calculate the Jacobian
   m_FAST%NextJacCalcTime = m_FAST%t_global                             ! We want to calculate the Jacobian on the first step
   
   
      ! ... Initialize NWTC Library (open console, set pi constants) ...
   CALL NWTC_Init( ProgNameIN=FAST_ver%Name, EchoLibVer=.FALSE. )       ! sets the pi constants, open console for output, etc...

   
      ! ... Open and read input files, initialize global parameters. ...
   IF (PRESENT(InFile)) THEN
      IF (PRESENT(ExternInitData)) THEN
         CALL FAST_Init( p_FAST, y_FAST, t_initial, ErrStat2, ErrMsg2, InFile, ExternInitData%TMax )  ! We have the name of the input file and the simulation length from somewhere else (e.g. Simulink)         
      ELSE         
         CALL FAST_Init( p_FAST, y_FAST, t_initial, ErrStat2, ErrMsg2, InFile )                       ! We have the name of the input file from somewhere else (e.g. Simulink)
      END IF
   ELSE
      CALL FAST_Init( p_FAST, y_FAST, t_initial, ErrStat2, ErrMsg2 )
   END IF
   
   CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF
      
      
   !...............................................................................................................................  
      
   p_FAST%dt_module = p_FAST%dt ! initialize time steps for each module   

   ! ........................
   ! initialize ElastoDyn (must be done first)
   ! ........................
   
   ALLOCATE( ED%Input( p_FAST%InterpOrder+1 ), ED%InputTimes( p_FAST%InterpOrder+1 ), ED%Output( p_FAST%InterpOrder+1 ),STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating ED%Input, ED%Output, and ED%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
      
   
   InitInData_ED%InputFile = p_FAST%EDFile
   IF ( p_FAST%CompAero == Module_AD14 ) THEN
      InitInData_ED%ADInputFile = p_FAST%AeroFile
   ELSE
      InitInData_ED%ADInputFile = ""
   END IF
   
   InitInData_ED%RootName      = p_FAST%OutFileRoot
   InitInData_ED%CompElast     = p_FAST%CompElast == Module_ED

   CALL ED_Init( InitInData_ED, ED%Input(1), ED%p, ED%x(STATE_CURR), ED%xd(STATE_CURR), ED%z(STATE_CURR), ED%OtherSt, &
                  ED%Output(1), p_FAST%dt_module( MODULE_ED ), InitOutData_ED, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
   p_FAST%ModuleInitialized(Module_ED) = .TRUE.
   CALL SetModuleSubstepTime(Module_ED, p_FAST, y_FAST, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      ! bjj: added this check per jmj; perhaps it would be better in ElastoDyn, but I'll leave it here for now:
   IF ( p_FAST%TurbineType == Type_Offshore_Floating ) THEN
      IF ( ED%p%TowerBsHt < 0.0_ReKi .AND. .NOT. EqualRealNos( ED%p%TowerBsHt, 0.0_ReKi ) ) THEN
         CALL SetErrStat(ErrID_Fatal,"ElastoDyn TowerBsHt must not be negative for floating offshore systems.",ErrStat,ErrMsg,RoutineName)
      END IF      
   END IF   
   
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF   
   
   ! ........................
   ! initialize BeamDyn 
   ! ........................
   IF ( p_FAST%CompElast == Module_BD ) THEN      
      p_FAST%nBeams = InitOutData_ED%NumBl          ! initialize number of BeamDyn instances = number of blades      
   ELSE
      p_FAST%nBeams = 0
   END IF

   ALLOCATE( BD%Input( p_FAST%InterpOrder+1, p_FAST%nBeams ), BD%InputTimes( p_FAST%InterpOrder+1, p_FAST%nBeams ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating BD%Input and BD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF  
                        
   ALLOCATE( BD%x(           p_FAST%nBeams,2), &
             BD%xd(          p_FAST%nBeams,2), &
             BD%z(           p_FAST%nBeams,2), &
             BD%OtherSt(     p_FAST%nBeams,2), &
             BD%p(           p_FAST%nBeams  ), &
             BD%u(           p_FAST%nBeams  ), &
             BD%y(           p_FAST%nBeams  ), &
             InitOutData_BD( p_FAST%nBeams  ), &
                                             STAT = ErrStat2 )                                                  
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating BeamDyn state, input, and output data.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF        
   
   IF (p_FAST%CompElast == Module_BD) THEN
      !bjj: FIX ME: TODO: verify that these are correct
      
      InitInData_BD%gravity      = (/ 0.0_ReKi, 0.0_ReKi, -InitOutData_ED%Gravity /)       ! "Gravitational acceleration" m/s^2
      
         ! now initialize IceD for additional legs (if necessary)
      dt_BD = p_FAST%dt_module( MODULE_BD )
                        
      DO k=1,p_FAST%nBeams  
         InitInData_BD%RootName     = TRIM(p_FAST%OutFileRoot)//'.BD'//TRIM( Num2LStr(k) )
         
         
         InitInData_BD%InputFile    = p_FAST%BDBldFile(k)
         
         InitInData_BD%GlbPos       = ED%Output(1)%BladeRootMotion(k)%Position(:,1)          ! {:}    - - "Initial Position Vector of the local blade coordinate system"
         InitInData_BD%GlbRot       = ED%Output(1)%BladeRootMotion(k)%RefOrientation(:,:,1)  ! {:}{:} - - "Initial direction cosine matrix of the local blade coordinate system"
         
         InitInData_BD%RootDisp     = ED%Output(1)%BladeRootMotion(k)%TranslationDisp(:,1)   ! {:}    - - "Initial root displacement"
         InitInData_BD%RootOri      = ED%Output(1)%BladeRootMotion(k)%Orientation(:,:,1)     ! {:}{:} - - "Initial root orientation"
         InitInData_BD%RootVel(1:3) = ED%Output(1)%BladeRootMotion(k)%TranslationVel(:,1)    ! {:}    - - "Initial root velocities and angular veolcities"                  
         InitInData_BD%RootVel(4:6) = ED%Output(1)%BladeRootMotion(k)%RotationVel(:,1)       ! {:}    - - "Initial root velocities and angular veolcities"                  
                           
         CALL BD_Init( InitInData_BD, BD%Input(1,k), BD%p(k),  BD%x(k,STATE_CURR), BD%xd(k,STATE_CURR), BD%z(k,STATE_CURR), &
                           BD%OtherSt(k,STATE_CURR), BD%y(k), dt_BD, InitOutData_BD(k), ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)                        
            
         !bjj: we're going to force this to have the same timestep because I don't want to have to deal with n BD modules with n timesteps.
         IF ( k == 1 ) THEN
            p_FAST%dt_module( MODULE_BD ) = dt_BD
            
            p_FAST%ModuleInitialized(Module_BD) = .TRUE. ! this really should be once per BD instance, but BD doesn't care so I won't go through the effort to track this
            CALL SetModuleSubstepTime(Module_BD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)            
         ELSEIF ( .NOT. EqualRealNos( p_FAST%dt_module( MODULE_BD ),dt_BD )) THEN
            CALL SetErrStat(ErrID_Fatal,"All instances of BeamDyn (one per blade) must have the same time step.",ErrStat,ErrMsg,RoutineName)
         END IF
                           
      END DO
               
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF           
      
   END IF         
      

   ! ........................
   ! initialize AeroDyn 
   ! ........................
   ALLOCATE( AD14%Input( p_FAST%InterpOrder+1 ), AD14%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating AD14%Input and AD14%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
        
   ALLOCATE( AD%Input( p_FAST%InterpOrder+1 ), AD%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating AD%Input and AD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
      
      
   IF ( p_FAST%CompAero == Module_AD14 ) THEN
               
      CALL AD_SetInitInput(InitInData_AD14, InitOutData_ED, ED%Output(1), p_FAST, ErrStat2, ErrMsg2)            ! set the values in InitInData_AD14
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
                                       
      CALL AD14_Init( InitInData_AD14, AD14%Input(1), AD14%p, AD14%x(STATE_CURR), AD14%xd(STATE_CURR), AD14%z(STATE_CURR), AD14%OtherSt, &
                     AD14%y, p_FAST%dt_module( MODULE_AD14 ), InitOutData_AD14, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      p_FAST%ModuleInitialized(Module_AD14) = .TRUE.            
      CALL SetModuleSubstepTime(Module_AD14, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
         ! bjj: this really shouldn't be in the FAST glue code, but I'm going to put this check here so people don't use an invalid model 
         !    and send me emails to debug numerical issues in their results.
      IF ( AD14%p%TwrProps%PJM_Version .AND. p_FAST%TurbineType == Type_Offshore_Floating ) THEN
         CALL SetErrStat(ErrID_Fatal,'AeroDyn v14 tower influence model "NEWTOWER" is invalid for models of floating offshore turbines.',ErrStat,ErrMsg,RoutineName)
      END IF         
            
      AirDens = InitOutData_AD14%AirDens
      
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF       
      
   ELSEIF ( p_FAST%CompAero == Module_AD ) THEN
      
      
         ! set initialization data for AD
      CALL AllocAry( InitInData_AD%BladeRootPosition,      3, InitOutData_ED%NumBl, 'InitInData_AD%BladeRootPosition', errStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL AllocAry( InitInData_AD%BladeRootOrientation,3, 3, InitOutData_ED%NumBl, 'InitInData_AD%BladeRootOrientation', errStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF
            
      InitInData_AD%InputFile          = p_FAST%AeroFile
      InitInData_AD%NumBlades          = InitOutData_ED%NumBl
      InitInData_AD%RootName           = p_FAST%OutFileRoot
      InitInData_AD%HubPosition        = ED%Output(1)%HubPtMotion%Position(:,1)
      InitInData_AD%HubOrientation     = ED%Output(1)%HubPtMotion%RefOrientation(:,:,1)
      
      do k=1,InitOutData_ED%NumBl
         InitInData_AD%BladeRootPosition(:,k)      = ED%Output(1)%BladeRootMotion(k)%Position(:,1)
         InitInData_AD%BladeRootOrientation(:,:,k) = ED%Output(1)%BladeRootMotion(k)%RefOrientation(:,:,1)
      end do
      
            
      CALL AD_Init( InitInData_AD, AD%Input(1), AD%p, AD%x(STATE_CURR), AD%xd(STATE_CURR), AD%z(STATE_CURR), AD%OtherSt, &
                     AD%y, p_FAST%dt_module( MODULE_AD ), InitOutData_AD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      p_FAST%ModuleInitialized(Module_AD) = .TRUE.            
      CALL SetModuleSubstepTime(Module_AD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
                               
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF       
      
      AirDens = InitOutData_AD%AirDens
      
   ELSE
      AirDens = 0.0_ReKi
   END IF ! CompAero
   
               
   ! ........................
   ! initialize InflowWind
   ! ........................   
   ALLOCATE( IfW%Input( p_FAST%InterpOrder+1 ), IfW%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating IfW%Input and IfW%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
              
   IF ( p_FAST%CompInflow == Module_IfW ) THEN
      
      InitInData_IfW%InputFileName    = p_FAST%InflowFile
      InitInData_IfW%RootName         = p_FAST%OutFileRoot
      InitInData_IfW%UseInputFile     = .TRUE.
   
      InitInData_IfW%NumWindPoints = 0      
      IF ( p_FAST%CompServo == Module_SrvD ) InitInData_IfW%NumWindPoints = InitInData_IfW%NumWindPoints + 1
      IF ( p_FAST%CompAero  == Module_AD14 ) THEN
         InitInData_IfW%NumWindPoints = InitInData_IfW%NumWindPoints + InitOutData_ED%NumBl * AD14%Input(1)%InputMarkers(1)%NNodes + AD14%Input(1)%Twr_InputMarkers%NNodes
      ELSEIF ( p_FAST%CompAero  == Module_AD ) THEN
         InitInData_IfW%NumWindPoints = InitInData_IfW%NumWindPoints + AD%Input(1)%TowerMotion%NNodes
         DO k=1,InitOutData_ED%NumBl
            InitInData_IfW%NumWindPoints = InitInData_IfW%NumWindPoints + AD%Input(1)%BladeMotion(k)%NNodes
         END DO
      END IF
      
      ! lidar        
      InitInData_IfW%lidar%Tmax                   = p_FAST%TMax
      InitInData_IfW%lidar%HubPosition            = ED%Output(1)%HubPtMotion%Position(:,1) 
         ! bjj: these should come from an InflowWind input file; I'm hard coding them here for now
      IF ( PRESENT(ExternInitData) ) THEN
         InitInData_IfW%lidar%SensorType          = ExternInitData%SensorType   
         InitInData_IfW%lidar%LidRadialVel        = ExternInitData%LidRadialVel   
         InitInData_IfW%lidar%RotorApexOffsetPos  = 0.0         
         InitInData_IfW%lidar%NumPulseGate        = 0
      ELSE
         InitInData_IfW%lidar%SensorType          = SensorType_None
      END IF
                                     
      CALL InflowWind_Init( InitInData_IfW, IfW%Input(1), IfW%p, IfW%x(STATE_CURR), IfW%xd(STATE_CURR), IfW%z(STATE_CURR), IfW%OtherSt, &
                     IfW%y, p_FAST%dt_module( MODULE_IfW ), InitOutData_IfW, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      p_FAST%ModuleInitialized(Module_IfW) = .TRUE.            
      CALL SetModuleSubstepTime(Module_IfW, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF       
      
   ELSEIF ( p_FAST%CompInflow == Module_OpFM ) THEN
      
         ! set up the data structures for integration with OpenFOAM
      CALL Init_OpFM( p_FAST, AirDens, AD14%Input(1), AD%Input(1), AD%y, ED%Output(1), OpFM, InitOutData_OpFM, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF       
                  
      !bjj: fix me!!! to do
      InitOutData_IfW%WindFileInfo%MWS = 0.0_ReKi
      
   ELSE
      InitOutData_IfW%WindFileInfo%MWS = 0.0_ReKi
   END IF   ! CompInflow
   

   ! ........................
   ! some checks for AeroDyn inputs with the high-speed shaft brake hack in ElastoDyn:
   ! (DO NOT COPY THIS CODE!)
   ! ........................   
   IF ( p_FAST%CompServo == Module_SrvD ) THEN
      
         ! bjj: this is a hack to get high-speed shaft braking in FAST v8
      
      IF ( InitOutData_SrvD%UseHSSBrake ) THEN
         IF ( p_FAST%CompAero == Module_AD14 ) THEN
            IF ( AD14%p%DYNINFL ) THEN
               CALL SetErrStat(ErrID_Fatal,'AeroDyn v14 "DYNINFL" InfModel is invalid for models with high-speed shaft braking.',ErrStat,ErrMsg,RoutineName)
            END IF
         END IF
         
            
         IF ( ED%p%method == 1 ) THEN ! bjj: should be using ElastoDyn's Method_ABM4 Method_AB4 parameters
            CALL SetErrStat(ErrID_Fatal,'ElastoDyn must use the AB4 or ABM4 integration method to implement high-speed shaft braking.',ErrStat,ErrMsg,RoutineName)
         END IF               
      END IF
      
   END IF  
   
   ! ........................
   ! some checks for AeroDyn14's Dynamic Inflow with Mean Wind Speed from InflowWind:
   ! (DO NOT COPY THIS CODE!)
   ! bjj: AeroDyn14 should not need this rule of thumb; it should check the instantaneous values when the code runs
   ! ........................   
   
   IF ( p_FAST%CompAero == Module_AD14 ) THEN
      IF (AD14%p%DynInfl) THEN               
         IF ( InitOutData_IfW%WindFileInfo%MWS  < 8.0 ) THEN
            CALL SetErrStat(ErrID_Fatal,'AeroDyn v14 "DYNINFL" InfModel is invalid for models with wind speeds less than 8 m/s.',ErrStat,ErrMsg,RoutineName)
            !CALL SetErrStat(ErrID_Info,'Estimated average inflow wind speed is less than 8 m/s. Dynamic Inflow will be turned off.',ErrStat,ErrMess,RoutineName )
         END IF
      END IF      
   END IF
   
   
   ! ........................
   ! initialize ServoDyn 
   ! ........................
   ALLOCATE( SrvD%Input( p_FAST%InterpOrder+1 ), SrvD%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating SrvD%Input and SrvD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
      
   IF ( p_FAST%CompServo == Module_SrvD ) THEN
      InitInData_SrvD%InputFile     = p_FAST%ServoFile
      InitInData_SrvD%RootName      = p_FAST%OutFileRoot
      InitInData_SrvD%NumBl         = InitOutData_ED%NumBl
      InitInData_SrvD%gravity       = InitOutData_ED%gravity
      InitInData_SrvD%r_N_O_G       = ED%Input(1)%NacelleLoads%Position(:,1) !bjj: check this!
      InitInData_SrvD%TMax          = p_FAST%TMax
      InitInData_SrvD%AirDens       = AirDens
      InitInData_SrvD%AvgWindSpeed  = InitOutData_IfW%WindFileInfo%MWS
      
      CALL AllocAry(InitInData_SrvD%BlPitchInit, InitOutData_ED%NumBl, 'BlPitchInit', ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      InitInData_SrvD%BlPitchInit   = InitOutData_ED%BlPitch
      CALL SrvD_Init( InitInData_SrvD, SrvD%Input(1), SrvD%p, SrvD%x(STATE_CURR), SrvD%xd(STATE_CURR), SrvD%z(STATE_CURR), &
                      SrvD%OtherSt, SrvD%y, p_FAST%dt_module( MODULE_SrvD ), InitOutData_SrvD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      p_FAST%ModuleInitialized(Module_SrvD) = .TRUE.

      !IF ( InitOutData_SrvD%CouplingScheme == ExplicitLoose ) THEN ...  bjj: abort if we're doing anything else!

      CALL SetModuleSubstepTime(Module_SrvD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      !! initialize y%ElecPwr and y%GenTq because they are one timestep different (used as input for the next step)
      !!bjj: perhaps this will require some better thought so that these two fields of y_SrvD_prev don't get set here in the glue code
      !CALL SrvD_CopyOutput( SrvD%y, SrvD%y_prev, MESH_NEWCOPY, ErrStat, ErrMsg)               
      !   
                  
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF          
   END IF

   
   
   ! ........................
   ! initialize HydroDyn 
   ! ........................
   ALLOCATE( HD%Input( p_FAST%InterpOrder+1 ), HD%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating HD%Input and HD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
      
   IF ( p_FAST%CompHydro == Module_HD ) THEN

      InitInData_HD%Gravity       = InitOutData_ED%Gravity
      InitInData_HD%UseInputFile  = .TRUE.
      InitInData_HD%InputFile     = p_FAST%HydroFile
      InitInData_HD%OutRootName   = p_FAST%OutFileRoot
      InitInData_HD%TMax          = p_FAST%TMax
      InitInData_HD%hasIce        = p_FAST%CompIce /= Module_None
      
         ! if wave field needs an offset, modify these values (added at request of SOWFA developers):
      InitInData_HD%PtfmLocationX = 0.0_ReKi  
      InitInData_HD%PtfmLocationY = 0.0_ReKi
      
      CALL HydroDyn_Init( InitInData_HD, HD%Input(1), HD%p,  HD%x(STATE_CURR), HD%xd(STATE_CURR), HD%z(STATE_CURR), HD%OtherSt, &
                          HD%y, p_FAST%dt_module( MODULE_HD ), InitOutData_HD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      p_FAST%ModuleInitialized(Module_HD) = .TRUE.
      CALL SetModuleSubstepTime(Module_HD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF       
   END IF   ! CompHydro

   ! ........................
   ! initialize SubDyn 
   ! ........................
   ALLOCATE( SD%Input( p_FAST%InterpOrder+1 ), SD%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating SD%Input and SD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF

   IF ( p_FAST%CompSub == Module_SD ) THEN
          
      IF ( p_FAST%CompHydro == Module_HD ) THEN
         InitInData_SD%WtrDpth = InitOutData_HD%WtrDpth
      ELSE
         InitInData_SD%WtrDpth = 0.0_ReKi
      END IF
            
      InitInData_SD%g             = InitOutData_ED%Gravity     
      !InitInData_SD%UseInputFile = .TRUE. 
      InitInData_SD%SDInputFile   = p_FAST%SubFile
      InitInData_SD%RootName      = p_FAST%OutFileRoot
      InitInData_SD%TP_RefPoint   = ED%Output(1)%PlatformPtMesh%Position(:,1)  ! bjj: not sure what this is supposed to be 
      InitInData_SD%SubRotateZ    = 0.0                                        ! bjj: not sure what this is supposed to be 
      
            
      CALL SD_Init( InitInData_SD, SD%Input(1), SD%p,  SD%x(STATE_CURR), SD%xd(STATE_CURR), SD%z(STATE_CURR), SD%OtherSt, &
                      SD%y, p_FAST%dt_module( MODULE_SD ), InitOutData_SD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      p_FAST%ModuleInitialized(Module_SD) = .TRUE.
      CALL SetModuleSubstepTime(Module_SD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF              
   END IF

   ! ------------------------------
   ! initialize CompMooring modules 
   ! ------------------------------
   ALLOCATE( MAPp%Input( p_FAST%InterpOrder+1 ), MAPp%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating MAPp%Input and MAPp%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
   ALLOCATE( MD%Input( p_FAST%InterpOrder+1 ), MD%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating MD%Input and MD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF   
   ALLOCATE( FEAM%Input( p_FAST%InterpOrder+1 ), FEAM%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating FEAM%Input and FEAM%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF   
      
   ! ........................
   ! initialize MAP 
   ! ........................
   IF (p_FAST%CompMooring == Module_MAP) THEN
      !bjj: until we modify this, MAP requires HydroDyn to be used. (perhaps we could send air density from AeroDyn or something...)
      
      CALL WrScr(NewLine) !bjj: I'm printing two blank lines here because MAP seems to be writing over the last line on the screen.
      

!      InitInData_MAP%rootname          =  p_FAST%OutFileRoot        ! Output file name 
      InitInData_MAP%gravity           =  InitOutData_ED%Gravity    ! This need to be according to g used in ElastoDyn
      InitInData_MAP%sea_density       =  InitOutData_HD%WtrDens    ! This needs to be set according to seawater density in HydroDyn
      InitInData_MAP%depth             =  InitOutData_HD%WtrDpth    ! This need to be set according to the water depth in HydroDyn
                  
   ! differences for MAP++
      InitInData_MAP%file_name         =  p_FAST%MooringFile        ! This needs to be set according to what is in the FAST input file. 
      InitInData_MAP%summary_file_name =  TRIM(p_FAST%OutFileRoot)//'.MAP.sum'        ! Output file name 
      InitInData_MAP%depth             = -InitOutData_HD%WtrDpth    ! This need to be set according to the water depth in HydroDyn
            

      
      CALL MAP_Init( InitInData_MAP, MAPp%Input(1), MAPp%p,  MAPp%x(STATE_CURR), MAPp%xd(STATE_CURR), MAPp%z(STATE_CURR), MAPp%OtherSt, &
                      MAPp%y, p_FAST%dt_module( MODULE_MAP ), InitOutData_MAP, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      p_FAST%ModuleInitialized(Module_MAP) = .TRUE.
      CALL SetModuleSubstepTime(Module_MAP, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
                  
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF              
   ! ........................
   ! initialize MoorDyn 
   ! ........................
   ELSEIF (p_FAST%CompMooring == Module_MD) THEN
                        
      InitInData_MD%FileName  = p_FAST%MooringFile         ! This needs to be set according to what is in the FAST input file. 
      InitInData_MD%RootName  = p_FAST%OutFileRoot
      
      InitInData_MD%PtfmInit  = InitOutData_ED%PlatformPos !ED%x(STATE_CURR)%QT(1:6)   ! initial position of the platform !bjj: this should come from InitOutData_ED, not x_ED
      InitInData_MD%g         = InitOutData_ED%Gravity     ! This need to be according to g used in ElastoDyn 
      InitInData_MD%rhoW      = InitOutData_HD%WtrDens     ! This needs to be set according to seawater density in HydroDyn      
      InitInData_MD%WtrDepth  = InitOutData_HD%WtrDpth    ! This need to be set according to the water depth in HydroDyn
            
      CALL MD_Init( InitInData_MD, MD%Input(1), MD%p,  MD%x(STATE_CURR), MD%xd(STATE_CURR), MD%z(STATE_CURR), MD%OtherSt, &
                      MD%y, p_FAST%dt_module( MODULE_MD ), InitOutData_MD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      p_FAST%ModuleInitialized(Module_MD) = .TRUE.
      CALL SetModuleSubstepTime(Module_MD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
   ! ........................
   ! initialize FEAM 
   ! ........................
   ELSEIF (p_FAST%CompMooring == Module_FEAM) THEN
            
      InitInData_FEAM%InputFile   = p_FAST%MooringFile         ! This needs to be set according to what is in the FAST input file. 
      InitInData_FEAM%RootName    = p_FAST%OutFileRoot
      
      InitInData_FEAM%PtfmInit    = InitOutData_ED%PlatformPos !ED%x(STATE_CURR)%QT(1:6)   ! initial position of the platform !bjj: this should come from InitOutData_ED, not x_ED
      InitInData_FEAM%NStepWave   = 1                          ! an arbitrary number > 0 (to set the size of the wave data, which currently contains all zero values)     
      InitInData_FEAM%gravity     = InitOutData_ED%Gravity     ! This need to be according to g used in ElastoDyn 
      InitInData_FEAM%WtrDens     = InitOutData_HD%WtrDens     ! This needs to be set according to seawater density in HydroDyn      
!      InitInData_FEAM%depth       =  InitOutData_HD%WtrDpth    ! This need to be set according to the water depth in HydroDyn
            
      CALL FEAM_Init( InitInData_FEAM, FEAM%Input(1), FEAM%p,  FEAM%x(STATE_CURR), FEAM%xd(STATE_CURR), FEAM%z(STATE_CURR), FEAM%OtherSt, &
                      FEAM%y, p_FAST%dt_module( MODULE_FEAM ), InitOutData_FEAM, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      p_FAST%ModuleInitialized(Module_FEAM) = .TRUE.
      CALL SetModuleSubstepTime(Module_FEAM, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF              
   END IF

   ! ------------------------------
   ! initialize CompIce modules 
   ! ------------------------------
   ALLOCATE( IceF%Input( p_FAST%InterpOrder+1 ), IceF%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating IceF%Input and IceF%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF  
      
      ! We need this to be allocated (else we have issues passing nonallocated arrays and using the first index of Input(),
      !   but we don't need the space of IceD_MaxLegs if we're not using it. 
   IF ( p_FAST%CompIce /= Module_IceD ) THEN   
      IceDim = 1
   ELSE
      IceDim = IceD_MaxLegs
   END IF
      
      ! because there may be multiple instances of IceDyn, we'll allocate arrays for that here
      ! we could allocate these after 
   ALLOCATE( IceD%Input( p_FAST%InterpOrder+1, IceDim ), IceD%InputTimes( p_FAST%InterpOrder+1, IceDim ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating IceD%Input and IceD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF  
      
     ALLOCATE( IceD%x(           IceDim,2), &
               IceD%xd(          IceDim,2), &
               IceD%z(           IceDim,2), &
               IceD%OtherSt(     IceDim  ), &
               IceD%p(           IceDim  ), &
               IceD%u(           IceDim  ), &
               IceD%y(           IceDim  ), &
               IceD%OtherSt_old( IceDim  ), &
                                             STAT = ErrStat2 )                                                  
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating IceD state, input, and output data.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF      
         
         
   ! ........................
   ! initialize IceFloe 
   ! ........................
   IF ( p_FAST%CompIce == Module_IceF ) THEN
                      
      InitInData_IceF%InputFile     = p_FAST%IceFile
      InitInData_IceF%RootName      = p_FAST%OutFileRoot     
      InitInData_IceF%simLength     = p_FAST%TMax  !bjj: IceFloe stores this as single-precision (ReKi) TMax is DbKi
      InitInData_IceF%MSL2SWL       = InitOutData_HD%MSL2SWL
      InitInData_IceF%gravity       = InitOutData_ED%Gravity
      
      CALL IceFloe_Init( InitInData_IceF, IceF%Input(1), IceF%p,  IceF%x(STATE_CURR), IceF%xd(STATE_CURR), IceF%z(STATE_CURR), IceF%OtherSt, &
                         IceF%y, p_FAST%dt_module( MODULE_IceF ), InitOutData_IceF, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      p_FAST%ModuleInitialized(Module_IceF) = .TRUE.
      CALL SetModuleSubstepTime(Module_IceF, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
              
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF              
   ! ........................
   ! initialize IceDyn 
   ! ........................
   ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN  
      
      InitInData_IceD%InputFile     = p_FAST%IceFile
      InitInData_IceD%RootName      = p_FAST%OutFileRoot     
      InitInData_IceD%MSL2SWL       = InitOutData_HD%MSL2SWL      
      InitInData_IceD%WtrDens       = InitOutData_HD%WtrDens    
      InitInData_IceD%gravity       = InitOutData_ED%Gravity
      InitInData_IceD%TMax          = p_FAST%TMax
      InitInData_IceD%LegNum        = 1
      
      CALL IceD_Init( InitInData_IceD, IceD%Input(1,1), IceD%p(1),  IceD%x(1,STATE_CURR), IceD%xd(1,STATE_CURR), IceD%z(1,STATE_CURR), &
                      IceD%OtherSt(1), IceD%y(1), p_FAST%dt_module( MODULE_IceD ), InitOutData_IceD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      p_FAST%ModuleInitialized(Module_IceD) = .TRUE.
      CALL SetModuleSubstepTime(Module_IceD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
         ! now initialize IceD for additional legs (if necessary)
      dt_IceD           = p_FAST%dt_module( MODULE_IceD )
      p_FAST%numIceLegs = InitOutData_IceD%numLegs     
      
      IF (p_FAST%numIceLegs > IceD_MaxLegs) THEN
         CALL SetErrStat(ErrID_Fatal,'IceDyn-FAST coupling is supported for up to '//TRIM(Num2LStr(IceD_MaxLegs))//' legs, but ' &
                           //TRIM(Num2LStr(p_FAST%numIceLegs))//' legs were specified.',ErrStat,ErrMsg,RoutineName)
      END IF
                  

      DO i=2,p_FAST%numIceLegs  ! basically, we just need IceDyn to set up its meshes for inputs/outputs and possibly initial values for states
         InitInData_IceD%LegNum = i
         
         CALL IceD_Init( InitInData_IceD, IceD%Input(1,i), IceD%p(i),  IceD%x(i,STATE_CURR), IceD%xd(i,STATE_CURR), IceD%z(i,STATE_CURR), &
                            IceD%OtherSt(i), IceD%y(i), dt_IceD, InitOutData_IceD, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            
         !bjj: we're going to force this to have the same timestep because I don't want to have to deal with n IceD modules with n timesteps.
         IF (.NOT. EqualRealNos( p_FAST%dt_module( MODULE_IceD ),dt_IceD )) THEN
            CALL SetErrStat(ErrID_Fatal,"All instances of IceDyn (one per support-structure leg) must be the same",ErrStat,ErrMsg,RoutineName)
         END IF
      END DO
            
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF           
      
   END IF   
   

   ! ........................
   ! Set up output for glue code (must be done after all modules are initialized so we have their WriteOutput information)
   ! ........................

   CALL FAST_InitOutput( p_FAST, y_FAST, InitOutData_ED, InitOutData_BD, InitOutData_SrvD, InitOutData_AD14, InitOutData_AD, &
                         InitOutData_IfW, InitOutData_OpFM, InitOutData_HD, InitOutData_SD, &
                         InitOutData_MAP, InitOutData_FEAM, InitOutData_MD, InitOutData_IceF, InitOutData_IceD, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)


   ! -------------------------------------------------------------------------
   ! Initialize mesh-mapping data
   ! -------------------------------------------------------------------------

   CALL InitModuleMappings(p_FAST, ED, BD, AD14, AD, HD, SD, SrvD, MAPp, FEAM, MD, IceF, IceD, MeshMapData, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF      
      
      
   ! -------------------------------------------------------------------------
   ! Write initialization data to FAST summary file:
   ! -------------------------------------------------------------------------
   
   CALL FAST_WrSum( p_FAST, y_FAST, MeshMapData, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
   
   ! -------------------------------------------------------------------------
   ! other misc variables initialized here:
   ! -------------------------------------------------------------------------
      
   m_FAST%t_global   = t_initial
         
   ! Initialize external inputs for first step  
   if ( p_FAST%CompServo == MODULE_SrvD ) then      
      m_FAST%ExternInput%GenTrq     = SrvD%Input(1)%ExternalGenTrq !0.0_ReKi
      m_FAST%ExternInput%ElecPwr    = SrvD%Input(1)%ExternalElecPwr
      m_FAST%ExternInput%YawPosCom  = SrvD%Input(1)%ExternalYawPosCom
      m_FAST%ExternInput%YawRateCom = SrvD%Input(1)%ExternalYawRateCom
      m_FAST%ExternInput%HSSBrFrac  = SrvD%Input(1)%ExternalHSSBrFrac
      
      do i=1,SIZE(SrvD%Input(1)%ExternalBlPitchCom)
         m_FAST%ExternInput%BlPitchCom(i) = SrvD%Input(1)%ExternalBlPitchCom(i)
      end do   
   end if
   
   m_FAST%ExternInput%LidarFocus = 1.0_ReKi  ! make this non-zero (until we add the initial position in the InflowWind input file)
         
   
   !...............................................................................................................................
   ! Destroy initializion data
   !...............................................................................................................................      
   CALL Cleanup()
   
CONTAINS
   SUBROUTINE Cleanup()
   !...............................................................................................................................
   ! Destroy initializion data
   !...............................................................................................................................
   
      CALL ED_DestroyInitInput(  InitInData_ED,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL ED_DestroyInitOutput( InitOutData_ED, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      CALL BD_DestroyInitInput(  InitInData_BD,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)         
      IF (ALLOCATED(InitOutData_BD)) THEN
         DO i=1,p_FAST%nBeams
            CALL BD_DestroyInitOutput( InitOutData_BD(i), ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         END DO
         DEALLOCATE(InitOutData_BD)
      END IF
                     
      CALL AD14_DestroyInitInput(  InitInData_AD14,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL AD14_DestroyInitOutput( InitOutData_AD14, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      CALL AD_DestroyInitInput(  InitInData_AD,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL AD_DestroyInitOutput( InitOutData_AD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
      CALL InflowWind_DestroyInitInput(  InitInData_IfW,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL InflowWind_DestroyInitOutput( InitOutData_IfW, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
      CALL OpFM_DestroyInitOutput( InitOutData_OpFM, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)         
         
      CALL SrvD_DestroyInitInput(  InitInData_SrvD,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL SrvD_DestroyInitOutput( InitOutData_SrvD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      CALL HydroDyn_DestroyInitInput(  InitInData_HD,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL HydroDyn_DestroyInitOutput( InitOutData_HD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      CALL SD_DestroyInitInput(  InitInData_SD,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL SD_DestroyInitOutput( InitOutData_SD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      CALL MAP_DestroyInitInput(  InitInData_MAP,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL MAP_DestroyInitOutput( InitOutData_MAP, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
      CALL FEAM_DestroyInitInput(  InitInData_FEAM,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL FEAM_DestroyInitOutput( InitOutData_FEAM, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      CALL MD_DestroyInitInput(  InitInData_MD,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL MD_DestroyInitOutput( InitOutData_MD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
                  
      CALL IceFloe_DestroyInitInput(  InitInData_IceF,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL IceFloe_DestroyInitOutput( InitOutData_IceF, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
      CALL IceD_DestroyInitInput(  InitInData_IceD,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL IceD_DestroyInitOutput( InitOutData_IceD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
   
   END SUBROUTINE Cleanup

END SUBROUTINE FAST_InitializeAll
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FAST_ExtrapInterpMods( t_global_next, p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, HD, SD, MAPp, FEAM, MD, IceF, IceD, ErrStat, ErrMsg )

   REAL(DbKi),               INTENT(IN   ) :: t_global_next       ! next global time step (t + dt), at which we're extrapolating inputs (and ED outputs)
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              ! Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(IN   ) :: y_FAST              ! Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(IN   ) :: m_FAST              ! Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  ! ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  ! BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                ! ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                ! AeroDyn14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  ! AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 ! InflowWind data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  ! HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  ! SubDyn data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                ! MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                ! FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  ! Data for the MoorDyn module
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                ! IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                ! All the IceDyn data used in time-step loop

   !TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         ! Data for mapping between modules
      
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              ! Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi)                          :: i, j, k             ! loop counters
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_ExtrapInterpMods'       
   
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Step 1.a: Extrapolate Inputs (gives predicted values at t+dt)
      ! 
      ! a) Extrapolate inputs (and outputs -- bjj: output extrapolation not necessary, yet) 
      !    to t + dt (i.e., t_global_next); will only be used by modules with an implicit dependence on input data.
      ! b) Shift "window" of the ModName_Input and ModName_Output
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
      ErrStat = ErrID_None
      ErrMsg  = ""
      
      ! ElastoDyn
      CALL ED_Input_ExtrapInterp(ED%Input, ED%InputTimes, ED%u, t_global_next, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
  
      CALL ED_Output_ExtrapInterp(ED%Output, ED%InputTimes, ED%y, t_global_next, ErrStat2, ErrMsg2) !this extrapolated value is used in the ED-HD coupling
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         
         
      DO j = p_FAST%InterpOrder, 1, -1
         CALL ED_CopyInput (ED%Input(j),  ED%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         CALL ED_CopyOutput(ED%Output(j), ED%Output(j+1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         ED%InputTimes(j+1) = ED%InputTimes(j)
         !ED_OutputTimes(j+1) = ED_OutputTimes(j)
      END DO
  
      CALL ED_CopyInput (ED%u,  ED%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
      CALL ED_CopyOutput (ED%y,  ED%Output(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
      ED%InputTimes(1)  = t_global_next
      !ED_OutputTimes(1) = t_global_next 
  
      
      ! BeamDyn
      IF (p_FAST%CompElast == Module_BD) THEN
         
         DO k = 1,p_FAST%nBeams
         
            CALL BD_Input_ExtrapInterp(BD%Input(:,k), BD%InputTimes(:,k), BD%u(k), t_global_next, ErrStat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
                        
            
            ! Shift "window" of BD%Input 
  
            DO j = p_FAST%InterpOrder, 1, -1
               CALL BD_CopyInput (BD%Input(j,k),  BD%Input(j+1,k),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
                  CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
               BD%InputTimes(j+1,k) = BD%InputTimes(j,k)
            END DO
  
            CALL BD_CopyInput (BD%u(k),  BD%Input(1,k),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            BD%InputTimes(1,k) = t_global_next          
            
         END DO ! k=p_FAST%nBeams
         
      END IF  ! BeamDyn      
      
      
      ! AeroDyn v14
      IF ( p_FAST%CompAero == Module_AD14 ) THEN
         
         CALL AD14_Input_ExtrapInterp(AD14%Input, AD14%InputTimes, AD14%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
         !CALL AD14_Output_ExtrapInterp(AD14_Output, AD14_OutputTimes, AD14%y, t_global_next, ErrStat2, ErrMsg2)
         !   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
            
         ! Shift "window" of AD14%Input and AD14_Output
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL AD14_CopyInput (AD14%Input(j),  AD14%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            AD14%InputTimes(j+1)  = AD14%InputTimes(j)
         END DO
  
         CALL AD14_CopyInput (AD14%u,  AD14%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         AD14%InputTimes(1)  = t_global_next          
            
      ELSEIF ( p_FAST%CompAero == Module_AD ) THEN
         
         CALL AD_Input_ExtrapInterp(AD%Input, AD%InputTimes, AD%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
                        
         ! Shift "window" of AD%Input 
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL AD_CopyInput (AD%Input(j),  AD%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            AD%InputTimes(j+1)  = AD%InputTimes(j)
         END DO
  
         CALL AD_CopyInput (AD%u,  AD%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         AD%InputTimes(1)  = t_global_next    
         
      END IF  ! CompAero      
      
         
      ! InflowWind
      IF ( p_FAST%CompInflow == Module_IfW ) THEN
         
         CALL InflowWind_Input_ExtrapInterp(IfW%Input, IfW%InputTimes, IfW%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
         !CALL InflowWind_Output_ExtrapInterp(IfW_Output, IfW_OutputTimes, IfW%y, t_global_next, ErrStat2, ErrMsg2)
         !   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
            
         ! Shift "window" of IfW%Input and IfW_Output
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL InflowWind_CopyInput (IfW%Input(j),  IfW%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
           !CALL InflowWind_CopyOutput(IfW_Output(j), IfW_Output(j+1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            IfW%InputTimes(j+1)  = IfW%InputTimes(j)
           !IfW_OutputTimes(j+1) = IfW_OutputTimes(j)
         END DO
  
         CALL InflowWind_CopyInput (IfW%u,  IfW%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
        !CALL InflowWind_CopyOutput(IfW%y,  IfW_Output(1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         IfW%InputTimes(1)  = t_global_next          
        !IfW_OutputTimes(1) = t_global_next 
            
      END IF  ! CompInflow          
      
      
      ! ServoDyn
      IF ( p_FAST%CompServo == Module_SrvD ) THEN
         
         CALL SrvD_Input_ExtrapInterp(SrvD%Input, SrvD%InputTimes, SrvD%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
                  
         !CALL SrvD_Output_ExtrapInterp(SrvD_Output, SrvD_OutputTimes, SrvD%y, t_global_next, ErrStat2, ErrMsg2)
         !   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
            
         ! Shift "window" of SrvD%Input and SrvD_Output
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL SrvD_CopyInput (SrvD%Input(j),  SrvD%Input(j+1),  MESH_UPDATECOPY, ErrStat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
           !CALL SrvD_CopyOutput(SrvD_Output(j), SrvD_Output(j+1), MESH_UPDATECOPY, ErrStat2, ErrMsg2)
            SrvD%InputTimes(j+1)  = SrvD%InputTimes(j)
           !SrvD_OutputTimes(j+1) = SrvD_OutputTimes(j)
         END DO
  
         CALL SrvD_CopyInput (SrvD%u,  SrvD%Input(1),  MESH_UPDATECOPY, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
        !CALL SrvD_CopyOutput(SrvD%y,  SrvD_Output(1), MESH_UPDATECOPY, ErrStat2, ErrMsg2)
         SrvD%InputTimes(1)  = t_global_next          
        !SrvD_OutputTimes(1) = t_global_next 
            
      END IF  ! ServoDyn       
      
      ! HydroDyn
      IF ( p_FAST%CompHydro == Module_HD ) THEN

         CALL HydroDyn_Input_ExtrapInterp(HD%Input, HD%InputTimes, HD%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
         !CALL HydroDyn_Output_ExtrapInterp(HD_Output, HD_OutputTimes, HD%y, t_global_next, ErrStat2, ErrMsg2)
         !   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
         ! Shift "window" of HD%Input and HD_Output
            
         DO j = p_FAST%InterpOrder, 1, -1

            CALL HydroDyn_CopyInput (HD%Input(j),  HD%Input(j+1),  MESH_UPDATECOPY, ErrStat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            !CALL HydroDyn_CopyOutput(HD_Output(j), HD_Output(j+1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            HD%InputTimes(j+1) = HD%InputTimes(j)
            !HD_OutputTimes(j+1)= HD_OutputTimes(j)
         END DO

         CALL HydroDyn_CopyInput (HD%u,  HD%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         !CALL HydroDyn_CopyOutput(HD%y,  HD_Output(1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         HD%InputTimes(1) = t_global_next          
         !HD_OutputTimes(1) = t_global_next
            
      END IF  ! HydroDyn

      
      ! SubDyn
      IF ( p_FAST%CompSub == Module_SD ) THEN

         CALL SD_Input_ExtrapInterp(SD%Input, SD%InputTimes, SD%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
                        
         !CALL SD_Output_ExtrapInterp(SD_Output, SD_OutputTimes, SD%y, t_global_next, ErrStat2, ErrMsg2)
         !   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
            
         ! Shift "window" of SD%Input and SD_Output
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL SD_CopyInput (SD%Input(j),  SD%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
           !CALL SD_CopyOutput(SD_Output(j), SD_Output(j+1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            SD%InputTimes(j+1) = SD%InputTimes(j)
            !SD_OutputTimes(j+1) = SD_OutputTimes(j)
         END DO
  
         CALL SD_CopyInput (SD%u,  SD%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         !CALL SD_CopyOutput(SD%y,  SD_Output(1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         SD%InputTimes(1) = t_global_next          
         !SD_OutputTimes(1) = t_global_next 
            
      END IF  ! SubDyn
      
      
      ! Mooring (MAP , FEAM , MoorDyn)
      ! MAP
      IF ( p_FAST%CompMooring == Module_MAP ) THEN
         
         CALL MAP_Input_ExtrapInterp(MAPp%Input, MAPp%InputTimes, MAPp%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
         !CALL MAP_Output_ExtrapInterp(MAP_Output, MAP_OutputTimes, MAPp%y, t_global_next, ErrStat2, ErrMsg2)
         !   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
            
         ! Shift "window" of MAPp%Input and MAP_Output
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL MAP_CopyInput (MAPp%Input(j),  MAPp%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
           !CALL MAP_CopyOutput(MAP_Output(j), MAP_Output(j+1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            MAPp%InputTimes(j+1) = MAPp%InputTimes(j)
            !MAP_OutputTimes(j+1) = MAP_OutputTimes(j)
         END DO
  
         CALL MAP_CopyInput (MAPp%u,  MAPp%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         !CALL MAP_CopyOutput(MAPp%y,  MAP_Output(1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         MAPp%InputTimes(1) = t_global_next          
         !MAP_OutputTimes(1) = t_global_next 
            
      ! MoorDyn
      ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
         
         CALL MD_Input_ExtrapInterp(MD%Input, MD%InputTimes, MD%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
         !CALL MD_Output_ExtrapInterp(MD_Output, MD_OutputTimes, MD%y, t_global_next, ErrStat2, ErrMsg2)
         !   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
            
         ! Shift "window" of MD%Input and MD_Output
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL MD_CopyInput (MD%Input(j),  MD%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
           !CALL MD_CopyOutput(MD_Output(j), MD_Output(j+1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            MD%InputTimes( j+1) = MD%InputTimes( j)
           !MD_OutputTimes(j+1) = MD_OutputTimes(j)
         END DO
  
         CALL MD_CopyInput (MD%u,  MD%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
        !CALL MD_CopyOutput(MD%y,  MD_Output(1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         MD%InputTimes(1)  = t_global_next          
        !MD_OutputTimes(1) = t_global_next 
         
      ! FEAM
      ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
         
         CALL FEAM_Input_ExtrapInterp(FEAM%Input, FEAM%InputTimes, FEAM%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
         !CALL FEAM_Output_ExtrapInterp(FEAM_Output, FEAM_OutputTimes, FEAM%y, t_global_next, ErrStat2, ErrMsg2)
         !   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
            
         ! Shift "window" of FEAM%Input and FEAM_Output
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL FEAM_CopyInput (FEAM%Input(j),  FEAM%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
           !CALL FEAM_CopyOutput(FEAM_Output(j), FEAM_Output(j+1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            FEAM%InputTimes( j+1) = FEAM%InputTimes( j)
           !FEAM_OutputTimes(j+1) = FEAM_OutputTimes(j)
         END DO
  
         CALL FEAM_CopyInput (FEAM%u,  FEAM%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
        !CALL FEAM_CopyOutput(FEAM%y,  FEAM_Output(1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         FEAM%InputTimes(1)  = t_global_next          
        !FEAM_OutputTimes(1) = t_global_next 
         
      END IF  ! MAP/FEAM/MoorDyn
           
            
      ! Ice (IceFloe or IceDyn)
      ! IceFloe
      IF ( p_FAST%CompIce == Module_IceF ) THEN
         
         CALL IceFloe_Input_ExtrapInterp(IceF%Input, IceF%InputTimes, IceF%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
                        
         !CALL IceFloe_Output_ExtrapInterp(IceF_Output, IceF_OutputTimes, IceF%y, t_global_next, ErrStat2, ErrMsg2)
         !   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
            
         ! Shift "window" of IceF%Input and IceF_Output
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL IceFloe_CopyInput (IceF%Input(j),  IceF%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
           !CALL IceFloe_CopyOutput(IceF_Output(j), IceF_Output(j+1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            IceF%InputTimes(j+1) = IceF%InputTimes(j)
            !IceF_OutputTimes(j+1) = IceF_OutputTimes(j)
         END DO
  
         CALL IceFloe_CopyInput (IceF%u,  IceF%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         !CALL IceFloe_CopyOutput(IceF%y,  IceF_Output(1), MESH_UPDATECOPY, Errstat, ErrMsg)
         IceF%InputTimes(1) = t_global_next          
         !IceF_OutputTimes(1) = t_global_next 
            
      ! IceDyn
      ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
         
         DO i = 1,p_FAST%numIceLegs
         
            CALL IceD_Input_ExtrapInterp(IceD%Input(:,i), IceD%InputTimes(:,i), IceD%u(i), t_global_next, ErrStat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
                        
            !CALL IceD_Output_ExtrapInterp(IceD%Output(:,i), IceD%OutputTimes(:,i), IceD%y(i), t_global_next, ErrStat2, ErrMsg2)
            !   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
            
            ! Shift "window" of IceD%Input and IceD%Output
  
            DO j = p_FAST%InterpOrder, 1, -1
               CALL IceD_CopyInput (IceD%Input(j,i),  IceD%Input(j+1,i),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
                  CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
              !CALL IceD_CopyOutput(IceD%Output(j,i), IceD%Output(j+1,i), MESH_UPDATECOPY, Errstat2, ErrMsg2)
               IceD%InputTimes(j+1,i) = IceD%InputTimes(j,i)
              !IceD%OutputTimes(j+1,i) = IceD%OutputTimes(j,i)
            END DO
  
            CALL IceD_CopyInput (IceD%u(i),  IceD%Input(1,i),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
           !CALL IceD_CopyOutput(IceD%y(i),  IceD%Output(1,i), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            IceD%InputTimes(1,i) = t_global_next          
           !IceD%OutputTimes(1,i) = t_global_next 
            
         END DO ! numIceLegs
         
      
      END IF  ! IceFloe/IceDyn



END SUBROUTINE FAST_ExtrapInterpMods
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FAST_InitIOarrays( t_initial, p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, HD, SD, MAPp, FEAM, MD, IceF, IceD, ErrStat, ErrMsg )

   REAL(DbKi),               INTENT(IN   ) :: t_initial           ! start time of the simulation 
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              ! Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(IN   ) :: y_FAST              ! Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(IN   ) :: m_FAST              ! Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  ! ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  ! BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                ! ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                ! AeroDyn v14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  ! AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 ! InflowWind data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  ! HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  ! SubDyn data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                ! MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                ! FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  ! MoorDyn data
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                ! IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                ! All the IceDyn data used in time-step loop
      
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              ! Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi)                          :: i, j, k             ! loop counters
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_InitIOarrays'       
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   ! We fill ED%InputTimes with negative times, but the ED%Input values are identical for each of those times; this allows
   ! us to use, e.g., quadratic interpolation that effectively acts as a zeroth-order extrapolation and first-order extrapolation
   ! for the first and second time steps.  (The interpolation order in the ExtrapInput routines are determined as
   ! order = SIZE(ED%Input)

   DO j = 1, p_FAST%InterpOrder + 1
      ED%InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
      !ED_OutputTimes(j) = t_initial - (j - 1) * dt
   END DO

   DO j = 2, p_FAST%InterpOrder + 1
      CALL ED_CopyInput (ED%Input(1),  ED%Input(j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      CALL ED_CopyOutput (ED%Output(1), ED%Output(j), MESH_NEWCOPY, Errstat2, ErrMsg2) !BJJ: THIS IS REALLY ONLY NECESSARY FOR ED-HD COUPLING AT THE MOMENT
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END DO
   CALL ED_CopyInput (ED%Input(1),  ED%u,  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL ED_CopyOutput (ED%Output(1), ED%y, MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
   
      
      ! Initialize predicted states for j_pc loop:
   CALL ED_CopyContState   (ED%x( STATE_CURR), ED%x( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL ED_CopyDiscState   (ED%xd(STATE_CURR), ED%xd(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)  
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL ED_CopyConstrState (ED%z( STATE_CURR), ED%z( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
   IF ( p_FAST%n_substeps( MODULE_ED ) > 1 ) THEN
      CALL ED_CopyOtherState( ED%OtherSt, ED%OtherSt_old, MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
   END IF   
      
   
   IF  (p_FAST%CompElast == Module_BD ) THEN      

      DO k = 1,p_FAST%nBeams
         
            ! Copy values for interpolation/extrapolation:
         DO j = 1, p_FAST%InterpOrder + 1
            BD%InputTimes(j,k) = t_initial - (j - 1) * p_FAST%dt
         END DO

         DO j = 2, p_FAST%InterpOrder + 1
            CALL BD_CopyInput (BD%Input(1,k),  BD%Input(j,k),  MESH_NEWCOPY, Errstat2, ErrMsg2)
               CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         END DO
         CALL BD_CopyInput (BD%Input(1,k),  BD%u(k),  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )      
                               
                     
            ! Initialize predicted states for j_pc loop:
         CALL BD_CopyContState   (BD%x( k,STATE_CURR), BD%x( k,STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL BD_CopyDiscState   (BD%xd(k,STATE_CURR), BD%xd(k,STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)  
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL BD_CopyConstrState (BD%z( k,STATE_CURR), BD%z( k,STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL BD_CopyOtherState (BD%OtherSt( k,STATE_CURR), BD%OtherSt( k,STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
      END DO ! nBeams
      
   END IF ! CompElast            
      
   
   IF ( p_FAST%CompServo == Module_SrvD ) THEN      
      ! Initialize Input-Output arrays for interpolation/extrapolation:
         
      DO j = 1, p_FAST%InterpOrder + 1
         SrvD%InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
         !SrvD_OutputTimes(j) = t_initial - (j - 1) * dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL SrvD_CopyInput (SrvD%Input(1),  SrvD%Input(j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
      CALL SrvD_CopyInput (SrvD%Input(1),  SrvD%u,  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
         ! Initialize predicted states for j_pc loop:
      CALL SrvD_CopyContState   (SrvD%x( STATE_CURR), SrvD%x( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SrvD_CopyDiscState   (SrvD%xd(STATE_CURR), SrvD%xd(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SrvD_CopyConstrState (SrvD%z( STATE_CURR), SrvD%z( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( p_FAST%n_substeps( MODULE_SrvD ) > 1 ) THEN
         CALL SrvD_CopyOtherState( SrvD%OtherSt, SrvD%OtherSt_old, MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
      END IF    
         
   END IF ! CompServo
   
   
   IF ( p_FAST%CompAero == Module_AD14 ) THEN      
         ! Copy values for interpolation/extrapolation:

      DO j = 1, p_FAST%InterpOrder + 1
         AD14%InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL AD14_CopyInput (AD14%Input(1),  AD14%Input(j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
      CALL AD14_CopyInput (AD14%Input(1),  AD14%u,  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


         ! Initialize predicted states for j_pc loop:
      CALL AD14_CopyContState   (AD14%x( STATE_CURR), AD14%x( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD14_CopyDiscState   (AD14%xd(STATE_CURR), AD14%xd(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD14_CopyConstrState (AD14%z( STATE_CURR), AD14%z( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )      
      IF ( p_FAST%n_substeps( MODULE_AD14 ) > 1 ) THEN
         CALL AD14_CopyOtherState( AD14%OtherSt, AD14%OtherSt_old, MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
      END IF         

   ELSEIF ( p_FAST%CompAero == Module_AD ) THEN      
         ! Copy values for interpolation/extrapolation:
   
      DO j = 1, p_FAST%InterpOrder + 1
         AD%InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
      END DO
   
      DO j = 2, p_FAST%InterpOrder + 1
         CALL AD_CopyInput (AD%Input(1),  AD%Input(j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
      CALL AD_CopyInput (AD%Input(1),  AD%u,  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
   
         ! Initialize predicted states for j_pc loop:
      CALL AD_CopyContState   (AD%x( STATE_CURR), AD%x( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD_CopyDiscState   (AD%xd(STATE_CURR), AD%xd(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD_CopyConstrState (AD%z( STATE_CURR), AD%z( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )      
      IF ( p_FAST%n_substeps( MODULE_AD ) > 1 ) THEN
         CALL AD_CopyOtherState( AD%OtherSt, AD%OtherSt_old, MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
      END IF         
            
   END IF ! CompAero == Module_AD    
   
   
   
   IF ( p_FAST%CompInflow == Module_IfW ) THEN      
         ! Copy values for interpolation/extrapolation:

      DO j = 1, p_FAST%InterpOrder + 1
         IfW%InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
         !IfW%OutputTimes(i) = t_initial - (j - 1) * dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL InflowWind_CopyInput (IfW%Input(1),  IfW%Input(j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
      CALL InflowWind_CopyInput (IfW%Input(1),  IfW%u,  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


         ! Initialize predicted states for j_pc loop:
      CALL InflowWind_CopyContState   (IfW%x( STATE_CURR), IfW%x( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL InflowWind_CopyDiscState   (IfW%xd(STATE_CURR), IfW%xd(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL InflowWind_CopyConstrState (IfW%z( STATE_CURR), IfW%z( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )      
      IF ( p_FAST%n_substeps( MODULE_IfW ) > 1 ) THEN
         CALL InflowWind_CopyOtherState( IfW%OtherSt, IfW%OtherSt_old, MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
      END IF         

   END IF ! CompInflow == Module_IfW 
      
   
   IF ( p_FAST%CompHydro == Module_HD ) THEN      
         ! Copy values for interpolation/extrapolation:
      DO j = 1, p_FAST%InterpOrder + 1
         HD%InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
         !HD_OutputTimes(i) = t_initial - (j - 1) * dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL HydroDyn_CopyInput (HD%Input(1),  HD%Input(j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
      CALL HydroDyn_CopyInput (HD%Input(1),  HD%u,  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


         ! Initialize predicted states for j_pc loop:
      CALL HydroDyn_CopyContState   (HD%x( STATE_CURR), HD%x( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL HydroDyn_CopyDiscState   (HD%xd(STATE_CURR), HD%xd(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL HydroDyn_CopyConstrState (HD%z( STATE_CURR), HD%z( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( p_FAST%n_substeps( MODULE_HD ) > 1 ) THEN
         CALL HydroDyn_CopyOtherState( HD%OtherSt, HD%OtherSt_old, MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
      END IF          
      
   END IF !CompHydro
         
   
   IF  (p_FAST%CompSub == Module_SD ) THEN      

         ! Copy values for interpolation/extrapolation:
      DO j = 1, p_FAST%InterpOrder + 1
         SD%InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
         !SD_OutputTimes(i) = t_initial - (j - 1) * dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL SD_CopyInput (SD%Input(1),  SD%Input(j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
      CALL SD_CopyInput (SD%Input(1),  SD%u,  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )      
                               
         
         ! Initialize predicted states for j_pc loop:
      CALL SD_CopyContState   (SD%x( STATE_CURR), SD%x( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SD_CopyDiscState   (SD%xd(STATE_CURR), SD%xd(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SD_CopyConstrState (SD%z( STATE_CURR), SD%z( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( p_FAST%n_substeps( MODULE_SD ) > 1 ) THEN
         CALL SD_CopyOtherState( SD%OtherSt_old, SD%OtherSt, MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
      END IF       
   END IF ! CompSub         
      
   
   IF (p_FAST%CompMooring == Module_MAP) THEN      
         ! Copy values for interpolation/extrapolation:

      DO j = 1, p_FAST%InterpOrder + 1
         MAPp%InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
         !MAP_OutputTimes(i) = t_initial - (j - 1) * dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL MAP_CopyInput (MAPp%Input(1),  MAPp%Input(j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
      CALL MAP_CopyInput (MAPp%Input(1),  MAPp%u,  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               
         ! Initialize predicted states for j_pc loop:
      CALL MAP_CopyContState   (MAPp%x( STATE_CURR), MAPp%x( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MAP_CopyDiscState   (MAPp%xd(STATE_CURR), MAPp%xd(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MAP_CopyConstrState (MAPp%z( STATE_CURR), MAPp%z( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( p_FAST%n_substeps( MODULE_MAP ) > 1 ) THEN
         CALL MAP_CopyOtherState( MAPp%OtherSt, MAPp%OtherSt_old, MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
      END IF  
      
   ELSEIF (p_FAST%CompMooring == Module_MD) THEN      
         ! Copy values for interpolation/extrapolation:

      DO j = 1, p_FAST%InterpOrder + 1
         MD%InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
         !MD_OutputTimes(i) = t_initial - (j - 1) * dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL MD_CopyInput (MD%Input(1),  MD%Input(j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
      CALL MD_CopyInput (MD%Input(1),  MD%u,  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               
         ! Initialize predicted states for j_pc loop:
      CALL MD_CopyContState   (MD%x( STATE_CURR), MD%x( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MD_CopyDiscState   (MD%xd(STATE_CURR), MD%xd(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MD_CopyConstrState (MD%z( STATE_CURR), MD%z( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( p_FAST%n_substeps( MODULE_MD ) > 1 ) THEN
         CALL MD_CopyOtherState( MD%OtherSt, MD%OtherSt_old, MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
      END IF        
      
   ELSEIF (p_FAST%CompMooring == Module_FEAM) THEN      
         ! Copy values for interpolation/extrapolation:

      DO j = 1, p_FAST%InterpOrder + 1
         FEAM%InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
         !FEAM_OutputTimes(i) = t_initial - (j - 1) * dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL FEAM_CopyInput (FEAM%Input(1),  FEAM%Input(j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
      CALL FEAM_CopyInput (FEAM%Input(1),  FEAM%u,  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               
         ! Initialize predicted states for j_pc loop:
      CALL FEAM_CopyContState   (FEAM%x( STATE_CURR), FEAM%x( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL FEAM_CopyDiscState   (FEAM%xd(STATE_CURR), FEAM%xd(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL FEAM_CopyConstrState (FEAM%z( STATE_CURR), FEAM%z( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( p_FAST%n_substeps( MODULE_FEAM ) > 1 ) THEN
         CALL FEAM_CopyOtherState( FEAM%OtherSt, FEAM%OtherSt_old, MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
      END IF           
   END IF ! CompMooring
                 
   
   IF  (p_FAST%CompIce == Module_IceF ) THEN      

         ! Copy values for interpolation/extrapolation:
      DO j = 1, p_FAST%InterpOrder + 1
         IceF%InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
         !IceF_OutputTimes(i) = t_initial - (j - 1) * dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL IceFloe_CopyInput (IceF%Input(1),  IceF%Input(j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
      CALL IceFloe_CopyInput (IceF%Input(1),  IceF%u,  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )      
                               
         
         ! Initialize predicted states for j_pc loop:
      CALL IceFloe_CopyContState   (IceF%x( STATE_CURR), IceF%x( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL IceFloe_CopyDiscState   (IceF%xd(STATE_CURR), IceF%xd(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL IceFloe_CopyConstrState (IceF%z( STATE_CURR), IceF%z( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( p_FAST%n_substeps( MODULE_IceF ) > 1 ) THEN
         CALL IceFloe_CopyOtherState( IceF%OtherSt_old, IceF%OtherSt, MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
      END IF       
      
   ELSEIF  (p_FAST%CompIce == Module_IceD ) THEN      

      DO i = 1,p_FAST%numIceLegs
         
            ! Copy values for interpolation/extrapolation:
         DO j = 1, p_FAST%InterpOrder + 1
            IceD%InputTimes(j,i) = t_initial - (j - 1) * p_FAST%dt
            !IceD%OutputTimes(j,i) = t_initial - (j - 1) * dt
         END DO

         DO j = 2, p_FAST%InterpOrder + 1
            CALL IceD_CopyInput (IceD%Input(1,i),  IceD%Input(j,i),  MESH_NEWCOPY, Errstat2, ErrMsg2)
               CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         END DO
         CALL IceD_CopyInput (IceD%Input(1,i),  IceD%u(i),  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )      
                               
         
            ! Initialize predicted states for j_pc loop:
         CALL IceD_CopyContState   (IceD%x( i,STATE_CURR), IceD%x( i,STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL IceD_CopyDiscState   (IceD%xd(i,STATE_CURR), IceD%xd(i,STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)  
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL IceD_CopyConstrState (IceD%z( i,STATE_CURR), IceD%z( i,STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF ( p_FAST%n_substeps( MODULE_IceD ) > 1 ) THEN
            CALL IceD_CopyOtherState( IceD%OtherSt_old(i), IceD%OtherSt(i), MESH_NEWCOPY, Errstat2, ErrMsg2)
               CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
         END IF       
         
      END DO ! numIceLegs
      
   END IF ! CompIce            
   
   
!bjj: TODO: FIX ME:   
      ! ServoDyn: copy current outputs to store as previous outputs for next step
      ! note that this is a violation of the framework as this is basically a state, but it's only used for the
      ! GH-Bladed DLL, which itself violates the framework....
   CALL SrvD_CopyOutput ( SrvD%y, SrvD%y_prev, MESH_NEWCOPY, ErrStat2, ErrMsg2)   
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
   
END SUBROUTINE FAST_InitIOarrays
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FAST_AdvanceStates( t_initial, n_t_global, p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, HD, SD, MAPp, FEAM, MD, IceF, IceD, ErrStat, ErrMsg )

   REAL(DbKi),               INTENT(IN   ) :: t_initial           ! initial simulation time (almost always 0)
   INTEGER(IntKi),           INTENT(IN   ) :: n_t_global          ! integer time step   
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              ! Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(IN   ) :: y_FAST              ! Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(IN   ) :: m_FAST              ! Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  ! ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  ! BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                ! ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                ! AeroDyn v14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  ! AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 ! InflowWind data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  ! HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  ! SubDyn data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                ! MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                ! FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  ! Data for the MoorDyn module
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                ! IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                ! All the IceDyn data used in time-step loop
      
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              ! Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi)                          :: i, k                ! loop counters
   
   REAL(DbKi)                              :: t_module            ! Current simulation time for module 
   INTEGER(IntKi)                          :: j_ss                ! substep loop counter 
   INTEGER(IntKi)                          :: n_t_module          ! simulation time step, loop counter for individual modules       
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_AdvanceStates'       
   
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""


   !----------------------------------------------------------------------------------------
   ! copy the states at step m_FAST%t_global and get prediction for step t_global_next
   ! (note that we need to copy the states because UpdateStates updates the values
   ! and we need to have the old values [at m_FAST%t_global] for the next j_pc step)
   !----------------------------------------------------------------------------------------
   ! ElastoDyn: get predicted states
   CALL ED_CopyContState   (ED%x( STATE_CURR), ED%x( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL ED_CopyDiscState   (ED%xd(STATE_CURR), ED%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL ED_CopyConstrState (ED%z( STATE_CURR), ED%z( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   IF ( p_FAST%n_substeps( MODULE_ED ) > 1 ) THEN
      CALL ED_CopyOtherState( ED%OtherSt, ED%OtherSt_old, MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF

   DO j_ss = 1, p_FAST%n_substeps( MODULE_ED )
      n_t_module = n_t_global*p_FAST%n_substeps( MODULE_ED ) + j_ss - 1
      t_module   = n_t_module*p_FAST%dt_module( MODULE_ED ) + t_initial
            
      CALL ED_UpdateStates( t_module, n_t_module, ED%Input, ED%InputTimes, ED%p, ED%x(STATE_PRED), ED%xd(STATE_PRED), ED%z(STATE_PRED), ED%OtherSt, ErrStat2, ErrMsg2 )
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               
   END DO !j_ss

         
   IF ( p_FAST%CompElast == Module_BD ) THEN
            
      DO k=1,p_FAST%nBeams
            
         CALL BD_CopyContState   (BD%x( k,STATE_CURR),BD%x( k,STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL BD_CopyDiscState   (BD%xd(k,STATE_CURR),BD%xd(k,STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL BD_CopyConstrState (BD%z( k,STATE_CURR),BD%z( k,STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL BD_CopyOtherState (BD%OtherSt( k,STATE_CURR),BD%OtherSt( k,STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
         DO j_ss = 1, p_FAST%n_substeps( Module_BD )
            n_t_module = n_t_global*p_FAST%n_substeps( Module_BD ) + j_ss - 1
            t_module   = n_t_module*p_FAST%dt_module( Module_BD ) + t_initial
                           
            CALL BD_UpdateStates( t_module, n_t_module, BD%Input(:,k), BD%InputTimes(:,k), BD%p(k), BD%x(k,STATE_PRED), &
                                       BD%xd(k,STATE_PRED), BD%z(k,STATE_PRED), BD%OtherSt(k,STATE_PRED), ErrStat2, ErrMsg2 )
               CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         END DO !j_ss
      END DO !nBeams
         
   END IF !CompElast
   
   
   ! AeroDyn: get predicted states
   IF ( p_FAST%CompAero == Module_AD14 ) THEN
      CALL AD14_CopyContState   (AD14%x( STATE_CURR), AD14%x( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD14_CopyDiscState   (AD14%xd(STATE_CURR), AD14%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD14_CopyConstrState (AD14%z( STATE_CURR), AD14%z( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
      IF ( p_FAST%n_substeps( Module_AD14 ) > 1 ) THEN
         CALL AD14_CopyOtherState( AD14%OtherSt, AD14%OtherSt_old, MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END IF
            
      DO j_ss = 1, p_FAST%n_substeps( MODULE_AD14 )
         n_t_module = n_t_global*p_FAST%n_substeps( MODULE_AD14 ) + j_ss - 1
         t_module   = n_t_module*p_FAST%dt_module( MODULE_AD14 ) + t_initial
            
         CALL AD14_UpdateStates( t_module, n_t_module, AD14%Input, AD14%InputTimes, AD14%p, AD14%x(STATE_PRED), &
                                AD14%xd(STATE_PRED), AD14%z(STATE_PRED), AD14%OtherSt, ErrStat2, ErrMsg2 )
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO !j_ss
   ELSEIF ( p_FAST%CompAero == Module_AD ) THEN
      CALL AD_CopyContState   (AD%x( STATE_CURR), AD%x( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD_CopyDiscState   (AD%xd(STATE_CURR), AD%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD_CopyConstrState (AD%z( STATE_CURR), AD%z( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
      IF ( p_FAST%n_substeps( Module_AD ) > 1 ) THEN
         CALL AD_CopyOtherState( AD%OtherSt, AD%OtherSt_old, MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END IF
            
      DO j_ss = 1, p_FAST%n_substeps( MODULE_AD )
         n_t_module = n_t_global*p_FAST%n_substeps( MODULE_AD ) + j_ss - 1
         t_module   = n_t_module*p_FAST%dt_module( MODULE_AD ) + t_initial
            
         CALL AD_UpdateStates( t_module, n_t_module, AD%Input, AD%InputTimes, AD%p, AD%x(STATE_PRED), AD%xd(STATE_PRED), AD%z(STATE_PRED), AD%OtherSt, ErrStat2, ErrMsg2 )
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO !j_ss
   END IF            

                        
   ! InflowWind: get predicted states
   IF ( p_FAST%CompInflow == Module_IfW ) THEN
      CALL InflowWind_CopyContState   (IfW%x( STATE_CURR), IfW%x( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL InflowWind_CopyDiscState   (IfW%xd(STATE_CURR), IfW%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL InflowWind_CopyConstrState (IfW%z( STATE_CURR), IfW%z( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
      IF ( p_FAST%n_substeps( Module_IfW ) > 1 ) THEN
         CALL InflowWind_CopyOtherState( IfW%OtherSt, IfW%OtherSt_old, MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END IF
            
      DO j_ss = 1, p_FAST%n_substeps( MODULE_IfW )
         n_t_module = n_t_global*p_FAST%n_substeps( MODULE_IfW ) + j_ss - 1
         t_module   = n_t_module*p_FAST%dt_module( MODULE_IfW ) + t_initial
            
         CALL InflowWind_UpdateStates( t_module, n_t_module, IfW%Input, IfW%InputTimes, IfW%p, IfW%x(STATE_PRED), IfW%xd(STATE_PRED), IfW%z(STATE_PRED), IfW%OtherSt, ErrStat2, ErrMsg2 )
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO !j_ss
   END IF          
   
   
   ! ServoDyn: get predicted states
   IF ( p_FAST%CompServo == Module_SrvD ) THEN
      CALL SrvD_CopyContState   (SrvD%x( STATE_CURR), SrvD%x( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SrvD_CopyDiscState   (SrvD%xd(STATE_CURR), SrvD%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SrvD_CopyConstrState (SrvD%z( STATE_CURR), SrvD%z( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
      IF ( p_FAST%n_substeps( Module_SrvD ) > 1 ) THEN
         CALL SrvD_CopyOtherState( SrvD%OtherSt, SrvD%OtherSt_old, MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END IF
         
      DO j_ss = 1, p_FAST%n_substeps( Module_SrvD )
         n_t_module = n_t_global*p_FAST%n_substeps( Module_SrvD ) + j_ss - 1
         t_module   = n_t_module*p_FAST%dt_module( Module_SrvD ) + t_initial
               
         CALL SrvD_UpdateStates( t_module, n_t_module, SrvD%Input, SrvD%InputTimes, SrvD%p, SrvD%x(STATE_PRED), SrvD%xd(STATE_PRED), SrvD%z(STATE_PRED), SrvD%OtherSt, ErrStat2, ErrMsg2 )
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO !j_ss
   END IF            
            

   ! HydroDyn: get predicted states
   IF ( p_FAST%CompHydro == Module_HD ) THEN
      CALL HydroDyn_CopyContState   (HD%x( STATE_CURR), HD%x( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL HydroDyn_CopyDiscState   (HD%xd(STATE_CURR), HD%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL HydroDyn_CopyConstrState (HD%z( STATE_CURR), HD%z( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
      IF ( p_FAST%n_substeps( Module_HD ) > 1 ) THEN
         CALL HydroDyn_CopyOtherState( HD%OtherSt, HD%OtherSt_old, MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END IF
         
      DO j_ss = 1, p_FAST%n_substeps( Module_HD )
         n_t_module = n_t_global*p_FAST%n_substeps( Module_HD ) + j_ss - 1
         t_module   = n_t_module*p_FAST%dt_module( Module_HD ) + t_initial
               
         CALL HydroDyn_UpdateStates( t_module, n_t_module, HD%Input, HD%InputTimes, HD%p, HD%x(STATE_PRED), HD%xd(STATE_PRED), HD%z(STATE_PRED), HD%OtherSt, ErrStat2, ErrMsg2 )
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO !j_ss
            
   END IF
            
         
   ! SubDyn: get predicted states
   IF ( p_FAST%CompSub == Module_SD ) THEN
      CALL SD_CopyContState   (SD%x( STATE_CURR), SD%x( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SD_CopyDiscState   (SD%xd(STATE_CURR), SD%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SD_CopyConstrState (SD%z( STATE_CURR), SD%z( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      IF ( p_FAST%n_substeps( Module_SD ) > 1 ) THEN
         CALL SD_CopyOtherState( SD%OtherSt, SD%OtherSt_old, MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END IF
            
      DO j_ss = 1, p_FAST%n_substeps( Module_SD )
         n_t_module = n_t_global*p_FAST%n_substeps( Module_SD ) + j_ss - 1
         t_module   = n_t_module*p_FAST%dt_module( Module_SD ) + t_initial
               
         CALL SD_UpdateStates( t_module, n_t_module, SD%Input, SD%InputTimes, SD%p, SD%x(STATE_PRED), SD%xd(STATE_PRED), SD%z(STATE_PRED), SD%OtherSt, ErrStat2, ErrMsg2 )
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO !j_ss
   END IF
            
            
   ! MAP/FEAM: get predicted states
   IF (p_FAST%CompMooring == Module_MAP) THEN
      CALL MAP_CopyContState   (MAPp%x( STATE_CURR), MAPp%x( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MAP_CopyDiscState   (MAPp%xd(STATE_CURR), MAPp%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MAP_CopyConstrState (MAPp%z( STATE_CURR), MAPp%z( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      IF ( p_FAST%n_substeps( Module_MAP ) > 1 ) THEN
         CALL MAP_CopyOtherState( MAPp%OtherSt, MAPp%OtherSt_old, MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END IF
         
      DO j_ss = 1, p_FAST%n_substeps( Module_MAP )
         n_t_module = n_t_global*p_FAST%n_substeps( Module_MAP ) + j_ss - 1
         t_module   = n_t_module*p_FAST%dt_module( Module_MAP ) + t_initial
               
         CALL MAP_UpdateStates( t_module, n_t_module, MAPp%Input, MAPp%InputTimes, MAPp%p, MAPp%x(STATE_PRED), MAPp%xd(STATE_PRED), MAPp%z(STATE_PRED), MAPp%OtherSt, ErrStat2, ErrMsg2 )
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO !j_ss
               
   ELSEIF (p_FAST%CompMooring == Module_MD) THEN
      CALL MD_CopyContState   (MD%x( STATE_CURR), MD%x( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MD_CopyDiscState   (MD%xd(STATE_CURR), MD%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MD_CopyConstrState (MD%z( STATE_CURR), MD%z( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
      IF ( p_FAST%n_substeps( Module_MD ) > 1 ) THEN
         CALL MD_CopyOtherState( MD%OtherSt, MD%OtherSt_old, MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END IF
            
      DO j_ss = 1, p_FAST%n_substeps( Module_MD )
         n_t_module = n_t_global*p_FAST%n_substeps( Module_MD ) + j_ss - 1
         t_module   = n_t_module*p_FAST%dt_module( Module_MD ) + t_initial
               
         CALL MD_UpdateStates( t_module, n_t_module, MD%Input, MD%InputTimes, MD%p, MD%x(STATE_PRED), MD%xd(STATE_PRED), MD%z(STATE_PRED), MD%OtherSt, ErrStat2, ErrMsg2 )
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO !j_ss
               
   ELSEIF (p_FAST%CompMooring == Module_FEAM) THEN
      CALL FEAM_CopyContState   (FEAM%x( STATE_CURR), FEAM%x( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL FEAM_CopyDiscState   (FEAM%xd(STATE_CURR), FEAM%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL FEAM_CopyConstrState (FEAM%z( STATE_CURR), FEAM%z( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
      IF ( p_FAST%n_substeps( Module_FEAM ) > 1 ) THEN
         CALL FEAM_CopyOtherState( FEAM%OtherSt, FEAM%OtherSt_old, MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END IF
            
      DO j_ss = 1, p_FAST%n_substeps( Module_FEAM )
         n_t_module = n_t_global*p_FAST%n_substeps( Module_FEAM ) + j_ss - 1
         t_module   = n_t_module*p_FAST%dt_module( Module_FEAM ) + t_initial
               
         CALL FEAM_UpdateStates( t_module, n_t_module, FEAM%Input, FEAM%InputTimes, FEAM%p, FEAM%x(STATE_PRED), FEAM%xd(STATE_PRED), FEAM%z(STATE_PRED), FEAM%OtherSt, ErrStat2, ErrMsg2 )
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO !j_ss
               
   END IF
             
         
   ! IceFloe/IceDyn: get predicted states
   IF ( p_FAST%CompIce == Module_IceF ) THEN
      CALL IceFloe_CopyContState   (IceF%x( STATE_CURR), IceF%x( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL IceFloe_CopyDiscState   (IceF%xd(STATE_CURR), IceF%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL IceFloe_CopyConstrState (IceF%z( STATE_CURR), IceF%z( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      IF ( p_FAST%n_substeps( Module_IceF ) > 1 ) THEN
         CALL IceFloe_CopyOtherState( IceF%OtherSt, IceF%OtherSt_old, MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END IF
            
      DO j_ss = 1, p_FAST%n_substeps( Module_IceF )
         n_t_module = n_t_global*p_FAST%n_substeps( Module_IceF ) + j_ss - 1
         t_module   = n_t_module*p_FAST%dt_module( Module_IceF ) + t_initial
               
         CALL IceFloe_UpdateStates( t_module, n_t_module, IceF%Input, IceF%InputTimes, IceF%p, IceF%x(STATE_PRED), IceF%xd(STATE_PRED), IceF%z(STATE_PRED), IceF%OtherSt, ErrStat2, ErrMsg2 )
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO !j_ss
   ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
            
      DO i=1,p_FAST%numIceLegs
            
         CALL IceD_CopyContState   (IceD%x( i,STATE_CURR),IceD%x( i,STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL IceD_CopyDiscState   (IceD%xd(i,STATE_CURR),IceD%xd(i,STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL IceD_CopyConstrState (IceD%z( i,STATE_CURR),IceD%z( i,STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         IF ( p_FAST%n_substeps( Module_IceD ) > 1 ) THEN
            CALL IceD_CopyOtherState( IceD%OtherSt(i), IceD%OtherSt_old(I), MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         END IF
            
         DO j_ss = 1, p_FAST%n_substeps( Module_IceD )
            n_t_module = n_t_global*p_FAST%n_substeps( Module_IceD ) + j_ss - 1
            t_module   = n_t_module*p_FAST%dt_module( Module_IceD ) + t_initial
               
            CALL IceD_UpdateStates( t_module, n_t_module, IceD%Input(:,i), IceD%InputTimes(:,i), IceD%p(i), IceD%x(i,STATE_PRED), &
                                       IceD%xd(i,STATE_PRED), IceD%z(i,STATE_PRED), IceD%OtherSt(i), ErrStat2, ErrMsg2 )
               CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         END DO !j_ss
      END DO
         
   END IF
         
END SUBROUTINE FAST_AdvanceStates
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FAST_EndMods( p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, HD, SD, MAPp, FEAM, MD, IceF, IceD, ErrStat, ErrMsg )
!...............................................................................................................................   

   TYPE(FAST_ParameterType), INTENT(INOUT) :: p_FAST              ! Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              ! Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              ! Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  ! ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  ! BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                ! ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                ! AeroDyn v14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  ! AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 ! InflowWind data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  ! HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  ! SubDyn data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                ! MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                ! FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  ! Data for the MoorDyn module
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                ! IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                ! All the IceDyn data used in time-step loop
      
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              ! Error message if ErrStat /= ErrID_None
   
   ! local variables
   INTEGER(IntKi)                          :: i, k                ! loop counter
   
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_EndMods'
                  
      !...............................................................................................................................
      ! End all modules (and write binary FAST output file)
      !...............................................................................................................................

   ErrStat = ErrID_None
   ErrMsg  = ""
            
      
   CALL FAST_EndOutput( p_FAST, y_FAST, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   IF ( p_FAST%ModuleInitialized(Module_ED) ) THEN
      CALL ED_End(   ED%Input(1),   ED%p,   ED%x(STATE_CURR),   ED%xd(STATE_CURR),   ED%z(STATE_CURR),   ED%OtherSt,   &
                     ED%Output(1),   ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   END IF

   IF ( p_FAST%ModuleInitialized(Module_BD) ) THEN
         
      DO k=1,p_FAST%nBeams                     
         CALL BD_End(BD%Input(1,k),  BD%p(k),  BD%x(k,STATE_CURR),  BD%xd(k,STATE_CURR),  BD%z(k,STATE_CURR), &
                        BD%OtherSt(k,STATE_CURR),  BD%y(k),  ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)            
      END DO
         
   END IF   
   
   
   IF ( p_FAST%ModuleInitialized(Module_AD14) ) THEN
      CALL AD14_End( AD14%Input(1), AD14%p, AD14%x(STATE_CURR), AD14%xd(STATE_CURR), AD14%z(STATE_CURR), AD14%OtherSt, &
                     AD14%y, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   ELSEIF ( p_FAST%ModuleInitialized(Module_AD) ) THEN
      CALL AD_End(   AD%Input(1),   AD%p,   AD%x(STATE_CURR),   AD%xd(STATE_CURR),   AD%z(STATE_CURR),   AD%OtherSt,   &
                     AD%y,   ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)   
   END IF
      
   IF ( p_FAST%ModuleInitialized(Module_IfW) ) THEN
      CALL InflowWind_End( IfW%Input(1), IfW%p, IfW%x(STATE_CURR), IfW%xd(STATE_CURR), IfW%z(STATE_CURR), IfW%OtherSt,   &
                           IfW%y,   ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   END IF   
   
   IF ( p_FAST%ModuleInitialized(Module_SrvD) ) THEN
      CALL SrvD_End( SrvD%Input(1), SrvD%p, SrvD%x(STATE_CURR), SrvD%xd(STATE_CURR), SrvD%z(STATE_CURR), SrvD%OtherSt, &
                     SrvD%y, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   END IF

   IF ( p_FAST%ModuleInitialized(Module_HD) ) THEN
      CALL HydroDyn_End(    HD%Input(1),   HD%p,   HD%x(STATE_CURR),   HD%xd(STATE_CURR),   HD%z(STATE_CURR),   HD%OtherSt,   &
                              HD%y,   ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   END IF

   IF ( p_FAST%ModuleInitialized(Module_SD) ) THEN
      CALL SD_End(    SD%Input(1),   SD%p,   SD%x(STATE_CURR),   SD%xd(STATE_CURR),   SD%z(STATE_CURR),   SD%OtherSt,   &
                        SD%y,   ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   END IF
      
   IF ( p_FAST%ModuleInitialized(Module_MAP) ) THEN
      CALL MAP_End(    MAPp%Input(1),   MAPp%p,   MAPp%x(STATE_CURR),   MAPp%xd(STATE_CURR),   MAPp%z(STATE_CURR),   MAPp%OtherSt,   &
                        MAPp%y,   ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   ELSEIF ( p_FAST%ModuleInitialized(Module_MD) ) THEN
      CALL MD_End(  MD%Input(1), MD%p, MD%x(STATE_CURR), MD%xd(STATE_CURR), MD%z(STATE_CURR), MD%OtherSt, MD%y, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   ELSEIF ( p_FAST%ModuleInitialized(Module_FEAM) ) THEN
      CALL FEAM_End(   FEAM%Input(1),  FEAM%p,  FEAM%x(STATE_CURR),  FEAM%xd(STATE_CURR),  FEAM%z(STATE_CURR),  FEAM%OtherSt,  &
                        FEAM%y,  ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   END IF
      
   IF ( p_FAST%ModuleInitialized(Module_IceF) ) THEN
      CALL IceFloe_End(IceF%Input(1),  IceF%p,  IceF%x(STATE_CURR),  IceF%xd(STATE_CURR),  IceF%z(STATE_CURR),  IceF%OtherSt, &
                        IceF%y,  ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   ELSEIF ( p_FAST%ModuleInitialized(Module_IceD) ) THEN
         
      DO i=1,p_FAST%numIceLegs                     
         CALL IceD_End(IceD%Input(1,i),  IceD%p(i),  IceD%x(i,STATE_CURR),  IceD%xd(i,STATE_CURR),  IceD%z(i,STATE_CURR), &
                        IceD%OtherSt(i),  IceD%y(i),  ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)            
      END DO
         
   END IF   
                     
END SUBROUTINE FAST_EndMods
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FAST_DestroyAll( p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, OpFM, HD, SD, MAPp, FEAM, MD, IceF, IceD, MeshMapData, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType), INTENT(INOUT) :: p_FAST              ! Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              ! Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              ! Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  ! ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  ! BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                ! ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                ! AeroDyn v14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  ! AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 ! InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                ! OpenFOAM data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  ! HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  ! SubDyn data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                ! MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                ! FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  ! Data for the MoorDyn module
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                ! IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                ! All the IceDyn data used in time-step loop

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         ! Data for mapping between modules
      
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              ! Error message if ErrStat /= ErrID_None
   
   ! local variables
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_DestroyAll'


   
   ! -------------------------------------------------------------------------
   ! Deallocate/Destroy structures associated with mesh mapping
   ! -------------------------------------------------------------------------

   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
   ! FAST
      
   CALL FAST_DestroyParam( p_FAST, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   
   CALL FAST_DestroyOutputFileType( y_FAST, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   
   CALL FAST_DestroyMiscVarType( m_FAST, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   
   ! ElastoDyn
   CALL FAST_DestroyElastoDyn_Data( ED, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)      
      
   ! BeamDyn
   CALL FAST_DestroyBeamDyn_Data( BD, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)            
      
   ! ServoDyn
   CALL FAST_DestroyServoDyn_Data( SrvD, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)      
   
   ! AeroDyn14
   CALL FAST_DestroyAeroDyn14_Data( AD14, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)      

   ! AeroDyn
   CALL FAST_DestroyAeroDyn_Data( AD, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)            
      
   ! InflowWind
   CALL FAST_DestroyInflowWind_Data( IfW, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)            
      
   ! OpenFOAM
   CALL FAST_DestroyOpenFOAM_Data( OpFM, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)            
            
   ! HydroDyn
   CALL FAST_DestroyHydroDyn_Data( HD, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)      
   
   ! SubDyn
   CALL FAST_DestroySubDyn_Data( SD, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)      
      
   ! MAP      
   CALL FAST_DestroyMAP_Data( MAPp, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)      
            
   ! FEAMooring 
   CALL FAST_DestroyFEAMooring_Data( FEAM, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)      

   ! MoorDyn 
   CALL FAST_DestroyMoorDyn_Data( MD, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)      

   ! IceFloe
   CALL FAST_DestroyIceFloe_Data( IceF, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)      
            
   ! IceDyn
   CALL FAST_DestroyIceDyn_Data( IceD, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)      

   ! Module (Mesh) Mapping data
   CALL FAST_DestroyModuleMapType( MeshMapData, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)            
       
      
   
END SUBROUTINE FAST_DestroyAll
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ExitThisProgram_T( Turbine, ErrLevel_in, ErrLocMsg )
   
   TYPE(FAST_TurbineType),   INTENT(INOUT) :: Turbine
   INTEGER(IntKi),           INTENT(IN)    :: ErrLevel_in         ! Error level when Error == .TRUE. (required when Error is .TRUE.)
   CHARACTER(*), OPTIONAL,   INTENT(IN)    :: ErrLocMsg           ! an optional message describing the location of the error
   
   
   IF (PRESENT(ErrLocMsg)) THEN
      
      CALL ExitThisProgram( Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
                     Turbine%ED, Turbine%BD, Turbine%SrvD, Turbine%AD14, Turbine%AD, Turbine%IfW, Turbine%OpFM, &
                     Turbine%HD, Turbine%SD, Turbine%MAP, Turbine%FEAM, Turbine%MD, &
                     Turbine%IceF, Turbine%IceD, Turbine%MeshMapData, ErrLevel_in, ErrLocMsg )
   
   ELSE     
      
      CALL ExitThisProgram( Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
                     Turbine%ED, Turbine%BD, Turbine%SrvD, Turbine%AD14, Turbine%AD, Turbine%IfW, Turbine%OpFM, &
                     Turbine%HD, Turbine%SD, Turbine%MAP, Turbine%FEAM, Turbine%MD, &
                     Turbine%IceF, Turbine%IceD, Turbine%MeshMapData, ErrLevel_in )
      
   END IF

END SUBROUTINE ExitThisProgram_T
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ExitThisProgram( p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, OpFM, HD, SD, MAPp, FEAM, MD, IceF, IceD, MeshMapData, ErrLevel_in, ErrLocMsg )
! This subroutine is called when FAST exits. It calls all the modules' end routines and cleans up variables declared in the
! main program. If there was an error, it also aborts. Otherwise, it prints the run times and performs a normal exit.
!...............................................................................................................................

      ! Passed arguments
   TYPE(FAST_ParameterType), INTENT(INOUT) :: p_FAST              ! Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              ! Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              ! Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  ! ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  ! BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                ! ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                ! AeroDyn v14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  ! AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 ! InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                ! OpenFOAM data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  ! HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  ! SubDyn data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                ! MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                ! FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  ! Data for the MoorDyn module
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                ! IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                ! All the IceDyn data used in time-step loop

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         ! Data for mapping between modules

   INTEGER(IntKi),           INTENT(IN)    :: ErrLevel_in         ! Error level when Error == .TRUE. (required when Error is .TRUE.)
   CHARACTER(*), OPTIONAL,   INTENT(IN)    :: ErrLocMsg           ! an optional message describing the location of the error


      ! Local variables:            
   INTEGER(IntKi)                          :: ErrorLevel
                                          
   INTEGER(IntKi)                          :: ErrStat2            ! Error status
   CHARACTER(1024)                         :: ErrMsg2             ! Error message
   CHARACTER(1024)                         :: SimMsg              ! optional message to print about where the error took place in the simulation
   
   CHARACTER(*), PARAMETER                 :: RoutineName = 'ExitThisProgram'       

      
   ErrorLevel = ErrLevel_in
      
      
      ! End all modules
   CALL FAST_EndMods( p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, HD, SD, MAPp, FEAM, MD, IceF, IceD, ErrStat2, ErrMsg2 )
      IF (ErrStat2 /= ErrID_None) THEN
         CALL WrScr( NewLine//RoutineName//':'//TRIM(ErrMsg2)//NewLine )
         ErrorLevel = MAX(ErrorLevel,ErrStat2)
      END IF
                  
      ! Destroy all data associated with FAST variables:

   CALL FAST_DestroyAll( p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, OpFM, HD, SD, MAPp, FEAM, MD, IceF, IceD, MeshMapData, ErrStat2, ErrMsg2 )
      IF (ErrStat2 /= ErrID_None) THEN
         CALL WrScr( NewLine//RoutineName//':'//TRIM(ErrMsg2)//NewLine )
         ErrorLevel = MAX(ErrorLevel,ErrStat2)
      END IF

      
   !............................................................................................................................
   ! Set exit error code if there was an error;
   !............................................................................................................................
   IF ( ErrorLevel >= AbortErrLev ) THEN
      
      IF (PRESENT(ErrLocMsg)) THEN
         SimMsg = ErrLocMsg
      ELSE
         SimMsg = 'after the simulation was complete'
      END IF
                                         
      CALL ProgAbort( 'FAST encountered an error '//TRIM(SimMsg)//'.'//NewLine//' Simulation error level: '&
                        //TRIM(GetErrStr(ErrorLevel)), TrapErrors=.FALSE., TimeWait=3._ReKi )  ! wait 3 seconds (in case they double-clicked and got an error)
   END IF
      
   !............................................................................................................................
   !  Write simulation times and stop
   !............................................................................................................................

   CALL RunTimes( m_FAST%StrtTime, m_FAST%UsrTime1, m_FAST%SimStrtTime, m_FAST%UsrTime2, m_FAST%t_global )

#if (defined COMPILE_SIMULINK || defined COMPILE_LABVIEW)
   ! for Simulink, this may not be a normal stop. It might call this after an error in the model.
   CALL WrScr( NewLine//' '//TRIM(FAST_Ver%Name)//' completed.'//NewLine )
#else   
   CALL NormStop( )
#endif   


END SUBROUTINE ExitThisProgram
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FAST_Solution0_T(Turbine, ErrStat, ErrMsg)

   TYPE(FAST_TurbineType),   INTENT(INOUT) :: Turbine             ! all data for one instance of a turbine
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              ! Error message if ErrStat /= ErrID_None


      CALL FAST_Solution0(Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
                     Turbine%ED, Turbine%BD, Turbine%SrvD, Turbine%AD14, Turbine%AD, Turbine%IfW, Turbine%OpFM, &
                     Turbine%HD, Turbine%SD, Turbine%MAP, Turbine%FEAM, Turbine%MD, &
                     Turbine%IceF, Turbine%IceD, Turbine%MeshMapData, ErrStat, ErrMsg )
      
END SUBROUTINE FAST_Solution0_T
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FAST_Solution0(p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, OpFM, HD, SD, MAPp, FEAM, MD, IceF, IceD, MeshMapData, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              ! Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              ! Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              ! Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  ! ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  ! BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                ! ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                ! AeroDyn14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  ! AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 ! InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                ! OpenFOAM data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  ! HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  ! SubDyn data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                ! MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                ! FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  ! Data for the MoorDyn module
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                ! IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                ! All the IceDyn data used in time-step loop

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         ! Data for mapping between modules
      
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              ! Error message if ErrStat /= ErrID_None
   
   ! local variables
   INTEGER(IntKi), PARAMETER               :: n_t_global = -1     ! loop counter
   
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_Solution0'

   
   !NOTE: m_FAST%t_global is t_initial in this routine
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
   CALL SimStatus_FirstTime( m_FAST%TiLstPrn, m_FAST%PrevClockTime, m_FAST%SimStrtTime, m_FAST%UsrTime2, m_FAST%t_global, p_FAST%TMax )

   ! Solve input-output relations; this section of code corresponds to Eq. (35) in Gasmi et al. (2013)
   ! This code will be specific to the underlying modules
   
#ifdef SOLVE_OPTION_1_BEFORE_2
! used for Option 1 before Option 2:

   IF ( p_FAST%CompSub == Module_SD .OR. p_FAST%CompHydro == Module_HD ) THEN
   ! Because SubDyn needs a better initial guess from ElastoDyn, we'll add an additional call to ED_CalcOutput to get them:
   ! (we'll do the same for HydroDyn, though I'm not sure it's as critical)
   
      CALL ED_CalcOutput( m_FAST%t_global, ED%Input(1), ED%p, ED%x(STATE_CURR), ED%xd(STATE_CURR), ED%z(STATE_CURR), ED%OtherSt, &
                          ED%Output(1), ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      CALL Transfer_ED_to_HD_SD_BD_Mooring( p_FAST, ED%Output(1), HD%Input(1), SD%Input(1), MAPp%Input(1), FEAM%Input(1), MD%Input(1), &
                                            BD%Input(1,:), MeshMapData, ErrStat2, ErrMsg2 )         
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               
   END IF   
#endif   

      ! the initial ServoDyn and IfW/Lidar inputs from Simulink:
   IF ( p_FAST%CompServo == Module_SrvD ) CALL SrvD_SetExternalInputs( p_FAST, m_FAST, SrvD%Input(1) )   
   IF ( p_FAST%CompInflow == Module_IfW ) CALL IfW_SetExternalInputs( IfW%p, m_FAST, ED%Output(1), IfW%Input(1) )  

   CALL CalcOutputs_And_SolveForInputs(  n_t_global, m_FAST%t_global,  STATE_CURR, m_FAST%calcJacobian, m_FAST%NextJacCalcTime, &
                        p_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, OpFM, HD, SD, MAPp, FEAM, MD, IceF, IceD, MeshMapData, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
      
      
! output some graphics. For now, WrGraphics is hard-coded to .FALSE.      
   IF (p_FAST%WrGraphics) THEN
      CALL WriteInputMeshesToFile( ED%Input(1), AD%Input(1), SD%Input(1), HD%Input(1), MAPp%Input(1), BD%Input(1,:), TRIM(p_FAST%OutFileRoot)//'.InputMeshes.bin', ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF 



   !----------------------------------------------------------------------------------------
   ! Check to see if we should output data this time step:
   !----------------------------------------------------------------------------------------

   CALL WriteOutputToFile(m_FAST%t_global, p_FAST, y_FAST, ED, BD, AD14, AD, IfW, OpFM, HD, SD, SrvD, MAPp, FEAM, MD, IceF, IceD, ErrStat2, ErrMsg2)   
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      
   !...............
   ! Copy values of these initial guesses for interpolation/extrapolation and 
   ! initialize predicted states for j_pc loop (use MESH_NEWCOPY here so we can use MESH_UPDATE copy later)
   !...............
         
   ! Initialize Input-Output arrays for interpolation/extrapolation:

   CALL FAST_InitIOarrays( m_FAST%t_global, p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, HD, SD, MAPp, FEAM, MD, IceF, IceD, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         

END SUBROUTINE FAST_Solution0
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FAST_Solution_T(t_initial, n_t_global, Turbine, ErrStat, ErrMsg )

   REAL(DbKi),               INTENT(IN   ) :: t_initial           ! initial time
   INTEGER(IntKi),           INTENT(IN   ) :: n_t_global          ! loop counter
   TYPE(FAST_TurbineType),   INTENT(INOUT) :: Turbine             ! all data for one instance of a turbine
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              ! Error message if ErrStat /= ErrID_None
   
      CALL FAST_Solution(t_initial, n_t_global, Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
                  Turbine%ED, Turbine%BD, Turbine%SrvD, Turbine%AD14, Turbine%AD, Turbine%IfW, Turbine%OpFM, &
                  Turbine%HD, Turbine%SD, Turbine%MAP, Turbine%FEAM, Turbine%MD, &
                  Turbine%IceF, Turbine%IceD, Turbine%MeshMapData, ErrStat, ErrMsg )                  
                  
END SUBROUTINE FAST_Solution_T
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FAST_Solution(t_initial, n_t_global, p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, OpFM, HD, SD, MAPp, FEAM, MD, IceF, IceD, MeshMapData, ErrStat, ErrMsg )
! this routine takes data from n_t_global and gets values at n_t_global + 1

   REAL(DbKi),               INTENT(IN   ) :: t_initial           ! initial time
   INTEGER(IntKi),           INTENT(IN   ) :: n_t_global          ! loop counter

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              ! Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              ! Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              ! Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  ! ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  ! BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                ! ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                ! AeroDyn14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  ! AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 ! InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                ! OpenFOAM data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  ! HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  ! SubDyn data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                ! MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                ! FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  ! Data for the MoorDyn module
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                ! IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                ! All the IceDyn data used in time-step loop

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         ! Data for mapping between modules
      
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              ! Error message if ErrStat /= ErrID_None
   
   ! local variables
   REAL(DbKi)                              :: t_global_next       ! next simulation time (m_FAST%t_global + p_FAST%dt)
   INTEGER(IntKi)                          :: j_pc                ! predictor-corrector loop counter 
   
   INTEGER(IntKi)                          :: I, k                ! generic loop counters
   
   
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_Solution'


   ErrStat = ErrID_None
   ErrMsg  = ""
   
   t_global_next = t_initial + (n_t_global+1)*p_FAST%DT  ! = m_FAST%t_global + p_FAST%dt
                       
      ! determine if the Jacobian should be calculated this time
   IF ( m_FAST%calcJacobian ) THEN ! this was true (possibly at initialization), so we'll advance the time for the next calculation of the Jacobian
      m_FAST%NextJacCalcTime = m_FAST%t_global + p_FAST%DT_UJac         
   END IF
      
      ! the ServoDyn inputs from Simulink are for t, not t+dt, so we're going to overwrite the inputs from
      ! the previous step before we extrapolate these inputs:
   IF ( p_FAST%CompServo == Module_SrvD ) CALL SrvD_SetExternalInputs( p_FAST, m_FAST, SrvD%Input(1) )   
   
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Step 1.a: Extrapolate Inputs (gives predicted values at t+dt)
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   CALL FAST_ExtrapInterpMods( t_global_next, p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, HD, SD, MAPp, FEAM, MD, IceF, IceD, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      
   ! predictor-corrector loop:
   DO j_pc = 0, p_FAST%NumCrctn
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Step 1.b: Advance states (yield state and constraint values at t_global_next)
   !
   ! x, xd, and z contain values at m_FAST%t_global;
   ! values at t_global_next are stored in the *_pred variables.
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      CALL FAST_AdvanceStates( t_initial, n_t_global, p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, HD, SD, MAPp, FEAM, MD, IceF, IceD, ErrStat2, ErrMsg2 )               
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN
         
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Step 1.c: Input-Output Solve      
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      CALL CalcOutputs_And_SolveForInputs( n_t_global, t_global_next,  STATE_PRED, m_FAST%calcJacobian, m_FAST%NextJacCalcTime, &
                  p_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, OpFM, HD, SD, MAPp, FEAM, MD, IceF, IceD, MeshMapData, ErrStat2, ErrMsg2 )            
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN
         
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Step 2: Correct (continue in loop) 
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF ( j_pc /= p_FAST%NumCrctn)  THEN          ! Don't copy these on the last loop iteration...
                  
         IF ( p_FAST%n_substeps( Module_ED ) > 1 ) THEN
            CALL ED_CopyOtherState( ED%OtherSt_old, ED%OtherSt, MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         END IF
            
         ! BD treats "OtherStates" as actual "other states", which are copied like the rest of the states (done now to avoid changing later when we add MiscVars)
         
         IF ( p_FAST%n_substeps( Module_AD14 ) > 1 ) THEN
            CALL AD14_CopyOtherState( AD14%OtherSt_old, AD14%OtherSt, MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         ELSEIF ( p_FAST%n_substeps( Module_AD ) > 1 ) THEN
            CALL AD_CopyOtherState( AD%OtherSt_old, AD%OtherSt, MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         END IF
            
         IF ( p_FAST%n_substeps( Module_IfW ) > 1 ) THEN
            CALL InflowWind_CopyOtherState( IfW%OtherSt_old, IfW%OtherSt, MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         END IF
                     
         IF ( p_FAST%n_substeps( Module_SrvD ) > 1 ) THEN
            CALL SrvD_CopyOtherState( SrvD%OtherSt_old, SrvD%OtherSt, MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         END IF
            
         IF ( p_FAST%n_substeps( Module_HD ) > 1 ) THEN
            CALL HydroDyn_CopyOtherState( HD%OtherSt_old, HD%OtherSt, MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         END IF
            
         IF ( p_FAST%n_substeps( Module_SD ) > 1 ) THEN
            CALL SD_CopyOtherState( SD%OtherSt_old, SD%OtherSt, MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         END IF

         IF ( p_FAST%n_substeps( Module_MAP ) > 1 ) THEN
            CALL MAP_CopyOtherState( MAPp%OtherSt_old, MAPp%OtherSt, MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         ELSEIF ( p_FAST%n_substeps( Module_MD ) > 1 ) THEN
            CALL MD_CopyOtherState( MD%OtherSt_old, MD%OtherSt, MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         ELSEIF ( p_FAST%n_substeps( Module_FEAM ) > 1 ) THEN
            CALL FEAM_CopyOtherState( FEAM%OtherSt_old, FEAM%OtherSt, MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         END IF
         
         IF ( p_FAST%n_substeps( Module_IceF ) > 1 ) THEN
            CALL IceFloe_CopyOtherState( IceF%OtherSt_old, IceF%OtherSt, MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         ELSEIF ( p_FAST%n_substeps( Module_IceD ) > 1 ) THEN
            DO i=1,p_FAST%numIceLegs
               CALL IceD_CopyOtherState( IceD%OtherSt_old(i), IceD%OtherSt(i), MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            END DO
         END IF
            
      END IF
      
      IF (ErrStat >= AbortErrLev) RETURN
                              
   enddo ! j_pc
      
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Step 3: Save all final variables (advance to next time)
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
   !----------------------------------------------------------------------------------------
   ! copy the final predicted states from step t_global_next to actual states for that step
   !----------------------------------------------------------------------------------------
      
   ! ElastoDyn: copy final predictions to actual states
   CALL ED_CopyContState   (ED%x( STATE_PRED), ED%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL ED_CopyDiscState   (ED%xd(STATE_PRED), ED%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL ED_CopyConstrState (ED%z( STATE_PRED), ED%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)      
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      
         ! BeamDyn: copy final predictions to actual states
   IF ( p_FAST%CompElast == Module_BD ) THEN
      DO k=1,p_FAST%nBeams
         CALL BD_CopyContState   (BD%x( k,STATE_PRED), BD%x( k,STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL BD_CopyDiscState   (BD%xd(k,STATE_PRED), BD%xd(k,STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL BD_CopyConstrState (BD%z( k,STATE_PRED), BD%z( k,STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL BD_CopyOtherState (BD%OtherSt( k,STATE_PRED), BD%OtherSt( k,STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
   END IF

   
   ! AeroDyn: copy final predictions to actual states; copy current outputs to next 
   IF ( p_FAST%CompAero == Module_AD14 ) THEN
      CALL AD14_CopyContState   (AD14%x( STATE_PRED), AD14%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD14_CopyDiscState   (AD14%xd(STATE_PRED), AD14%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD14_CopyConstrState (AD14%z( STATE_PRED), AD14%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)      
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ELSEIF ( p_FAST%CompAero == Module_AD ) THEN
      CALL AD_CopyContState   (AD%x( STATE_PRED), AD%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD_CopyDiscState   (AD%xd(STATE_PRED), AD%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD_CopyConstrState (AD%z( STATE_PRED), AD%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)      
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF
            
      
   ! InflowWind: copy final predictions to actual states; copy current outputs to next 
   IF ( p_FAST%CompInflow == Module_IfW ) THEN
      CALL InflowWind_CopyContState   (IfW%x( STATE_PRED), IfW%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL InflowWind_CopyDiscState   (IfW%xd(STATE_PRED), IfW%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL InflowWind_CopyConstrState (IfW%z( STATE_PRED), IfW%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)      
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF
            
      
   ! ServoDyn: copy final predictions to actual states; copy current outputs to next 
   IF ( p_FAST%CompServo == Module_SrvD ) THEN
      CALL SrvD_CopyContState   (SrvD%x( STATE_PRED), SrvD%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SrvD_CopyDiscState   (SrvD%xd(STATE_PRED), SrvD%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SrvD_CopyConstrState (SrvD%z( STATE_PRED), SrvD%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)      
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF
      
      
   ! HydroDyn: copy final predictions to actual states
   IF ( p_FAST%CompHydro == Module_HD ) THEN         
      CALL HydroDyn_CopyContState   (HD%x( STATE_PRED), HD%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL HydroDyn_CopyDiscState   (HD%xd(STATE_PRED), HD%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL HydroDyn_CopyConstrState (HD%z( STATE_PRED), HD%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF
            
            
   ! SubDyn: copy final predictions to actual states
   IF ( p_FAST%CompSub == Module_SD ) THEN
      CALL SD_CopyContState   (SD%x( STATE_PRED), SD%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SD_CopyDiscState   (SD%xd(STATE_PRED), SD%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SD_CopyConstrState (SD%z( STATE_PRED), SD%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF
         
      
   ! MAP: copy final predictions to actual states
   IF (p_FAST%CompMooring == Module_MAP) THEN
      CALL MAP_CopyContState   (MAPp%x( STATE_PRED), MAPp%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MAP_CopyDiscState   (MAPp%xd(STATE_PRED), MAPp%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MAP_CopyConstrState (MAPp%z( STATE_PRED), MAPp%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ELSEIF (p_FAST%CompMooring == Module_MD) THEN
      CALL MD_CopyContState   (MD%x( STATE_PRED), MD%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MD_CopyDiscState   (MD%xd(STATE_PRED), MD%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MD_CopyConstrState (MD%z( STATE_PRED), MD%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ELSEIF (p_FAST%CompMooring == Module_FEAM) THEN
      CALL FEAM_CopyContState   (FEAM%x( STATE_PRED), FEAM%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL FEAM_CopyDiscState   (FEAM%xd(STATE_PRED), FEAM%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL FEAM_CopyConstrState (FEAM%z( STATE_PRED), FEAM%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF
             
         ! IceFloe: copy final predictions to actual states
   IF ( p_FAST%CompIce == Module_IceF ) THEN
      CALL IceFloe_CopyContState   (IceF%x( STATE_PRED), IceF%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL IceFloe_CopyDiscState   (IceF%xd(STATE_PRED), IceF%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL IceFloe_CopyConstrState (IceF%z( STATE_PRED), IceF%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
      DO i=1,p_FAST%numIceLegs
         CALL IceD_CopyContState   (IceD%x( i,STATE_PRED), IceD%x( i,STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL IceD_CopyDiscState   (IceD%xd(i,STATE_PRED), IceD%xd(i,STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL IceD_CopyConstrState (IceD%z( i,STATE_PRED), IceD%z( i,STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
   END IF

            
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! We've advanced everything to the next time step: 
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                       
      
   ! update the global time 
  
   m_FAST%t_global = t_global_next 
      
      
   !----------------------------------------------------------------------------------------
   ! Check to see if we should output data this time step:
   !----------------------------------------------------------------------------------------

   CALL WriteOutputToFile(m_FAST%t_global, p_FAST, y_FAST, ED, BD, AD14, AD, IfW, OpFM, HD, SD, SrvD, MAPp, FEAM, MD, IceF, IceD, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   !----------------------------------------------------------------------------------------
   ! Display simulation status every SttsTime-seconds (i.e., n_SttsTime steps):
   !----------------------------------------------------------------------------------------   
      
   IF ( MOD( n_t_global + 1, p_FAST%n_SttsTime ) == 0 ) THEN
      
      if (.not. Cmpl4SFun) then   
         CALL SimStatus( m_FAST%TiLstPrn, m_FAST%PrevClockTime, m_FAST%t_global, p_FAST%TMax )
      end if
      
   ENDIF   
     
END SUBROUTINE FAST_Solution
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FAST_CreateCheckpoint_Tary(t_initial, n_t_global, Turbine, CheckpointRoot, ErrStat, ErrMsg)

   REAL(DbKi),               INTENT(IN   ) :: t_initial           ! initial time
   INTEGER(IntKi),           INTENT(IN   ) :: n_t_global          ! loop counter
   TYPE(FAST_TurbineType),   INTENT(INOUT) :: Turbine(:)          ! all data for all turbines
   CHARACTER(*),             INTENT(IN   ) :: CheckpointRoot      ! Rootname of checkpoint file
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              ! Error message if ErrStat /= ErrID_None

      ! local variables
   INTEGER(IntKi)                          :: NumTurbines         ! Number of turbines in this simulation
   INTEGER(IntKi)                          :: i_turb
   INTEGER                                 :: Unit
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(1024)                         :: ErrMsg2             ! local error message
   CHARACTER(*),             PARAMETER     :: RoutineName = 'FAST_CreateCheckpoint_Tary' 
   
   
   NumTurbines = SIZE(Turbine)   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   ! TRIM(CheckpointRoot)//'.'//TRIM(Num2LStr(Turbine%TurbID))//
   
      ! This allows us to put all the turbine data in one file.
   Unit = -1         
   DO i_turb = 1,NumTurbines
      CALL FAST_CreateCheckpoint_T(t_initial, n_t_global, NumTurbines, Turbine(i_turb), CheckpointRoot, ErrStat2, ErrMsg2, Unit )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev ) RETURN
   END DO
               
   
END SUBROUTINE FAST_CreateCheckpoint_Tary
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FAST_CreateCheckpoint_T(t_initial, n_t_global, NumTurbines, Turbine, CheckpointRoot, ErrStat, ErrMsg, Unit )

   USE BladedInterface, ONLY: CallBladedDLL  ! Hack for Bladed-style DLL

   REAL(DbKi),               INTENT(IN   ) :: t_initial           ! initial time
   INTEGER(IntKi),           INTENT(IN   ) :: n_t_global          ! loop counter
   INTEGER(IntKi),           INTENT(IN   ) :: NumTurbines         ! Number of turbines in this simulation
   TYPE(FAST_TurbineType),   INTENT(INOUT) :: Turbine             ! all data for one instance of a turbine (INTENT(OUT) only because of hack for Bladed DLL)
   CHARACTER(*),             INTENT(IN   ) :: CheckpointRoot      ! Rootname of checkpoint file
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              ! Error message if ErrStat /= ErrID_None
   INTEGER(IntKi), OPTIONAL, INTENT(INOUT) :: Unit                ! unit number for output file 
   
      ! local variables:
   REAL(ReKi),               ALLOCATABLE   :: ReKiBuf(:)
   REAL(DbKi),               ALLOCATABLE   :: DbKiBuf(:)
   INTEGER(IntKi),           ALLOCATABLE   :: IntKiBuf(:)
   
   INTEGER(B4Ki)                           :: ArraySizes(3) 
   
   INTEGER(IntKi)                          :: unOut               ! unit number for output file 
   INTEGER(IntKi)                          :: old_avrSwap1        ! previous value of avrSwap(1) !hack for Bladed DLL checkpoint/restore
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(1024)                         :: ErrMsg2             ! local error message
   CHARACTER(*),             PARAMETER     :: RoutineName = 'FAST_CreateCheckpoint_T' 
  
   CHARACTER(1024)                         :: FileName            ! Name of the (output) checkpoint file
   CHARACTER(1024)                         :: DLLFileName         ! Name of the (output) checkpoint file
   
      ! init error status
   ErrStat = ErrID_None
   ErrMsg  = ""
   
      ! Get the arrays of data to be stored in the output file
   CALL FAST_PackTurbineType( ReKiBuf, DbKiBuf, IntKiBuf, Turbine, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev ) RETURN
      
   ArraySizes = 0   
   IF ( ALLOCATED(ReKiBuf)  ) ArraySizes(1) = SIZE(ReKiBuf)
   IF ( ALLOCATED(DbKiBuf)  ) ArraySizes(2) = SIZE(DbKiBuf)
   IF ( ALLOCATED(IntKiBuf) ) ArraySizes(3) = SIZE(IntKiBuf)

   FileName    = TRIM(CheckpointRoot)//'.chkp'
   DLLFileName = TRIM(CheckpointRoot)//'.dll.chkp'

   unOut=-1      
   IF (PRESENT(Unit)) unOut = Unit
         
   IF ( unOut < 0 ) THEN

      CALL GetNewUnit( unOut, ErrStat2, ErrMsg2 )      
      CALL OpenBOutFile ( unOut, FileName, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev ) RETURN
  
         ! checkpoint file header:
      WRITE (unOut, IOSTAT=ErrStat2)   INT(ReKi              ,B4Ki)     ! let's make sure we've got the correct number of bytes for reals on restart.
      WRITE (unOut, IOSTAT=ErrStat2)   INT(DbKi              ,B4Ki)     ! let's make sure we've got the correct number of bytes for doubles on restart.
      WRITE (unOut, IOSTAT=ErrStat2)   INT(IntKi             ,B4Ki)     ! let's make sure we've got the correct number of bytes for integers on restart.
      WRITE (unOut, IOSTAT=ErrStat2)   AbortErrLev
      WRITE (unOut, IOSTAT=ErrStat2)   NumTurbines                      ! Number of turbines
      WRITE (unOut, IOSTAT=ErrStat2)   t_initial                        ! initial time
      WRITE (unOut, IOSTAT=ErrStat2)   n_t_global                       ! current time step
   
   END IF
      
      
      ! data from current turbine at time step:
   WRITE (unOut, IOSTAT=ErrStat2)   ArraySizes                       ! Number of reals, doubles, and integers written to file
   WRITE (unOut, IOSTAT=ErrStat2)   ReKiBuf                          ! Packed reals
   WRITE (unOut, IOSTAT=ErrStat2)   DbKiBuf                          ! Packed doubles
   WRITE (unOut, IOSTAT=ErrStat2)   IntKiBuf                         ! Packed integers
      
      
   IF ( ALLOCATED(ReKiBuf)  ) DEALLOCATE(ReKiBuf)
   IF ( ALLOCATED(DbKiBuf)  ) DEALLOCATE(DbKiBuf)
   IF ( ALLOCATED(IntKiBuf) ) DEALLOCATE(IntKiBuf)
   
      !CALL FAST_CreateCheckpoint(t_initial, n_t_global, Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
      !            Turbine%ED, Turbine%SrvD, Turbine%AD, Turbine%IfW, &
      !            Turbine%HD, Turbine%SD, Turbine%MAP, Turbine%FEAM, Turbine%MD, &
      !            Turbine%IceF, Turbine%IceD, Turbine%MeshMapData, ErrStat, ErrMsg )              
   
      
   IF (Turbine%TurbID == NumTurbines .OR. .NOT. PRESENT(Unit)) THEN
      CLOSE(unOut)
      unOut = -1
   END IF
   
   IF (PRESENT(Unit)) Unit = unOut
      
      ! A hack to pack Bladed-style DLL data
   IF (Turbine%SrvD%p%UseBladedInterface .AND. Turbine%SrvD%OtherSt%dll_data%avrSWAP( 1) > 0   ) THEN
         ! store value to be overwritten
      old_avrSwap1 = Turbine%SrvD%OtherSt%dll_data%avrSWAP( 1) 
      FileName     = Turbine%SrvD%p%DLL_InFile
         ! overwrite values:
      Turbine%SrvD%p%DLL_InFile = DLLFileName
      Turbine%SrvD%OtherSt%dll_data%avrSWAP(50) = REAL( LEN_TRIM(DLLFileName) ) +1 ! No. of characters in the "INFILE"  argument (-) (we add one for the C NULL CHARACTER)
      Turbine%SrvD%OtherSt%dll_data%avrSWAP( 1) = -8
      CALL CallBladedDLL(Turbine%SrvD%p%DLL_Trgt,  Turbine%SrvD%OtherSt%dll_data, Turbine%SrvD%p, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         ! put values back:
      Turbine%SrvD%p%DLL_InFile = FileName
      Turbine%SrvD%OtherSt%dll_data%avrSWAP(50) = REAL( LEN_TRIM(FileName) ) +1 ! No. of characters in the "INFILE"  argument (-) (we add one for the C NULL CHARACTER)
      Turbine%SrvD%OtherSt%dll_data%avrSWAP( 1) = old_avrSwap1
   END IF
   
   
                  
END SUBROUTINE FAST_CreateCheckpoint_T
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FAST_RestoreFromCheckpoint_Tary(t_initial, n_t_global, Turbine, CheckpointRoot, ErrStat, ErrMsg  )

   REAL(DbKi),               INTENT(IN   ) :: t_initial           ! initial time
   INTEGER(IntKi),           INTENT(  OUT) :: n_t_global          ! loop counter
   TYPE(FAST_TurbineType),   INTENT(  OUT) :: Turbine(:)          ! all data for one instance of a turbine
   CHARACTER(*),             INTENT(IN   ) :: CheckpointRoot      ! Rootname of checkpoint file
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              ! Error message if ErrStat /= ErrID_None

      ! local variables
   REAL(DbKi)                              :: t_initial_out
   INTEGER(IntKi)                          :: NumTurbines_out
   INTEGER(IntKi)                          :: NumTurbines         ! Number of turbines in this simulation
   INTEGER(IntKi)                          :: i_turb
   INTEGER                                 :: Unit
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(1024)                         :: ErrMsg2             ! local error message
   CHARACTER(*),             PARAMETER     :: RoutineName = 'FAST_RestoreFromCheckpoint_Tary' 
   
   
   NumTurbines = SIZE(Turbine)   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
      ! Display the copyright notice
   CALL DispCopyrightLicense( FAST_Ver )

      ! Tell our users what they're running
   CALL WrScr( ' Running '//GetVersion()//NewLine//' linked with '//TRIM( GetNVD( NWTC_Ver ))//NewLine )
   
   Unit = -1         
   DO i_turb = 1,NumTurbines
      CALL FAST_RestoreFromCheckpoint_T(t_initial_out, n_t_global, NumTurbines_out, Turbine(i_turb), CheckpointRoot, ErrStat2, ErrMsg2, Unit )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         IF (t_initial_out /= t_initial) CALL SetErrStat(ErrID_Fatal, "invalid value of t_initial.", ErrStat, ErrMsg, RoutineName )
         IF (NumTurbines_out /= NumTurbines) CALL SetErrStat(ErrID_Fatal, "invalid value of NumTurbines.", ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN
   END DO

   CALL WrScr( ' Restarting simulation at '//TRIM(Num2LStr(n_t_global*Turbine(1)%p_FAST%DT))//' seconds.' )
   
   
END SUBROUTINE FAST_RestoreFromCheckpoint_Tary
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FAST_RestoreFromCheckpoint_T(t_initial, n_t_global, NumTurbines, Turbine, CheckpointRoot, ErrStat, ErrMsg, Unit )
! The inverse of FAST_CreateCheckpoint_T
   USE BladedInterface, ONLY: CallBladedDLL  ! Hack for Bladed-style DLL

   REAL(DbKi),               INTENT(INOUT) :: t_initial           ! initial time
   INTEGER(IntKi),           INTENT(INOUT) :: n_t_global          ! loop counter
   INTEGER(IntKi),           INTENT(INOUT) :: NumTurbines         ! Number of turbines in this simulation
   TYPE(FAST_TurbineType),   INTENT(  OUT) :: Turbine             ! all data for one instance of a turbine
   CHARACTER(*),             INTENT(IN   ) :: CheckpointRoot      ! Rootname of checkpoint file
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              ! Error message if ErrStat /= ErrID_None
   INTEGER(IntKi), OPTIONAL, INTENT(INOUT) :: Unit                ! unit number for output file 
   
      ! local variables:
   REAL(ReKi),               ALLOCATABLE   :: ReKiBuf(:)
   REAL(DbKi),               ALLOCATABLE   :: DbKiBuf(:)
   INTEGER(IntKi),           ALLOCATABLE   :: IntKiBuf(:)
   
   INTEGER(B4Ki)                           :: ArraySizes(3) 
   
   INTEGER(IntKi)                          :: unIn                ! unit number for input file 
   INTEGER(IntKi)                          :: old_avrSwap1        ! previous value of avrSwap(1) !hack for Bladed DLL checkpoint/restore
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(1024)                         :: ErrMsg2             ! local error message
   CHARACTER(*),             PARAMETER     :: RoutineName = 'FAST_RestoreFromCheckpoint_T' 
  
   CHARACTER(1024)                         :: FileName            ! Name of the (input) checkpoint file
   CHARACTER(1024)                         :: DLLFileName         ! Name of the (input) checkpoint file

            
   ErrStat=ErrID_None
   ErrMsg=""
   
   FileName    = TRIM(CheckpointRoot)//'.chkp'
   DLLFileName = TRIM(CheckpointRoot)//'.dll.chkp'
   ! FileName = TRIM(CheckpointRoot)//'.cp'   
   unIn=-1      
   IF (PRESENT(Unit)) unIn = Unit
         
   IF ( unIn < 0 ) THEN

      CALL GetNewUnit( unIn, ErrStat2, ErrMsg2 )
      
      CALL OpenBInpFile ( unIn, FileName, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev ) RETURN
  
         ! checkpoint file header:
      READ (unIn, IOSTAT=ErrStat2)   ArraySizes     ! let's make sure we've got the correct number of bytes for reals, doubles, and integers on restart.
      
      IF ( ArraySizes(1) /= ReKi  ) CALL SetErrStat(ErrID_Fatal,"ReKi on restart is different than when checkpoint file was created.",ErrStat,ErrMsg,RoutineName)
      IF ( ArraySizes(2) /= DbKi  ) CALL SetErrStat(ErrID_Fatal,"DbKi on restart is different than when checkpoint file was created.",ErrStat,ErrMsg,RoutineName)
      IF ( ArraySizes(3) /= IntKi ) CALL SetErrStat(ErrID_Fatal,"IntKi on restart is different than when checkpoint file was created.",ErrStat,ErrMsg,RoutineName)
      IF (ErrStat >= AbortErrLev) THEN
         CLOSE(unIn)
         unIn = -1
         IF (PRESENT(Unit)) Unit = unIn
         RETURN
      END IF
      
      READ (unIn, IOSTAT=ErrStat2)   AbortErrLev
      READ (unIn, IOSTAT=ErrStat2)   NumTurbines                      ! Number of turbines
      READ (unIn, IOSTAT=ErrStat2)   t_initial                        ! initial time
      READ (unIn, IOSTAT=ErrStat2)   n_t_global                       ! current time step
   
   END IF
      
      
      ! data from current time step:
   READ (unIn, IOSTAT=ErrStat2)   ArraySizes                       ! Number of reals, doubles, and integers written to file
   
   ALLOCATE(ReKiBuf( ArraySizes(1)), STAT=ErrStat2)
      IF (ErrStat2 /=0) CALL SetErrStat(ErrID_Fatal, "Could not allocate ReKiBuf", ErrStat, ErrMsg, RoutineName )
   ALLOCATE(DbKiBuf( ArraySizes(2)), STAT=ErrStat2)
      IF (ErrStat2 /=0) CALL SetErrStat(ErrID_Fatal, "Could not allocate DbKiBuf", ErrStat, ErrMsg, RoutineName )
   ALLOCATE(IntKiBuf(ArraySizes(3)), STAT=ErrStat2)
      IF (ErrStat2 /=0) CALL SetErrStat(ErrID_Fatal, "Could not allocate IntKiBuf", ErrStat, ErrMsg, RoutineName )
      IF (ErrStat < AbortErrLev) THEN
      
         READ (unIn, IOSTAT=ErrStat2)   ReKiBuf                          ! Packed reals
         READ (unIn, IOSTAT=ErrStat2)   DbKiBuf                          ! Packed doubles
         READ (unIn, IOSTAT=ErrStat2)   IntKiBuf                         ! Packed integers
         
      END IF

      ! close file if necessary      
   IF (Turbine%TurbID == NumTurbines .OR. .NOT. PRESENT(Unit)) THEN
      CLOSE(unIn)
      unIn = -1
   END IF
   
   IF (PRESENT(Unit)) Unit = unIn

   IF (ErrStat < AbortErrLev) THEN
         ! Put the arrays back in the data types
      CALL FAST_UnpackTurbineType( ReKiBuf, DbKiBuf, IntKiBuf, Turbine, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF
   
   
   IF ( ALLOCATED(ReKiBuf)  ) DEALLOCATE(ReKiBuf)
   IF ( ALLOCATED(DbKiBuf)  ) DEALLOCATE(DbKiBuf)
   IF ( ALLOCATED(IntKiBuf) ) DEALLOCATE(IntKiBuf)    
   
   
      ! A sort-of hack to restore MAP DLL data (in particular Turbine%MAP%OtherSt%C_Obj%object)
    ! these must be the same variables that are used in MAP_Init because they get allocated in the DLL and
    ! destroyed in MAP_End
   IF (Turbine%p_FAST%CompMooring == Module_MAP) THEN
      CALL MAP_Restart( Turbine%MAP%Input(1), Turbine%MAP%p, Turbine%MAP%x(STATE_CURR), Turbine%MAP%xd(STATE_CURR), &
                        Turbine%MAP%z(STATE_CURR), Turbine%MAP%OtherSt, Turbine%MAP%y, ErrStat2, ErrMsg2 )   
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                           
   END IF
   
   
      ! A hack to restore Bladed-style DLL data
   IF (Turbine%SrvD%p%UseBladedInterface .AND. Turbine%SrvD%OtherSt%dll_data%avrSWAP( 1) > 0   ) THEN
         ! store value to be overwritten
      old_avrSwap1 = Turbine%SrvD%OtherSt%dll_data%avrSWAP( 1) 
      FileName     = Turbine%SrvD%p%DLL_InFile
         ! overwrite values before calling DLL:
      Turbine%SrvD%p%DLL_InFile = DLLFileName
      Turbine%SrvD%OtherSt%dll_data%avrSWAP(50) = REAL( LEN_TRIM(DLLFileName) ) +1 ! No. of characters in the "INFILE"  argument (-) (we add one for the C NULL CHARACTER)
      Turbine%SrvD%OtherSt%dll_data%avrSWAP( 1) = -9
      CALL CallBladedDLL(Turbine%SrvD%p%DLL_Trgt,  Turbine%SrvD%OtherSt%dll_data, Turbine%SrvD%p, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                           
         ! put values back:
      Turbine%SrvD%p%DLL_InFile = FileName
      Turbine%SrvD%OtherSt%dll_data%avrSWAP(50) = REAL( LEN_TRIM(FileName) ) +1 ! No. of characters in the "INFILE"  argument (-) (we add one for the C NULL CHARACTER)
      Turbine%SrvD%OtherSt%dll_data%avrSWAP( 1) = old_avrSwap1            
   END IF   
   
      ! deal with sibling meshes here:
      

END SUBROUTINE FAST_RestoreFromCheckpoint_T
!----------------------------------------------------------------------------------------------------------------------------------

END MODULE FAST_Subs
