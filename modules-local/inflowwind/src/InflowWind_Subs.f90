!**********************************************************************************************************************************
! $Id: InflowWind.f90 125 2014-10-29 22:28:35Z aplatt $
!
! This module is used to read and process the (undisturbed) inflow winds.  It must be initialized
! using InflowWind_Init() with the name of the file, the file type, and possibly reference height and
! width (depending on the type of wind file being used).  This module calls appropriate routines
! in the wind modules so that the type of wind becomes seamless to the user.  InflowWind_End()
! should be called when the program has finshed.
!
! Data are assumed to be in units of meters and seconds.  Z is measured from the ground (NOT the hub!).
!
!  7 Oct 2009    Initial Release with AeroDyn 13.00.00         B. Jonkman, NREL/NWTC
! 14 Nov 2011    v1.00.01b-bjj                                 B. Jonkman
!  1 Aug 2012    v1.01.00a-bjj                                 B. Jonkman
! 10 Aug 2012    v1.01.00b-bjj                                 B. Jonkman
!    Feb 2013    v2.00.00a-adp   conversion to Framework       A. Platt
!
!..................................................................................................................................
! Files with this module:
!!!!!  InflowWind_Subs.f90
!  InflowWind.txt       -- InflowWind_Types will be auto-generated based on the descriptions found in this file.
!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
!
!    This file is part of InflowWind.
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
! File last committed: $Date: 2014-10-29 16:28:35 -0600 (Wed, 29 Oct 2014) $
! (File) Revision #: $Rev: 125 $
! URL: $HeadURL: https://windsvn.nrel.gov/InflowWind/branches/modularization2/Source/InflowWind.f90 $
!**********************************************************************************************************************************
MODULE InflowWind_Subs


   USE                              InflowWind_Types
   USE                              NWTC_Library

      !-------------------------------------------------------------------------------------------------
      ! The included wind modules
      !-------------------------------------------------------------------------------------------------

   USE                              IfW_UniformWind_Types      ! Types for IfW_UniformWind
   USE                              IfW_UniformWind            ! uniform wind files (text files)
!!!   USE                              IfW_FFWind_Types           ! Types for IfW_FFWind
!!!   USE                              IfW_FFWind                 ! full-field binary wind files
!!!   USE                              HAWCWind                   ! full-field binary wind files in HAWC format
!!!   USE                              FDWind                     ! 4-D binary wind files
!!!   USE                              CTWind                     ! coherent turbulence from KH billow - binary file superimposed on another wind type
!!!   USE                              UserWind                   ! user-defined wind module




   IMPLICIT NONE
   PRIVATE


      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: InflowWind_ReadInput            ! Read IfW input file
   PUBLIC :: InflowWind_ValidateInput        ! Check data read from file
   PUBLIC :: InflowWind_SetParameters        ! Store input file data in Parameters
   PUBLIC :: PrintBadChannelWarning          ! For requested output channels


CONTAINS


!====================================================================================================
SUBROUTINE InflowWind_ReadInput( InputFileName, EchoFileName, InputFileData, ErrStat, ErrMsg )
!     This public subroutine reads the input required for InflowWind from the file whose name is an
!     input parameter.
!----------------------------------------------------------------------------------------------------


      ! Passed variables
   CHARACTER(*),                       INTENT(IN   )  :: InputFileName
   CHARACTER(*),                       INTENT(IN   )  :: EchoFileName
   TYPE(InflowWind_InputFile),         INTENT(INOUT)  :: InputFileData            !< The data for initialization
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat              !< Returned error status  from this subroutine
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg               !< Returned error message from this subroutine


      ! Local variables
   INTEGER(IntKi)                                     :: UnitInput            !< Unit number for the input file
   INTEGER(IntKi)                                     :: UnitEcho             !< The local unit number for this module's echo file
   CHARACTER(1024)                                    :: TmpPath              !< Temporary storage for relative path name
   CHARACTER(1024)                                    :: TmpFmt               !< Temporary storage for format statement
   CHARACTER(35)                                      :: Frmt                 !< Output format for logical parameters. (matches NWTC Subroutine Library format)


      ! Temoporary messages
   INTEGER(IntKi)                                     :: TmpErrStat
   CHARACTER(LEN(ErrMsg))                             :: TmpErrMsg


      ! Initialize local data

   UnitEcho                = -1
   Frmt                    = "( 2X, L11, 2X, A, T30, ' - ', A )"
   ErrStat                 = ErrID_None
   ErrMsg                  = ""
   InputFileData%EchoFlag  = .FALSE.  ! initialize for error handling (cleanup() routine)



   !-------------------------------------------------------------------------------------------------
   ! Open the file
   !-------------------------------------------------------------------------------------------------

   CALL GetNewUnit( UnitInput, TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )

   CALL OpenFInpFile( UnitInput, TRIM(InputFileName), TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF


      ! Tell the users what we are doing... though this part should be rather quick.
   CALL WrScr( ' Opening InflowWind input file:  '//InputFileName )


   !-------------------------------------------------------------------------------------------------
   ! File header
   !-------------------------------------------------------------------------------------------------

   CALL ReadCom( UnitInput, InputFileName, 'InflowWind input file header line 1', TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

   CALL ReadCom( UnitInput, InputFileName, 'InflowWind input file header line 2', TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

   CALL ReadCom( UnitInput, InputFileName, 'InflowWind input file separator line', TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF


     ! Echo Input Files.

   CALL ReadVar ( UnitInput, InputFileName, InputFileData%EchoFlag, 'Echo', 'Echo Input', TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

      ! If we are Echoing the input then we should re-read the first three lines so that we can echo them
      ! using the NWTC_Library routines.  The echoing is done inside those routines via a global variable
      ! which we must store, set, and then replace on error or completion.

   IF ( InputFileData%EchoFlag ) THEN

      CALL OpenEcho ( UnitEcho, TRIM(EchoFileName), TmpErrStat, TmpErrMsg )
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

      REWIND(UnitInput)


         ! The input file was already successfully read through up to this point, so we shouldn't have any read
         ! errors in the first four lines.  So, we won't worry about checking the error status here.

      CALL ReadCom( UnitInput, InputFileName, 'InflowWind input file header line 1', TmpErrStat, TmpErrMsg, UnitEcho )
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )

      CALL ReadCom( UnitInput, InputFileName, 'InflowWind input file header line 2', TmpErrStat, TmpErrMsg, UnitEcho )
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )

      CALL ReadCom( UnitInput, InputFileName, 'InflowWind input file separator line', TmpErrStat, TmpErrMsg, UnitEcho )
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )

         ! Echo Input Files.

      CALL ReadVar ( UnitInput, InputFileName, InputFileData%EchoFlag, 'Echo', 'Echo the input file data', TmpErrStat, TmpErrMsg, UnitEcho )
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )

   END IF



   !-------------------------------------------------------------------------------------------------
   !> Read general section with wind type, direction, and output point list (applies to all wind types)
   !-------------------------------------------------------------------------------------------------


      ! Read WindType
   CALL ReadVar( UnitInput, InputFileName, InputFileData%WindType, 'WindType', &
               'switch for wind file type (1=steady; 2=uniform; 3=binary TurbSim FF; '//&
               '4=binary Bladed-style FF; 5=HAWC format; 6=User defined)', &
               TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput')
   IF (ErrStat >= AbortErrLev) THEN
      CALL CleanUp()
      RETURN
   ENDIF


      ! Read PropogationDir
   CALL ReadVar( UnitInput, InputFileName, InputFileData%PropogationDir, 'PropogationDir', &
               'Direction of wind propogation (meteoroligical direction)', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput')
   IF (ErrStat >= AbortErrLev) THEN
      CALL CleanUp()
      RETURN
   ENDIF


      ! Read the number of points for the wind velocity output
   CALL ReadVar( UnitInput, InputFileName, InputFileData%NWindVel, 'NWindVel', &
               'Number of points to output the wind velocity (0 to 9)', &
               TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput')
   IF (ErrStat >= AbortErrLev) THEN
      CALL CleanUp()
      RETURN
   ENDIF

      ! Before proceeding, make sure that NWindVel makes sense
   IF ( InputFileData%NWindVel < 0 .OR. InputFileData%NwindVel > 9 ) THEN
      CALL SetErrStat( ErrID_Fatal, 'NWindVel must be greater than or equal to zero and less than 10.', &
                        ErrStat, ErrMsg, 'InflowWind_ReadInput' )
      CALL CleanUp()
      RETURN
   ELSE

      ! Allocate space for the output location arrays:
      CALL AllocAry( InputFileData%WindVxiList, InputFileData%NWindVel, 'WindVxiList', TmpErrStat, TmpErrMsg )
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
      CALL AllocAry( InputFileData%WindVyiList, InputFileData%NWindVel, 'WindVyiList', TmpErrStat, TmpErrMsg )
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
      CALL AllocAry( InputFileData%WindVziList, InputFileData%NWindVel, 'WindVziList', TmpErrStat, TmpErrMsg )
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      ENDIF
   ENDIF

      ! Read in the values of WindVxiList
   CALL ReadAry( UnitInput, InputFileName, InputFileData%WindVxiList, InputFileData%NWindVel, 'WindVxiList', &
               'List of coordinates in the inertial X direction (m)', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput')
   IF (ErrStat >= AbortErrLev) THEN
      CALL CleanUp()
      RETURN
   ENDIF

      ! Read in the values of WindVxiList
   CALL ReadAry( UnitInput, InputFileName, InputFileData%WindVyiList, InputFileData%NWindVel, 'WindVyiList', &
               'List of coordinates in the inertial Y direction (m)', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput')
   IF (ErrStat >= AbortErrLev) THEN
      CALL CleanUp()
      RETURN
   ENDIF

      ! Read in the values of WindVziList
   CALL ReadAry( UnitInput, InputFileName, InputFileData%WindVziList, InputFileData%NWindVel, 'WindVziList', &
               'List of coordinates in the inertial Z direction (m)', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput')
   IF (ErrStat >= AbortErrLev) THEN
      CALL CleanUp()
      RETURN
   ENDIF


   !-------------------------------------------------------------------------------------------------
   !> Read the _Parameters for Steady Wind Conditions [used only for WindType = 1]_ section
   !-------------------------------------------------------------------------------------------------

      ! Section separator line
   CALL ReadCom( UnitInput, InputFileName, 'InflowWind input file separator line', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF


      ! Read HWindSpeed
   CALL ReadVar( UnitInput, InputFileName, InputFileData%Steady_HWindSpeed, 'HWindSpeed', &
                  'Horizontal windspeed for steady wind', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput')
   IF (ErrStat >= AbortErrLev) THEN
      CALL CleanUp()
      RETURN
   ENDIF

      ! Read RefHt
   CALL ReadVar( UnitInput, InputFileName, InputFileData%Steady_RefHt, 'RefHt', &
                  'Reference height for horizontal wind speed for steady wind', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput')
   IF (ErrStat >= AbortErrLev) THEN
      CALL CleanUp()
      RETURN
   ENDIF

      ! Read PLexp
   CALL ReadVar( UnitInput, InputFileName, InputFileData%Steady_PLexp, 'PLexp', &
                  'Power law exponent for steady wind', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput')
   IF (ErrStat >= AbortErrLev) THEN
      CALL CleanUp()
      RETURN
   ENDIF


   !-------------------------------------------------------------------------------------------------
   !> Read the _Parameters for Uniform wind file [used only for WindType = 2]_ section
   !-------------------------------------------------------------------------------------------------

      ! Section separator line
   CALL ReadCom( UnitInput, InputFileName, 'InflowWind input file separator line', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

      ! Read UniformWindFile
   CALL ReadVar( UnitInput, InputFileName, InputFileData%Uniform_FileName, 'WindFileName', &
                  'Filename of time series data for uniform wind field', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput')
   IF (ErrStat >= AbortErrLev) THEN
      CALL CleanUp()
      RETURN
   ENDIF

      ! Read RefHt
   CALL ReadVar( UnitInput, InputFileName, InputFileData%Uniform_RefHt, 'RefHt', &
                  'Reference height for uniform wind file', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput')
   IF (ErrStat >= AbortErrLev) THEN
      CALL CleanUp()
      RETURN
   ENDIF

      ! Read RefLength
   CALL ReadVar( UnitInput, InputFileName, InputFileData%Uniform_RefLength, 'RefLength', &
                  'Reference length for uniform wind file', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput')
   IF (ErrStat >= AbortErrLev) THEN
      CALL CleanUp()
      RETURN
   ENDIF


   !-------------------------------------------------------------------------------------------------
   !> Read the _Parameters for Binary TurbSim Full-Field files [used only for WindType = 3]_ section
   !-------------------------------------------------------------------------------------------------

      ! Section separator line
   CALL ReadCom( UnitInput, InputFileName, 'InflowWind input file separator line', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

      ! Read TSFFWind info
   CALL ReadVar( UnitInput, InputFileName, InputFileData%TSFF_FileName, 'FileName', &
               'Name of the TurbSim full field wind file to use (.bts)', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput')
   IF (ErrStat >= AbortErrLev) THEN
      CALL CleanUp()
      RETURN
   ENDIF


   !-------------------------------------------------------------------------------------------------
   !> Read the _Parameters for Binary Bladed-style Full-Field files [used only for WindType = 4]_ section
   !-------------------------------------------------------------------------------------------------

      ! Section separator line
   CALL ReadCom( UnitInput, InputFileName, 'InflowWind input file separator line', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

      ! Read BladedStyle%WindFileName
   CALL ReadVar( UnitInput, InputFileName, InputFileData%Bladed_FileName, 'FileName', &
               'Name of the TurbSim full field wind file to use (.bts)', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput')
   IF (ErrStat >= AbortErrLev) THEN
      CALL CleanUp()
      RETURN
   ENDIF

      ! Read TowerFileFlag
   CALL ReadVar( UnitInput, InputFileName, InputFileData%Bladed_TowerFile, 'TowerFileFlag', &
               'Have tower file (.twr) [flag]', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput')
   IF (ErrStat >= AbortErrLev) THEN
      CALL CleanUp()
      RETURN
   ENDIF


   !-------------------------------------------------------------------------------------------------
   !> Read the _Parameters for coherent turbulence [used only for WindType = 3 or 4]_ section
   !-------------------------------------------------------------------------------------------------

      ! Section separator line
   CALL ReadCom( UnitInput, InputFileName, 'InflowWind input file separator line', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

      ! Read CTTS_Flag
   CALL ReadVar( UnitInput, InputFileName, InputFileData%CTTS_CoherentTurb, 'CTTS_CoherentTurbFlag', &
               'Flag to coherent turbulence', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput')
   IF (ErrStat >= AbortErrLev) THEN
      CALL CleanUp()
      RETURN
   ENDIF

      ! Read CTWind%WindFileName
   CALL ReadVar( UnitInput, InputFileName, InputFileData%CTTS_FileName, 'CTTS_FileName', &
               'Name of coherent turbulence file', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput')
   IF (ErrStat >= AbortErrLev) THEN
      CALL CleanUp()
      RETURN
   ENDIF

      ! Read CTWind%PathName
   CALL ReadVar( UnitInput, InputFileName, InputFileData%CTTS_Path, 'CTTS_Path', &
               'Path to coherent turbulence binary files', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput')
   IF (ErrStat >= AbortErrLev) THEN
      CALL CleanUp()
      RETURN
   ENDIF



   !-------------------------------------------------------------------------------------------------
   !> Read the _Parameters for HAWC-formatted binary files [used only for WindType = 5]_ section
   !-------------------------------------------------------------------------------------------------

      ! Section separator line
   CALL ReadCom( UnitInput, InputFileName, 'InflowWind input file separator line', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

      ! Read HAWC_FileName_u
   CALL ReadVar( UnitInput, InputFileName, InputFileData%HAWC_FileName_u, 'HAWC_FileName_u', &
               'Name of the file containing the u-component fluctuating wind', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

      ! Read HAWC_FileName_v
   CALL ReadVar( UnitInput, InputFileName, InputFileData%HAWC_FileName_v, 'HAWC_FileName_v', &
               'Name of the file containing the v-component fluctuating wind', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

      ! Read HAWC_FileName_w
   CALL ReadVar( UnitInput, InputFileName, InputFileData%HAWC_FileName_w, 'HAWC_FileName_w', &
               'Name of the file containing the w-component fluctuating wind', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

      ! Read HAWC_nx
   CALL ReadVar( UnitInput, InputFileName, InputFileData%HAWC_nx, 'HAWC_nx', &
               'Number of grids in the x direction (in the 3 files above)', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

      ! Read HAWC_ny
   CALL ReadVar( UnitInput, InputFileName, InputFileData%HAWC_ny, 'HAWC_ny', &
               'Number of grids in the y direction (in the 3 files above)', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

      ! Read HAWC_nz
   CALL ReadVar( UnitInput, InputFileName, InputFileData%HAWC_nz, 'HAWC_nz', &
               'Number of grids in the z direction (in the 3 files above)', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

      ! Read HAWC_dx
   CALL ReadVar( UnitInput, InputFileName, InputFileData%HAWC_dx, 'HAWC_dx', &
               'Number of grids in the x direction (in the 3 files above)', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

      ! Read HAWC_dy
   CALL ReadVar( UnitInput, InputFileName, InputFileData%HAWC_dy, 'HAWC_dy', &
               'Number of grids in the y direction (in the 3 files above)', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

      ! Read HAWC_dz
   CALL ReadVar( UnitInput, InputFileName, InputFileData%HAWC_dz, 'HAWC_dz', &
               'Number of grids in the z direction (in the 3 files above)', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

      ! Read HAWC_RefHt
   CALL ReadVar( UnitInput, InputFileName, InputFileData%HAWC_RefHt, 'HAWC_RefHt', &
               'Reference (hub) height of the grid', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF



      !----------------------------------------------------------------------------------------------
      !> Read the _Scaling parameters for turbulence (HAWC-format files) [used only for WindType = 5]_ subsection
      !----------------------------------------------------------------------------------------------

      ! Section separator line
   CALL ReadCom( UnitInput, InputFileName, 'InflowWind input file separator line', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

      ! Read HAWC_ScaleMethod
   CALL ReadVar( UnitInput, InputFileName, InputFileData%HAWC_ScaleMethod, 'HAWC_ScaleMethod', &
               'Turbulence scaling method [0=none, 1=direct scaling, 2= calculate scaling '// &
               'factor based on a desired standard deviation]', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

      ! Read HAWC_SFx
   CALL ReadVar( UnitInput, InputFileName, InputFileData%HAWC_SFx, 'HAWC_SFx', &
               'Turbulence scaling factor for the x direction [ScaleMethod=1]', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

      ! Read HAWC_SFy
   CALL ReadVar( UnitInput, InputFileName, InputFileData%HAWC_SFy, 'HAWC_SFy', &
               'Turbulence scaling factor for the y direction [ScaleMethod=1]', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

      ! Read HAWC_SFz
   CALL ReadVar( UnitInput, InputFileName, InputFileData%HAWC_SFz, 'HAWC_SFz', &
               'Turbulence scaling factor for the z direction [ScaleMethod=1]', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

      ! Read HAWC_SigmaFx
   CALL ReadVar( UnitInput, InputFileName, InputFileData%HAWC_SigmaFx, 'HAWC_SigmaFx', &
               'Turbulence standard deviation to calculate scaling from in x direction [ScaleMethod=2]', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

      ! Read HAWC_SigmaFy
   CALL ReadVar( UnitInput, InputFileName, InputFileData%HAWC_SigmaFy, 'HAWC_SigmaFy', &
               'Turbulence standard deviation to calculate scaling from in y direction [ScaleMethod=2]', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

      ! Read HAWC_SigmaFz
   CALL ReadVar( UnitInput, InputFileName, InputFileData%HAWC_SigmaFz, 'HAWC_SigmaFz', &
               'Turbulence standard deviation to calculate scaling from in z direction [ScaleMethod=2]', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

!FIXME:  TStart has no comment
      ! Read HAWC_TStart
   CALL ReadVar( UnitInput, InputFileName, InputFileData%HAWC_TStart, 'HAWC_TStart', &
               '', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

!FIXME:  TEnd has no comment
      ! Read HAWC_TEnd
   CALL ReadVar( UnitInput, InputFileName, InputFileData%HAWC_TEnd, 'HAWC_TEnd', &
               '', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF


      !----------------------------------------------------------------------------------------------
      !> Read the _Mean wind profile paramters (added to HAWC-format files) [used only for WindType = 5]_ subsection
      !----------------------------------------------------------------------------------------------

      ! Section separator line
   CALL ReadCom( UnitInput, InputFileName, 'InflowWind input file separator line', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

      ! Read HAWC_URef
   CALL ReadVar( UnitInput, InputFileName, InputFileData%HAWC_URef, 'HAWC_URef', &
               'Mean u-component wind speed at the reference height', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

      ! Read HAWC_ProfileType
   CALL ReadVar( UnitInput, InputFileName, InputFileData%HAWC_ProfileType, 'HAWC_ProfileType', &
               'Wind profile type ("LOG"=logarithmic, "PL"=power law, or "UD"=user defined)', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

      ! Read HAWC_PLExp
   CALL ReadVar( UnitInput, InputFileName, InputFileData%HAWC_PLExp, 'HAWC_PLExp', &
               'Power law exponent (used for PL wind profile type only)', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

      ! Read HAWC_Z0
   CALL ReadVar( UnitInput, InputFileName, InputFileData%HAWC_Z0, 'HAWC_Z0', &
               'Surface roughness length (used for LOG wind profile type only)', TmpErrStat, TmpErrMsg, UnitEcho )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'InflowWind_ReadInput' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF




!FIXME: read the outlist


   !-------------------------------------------------------------------------------------------------
   ! This is the end of the input file
   !-------------------------------------------------------------------------------------------------

   CALL Cleanup()

   RETURN

      CONTAINS
         !..............................
         SUBROUTINE Cleanup()


               ! Close input file
            CLOSE ( UnitInput )

               ! Cleanup the Echo file and global variables
            IF ( InputFileData%EchoFlag ) THEN
               CLOSE(UnitEcho)
            END IF


         END SUBROUTINE Cleanup

END SUBROUTINE InflowWind_ReadInput


!====================================================================================================
!> This private subroutine verifies the input required for InflowWind is correctly specified.  This
!! routine checkes all the parameters that are common with all the wind types, then calls subroutines
!! that check the parameters specific to each wind type.  Only the parameters corresponding to the
!! desired wind type are evaluated; the rest are ignored.  Additional checks will be performed after
!! the respective wind file has been read in, but these checks will be performed within the respective
!! wind module.
!
!  The reason for structuring it this way is to allow for relocating the validation routines for the
!  wind type into their respective modules. It might also prove useful later if we change languages
!  but retain the fortran wind modules.
SUBROUTINE InflowWind_ValidateInput( InputFileData, ErrStat, ErrMsg )


      ! Passed variables

   TYPE(InflowWind_InputFile),         INTENT(INOUT)  :: InputFileData        !< The data for initialization
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat              !< Error status  from this subroutine
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg               !< Error message from this subroutine


      ! Temporary variables
   INTEGER(IntKi)                                     :: TmpErrStat           !< Temporary error status  for subroutine and function calls
   CHARACTER(LEN(ErrMsg))                             :: TmpErrMsg            !< Temporary error message for subroutine and function calls

      ! Local variables



      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""


      !-----------------------------------------------
      ! Data that applies to all wind types
      !-----------------------------------------------

      ! WindType
   IF ( InputFileData%WindType <= Undef_WindNumber .OR. InputFileData%WindType > Highest_WindNumber ) THEN
      CALL SetErrStat(ErrID_Fatal,' Invalid WindType specified.  Only values of '//   &
            TRIM(Num2LStr( Steady_WindNumber    ))//','// &
            TRIM(Num2LStr( Uniform_WindNumber   ))//','// &
            TRIM(Num2LStr( TSFF_WindNumber      ))//','// &
            TRIM(Num2LStr( BladedFF_WindNumber  ))//','// &
            TRIM(Num2LStr( HAWC_WindNumber      ))//', or '// &
            TRIM(Num2LStr( User_WindNumber      ))//' are supported.', &
            ErrStat,ErrMsg,'InflowWind_ValidateInput')
      RETURN
   ENDIF




      !-----------------------------------------------
      ! Data specific to a WindType
      !-----------------------------------------------

   SELECT CASE ( InputFileData%WindType )

      CASE ( Steady_WindNumber )
         CALL Steady_ValidateInput()

      CASE ( Uniform_WindNumber )
         CALL Uniform_ValidateInput()

      CASE ( TSFF_WindNumber )
         CALL TSFF_ValidateInput()

      CASE ( BladedFF_WindNumber )
         CALL BladedFF_ValidateInput()

      CASE ( HAWC_WindNumber )
         CALL HAWC_ValidateInput()

      CASE ( User_WindNumber )
         CALL User_ValidateInput()

      CASE DEFAULT  ! keep this check to make sure that all new wind types have been accounted for
         CALL SetErrStat(ErrID_Fatal,' Undefined wind type.',ErrStat,ErrMsg,'InflowWind_ValidateInput')

   END SELECT

      ! NOTE: we are not checking the error status yet, but rather will run through the Coherent turbulence section
      !        before returning.  This way a user will know of errors there as well before having to fix the first
      !        section.

      ! Coherent turbulence
   IF ( InputFileData%CTTS_CoherentTurb ) THEN
      CALL CTTS_ValidateInput()
   ENDIF


   RETURN

CONTAINS

   !> This subroutine checks that the values for the steady wind make sense.  There aren't many to check.
   SUBROUTINE Steady_ValidateInput()

         ! Check that HWindSpeed is positive
      IF ( InputFileData%Steady_HWindSpeed <= 0.0_ReKi ) THEN
         CALL SetErrStat(ErrID_Fatal,' Horizontal wind speed (HWindSpeed) for steady winds must be greater than zero.',  &
               ErrStat,ErrMsg,'Steady_ValidateInput')
      ENDIF

         ! Check that RefHt is positive
      IF ( InputFileData%Steady_RefHt <= 0.0_ReKi ) THEN
         CALL SetErrStat(ErrID_Fatal,' Reference height (RefHt) for steady winds must be greater than zero.',  &
               ErrStat,ErrMsg,'Steady_ValidateInput')
      ENDIF

         ! Check that PLexp is positive
      IF ( InputFileData%Steady_PLexp < 0.0_ReKi ) THEN
         CALL SetErrStat(ErrID_Fatal,' Power law exponent (PLexp) for steady winds must be positive or zero.',  &
               ErrStat,ErrMsg,'Steady_ValidateInput')
      ENDIF

      RETURN
   END SUBROUTINE Steady_ValidateInput



   !> This subroutine checks that the values for the uniform wind make sense.  There isn't much to check.
   SUBROUTINE Uniform_ValidateInput()

         ! Local variables
      LOGICAL                                         :: TmpFileExist

         ! Check that the filename requested actually exists
      INQUIRE( file=InputFileData%Uniform_FileName, exist=TmpFileExist )
      IF ( .NOT. TmpFileExist ) THEN
         CALL SetErrStat( ErrID_Fatal," Cannot find uniform wind input file: '"//TRIM(InputFileData%Uniform_FileName)//"'", &
               ErrStat,ErrMsg,'Uniform_ValidateInput')
      ENDIF

         ! Check that RefHt is positive
      IF ( InputFileData%Uniform_RefHt <= 0.0_ReKi ) THEN
         CALL SetErrStat(ErrID_Fatal,' Reference height (RefHt) for uniform winds must be greater than zero.',  &
               ErrStat,ErrMsg,'Uniform_ValidateInput')
      ENDIF

!FIXME: Do I need to check the RefLength?
      RETURN

   END SUBROUTINE Uniform_ValidateInput



   !> There isn't much to check for the TurbSim full-field section of the input file.
   SUBROUTINE TSFF_ValidateInput()

         ! Local variables
      LOGICAL                                         :: TmpFileExist

         ! Check that the filename requested actually exists
      INQUIRE( file=InputFileData%TSFF_FileName, exist=TmpFileExist )
      IF ( .NOT. TmpFileExist ) THEN
         CALL SetErrStat( ErrID_Fatal," Cannot find TurbSim full-field wind input file: '"//TRIM(InputFileData%TSFF_FileName)//"'", &
               ErrStat,ErrMsg,'TSFF_ValidateInput')
      ENDIF

      RETURN

   END SUBROUTINE TSFF_ValidateInput



   !> We will only check the main bladed input file here.  We will check the existance of the towerfile (if requested) in the SetParameters routine.
   SUBROUTINE BladedFF_ValidateInput()

         ! Local variables
      LOGICAL                                         :: TmpFileExist

         ! Check that the filename requested actually exists
      INQUIRE( file=InputFileData%Bladed_FileName, exist=TmpFileExist )
      IF ( .NOT. TmpFileExist ) THEN
         CALL SetErrStat( ErrID_Fatal," Cannot find Bladed-style full-field wind input file: '"//TRIM(InputFileData%Bladed_FileName)//"'", &
               ErrStat,ErrMsg,'BladedFF_ValidateInput')
      ENDIF

      RETURN
   END SUBROUTINE BladedFF_ValidateInput



   !> This routine checks to see if coherent turbulence can be used or not.  We check the path and file
   !! information before proceeding.  The coherent turbulence info is loaded after the full-field files have
   !! been read, and those can be very slow to load.  So we preemptively check the existence of the coherent
   !! turbulence file.  Further checking that the coherent turbulence can actually be used is check by the
   !! CTTS_Wind submodule (after the loading of the other stuff).
   SUBROUTINE CTTS_ValidateInput()

         ! Local variables
      LOGICAL                                         :: TmpFileExist

         ! Check if coherent turbulence can be used or not.  It is only applicable to the TSFF and BladedFF wind
         ! file types.
      IF ( InputFileData%CTTS_CoherentTurb .AND. ( .NOT. InputFileData%WindType == TSFF_WindNumber ) &
               .AND. ( .NOT. InputFileData%WindType == BladedFF_WindNumber ) ) THEN
         CALL SetErrStat(ErrID_Fatal,' Coherent turbulence may only be used with full field wind data (WindType = '// &
               TRIM(Num2LStr(TSFF_WindNumber))//' or '//TRIM(Num2LStr(BladedFF_WindNumber))//').', &
               ErrStat,ErrMsg,'CTTS_ValidateInput')
         RETURN
      ENDIF

         ! Check that the CT input file exists.
      INQUIRE( file=InputFileData%CTTS_FileName, exist=TmpFileExist )
      IF ( .NOT. TmpFileExist ) THEN
         CALL SetErrStat( ErrID_Fatal," Cannot find the coherent turbulence input file: '"//TRIM(InputFileData%CTTS_FileName)//"'", &
               ErrStat,ErrMsg,'CTTS_ValidateInput')
      ENDIF

         ! Check that the CT Event path exists.
      ! Unfortunately there is no good portable method to do this, so we won't try.
!FIXME: See if we have anything in the library.  -- Nothing I found...

      RETURN

   END SUBROUTINE CTTS_ValidateInput



   !> This routine checks the HAWC wind file type section of the input file.
   SUBROUTINE HAWC_ValidateInput()

         ! Local variables
      LOGICAL                                         :: TmpFileExist

         ! Check that the wind files exist.
      INQUIRE( file=InputFileData%HAWC_FileName_u, exist=TmpFileExist )
      IF ( .NOT. TmpFileExist ) THEN
         CALL SetErrStat( ErrID_Fatal," Cannot find the HAWC u-component wind file: '"//TRIM(InputFileData%HAWC_FileName_u)//"'", &
               ErrStat,ErrMsg,'HAWC_ValidateInput')
      ENDIF

      INQUIRE( file=InputFileData%HAWC_FileName_v, exist=TmpFileExist )
      IF ( .NOT. TmpFileExist ) THEN
         CALL SetErrStat( ErrID_Fatal," Cannot find the HAWC v-component wind file: '"//TRIM(InputFileData%HAWC_FileName_v)//"'", &
               ErrStat,ErrMsg,'HAWC_ValidateInput')
      ENDIF
      INQUIRE( file=InputFileData%HAWC_FileName_w, exist=TmpFileExist )
      IF ( .NOT. TmpFileExist ) THEN
         CALL SetErrStat( ErrID_Fatal," Cannot find the HAWC w-component wind file: '"//TRIM(InputFileData%HAWC_FileName_w)//"'", &
               ErrStat,ErrMsg,'HAWC_ValidateInput')
      ENDIF


         ! Check that the number of grids make some sense
      IF ( InputFileData%HAWC_nx <= 0_IntKi ) THEN
         CALL SetErrStat( ErrID_Fatal,' Number of grids in the x-direction of HAWC wind files (nx) must be greater than zero.',  &
               ErrStat,ErrMsg,'HAWC_ValidateInput')
      ENDIF

      IF ( InputFileData%HAWC_ny <= 0_IntKi ) THEN
         CALL SetErrStat( ErrID_Fatal,' Number of grids in the y-direction of HAWC wind files (ny) must be greater than zero.',  &
               ErrStat,ErrMsg,'HAWC_ValidateInput')
      ENDIF

      IF ( InputFileData%HAWC_nz <= 0_IntKi ) THEN
         CALL SetErrStat( ErrID_Fatal,' Number of grids in the z-direction of HAWC wind files (nz) must be greater than zero.',  &
               ErrStat,ErrMsg,'HAWC_ValidateInput')
      ENDIF


         ! Check that the distance between points is positive
      IF ( InputFileData%HAWC_dx <= 0.0_ReKi ) THEN
         CALL SetErrStat( ErrID_Fatal,' Distance between points in the x-direction of HAWC wind files (dx) must be greater than zero.',  &
               ErrStat,ErrMsg,'HAWC_ValidateInput')
      ENDIF

      IF ( InputFileData%HAWC_dy <= 0.0_ReKi ) THEN
         CALL SetErrStat( ErrID_Fatal,' Distance between points in the y-direction of HAWC wind files (dy) must be greater than zero.',  &
               ErrStat,ErrMsg,'HAWC_ValidateInput')
      ENDIF

      IF ( InputFileData%HAWC_dz <= 0.0_ReKi ) THEN
         CALL SetErrStat( ErrID_Fatal,' Distance between points in the z-direction of HAWC wind files (dz) must be greater than zero.',  &
               ErrStat,ErrMsg,'HAWC_ValidateInput')
      ENDIF


         ! Check that RefHt is positive
      IF ( InputFileData%HAWC_RefHt <= 0.0_ReKi ) THEN
         CALL SetErrStat( ErrID_Fatal,' Reference height (RefHt) for HAWC winds must be greater than zero.',  &
               ErrStat,ErrMsg,'HAWC_ValidateInput')
      ENDIF


         !----------------------
         ! Scaling method info
         !----------------------

      IF ( ( InputFileData%HAWC_ScaleMethod < 0_IntKi ) .OR. ( InputFileData%HAWC_ScaleMethod > 3_IntKi ) ) THEN
         CALL SetErrStat( ErrID_Fatal,' The scaling method (ScaleMethod) for HAWC winds can only be 0, 1, or 2.', &
               ErrStat,ErrMsg,'HAWC_ValidateInput')
      ENDIF

         ! Check the other scaling parameters
      SELECT CASE( InputFileData%HAWC_ScaleMethod )

      CASE ( 1 )

         IF ( InputFileData%HAWC_SFx <= 0.0_ReKi ) THEN
            CALL SetErrStat( ErrID_Fatal,' HAWC wind scaling paramter SFx must be greater than zero.', &
                  ErrStat,ErrMsg,'HAWC_ValidateInput')
         ENDIF

         IF ( InputFileData%HAWC_SFy <= 0.0_ReKi ) THEN
            CALL SetErrStat( ErrID_Fatal,' HAWC wind scaling paramter SFy must be greater than zero.', &
                  ErrStat,ErrMsg,'HAWC_ValidateInput')
         ENDIF

         IF ( InputFileData%HAWC_SFz <= 0.0_ReKi ) THEN
            CALL SetErrStat( ErrID_Fatal,' HAWC wind scaling paramter SFz must be greater than zero.', &
                  ErrStat,ErrMsg,'HAWC_ValidateInput')
         ENDIF

      CASE ( 2 )

         IF ( InputFileData%HAWC_SigmaFx <= 0.0_ReKi ) THEN
            CALL SetErrStat( ErrID_Fatal,' HAWC wind scaling paramter SigmaFx must be greater than zero.', &
                  ErrStat,ErrMsg,'HAWC_ValidateInput')
         ENDIF

         IF ( InputFileData%HAWC_SigmaFy <= 0.0_ReKi ) THEN
            CALL SetErrStat( ErrID_Fatal,' HAWC wind scaling paramter SigmaFy must be greater than zero.', &
                  ErrStat,ErrMsg,'HAWC_ValidateInput')
         ENDIF

         IF ( InputFileData%HAWC_SigmaFz <= 0.0_ReKi ) THEN
            CALL SetErrStat( ErrID_Fatal,' HAWC wind scaling paramter SigmaFz must be greater than zero.', &
                  ErrStat,ErrMsg,'HAWC_ValidateInput')
         ENDIF

      END SELECT


         ! Start and end times of the scaling, but only if ScaleMethod is set.
      IF ( ( InputFileData%HAWC_ScaleMethod == 1_IntKi ) .OR. (InputFileData%HAWC_ScaleMethod == 2_IntKi ) ) THEN
            ! Check that the start time is >= 0.  NOTE: this may be an invalid test.
         IF ( InputFileData%HAWC_TStart < 0.0_ReKi ) THEN
            CALL SetErrStat( ErrID_Fatal,' HAWC wind scaling TStart must be zero or positive.',  &
                  ErrStat,ErrMsg,'HAWC_ValidateInput')
         ENDIF

         IF ( InputFileData%HAWC_TStart < 0.0_ReKi ) THEN
            CALL SetErrStat( ErrID_Fatal,' HAWC wind scaling TStart must be zero or positive.',  &
                  ErrStat,ErrMsg,'HAWC_ValidateInput')
         ENDIF

         IF ( InputFileData%HAWC_TStart > InputFileData%HAWC_TEnd ) THEN
            CALL SetErrStat( ErrID_Fatal,' HAWC wind scaling start time must occur before the end time.',   &
                  ErrStat,ErrMsg,'HAWC_ValidateInput')
         ENDIF

!FIXME:  How do we want to handle having the start and end times both being exactly zero?  Is that a do not apply, or apply for all time?
!FIXME:  What about TStart == TEnd ?
         IF ( EqualRealNos( InputFileData%HAWC_TStart, InputFileData%HAWC_TEnd ) ) THEN
            CALL SetErrStat( ErrID_Severe,' Start and end time for HAWC wind scaling are identical.  No time elapses.', &
                  ErrStat,ErrMsg,'HAWC_ValidateInput')
         ENDIF

      ENDIF

      RETURN

   END SUBROUTINE HAWC_ValidateInput



   !> This routine is a placeholder for later development by the end user.  It should only be used if additional inputs
   !! for the user defined routines are added to the input file. Otherwise, this routine may be ignored by the user.
   SUBROUTINE User_ValidateInput()
      RETURN
   END SUBROUTINE User_ValidateInput


END SUBROUTINE InflowWind_ValidateInput





!====================================================================================================
!> This private subroutine copies the info from the input file over to the parameters for InflowWind.
SUBROUTINE InflowWind_SetParameters( InputFileData, ParamData, ErrStat, ErrMsg )


      ! Passed variables

   TYPE(InflowWind_InputFile),         INTENT(INOUT)  :: InputFileData        !< The data for initialization
   TYPE(InflowWind_ParameterType),     INTENT(INOUT)  :: ParamData            !< The parameters for InflowWind
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat              !< Error status  from this subroutine
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg               !< Error message from this subroutine


      ! Temporary variables
   INTEGER(IntKi)                                     :: TmpErrStat           !< Temporary error status  for subroutine and function calls
   CHARACTER(LEN(ErrMsg))                             :: TmpErrMsg            !< Temporary error message for subroutine and function calls

      ! Local variables



      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""


   !-----------------------------------------------------------------
   ! Copy over the general information that applies to everything
   !-----------------------------------------------------------------

      ! Copy the WindType over.

   ParamData%WindType   =  InputFileData%WindType


      ! Convert the PropogationDir to radians and store this.  For simplicity, we will shift it to be between -pi and pi

   ParamData%PropogationDir   = D2R * InputFileData%PropogationDir
   CALL MPi2Pi( ParamData%PropogationDir )         ! Shift if necessary so that the value is between -pi and pi


      ! Copy over the list of wind coordinates.  Move the arrays to the new one.

   ParamData%NWindVel   =  InputFileData%NWindVel
   CALL MOVE_ALLOC( InputFileData%WindVxiList, ParamData%WindVxiList )
   CALL MOVE_ALLOC( InputFileData%WindVyiList, ParamData%WindVyiList )
   CALL MOVE_ALLOC( InputFileData%WindVziList, ParamData%WindVziList )



   SELECT CASE ( ParamData%WindType )


      CASE ( Steady_WindNumber )
         CALL Steady_SetParameters()

      CASE ( Uniform_WindNumber )
         CALL Uniform_SetParameters()

      CASE ( TSFF_WindNumber )
         CALL TSFF_SetParameters()

      CASE ( BladedFF_WindNumber )
         CALL BladedFF_SetParameters()

      CASE ( HAWC_WindNumber )
         CALL HAWC_SetParameters()

      CASE ( User_WindNumber )
         CALL User_SetParameters()

      CASE DEFAULT  ! keep this check to make sure that all new wind types have been accounted for
         CALL SetErrStat(ErrID_Fatal,' Undefined wind type.',ErrStat,ErrMsg,'InflowWind_SetParameters')

   END SELECT


   IF ( InputFileData%CTTS_CoherentTurb ) THEN
      CALL CTTS_SetParameters()
   ENDIF



CONTAINS
   SUBROUTINE Steady_SetParameters()
!      Uniform_InitData%ReferenceHeight =  InputFileData%Uniform_RefHt
!      Uniform_InitData%Width           =  InputFileData%Uniform_RefLength 
!      Uniform_InitData%WindFileName    =  InputFileData%Uniform_FileName
!      CALL SetErrStat(ErrID_Warn,' This subroutine has not been written yet.',ErrStat,ErrMsg,'Steady_SetParameters')
   END SUBROUTINE Steady_SetParameters

   SUBROUTINE Uniform_SetParameters()
      ! No parameters to set here.  This is done in the Uniform_Init routine.
   END SUBROUTINE Uniform_SetParameters

   SUBROUTINE TSFF_SetParameters()
      CALL SetErrStat(ErrID_Warn,' This subroutine has not been written yet.',ErrStat,ErrMsg,'TSFF_SetParameters')
   END SUBROUTINE TSFF_SetParameters

   SUBROUTINE BladedFF_SetParameters()
!Check the tower file status here.
      CALL SetErrStat(ErrID_Warn,' This subroutine has not been written yet.',ErrStat,ErrMsg,'BladedFF_SetParameters')
   END SUBROUTINE BladedFF_SetParameters

   SUBROUTINE HAWC_SetParameters()
      CALL SetErrStat(ErrID_Warn,' This subroutine has not been written yet.',ErrStat,ErrMsg,'HAWC_SetParameters')
   END SUBROUTINE HAWC_SetParameters

   SUBROUTINE User_SetParameters()
      CALL SetErrStat(ErrID_Warn,' This subroutine has not been written yet.',ErrStat,ErrMsg,'User_SetParameters')
   END SUBROUTINE User_SetParameters

   SUBROUTINE CTTS_SetParameters()
      CALL SetErrStat(ErrID_Warn,' This subroutine has not been written yet.',ErrStat,ErrMsg,'CTTS_SetParameters')
   END SUBROUTINE CTTS_SetParameters

END SUBROUTINE InflowWind_SetParameters




!====================================================================================================
!> The routine prints out warning messages if the user has requested invalid output channel names
!! The errstat is set to ErrID_Warning if any element in foundMask is .FALSE.
!FIXME: revamp this routine.
SUBROUTINE PrintBadChannelWarning(NUserOutputs, UserOutputs , foundMask, ErrStat, ErrMsg )
!----------------------------------------------------------------------------------------------------
   INTEGER(IntKi),                     INTENT(IN   )  :: NUserOutputs         !< Number of user-specified output channels
   CHARACTER(10),                      INTENT(IN   )  :: UserOutputs(:)       !< An array holding the names of the requested output channels.
   LOGICAL,                            INTENT(IN   )  :: foundMask(:)         !< A mask indicating whether a user requested channel belongs to a module's output channels.
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat              !< returns a non-zero value when an error occurs
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg               !< Error message if ErrStat /= ErrID_None


  INTEGER(IntKi)                                      :: I                    !< Generic loop counter

   ErrStat = ErrID_None
   ErrMsg  = ''

   DO I = 1, NUserOutputs
      IF (.NOT. foundMask(I)) THEN
         ErrMsg  = ' A requested output channel is invalid'
         CALL ProgWarn( 'The requested output channel is invalid: ' // UserOutputs(I) )
         ErrStat = ErrID_Warn
      END IF
   END DO



END SUBROUTINE PrintBadChannelWarning






!====================================================================================================
END MODULE InflowWind_Subs
