!>  This module uses full-field binary wind files to determine the wind inflow.
!!  This module assumes that the origin, (0,0,0), is located at the tower centerline at ground level,
!!  and that all units are specified in the metric system (using meters and seconds).
!!  Data is shifted by half the grid width to account for turbine yaw (so that data in the X
!!  direction actually starts at -1*ParamData%FFYHWid meters).
MODULE IfW_BladedFFWind
!!
!!  Created 25-Sep-2009 by B. Jonkman, National Renewable Energy Laboratory
!!     using subroutines and modules from AeroDyn v12.58
!!
!!----------------------------------------------------------------------------------------------------
!!  Feb 2013    v2.00.00          A. Platt
!!     -- updated to the new framework
!!     -- Modified to use NWTC_Library v. 2.0
!!     -- Note:  Jacobians are not included in this version.
!!
!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
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

   USE                                          NWTC_Library
   USE                                          IfW_BladedFFWind_Types
   USE                                          IfW_FFWind_Base

   IMPLICIT                                     NONE
   PRIVATE

   TYPE(ProgDesc),   PARAMETER               :: IfW_BladedFFWind_Ver = ProgDesc( 'IfW_BladedFFWind', '', '' )

   PUBLIC                                    :: IfW_BladedFFWind_Init
   PUBLIC                                    :: IfW_BladedFFWind_End
   PUBLIC                                    :: IfW_BladedFFWind_CalcOutput




CONTAINS
!====================================================================================================
!>  This routine is used read the full-field turbulence data.
!!  09/25/1997  - Created by M. Buhl from GETFILES in ViewWind.
!!  09/23/2009  - modified by B. Jonkman: this subroutine was split into several subroutines (was ReadFF)
!!  16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
SUBROUTINE IfW_BladedFFWind_Init(InitInp, ParamData, MiscVars, InitOutData, ErrStat, ErrMsg)

   CHARACTER(*),              PARAMETER                        :: RoutineName="IfW_BladedFFWind_Init"

      ! Passed Variables
   TYPE(IfW_BladedFFWind_InitInputType),        INTENT(IN   )  :: InitInp        !< Initialization data passed to the module
   TYPE(IfW_BladedFFWind_ParameterType),        INTENT(  OUT)  :: ParamData      !< Parameters
   TYPE(IfW_BladedFFWind_MiscVarType),          INTENT(  OUT)  :: MiscVars       !< misc/optimization data   (storage for the main data)
   TYPE(IfW_BladedFFWind_InitOutputType),       INTENT(  OUT)  :: InitOutData    !< Initial output


      ! Error Handling
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat        !< determines if an error has been encountered
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg         !< Message about errors


      ! Temporary variables for error handling
   REAL(ReKi)                                            :: TI(3)
   REAL(ReKi)                                            :: ScaleFactors(3)
   INTEGER(IntKi)                                        :: TmpErrStat     ! temporary error status
   CHARACTER(ErrMsgLen)                                  :: TmpErrMsg      ! temporary error message
   TYPE(IfW_FFWind_InitInputType)                        :: FF_InitInp     ! Initialization input data for FF scaling

   

   ErrMsg      = ''
   ErrStat     = ErrID_None
   
   
   CALL ReadFiles(InitInp, FF_InitInp, InitOutData, ParamData, TI, TmpErrStat, TmpErrMsg)
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) RETURN
   
      !-------------------------------------------------------------------------------------------------
      ! If the wind file has zero-mean and unit standard deviation (native Bladed format), scale the data:
      !-------------------------------------------------------------------------------------------------
   ParamData%FF%AddMeanAfterInterp = .false.
   ParamData%FF%WindProfileType    = FF_InitInp%WindProfileType
   ParamData%FF%Z0                 = FF_InitInp%Z0
   ParamData%FF%PLExp              = FF_InitInp%PLExp

   if (InitInp%NativeBladedFmt) then
      ParamData%FF%InterpTower = .true.
      
         ! Validate scaling data if we've got native-Bladed format
      CALL FFWind_ValidateInput(FF_InitInp, ParamData%FF%NFFComp, TmpErrStat, TmpErrMsg)
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN
         
         ! scale to requested TI (or use requested scale factors)
      call ScaleTurbulence(FF_InitInp, ParamData%FF%FFData(:,:,:,1:ParamData%FF%NFFSteps), ScaleFactors, TmpErrStat, TmpErrMsg)      
         CALL SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName) 
         IF (ErrStat >= AbortErrLev) RETURN      
      
         ! Add the mean wind speed to the u component.
      call AddMeanVelocity(FF_InitInp, ParamData%FF%GridBase, 1.0_ReKi/ParamData%FF%InvFFZD, ParamData%FF%FFData)
   else
      ParamData%FF%InterpTower = .false.
   end if
   

   IF (ParamData%FF%Periodic) THEN
      ParamData%FF%InitXPosition = 0                ! start at the hub
      ParamData%FF%TotalTime     = ParamData%FF%NFFSteps*ParamData%FF%FFDTime
   ELSE
      ParamData%FF%InitXPosition = ParamData%FF%FFYHWid          ! start half the grid width ahead of the turbine
      ParamData%FF%TotalTime     = (ParamData%FF%NFFSteps-1)*ParamData%FF%FFDTime
   ENDIF

      ! overwrite the offset
   IF (InitInp%NativeBladedFmt) THEN
      ParamData%FF%InitXPosition = FF_InitInp%XOffset
   END IF
   
   
      !-------------------------------------------------------------------------------------------------
      ! Set the InitOutput information
      !-------------------------------------------------------------------------------------------------

   InitOutdata%Ver         =  IfW_BladedFFWind_Ver
   InitOutdata%TI          =  TI



      !-------------------------------------------------------------------------------------------------
      ! Write to the summary file
      !-------------------------------------------------------------------------------------------------

   IF ( InitInp%SumFileUnit > 0 ) THEN
      WRITE(InitInp%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)
      WRITE(InitInp%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)    'Bladed-style wind type.  Read by InflowWind sub-module '//   &
                                                                     TRIM(IfW_BladedFFWind_Ver%Name)//' '//TRIM(IfW_BladedFFWind_Ver%Ver)
      WRITE(InitInp%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)    TRIM(TmpErrMsg)
      WRITE(InitInp%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)    '     FileName:                    '//TRIM(InitInp%WindFileName)
      WRITE(InitInp%SumFileUnit,'(A34,I3)',   IOSTAT=TmpErrStat)    '     Binary file format id:       ',ParamData%FF%WindFileFormat
      WRITE(InitInp%SumFileUnit,'(A34,G12.4)',IOSTAT=TmpErrStat)    '     Reference height (m):        ',ParamData%FF%RefHt
      WRITE(InitInp%SumFileUnit,'(A34,G12.4)',IOSTAT=TmpErrStat)    '     Timestep (s):                ',ParamData%FF%FFDTime
      WRITE(InitInp%SumFileUnit,'(A34,I12)',  IOSTAT=TmpErrStat)    '     Number of timesteps:         ',ParamData%FF%NFFSteps
      WRITE(InitInp%SumFileUnit,'(A34,G12.4)',IOSTAT=TmpErrStat)    '     Mean windspeed (m/s):        ',ParamData%FF%MeanFFWS
      WRITE(InitInp%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)    '     Characteristic TI:            [ '// &
                                                TRIM(Num2LStr(TI(1)))//', '//TRIM(Num2LStr(TI(2)))//', '//TRIM(Num2LStr(TI(3)))//' ] '
      WRITE(InitInp%SumFileUnit,'(A34,L1)',   IOSTAT=TmpErrStat)    '     Windfile is periodic:        ',ParamData%FF%Periodic
      WRITE(InitInp%SumFileUnit,'(A34,L1)',   IOSTAT=TmpErrStat)    '     Windfile includes tower:     ',ParamData%FF%NTGrids > 0

      IF ( ParamData%FF%Periodic ) THEN
         WRITE(InitInp%SumFileUnit,'(A)',     IOSTAT=TmpErrStat)    '     Time range (s):              [ '// &
                     TRIM(Num2LStr(0.0_ReKi))//' : '//TRIM(Num2LStr(ParamData%FF%TotalTime))//' ]'
      ELSE  ! Shift the time range to compensate for the shifting of the wind grid
         WRITE(InitInp%SumFileUnit,'(A)',     IOSTAT=TmpErrStat)    '     Time range (s):              [ '// &
                     TRIM(Num2LStr(-ParamData%FF%InitXPosition*ParamData%FF%InvMFFWS))//' : '// &
                     TRIM(Num2LStr(ParamData%FF%TotalTime-ParamData%FF%InitXPosition*ParamData%FF%InvMFFWS))//' ]'
      ENDIF

      WRITE(InitInp%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)    '     Y range (m):                 [ '// &
                     TRIM(Num2LStr(-ParamData%FF%FFYHWid))//' : '//TRIM(Num2LStr(ParamData%FF%FFYHWid))//' ]'

      IF ( ParamData%FF%NTGrids > 0 ) THEN
         WRITE(InitInp%SumFileUnit,'(A)',     IOSTAT=TmpErrStat)    '     Z range (m):                 [ '// &
                     TRIM(Num2LStr(0.0_ReKi))//' : '//TRIM(Num2LStr(ParamData%FF%RefHt + ParamData%FF%FFZHWid))//' ]'
      ELSE
         WRITE(InitInp%SumFileUnit,'(A)',     IOSTAT=TmpErrStat)    '     Z range (m):                 [ '// &
                     TRIM(Num2LStr(ParamData%FF%RefHt - ParamData%FF%FFZHWid))//' : '//TRIM(Num2LStr(ParamData%FF%RefHt + ParamData%FF%FFZHWid))//' ]'
      ENDIF


         ! We are assuming that if the last line was written ok, then all of them were.
      IF (TmpErrStat /= 0_IntKi) THEN
         CALL SetErrStat(ErrID_Fatal,'Error writing to summary file.',ErrStat,ErrMsg,RoutineName)
         RETURN
      ENDIF   
   ENDIF 


      
 
   RETURN

END SUBROUTINE IfW_BladedFFWind_Init
!========================================================================================================
SUBROUTINE ReadFiles(InitInp, FF_InitInp, InitOut, ParamData, TI, ErrStat, ErrMsg)

   CHARACTER(*),              PARAMETER                        :: RoutineName="ReadFiles"

      ! Passed Variables
   TYPE(IfW_BladedFFWind_InitInputType),        INTENT(IN   )  :: InitInp        !< Initialization data passed to the module
   TYPE(IfW_FFWind_InitInputType),              INTENT(  OUT)  :: FF_InitInp     !< Initialization data for scaling
   TYPE(IfW_BladedFFWind_InitOutputType),       INTENT(INOUT)  :: InitOut        !< Initial output
   TYPE(IfW_BladedFFWind_ParameterType),        INTENT(  OUT)  :: ParamData      !< Parameters
   REAL(ReKi)                           ,       INTENT(  OUT)  :: TI      (3)    !< turbulence intensities of the wind components as defined in the FF file, not necessarially the actual TI


      ! Error Handling
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat        !< determines if an error has been encountered
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg         !< Message about errors


      ! Temporary variables for error handling
   INTEGER(IntKi)                                        :: TmpErrStat     ! temporary error status
   CHARACTER(ErrMsgLen)                                  :: TmpErrMsg      ! temporary error message


      ! Local Variables:

   REAL(ReKi)                                            :: BinTI   (3)    ! turbulence intensities of the wind components as defined in the FF binary file, not necessarially the actual TI
   REAL(ReKi)                                            :: NatTI   (3)    ! turbulence intensities of the wind components as defined in the native FF summary file
   REAL(ReKi)                                            :: UBar
   REAL(ReKi)                                            :: ZCenter
 
   INTEGER(IntKi)                                        :: UnitWind      ! Unit number for the InflowWind input file
   INTEGER(B2Ki)                                         :: Dum_Int2
   INTEGER(IntKi)                                        :: I
   LOGICAL                                               :: CWise
   LOGICAL                                               :: LHR            ! Left-hand rule for Bladed files (is the v component aligned along *negative* Y?)
   
   LOGICAL                                               :: Exists
   CHARACTER( 1028 )                                     :: SumFile        ! length is LEN(ParamData%WindFileName) + the 4-character extension.
   CHARACTER( 1028 )                                     :: TwrFile        ! length is LEN(ParamData%WindFileName) + the 4-character extension.

   CHARACTER(1024)                                       :: BinFileName 
   CHARACTER(1024)                                       :: PriPath 

   
   ErrMsg      = ''
   ErrStat     = ErrID_None
   
   
   if (InitInp%NativeBladedFmt) then
      call Read_NativeBladedSummary(InitInp%WindFileName, FF_InitInp%PLExp, NatTI, ParamData%FF%MeanFFWS, ParamData%FF%RefHt, InitOut%PropagationDir, InitOut%VFlowAngle, BinFileName, FF_InitInp%XOffset, TmpErrStat, TmpErrMsg)
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN
      
      
      if (pathIsRelative(BinFileName)) then
         CALL GetPath( InitInp%WindFileName, PriPath )     ! Binary file will be relative to the path where the primary input file is located.
         BinFileName = TRIM(PriPath)//TRIM(BinFileName)
      end if
      
      IF ( InitInp%FixedWindFileRootName ) THEN ! .TRUE. when FAST.Farm uses multiple instances of InflowWind for ambient wind data
         IF ( InitInp%TurbineID == 0 ) THEN     ! .TRUE. for the FAST.Farm low-resolution domain
            BinFileName = TRIM(BinFileName)//TRIM(PathSep)//'Low'
         ELSE                                   ! FAST.Farm high-resolution domain(s)
            BinFileName = TRIM(BinFileName)//TRIM(PathSep)//'HighT'//TRIM(Num2Lstr(InitInp%TurbineID))
         ENDIF
      ENDIF
         
         ! default values for Bladed Format
      CWise = .false.
      ZCenter = ParamData%FF%RefHt
      ParamData%FF%Periodic = .true.
      
      FF_InitInp%ScaleMethod = ScaleMethod_StdDev
      FF_InitInp%SigmaF = NatTI * ParamData%FF%MeanFFWS
      FF_InitInp%sf     = FF_InitInp%SigmaF
!      FF_InitInp%ScaleMethod = ScaleMethod_Direct ! Bladed files should have std of 1, so we'll just multiply (closer to what Bladed does)
      
      FF_InitInp%RefHt = ParamData%FF%RefHt
      FF_InitInp%URef = ParamData%FF%MeanFFWS
      FF_InitInp%WindProfileType = WindProfileType_PL ! it could also have logarithmic, but I'm going to leave that off for now
      
      TI   = 100.0_ReKi
      UBar =   0.0_ReKi
      LHR  = .true.
            
   else
      InitOut%PropagationDir = 0.0_ReKi
      InitOut%VFlowAngle = 0.0_ReKi
      BinFileName = InitInp%WindFileName
   end if
      
   
      ! Get a unit number to use

   CALL GetNewUnit(UnitWind, TmpErrStat, TmpErrMsg)
   CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
   IF ( ErrStat >= AbortErrLev ) RETURN


      !----------------------------------------------------------------------------------------------
      ! Open the binary file, read its "header" (first 2-byte integer) to determine what format
      ! binary file it is, and close it.
      !----------------------------------------------------------------------------------------------

   CALL OpenBInpFile (UnitWind, TRIM(BinFileName), TmpErrStat, TmpErrMsg)
   CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
   IF ( ErrStat >= AbortErrLev ) RETURN

      ! Read the first binary integer from the file to get info on the type.
      ! Cannot use library read routines since this is a 2-byte integer.
   READ ( UnitWind, IOSTAT=TmpErrStat )  Dum_Int2
   CLOSE( UnitWind )

   IF (TmpErrStat /= 0) THEN
      CALL SetErrStat(ErrID_Fatal,' Error reading first binary integer from file "'//TRIM(BinFileName)//'."',   &
               ErrStat,ErrMsg,RoutineName)
      RETURN
   ENDIF


      !----------------------------------------------------------------------------------------------
      ! Read the files to get the required FF data.
      !----------------------------------------------------------------------------------------------

      ! Store the binary format information so the InflowWind code can use it.
      ! Also changes to IntKi from INT(2) to compare in the SELECT below
   ParamData%FF%WindFileFormat = Dum_Int2


   SELECT CASE (ParamData%FF%WindFileFormat)

      CASE ( -1, -2, -3, -99 )                                         ! Bladed-style binary format

         IF (.not. InitInp%NativeBladedFmt) THEN
                  
               !...........................................................................................
               ! Create full-field summary file name from binary file root name.  Also get tower file
               ! name.
               !...........................................................................................

            CALL GetRoot(BinFileName, SumFile)

            TwrFile = TRIM(SumFile)//'.twr'
            SumFile = TRIM(SumFile)//'.sum'


               !...........................................................................................
               ! Read the summary file to get necessary scaling information
               !...........................................................................................
            
            CALL Read_Summary_FF (UnitWind, TRIM(SumFile), CWise, ZCenter, TI, UBar, ParamData%FF%RefHt, ParamData%FF%Periodic, LHR, TmpErrStat, TmpErrMsg )
               CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )                  
               IF ( ErrStat >= AbortErrLev ) THEN
                  CLOSE ( UnitWind )
                  RETURN
               END IF  
               
         END IF
         

         !...........................................................................................
         ! Open the binary file and read its header
         !...........................................................................................

            CALL OpenBInpFile (UnitWind, TRIM(BinFileName), TmpErrStat, TmpErrMsg )
               CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )                  
               IF ( ErrStat >= AbortErrLev ) THEN
                  CLOSE ( UnitWind )
                  RETURN
               END IF
            
            IF ( Dum_Int2 == -99 ) THEN                                                      ! Newer-style BLADED format
               CALL Read_Bladed_FF_Header1 (UnitWind, BinTI, ParamData%FF, InitInp%NativeBladedFmt, TmpErrStat, TmpErrMsg)
                  CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )                  
                  IF ( ErrStat >= AbortErrLev ) THEN
                     CLOSE ( UnitWind )
                     RETURN
                  END IF

                  ! If the TIs are also in the binary file (BinTI > 0),
                  ! use those numbers instead of ones from the summary file

               if (.not. InitInp%NativeBladedFmt) then
                  DO I =1,ParamData%FF%NFFComp
                     IF ( BinTI(I) > 0 ) TI(I) = BinTI(I)
                  ENDDO
               end if
               

            ELSE
               CALL Read_Bladed_FF_Header0 (UnitWind, ParamData%FF, InitInp%NativeBladedFmt, TmpErrStat, TmpErrMsg)     ! Older-style BLADED format
                  CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )                  
                  IF ( ErrStat >= AbortErrLev ) THEN
                     CLOSE ( UnitWind )
                     RETURN
                  END IF

            ENDIF



         !...........................................................................................
         ! Let's see if the summary and binary FF wind files go together before continuing.
         !...........................................................................................

            IF (.not. InitInp%NativeBladedFmt) THEN            
               IF ( ABS( UBar - ParamData%FF%MeanFFWS ) > 0.1 )  THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error: Incompatible mean hub-height wind speeds in FF wind files. '//&
                              '(Check that the .sum and .wnd files were generated together.)', ErrStat, ErrMsg, RoutineName )                  
                  CLOSE ( UnitWind )
                  RETURN
               ENDIF

            END IF
            
         !...........................................................................................
         ! Calculate the height of the bottom of the grid
         !...........................................................................................

            ParamData%FF%GridBase = ZCenter - ParamData%FF%FFZHWid         ! the location, in meters, of the bottom of the grid
            IF ( ParamData%FF%GridBase < 0.0_ReKi ) THEN
               call SetErrStat( ErrID_Severe, 'WARNING: The bottom of the grid is located at a height of '//&
                               TRIM( Num2LStr(ParamData%FF%GridBase) )//' meters, which is below the ground.'//&
                      ' Winds below the ground will be set to 0.', ErrStat,ErrMsg, RoutineName)
            END IF

         !...........................................................................................
         ! Read the binary grids (converted to m/s) and close the file
         !...........................................................................................

            CALL Read_Bladed_Grids( UnitWind, InitInp%NativeBladedFmt, CWise, LHR, TI, ParamData%FF, TmpErrStat, TmpErrMsg)
               CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )                  
               
            CLOSE ( UnitWind )
            if (InitInp%NativeBladedFmt) TI = NatTI*100.0_ReKi  ! report these TI for the native Bladed format in percent

            IF ( ErrStat >= AbortErrLev ) RETURN
            
         !...........................................................................................
         ! Read the tower points file
         !...........................................................................................

            IF ( InitInp%TowerFileExist .AND. .NOT. InitInp%NativeBladedFmt) THEN      ! If we specified a tower file
               INQUIRE ( FILE=TRIM(TwrFile) , EXIST=Exists )

                  ! Double check that the tower file exists and read it.  If it was requested but doesn't exist,
                  ! throw fatal error and exit.
               IF (  Exists )  THEN
                  CALL Read_FF_Tower( UnitWind, ParamData%FF, TwrFile, TmpErrStat, TmpErrMsg )
                  CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )                  
                  IF ( ErrStat >= AbortErrLev ) THEN
                     CLOSE ( UnitWind )
                     RETURN
                  END IF
               ELSE
                  CALL SetErrStat( ErrID_Fatal, ' Tower file '//TRIM(TwrFile)//' specified for Bladed full-field '// &
                           'wind files does not exist.', ErrStat, ErrMsg, RoutineName)
                  CLOSE ( UnitWind )
                  RETURN
               ENDIF
            ELSE
               ParamData%FF%NTGrids  =  0_IntKi
            ENDIF


      CASE DEFAULT
         CALL SetErrStat( ErrID_Fatal, ' This is not a bladed-style binary wind file (binary format identifier: '//  &
                  TRIM(Num2LStr(ParamData%FF%WindFileFormat))//'.  This might be a TurbSim binary wind file.', &
                  ErrStat, ErrMsg, RoutineName )
         RETURN

      END SELECT
      
END SUBROUTINE ReadFiles

   !====================================================================================================
   !> This subroutine reads the text summary file to get normalizing parameters, the location of the
   !! grid, and the direction the grid was written to the binary file
   !!
   !!   16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   SUBROUTINE Read_Summary_FF ( UnitWind, FileName, CWise, ZCenter, TI, UBar, RefHt, Periodic, LHR, ErrStat, ErrMsg )

      IMPLICIT                                              NONE

      CHARACTER(*),        PARAMETER                     :: RoutineName="Read_Summary_FF"


         ! Passed variables
      INTEGER(IntKi),                     INTENT(IN   )  :: UnitWind       !< unit number for the file to open
      CHARACTER(*),                       INTENT(IN   )  :: FileName       !< name of the summary file
      LOGICAL,                            INTENT(  OUT)  :: CWise          !< rotation (for reading the order of the binary data)
      REAL(ReKi),                         INTENT(  OUT)  :: ZCenter        !< the height at the center of the grid
      REAL(ReKi),                         INTENT(  OUT)  :: TI      (3)    !< turbulence intensities of the wind components as defined in the FF file, not necessarially the actual TI
      REAL(ReKi),                         INTENT(  OUT)  :: UBar           !< mean (advection) wind speed
      REAL(ReKi),                         INTENT(  OUT)  :: RefHt          !< Reference height
      LOGICAL,                            INTENT(  OUT)  :: Periodic       !< rotation (for reading the order of the binary data)
      LOGICAL,                            INTENT(  OUT)  :: LHR            !< Left-hand rule for Bladed files (is the v component aligned along *negative* Y?)
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat        !< returns 0 if no error encountered in the subroutine
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg         !< holds the error messages

        ! Local variables
      REAL(ReKi)                                         :: ZGOffset       ! The vertical offset of the turbine on rectangular grid (allows turbulence not centered on turbine hub)

      INTEGER, PARAMETER                                 :: NumStrings = 7 ! number of strings to be looking for in the file

      INTEGER(IntKi)                                     :: FirstIndx      ! The first character of a line where data is located
      INTEGER(IntKi)                                     :: I              ! A loop counter
      INTEGER(IntKi)                                     :: LastIndx       ! The last  character of a line where data is located
      INTEGER(IntKi)                                     :: LineCount      ! Number of lines that have been read in the file

      LOGICAL                                            :: StrNeeded(NumStrings)   ! if the string has been found

      CHARACTER(1024)                                    :: LINE           ! temporary storage for reading a line from the file

         ! Temporary variables for error handling
      INTEGER(IntKi)                                     :: TmpErrStat     ! temporary error status
      CHARACTER(ErrMsgLen)                               :: TmpErrMsg      ! temporary error message

         !----------------------------------------------------------------------------------------------
         ! Initialize some variables
         !----------------------------------------------------------------------------------------------

      ErrStat              = ErrID_None
      ErrMsg               = ''

      LineCount            = 0
      StrNeeded(:)         = .TRUE.
      ZGOffset             = 0.0
      RefHt                = 0.0
      Periodic             = .FALSE.
      LHR                  = .FALSE.

         !----------------------------------------------------------------------------------------------
         ! Open summary file.
         !----------------------------------------------------------------------------------------------

      CALL OpenFInpFile ( UnitWind, TRIM( FileName ), TmpErrStat, TmpErrMsg )
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )                  
      IF ( ErrStat >= AbortErrLev ) RETURN


         !----------------------------------------------------------------------------------------------
         ! Read the summary file.
         !----------------------------------------------------------------------------------------------

      ! Here are the strings we're looking for, in this order:
      ! 1) 'CLOCKWISE'
      ! 2) 'HUB HEIGHT'
      ! 3)     (unused; decided we didn't need to read data also stored in the binary file)
      ! 4) 'UBAR'
      ! 5) 'HEIGHT OFFSET' (optional)
      ! 6) 'PERIODIC' (optional)
      ! 7) 'BLADED LEFT-HAND RULE' (optional)
         
         
      DO WHILE ( ( ErrStat == ErrID_None ) .AND. StrNeeded(NumStrings) )

         LineCount = LineCount + 1

         READ ( UnitWind, '(A)', IOSTAT=TmpErrStat ) LINE
         IF ( TmpErrStat /= 0 ) THEN

               ! the "HEIGHT OFFSET", "PERIODIC", and "BLADED LEFT-HAND RULE" parameters are not necessary.  We'll assume they are zero/false if we didn't find it.
            IF ( StrNeeded(1) .OR. StrNeeded(2) .OR. StrNeeded(4)  ) THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading line #'//TRIM(Num2LStr(LineCount))//' of the summary file, "'// &
                           TRIM(FileName)//'". Could not find all of the required parameters.', ErrStat, ErrMsg, RoutineName )                  
               RETURN
            ELSE
               EXIT
            ENDIF

         ENDIF

         CALL Conv2UC ( LINE )


         IF ( StrNeeded(1) ) THEN

            !-------------------------------------------------------------------------------------------
            ! #1: Get the rotation direction, using the string "CLOCKWISE"
            !-------------------------------------------------------------------------------------------

            IF ( INDEX( LINE, 'CLOCKWISE' ) > 0 ) THEN

               READ (LINE, *, IOSTAT = TmpErrStat)  CWise          ! Look for True/False values

               IF ( TmpErrStat /= 0 ) THEN                         ! Look for Yes/No values instead

                  LINE = ADJUSTL ( LINE )                      ! Remove leading spaces from input line

                  SELECT CASE (LINE(1:1) )
                     CASE ('Y')
                        CWise = .TRUE.
                     CASE ('N')
                        CWise = .FALSE.
                     CASE DEFAULT
                        CALL SetErrStat( ErrID_Fatal, ' Error reading rotation direction (CLOCKWISE) from FF summary file.', ErrStat, ErrMsg, RoutineName )                  
                        RETURN
                  END SELECT

               ENDIF ! TmpErrStat /= 0
               StrNeeded(1) = .FALSE.

            ENDIF   ! INDEX for "CLOCKWISE"

         ELSEIF ( StrNeeded(2) ) THEN

            !-------------------------------------------------------------------------------------------
            ! #2: Get the hub height, using the strings "HUB HEIGHT" or "ZHUB"
            !-------------------------------------------------------------------------------------------

            IF ( INDEX( LINE, 'HUB HEIGHT' ) > 0 .OR. INDEX( LINE, 'ZHUB' ) > 0 ) THEN

               READ (LINE, *, IOSTAT = TmpErrStat) RefHt

               IF ( TmpErrStat /= 0 ) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading hub height from FF summary file.', ErrStat, ErrMsg, RoutineName )                  
                  RETURN
               ENDIF
               StrNeeded(2) = .FALSE.

            ENDIF !INDEX for "HUB HEIGHT" or "ZHUB"

   !      ELSEIF ( StrNeeded(3) ) THEN
   !
   !         !-------------------------------------------------------------------------------------------
   !         ! #3: Get the grid width (& height, if available), using the strings "GRID WIDTH" or "RDIAM"
   !         !    If GRID HEIGHT is specified, use it, too. -- THIS IS UNNECESSARY AS IT'S STORED IN THE BINARY FILE
   !         !-------------------------------------------------------------------------------------------

         ELSEIF ( StrNeeded(4) ) THEN

            !-------------------------------------------------------------------------------------------
            ! #4: Get the mean wind speed "UBAR" and turbulence intensities from following lines for
            !     scaling Bladed-style FF binary files
            !-------------------------------------------------------------------------------------------

            IF ( INDEX( LINE, 'UBAR') > 0 ) THEN

               FirstIndx = INDEX( LINE, '=' ) + 1        ! Look for the equal siqn to find the number we're looking for

               READ ( LINE( FirstIndx:LEN(LINE) ), *, IOSTAT=TmpErrStat ) UBar

               IF ( TmpErrStat /= 0 ) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading UBar binary data normalizing parameter from FF summary file.', ErrStat, ErrMsg, RoutineName )                  
                  RETURN
               ENDIF

               DO I = 1,3

                  LineCount = LineCount + 1

                  READ ( UnitWind, '(A)', IOSTAT=TmpErrStat ) LINE
                  IF ( TmpErrStat /= 0 ) THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading line #'//TRIM(Num2LStr(LineCount))//' of the summary file, "'//TRIM(FileName)//&
                                          '". Could not find all of the required parameters.', ErrStat, ErrMsg, RoutineName )                  
                     RETURN
                  ENDIF

                  FirstIndx = INDEX( LINE, '=' ) + 1     ! Read the number between the = and % signs
                  LastIndx  = INDEX( LINE, '%' ) - 1

                  IF ( LastIndx <= FirstIndx ) LastIndx = LEN( LINE )   ! If there's no % sign, read to the end of the line

                  READ ( LINE( FirstIndx:LastIndx ), *, IOSTAT=TmpErrStat ) TI(I)
                  IF ( TmpErrStat /= 0 ) THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading TI('//TRIM(Num2LStr(I))// &
                                 ') binary data normalizing parameter from FF summary file.', ErrStat, ErrMsg, RoutineName )                  
                     RETURN
                  ENDIF

               ENDDO !I

               StrNeeded(4) = .FALSE.

             ENDIF

         ELSEIF ( StrNeeded(5) ) THEN

            !-------------------------------------------------------------------------------------------
            ! #5: Get the grid "HEIGHT OFFSET", if it exists (in TurbSim). Otherwise, assume it's zero
            !           ZGOffset = HH - GridBase - ParamData%FF%FFZHWid
            !-------------------------------------------------------------------------------------------
            IF ( INDEX( LINE, 'HEIGHT OFFSET' ) > 0  ) THEN

               FirstIndx = INDEX ( LINE, '=' ) + 1

               READ ( LINE( FirstIndx:LEN(LINE) ), *, IOSTAT=TmpErrStat ) ZGOffset

               IF ( TmpErrStat /= 0 ) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading height offset from FF summary file.', ErrStat, ErrMsg, RoutineName )                  
                  RETURN
               ENDIF

               StrNeeded(5) = .FALSE.

            ENDIF !INDEX for "HEIGHT OFFSET"

         ELSE
         
            IF ( StrNeeded(6) ) THEN

               !-------------------------------------------------------------------------------------------
               ! #6: Get the grid "PERIODIC", if it exists (in TurbSim). Otherwise, assume it's
               !        not a periodic file (would only show up if the HEIGHT OFFSET is in the file)
               !-------------------------------------------------------------------------------------------
               IF ( INDEX( LINE, 'PERIODIC' ) > 0  ) THEN

                  Periodic      = .TRUE.
                  StrNeeded(6)  = .FALSE.
                  CYCLE
               ENDIF !INDEX for "PERIODIC"
            END IF

            IF ( StrNeeded(7) ) THEN

               IF ( INDEX( LINE, 'BLADED LEFT-HAND RULE') > 0 ) THEN
                  LHR           = .TRUE.
                  StrNeeded(7)  = .FALSE.
               END IF ! INDEX for "BLADED LEFT-HAND RULE"

            END IF

         ENDIF ! StrNeeded

      ENDDO !WHILE

      !-------------------------------------------------------------------------------------------------
      ! Close the summary file
      !-------------------------------------------------------------------------------------------------

      CLOSE ( UnitWind )


      !-------------------------------------------------------------------------------------------------
      ! Calculate the height of the grid center
      !-------------------------------------------------------------------------------------------------

       ZCenter  = RefHt - ZGOffset


   END SUBROUTINE Read_Summary_FF

   !====================================================================================================
   !>   Reads the binary headers from the turbulence files of the old Bladed variety.  Note that
   !!   because of the normalization, neither ParamData%FF%NZGrids or ParamData%FF%NYGrids are larger than 32 points.
   !!   21-Sep-2009 - B. Jonkman, NREL/NWTC.
   !!   16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   SUBROUTINE Read_Bladed_FF_Header0 (UnitWind, p, NativeBladedFmt, ErrStat, ErrMsg)

      IMPLICIT                                              NONE

      CHARACTER(*),        PARAMETER                     :: RoutineName="Read_Bladed_FF_Header0"

         ! Passed Variables:

      INTEGER(IntKi),                     INTENT(IN   )  :: UnitWind  !< unit number of already-opened wind file
      TYPE(IfW_FFWind_ParameterType),     INTENT(INOUT)  :: p         !< Parameters
      LOGICAL,                            INTENT(IN   )  :: NativeBladedFmt !< Whether this should ignore the advection speed in the binary file
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat   !< error status 
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg    !< error message 


         ! Local Variables:
      REAL(ReKi)                                         :: FFXDelt
      REAL(ReKi)                                         :: FFYDelt
      REAL(ReKi)                                         :: FFZDelt

      INTEGER(B2Ki)                                      :: Dum_Int2
      INTEGER(IntKi)                                     :: I


         ! Temporary Error Handling
      INTEGER(IntKi)                                     :: TmpErrStat     ! for checking the IOSTAT from a READ or Open statement


      !-------------------------------------------------------------------------------------------------
      ! Initializations
      !-------------------------------------------------------------------------------------------------

      ErrStat  = ErrID_None
      ErrMsg   = ''


      !-------------------------------------------------------------------------------------------------
      ! Read the header (file has just been opened)
      !-------------------------------------------------------------------------------------------------

         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! -NFFC (file ID)

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading number of wind components from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         p%NFFComp = -1*Dum_Int2


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! delta z (mm)

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading dz from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         FFZDelt = 0.001*Dum_Int2
         p%InvFFZD = 1.0/FFZDelt


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! delta y (mm)

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading dy from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         FFYDelt = 0.001*Dum_Int2
         p%InvFFYD = 1.0/FFYDelt


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! delta x (mm)

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading dx from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         FFXDelt = 0.001*Dum_Int2


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! half the number of time steps

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading number of time steps from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         p%NFFSteps = 2*Dum_Int2


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! 10 times the mean full-field wind speed

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading mean full-field wind speed from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         if (.not. NativeBladedFmt) p%MeanFFWS = 0.1*Dum_Int2
         p%InvMFFWS = 1.0/p%MeanFFWS
         p%FFDTime  = FFXDelt/p%MeanFFWS
         p%FFRate   = 1.0/p%FFDTime


      DO I = 1,5

         ! Read 2-byte integer. Can't use library routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                              ! unused variables: zLu, yLu, xLu, dummy, random seed

            IF (TmpErrStat /= 0) THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading 2-byte integers from binary FF file.', ErrStat, ErrMsg, RoutineName)
               RETURN
            ENDIF

      END DO


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! 1000*nz

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading nz from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         p%NZGrids  = Dum_Int2/1000
         p%FFZHWid  = 0.5*FFZDelt*( p%NZGrids - 1 )    ! half the vertical size of the grid


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! 1000*ny

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading ny from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         p%NYGrids  = Dum_Int2/1000
         p%FFYHWid  = 0.5*FFYDelt*( p%NYGrids - 1 )


      IF (p%NFFComp == 3) THEN

         DO I=1,6

               ! Read 2-byte integer. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                           ! unused variables: zLv, yLv, xLv, zLw, yLw, xLw

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading 2-byte length scales from binary FF file.', ErrStat, ErrMsg, RoutineName)
                  RETURN
               ENDIF

         ENDDO !I

      ENDIF !NFFComp


      RETURN

   END SUBROUTINE Read_Bladed_FF_Header0
   !====================================================================================================
   !>   Reads the binary headers from the turbulence files of the new Bladed variety.
   !!   16-May-2002 - Windward Engineering.
   !!   21-Sep-2009 - B. Jonkman, NREL.  updated to trap errors and add extra parameters for MANN model
   !!   16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   SUBROUTINE Read_Bladed_FF_Header1 (UnitWind, TI, p, NativeBladedFmt, ErrStat, ErrMsg)

      IMPLICIT                                              NONE

      CHARACTER(*),        PARAMETER                     :: RoutineName="Read_Bladed_FF_Header1"


         ! Passed Variables:

      INTEGER(IntKi),                     INTENT(IN   )  :: UnitWind  !< unit number of already-opened wind file
      REAL(ReKi),                         INTENT(  OUT)  :: TI(3)     !< turbulence intensity contained in file header 
      TYPE(IfW_FFWind_ParameterType),     INTENT(INOUT)  :: p         !< Parameters
      LOGICAL,                            INTENT(IN   )  :: NativeBladedFmt !< Whether this should ignore the advection speed in the binary file
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat   !< error status
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg    !< error message


         ! Local Variables:

      REAL(ReKi)                                         :: FFXDelt
      REAL(ReKi)                                         :: FFYDelt
      REAL(ReKi)                                         :: FFZDelt

      REAL(SiKi)                                         :: Dum_Real4
      INTEGER(B2Ki)                                      :: Dum_Int2
      INTEGER(B4Ki)                                      :: Dum_Int4

      INTEGER(IntKi)                                     :: I
      INTEGER(IntKi)                                     :: TurbType


         ! Temporary Error Handling
      INTEGER(IntKi)                                     :: TmpErrStat


      !-------------------------------------------------------------------------------------------------
      ! Initializations
      !-------------------------------------------------------------------------------------------------

      ErrStat  = ErrID_None
      ErrMsg   = ''

      TI(:) = -1                                                                                !Initialize to -1 (not all models contain TI)

      !-------------------------------------------------------------------------------------------------
      ! File reading
      !-------------------------------------------------------------------------------------------------

         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! -99 (file ID)

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading integer from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! turbulence type

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading turbulence type from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         TurbType = Dum_Int2


      SELECT CASE (TurbType)
         CASE(1, 2)
            !----------------------------------------
            !1-component Von Karman (1) or Kaimal (2)
            !----------------------------------------
               p%NFFComp = 1

         CASE(3, 5)
            !----------------------------------------
            !3-component Von Karman (3) or IEC-2
            ! Kaimal (5)
            !----------------------------------------
               p%NFFComp = 3

         CASE(4)
            !----------------------------------------
            !improved Von Karman
            !----------------------------------------

                  ! Read 2-byte integer. Can't use library routines for this.
               READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                        ! number of components (should be 3)

                  IF (TmpErrStat /= 0) THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading number of components from binary FF file.', ErrStat, ErrMsg, RoutineName)
                     RETURN
                  ENDIF
                  p%NFFComp = Dum_Int4

                  ! Read 4-byte real. Can't use library routines for this.
               READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                       ! Latitude (deg)

                  IF (TmpErrStat /= 0) THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading latitude from binary FF file.', ErrStat, ErrMsg, RoutineName)
                     RETURN
                  ENDIF

                  ! Read 4-byte real. Can't use library routines for this.
               READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                       ! Roughness length (m)

                  IF (TmpErrStat /= 0) THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading roughness length from binary FF file.', ErrStat, ErrMsg, RoutineName)
                     RETURN
                  ENDIF

                  ! Read 4-byte real. Can't use library routines for this.
               READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                       ! Reference height (m) = Z(1) + GridHeight / 2.0

                  IF (TmpErrStat /= 0) THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading reference height from binary FF file.', ErrStat, ErrMsg, RoutineName)
                     RETURN
                  ENDIF


               DO I = 1,3
                     ! Read 4-byte real. Can't use library routines for this.
                  READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                    ! TI(u, v, w) (%)

                     IF (TmpErrStat /= 0) THEN
                        CALL SetErrStat( ErrID_Fatal, ' Error reading TI('//'TRIM(Num2LStr(I))'//') from binary FF file.', ErrStat, ErrMsg, RoutineName)
                        RETURN
                     ENDIF
                     TI(I) = Dum_Real4                                                          ! This overwrites the TI read in the summary file

               END DO !I


         CASE (7, 8)
            !----------------------------------------
            ! General Kaimal (7) or  Mann model (8)
            !----------------------------------------

                  ! Read 4-byte integer. Can't use library routines for this.
               READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                        ! number of bytes in header

                  IF (TmpErrStat /= 0) THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading number of header records from binary FF file.', ErrStat, ErrMsg, RoutineName)
                     RETURN
                  ENDIF

                  ! Read 4-byte integer. Can't use library routines for this.
               READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                        ! number of components

                  IF (TmpErrStat /= 0) THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading number of data from binary FF file.', ErrStat, ErrMsg, RoutineName)
                     RETURN
                  ENDIF
                  p%NFFComp = Dum_Int4


         CASE DEFAULT

            CALL SetErrStat( ErrID_Warn, ' InflowWind does not recognize the full-field turbulence file type ='// &
                        TRIM(Num2LStr(TurbType))//'.', ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN

      END SELECT !TurbType


         ! Read 4-byte real. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                                ! delta z (m)

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading dz from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         FFZDelt = Dum_Real4
         p%InvFFZD = 1.0/FFZDelt


         ! Read 4-byte real. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                               ! delta y (m)

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading dy from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         FFYDelt = Dum_Real4
         p%InvFFYD = 1.0/FFYDelt

         ! Read 4-byte real. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                               ! delta x (m)

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading dx from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         FFXDelt = Dum_Real4


         ! Read 4-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                                ! half the number of time steps

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading number of time steps from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         p%NFFSteps = 2*Dum_Int4


         ! Read 4-byte real. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                               ! mean full-field wind speed

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading mean full-field wind speed from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         if (.not. NativeBladedFmt) p%MeanFFWS = Dum_Real4
         p%InvMFFWS = 1.0/p%MeanFFWS
         p%FFDTime  = FFXDelt/p%MeanFFWS
         p%FFRate   = 1.0/p%FFDTime


      DO I = 1,3

            ! Read 4-byte real. Can't use library routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                            ! unused variables: zLu, yLu, xLu

            IF (TmpErrStat /= 0) THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading 4-byte length scales from binary FF file.', ErrStat, ErrMsg, RoutineName)
               RETURN
            ENDIF

      END DO


      DO I = 1,2

         ! Read 4-byte integer. Can't use library routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                             ! unused variables: dummy, random seed

            IF (TmpErrStat /= 0) THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading 4-byte integers from binary FF file.', ErrStat, ErrMsg, RoutineName)
               RETURN
            ENDIF

      END DO


         ! Read 4-integer real. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                                ! nz

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading nz from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         p%NZGrids  = Dum_Int4
         p%FFZHWid  = 0.5*FFZDelt*( p%NZGrids - 1 )    ! half the vertical size of the grid


         ! Read 4-integer real. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                                ! ny

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading ny from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         p%NYGrids  = Dum_Int4
         p%FFYHWid  = 0.5*FFYDelt*( p%NYGrids - 1 )


      IF (p%NFFComp == 3) THEN

         DO I=1,6

               ! Read 4-real real. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                         ! unused variables: zLv, yLv, xLv, zLw, yLw, xLw

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading 4-byte length scales from binary FF file.', ErrStat, ErrMsg, RoutineName)
                  RETURN
               ENDIF

         ENDDO !I

      ENDIF !NFFComp



      IF ( TurbType == 7 ) THEN     ! General Kaimal model

               ! Read 4-real real. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                         ! unused variable: coherence decay constant

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading coherence decay constant from binary FF file.', ErrStat, ErrMsg, RoutineName)
                  RETURN
               ENDIF

               ! Read 4-real real. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                         ! unused variables: coherence scale parameter in m

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading coherence scale parameter from binary FF file.', ErrStat, ErrMsg, RoutineName)
                  RETURN
               ENDIF

      ELSE IF ( TurbType == 8 ) THEN     ! Mann model

         DO I=1,2

               ! Read 4-real real. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                         ! unused variables: shear parameter (gamma), scale length

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading 4-byte parameters from binary FF file.', ErrStat, ErrMsg, RoutineName)
                  RETURN
               ENDIF

         ENDDO !I

         DO I=1,4

               ! Read 4-real real. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                         ! unused variables

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading 4-byte parameters from binary FF file.', ErrStat, ErrMsg, RoutineName)
                  RETURN
               ENDIF

         ENDDO !I

         DO I=1,3

               ! Read 4-integer real. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                          ! unused variables

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading 4-byte parameters from binary FF file.', ErrStat, ErrMsg, RoutineName)
                  RETURN
               ENDIF

         ENDDO !I

         DO I=1,2

               ! Read 4-real real. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                         ! unused variables

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading 4-byte parameters from binary FF file.', ErrStat, ErrMsg, RoutineName)
                  RETURN
               ENDIF

         ENDDO !I

         DO I=1,3

               ! Read 4-integer real. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                          ! unused variables

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading 4-byte parameters from binary FF file.', ErrStat, ErrMsg, RoutineName)
                  RETURN
               ENDIF

         ENDDO !I

         DO I=1,2

               ! Read 4-real real. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                         ! unused variables

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading 4-byte parameters from binary FF file.', ErrStat, ErrMsg, RoutineName)
                  RETURN
               ENDIF

         ENDDO !I


      ENDIF !TurbType


      RETURN

   END SUBROUTINE Read_Bladed_FF_Header1
   !====================================================================================================
   !> This subroutine continues reading UnitWind, starting after the headers have been read.
   !! It reads the Grids and converts the data to un-normalized wind speeds in m/s.
   !!   16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   SUBROUTINE Read_Bladed_Grids ( UnitWind, NativeBladedFmt, CWise, LHR, TI, p, ErrStat, ErrMsg )
      IMPLICIT                                              NONE

      CHARACTER(*),        PARAMETER                     :: RoutineName="Read_Bladed_Grids"

         ! Passed variables

      INTEGER(IntKi),                     INTENT(IN   )  :: UnitWind          !< unit number of already-opened wind file
      LOGICAL,                            INTENT(IN   )  :: NativeBladedFmt   !< whether this data is in native Bladed format (scale to zero mean and unit standard deviation)
      LOGICAL,                            INTENT(IN   )  :: CWise             !< clockwise flag (determines if y is increasing or decreasing in file)
      LOGICAL,                            INTENT(IN   )  :: LHR               !< Left-hand rule for Bladed files (is the v component aligned along *negative* Y?)
      REAL(ReKi),                         INTENT(IN   )  :: TI      (3)       !< turbulence intensities of the wind components as defined in the FF file, not necessarially the actual TI
      TYPE(IfW_FFWind_ParameterType),     INTENT(INOUT)  :: p                 !< Parameters
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat           !< error status
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg            !< error message 

      REAL(ReKi)                                         :: FF_Scale(3)       !< used for "un-normalizing" the data
      REAL(ReKi)                                         :: FF_Offset(3)      !< used for "un-normalizing" the data

      INTEGER(IntKi)                                     :: CFirst
      INTEGER(IntKi)                                     :: CLast
      INTEGER(IntKi)                                     :: CStep
      INTEGER(B2Ki)                                      :: Dum_Int2
      INTEGER(IntKi)                                     :: I
      INTEGER(IntKi)                                     :: IC
      INTEGER(IntKi)                                     :: IR
      INTEGER(IntKi)                                     :: IT

      INTEGER(IntKi)                                     :: TmpNumSteps

         ! Temporary variables for error handling

      INTEGER(IntKi)                                     :: TmpErrStat     ! for checking the result of IOSTAT on READ or Open statements
      CHARACTER(ErrMsgLen)                               :: TmpErrMsg

      IF (NativeBladedFmt) THEN
         FF_Scale  = 0.001_ReKi
         FF_Offset = 0.0_ReKi
      ELSE
         FF_Scale = 0.001_ReKi*p%MeanFFWS*TI/100.0_ReKi
         FF_Offset= (/ p%MeanFFWS, 0.0_ReKi, 0.0_ReKi /)  ! used for "un-normalizing" the data
      END IF

         ! Bladed convention has positive V pointed along negative Y
      IF (LHR) THEN ! left-hand rule
         FF_Scale(2) = -FF_Scale(2)
      END IF
         
            
      !-------------------------------------------------------------------------------------------------
      ! Generate an informative message. Initialize the ErrStat.
      !-------------------------------------------------------------------------------------------------
         ! This could take a while, so we'll write a message to tell users what's going on:
         
      CALL WrScr( NewLine//'   Reading a '//TRIM( Num2LStr(p%NYGrids) )//'x'//TRIM( Num2LStr(p%NZGrids) )//  &
                  ' grid ('//TRIM( Num2LStr(p%FFYHWid*2) )//' m wide, '// &
                  TRIM( Num2LStr(p%GridBase) )//' m to '// &
                  TRIM( Num2LStr(p%GridBase+p%FFZHWid*2) )//&
                  ' m above ground) with a characteristic wind speed of '//TRIM( Num2LStr(p%MeanFFWS) )//' m/s. ' )
      ErrMsg   = ""
      ErrStat  =  ErrID_None


      !-------------------------------------------------------------------------------------------------
      ! Allocate space for the FF array
      !-------------------------------------------------------------------------------------------------

      TmpNumSteps = p%NFFSteps + 1       ! add another step, just in case there is an odd number of steps.

   !bjj: should we reorganize this FFData array so we access the data faster?

      IF ( .NOT. ALLOCATED( p%FFData ) ) THEN
         CALL AllocAry( p%FFData, p%NZGrids,p%NYGrids,p%NFFComp,TmpNumSteps, &
                  'Full-field wind data array.', TmpErrStat, TmpErrMsg )
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN

      ELSE
         IF (SIZE(p%FFData,1) /= p%NZGrids .OR. SIZE(p%FFData,2) /= p%NYGrids .OR. &
             SIZE(p%FFData,3) /= p%NFFComp .OR. SIZE(p%FFData,3) /= TmpNumSteps ) THEN

               ! Let's make the array the correct size (we should never get here, but you never know)

            DEALLOCATE( p%FFData )

            CALL AllocAry( p%FFData, p%NZGrids,p%NYGrids,p%NFFComp,TmpNumSteps, &
                  'Full-field wind data array.', TmpErrStat, TmpErrMsg )
               CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
               IF ( ErrStat >= AbortErrLev ) RETURN            

         ENDIF !Incorrect size
      ENDIF ! allocated

      !-------------------------------------------------------------------------------------------------
      ! Initialize the data and set column indexing to account for direction of turbine rotation (CWise)
      !-------------------------------------------------------------------------------------------------

      p%FFData(:,:,:,:) = 0.0                        ! we may have only one component

      IF ( CWise )  THEN
         CFirst    = p%NYGrids
         CLast     = 1
         CStep     = -1
      ELSE
         CFirst    = 1
         CLast     = p%NYGrids
         CStep     = 1
      ENDIF


      !-------------------------------------------------------------------------------------------------
      ! Loop through all the time steps, reading the data and converting to m/s
      !-------------------------------------------------------------------------------------------------
   !bjj: should we reorganize this FFData array so we access the data faster?

      p%NFFSteps = TmpNumSteps

   TIME_LOOP:  DO IT=1,TmpNumSteps     ! time (add 1 to see if there is an odd number of grids)

         DO IR=1,p%NZGrids               ! the rows (vertical)

            DO IC=CFirst,CLast,CStep   ! the columns (lateral)

               DO I=1,p%NFFComp          ! wind components (U, V, W)

                     ! Get the next integer from the file.
                     ! This is a 2-byte integer, so we can't use the library read routines.
                  READ (UnitWind,IOStat=TmpErrStat)  Dum_Int2
                  IF (TmpErrStat /= 0) THEN
                     IF ( IT == TmpNumSteps ) THEN ! There really were an even number of steps
                        p%NFFSteps = TmpNumSteps - 1
                        ErrStat  = 0
                        EXIT TIME_LOOP
                     ELSE
                        CALL SetErrStat( ErrID_Fatal, ' Error reading binary data file. '// &
                                    'ic = '//TRIM(Num2LStr(ic))// &
                                    ', ir = '//TRIM(Num2LStr(ir))// &
                                    ', it = '//TRIM(Num2LStr(it))// &
                                    ', nffsteps = '//TRIM(Num2LStr(p%NFFSteps)), ErrStat, ErrMsg, RoutineName)
                        RETURN
                     ENDIF
                  ELSE
                     p%FFData(IR,IC,I,IT) = FF_Offset(I)+FF_Scale(I)*Dum_Int2
                  ENDIF

               END DO !I

            END DO !IC

         END DO !IR

      END DO TIME_LOOP !IT

      IF ( p%Periodic ) THEN
         TmpErrMsg = '   Processed '//TRIM( Num2LStr( p%NFFSteps ) )//' time steps of '// &
                    TRIM( Num2LStr ( p%FFRate ) )//'-Hz full-field data (period of '// &
                    TRIM( Num2LStr( p%FFDTime*p%NFFSteps ) )//' seconds).'
         
      ELSE
         TmpErrMsg= '   Processed '//TRIM( Num2LStr( p%NFFSteps ) )//' time steps of '// &
                    TRIM( Num2LStr ( p%FFRate ) )//'-Hz full-field data ('// &
                    TRIM( Num2LStr( p%FFDTime*( p%NFFSteps - 1 ) ) )//' seconds).'
      ENDIF
      CALL WrScr( NewLine//TRIM(TmpErrMsg) )
      !CALL SetErrStat( ErrID_Info, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
      
      

   END SUBROUTINE Read_Bladed_Grids
   !====================================================================================================
   !> This subroutine reads the binary tower file that corresponds with the Bladed-style FF binary file.
   !! The FF grid must be read before this subroutine is called! (many checks are made to ensure the
   !! files belong together)
   !!   16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   SUBROUTINE Read_FF_Tower( UnitWind, p, TwrFileName,  ErrStat, ErrMsg )
      IMPLICIT                                              NONE

      CHARACTER(*),           PARAMETER                  :: RoutineName="Read_FF_Tower"


         ! Passed Variables:
      INTEGER(IntKi)                                     :: UnitWind          !< unit number of wind file to be opened
      TYPE(IfW_FFWind_ParameterType),     INTENT(INOUT)  :: p                 !< Parameters
      CHARACTER(*),                       INTENT(IN   )  :: TwrFileName
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat        !< error status return value (0=no error; non-zero is error)
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg         !< a message for errors that occur

         ! Local Variables:

      REAL(SiKi)                                         :: Dum_Real4      ! dummy 4-byte real number
      INTEGER(B2Ki)                                      :: Dum_Int2       ! dummy 2-byte integer
      INTEGER(B4Ki)                                      :: Dum_Int4       ! dummy 4-byte integer

      INTEGER(IntKi)                                     :: IC             ! loop counter for wind components
      INTEGER(IntKi)                                     :: IT             ! loop counter for time
      INTEGER(IntKi)                                     :: IZ             ! loop counter for z

      REAL(ReKi),    PARAMETER                           :: TOL = 1E-4     ! tolerence for wind file comparisons

      REAL(ReKi),    PARAMETER                           :: FF_Offset(3) = (/ 1.0, 0.0, 0.0 /)  ! used for "un-normalizing" the data
      REAL(SiKi)                                         :: TI       (3)   ! scaling values for "un-normalizing the data" [approx. turbulence intensities of the wind components]


         ! Temporary Error Handling

      INTEGER(IntKi)                                     :: TmpErrStat     ! IOSTAT value.
      CHARACTER(ErrMsgLen)                               :: TmpErrMsg

      !-------------------------------------------------------------------------------------------------
      ! Initialization
      !-------------------------------------------------------------------------------------------------

      ErrMsg   = ''
      ErrStat  = ErrID_None

      p%NTGrids = 0

      IF ( p%NFFComp /= 3 ) THEN
         CALL SetErrStat( ErrID_Fatal, ' Error: Tower binary files require 3 wind components.', ErrStat, ErrMsg, RoutineName )
         RETURN
      ENDIF

      !-------------------------------------------------------------------------------------------------
      ! Open the file
      !-------------------------------------------------------------------------------------------------

      CALL OpenBInpFile (UnitWind, TRIM(TwrFileName), TmpErrStat, TmpErrMsg)
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN

      !-------------------------------------------------------------------------------------------------
      ! Read the header information and check that it's compatible with the FF Bladed-style binary
      ! parameters already read.
      !-------------------------------------------------------------------------------------------------
            ! This is a 4-byte real, so we can't use the library read routines.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4               ! dz, in meters [4-byte REAL]
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading dz in the binary tower file "'//TRIM( TwrFileName )//'."', ErrStat, ErrMsg, RoutineName )
               RETURN
            ENDIF

            IF ( ABS(Dum_Real4*p%InvFFZD-1) > TOL ) THEN
               CALL SetErrStat( ErrID_Fatal, ' Resolution in the FF binary file does not match the tower file.', ErrStat, ErrMsg, RoutineName )
               RETURN
            ENDIF


            ! This is a 4-byte real, so we can't use the library read routines.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4               ! dx, in meters [4-byte REAL]
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading dx in the binary tower file "'//TRIM( TwrFileName )//'."', ErrStat, ErrMsg, RoutineName )
               RETURN
            ENDIF

            IF ( ABS(Dum_Real4*p%InvMFFWS/p%FFDTime-1) > TOL ) THEN
               CALL SetErrStat( ErrID_Fatal, ' Time resolution in the FF binary file does not match the tower file.', ErrStat, ErrMsg, RoutineName )
               RETURN
            ENDIF


            ! This is a 4-byte real, so we can't use the library read routines.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4               ! Zmax, in meters [4-byte REAL]
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading Zmax in the binary tower file "'//TRIM( TwrFileName )//'."', ErrStat, ErrMsg, RoutineName )
               RETURN
            ENDIF

            IF ( ABS(Dum_Real4/p%GridBase-1) > TOL ) THEN
               CALL SetErrStat( ErrID_Fatal, ' Height in the FF binary file does not match the tower file "'//TRIM( TwrFileName )//'."', ErrStat, ErrMsg, RoutineName )
               RETURN
            ENDIF


            ! This is a 4-byte integer, so we can't use the library read routines.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                ! NumOutSteps [4-byte INTEGER]
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading NumOutSteps in the binary tower file "'//TRIM( TwrFileName )//'."', ErrStat, ErrMsg, RoutineName )
               RETURN
            ENDIF

            IF ( Dum_Int4 /= p%NFFSteps ) THEN
               CALL SetErrStat( ErrID_Fatal, ' Number of time steps in the FF binary file does not match the tower file.', ErrStat, ErrMsg, RoutineName )
               RETURN
            ENDIF


            ! This is a 4-byte integer, so we can't use the library read routines.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                ! NumZ      [4-byte INTEGER]
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading NumZ in the binary tower file "'//TRIM( TwrFileName )//'."', ErrStat, ErrMsg, RoutineName )
               RETURN
            ENDIF
            p%NTGrids = Dum_Int4


            ! This is a 4-byte real, so we can't use the library read routines.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4               ! UHub      [4-byte REAL]
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading UHub in the binary tower file "'//TRIM( TwrFileName )//'."', ErrStat, ErrMsg, RoutineName )
               RETURN
            ENDIF

            IF ( ABS(Dum_Real4*p%InvMFFWS - 1) > TOL ) THEN
               CALL SetErrStat( ErrID_Fatal, ' Mean wind speed in the FF binary file does not match the tower file.', ErrStat, ErrMsg, RoutineName )
               p%NTGrids  = 0
               RETURN
            ENDIF


         DO IC=1,3
               ! Read the TI values fromthe tower file: 4-byte reals.
               
               !bjj: not sure you can call this routine to read from a binary file...
            !CALL ReadVar( UnitWind, TRIM(InitInp%WindFileName), TI(IC), 'TI('//TRIM(Num2LStr(IC))//')', 'TI value for u,v, or w', TmpErrStat, TmpErrMsg )
            !IF (TmpErrStat /= ErrID_None) THEN
            !   p%NTGrids  = 0
            !   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
            !   IF (ErrStat >= AbortErrLev) RETURN
            !ENDIF
            !
            READ (UnitWind, IOSTAT=TmpErrStat)   TI(IC)               ! TI(u), TI(v), TI(w)  [4-byte REAL]
            
            IF (TmpErrStat /= 0) THEN
               p%NTGrids  = 0
               CALL SetErrStat( ErrID_Fatal, ' Error reading TI('//TRIM(Num2LStr(IC))//') in the binary tower file "' &
                               //TRIM( TwrFileName )//'."', ErrStat, ErrMsg, RoutineName )
               RETURN
            ENDIF
         
         END DO

      !----------------------------------------------------------------------------------------------
      ! Allocate arrays for the tower points
      !----------------------------------------------------------------------------------------------

         IF ( p%NTGrids > 0 ) THEN

            IF ( .NOT. ALLOCATED( p%FFTower ) ) THEN
               CALL AllocAry( p%FFTower, p%NFFComp, p%NTGrids, p%NFFSteps, &
                  'Tower wind data array.', TmpErrStat, TmpErrMsg )
               CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
               IF (ErrStat >= AbortErrLev) RETURN
               
            ELSE
               ! Check sizes here!
            ENDIF

         ENDIF

      !-------------------------------------------------------------------------------------------------
      ! Read the 16-bit time-series data and scale it to 32-bit reals
      !-------------------------------------------------------------------------------------------------

         ! Loop through time.

         DO IT=1,p%NFFSteps

            DO IZ=1,p%NTGrids         ! If NTGrids<1, there are no tower points & FFTower is not allocated

               ! Ytower     = 0               ! Lateral location of the tower data point, in m relative to tower centerline
               ! Ztower(IZ) = Z1 - (IZ-1)*dz  ! Vertical location of tower data point, in m relative to ground

               DO IC=1,p%NFFComp   ! number of wind components

                     ! Read in the 2-byte integer. Can't use library read routines for this.
                  READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2       ! normalized wind-component, INT(2)
                  IF ( TmpErrStat /= 0 )  THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading binary tower data file. it = '//TRIM(Num2LStr(it))// &
                                    ', nffsteps = '//TRIM(Num2LStr(p%NFFSteps)), ErrStat, ErrMsg, RoutineName )
                     p%NTGrids  = 0
                     RETURN
                  ENDIF

                  p%FFTower(IC,IZ,IT) = p%MeanFFWS*(FF_Offset(IC)+0.00001*TI(IC)*Dum_Int2)   ! wind-component scaled to m/s

               ENDDO !IC

            ENDDO ! IZ


         ENDDO ! IT

      !-------------------------------------------------------------------------------------------------
      ! Close the file
      !-------------------------------------------------------------------------------------------------
      CLOSE ( UnitWind )

      TmpErrMsg = '   Processed '//TRIM( Num2LStr(p%NFFSteps) )//' time steps of '// &
            TRIM( Num2LStr(p%NTGrids) )//'x1 tower data grids.'
      
      !CALL SetErrStat( ErrID_Info, ErrMsgLcl, ErrStat, ErrMsg, RoutineName )
      CALL WrScr( NewLine//TRIM(TmpErrMsg) )
      
      RETURN

   END SUBROUTINE Read_FF_Tower
!====================================================================================================
!> This subroutine reads the text summary file to get normalizing parameters, the location of the
!! grid, and the direction the grid was written to the binary file
SUBROUTINE Read_NativeBladedSummary ( FileName, PLExp, TI, UBar, RefHt, PropagationDir, VFlowAngle, BinFileName, XOffset, ErrStat, ErrMsg )

   IMPLICIT                                              NONE

   CHARACTER(*),        PARAMETER                     :: RoutineName="Read_NativeBladedSummary"


      ! Passed variables
   CHARACTER(*),                       INTENT(IN   )  :: FileName       !< name of the summary file
   REAL(ReKi),                         INTENT(  OUT)  :: PLExp          !< the power-law exponent for vertical wind shear
   REAL(ReKi),                         INTENT(  OUT)  :: TI      (3)    !< turbulence intensities of the wind components as defined in the FF file, not necessarially the actual TI
   REAL(ReKi),                         INTENT(  OUT)  :: UBar           !< mean (advection) wind speed
   REAL(ReKi),                         INTENT(  OUT)  :: RefHt          !< Reference height
   REAL(ReKi),                         INTENT(  OUT)  :: PropagationDir !< propagation direction
   REAL(ReKi),                         INTENT(  OUT)  :: VFlowAngle     !< vertical flow angle
   CHARACTER(*),                       INTENT(  OUT)  :: BinFileName    !< name of the binary file containing wind data
   REAL(ReKi),                         INTENT(  OUT)  :: XOffset        !< distance offset for start of wind files
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat        !< returns 0 if no error encountered in the subroutine
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg         !< holds the error messages

      ! Local variables
   INTEGER(IntKi), PARAMETER                          :: UnEc= -1       ! echo file unit number (set to something else > 0 for debugging)
   INTEGER(IntKi)                                     :: CurLine        ! Current line to parse in FileInfo data structure
   INTEGER(IntKi)                                     :: ErrStat2       ! temporary error status
   CHARACTER(ErrMsgLen)                               :: ErrMsg2        ! temporary error message

   TYPE (FileInfoType)                                :: FileInfo       ! The derived type for holding the file information.
   
   
      !----------------------------------------------------------------------------------------------
      ! Initialize some variables
      !----------------------------------------------------------------------------------------------

   ErrStat      = ErrID_None
   ErrMsg       = ''
   
      !----------------------------------------------------------------------------------------------
      ! Open and read the summary file; store data in FileInfo structure.
      !----------------------------------------------------------------------------------------------

   CALL ProcessComFile ( FileName, FileInfo, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
      
      !-------------------------------------------------------------------------------------------------
      ! Process the lines stored in FileInfo
      !-------------------------------------------------------------------------------------------------
      
   CurLine = 1   

   CALL ParseVar ( FileInfo, CurLine, 'UBAR', UBar, ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL ParseVar ( FileInfo, CurLine, 'REFHT', RefHt, ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   CALL ParseVar ( FileInfo, CurLine, 'TI', TI(1), ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL ParseVar ( FileInfo, CurLine, 'TI_V', TI(2), ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   CALL ParseVar ( FileInfo, CurLine, 'TI_W', TI(3), ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL ParseVar ( FileInfo, CurLine, 'WDIR', PropagationDir, ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      PropagationDir = R2D*PropagationDir

   CALL ParseVar ( FileInfo, CurLine, 'FLINC', VFlowAngle, ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      VFlowAngle = R2D*VFlowAngle ! convert to degrees 

   CALL ParseVar ( FileInfo, CurLine, 'WINDF', BinFileName, ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   CALL ParseVar ( FileInfo, CurLine, 'WSHEAR', PLExp, ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF

   CALL ParseVar ( FileInfo, CurLine, 'XOffset', XOffset, ErrStat2, ErrMsg2, UnEc )
   if (ErrStat2/=ErrID_None) then
      XOffset  = 0.0_ReKi ! this will be the default if offset is not in the file
   end if
            
      
      !-------------------------------------------------------------------------------------------------
      ! Get rid of the FileInfo data structure (including pointers and allocatable array):
      !-------------------------------------------------------------------------------------------------

   call Cleanup ( )


CONTAINS
                  
   SUBROUTINE Cleanup ()
      CALL NWTC_Library_DestroyFileInfoType (FileInfo, ErrStat2, ErrMsg2)            
   END SUBROUTINE Cleanup

END SUBROUTINE Read_NativeBladedSummary
!====================================================================================================


!====================================================================================================
!> This routine acts as a wrapper for the GetWindSpeed routine. It steps through the array of input
!! positions and calls the GetWindSpeed routine to calculate the velocities at each point.
!!
!! There are inefficiencies in how this set of routines is coded, but that is a problem for another
!! day. For now, it merely needs to be functional. It can be fixed up and made all pretty later.
!!
!!   16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
SUBROUTINE IfW_BladedFFWind_CalcOutput(Time, PositionXYZ, ParamData, Velocity, DiskVel, MiscVars, ErrStat, ErrMsg)

   IMPLICIT                                                       NONE

   CHARACTER(*),        PARAMETER                              :: RoutineName="IfW_BladedFFWind_CalcOutput"

      ! Passed Variables
   REAL(DbKi),                                  INTENT(IN   )  :: Time              !< time from the start of the simulation
   REAL(ReKi),                                  INTENT(IN   )  :: PositionXYZ(:,:)  !< Array of XYZ coordinates, 3xN
   TYPE(IfW_BladedFFWind_ParameterType),        INTENT(IN   )  :: ParamData         !< Parameters
   REAL(ReKi),                                  INTENT(INOUT)  :: Velocity(:,:)     !< Velocity output at Time    (Set to INOUT so that array does not get deallocated)
   REAL(ReKi),                                  INTENT(  OUT)  :: DiskVel(3)        !< HACK for AD14: disk velocity output at Time
   TYPE(IfW_BladedFFWind_MiscVarType),          INTENT(INOUT)  :: MiscVars          !< misc/optimization data (storage for the main data)

      ! Error handling
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           !< error status
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            !< The error message



   CALL IfW_FFWind_CalcOutput(Time, PositionXYZ, ParamData%FF, Velocity, DiskVel, ErrStat, ErrMsg)


   RETURN

END SUBROUTINE IfW_BladedFFWind_CalcOutput

!====================================================================================================
!!  This subroutine cleans up any data that is still allocated.  The (possibly) open files are
!!  closed in InflowWindMod.
!!
!!  16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
SUBROUTINE IfW_BladedFFWind_End( ParamData, MiscVars, ErrStat, ErrMsg)

   IMPLICIT                                                       NONE

   CHARACTER(*),                 PARAMETER                     :: RoutineName="IfW_BladedFFWind_End"
      ! Passed Variables
   TYPE(IfW_BladedFFWind_ParameterType),        INTENT(INOUT)  :: ParamData         !< Parameters
   TYPE(IfW_BladedFFWind_MiscVarType),          INTENT(INOUT)  :: MiscVars          !< misc/optimization data (storage for the main data)


      ! Error Handling
   INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat        !< determines if an error has been encountered
   CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg         !< Message about errors


      ! Local Variables
   INTEGER(IntKi)                                        :: TmpErrStat     ! temporary error status
   CHARACTER(ErrMsgLen)                                  :: TmpErrMsg      ! temporary error message


      !-=- Initialize the routine -=-

   ErrMsg   = ''
   ErrStat  = ErrID_None



      ! Destroy parameter data

   CALL IfW_BladedFFWind_DestroyParam(       ParamData,     TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )


      ! Destroy the state data

   CALL IfW_BladedFFWind_DestroyMisc(  MiscVars,   TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )



END SUBROUTINE IfW_BladedFFWind_End

!====================================================================================================
END MODULE IfW_BladedFFWind
