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
SUBROUTINE IfW_BladedFFWind_Init(InitData, ParamData, MiscVars, Interval, InitOutData, ErrStat, ErrMsg)
   IMPLICIT                                                       NONE

   CHARACTER(*),              PARAMETER                        :: RoutineName="IfW_BladedFFWind_Init"

      ! Passed Variables
   TYPE(IfW_BladedFFWind_InitInputType),        INTENT(IN   )  :: InitData       !< Initialization data passed to the module
   TYPE(IfW_BladedFFWind_ParameterType),        INTENT(  OUT)  :: ParamData      !< Parameters
   TYPE(IfW_BladedFFWind_MiscVarType),          INTENT(  OUT)  :: MiscVars       !< misc/optimization data   (storage for the main data)
   TYPE(IfW_BladedFFWind_InitOutputType),       INTENT(  OUT)  :: InitOutData    !< Initial output

   REAL(DbKi),                                  INTENT(IN   )  :: Interval       !< Time Interval to use (passed through here)


      ! Error Handling
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat        !< determines if an error has been encountered
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg         !< Message about errors


      ! Temporary variables for error handling
   INTEGER(IntKi)                                        :: TmpErrStat     ! temporary error status
   CHARACTER(ErrMsgLen)                                  :: TmpErrMsg      ! temporary error message


      ! Local Variables:

   REAL(ReKi)                                            :: TI      (3)    ! turbulence intensities of the wind components as defined in the FF file, not necessarially the actual TI
   REAL(ReKi)                                            :: BinTI   (3)    ! turbulence intensities of the wind components as defined in the FF binary file, not necessarially the actual TI
   REAL(ReKi)                                            :: UBar
   REAL(ReKi)                                            :: ZCenter
 
   INTEGER(IntKi)                                        :: UnitWind     ! Unit number for the InflowWind input file
   INTEGER(B2Ki)                                         :: Dum_Int2
   INTEGER(IntKi)                                        :: I
   LOGICAL                                               :: CWise
   LOGICAL                                               :: Exists
   CHARACTER( 1028 )                                     :: SumFile        ! length is LEN(ParamData%WindFileName) + the 4-character extension.
   CHARACTER( 1028 )                                     :: TwrFile        ! length is LEN(ParamData%WindFileName) + the 4-character extension.



      !-------------------------------------------------------------------------------------------------
      ! Initialize temporary variables
      !-------------------------------------------------------------------------------------------------

   ErrMsg      = ''
   ErrStat     = ErrID_None

   TmpErrMsg   = ''
   TmpErrStat  = ErrID_None


      ! Get a unit number to use

   CALL GetNewUnit(UnitWind, TmpErrStat, TmpErrMsg)
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'IfW_BladedFFWind')
   IF ( ErrStat >= AbortErrLev ) RETURN


      !-------------------------------------------------------------------------------------------------
      ! Copy things from the InitData to the ParamData
      !-------------------------------------------------------------------------------------------------


      !----------------------------------------------------------------------------------------------
      ! Open the binary file, read its "header" (first 2-byte integer) to determine what format
      ! binary file it is, and close it.
      !----------------------------------------------------------------------------------------------

   CALL OpenBInpFile (UnitWind, TRIM(InitData%WindFileName), TmpErrStat, TmpErrMsg)
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'IfW_BladedFFWind')
   IF ( ErrStat >= AbortErrLev ) RETURN

      ! Read the first binary integer from the file to get info on the type.
      ! Cannot use library read routines since this is a 2-byte integer.
   READ ( UnitWind, IOSTAT=TmpErrStat )  Dum_Int2
   CLOSE( UnitWind )

   IF (TmpErrStat /= 0) THEN
      CALL SetErrStat( ErrID_Fatal, ' Error reading first binary integer from file "'//TRIM(InitData%WindFileName)//'."', &
            ErrStat, ErrMsg, 'IfW_BladedFFWind')
      RETURN
   ENDIF


      !----------------------------------------------------------------------------------------------
      ! Read the files to get the required FF data.
      !----------------------------------------------------------------------------------------------

      ! Store the binary format information so the InflowWind code can use it.
      ! Also changes to IntKi from INT(2) to compare in the SELECT below
   ParamData%WindFileFormat = Dum_Int2


   SELECT CASE (ParamData%WindFileFormat)

      CASE ( -1, -2, -3, -99 )                                         ! Bladed-style binary format

         !...........................................................................................
         ! Create full-field summary file name from binary file root name.  Also get tower file
         ! name.
         !...........................................................................................

            CALL GetRoot(InitData%WindFileName, SumFile)

            TwrFile = TRIM(SumFile)//'.twr'
            SumFile = TRIM(SumFile)//'.sum'


         !...........................................................................................
         ! Read the summary file to get necessary scaling information
         !...........................................................................................

            CALL Read_Summary_FF (UnitWind, TRIM(SumFile), CWise, ZCenter, TI, TmpErrStat, TmpErrMsg )
               CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )                  
               IF ( ErrStat >= AbortErrLev ) THEN
                  CLOSE ( UnitWind )
                  RETURN
               END IF
            
            UBar = ParamData%MeanFFWS      ! temporary storage .... this is our only check to see if the summary and binary files "match"


         !...........................................................................................
         ! Open the binary file and read its header
         !...........................................................................................

            CALL OpenBInpFile (UnitWind, TRIM(InitData%WindFileName), TmpErrStat, TmpErrMsg )
               CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )                  
               IF ( ErrStat >= AbortErrLev ) THEN
                  CLOSE ( UnitWind )
                  RETURN
               END IF
            
            IF ( Dum_Int2 == -99 ) THEN                                                      ! Newer-style BLADED format
               CALL Read_Bladed_FF_Header1 (UnitWind, BinTI, TmpErrStat, TmpErrMsg)
                  CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )                  
                  IF ( ErrStat >= AbortErrLev ) THEN
                     CLOSE ( UnitWind )
                     RETURN
                  END IF

                  ! If the TIs are also in the binary file (BinTI > 0),
                  ! use those numbers instead of ones from the summary file

               DO I =1,ParamData%NFFComp
                  IF ( BinTI(I) > 0 ) TI(I) = BinTI(I)
               ENDDO

            ELSE
               CALL Read_Bladed_FF_Header0 (UnitWind, TmpErrStat, TmpErrMsg)     ! Older-style BLADED format
                  CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )                  
                  IF ( ErrStat >= AbortErrLev ) THEN
                     CLOSE ( UnitWind )
                     RETURN
                  END IF

            ENDIF



         !...........................................................................................
         ! Let's see if the summary and binary FF wind files go together before continuing.
         !...........................................................................................

            IF ( ABS( UBar - ParamData%MeanFFWS ) > 0.1 )  THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error: Incompatible mean hub-height wind speeds in FF wind files. '//&
                           '(Check that the .sum and .wnd files were generated together.)', ErrStat, ErrMsg, RoutineName )                  
               RETURN
            ENDIF


         !...........................................................................................
         ! Calculate the height of the bottom of the grid
         !...........................................................................................

            ParamData%GridBase = ZCenter - ParamData%FFZHWid         ! the location, in meters, of the bottom of the grid


         !...........................................................................................
         ! Read the binary grids (converted to m/s) and close the file
         !...........................................................................................

            CALL Read_Bladed_Grids( MiscVars, CWise, TI, TmpErrStat, TmpErrMsg)
            CLOSE ( UnitWind )

               CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )                  
               IF ( ErrStat >= AbortErrLev ) RETURN
            
         !...........................................................................................
         ! Read the tower points file
         !...........................................................................................

            IF ( InitData%TowerFileExist ) THEN      ! If we specified a tower file
               INQUIRE ( FILE=TRIM(TwrFile) , EXIST=Exists )

                  ! Double check that the tower file exists and read it.  If it was requested but doesn't exist,
                  ! throw fatal error and exit.
               IF (  Exists )  THEN
                  CALL Read_FF_Tower( MiscVars, TmpErrStat, TmpErrMsg )
                  CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )                  
                  IF ( ErrStat >= AbortErrLev ) THEN
                     CLOSE ( UnitWind )
                     RETURN
                  END IF
                  ParamData%TowerDataExist   =  .TRUE.
               ELSE
                  CALL SetErrStat( ErrID_Fatal, ' Tower file '//TRIM(TwrFile)//' specified for Bladed full-field '// &
                           'wind files does not exist.', ErrStat, ErrMsg, RoutineName)
                  CLOSE ( UnitWind )
                  RETURN
               ENDIF
            ELSE
               ParamData%NTGrids          =  0_IntKi
               ParamData%TowerDataExist   =  .FALSE.
            ENDIF


      CASE DEFAULT
         CALL SetErrStat( ErrID_Fatal, ' This is not a bladed-style binary wind file (binary format identifier: '//  &
                  TRIM(Num2LStr(ParamData%WindFileFormat))//'.  This might be a TurbSim binary wind file.', &
                  ErrStat, ErrMsg, RoutineName )
         RETURN

   END SELECT


   IF (ParamData%Periodic) THEN
      ParamData%InitXPosition = 0                ! start at the hub
      ParamData%TotalTime     = ParamData%NFFSteps*ParamData%FFDTime
   ELSE
      ParamData%InitXPosition = ParamData%FFYHWid          ! start half the grid width ahead of the turbine
      ParamData%TotalTime     = (ParamData%NFFSteps-1)*ParamData%FFDTime
   ENDIF



      !-------------------------------------------------------------------------------------------------
      ! Set the InitOutput information
      !-------------------------------------------------------------------------------------------------

   InitOutdata%Ver         =  IfW_BladedFFWind_Ver
   InitOutdata%TI          =  TI



      !-------------------------------------------------------------------------------------------------
      ! Write to the summary file
      !-------------------------------------------------------------------------------------------------

   IF ( InitData%SumFileUnit > 0 ) THEN
      WRITE(InitData%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)
      WRITE(InitData%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)    'Bladed-style wind type.  Read by InflowWind sub-module '//   &
                                                                     TRIM(IfW_BladedFFWind_Ver%Name)//' '//TRIM(IfW_BladedFFWind_Ver%Ver)
      WRITE(InitData%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)    TRIM(TmpErrMsg)
      WRITE(InitData%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)    '     FileName:                    '//TRIM(InitData%WindFileName)
      WRITE(InitData%SumFileUnit,'(A34,I3)',   IOSTAT=TmpErrStat)    '     Binary file format id:       ',ParamData%WindFileFormat
      WRITE(InitData%SumFileUnit,'(A34,G12.4)',IOSTAT=TmpErrStat)    '     Reference height (m):        ',ParamData%RefHt
      WRITE(InitData%SumFileUnit,'(A34,G12.4)',IOSTAT=TmpErrStat)    '     Timestep (s):                ',ParamData%FFDTime
      WRITE(InitData%SumFileUnit,'(A34,I12)',  IOSTAT=TmpErrStat)    '     Number of timesteps:         ',ParamData%NFFSteps
      WRITE(InitData%SumFileUnit,'(A34,G12.4)',IOSTAT=TmpErrStat)    '     Mean windspeed (m/s):        ',ParamData%MeanFFWS
      WRITE(InitData%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)    '     Characteristic TI:            [ '// &
                                                TRIM(Num2LStr(TI(1)))//', '//TRIM(Num2LStr(TI(2)))//', '//TRIM(Num2LStr(TI(3)))//' ] '
      WRITE(InitData%SumFileUnit,'(A34,L1)',   IOSTAT=TmpErrStat)    '     Windfile is periodic:        ',ParamData%Periodic
      WRITE(InitData%SumFileUnit,'(A34,L1)',   IOSTAT=TmpErrStat)    '     Windfile includes tower:     ',ParamData%TowerDataExist

      IF ( ParamData%Periodic ) THEN
         WRITE(InitData%SumFileUnit,'(A)',     IOSTAT=TmpErrStat)    '     Time range (s):              [ '// &
                     TRIM(Num2LStr(0.0_ReKi))//' : '//TRIM(Num2LStr(ParamData%TotalTime))//' ]'
      ELSE  ! Shift the time range to compensate for the shifting of the wind grid
         WRITE(InitData%SumFileUnit,'(A)',     IOSTAT=TmpErrStat)    '     Time range (s):              [ '// &
                     TRIM(Num2LStr(-ParamData%InitXPosition*ParamData%InvMFFWS))//' : '// &
                     TRIM(Num2LStr(ParamData%TotalTime-ParamData%InitXPosition*ParamData%InvMFFWS))//' ]'
      ENDIF

      WRITE(InitData%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)    '     Y range (m):                 [ '// &
                     TRIM(Num2LStr(-ParamData%FFYHWid))//' : '//TRIM(Num2LStr(ParamData%FFYHWid))//' ]'

      IF ( ParamData%TowerDataExist ) THEN
         WRITE(InitData%SumFileUnit,'(A)',     IOSTAT=TmpErrStat)    '     Z range (m):                 [ '// &
                     TRIM(Num2LStr(0.0_ReKi))//' : '//TRIM(Num2LStr(ParamData%RefHt + ParamData%FFZHWid))//' ]'
      ELSE
         WRITE(InitData%SumFileUnit,'(A)',     IOSTAT=TmpErrStat)    '     Z range (m):                 [ '// &
                     TRIM(Num2LStr(ParamData%RefHt - ParamData%FFZHWid))//' : '//TRIM(Num2LStr(ParamData%RefHt + ParamData%FFZHWid))//' ]'
      ENDIF


         ! We are assuming that if the last line was written ok, then all of them were.
      IF (TmpErrStat /= 0_IntKi) THEN
         CALL SetErrStat(ErrID_Fatal,'Error writing to summary file.',ErrStat,ErrMsg,RoutineName)
         RETURN
      ENDIF   
   ENDIF 



 
   RETURN


   CONTAINS

   !====================================================================================================
   !> This subroutine reads the text summary file to get normalizing parameters, the location of the
   !! grid, and the direction the grid was written to the binary file
   !!
   !!   16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   SUBROUTINE Read_Summary_FF ( UnitWind, FileName, CWise, ZCenter, TI, ErrStat, ErrMsg )

      IMPLICIT                                              NONE

      CHARACTER(*),        PARAMETER                     :: RoutineName="Read_Summary_FF"


         ! Passed variables
      INTEGER(IntKi),                     INTENT(IN   )  :: UnitWind       !< unit number for the file to open
      CHARACTER(*),                       INTENT(IN   )  :: FileName       !< name of the summary file
      LOGICAL,                            INTENT(  OUT)  :: CWise          !< rotation (for reading the order of the binary data)
      REAL(ReKi),                         INTENT(  OUT)  :: ZCenter        !< the height at the center of the grid
      REAL(ReKi),                         INTENT(  OUT)  :: TI      (3)    !< turbulence intensities of the wind components as defined in the FF file, not necessarially the actual TI
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat        !< returns 0 if no error encountered in the subroutine
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg         !< holds the error messages

        ! Local variables
      REAL(ReKi)                                         :: ZGOffset       ! The vertical offset of the turbine on rectangular grid (allows turbulence not centered on turbine hub)

      INTEGER, PARAMETER                                 :: NumStrings = 6 ! number of strings to be looking for in the file

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
      ParamData%RefHt      = 0.0
      ParamData%Periodic   = .FALSE.

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
         
         
      DO WHILE ( ( ErrStat == ErrID_None ) .AND. StrNeeded(NumStrings) )

         LineCount = LineCount + 1

         READ ( UnitWind, '(A)', IOSTAT=TmpErrStat ) LINE
         IF ( TmpErrStat /= 0 ) THEN

               ! the "HEIGHT OFFSET" and "PERIODIC" parameters are not necessary.  We'll assume they are zero/false if we didn't find it.
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

               READ (LINE, *, IOSTAT = TmpErrStat) ParamData%RefHt

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

               READ ( LINE( FirstIndx:LEN(LINE) ), *, IOSTAT=TmpErrStat ) ParamData%MeanFFWS

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
            !           ZGOffset = HH - GridBase - ParamData%FFZHWid
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

         ELSEIF ( StrNeeded(6) ) THEN

            !-------------------------------------------------------------------------------------------
            ! #5: Get the grid "PERIODIC", if it exists (in TurbSim). Otherwise, assume it's
            !        not a periodic file
            !-------------------------------------------------------------------------------------------
            IF ( INDEX( LINE, 'PERIODIC' ) > 0  ) THEN

               ParamData%Periodic   = .TRUE.
               StrNeeded(6)         = .FALSE.

            ENDIF !INDEX for "PERIODIC"

         ENDIF ! StrNeeded


      ENDDO !WHILE


      !-------------------------------------------------------------------------------------------------
      ! Close the summary file
      !-------------------------------------------------------------------------------------------------

      CLOSE ( UnitWind )


      !-------------------------------------------------------------------------------------------------
      ! Calculate the height of the grid center
      !-------------------------------------------------------------------------------------------------

       ZCenter  = ParamData%RefHt - ZGOffset


   END SUBROUTINE Read_Summary_FF

   !====================================================================================================
   !>   Reads the binary headers from the turbulence files of the old Bladed variety.  Note that
   !!   because of the normalization, neither ParamData%NZGrids or ParamData%NYGrids are larger than 32 points.
   !!   21-Sep-2009 - B. Jonkman, NREL/NWTC.
   !!   16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   SUBROUTINE Read_Bladed_FF_Header0 (UnitWind, ErrStat, ErrMsg)

      IMPLICIT                                              NONE

      CHARACTER(*),        PARAMETER                     :: RoutineName="Read_Bladed_FF_Header0"

         ! Passed Variables:

      INTEGER(IntKi),                     INTENT(IN   )  :: UnitWind  !< unit number of already-opened wind file
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
      CHARACTER(ErrMsgLen)                               :: TmpErrMsg      ! Temporary ErrMsg


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
         ParamData%NFFComp = -1*Dum_Int2


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! delta z (mm)

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading dz from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         FFZDelt = 0.001*Dum_Int2
         ParamData%InvFFZD = 1.0/FFZDelt


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! delta y (mm)

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading dy from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         FFYDelt = 0.001*Dum_Int2
         ParamData%InvFFYD = 1.0/FFYDelt


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
         ParamData%NFFSteps = 2*Dum_Int2


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! 10 times the mean full-field wind speed

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading mean full-field wind speed from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         ParamData%MeanFFWS = 0.1*Dum_Int2
         ParamData%InvMFFWS = 1.0/ParamData%MeanFFWS
         ParamData%FFDTime  = FFXDelt/ParamData%MeanFFWS
         ParamData%FFRate   = 1.0/ParamData%FFDTime


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
         ParamData%NZGrids  = Dum_Int2/1000
         ParamData%FFZHWid  = 0.5*FFZDelt*( ParamData%NZGrids - 1 )    ! half the vertical size of the grid


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! 1000*ny

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading ny from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         ParamData%NYGrids  = Dum_Int2/1000
         ParamData%FFYHWid  = 0.5*FFYDelt*( ParamData%NYGrids - 1 )


      IF (ParamData%NFFComp == 3) THEN

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
   SUBROUTINE Read_Bladed_FF_Header1 (UnitWind, TI, ErrStat, ErrMsg)

      IMPLICIT                                              NONE

      CHARACTER(*),        PARAMETER                     :: RoutineName="Read_Bladed_FF_Header1"


         ! Passed Variables:

      INTEGER(IntKi),                     INTENT(IN   )  :: UnitWind  !< unit number of already-opened wind file
      REAL(ReKi),                         INTENT(  OUT)  :: TI(3)     !< turbulence intensity contained in file header 
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
      CHARACTER(ErrMsgLen)                               :: TmpErrMsg


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
               ParamData%NFFComp = 1

         CASE(3, 5)
            !----------------------------------------
            !3-component Von Karman (3) or IEC-2
            ! Kaimal (5)
            !----------------------------------------
               ParamData%NFFComp = 3

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
                  ParamData%NFFComp = Dum_Int4

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
                  ParamData%NFFComp = Dum_Int4


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
         ParamData%InvFFZD = 1.0/FFZDelt


         ! Read 4-byte real. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                               ! delta y (m)

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading dy from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         FFYDelt = Dum_Real4
         ParamData%InvFFYD = 1.0/FFYDelt

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
         ParamData%NFFSteps = 2*Dum_Int4


         ! Read 4-byte real. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                               ! mean full-field wind speed

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading mean full-field wind speed from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         ParamData%MeanFFWS = Dum_Real4
         ParamData%InvMFFWS = 1.0/ParamData%MeanFFWS
         ParamData%FFDTime  = FFXDelt/ParamData%MeanFFWS
         ParamData%FFRate   = 1.0/ParamData%FFDTime


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
         ParamData%NZGrids  = Dum_Int4
         ParamData%FFZHWid  = 0.5*FFZDelt*( ParamData%NZGrids - 1 )    ! half the vertical size of the grid


         ! Read 4-integer real. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                                ! ny

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading ny from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         ParamData%NYGrids  = Dum_Int4
         ParamData%FFYHWid  = 0.5*FFYDelt*( ParamData%NYGrids - 1 )


      IF (ParamData%NFFComp == 3) THEN

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
   SUBROUTINE Read_Bladed_Grids ( MiscVars, CWise, TI, ErrStat, ErrMsg )
      IMPLICIT                                              NONE

      CHARACTER(*),        PARAMETER                     :: RoutineName="Read_Bladed_Grids"

         ! Passed variables

   TYPE(IfW_BladedFFWind_MiscVarType),    INTENT(INOUT)  :: MiscVars       !< misc/optimization data (storage for the main data)
      LOGICAL,                            INTENT(IN   )  :: CWise          !< clockwise flag (determines if y is increasing or decreasing in file)
      REAL(ReKi),                         INTENT(IN   )  :: TI      (3)    !< turbulence intensities of the wind components as defined in the FF file, not necessarially the actual TI
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat        !< error status
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg         !< error message 

      REAL(ReKi),    PARAMETER                           :: FF_Offset(3) = (/ 1.0, 0.0, 0.0 /)  ! used for "un-normalizing" the data

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


      !-------------------------------------------------------------------------------------------------
      ! Generate an informative message. Initialize the ErrStat.
      !-------------------------------------------------------------------------------------------------
         ! This could take a while, so we'll write a message to tell users what's going on:
         
      CALL WrScr( NewLine//'   Reading a '//TRIM( Num2LStr(ParamData%NYGrids) )//'x'//TRIM( Num2LStr(ParamData%NZGrids) )//  &
                  ' grid ('//TRIM( Num2LStr(ParamData%FFYHWid*2) )//' m wide, '// &
                  TRIM( Num2LStr(ParamData%GridBase) )//' m to '// &
                  TRIM( Num2LStr(ParamData%GridBase+ParamData%FFZHWid*2) )//&
                  ' m above ground) with a characteristic wind speed of '//TRIM( Num2LStr(ParamData%MeanFFWS) )//' m/s. ' )
      ErrMsg   = ""
      ErrStat  =  ErrID_None


      !-------------------------------------------------------------------------------------------------
      ! Allocate space for the FF array
      !-------------------------------------------------------------------------------------------------

      TmpNumSteps = ParamData%NFFSteps + 1       ! add another step, just in case there is an odd number of steps.

   !bjj: should we reorganize this FFData array so we access the data faster?

      IF ( .NOT. ALLOCATED( ParamData%FFData ) ) THEN
         CALL AllocAry( ParamData%FFData, ParamData%NZGrids,ParamData%NYGrids,ParamData%NFFComp,TmpNumSteps, &
                  'Full-field wind data array.', TmpErrStat, TmpErrMsg )
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN

      ELSE
         IF (SIZE(ParamData%FFData,1) /= ParamData%NZGrids .OR. SIZE(ParamData%FFData,2) /= ParamData%NYGrids .OR. &
             SIZE(ParamData%FFData,3) /= ParamData%NFFComp .OR. SIZE(ParamData%FFData,3) /= TmpNumSteps ) THEN

               ! Let's make the array the correct size (we should never get here, but you never know)

            DEALLOCATE( ParamData%FFData )

            CALL AllocAry( ParamData%FFData, ParamData%NZGrids,ParamData%NYGrids,ParamData%NFFComp,TmpNumSteps, &
                  'Full-field wind data array.', TmpErrStat, TmpErrMsg )
               CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
               IF ( ErrStat >= AbortErrLev ) RETURN            

         ENDIF !Incorrect size
      ENDIF ! allocated

      !-------------------------------------------------------------------------------------------------
      ! Initialize the data and set column indexing to account for direction of turbine rotation (CWise)
      !-------------------------------------------------------------------------------------------------

      ParamData%FFData(:,:,:,:) = 0.0                        ! we may have only one component

      IF ( CWise )  THEN
         CFirst    = ParamData%NYGrids
         CLast     = 1
         CStep     = -1
      ELSE
         CFirst    = 1
         CLast     = ParamData%NYGrids
         CStep     = 1
      ENDIF


      !-------------------------------------------------------------------------------------------------
      ! Loop through all the time steps, reading the data and converting to m/s
      !-------------------------------------------------------------------------------------------------
   !bjj: should we reorganize this FFData array so we access the data faster?

      ParamData%NFFSteps = TmpNumSteps

   TIME_LOOP:  DO IT=1,TmpNumSteps     ! time (add 1 to see if there is an odd number of grids)

         DO IR=1,ParamData%NZGrids               ! the rows (vertical)

            DO IC=CFirst,CLast,CStep   ! the columns (lateral)

               DO I=1,ParamData%NFFComp          ! wind components (U, V, W)

                     ! Get the next integer from the file.
                     ! This is a 2-byte integer, so we can't use the library read routines.
                  READ (UnitWind,IOStat=TmpErrStat)  Dum_Int2
                  IF (TmpErrStat /= 0) THEN
                     IF ( IT == TmpNumSteps ) THEN ! There really were an even number of steps
                        ParamData%NFFSteps = TmpNumSteps - 1
                        ErrStat  = 0
                        EXIT TIME_LOOP
                     ELSE
                        CALL SetErrStat( ErrID_Fatal, ' Error reading binary data file. '// &
                                    'ic = '//TRIM(Num2LStr(ic))// &
                                    ', ir = '//TRIM(Num2LStr(ir))// &
                                    ', it = '//TRIM(Num2LStr(it))// &
                                    ', nffsteps = '//TRIM(Num2LStr(ParamData%NFFSteps)), ErrStat, ErrMsg, RoutineName)
                        RETURN
                     ENDIF
                  ELSE
                     ParamData%FFData(IR,IC,I,IT) = ParamData%MeanFFWS*(FF_Offset(I)+0.00001*TI(I)*Dum_Int2)
                  ENDIF

               END DO !I

            END DO !IC

         END DO !IR

      END DO TIME_LOOP !IT

      IF ( ParamData%Periodic ) THEN
         TmpErrMsg = '   Processed '//TRIM( Num2LStr( ParamData%NFFSteps ) )//' time steps of '// &
                    TRIM( Num2LStr ( ParamData%FFRate ) )//'-Hz full-field data (period of '// &
                    TRIM( Num2LStr( ParamData%FFDTime*ParamData%NFFSteps ) )//' seconds).'
         
      ELSE
         TmpErrMsg= '   Processed '//TRIM( Num2LStr( ParamData%NFFSteps ) )//' time steps of '// &
                    TRIM( Num2LStr ( ParamData%FFRate ) )//'-Hz full-field data ('// &
                    TRIM( Num2LStr( ParamData%FFDTime*( ParamData%NFFSteps - 1 ) ) )//' seconds).'
      ENDIF
      CALL WrScr( NewLine//TRIM(TmpErrMsg) )
      !CALL SetErrStat( ErrID_Info, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
      
      

   END SUBROUTINE Read_Bladed_Grids
   !====================================================================================================
   !> This subroutine reads the binary tower file that corresponds with the Bladed-style FF binary file.
   !! The FF grid must be read before this subroutine is called! (many checks are made to ensure the
   !! files belong together)
   !!   16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   SUBROUTINE Read_FF_Tower( MiscVars, ErrStat, ErrMsg )
      IMPLICIT                                              NONE

      CHARACTER(*),           PARAMETER                  :: RoutineName="Read_FF_Tower"


         ! Passed Variables:

      TYPE(IfW_BladedFFWind_MiscVarType), INTENT(INOUT)  :: MiscVars       !< misc/optimization data (storage for the main data)
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

      ParamData%NTGrids = 0

      IF ( ParamData%NFFComp /= 3 ) THEN
         CALL SetErrStat( ErrID_Fatal, ' Error: Tower binary files require 3 wind components.', ErrStat, ErrMsg, RoutineName )
         RETURN
      ENDIF

      !-------------------------------------------------------------------------------------------------
      ! Open the file
      !-------------------------------------------------------------------------------------------------

      CALL OpenBInpFile (UnitWind, TRIM(TwrFile), TmpErrStat, TmpErrMsg)
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN

      !-------------------------------------------------------------------------------------------------
      ! Read the header information and check that it's compatible with the FF Bladed-style binary
      ! parameters already read.
      !-------------------------------------------------------------------------------------------------
            ! This is a 4-byte real, so we can't use the library read routines.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4               ! dz, in meters [4-byte REAL]
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading dz in the binary tower file "'//TRIM( TwrFile )//'."', ErrStat, ErrMsg, RoutineName )
               RETURN
            ENDIF

            IF ( ABS(Dum_Real4*ParamData%InvFFZD-1) > TOL ) THEN
               CALL SetErrStat( ErrID_Fatal, ' Resolution in the FF binary file does not match the tower file.', ErrStat, ErrMsg, RoutineName )
               RETURN
            ENDIF


            ! This is a 4-byte real, so we can't use the library read routines.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4               ! dx, in meters [4-byte REAL]
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading dx in the binary tower file "'//TRIM( TwrFile )//'."', ErrStat, ErrMsg, RoutineName )
               RETURN
            ENDIF

            IF ( ABS(Dum_Real4*ParamData%InvMFFWS/ParamData%FFDTime-1) > TOL ) THEN
               CALL SetErrStat( ErrID_Fatal, ' Time resolution in the FF binary file does not match the tower file.', ErrStat, ErrMsg, RoutineName )
               RETURN
            ENDIF


            ! This is a 4-byte real, so we can't use the library read routines.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4               ! Zmax, in meters [4-byte REAL]
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading Zmax in the binary tower file "'//TRIM( TwrFile )//'."', ErrStat, ErrMsg, RoutineName )
               RETURN
            ENDIF

            IF ( ABS(Dum_Real4/ParamData%GridBase-1) > TOL ) THEN
               CALL SetErrStat( ErrID_Fatal, ' Height in the FF binary file does not match the tower file "'//TRIM( TwrFile )//'."', ErrStat, ErrMsg, RoutineName )
               RETURN
            ENDIF


            ! This is a 4-byte integer, so we can't use the library read routines.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                ! NumOutSteps [4-byte INTEGER]
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading NumOutSteps in the binary tower file "'//TRIM( TwrFile )//'."', ErrStat, ErrMsg, RoutineName )
               RETURN
            ENDIF

            IF ( Dum_Int4 /= ParamData%NFFSteps ) THEN
               CALL SetErrStat( ErrID_Fatal, ' Number of time steps in the FF binary file does not match the tower file.', ErrStat, ErrMsg, RoutineName )
               RETURN
            ENDIF


            ! This is a 4-byte integer, so we can't use the library read routines.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                ! NumZ      [4-byte INTEGER]
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading NumZ in the binary tower file "'//TRIM( TwrFile )//'."', ErrStat, ErrMsg, RoutineName )
               RETURN
            ENDIF
            ParamData%NTGrids = Dum_Int4


            ! This is a 4-byte real, so we can't use the library read routines.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4               ! UHub      [4-byte REAL]
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading UHub in the binary tower file "'//TRIM( TwrFile )//'."', ErrStat, ErrMsg, RoutineName )
               RETURN
            ENDIF

            IF ( ABS(Dum_Real4*ParamData%InvMFFWS - 1) > TOL ) THEN
               CALL SetErrStat( ErrID_Fatal, ' Mean wind speed in the FF binary file does not match the tower file.', ErrStat, ErrMsg, RoutineName )
               ParamData%NTGrids  = 0
               RETURN
            ENDIF


         DO IC=1,3
               ! Read the TI values fromthe tower file: 4-byte reals.
               
               !bjj: not sure you can call this routine to read from a binary file...
            !CALL ReadVar( UnitWind, TRIM(TwrFile), TI(IC), 'TI('//TRIM(Num2LStr(IC))//')', 'TI value for u,v, or w', TmpErrStat, TmpErrMsg )
            !IF (TmpErrStat /= ErrID_None) THEN
            !   ParamData%NTGrids  = 0
            !   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
            !   IF (ErrStat >= AbortErrLev) RETURN
            !ENDIF
            !
            READ (UnitWind, IOSTAT=TmpErrStat)   TI(IC)               ! TI(u), TI(v), TI(w)  [4-byte REAL]
            
            IF (TmpErrStat /= 0) THEN
               ParamData%NTGrids  = 0
               CALL SetErrStat( ErrID_Fatal, ' Error reading TI('//TRIM(Num2LStr(IC))//') in the binary tower file "' &
                               //TRIM( TwrFile )//'."', ErrStat, ErrMsg, RoutineName )
               RETURN
            ENDIF
         
         END DO

      !----------------------------------------------------------------------------------------------
      ! Allocate arrays for the tower points
      !----------------------------------------------------------------------------------------------

         IF ( ParamData%NTGrids > 0 ) THEN

            IF ( .NOT. ALLOCATED( ParamData%FFTower ) ) THEN
               CALL AllocAry( ParamData%FFTower, ParamData%NFFComp, ParamData%NTGrids, ParamData%NFFSteps, &
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

         DO IT=1,ParamData%NFFSteps

            DO IZ=1,ParamData%NTGrids         ! If NTGrids<1, there are no tower points & FFTower is not allocated

               ! Ytower     = 0               ! Lateral location of the tower data point, in m relative to tower centerline
               ! Ztower(IZ) = Z1 - (IZ-1)*dz  ! Vertical location of tower data point, in m relative to ground

               DO IC=1,ParamData%NFFComp   ! number of wind components

                     ! Read in the 2-byte integer. Can't use library read routines for this.
                  READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2       ! normalized wind-component, INT(2)
                  IF ( TmpErrStat /= 0 )  THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading binary tower data file. it = '//TRIM(Num2LStr(it))// &
                                    ', nffsteps = '//TRIM(Num2LStr(ParamData%NFFSteps)), ErrStat, ErrMsg, RoutineName )
                     ParamData%NTGrids  = 0
                     RETURN
                  ENDIF

                  ParamData%FFTower(IC,IZ,IT) = ParamData%MeanFFWS*(FF_Offset(IC)+0.00001*TI(IC)*Dum_Int2)   ! wind-component scaled to m/s

               ENDDO !IC

            ENDDO ! IZ


         ENDDO ! IT

      !-------------------------------------------------------------------------------------------------
      ! Close the file
      !-------------------------------------------------------------------------------------------------
      CLOSE ( UnitWind )

      TmpErrMsg = '   Processed '//TRIM( Num2LStr(ParamData%NFFSteps) )//' time steps of '// &
            TRIM( Num2LStr(ParamData%NTGrids) )//'x1 tower data grids.'
      
      !CALL SetErrStat( ErrID_Info, ErrMsgLcl, ErrStat, ErrMsg, RoutineName )
      CALL WrScr( NewLine//TRIM(TmpErrMsg) )
      
      RETURN

   END SUBROUTINE Read_FF_Tower
END SUBROUTINE IfW_BladedFFWind_Init
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

      ! local variables
   INTEGER(IntKi)                                              :: NumPoints         ! Number of points specified by the PositionXYZ array

      ! local counters
   INTEGER(IntKi)                                              :: PointNum          ! a loop counter for the current point

      ! temporary variables
   INTEGER(IntKi)                                              :: TmpErrStat        ! temporary error status
   CHARACTER(ErrMsgLen)                                        :: TmpErrMsg         ! temporary error message


      !-------------------------------------------------------------------------------------------------
      ! Check that the module has been initialized.
      !-------------------------------------------------------------------------------------------------

   ErrStat     = ErrID_None
   ErrMsg      = ''

      !-------------------------------------------------------------------------------------------------
      ! Initialize some things
      !-------------------------------------------------------------------------------------------------


      ! The array is transposed so that the number of points is the second index, x/y/z is the first.
      ! This is just in case we only have a single point, the SIZE command returns the correct number of points.
   NumPoints   =  SIZE(PositionXYZ,2)


      ! Step through all the positions and get the velocities
   DO PointNum = 1, NumPoints

         ! Calculate the velocity for the position
      Velocity(:,PointNum) = FF_Interp(Time,PositionXYZ(:,PointNum),ParamData,MiscVars,TmpErrStat,TmpErrMsg)

         ! Error handling
      IF (TmpErrStat /= ErrID_None) THEN  !  adding this so we don't have to convert numbers to strings every time
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, "IfW_BladedFFWind:CalcOutput [position=("//   &
                                                      TRIM(Num2LStr(PositionXYZ(1,PointNum)))//", "// &
                                                      TRIM(Num2LStr(PositionXYZ(2,PointNum)))//", "// &
                                                      TRIM(Num2LStr(PositionXYZ(3,PointNum)))//")]" )
         IF (ErrStat >= AbortErrLev) RETURN
      END IF

   ENDDO



      !REMOVE THIS for AeroDyn 15
      ! Return the average disk velocity values needed by AeroDyn 14.  This is the WindInf_ADhack_diskVel routine.
   DiskVel(1)   =  ParamData%MeanFFWS
   DiskVel(2:3) =  0.0_ReKi


   RETURN

END SUBROUTINE IfW_BladedFFWind_CalcOutput
   !+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
   !>    This function is used to interpolate into the full-field wind array or tower array if it has
   !!    been defined and is necessary for the given inputs.  It receives X, Y, Z and
   !!    TIME from the calling routine.  It then computes a time shift due to a nonzero X based upon
   !!    the average windspeed.  The modified time is used to decide which pair of time slices to interpolate
   !!    within and between.  After finding the two time slices, it decides which four grid points bound the
   !!    (Y,Z) pair.  It does a bilinear interpolation for each time slice. Linear interpolation is then used
   !!    to interpolate between time slices.  This routine assumes that X is downwind, Y is to the left when
   !!    looking downwind and Z is up.  It also assumes that no extrapolation will be needed.
   !!
   !!    If tower points are used, it assumes the velocity at the ground is 0.  It interpolates between
   !!    heights and between time slices, but ignores the Y input.
   !!
   !!    11/07/1994 - Created by M. Buhl from the original TURBINT.
   !!    09/25/1997 - Modified by M. Buhl to use f90 constructs and new variable names.  Renamed to FF_Interp.
   !!    09/23/2009 - Modified by B. Jonkman to use arguments instead of modules to determine time and position.
   !!                 Height is now relative to the ground
   !!   16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   FUNCTION FF_Interp(Time, Position, ParamData, MiscVars, ErrStat, ErrMsg)

      IMPLICIT                                              NONE

      CHARACTER(*),           PARAMETER                     :: RoutineName="FF_Interp"

      REAL(DbKi),                            INTENT(IN   )  :: Time           !< time (s)
      REAL(ReKi),                            INTENT(IN   )  :: Position(3)    !< takes the place of XGrnd, YGrnd, ZGrnd
      TYPE(IfW_BladedFFWind_ParameterType),  INTENT(IN   )  :: ParamData      !< Parameters
      TYPE(IfW_BladedFFWind_MiscVarType),    INTENT(INOUT)  :: MiscVars       !< misc/optimization data (storage for the main data)
      REAL(ReKi)                                            :: FF_Interp(3)   !< The U, V, W velocities

      INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat        !< error status
      CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg         !< error message 

         ! Local Variables:

      REAL(ReKi)                                            :: TimeShifted
      REAL(ReKi),PARAMETER                                  :: Tol = 1.0E-3   ! a tolerance for determining if two reals are the same (for extrapolation)
      REAL(ReKi)                                            :: T
      REAL(ReKi)                                            :: TGRID
      REAL(ReKi)                                            :: Y
      REAL(ReKi)                                            :: YGRID
      REAL(ReKi)                                            :: Z
      REAL(ReKi)                                            :: ZGRID
      REAL(ReKi)                                            :: N(8)           ! array for holding scaling factors for the interpolation algorithm
      REAL(ReKi)                                            :: u(8)           ! array for holding the corner values for the interpolation algorithm across a cubic volume
      REAL(ReKi)                                            :: M(4)           ! array for holding scaling factors for the interpolation algorithm
      REAL(ReKi)                                            :: v(4)           ! array for holding the corner values for the interpolation algorithm across an area

      INTEGER(IntKi)                                        :: IDIM
      INTEGER(IntKi)                                        :: ITHI
      INTEGER(IntKi)                                        :: ITLO
      INTEGER(IntKi)                                        :: IYHI
      INTEGER(IntKi)                                        :: IYLO
      INTEGER(IntKi)                                        :: IZHI
      INTEGER(IntKi)                                        :: IZLO

      LOGICAL                                               :: OnGrid

      !-------------------------------------------------------------------------------------------------
      ! Initialize variables
      !-------------------------------------------------------------------------------------------------

      FF_Interp(:)        = 0.0_ReKi                         ! the output velocities (in case ParamData%NFFComp /= 3)

      ErrStat              = ErrID_None
      ErrMsg               = ""
      
      !-------------------------------------------------------------------------------------------------
      ! Find the bounding time slices.
      !-------------------------------------------------------------------------------------------------

      ! Perform the time shift.  At time=0, a point half the grid width downstream (ParamData%FFYHWid) will index into the zero time slice.
      ! If we did not do this, any point downstream of the tower at the beginning of the run would index outside of the array.
      ! This all assumes the grid width is at least as large as the rotor.  If it isn't, then the interpolation will not work.


      TimeShifted = TIME + ( ParamData%InitXPosition - Position(1) )*ParamData%InvMFFWS    ! in distance, X: InputInfo%Position(1) - ParamData%InitXPosition - TIME*ParamData%MeanFFWS


      IF ( ParamData%Periodic ) THEN ! translate TimeShifted to ( 0 <= TimeShifted < ParamData%TotalTime )

         TimeShifted = MODULO( TimeShifted, ParamData%TotalTime )
             ! If TimeShifted is a very small negative number, modulo returns the incorrect value due to internal rounding errors.
             ! See bug report #471
         IF (TimeShifted == ParamData%TotalTime) TimeShifted = 0.0_ReKi

         TGRID = TimeShifted*ParamData%FFRate
         ITLO  = INT( TGRID )             ! convert REAL to INTEGER (add 1 later because our grids start at 1, not 0)
         T     = 2.0_ReKi * ( TGRID - REAL(ITLO, ReKi) ) - 1.0_ReKi     ! a value between -1 and 1 that indicates a relative position between ITLO and ITHI

         ITLO = ITLO + 1
         IF ( ITLO == ParamData%NFFSteps ) THEN
            ITHI = 1
         ELSE
            ITHI = ITLO + 1
         ENDIF


      ELSE

         TGRID = TimeShifted*ParamData%FFRate
         ITLO  = INT( TGRID )             ! convert REAL to INTEGER (add 1 later because our grids start at 1, not 0)
         T     = 2.0_ReKi * ( TGRID - REAL(ITLO, ReKi) ) - 1.0_ReKi     ! a value between -1 and 1 that indicates a relative position between ITLO and ITHI

         ITLO = ITLO + 1                  ! add one since our grids start at 1, not 0
         ITHI = ITLO + 1

         IF ( ITLO >= ParamData%NFFSteps .OR. ITLO < 1 ) THEN
            IF ( ITLO == ParamData%NFFSteps  ) THEN
               ITHI = ITLO
               IF ( T <= TOL ) THEN ! we're on the last point
                  T = -1.0_ReKi
               ELSE  ! We'll extrapolate one dt past the last value in the file
                  ITLO = ITHI - 1
               ENDIF
            ELSE
               ErrMsg   = ' Error: FF wind array was exhausted at '//TRIM( Num2LStr( REAL( TIME,   ReKi ) ) )// &
                          ' seconds (trying to access data at '//TRIM( Num2LStr( REAL( TimeShifted, ReKi ) ) )//' seconds).'
               ErrStat  = ErrID_Fatal
               RETURN
            ENDIF
         ENDIF

      ENDIF


      !-------------------------------------------------------------------------------------------------
      ! Find the bounding rows for the Z position. [The lower-left corner is (1,1) when looking upwind.]
      !-------------------------------------------------------------------------------------------------

      ZGRID = ( Position(3) - ParamData%GridBase )*ParamData%InvFFZD

      IF (ZGRID > -1*TOL) THEN
         OnGrid = .TRUE.

            ! Index for start and end slices
         IZLO = INT( ZGRID ) + 1             ! convert REAL to INTEGER, then add one since our grids start at 1, not 0
         IZHI = IZLO + 1

            ! Set Z as a value between -1 and 1 for the relative location between IZLO and IZHI.
            ! Subtract 1_IntKi from Z since the indices are starting at 1, not 0
         Z = 2.0_ReKi * (ZGRID - REAL(IZLO - 1_IntKi, ReKi)) - 1.0_ReKi

         IF ( IZLO < 1 ) THEN
            IF ( IZLO == 0 .AND. Z >= 1.0-TOL ) THEN
               Z    = -1.0_ReKi
               IZLO = 1
            ELSE
               ErrMsg   = ' FF wind array boundaries violated. Grid too small in Z direction (Z='//&
                          TRIM(Num2LStr(Position(3)))//' m is below the grid).'
               ErrStat  = ErrID_Fatal
               RETURN
            ENDIF
         ELSEIF ( IZLO >= ParamData%NZGrids ) THEN
            IF ( IZLO == ParamData%NZGrids .AND. Z <= TOL ) THEN
               Z    = -1.0_ReKi
               IZHI = IZLO                   ! We're right on the last point, which is still okay
            ELSE
               ErrMsg   = ' FF wind array boundaries violated. Grid too small in Z direction (Z='//&
                          TRIM(Num2LStr(Position(3)))//' m is above the grid).'
               ErrStat  = ErrID_Fatal
               RETURN
            ENDIF
         ENDIF

      ELSE

         OnGrid = .FALSE.  ! this is on the tower

         IF ( ParamData%NTGrids < 1 ) THEN
            ErrMsg   = ' FF wind array boundaries violated. Grid too small in Z direction '// &
                       '(height (Z='//TRIM(Num2LStr(Position(3)))//' m) is below the grid and no tower points are defined).'
            ErrStat  = ErrID_Fatal
            RETURN
         ENDIF

         IZLO = INT( -1.0*ZGRID ) + 1            ! convert REAL to INTEGER, then add one since our grids start at 1, not 0


         IF ( IZLO >= ParamData%NTGrids ) THEN  !our dz is the difference between the bottom tower point and the ground
            IZLO  = ParamData%NTGrids

               ! Check that this isn't zero.  Value between -1 and 1 corresponding to the relative position.
            Z = 1.0_ReKi - 2.0_ReKi * (Position(3) / (ParamData%GridBase - REAL(IZLO - 1_IntKi, ReKi)/ParamData%InvFFZD))

         ELSE

               ! Set Z as a value between -1 and 1 for the relative location between IZLO and IZHI.  Used in the interpolation.
            Z = 2.0_ReKi * (ABS(ZGRID) - REAL(IZLO - 1_IntKi, ReKi)) - 1.0_ReKi

         ENDIF
         IZHI = IZLO + 1

      ENDIF


      IF ( OnGrid ) THEN      ! The tower points don't use this

         !-------------------------------------------------------------------------------------------------
         ! Find the bounding columns for the Y position. [The lower-left corner is (1,1) when looking upwind.]
         !-------------------------------------------------------------------------------------------------

            YGRID = ( Position(2) + ParamData%FFYHWid )*ParamData%InvFFYD    ! really, it's (Position(2) - -1.0*ParamData%FFYHWid)

            IYLO = INT( YGRID ) + 1             ! convert REAL to INTEGER, then add one since our grids start at 1, not 0
            IYHI = IYLO + 1

               ! Set Y as a value between -1 and 1 for the relative location between IYLO and IYHI.  Used in the interpolation.
               ! Subtract 1_IntKi from IYLO since grids start at index 1, not 0
            Y = 2.0_ReKi * (YGRID - REAL(IYLO - 1_IntKi, ReKi)) - 1.0_ReKi

            IF ( IYLO >= ParamData%NYGrids .OR. IYLO < 1 ) THEN
               IF ( IYLO == 0 .AND. Y >= 1.0-TOL ) THEN
                  Y    = -1.0_ReKi
                  IYLO = 1
               ELSE IF ( IYLO == ParamData%NYGrids .AND. Y <= TOL ) THEN
                  Y    = -1.0_ReKi
                  IYHI = IYLO                   ! We're right on the last point, which is still okay
               ELSE
                  ErrMsg   = ' FF wind array boundaries violated: Grid too small in Y direction. Y='// &
                             TRIM(Num2LStr(Position(2)))//'; Y boundaries = ['//TRIM(Num2LStr(-1.0*ParamData%FFYHWid))// &
                             ', '//TRIM(Num2LStr(ParamData%FFYHWid))//']'
                  ErrStat = ErrID_Fatal         ! we don't return anything
                  RETURN
               ENDIF
            ENDIF

         !-------------------------------------------------------------------------------------------------
         ! Interpolate on the grid
         !-------------------------------------------------------------------------------------------------

         DO IDIM=1,ParamData%NFFComp       ! all the components


!New Algorithm here
            N(1)  = ( 1.0_ReKi + Z )*( 1.0_ReKi - Y )*( 1.0_ReKi - T )
            N(2)  = ( 1.0_ReKi + Z )*( 1.0_ReKi + Y )*( 1.0_ReKi - T )
            N(3)  = ( 1.0_ReKi - Z )*( 1.0_ReKi + Y )*( 1.0_ReKi - T )
            N(4)  = ( 1.0_ReKi - Z )*( 1.0_ReKi - Y )*( 1.0_ReKi - T )
            N(5)  = ( 1.0_ReKi + Z )*( 1.0_ReKi - Y )*( 1.0_ReKi + T )
            N(6)  = ( 1.0_ReKi + Z )*( 1.0_ReKi + Y )*( 1.0_ReKi + T )
            N(7)  = ( 1.0_ReKi - Z )*( 1.0_ReKi + Y )*( 1.0_ReKi + T )
            N(8)  = ( 1.0_ReKi - Z )*( 1.0_ReKi - Y )*( 1.0_ReKi + T )
            N     = N / REAL( SIZE(N), ReKi )  ! normalize


            u(1)  = ParamData%FFData( IZHI, IYLO, IDIM, ITLO )
            u(2)  = ParamData%FFData( IZHI, IYHI, IDIM, ITLO )
            u(3)  = ParamData%FFData( IZLO, IYHI, IDIM, ITLO )
            u(4)  = ParamData%FFData( IZLO, IYLO, IDIM, ITLO )
            u(5)  = ParamData%FFData( IZHI, IYLO, IDIM, ITHI )
            u(6)  = ParamData%FFData( IZHI, IYHI, IDIM, ITHI )
            u(7)  = ParamData%FFData( IZLO, IYHI, IDIM, ITHI )
            u(8)  = ParamData%FFData( IZLO, IYLO, IDIM, ITHI )
            
            FF_Interp(IDIM)  =  SUM ( N * u ) 


         END DO !IDIM

      ELSE

      !-------------------------------------------------------------------------------------------------
      ! Interpolate on the tower array
      !-------------------------------------------------------------------------------------------------

         DO IDIM=1,ParamData%NFFComp    ! all the components

            !----------------------------------------------------------------------------------------------
            ! Interpolate between the two times using an area interpolation.
            !----------------------------------------------------------------------------------------------

               ! Setup the scaling factors.  Set the unused portion of the array to zero
            M(1)  =  ( 1.0_ReKi + Z )*( 1.0_ReKi - T )
            M(2)  =  ( 1.0_ReKi + Z )*( 1.0_ReKi + T )
            M(3)  =  ( 1.0_ReKi - Z )*( 1.0_ReKi - T )
            M(4)  =  ( 1.0_ReKi - Z )*( 1.0_ReKi + T )
            M     =  M / 4.0_ReKi               ! normalize

            IF (IZHI > ParamData%NTGrids) THEN
               v(1)  =  0.0_ReKi  ! on the ground
               v(2)  =  0.0_ReKi  ! on the ground
            ELSE
               v(1)  =  ParamData%FFTower( IDIM, IZHI, ITLO )
               v(2)  =  ParamData%FFTower( IDIM, IZHI, ITHI )
            END IF
            
            v(3)  =  ParamData%FFTower( IDIM, IZLO, ITLO )
            v(4)  =  ParamData%FFTower( IDIM, IZLO, ITHI )
            
            FF_Interp(IDIM)  =  SUM ( M * v ) 


         END DO !IDIM

      ENDIF ! OnGrid
      RETURN

   END FUNCTION FF_Interp


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
