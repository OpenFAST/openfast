!>  This module uses full-field binary wind files to determine the wind inflow.
!!  This module assumes that the origin, (0,0,0), is located at the tower centerline at ground level,
!!  and that all units are specified in the metric system (using meters and seconds).
!!  Data is shifted by half the grid width to account for turbine yaw (so that data in the X
!!  direction actually starts at -1*ParamData%FFYHWid meters).
MODULE IfW_TSFFWind
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
   USE                                          IfW_TSFFWind_Types

   IMPLICIT                                     NONE
   PRIVATE

   TYPE(ProgDesc),   PARAMETER               :: IfW_TSFFWind_Ver = ProgDesc( 'IfW_TSFFWind', '', '' )

   PUBLIC                                    :: IfW_TSFFWind_Init
   PUBLIC                                    :: IfW_TSFFWind_End
   PUBLIC                                    :: IfW_TSFFWind_CalcOutput




CONTAINS
!====================================================================================================
!>  This routine is used read the full-field turbulence data.
!!  09/25/1997  - Created by M. Buhl from GETFILES in ViewWind.
!!  09/23/2009  - modified by B. Jonkman: this subroutine was split into several subroutines (was ReadFF)
!!  16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
SUBROUTINE IfW_TSFFWind_Init(InitData, ParamData, MiscVars, Interval, InitOutData, ErrStat, ErrMsg)

   IMPLICIT                       NONE

   CHARACTER(*),           PARAMETER                        :: RoutineName="IfW_TSFFWind_Init"


      ! Passed Variables
   TYPE(IfW_TSFFWind_InitInputType),         INTENT(IN   )  :: InitData          !< Initialization data passed to the module
   TYPE(IfW_TSFFWind_ParameterType),         INTENT(  OUT)  :: ParamData         !< Parameters
   TYPE(IfW_TSFFWind_MiscVarType),           INTENT(  OUT)  :: MiscVars          !< Misc variables for optimization (not copied in glue code)
   TYPE(IfW_TSFFWind_InitOutputType),        INTENT(  OUT)  :: InitOutData       !< Initial output

   REAL(DbKi),                               INTENT(IN   )  :: Interval          !< Time Interval to use (passed through here)


      ! Error Handling
   INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat           !< determines if an error has been encountered
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg            !< Message about errors


      ! Temporary variables for error handling
   INTEGER(IntKi)                                           :: TmpErrStat        ! temporary error status
   CHARACTER(ErrMsgLen)                                     :: TmpErrMsg         ! temporary error message


      ! Local Variables:

   INTEGER(IntKi)                                           :: UnitWind     ! Unit number for the InflowWind input file
   INTEGER(B2Ki)                                            :: Dum_Int2

      !-------------------------------------------------------------------------------------------------
      ! Initialize temporary variables
      !-------------------------------------------------------------------------------------------------

   ErrMsg      = ''
   ErrStat     = ErrID_None

   TmpErrMsg   = ''
   TmpErrStat  = ErrID_None



      ! Get a unit number to use

   CALL GetNewUnit(UnitWind, TmpErrStat, TmpErrMsg)
   CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
   IF ( ErrStat >= AbortErrLev ) RETURN


      !-------------------------------------------------------------------------------------------------
      ! Copy things from the InitData to the ParamData
      !-------------------------------------------------------------------------------------------------


      !----------------------------------------------------------------------------------------------
      ! Open the binary file, read its "header" (first 2-byte integer) to determine what format
      ! binary file it is, and close it.
      !----------------------------------------------------------------------------------------------

   CALL OpenBInpFile (UnitWind, TRIM(InitData%WindFileName), TmpErrStat, TmpErrMsg)
   CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
   IF ( ErrStat >= AbortErrLev ) RETURN

      ! Read the first binary integer from the file to get info on the type.
      ! Cannot use library read routines since this is a 2-byte integer.
   READ ( UnitWind, IOSTAT=TmpErrStat )  Dum_Int2
   CLOSE( UnitWind )

   IF (TmpErrStat /= 0) THEN
      CALL SetErrStat(ErrID_Fatal,' Error reading first binary integer from file "'//TRIM(InitData%WindFileName)//'."',   &
               ErrStat,ErrMsg,RoutineName)
      RETURN
   ENDIF


      !----------------------------------------------------------------------------------------------
      ! Read the files to get the required FF data.
      !----------------------------------------------------------------------------------------------
      ! Store the binary format information so the InflowWind code can use it.
      ! Also changes to IntKi from INT(2) to compare in the SELECT below
   ParamData%WindFileFormat = Dum_Int2

   SELECT CASE (ParamData%WindFileFormat)

      CASE ( 7, 8 )                                                    ! TurbSim binary format

         CALL Read_TurbSim_FF(UnitWind, TmpErrStat, TmpErrMsg)
            CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )                  
            IF ( ErrStat >= AbortErrLev ) THEN
               CLOSE ( UnitWind )
               RETURN
            END IF
         

      CASE DEFAULT
         CALL SetErrStat( ErrID_Fatal, ' Error: This is not a TurbSim binary wind file type (binary format identifier: '//  &
                  TRIM(Num2LStr(ParamData%WindFileFormat))//'.  This might be a Bladed style binary wind file.', &
                  ErrStat, ErrMsg, RoutineName )
         RETURN

   END SELECT


   IF (ParamData%Periodic) THEN
      ParamData%InitXPosition = 0                ! start at the hub
      ParamData%TotalTime     = ParamData%NFFSteps*ParamData%FFDTime
   ELSE
      ParamData%InitXPosition = ParamData%FFYHWid          ! start half the grid with ahead of the turbine
      ParamData%TotalTime     = (ParamData%NFFSteps-1)*ParamData%FFDTime
   ENDIF


      !-------------------------------------------------------------------------------------------------
      ! Set the InitOutput information
      !-------------------------------------------------------------------------------------------------

   InitOutdata%Ver         = IfW_TSFFWind_Ver



   RETURN


   CONTAINS

   !====================================================================================================
   !> This subroutine reads the binary TurbSim-format FF file (.bts).  It fills the FFData array with
   !! velocity data for the grids and fills the FFTower array with velocities at points on the tower
   !! (if data exists).
   !!   16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   SUBROUTINE Read_TurbSim_FF(UnitWind, ErrStat, ErrMsg)
   !----------------------------------------------------------------------------------------------------

      CHARACTER(*),           PARAMETER                  :: RoutineName="READ_TurbSim_FF"


         ! Passed Variables:

      INTEGER(IntKi),                     INTENT(IN   )  :: UnitWind       !< unit number for the wind file
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat        !< error status return value (0=no error; non-zero is error)
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg         !< message about the error encountered

         ! Local Variables:

      REAL(SiKi)                                         :: Dum_Real4      ! dummy 4-byte real number
      INTEGER(B1Ki)                                      :: Dum_Int1       ! dummy 1-byte integer
      INTEGER(B2Ki)                                      :: Dum_Int2       ! dummy 2-byte integer
      INTEGER(B4Ki)                                      :: Dum_Int4       ! dummy 4-byte integer

      INTEGER(IntKi)                                     :: IC             ! loop counter for wind components
      INTEGER(IntKi)                                     :: IT             ! loop counter for time
      INTEGER(IntKi)                                     :: IY             ! loop counter for y
      INTEGER(IntKi)                                     :: IZ             ! loop counter for z
      INTEGER(IntKi)                                     :: NChar          ! number of characters in the description string

      REAL(SiKi)                                         :: Vslope(3)      ! slope  for "un-normalizing" data
      REAL(SiKi)                                         :: Voffset(3)     ! offset for "un-normalizing" data

      LOGICAL                                            :: FirstWarn      ! we don't need to print warning for each character that exceeds the 
      CHARACTER(1024)                                    :: DescStr        ! description string contained in the file

      
         ! Temporary variables for error handling
      INTEGER(IntKi)                                     :: TmpErrStat     ! temporary error status
      CHARACTER(ErrMsgLen)                               :: TmpErrMsg      ! temporary error message


      ParamData%NFFComp = 3                                              ! this file contains 3 wind components
      ErrStat = ErrID_None
      ErrMsg  = ""
      
   !-------------------------------------------------------------------------------------------------
   ! Open the file
   !-------------------------------------------------------------------------------------------------

      CALL OpenBInpFile (UnitWind, TRIM(InitData%WindFileName), TmpErrStat, TmpErrMsg)
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )                  
      IF ( ErrStat >= AbortErrLev ) RETURN

      !-------------------------------------------------------------------------------------------------
      ! Read the header information
      !-------------------------------------------------------------------------------------------------
            ! Read in the 2-byte integer. Can't use library read routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2             ! the file identifier, INT(2)
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading the file identifier in the FF binary file "'//TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName )                  
               RETURN
            ENDIF
            ParamData%Periodic = Dum_Int2 == INT( 8, B2Ki) ! the number 7 is used for non-periodic wind files; 8 is periodic wind


            ! Read in the 4-byte integer. Can't use library read routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4             ! the number of grid points vertically, INT(4)
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading the number of z grid points in the FF binary file "'//TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName )                  
               RETURN
            ENDIF
            ParamData%NZGrids = Dum_Int4


            ! Read in the 4-byte integer. Can't use library read routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4             ! the number of grid points laterally, INT(4)
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading the number of y grid points in the FF binary file "'//TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName )                  
               RETURN
            ENDIF
            ParamData%NYGrids = Dum_Int4


            ! Read in the 4-byte integer. Can't use library read routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4             ! the number of tower points, INT(4)
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading the number of tower points in the FF binary file "'//TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName )                  
               RETURN
            ENDIF
            ParamData%NTGrids = Dum_Int4


            ! Read in the 4-byte integer. Can't use library read routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4             ! the number of time steps, INT(4)
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading the number of time steps in the FF binary file "'//TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName )                  
               RETURN
            ENDIF
            ParamData%NFFSteps = Dum_Int4


            ! Read in the 4-byte real. Can't use library read routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4            ! grid spacing in vertical direction (dz), REAL(4), in m
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading dz in the FF binary file "'//TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName )                  
               RETURN
            ENDIF
            ParamData%InvFFZD = 1.0/Dum_Real4                            ! 1/dz
            ParamData%FFZHWid = 0.5*(ParamData%NZGrids-1)*Dum_Real4                ! half the grid height


            ! Read in the 4-byte real. Can't use library read routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4            ! grid spacing in lateral direction (dy), REAL(4), in m
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading dy in the FF binary file "'//TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName )                  
               RETURN
            ENDIF
            ParamData%InvFFYD = 1.0 / Dum_Real4                          ! 1/dy
            ParamData%FFYHWid = 0.5*(ParamData%NYGrids-1)*Dum_Real4                ! half grid grid width


            ! Read in the 4-byte real. Can't use library read routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4            ! grid spacing in time (dt), REAL(4), in m/s
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading dt in the FF binary file "'//TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName )                  
               RETURN
            ENDIF
            ParamData%FFDTime = Dum_Real4
            ParamData%FFRate  = 1.0/ParamData%FFDTime


            ! Read in the 4-byte real. Can't use library read routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4            ! the mean wind speed at hub height, REAL(4), in m/s
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading mean wind speed in the FF binary file "'//TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName )                  
               RETURN
            ENDIF
            ParamData%MeanFFWS = Dum_Real4
            ParamData%InvMFFWS = 1.0 / ParamData%MeanFFWS


            ! Read in the 4-byte real. Can't use library read routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4            ! height of the hub, REAL(4), in m
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading zHub in the FF binary file "'//TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName )                  
               RETURN
            ENDIF
            ParamData%RefHt = Dum_Real4


            ! Read in the 4-byte real. Can't use library read routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4            ! height of the bottom of the grid, REAL(4), in m
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading GridBase in the FF binary file "'//TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName )                  
               RETURN
            ENDIF
            ParamData%GridBase = Dum_Real4

    !        ZGOffset = ParamData%RefHt - ParamData%GridBase  - ParamData%FFZHWid


         !----------------------------------------------------------------------------------------------
         ! Read the binary scaling factors
         !----------------------------------------------------------------------------------------------

            DO IC = 1,ParamData%NFFComp
                  ! Read in the 4-byte real. Can't use library read routines for this.
               READ (UnitWind, IOSTAT=TmpErrStat)   Vslope(IC)     ! the IC-component slope for scaling, REAL(4)
                  IF ( TmpErrStat /= 0 )  THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading Vslope('//Num2LStr(IC)//') in the FF binary file "'//TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName )                  
                     RETURN
                  ENDIF


                  ! Read in the 4-byte real. Can't use library read routines for this.
               READ (UnitWind, IOSTAT=TmpErrStat)   Voffset(IC)    ! the IC-component offset for scaling, REAL(4)
                  IF ( TmpErrStat /= 0 )  THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading Voffset('//Num2LStr(IC)//') in the FF binary file "'//TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName )                  
                     RETURN
                  ENDIF

            ENDDO !IC


         !----------------------------------------------------------------------------------------------
         ! Read the description string: "Generated by TurbSim (vx.xx, dd-mmm-yyyy) on dd-mmm-yyyy at hh:mm:ss."
         !----------------------------------------------------------------------------------------------

            ! Read in the 4-byte integer. Can't use library read routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4          ! the number of characters in the description string, max 200, INT(4)
               IF ( TmpErrStat /= 0 )  THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading NCHAR in the FF binary file "'//TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName )                  
                  RETURN
               ENDIF
               nchar = Dum_Int4

            DescStr = ''                                       ! Initialize the description string
            FirstWarn = .true.
            
            DO IC=1,nchar

                  ! Read in the 1-byte integer. Can't use library read routines for this.
               READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int1       ! the ASCII integer representation of the character, INT(1)
               IF ( TmpErrStat /= 0 )  THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading description line in the FF binary file "'//TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName ) 
                  RETURN
               ENDIF

               IF ( LEN(DescStr) >= IC ) THEN
                  DescStr(IC:IC) = ACHAR( Dum_Int1 )              ! converted ASCII characters
               ELSEIF ( FirstWarn ) THEN
                  FirstWarn = .FALSE.
                  CALL SetErrStat( ErrID_Info, ' Description string was too long for variable.'//TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName ) 
               ENDIF

            ENDDO !IC


      !-------------------------------------------------------------------------------------------------
      ! Get the grid and tower velocities
      !-------------------------------------------------------------------------------------------------

         ! this could take a while, so we'll write a message indicating what's going on:
         
         CALL WrScr( NewLine//'   Reading a '//TRIM( Num2LStr(ParamData%NYGrids) )//'x'//TRIM( Num2LStr(ParamData%NZGrids) )//  &
                    ' grid ('//TRIM( Num2LStr(ParamData%FFYHWid*2) )//' m wide, '// &
                    TRIM( Num2LStr(ParamData%GridBase) )//' m to '// &
                    TRIM( Num2LStr(ParamData%GridBase+ParamData%FFZHWid*2) )//&
                    ' m above ground) with a characteristic wind speed of '// &
                    TRIM( Num2LStr(ParamData%MeanFFWS) )//' m/s. '//TRIM(DescStr) )


      !----------------------------------------------------------------------------------------------
      ! Allocate arrays for the FF grid as well as the tower points, if they exist
      !----------------------------------------------------------------------------------------------

         IF ( .NOT. ALLOCATED( ParamData%FFData ) ) THEN
            CALL AllocAry( ParamData%FFData, ParamData%NZGrids, ParamData%NYGrids, ParamData%NFFComp, ParamData%NFFSteps, &
                  'Full-field wind data array.', TmpErrStat, TmpErrMsg )
            CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName ) 
            IF ( ErrStat >= AbortErrLev ) RETURN
         ENDIF


         IF ( ParamData%NTGrids > 0 ) THEN

            ParamData%TowerDataExist   =  .TRUE.

            IF ( .NOT. ALLOCATED( ParamData%FFTower ) ) THEN
               CALL AllocAry( ParamData%FFTower, ParamData%NFFComp, ParamData%NTGrids, ParamData%NFFSteps, &
                     'Tower wind file data array.', TmpErrStat, TmpErrMsg )
               CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName ) 
               IF ( ErrStat >= AbortErrLev ) RETURN
            ENDIF

         ENDIF

      !-------------------------------------------------------------------------------------------------
      ! Read the 16-bit data and scale it to 32-bit reals
      !-------------------------------------------------------------------------------------------------

         ! Loop through time.

         DO IT=1,ParamData%NFFSteps

            !...........................................................................................
            ! Read grid data at this time step.
            !...........................................................................................

            DO IZ=1,ParamData%NZGrids
               ! Zgrid(IZ) = Z1 + (IZ-1)*dz                 ! Vertical location of grid data point, in m relative to ground

               DO IY=1,ParamData%NYGrids
                  ! Ygrid(IY) = -0.5*(ny-1)*dy + (IY-1)*dy  ! Horizontal location of grid data point, in m relative to tower centerline

                  DO IC=1,ParamData%NFFComp                           ! number of wind components (U, V, W)

                        ! Read in the 2-byte integer. Can't use library read routines for this.
                     READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2       ! normalized wind-component, INT(2)
                     IF ( TmpErrStat /= 0 )  THEN                        
                        CALL SetErrStat( ErrID_Fatal, ' Error reading grid wind components in the FF binary file "'// &
                                    TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName ) 
                        RETURN
                     ENDIF

                     ParamData%FFData(IZ,IY,IC,IT) = ( Dum_Int2 - Voffset(IC) ) / VSlope(IC)

                  ENDDO !IC

               ENDDO !IY

            ENDDO ! IZ


            !...........................................................................................
            ! Read the tower data at this time step.
            !...........................................................................................

            DO IZ=1,ParamData%NTGrids         ! If NTGrids<1, there are no tower points & FFTower is not allocated

               ! Ytower     = 0               ! Lateral location of the tower data point, in m relative to tower centerline
               ! Ztower(IZ) = Z1 - (IZ-1)*dz  ! Vertical location of tower data point, in m relative to ground

               DO IC=1,ParamData%NFFComp   ! number of wind components

                     ! Read in a 2-byte integer. Can't use library routines for this.
                  READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2       ! normalized wind-component, INT(2)
                  IF ( TmpErrStat /= 0 )  THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading tower wind components in the FF binary file "'//TRIM(InitData%WindFileName)//'."'&
                                      , ErrStat, ErrMsg, RoutineName )                    
                     RETURN
                  ENDIF

                  ParamData%FFTower(IC,IZ,IT) = ( Dum_Int2 - Voffset(IC) ) / VSlope(IC)  ! wind-component scaled to m/s

               ENDDO !IC

            ENDDO ! IZ


         ENDDO ! IT

      !-------------------------------------------------------------------------------------------------
      ! close the file and return
      !-------------------------------------------------------------------------------------------------

      CLOSE ( UnitWind )

      IF ( ParamData%Periodic ) THEN
         TmpErrMsg   = '   Processed '//TRIM( Num2LStr( ParamData%NFFSteps ) )//' time steps of '// &
                        TRIM( Num2LStr ( ParamData%FFRate ) )//'-Hz full-field data (period of '// &
                        TRIM( Num2LStr( ParamData%FFDTime*( ParamData%NFFSteps ) ) )//' seconds).'
      ELSE
         TmpErrMsg   = '   Processed '//TRIM( Num2LStr( ParamData%NFFSteps ) )//' time steps of '// &
                        TRIM( Num2LStr ( ParamData%FFRate ) )//'-Hz full-field data ('// &
                        TRIM( Num2LStr( ParamData%FFDTime*( ParamData%NFFSteps - 1 ) ) )//' seconds).'
      ENDIF
  
      CALL WrScr( NewLine//TRIM(TmpErrMsg) )    ! Note: the TmpErrMsg gets used below for the summary file 



      !-------------------------------------------------------------------------------------------------
      ! Write to the summary file
      !-------------------------------------------------------------------------------------------------

   IF ( InitData%SumFileUnit > 0 ) THEN
      WRITE(InitData%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)
      WRITE(InitData%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)    'TurbSim wind type.  Read by InflowWind sub-module '//TRIM(IfW_TSFFWind_Ver%Name)//  &
                                                                                 ' '//TRIM(IfW_TSFFWind_Ver%Ver)
      WRITE(InitData%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)    TRIM(TmpErrMsg)
      WRITE(InitData%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)    '     FileName:                    '//TRIM(InitData%WindFileName)
      WRITE(InitData%SumFileUnit,'(A34,I3)',   IOSTAT=TmpErrStat)    '     Binary file format id:       ',ParamData%WindFileFormat
      WRITE(InitData%SumFileUnit,'(A34,G12.4)',IOSTAT=TmpErrStat)    '     Reference height (m):        ',ParamData%RefHt
      WRITE(InitData%SumFileUnit,'(A34,G12.4)',IOSTAT=TmpErrStat)    '     Timestep (s):                ',ParamData%FFDTime
      WRITE(InitData%SumFileUnit,'(A34,I12)',  IOSTAT=TmpErrStat)    '     Number of timesteps:         ',ParamData%NFFSteps
      WRITE(InitData%SumFileUnit,'(A34,G12.4)',IOSTAT=TmpErrStat)    '     Mean windspeed (m/s):        ',ParamData%MeanFFWS
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

   END SUBROUTINE READ_TurbSim_FF

END SUBROUTINE IfW_TSFFWind_Init
!====================================================================================================


!====================================================================================================
!> This routine acts as a wrapper for the GetWindSpeed routine. It steps through the array of input
!! positions and calls the GetWindSpeed routine to calculate the velocities at each point.
!!
!! There are inefficiencies in how this set of routines is coded, but that is a problem for another
!! day. For now, it merely needs to be functional. It can be fixed up and made all pretty later.
!!
!!   16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
SUBROUTINE IfW_TSFFWind_CalcOutput(Time, PositionXYZ, ParamData,  Velocity, DiskVel, MiscVars, ErrStat, ErrMsg)

   IMPLICIT                                                 NONE


   CHARACTER(*),           PARAMETER                  :: RoutineName="IfW_TSFFWind_CalcOutput"


      ! Passed Variables
   REAL(DbKi),                               INTENT(IN   )  :: Time              !< time from the start of the simulation
   REAL(ReKi),                               INTENT(IN   )  :: PositionXYZ(:,:)  !< Array of XYZ coordinates, 3xN
   TYPE(IfW_TSFFWind_ParameterType),         INTENT(IN   )  :: ParamData         !< Parameters
   REAL(ReKi),                               INTENT(INOUT)  :: Velocity(:,:)     !< Velocity output at Time    (Set to INOUT so that array does not get deallocated)
   REAL(ReKi),                               INTENT(  OUT)  :: DiskVel(3)        !< HACK for AD14: disk velocity output at Time
   TYPE(IfW_TSFFWind_MiscVarType),           INTENT(INOUT)  :: MiscVars          !< Misc variables for optimization (not copied in glue code)

      ! Error handling
   INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat           !< error status
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg            !< The error message

      ! local variables
   INTEGER(IntKi)                                           :: NumPoints      ! Number of points specified by the PositionXYZ array

      ! local counters
   INTEGER(IntKi)                                           :: PointNum       ! a loop counter for the current point

      ! temporary variables
   INTEGER(IntKi)                                           :: TmpErrStat     ! temporary error status
   CHARACTER(ErrMsgLen)                                     :: TmpErrMsg      ! temporary error message



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
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName//" [position=("//   &
                                                      TRIM(Num2LStr(PositionXYZ(1,PointNum)))//", "// &
                                                      TRIM(Num2LStr(PositionXYZ(2,PointNum)))//", "// &
                                                      TRIM(Num2LStr(PositionXYZ(3,PointNum)))//") in wind-file coordinates]" )
         IF (ErrStat >= AbortErrLev) RETURN
      END IF

   ENDDO



      !REMOVE THIS for AeroDyn 15
      ! Return the average disk velocity values needed by AeroDyn 14.  This is the WindInf_ADhack_diskVel routine.
   DiskVel(1)   =  ParamData%MeanFFWS
   DiskVel(2:3) =  0.0_ReKi


   RETURN

END SUBROUTINE IfW_TSFFWind_CalcOutput
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
   FUNCTION FF_Interp(t_in, Position, ParamData, MiscVars, ErrStat, ErrMsg)
   !----------------------------------------------------------------------------------------------------

      IMPLICIT                                              NONE

      CHARACTER(*),           PARAMETER                  :: RoutineName="FF_Interp"

      REAL(DbKi),                         INTENT(IN   )  :: t_in           !< input time
      REAL(ReKi),                         INTENT(IN   )  :: Position(3)    !< takes the place of XGrnd, YGrnd, ZGrnd
      TYPE(IfW_TSFFWind_ParameterType),   INTENT(IN   )  :: ParamData      !< Parameters
      TYPE(IfW_TSFFWind_MiscVarType),     INTENT(INOUT)  :: MiscVars       !< Misc variables for optimization (not copied in glue code)
      REAL(ReKi)                                         :: FF_Interp(3)   !< The U, V, W velocities

      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat        !< error status
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg         !< error message

         ! Local Variables:

      REAL(ReKi)                                         :: TimeShifted
      REAL(ReKi),PARAMETER                               :: Tol = 1.0E-3   ! a tolerance for determining if two reals are the same (for extrapolation)
      REAL(ReKi)                                         :: T
      REAL(ReKi)                                         :: TGRID
      REAL(ReKi)                                         :: Y
      REAL(ReKi)                                         :: YGRID
      REAL(ReKi)                                         :: Z
      REAL(ReKi)                                         :: ZGRID
      REAL(ReKi)                                         :: N(8)           ! array for holding scaling factors for the interpolation algorithm
      REAL(ReKi)                                         :: u(8)           ! array for holding the corner values for the interpolation algorithm across a cubic volume
      REAL(ReKi)                                         :: M(4)           ! array for holding scaling factors for the interpolation algorithm
      REAL(ReKi)                                         :: v(4)           ! array for holding the corner values for the interpolation algorithm across an area

      INTEGER(IntKi)                                     :: IDIM
      INTEGER(IntKi)                                     :: ITHI
      INTEGER(IntKi)                                     :: ITLO
      INTEGER(IntKi)                                     :: IYHI
      INTEGER(IntKi)                                     :: IYLO
      INTEGER(IntKi)                                     :: IZHI
      INTEGER(IntKi)                                     :: IZLO

      LOGICAL                                            :: OnGrid

      !-------------------------------------------------------------------------------------------------
      ! Initialize variables
      !-------------------------------------------------------------------------------------------------

      FF_Interp(:)        = 0.0_ReKi                         ! the output velocities (in case ParamData%NFFComp /= 3)

      ErrStat              = ErrID_None
      ErrMsg               = ""
      
      !-------------------------------------------------------------------------------------------------
      ! Find the bounding time slices.
      !-------------------------------------------------------------------------------------------------

      ! Perform the time shift.  At t_in=0, a point half the grid width downstream (ParamData%FFYHWid) will index into the zero time slice.
      ! If we did not do this, any point downstream of the tower at the beginning of the run would index outside of the array.
      ! This all assumes the grid width is at least as large as the rotor.  If it isn't, then the interpolation will not work.


      TimeShifted = t_in + ( ParamData%InitXPosition - Position(1) )*ParamData%InvMFFWS    ! in distance, X: InputInfo%Position(1) - ParamData%InitXPosition - t*ParamData%MeanFFWS


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
               ErrMsg   = ' Error: FF wind array was exhausted at '//TRIM( Num2LStr( REAL( t_in,   ReKi ) ) )// &
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
!>  This subroutine cleans up any data that is still allocated.  The (possibly) open files are
!!  closed in InflowWindMod.
!!
!!  16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
SUBROUTINE IfW_TSFFWind_End( ParamData, MiscVars, ErrStat, ErrMsg)

   IMPLICIT                                                 NONE

   CHARACTER(*),           PARAMETER                     :: RoutineName="IfW_TSFFWind_End"



      ! Passed Variables
   TYPE(IfW_TSFFWind_ParameterType),      INTENT(INOUT)  :: ParamData         !< Parameters
   TYPE(IfW_TSFFWind_MiscVarType),        INTENT(INOUT)  :: MiscVars          !< Misc variables for optimization (not copied in glue code)


      ! Error Handling
   INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat           !< determines if an error has been encountered
   CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg            !< Message about errors


      ! Local Variables
   INTEGER(IntKi)                                        :: TmpErrStat     ! temporary error status
   CHARACTER(ErrMsgLen)                                  :: TmpErrMsg      ! temporary error message


      !-=- Initialize the routine -=-

   ErrMsg   = ''
   ErrStat  = ErrID_None



      ! Destroy parameter data

   CALL IfW_TSFFWind_DestroyParam(       ParamData,     TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)


      ! Destroy the misc data

   CALL IfW_TSFFWind_DestroyMisc(  MiscVars,   TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)



END SUBROUTINE IfW_TSFFWind_End

!====================================================================================================
END MODULE IfW_TSFFWind
