!>  This module uses full-field binary wind files to determine the wind inflow.
!!  This module assumes that the origin, (0,0,0), is located at the tower centerline at ground level,
!!  and that all units are specified in the metric system (using meters and seconds).
!!  Data is shifted by half the grid width to account for turbine yaw (so that data in the X
!!  direction actually starts at -1*ParamData%FF%FFYHWid meters).
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
   USE                                          IfW_FFWind_Base

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
SUBROUTINE IfW_TSFFWind_Init(InitData, ParamData, MiscVars, InitOutData, ErrStat, ErrMsg)

   IMPLICIT                       NONE

   CHARACTER(*),           PARAMETER                        :: RoutineName="IfW_TSFFWind_Init"


      ! Passed Variables
   TYPE(IfW_TSFFWind_InitInputType),         INTENT(IN   )  :: InitData          !< Initialization data passed to the module
   TYPE(IfW_TSFFWind_ParameterType),         INTENT(  OUT)  :: ParamData         !< Parameters
   TYPE(IfW_TSFFWind_MiscVarType),           INTENT(  OUT)  :: MiscVars          !< Misc variables for optimization (not copied in glue code)
   TYPE(IfW_TSFFWind_InitOutputType),        INTENT(  OUT)  :: InitOutData       !< Initial output


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

   ParamData%FF%InterpTower = .false.
   ParamData%FF%AddMeanAfterInterp = .false.


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
   ParamData%FF%WindFileFormat = Dum_Int2

   SELECT CASE (ParamData%FF%WindFileFormat)

      CASE ( 7, 8 )                                                    ! TurbSim binary format

         CALL Read_TurbSim_FF(UnitWind, TmpErrStat, TmpErrMsg)
            CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )                  
            IF ( ErrStat >= AbortErrLev ) THEN
               CLOSE ( UnitWind )
               RETURN
            END IF
         

      CASE DEFAULT
         CALL SetErrStat( ErrID_Fatal, ' Error: This is not a TurbSim binary wind file type (binary format identifier: '//  &
                  TRIM(Num2LStr(ParamData%FF%WindFileFormat))//'.  This might be a Bladed style binary wind file.', &
                  ErrStat, ErrMsg, RoutineName )
         RETURN

   END SELECT


   IF (ParamData%FF%Periodic) THEN
      ParamData%FF%InitXPosition = 0                ! start at the hub
      ParamData%FF%TotalTime     = ParamData%FF%NFFSteps*ParamData%FF%FFDTime
   ELSE
      ParamData%FF%InitXPosition = ParamData%FF%FFYHWid          ! start half the grid width ahead of the turbine
      ParamData%FF%TotalTime     = (ParamData%FF%NFFSteps-1)*ParamData%FF%FFDTime
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


      ParamData%FF%NFFComp = 3                                              ! this file contains 3 wind components
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
            ParamData%FF%Periodic = Dum_Int2 == INT( 8, B2Ki) ! the number 7 is used for non-periodic wind files; 8 is periodic wind


            ! Read in the 4-byte integer. Can't use library read routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4             ! the number of grid points vertically, INT(4)
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading the number of z grid points in the FF binary file "'//TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName )                  
               RETURN
            ENDIF
            ParamData%FF%NZGrids = Dum_Int4


            ! Read in the 4-byte integer. Can't use library read routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4             ! the number of grid points laterally, INT(4)
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading the number of y grid points in the FF binary file "'//TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName )                  
               RETURN
            ENDIF
            ParamData%FF%NYGrids = Dum_Int4


            ! Read in the 4-byte integer. Can't use library read routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4             ! the number of tower points, INT(4)
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading the number of tower points in the FF binary file "'//TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName )                  
               RETURN
            ENDIF
            ParamData%FF%NTGrids = Dum_Int4


            ! Read in the 4-byte integer. Can't use library read routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4             ! the number of time steps, INT(4)
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading the number of time steps in the FF binary file "'//TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName )                  
               RETURN
            ENDIF
            ParamData%FF%NFFSteps = Dum_Int4


            ! Read in the 4-byte real. Can't use library read routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4            ! grid spacing in vertical direction (dz), REAL(4), in m
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading dz in the FF binary file "'//TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName )                  
               RETURN
            ENDIF
            ParamData%FF%InvFFZD = 1.0/Dum_Real4                            ! 1/dz
            ParamData%FF%FFZHWid = 0.5*(ParamData%FF%NZGrids-1)*Dum_Real4                ! half the grid height


            ! Read in the 4-byte real. Can't use library read routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4            ! grid spacing in lateral direction (dy), REAL(4), in m
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading dy in the FF binary file "'//TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName )                  
               RETURN
            ENDIF
            ParamData%FF%InvFFYD = 1.0 / Dum_Real4                          ! 1/dy
            ParamData%FF%FFYHWid = 0.5*(ParamData%FF%NYGrids-1)*Dum_Real4                ! half grid grid width


            ! Read in the 4-byte real. Can't use library read routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4            ! grid spacing in time (dt), REAL(4), in m/s
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading dt in the FF binary file "'//TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName )                  
               RETURN
            ENDIF
            ParamData%FF%FFDTime = Dum_Real4
            ParamData%FF%FFRate  = 1.0/ParamData%FF%FFDTime


            ! Read in the 4-byte real. Can't use library read routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4            ! the mean wind speed at hub height, REAL(4), in m/s
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading mean wind speed in the FF binary file "'//TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName )                  
               RETURN
            ENDIF
            ParamData%FF%MeanFFWS = Dum_Real4
            ParamData%FF%InvMFFWS = 1.0 / ParamData%FF%MeanFFWS


            ! Read in the 4-byte real. Can't use library read routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4            ! height of the hub, REAL(4), in m
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading zHub in the FF binary file "'//TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName )                  
               RETURN
            ENDIF
            ParamData%FF%RefHt = Dum_Real4


            ! Read in the 4-byte real. Can't use library read routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4            ! height of the bottom of the grid, REAL(4), in m
            IF ( TmpErrStat /= 0 )  THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading GridBase in the FF binary file "'//TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName )                  
               RETURN
            ENDIF
            ParamData%FF%GridBase = Dum_Real4

    !        ZGOffset = ParamData%FF%RefHt - ParamData%FF%GridBase  - ParamData%FF%FFZHWid


         !----------------------------------------------------------------------------------------------
         ! Read the binary scaling factors
         !----------------------------------------------------------------------------------------------

            DO IC = 1,ParamData%FF%NFFComp
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
         
         CALL WrScr( NewLine//'   Reading a '//TRIM( Num2LStr(ParamData%FF%NYGrids) )//'x'//TRIM( Num2LStr(ParamData%FF%NZGrids) )//  &
                    ' grid ('//TRIM( Num2LStr(ParamData%FF%FFYHWid*2) )//' m wide, '// &
                    TRIM( Num2LStr(ParamData%FF%GridBase) )//' m to '// &
                    TRIM( Num2LStr(ParamData%FF%GridBase+ParamData%FF%FFZHWid*2) )//&
                    ' m above ground) with a characteristic wind speed of '// &
                    TRIM( Num2LStr(ParamData%FF%MeanFFWS) )//' m/s. '//TRIM(DescStr) )


      !----------------------------------------------------------------------------------------------
      ! Allocate arrays for the FF grid as well as the tower points, if they exist
      !----------------------------------------------------------------------------------------------

         IF ( .NOT. ALLOCATED( ParamData%FF%FFData ) ) THEN
            CALL AllocAry( ParamData%FF%FFData, ParamData%FF%NZGrids, ParamData%FF%NYGrids, ParamData%FF%NFFComp, ParamData%FF%NFFSteps, &
                  'Full-field wind data array.', TmpErrStat, TmpErrMsg )
            CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName ) 
            IF ( ErrStat >= AbortErrLev ) RETURN
         ENDIF


         IF ( ParamData%FF%NTGrids > 0 ) THEN

            IF ( .NOT. ALLOCATED( ParamData%FF%FFTower ) ) THEN
               CALL AllocAry( ParamData%FF%FFTower, ParamData%FF%NFFComp, ParamData%FF%NTGrids, ParamData%FF%NFFSteps, &
                     'Tower wind file data array.', TmpErrStat, TmpErrMsg )
               CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName ) 
               IF ( ErrStat >= AbortErrLev ) RETURN
            ENDIF

         ENDIF

      !-------------------------------------------------------------------------------------------------
      ! Read the 16-bit data and scale it to 32-bit reals
      !-------------------------------------------------------------------------------------------------

         ! Loop through time.

         DO IT=1,ParamData%FF%NFFSteps

            !...........................................................................................
            ! Read grid data at this time step.
            !...........................................................................................

            DO IZ=1,ParamData%FF%NZGrids
               ! Zgrid(IZ) = Z1 + (IZ-1)*dz                 ! Vertical location of grid data point, in m relative to ground

               DO IY=1,ParamData%FF%NYGrids
                  ! Ygrid(IY) = -0.5*(ny-1)*dy + (IY-1)*dy  ! Horizontal location of grid data point, in m relative to tower centerline

                  DO IC=1,ParamData%FF%NFFComp                           ! number of wind components (U, V, W)

                        ! Read in the 2-byte integer. Can't use library read routines for this.
                     READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2       ! normalized wind-component, INT(2)
                     IF ( TmpErrStat /= 0 )  THEN                        
                        CALL SetErrStat( ErrID_Fatal, ' Error reading grid wind components in the FF binary file "'// &
                                    TRIM( InitData%WindFileName )//'."', ErrStat, ErrMsg, RoutineName ) 
                        RETURN
                     ENDIF

                     ParamData%FF%FFData(IZ,IY,IC,IT) = ( Dum_Int2 - Voffset(IC) ) / VSlope(IC)

                  ENDDO !IC

               ENDDO !IY

            ENDDO ! IZ


            !...........................................................................................
            ! Read the tower data at this time step.
            !...........................................................................................

            DO IZ=1,ParamData%FF%NTGrids         ! If NTGrids<1, there are no tower points & FFTower is not allocated

               ! Ytower     = 0               ! Lateral location of the tower data point, in m relative to tower centerline
               ! Ztower(IZ) = Z1 - (IZ-1)*dz  ! Vertical location of tower data point, in m relative to ground

               DO IC=1,ParamData%FF%NFFComp   ! number of wind components

                     ! Read in a 2-byte integer. Can't use library routines for this.
                  READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2       ! normalized wind-component, INT(2)
                  IF ( TmpErrStat /= 0 )  THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading tower wind components in the FF binary file "'//TRIM(InitData%WindFileName)//'."'&
                                      , ErrStat, ErrMsg, RoutineName )                    
                     RETURN
                  ENDIF

                  ParamData%FF%FFTower(IC,IZ,IT) = ( Dum_Int2 - Voffset(IC) ) / VSlope(IC)  ! wind-component scaled to m/s

               ENDDO !IC

            ENDDO ! IZ


         ENDDO ! IT

      !-------------------------------------------------------------------------------------------------
      ! close the file and return
      !-------------------------------------------------------------------------------------------------

      CLOSE ( UnitWind )

      IF ( ParamData%FF%Periodic ) THEN
         TmpErrMsg   = '   Processed '//TRIM( Num2LStr( ParamData%FF%NFFSteps ) )//' time steps of '// &
                        TRIM( Num2LStr ( ParamData%FF%FFRate ) )//'-Hz full-field data (period of '// &
                        TRIM( Num2LStr( ParamData%FF%FFDTime*( ParamData%FF%NFFSteps ) ) )//' seconds).'
      ELSE
         TmpErrMsg   = '   Processed '//TRIM( Num2LStr( ParamData%FF%NFFSteps ) )//' time steps of '// &
                        TRIM( Num2LStr ( ParamData%FF%FFRate ) )//'-Hz full-field data ('// &
                        TRIM( Num2LStr( ParamData%FF%FFDTime*( ParamData%FF%NFFSteps - 1 ) ) )//' seconds).'
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
      WRITE(InitData%SumFileUnit,'(A34,I3)',   IOSTAT=TmpErrStat)    '     Binary file format id:       ',ParamData%FF%WindFileFormat
      WRITE(InitData%SumFileUnit,'(A34,G12.4)',IOSTAT=TmpErrStat)    '     Reference height (m):        ',ParamData%FF%RefHt
      WRITE(InitData%SumFileUnit,'(A34,G12.4)',IOSTAT=TmpErrStat)    '     Timestep (s):                ',ParamData%FF%FFDTime
      WRITE(InitData%SumFileUnit,'(A34,I12)',  IOSTAT=TmpErrStat)    '     Number of timesteps:         ',ParamData%FF%NFFSteps
      WRITE(InitData%SumFileUnit,'(A34,G12.4)',IOSTAT=TmpErrStat)    '     Mean windspeed (m/s):        ',ParamData%FF%MeanFFWS
      WRITE(InitData%SumFileUnit,'(A34,L1)',   IOSTAT=TmpErrStat)    '     Windfile is periodic:        ',ParamData%FF%Periodic
      WRITE(InitData%SumFileUnit,'(A34,L1)',   IOSTAT=TmpErrStat)    '     Windfile includes tower:     ',ParamData%FF%NTGrids > 0

      IF ( ParamData%FF%Periodic ) THEN
         WRITE(InitData%SumFileUnit,'(A)',     IOSTAT=TmpErrStat)    '     Time range (s):              [ '// &
                     TRIM(Num2LStr(0.0_ReKi))//' : '//TRIM(Num2LStr(ParamData%FF%TotalTime))//' ]'
      ELSE  ! Shift the time range to compensate for the shifting of the wind grid
         WRITE(InitData%SumFileUnit,'(A)',     IOSTAT=TmpErrStat)    '     Time range (s):              [ '// &
                     TRIM(Num2LStr(-ParamData%FF%InitXPosition*ParamData%FF%InvMFFWS))//' : '// &
                     TRIM(Num2LStr(ParamData%FF%TotalTime-ParamData%FF%InitXPosition*ParamData%FF%InvMFFWS))//' ]'
      ENDIF

      WRITE(InitData%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)    '     Y range (m):                 [ '// &
                     TRIM(Num2LStr(-ParamData%FF%FFYHWid))//' : '//TRIM(Num2LStr(ParamData%FF%FFYHWid))//' ]'

      IF ( ParamData%FF%NTGrids > 0 ) THEN
         WRITE(InitData%SumFileUnit,'(A)',     IOSTAT=TmpErrStat)    '     Z range (m):                 [ '// &
                     TRIM(Num2LStr(0.0_ReKi))//' : '//TRIM(Num2LStr(ParamData%FF%RefHt + ParamData%FF%FFZHWid))//' ]'
      ELSE
         WRITE(InitData%SumFileUnit,'(A)',     IOSTAT=TmpErrStat)    '     Z range (m):                 [ '// &
                     TRIM(Num2LStr(ParamData%FF%RefHt - ParamData%FF%FFZHWid))//' : '//TRIM(Num2LStr(ParamData%FF%RefHt + ParamData%FF%FFZHWid))//' ]'
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
!> This routine computes the wind speed at each of the PositionXYZ points.
SUBROUTINE IfW_TSFFWind_CalcOutput(Time, PositionXYZ, p,  Velocity, DiskVel, MiscVars, ErrStat, ErrMsg)

   IMPLICIT                                                 NONE


   CHARACTER(*),           PARAMETER                  :: RoutineName="IfW_TSFFWind_CalcOutput"


      ! Passed Variables
   REAL(DbKi),                               INTENT(IN   )  :: Time              !< time from the start of the simulation
   REAL(ReKi),                               INTENT(IN   )  :: PositionXYZ(:,:)  !< Array of XYZ coordinates, 3xN
   TYPE(IfW_TSFFWind_ParameterType),         INTENT(IN   )  :: p                 !< Parameters
   REAL(ReKi),                               INTENT(INOUT)  :: Velocity(:,:)     !< Velocity output at Time    (Set to INOUT so that array does not get deallocated)
   REAL(ReKi),                               INTENT(  OUT)  :: DiskVel(3)        !< HACK for AD14: disk velocity output at Time
   TYPE(IfW_TSFFWind_MiscVarType),           INTENT(INOUT)  :: MiscVars          !< Misc variables for optimization (not copied in glue code)

      ! Error handling
   INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat           !< error status
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg            !< The error message



   CALL IfW_FFWind_CalcOutput(Time, PositionXYZ, p%FF, Velocity, DiskVel, ErrStat, ErrMsg)


   RETURN

END SUBROUTINE IfW_TSFFWind_CalcOutput
!====================================================================================================
!>  This subroutine cleans up any data that is still allocated.  The (possibly) open files are
!!  closed in InflowWindMod.
!!
!!  16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
SUBROUTINE IfW_TSFFWind_End( p, m, ErrStat, ErrMsg)

   IMPLICIT                                                 NONE

   CHARACTER(*),           PARAMETER                     :: RoutineName="IfW_TSFFWind_End"



      ! Passed Variables
   TYPE(IfW_TSFFWind_ParameterType),      INTENT(INOUT)  :: p             !< Parameters
   TYPE(IfW_TSFFWind_MiscVarType),        INTENT(INOUT)  :: m             !< Misc variables for optimization (not copied in glue code)


      ! Error Handling
   INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat       !< determines if an error has been encountered
   CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg        !< Message about errors


      ! Local Variables
   INTEGER(IntKi)                                        :: TmpErrStat     ! temporary error status
   CHARACTER(ErrMsgLen)                                  :: TmpErrMsg      ! temporary error message


      !-=- Initialize the routine -=-

   ErrMsg   = ''
   ErrStat  = ErrID_None



      ! Destroy parameter data

   CALL IfW_TSFFWind_DestroyParam( p, TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)


      ! Destroy the misc data

   CALL IfW_TSFFWind_DestroyMisc( m, TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)



END SUBROUTINE IfW_TSFFWind_End

!====================================================================================================
END MODULE IfW_TSFFWind
