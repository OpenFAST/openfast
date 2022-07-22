MODULE FDWind
! This module reads and processes 4-dimensional wind fields.
! The subroutines were originally created by Marshall Buhl to read LES data provided by researchers
! at NCAR. It was later updated by Bonnie Jonkman to read DNS data provided by researchers at CoRA.
!
! Data are assumed to be in units of meters and seconds.
!
!  7 Oct 2009    B. Jonkman, NREL/NWTC using subroutines from AeroDyn 12.57
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
! File last committed: $Date: 2014-07-29 13:30:04 -0600 (Tue, 29 Jul 2014) $
! (File) Revision #: $Rev: 125 $
! URL: $HeadURL$
!**********************************************************************************************************************************

   USE                     NWTC_Library
   USE                     SharedInflowDefs
   USE                     InflowWind_Module_Types

   IMPLICIT                NONE
   PRIVATE

      ! FD_Wind

   REAL(ReKi)                   :: DelXgrid                                   ! The nondimensional distance between grid points in the x direction.
   REAL(ReKi)                   :: DelYgrid                                   ! The nondimensional distance between grid points in the y direction.
   REAL(ReKi)                   :: DelZgrid                                   ! The nondimensional distance between grid points in the z direction.
   REAL(ReKi)                   :: FDper                                      ! Total time in dataset.
   REAL(ReKi)                   :: FDTime   (2)                               ! Times for the 4D wind files.
   REAL(ReKi), ALLOCATABLE      :: FDu      (:,:,:,:)                         ! The u-component array of 4D wind data.
   REAL(ReKi), ALLOCATABLE      :: FDv      (:,:,:,:)                         ! The v-component array of 4D wind data.
   REAL(ReKi), ALLOCATABLE      :: FDw      (:,:,:,:)                         ! The w-component array of 4D wind data.
   REAL(ReKi), ALLOCATABLE      :: FDuData  (:,:,:,:)                         ! The u-component array of all 4D wind data when used with advection.
   REAL(ReKi), ALLOCATABLE      :: FDvData  (:,:,:,:)                         ! The v-component array of all 4D wind data when used with advection.
   REAL(ReKi), ALLOCATABLE      :: FDwData  (:,:,:,:)                         ! The w-component array of all 4D wind data when used with advection.
   REAL(ReKi)                   :: Lx                                         ! Fractional location of tower centerline from upwind end to downwind end of the dataset.
   REAL(ReKi)                   :: Ly                                         ! Fractional location of tower centerline from right (looking downwind) to left side of the dataset.
   REAL(ReKi)                   :: Lz                                         ! Fractional location of hub height from bottom to top of dataset.
   REAL(ReKi)                   :: Offsets  (3)                               ! Offsets to convert integer data to actual wind speeds.
!FIXME: move to otherstates
   REAL(ReKi), SAVE             :: PrevTime                                   ! The previous time this was called -- so we can go back in time if necessary
   REAL(ReKi)                   :: RotDiam                                    ! Rotor diameter.
   REAL(ReKi)                   :: ScalFact (3)                               ! Scaling factors to convert integer data to actual wind speeds.
   REAL(ReKi)                   :: ScaleVel                                   ! Scaling velocity, U0.  2*U0 is the difference in wind speed between the top and bottom of the wave.
   REAL(ReKi), ALLOCATABLE      :: Times4D  (:)                               ! The list of times for the 4D-wind input files.
   REAL(ReKi)                   :: Tm_max                                     ! The total nondimensional time of the dataset.
   REAL(ReKi)                   :: TSclFact                                   ! Scale factor for time (h/U0).
   REAL(ReKi)                   :: T_4D_En                                    ! Time at which the wave event ends.
   REAL(ReKi)                   :: T_4D_St                                    ! Time at which the wave event starts.
   REAL(ReKi)                   :: Xmax                                       ! The dimensional downwind length of the dataset.
   REAL(ReKi)                   :: Xt                                         ! Distance of the tower from the upwind end of the dataset.
   REAL(ReKi)                   :: Ymax                                       ! The dimensional lateral width of the dataset.
   REAL(ReKi)                   :: Yt                                         ! Distance of the tower from the right side of the dataset (looking downwind).
   REAL(ReKi)                   :: Zmax                                       ! The dimensional vertical height of the dataset.
   REAL(ReKi)                   :: Zt                                         ! Distance of the hub from the bottom of the dataset.
   REAL(ReKi)                   :: Zref                                       ! The reference height (hub height)

   INTEGER                      :: FD_DF_X                                    ! The decimation factor for the 4D wind data in the x direction.
   INTEGER                      :: FD_DF_Y                                    ! The decimation factor for the 4D wind data in the y direction.
   INTEGER                      :: FD_DF_Z                                    ! The decimation factor for the 4D wind data in the z direction.
   INTEGER                      :: FDFileNo                                   ! The 4D wind file number.
   INTEGER                      :: FDRecL                                     ! The length, in bytes, of the LE binary records.
   INTEGER                      :: Ind4DAdv                                   ! Index of the file to be used in advection
   INTEGER                      :: Ind4Dnew                                   ! Index of the newest 4D wind file.
   INTEGER                      :: Ind4Dold                                   ! Index of the older 4D wind file.
   INTEGER                      :: Num4Dt                                     ! The number of 4D wind grids, one grid per time step.
   INTEGER, PARAMETER           :: Num4DtD = 2                                ! The number of 4D wind grids stored in memory, normally 2
   INTEGER                      :: Num4Dx                                     ! The number of 4D wind grid points in the x direction.
   INTEGER                      :: Num4DxD                                    ! The decimated number of 4D wind grid points in the x direction.
   INTEGER                      :: Num4DxD1                                   ! The decimated number of 4D wind grid points in the x direction minus 1.
   INTEGER                      :: Num4Dy                                     ! The number of 4D wind grid points in the y direction.
   INTEGER                      :: Num4DyD                                    ! The decimated number of 4D wind grid points in the y direction.
   INTEGER                      :: Num4DyD1                                   ! The decimated number of 4D wind grid points in the y direction minus 1.
   INTEGER                      :: Num4Dz                                     ! The number of 4D wind grid points in the z direction.
   INTEGER                      :: Num4DzD                                    ! The decimated number of 4D wind grid points in the z direction.
   INTEGER                      :: Num4DzD1                                   ! The decimated number of 4D wind grid points in the z direction minus 1.
   INTEGER                      :: NumAdvect                                  ! Number of frozen timesteps to advect past the turbine
   INTEGER                      :: Shft4Dnew                                  ! Number of times the x-data needs to be shifted for advection
   INTEGER, ALLOCATABLE         :: Times4DIx (:)                              ! Index number of the 4D time files (used for advection)

   INTEGER                      :: FDUnit                                     ! Unit number for reading wind files

   LOGICAL                      :: Advect                                     ! Flag to indicate whether or not to advect a given data set or to just use the time step files
   LOGICAL                      :: VertShft                                   ! Flag to indicate whether or not to shift the z values for the w component.

!FIXME: move to parameters or otherstates
   LOGICAL, SAVE                :: Initialized = .FALSE.

   CHARACTER(5), ALLOCATABLE    :: AdvFiles (:)
   CHARACTER(1024)              :: FDSpath                                    ! The path to the 4D wind files.


   PUBLIC                       :: FD_Init
   PUBLIC                       :: FD_GetWindSpeed
   PUBLIC                       :: FD_Terminate
   PUBLIC                       :: FD_GetValue


CONTAINS
!====================================================================================================
SUBROUTINE FD_Init(UnWind, WindFile, RefHt, ErrStat)
!  This subroutine is called at the beginning of a simulation to initialize the module.
!----------------------------------------------------------------------------------------------------

      ! Passed variables

   INTEGER,         INTENT(IN)    :: UnWind                       ! unit number for reading wind files
   CHARACTER(*),    INTENT(IN)    :: WindFile                     ! Name of the 4D wind parameter file (.fdp)
   REAL(ReKi),      INTENT(IN)    :: RefHt                        ! The reference height for the billow (should be hub height)
   INTEGER,         INTENT(OUT)   :: ErrStat                      ! return 0 if no errors; non-zero otherwise

      ! Local variables

   CHARACTER(1024)                :: FDTSfile                     ! name of the 4D time step file
   REAL(ReKi)                     :: FDTimStp                     ! Average time step for 4D wind data.
   INTEGER                        :: IT

   !-------------------------------------------------------------------------------------------------
   ! Check that the module hasn't already been initialized.
   !-------------------------------------------------------------------------------------------------

   IF ( Initialized ) THEN
      CALL WrScr( ' FDWind has already been initialized.' )
      ErrStat = 1
      RETURN
   ELSE
      ErrStat = 0
!      CALL NWTC_Init()    ! Initialized in IfW_Init
   END IF


   !-------------------------------------------------------------------------------------------------
   ! Set the reference height for the wind file (this takes the place of HH that was used earlier)
   !-------------------------------------------------------------------------------------------------

   ZRef = RefHt

   !-------------------------------------------------------------------------------------------------
   ! Read the main 4D input file
   !-------------------------------------------------------------------------------------------------

   CALL ReadFDP( UnWind, WindFile, FDTSfile, ErrStat )
   IF ( ErrStat /= 0 ) RETURN

   !-------------------------------------------------------------------------------------------------
   ! Get the times array, which must be scaled and shifted later using TSclFact and T_4D_St
   !-------------------------------------------------------------------------------------------------

   CALL Read4Dtimes ( UnWind, FDTSfile, ErrStat )
   IF ( ErrStat /= 0 ) RETURN


   !-------------------------------------------------------------------------------------------------
   ! Calculate some values that don't change during the run.
   !-------------------------------------------------------------------------------------------------

   FDRecL      = 2*Num4Dx*Num4Dy                                           ! The length, in bytes, of the 4D binary records.
   Num4DxD     = ( Num4Dx + FD_DF_X - 1 )/FD_DF_X                          ! The decimated number of 4D wind grid points in the x direction.
   Num4DyD     = ( Num4Dy + FD_DF_Y - 1 )/FD_DF_Y                          ! The decimated number of 4D wind grid points in the y direction.
   Num4DzD     = ( Num4Dz + FD_DF_Z - 1 )/FD_DF_Z                          ! The decimated number of 4D wind grid points in the z direction.
   Num4DxD1    = Num4DxD - 1                                               ! The decimated number of 4D wind grid points in the x direction minus 1.
   Num4DyD1    = Num4DyD - 1                                               ! The decimated number of 4D wind grid points in the y direction minus 1.
   Num4DzD1    = Num4DzD - 1                                               ! The decimated number of 4D wind grid points in the z direction minus 1.

   Tm_max      = Times4D(Num4Dt)                                           ! Time of end of dataset.
   IF ( ADVECT ) THEN
      FDTimStp   = Xmax / ( ( Num4Dx - 1 )*( ScaleVel )*Num4Dt )           ! The timestep is calculated by the approximation dx/dt ~= U0 (divide by num4dt to get delta for a full timestep).
      FDper      = FDTimStp * Num4Dt                                       ! Total time in dataset. (We have periodic time, so multiply by number of time steps, without subtracting 1)
      TSclFact   = FDper / Tm_max                                          ! Equivalent scale factor for time.
   ELSE
      FDper       = TSclFact*Tm_max                                        ! Total time in dataset.
      FDTimStp    = FDper/( Num4Dt - 1 )                                   ! Average time step.
   ENDIF

   T_4D_En     = T_4D_St + FDper                                           ! Time for the end of the dataset.
   Xt          = Xmax*Lx                                                   ! Distance of the tower from the upwind end of the dataset.
   Yt          = Ymax*Ly                                                   ! Distance of the tower from the right side of the dataset (looking downwind).
   Zt          = Zmax*Lz                                                   ! Distance of the hub from the bottom of the dataset.
   DelXgrid    = 1.0/Num4DxD1                                              ! The nondimensional distance between grid points in the x direction.
   DelYgrid    = 1.0/Num4DyD1                                              ! The nondimensional distance between grid points in the y direction.
   DelZgrid    = 1.0/Num4DzD1                                              ! The nondimensional distance between grid points in the z direction.


   !-------------------------------------------------------------------------------------------------
   ! Scale and shift the times array using TSclFact and T_4D_St
   !-------------------------------------------------------------------------------------------------

   DO IT=1,Num4Dt

      Times4D(IT) = TSclFact*Times4D(IT) + T_4D_St

   ENDDO ! IT


   !-------------------------------------------------------------------------------------------------
   ! Allocate velocity arrays and fill Data arrays for advection (DNS files)
   !-------------------------------------------------------------------------------------------------

   IF (.NOT. ALLOCATED(FDu) ) THEN
!      CALL AllocAry ( FDu, Num4DxD, Num4DyD, Num4DzD, 2, 'U-component velocity array (FDu)', ErrStat)
      ALLOCATE ( FDu(Num4DxD,Num4DyD,Num4DzD,2), STAT=ErrStat )

      IF ( ErrStat /= 0 )  THEN
         CALL WrScr ( ' Error allocating memory for the U-component velocity array (FDu) array.' )
         RETURN
      END IF
   END IF

   IF (.NOT. ALLOCATED(FDv) ) THEN
!      CALL AllocAry ( FDv, Num4DxD, Num4DyD, Num4DzD, 2, 'V-component velocity array (FDv)', ErrStat)
      ALLOCATE ( FDv(Num4DxD,Num4DyD,Num4DzD,2), STAT=ErrStat )

      IF ( ErrStat /= 0 )  THEN
         CALL WrScr ( ' Error allocating memory for the V-component velocity array (FDv) array.' )
         RETURN
      END IF
   END IF

   IF (.NOT. ALLOCATED(FDw) ) THEN
!      CALL AllocAry ( FDw, Num4DxD, Num4DyD, Num4DzD, 2, 'W-component velocity array (FDw)', ErrStat)
      ALLOCATE ( FDw(Num4DxD,Num4DyD,Num4DzD,2), STAT=ErrStat )

      IF ( ErrStat /= 0 )  THEN
         CALL WrScr ( ' Error allocating memory for the W-component velocity array (FDw) array.' )
         RETURN
      END IF
   END IF

   IF ( ADVECT ) THEN

      IF (.NOT. ALLOCATED(FDuData) ) THEN
!         CALL AllocAry ( FDuData, Num4DxD, Num4DyD, Num4DzD, Num4Dt, 'U-component velocity array (FDuData)', ErrStat)
         ALLOCATE ( FDuData(Num4DxD,Num4DyD,Num4DzD,Num4Dt), STAT=ErrStat )

         IF ( ErrStat /= 0 )  THEN
            CALL WrScr ( ' Error allocating memory for the U-component velocity array (FDuData) array.' )
            RETURN
         END IF
      END IF

      IF (.NOT. ALLOCATED(FDvData) ) THEN
!         CALL AllocAry ( FDvData, Num4DxD, Num4DyD, Num4DzD, Num4Dt, 'V-component velocity array (FDvData)', ErrStat)
         ALLOCATE ( FDvData(Num4DxD,Num4DyD,Num4DzD,Num4Dt), STAT=ErrStat )

         IF ( ErrStat /= 0 )  THEN
            CALL WrScr ( ' Error allocating memory for the V-component velocity array (FDvData) array.' )
            RETURN
         END IF
      END IF

      IF (.NOT. ALLOCATED(FDwData) ) THEN
!         CALL AllocAry ( FDwData, Num4DxD, Num4DyD, Num4DzD, Num4Dt, 'W-component velocity array (FDwData)', ErrStat)
         ALLOCATE ( FDwData(Num4DxD,Num4DyD,Num4DzD,Num4Dt), STAT=ErrStat )

         IF ( ErrStat /= 0 )  THEN
            CALL WrScr ( ' Error allocating memory for the W-component velocity array (FDwData) array.' )
            RETURN
         END IF
      END IF

      CALL ReadAll4DData(UnWind, ErrStat) !This needs AdvFiles(:), which was is read in ReadFDP()
      IF ( ErrStat /= 0 ) RETURN

   ENDIF


   !-------------------------------------------------------------------------------------------------
   ! Determine the first file needed for this simulation.
   !-------------------------------------------------------------------------------------------------
   Ind4Dold  = 1                                                           ! Put the old stuff in the first part of the array.
   Ind4Dnew  = 2                                                           ! Put the new stuff in the second part of the array.

   Shft4Dnew = 0


   IF ( T_4D_St >= 0.0 )  THEN
      FDFileNo = 1
   ELSE
      FDFileNo = Num4Dt
      DO IT=1,Num4Dt
         IF ( Times4D(IT) > 0.0 )  THEN
            FDFileNo = IT - 1
            EXIT
         END IF
      END DO ! IT
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Open, read, and close the first set of files.
   !-------------------------------------------------------------------------------------------------
   FDTime(Ind4Dold) = Times4D(FDFileNo)                                 ! Set the time for this file.

   IF ( ADVECT ) THEN
      CALL Load4DData(Ind4Dold)     ! load data stored in FDuData, FDvData, and FDwData arrays
   ELSE
      CALL LoadLESData( UnWind, FDFileNo, Ind4Dold, ErrStat )
   END IF


   !-------------------------------------------------------------------------------------------------
   ! Open, read, and close the second set of files.
   !-------------------------------------------------------------------------------------------------
   FDFileNo  = FDFileNo + 1


   IF ( ADVECT ) THEN
      FDFileNo = MOD(FDFileNo-1,Num4Dt) + 1

      IF (FDFileNo == 1) THEN
         Shft4Dnew = Shft4Dnew + 1

         IF (Ind4DAdv <= NumAdvect) THEN                             ! Ind4DAdv was set in ReadFDP
            IF ( MOD( Shft4Dnew, Num4Dx ) == 0 ) THEN
               CALL ReadAll4DData(UnWind, ErrStat)
               IF ( ErrStat /= 0 ) RETURN
            END IF
         END IF

      ENDIF

      FDTime(Ind4Dnew) = Times4D(FDFileNo) + Shft4Dnew*FDPer         ! Set the time for this file.

      CALL Load4DData( Ind4Dnew )    ! shift the data

   ELSE
      FDTime(Ind4Dnew) = Times4D(FDFileNo)                                           ! Set the time for this file.

      CALL LoadLESData( UnWind, FDFileNo, Ind4Dnew, ErrStat )
   ENDIF


   !-------------------------------------------------------------------------------------------------
   ! Set the initialization flag
   !-------------------------------------------------------------------------------------------------
   FDUnit      = UnWind
   PrevTime    = 0.0

   Initialized = .TRUE.

   RETURN

END SUBROUTINE FD_Init
!====================================================================================================
SUBROUTINE ReadFDP ( UnWind, FileName, FDTSfile, ErrStat )
!  This subroutine is used to read the input parameters for the large-eddy simulation.
!----------------------------------------------------------------------------------------------------

      ! Passed variables

   INTEGER,         INTENT(IN)    :: UnWind                       ! unit number for reading wind files
   CHARACTER(*),    INTENT(IN)    :: FileName                     ! Then name of the LE data file.
   CHARACTER(*),    INTENT(OUT)   :: FDTSfile                     ! The name of the file containing the time-step history of the wind files.
   INTEGER,         INTENT(OUT)   :: ErrStat                      ! return 0 if no errors encountered; non-zero otherwise


      ! Local variables

   CHARACTER(1024)                :: HeaderLine
!FIXME: this invokes SAVE
   CHARACTER(1),PARAMETER         :: Comp(3) = (/'U', 'V', 'W' /) ! the wind components

   REAL(ReKi)                     :: CoefTE                       ! Coefficient of thermal expansion.
   REAL(ReKi)                     :: DistScal                     ! Disturbance scale (ratio of wave height to rotor diameter) from input file.
   REAL(ReKi)                     :: Grav                         ! Gravitational acceleration.
   REAL(ReKi)                     :: LenScale                     ! Length scale (h).
   REAL(ReKi)                     :: Ri                           ! Richardson number.
   REAL(ReKi)                     :: Ubot                         ! Steady u-component wind speed at the bottom of the wave.
   REAL(ReKi)                     :: Zm_maxo                      ! The nondimensional vertical height of the untrimmed dataset.

   REAL(ReKi)                     :: Xm_max                       ! The nondimensional downwind length of the dataset.
   REAL(ReKi)                     :: Ym_max                       ! The nondimensional lateral width of the dataset.
   REAL(ReKi)                     :: Zm_max                       ! The nondimensional vertical height of the dataset.

   INTEGER                        :: I

   !-------------------------------------------------------------------------------------------------
   ! Open the 4D parameter file for reading
   !-------------------------------------------------------------------------------------------------
   CALL OpenFInpFile ( UnWind, TRIM( FileName ), ErrStat)
   IF (ErrStat /= 0) RETURN


   !-------------------------------------------------------------------------------------------------
   ! Read the 4D parameter input file
   !-------------------------------------------------------------------------------------------------

      !..............................................................................................
      ! Read the 4D wind parameters specific to this turbine simulation.
      !..............................................................................................

   CALL ReadStr( UnWind, TRIM( FileName ), HeaderLine, 'Header line', 'The header line in the FTP file', ErrStat )
   IF (ErrStat /= 0) RETURN
   CALL WrScr ( ' Heading of the 4D-wind-parameter file: "'//TRIM(HeaderLine)//'"' )


   CALL ReadCom( UnWind, TRIM( FileName ), 'Header line', ErrStat )
   IF (ErrStat /= 0) RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), FDSpath,  'FDSpath', 'Location (path) of the binary dataset', ErrStat )
   IF (ErrStat /= 0) RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), FDTSfile,  'FDTSfile', &
                                  'Name of the file containing the time-step history of the wind files', ErrStat )
   IF (ErrStat /= 0) RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), Ubot,  'Ubot', 'Steady u-component wind speed at the bottom of the wave', ErrStat )
   IF (ErrStat /= 0) RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), DistScal,  'DistScal', 'Disturbance scale', ErrStat )
   IF (ErrStat /= 0) RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), Lx,  'Lx', &
                            'Fractional location of tower centerline from upwind end to downwind end of the dataset', ErrStat )
   IF (ErrStat /= 0) RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), Ly,  'Ly', &
                 'Fractional location of tower centerline from right (looking downwind) to left side of the dataset', ErrStat )
   IF (ErrStat /= 0) RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), Lz,  'Lz', &
                                          'Fractional location of hub height from bottom to top of dataset', ErrStat )
   IF (ErrStat /= 0) RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), T_4D_St,  'T_4D_St', 'Time at which the wave event starts', ErrStat )
   IF (ErrStat /= 0) RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), ScaleVel,  'ScaleVel', &
                 'Scaling velocity, U0: half the difference in wind speed between the top and bottom of the billow.', ErrStat )
   IF (ErrStat /= 0) RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), RotDiam,  'RotDiam', 'Rotor diameter', ErrStat )
   IF (ErrStat /= 0) RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), FD_DF_X,  'FD_DF_X', 'Decimation factor in X direction', ErrStat )
   IF (ErrStat /= 0) RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), FD_DF_Y,  'FD_DF_Y', 'Decimation factor in Y direction', ErrStat )
   IF (ErrStat /= 0) RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), FD_DF_Z,  'FD_DF_Z', 'Decimation factor in Z direction', ErrStat )
   IF (ErrStat /= 0) RETURN


   CALL ReadCom( UnWind, TRIM( FileName ), 'blank line', ErrStat )
   IF (ErrStat /= 0) RETURN

      !..............................................................................................
      ! Read the 4D wind parameters specific to the K-H billow simulation being used.
      !..............................................................................................

   CALL ReadCom( UnWind, TRIM( FileName ), 'LES parameters specific to the K-H billow simulation being used', ErrStat )
   IF (ErrStat /= 0) RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), VertShft,  'VertShft', &
                           'Flag to indicate whether or not to shift the z values for the w component', ErrStat )
   IF (ErrStat /= 0) RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), Xm_max,  'Xm_max', &
                           'Maximum nondimensional downwind distance from center of dataset', ErrStat )
   IF (ErrStat /= 0) RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), Ym_max,  'Ym_max', &
                           'Maximum nondimensional lateral distance from center of dataset', ErrStat )
   IF (ErrStat /= 0) RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), Zm_max,  'Zm_max', &
                           'Maximum nondimensional vertical distance from center of dataset', ErrStat )
   IF (ErrStat /= 0) RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), Zm_maxo,  'Zm_maxo', &
                 'Maximum nondimensional vertical distance from center of untrimmed dataset', ErrStat )
   IF (ErrStat /= 0) RETURN


   DO I = 1,3

      CALL ReadVar( UnWind, TRIM( FileName ), ScalFact(I),  Comp(I)//'Scl', &
                    Comp(I)//'-component scale factor for converting from integers to reals', ErrStat )
      IF (ErrStat /= 0) RETURN
      ScalFact(I) = ScalFact(I) * ScaleVel


      CALL ReadVar( UnWind, TRIM( FileName ), Offsets(I), Comp(I)//'Off', &
                    Comp(I)//'-component offset for converting from integers to reals', ErrStat )
      IF (ErrStat /= 0) RETURN
      Offsets(I) = Offsets(I) * ScaleVel

   END DO
   Offsets (1) = Offsets (1) + ScaleVel + Ubot                           ! u-component offset to convert integer data to actual wind speeds.


   CALL ReadVar( UnWind, TRIM( FileName ), Num4Dt, 'Num4Dt', 'The number of LE grids, one grid per time step', ErrStat )
   IF (ErrStat /= 0) RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), Num4Dx, 'Num4Dx', 'The number of LE grid points in the x direction', ErrStat )
   IF (ErrStat /= 0) RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), Num4Dy, 'Num4Dy', 'The number of LE grid points in the y direction', ErrStat )
   IF (ErrStat /= 0) RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), Num4Dz, 'Num4Dz', 'The number of LE grid points in the z direction', ErrStat )
   IF (ErrStat /= 0) RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), Ri, 'Ri', 'Richardson number', ErrStat )
   IF (ErrStat /= 0) RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), CoefTE, 'CoefTE', 'Coefficient of thermal expansion', ErrStat )
   IF (ErrStat /= 0) RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), Grav, 'Grav', 'Gravitational acceleration', ErrStat )
   IF (ErrStat /= 0) RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), Advect, 'Advect', 'Advection flag', ErrStat )

   IF (ErrStat /= 0) THEN

      Advect   = .FALSE.
      Ind4DAdv = 0
      ErrStat  = 0
      CALL WrScr( ' Advection will not be used.')

   ELSE

      IF (Advect) THEN
         IF ( FD_DF_X /= 1 ) THEN
            CALL WrScr( ' FD_DF_X must be 1 when using advection. ' )
            FD_DF_X = 1
         ENDIF

         CALL ReadVar( UnWind, TRIM( FileName ), NumAdvect, 'NumAdvect', 'Number of 4D files for advection', ErrStat )
         IF (ErrStat /= 0) RETURN


         IF ( NumAdvect < 1 ) THEN
            CALL WrScr( ' NumAdvect in 4D-wind-parameter file, "'//TRIM( FileName )//'," must be at least 1.' )
            ErrStat = 1
            RETURN
         ENDIF

         IF ( .NOT. ALLOCATED( AdvFiles ) ) THEN
!            CALL AllocAry( AdvFiles, NumAdvect, 'AdvFiles array', ErrStat )
            ALLOCATE ( AdvFiles(NumAdvect), STAT=ErrStat )

            IF ( ErrStat /= 0 )  THEN
               CALL WrScr ( ' Error allocating memory for the AdvFiles array.' )
               RETURN
            END IF
         ENDIF

         CALL ReadAryLines( UnWind, TRIM( FileName ), AdvFiles, NumAdvect, 'AdvFiles', 'Advection file names', ErrStat )
         IF (ErrStat /= 0) RETURN
         Ind4DAdv = 1

      ELSE
         Ind4DAdv = 0
      ENDIF !Advect == .TRUE.

   END IF

   !-------------------------------------------------------------------------------------------------
   ! Close the 4D parameter input file
   !-------------------------------------------------------------------------------------------------
   CLOSE ( UnWind )

   !-------------------------------------------------------------------------------------------------
   ! Close the 4D parameter input file
   !-------------------------------------------------------------------------------------------------

   LenScale    = RotDiam*DistScal/Zm_max                             ! Length scale (h).
   Xmax        = Xm_max*LenScale                                     ! The dimensional length of the dataset.
   Ymax        = Ym_max*LenScale                                     ! The dimensional width of the dataset
   Zmax        = Zm_max*LenScale                                     ! The dimensional vertical height of the dataset.
   TSclFact    = LenScale/ScaleVel                                   ! Scale factor for time (h/U0).



   RETURN

END SUBROUTINE ReadFDP
!====================================================================================================
SUBROUTINE Read4Dtimes ( UnWind, FileName, ErrStat )
!  This subroutine is used to read the time array for the 4D data.  The times in the file are
!  non-dimensional and non-uniformly spaced. They are scaled using TSclFact to obtain units of seconds
!  and T_4D_St is added to allow the billow to start at non-zero time.
!----------------------------------------------------------------------------------------------------

      ! Passed variables

   INTEGER,         INTENT(IN)    :: UnWind                       ! unit number for reading wind files
   CHARACTER(*),    INTENT(IN)    :: FileName                     ! Then name of the LE data file.
   INTEGER,         INTENT(OUT)   :: ErrStat                      ! return 0 if no errors encountered; non-zero otherwise


      ! Local variables

   INTEGER                        :: I                            ! Loop counter

   !-------------------------------------------------------------------------------------------------
   ! Allocate arrays to store the data in
   !-------------------------------------------------------------------------------------------------

   IF (.NOT. ALLOCATED( Times4D) ) THEN
!      CALL AllocAry( Times4D, Num4Dt, '4D time array', ErrStat)
      ALLOCATE ( Times4D(Num4Dt), STAT=ErrStat )

      IF ( ErrStat /= 0 )  THEN
         CALL WrScr ( ' Error allocating memory for the Times4D array.' )
         RETURN
      END IF
   END IF

   IF (.NOT. ALLOCATED( Times4DIx) ) THEN
!      CALL AllocAry( Times4DIx, Num4Dt, '4D time array', ErrStat)
      ALLOCATE ( Times4DIx(Num4Dt), STAT=ErrStat )

      IF ( ErrStat /= 0 )  THEN
         CALL WrScr ( ' Error allocating memory for the Times4DIx array.' )
         RETURN
      END IF
   END IF


   !-------------------------------------------------------------------------------------------------
   ! Open the 4D times file
   !-------------------------------------------------------------------------------------------------
   CALL OpenFInpFile ( UnWind, TRIM( FileName ), ErrStat)
   IF ( ErrStat /= 0 ) RETURN

   !-------------------------------------------------------------------------------------------------
   ! Read the 4D times file
   !-------------------------------------------------------------------------------------------------
   CALL ReadCom( UnWind, TRIM( FileName ), 'first line', ErrStat)
   IF ( ErrStat /= 0 ) RETURN

   DO I=1,Num4Dt

      READ (UnWind,*,IOSTAT=ErrStat)  Times4DIx(I), Times4D(I)

      IF ( ErrStat /= 0 )  THEN

         CALL WrScr( ' Error reading line '//TRIM( Num2LStr( I+1 ) )// &
                        ' of the 4D-wind time-steps file, "'//TRIM( FileName )//'."')
         RETURN

      ENDIF

   ENDDO ! I


   !-------------------------------------------------------------------------------------------------
   ! Close the 4D times file
   !-------------------------------------------------------------------------------------------------

   CLOSE ( UnWind )

   RETURN

END SUBROUTINE Read4Dtimes
!====================================================================================================
SUBROUTINE ReadAll4DData(UnWind, ErrStat)
! This subroutine reads the data into one array to be accessed later when ADVECT=.TRUE. Since there
! are just a few time steps, we'll load them into memory to (hopefully) save I/O time.
!----------------------------------------------------------------------------------------------------

   INTEGER, INTENT(IN)        :: UnWind
   INTEGER, INTENT(OUT)       :: ErrStat                            !
   INTEGER                    :: IT

   CHARACTER(1)               :: FDNum
   CHARACTER(20)              :: DNSFileName                        ! String containing part of the current file name.


   DO IT = 1,Num4Dt

      WRITE(FDNum,'(I1.1)') Times4DIx(IT)
      DNSFileName = TRIM(AdvFiles(Ind4DAdv))//'_'//TRIM(FDNum)//'.dns'

      CALL Read4DData ( UnWind, TRIM( FDSpath )//'\u\u_16i_'//TRIM(DNSFileName), FDuData, IT, ScalFact(1), Offsets(1), ErrStat )
      IF ( ErrStat /= 0 ) RETURN

      CALL Read4DData ( UnWind, TRIM( FDSpath )//'\v\v_16i_'//TRIM(DNSFileName), FDvData, IT, ScalFact(2), Offsets(2), ErrStat )
      IF ( ErrStat /= 0 ) RETURN

      CALL Read4DData ( UnWind, TRIM( FDSpath )//'\w\w_16i_'//TRIM(DNSFileName), FDwData, IT, ScalFact(3), Offsets(3), ErrStat )
      IF ( ErrStat /= 0 ) RETURN

   ENDDO ! IT

   Ind4DAdv = Ind4DAdv + 1

   RETURN

END SUBROUTINE ReadAll4DData
!====================================================================================================
SUBROUTINE LoadLESData( UnWind, FileNo, Indx, ErrStat )
! This subroutine reads binary data from the U, V, and W files and stores them in the arrays FDu,
! FDv, and FDw (by calling Read4DData).
!----------------------------------------------------------------------------------------------------
      ! Passed variables

   INTEGER,         INTENT(IN)    :: UnWind                       ! unit number for reading wind files
   INTEGER,         INTENT(IN)    :: FileNo                       ! current file number to read
   INTEGER,         INTENT(IN)    :: Indx                         ! index into the data arrays
   INTEGER,         INTENT(OUT)   :: ErrStat                      ! return 0 if no errors encountered; non-zero otherwise

      ! local variables
   CHARACTER(5)                   :: FDNum
   CHARACTER(20)                  :: LESFileName                  ! String containing part of the current file name.


      ! get the file name for the file number

   WRITE(FDNum,'(I5.5)', IOStat=ErrStat) FileNo
   IF ( ErrStat /= 0 ) RETURN

   LESFileName = TRIM(FDNum)//'.les'


      ! set the paths and read the data for each component

   CALL Read4DData ( UnWind, TRIM( FDSpath )//'\u\u_16i_'//TRIM(LESFileName), FDu, Indx, ScalFact(1), Offsets(1), ErrStat )
   IF ( ErrStat /= 0 ) RETURN

   CALL Read4DData ( UnWind, TRIM( FDSpath )//'\v\v_16i_'//TRIM(LESFileName), FDv, Indx, ScalFact(2), Offsets(2), ErrStat )
   IF ( ErrStat /= 0 ) RETURN

   CALL Read4DData ( UnWind, TRIM( FDSpath )//'\w\w_16i_'//TRIM(LESFileName), FDw, Indx, ScalFact(3), Offsets(3), ErrStat )


END SUBROUTINE LoadLESData
!====================================================================================================
SUBROUTINE Read4DData ( UnWind, FileName, Comp, Indx4, Scale, Offset,  ErrStat)
! This subroutine is used to read one time-step's worth of large-eddy wind data for one component
! from a file.
!----------------------------------------------------------------------------------------------------

      ! Passed variables

   INTEGER,     INTENT(IN)    :: UnWind               ! The I/O unit of the LE file.
   CHARACTER(*),INTENT(IN)    :: FileName             ! Then name of the LE data file.

   REAL(ReKi),  INTENT(INOUT) :: Comp (:,:,:,:)       ! The velocity array [do NOT make this INTENT(OUT): other parts of the array may become undefined]
   INTEGER,     INTENT(IN)    :: Indx4                ! The index of the 4th dimension of Comp, which is to be read.
   REAL(ReKi),  INTENT(IN)    :: Scale                ! The scale factor for converting from intergers to non-normalized reals.
   REAL(ReKi),  INTENT(IN)    :: Offset               ! The offset for converting from intergers to non-normalized reals.

   INTEGER,     INTENT(OUT)   :: ErrStat              ! The returned status of a READ.

      ! Local variables

   INTEGER                    :: IX                   ! A DO index for indexing the arrays in the x direction.
   INTEGER                    :: IXK                  ! An index for the decimated arrays in the x direction.
   INTEGER                    :: IY                   ! A DO index for indexing the arrays in the y direction.
   INTEGER                    :: IYK                  ! An index for the decimated arrays in the y direction.
   INTEGER                    :: IZ                   ! A DO index for indexing the arrays in the z direction.
   INTEGER                    :: IZK                  ! An index for the decimated arrays in the z direction.

   INTEGER(B2Ki)              :: Com (Num4Dx,Num4Dy)  ! Temporary array to hold component's integer values for a given Z.


   !-------------------------------------------------------------------------------------------------
   ! Open the binary input file
   !-------------------------------------------------------------------------------------------------
   CALL OpenUInBEFile( UnWind, TRIM( FileName ), FDRecL, ErrStat )
   IF ( ErrStat /= 0 ) RETURN


   !-------------------------------------------------------------------------------------------------
   ! Read the input file
   !-------------------------------------------------------------------------------------------------

   IZK = 0
   DO IZ=1,Num4Dz,FD_DF_Z

      READ (UnWind,REC=IZ,IOSTAT=ErrStat)  Com

      IF ( ErrStat /= 0 )  THEN

         CALL WrScr( ' Error reading record '//TRIM( Num2LStr( IZ ) )// &
                                            ' of the binary 4D wind file, "'//TRIM( FileName )//'".')
         RETURN

      ENDIF

      IZK = IZK + 1                                ! IZK = ( IZ - 1 + FD_DF_Z )/FD_DF_Z
      IYK = 0

      DO IY=1,Num4Dy,FD_DF_Y

         IYK = IYK + 1                             ! IYK = ( IY - 1 + FD_DF_Y )/FD_DF_Y

         DO IX=1,Num4Dx,FD_DF_X

               ! shift the x-index, if necessary, to perform Advection

            !IXK = ( IX + FD_DF_X - 1 )/FD_DF_X
            IXK = ( MOD(IX+Shft4Dnew-1,Num4Dx) + FD_DF_X )/FD_DF_X

            Comp(IXK,IYK,IZK,Indx4) = Scale*Com(IX,IY) + Offset

         ENDDO ! IX

      ENDDO ! IY

   ENDDO ! IZ


   !-------------------------------------------------------------------------------------------------
   ! Close the file
   !-------------------------------------------------------------------------------------------------
   CLOSE ( UnWind )

   RETURN

END SUBROUTINE Read4DData
!====================================================================================================
SUBROUTINE Load4DData( InpIndx )
! This subroutine takes the data from the storage array (used when ADVECT=.TRUE., shifts it if necessary,
! and loads it into the array for the time slice indexed by InpIndx.
!----------------------------------------------------------------------------------------------------

   INTEGER, INTENT(IN) :: InpIndx

   INTEGER             :: IX
   INTEGER             :: IXK


   DO IX=1,Num4Dx,FD_DF_X

         ! shift the x-index, if necessary, to perform Advection
      IXK = ( MOD(IX+Shft4Dnew-1,Num4Dx) + FD_DF_X )/FD_DF_X

      FDu(IXK,:,:,InpIndx) = FDuData(IX,:,:,FDFileNo)
      FDv(IXK,:,:,InpIndx) = FDvData(IX,:,:,FDFileNo)
      FDw(IXK,:,:,InpIndx) = FDwData(IX,:,:,FDFileNo)

   ENDDO ! IX


   RETURN

END SUBROUTINE Load4DData
!====================================================================================================
FUNCTION FD_GetValue(RVarName, ErrStat)
!  This function returns a real scalar value whose name is listed in the RVarName input argument.
!  If the name is not recognized, an error is returned in ErrStat.
!----------------------------------------------------------------------------------------------------

   CHARACTER(*),   INTENT(IN)    :: RVarName
   INTEGER,        INTENT(OUT)   :: ErrStat
   REAL(ReKi)                    :: FD_GetValue


   CHARACTER(20)                 :: VarNameUC


   !-------------------------------------------------------------------------------------------------
   ! Check that the module has been initialized.
   !-------------------------------------------------------------------------------------------------

   IF ( .NOT. Initialized ) THEN
      CALL WrScr( ' Initialialize the FDWind module before calling its subroutines.' )
      ErrStat = 1
      RETURN
   ELSE
      ErrStat = 0
   END IF


   !-------------------------------------------------------------------------------------------------
   ! Return the requested values.
   !-------------------------------------------------------------------------------------------------

   VarNameUC = RVarName
   CALL Conv2UC( VarNameUC )

   SELECT CASE ( TRIM(VarNameUC) )

      CASE ('ROTDIAM' )
         FD_GetValue = RotDiam

      CASE DEFAULT
         CALL WrScr( ' Invalid variable name in FD_GetRValue().' )
         ErrStat = 1

   END SELECT

END FUNCTION FD_GetValue
!====================================================================================================
FUNCTION FD_GetWindSpeed(Time, InputPosition, ErrStat)
! This function is used to interpolate into the 4D wind arrays.  It receives X, Y, Z and TIME from the
! calling routine.  The time since the start of the 4D data is used to decide which pair of time slices
! to interpolate within and between.  After finding the two time slices, it decides which eight grid
! points bound the (X,Y,Z) pair. It does a trilinear interpolation for each time slice. Linear
! interpolation is then used to interpolate between time slices.  This routine assumes that X is
! downwind, Y is to the left when looking downwind and Z is up.  It also assumes that no
! extrapolation will be needed except in time and the Z direction.  In those cases, the appropriate
! steady winds are used.
!----------------------------------------------------------------------------------------------------

      ! Passed variables:

   REAL(DbKi),        INTENT(IN) :: Time                                   ! the time
   REAL(ReKi),        INTENT(IN) :: InputPosition(3)                       ! structure that contains the position
   INTEGER,           INTENT(OUT):: ErrStat                                ! returns 0 if no error; non-zero otherwise
!FIXME:delete
!   TYPE(InflIntrpOut)            :: FD_GetWindSpeed                        ! the resultant wind speed
   REAL(ReKi)                 :: FD_GetWindSpeed(3)                        ! the resultant wind speed


      ! Local Variables:

   REAL(ReKi)                 :: Ixhyz                                     ! Temporary interpolated value.
   REAL(ReKi)                 :: Ixlyz                                     ! Temporary interpolated value.
   REAL(ReKi)                 :: Ixyzo                                     ! Temporary interpolated value.
   REAL(ReKi)                 :: Iyhz                                      ! Temporary interpolated value.
   REAL(ReKi)                 :: Iylz                                      ! Temporary interpolated value.
   REAL(ReKi)                 :: Ixyzn                                     ! Temporary interpolated value.
   REAL(ReKi)                 :: Tgrid                                     ! Fractional distance between time grids.
   REAL(ReKi)                 :: Xgrid                                     ! Fractional distance between grids in the x direction.
   REAL(ReKi)                 :: Xnorm                                     ! Nondimensional downwind distance of the analysis point from upwind end of dataset.
   REAL(ReKi)                 :: Ygrid                                     ! Fractional distance between grids in the y direction.
   REAL(ReKi)                 :: Ynorm                                     ! Nondimensional lateral distance of the analysis point from right side of dataset (looking downwind).
   REAL(ReKi)                 :: Zgrid                                     ! Fractional distance between grids in the z direction.
   REAL(ReKi)                 :: Zgrid_w                                   ! Fractional distance between grids in the z direction for the w component.
   REAL(ReKi)                 :: Znorm                                     ! Nondimensional vertical distance of the analysis point from bottom of dataset.
   REAL(ReKi)                 :: Znorm_w                                   ! Nondimensional vertical distance of the analysis point from bottom of dataset for the w component.

   INTEGER                    :: IT                                        ! Index for do loop
   INTEGER                    :: IXHI                                      ! Index for the more-positive x value.
   INTEGER                    :: IXLO                                      ! Index for the more-negative x value.
   INTEGER                    :: IYHI                                      ! Index for the more-positive y value.
   INTEGER                    :: IYLO                                      ! Index for the more-negative y value.
   INTEGER                    :: IZHI                                      ! Index for the more-positive z value.
   INTEGER                    :: IZHI_w                                    ! Index for the more-positive z value for the w component.
   INTEGER                    :: IZLO                                      ! Index for the more-negative z value.
   INTEGER                    :: IZLO_w                                    ! Index for the more-negative z value for the w component.

   REAL(ReKi)                    :: TempWindSpeed(3)                       ! Temporary variable to hold the windspeed before returning

   !-------------------------------------------------------------------------------------------------
   ! Check that we've initialized everything first
   !-------------------------------------------------------------------------------------------------

   IF ( .NOT. Initialized ) THEN
      CALL WrScr( ' Initialialize the FDWind module before calling its subroutines.' )
      ErrStat = 1
      RETURN
   ELSE
      ErrStat = 0
   END IF

   !-------------------------------------------------------------------------------------------------
   ! If the TIME is greater than the time for the last file read, read another set of files until we straddle the current time.
   ! Stick with the last file if we've exhausted the data.
   ! We're assuming here that the simulation time step is smaller than the wind-file time step.
   !-------------------------------------------------------------------------------------------------

   IF ( Time < PrevTime .AND. Time < FDTime(Ind4Dold) ) THEN  ! bjj: GET THE CORRECT TIME if we're going backward!

      !----------------------------------------------------------------------------------------------
      ! Determine the first file needed for this simulation.
      !----------------------------------------------------------------------------------------------
      Ind4Dold  = 1                                                           ! Put the old stuff in the first part of the array.
      Ind4Dnew  = 2                                                           ! Put the new stuff in the second part of the array.

      FDFileNo  = Num4Dt
      DO IT=1,Num4Dt
         IF ( Times4D(IT) > Time )  THEN
            FDFileNo = IT - 1
            EXIT
         END IF
      END DO ! IT

      !----------------------------------------------------------------------------------------------
      ! Open, read, and close the first set of files.
      !----------------------------------------------------------------------------------------------
      FDTime(Ind4Dold) = Times4D(FDFileNo)                                 ! Set the time for this file.

      IF ( ADVECT ) THEN
         CALL Load4DData(Ind4Dold)     ! load data stored in FDuData, FDvData, and FDwData arrays
      ELSE
         CALL LoadLESData( FDUnit, FDFileNo, Ind4Dold, ErrStat )
      END IF

      !----------------------------------------------------------------------------------------------
      ! Open, read, and close the second set of files.
      !----------------------------------------------------------------------------------------------
      FDFileNo  = MIN(FDFileNo + 1, Num4Dt)
      Shft4Dnew = 0

      IF ( ADVECT ) THEN
         FDFileNo = MOD(FDFileNo-1,Num4Dt) + 1

         IF (FDFileNo == 1) THEN
            Shft4Dnew = Shft4Dnew + 1

            IF (Ind4DAdv <= NumAdvect) THEN                             ! Ind4DAdv was set in ReadFDP
               IF ( MOD( Shft4Dnew, Num4Dx ) == 0 ) THEN
                  CALL ReadAll4DData(FDUnit, ErrStat)
                  IF ( ErrStat /= 0 ) RETURN
               END IF
            END IF

         ENDIF

         FDTime(Ind4Dnew) = Times4D(FDFileNo) + Shft4Dnew*FDPer         ! Set the time for this file.

         CALL Load4DData( Ind4Dnew )    ! shift the data

      ELSE
         FDTime(Ind4Dnew) = Times4D(FDFileNo)                                           ! Set the time for this file.
!
         CALL LoadLESData( FDUnit, FDFileNo, Ind4Dnew, ErrStat )
      ENDIF

   END IF

   !-------------------------------------------------------------------------------------------------
   ! Move forward in time
   !-------------------------------------------------------------------------------------------------

   DO WHILE ( Time > FDTime(Ind4Dnew) .AND. ( Time < T_4D_En .OR. ADVECT ) )

      Ind4Dnew         = Ind4Dold                                          ! Reverse array indices (1 or 2).
      Ind4Dold         = 3 - Ind4Dnew
      FDFileNo         = FDFileNo + 1                                      ! Increment file number.


      IF ( ADVECT ) THEN
         FDFileNo = MOD(FDFileNo-1,Num4Dt) + 1

         IF (FDFileNo == 1) THEN
               Shft4Dnew = Shft4Dnew + 1

               IF (Ind4DAdv <= NumAdvect) THEN
                  IF ( MOD( Shft4Dnew, Num4Dx ) == 0 ) THEN
                     CALL ReadAll4DData(FDUnit, ErrStat)
                     IF ( ErrStat /= 0 ) RETURN
                  END IF
               ENDIF

         ENDIF

         FDTime(Ind4Dnew) = Times4D(FDFileNo) + Shft4Dnew*FDPer

         CALL Load4DData( Ind4Dnew )  ! shift the data
      ELSE
         FDTime(Ind4Dnew) = Times4D(FDFileNo)

         CALL LoadLESData( FDUnit, FDFileNo, Ind4Dnew, ErrStat )
      ENDIF

   ENDDO


   !.................................................................................................
   ! Find the bounding rows, columns, and planes for the X,Y,Z position.  The near, lower-right
   ! corner is (1,1,1) when looking downwind. Make sure the lowest possible value is 1.
   !.................................................................................................


   !-------------------------------------------------------------------------------------------------
   ! get values of Time for interpolation. Linear interpolation; Nearest-neighbor extrapolation.
   !-------------------------------------------------------------------------------------------------

      ! Find out fractionally how far we are between grids in time and between grid points in each direction.
      !  Limit values to avoid extrapolation.  We need this for interpolation later on.

   Tgrid = MIN( MAX( ( Time - FDTime(Ind4Dold) )/( FDTime(Ind4Dnew) - FDTime(Ind4Dold) ), 0.0 ), 1.0 )


   !-------------------------------------------------------------------------------------------------
   ! get values of X for interpolation. Grid is periodic in X.
   !-------------------------------------------------------------------------------------------------
   Xnorm = ( Xt + InputPosition(1) )/Xmax

   DO WHILE ( Xnorm < 0.0 )   ! Ensure Xnorm is not negative.  The wave is periodic in x.
      Xnorm = Xnorm + 1.0
   ENDDO

   Xgrid = MIN( MAX( MOD( Xnorm, DelXgrid ), 0.0 ), 1.0 )
   IXLo  = MAX( MOD( INT( Xnorm*Num4DxD1 ) + 1, Num4DxD1 ), 1 )
   IXHi  = MOD( IXLo, Num4DxD ) + 1

   !-------------------------------------------------------------------------------------------------
   ! get values of Y for interpolation. Grid is periodic in Y.
   !-------------------------------------------------------------------------------------------------
   Ynorm = ( Yt + InputPosition(2) )/Ymax

   DO WHILE ( Ynorm < 0.0 )  ! Ensure Ynorm is not negative.  The wave is periodic in y.
      Ynorm = Ynorm + 1.0
   ENDDO

   Ygrid = MIN( MAX( MOD( Ynorm, DelYgrid ), 0.0 ), 1.0 )
   IYLo  = MAX( MOD( INT( Ynorm*Num4DyD1 ) + 1, Num4DyD1 ), 1 )
   IYHi  = MOD( IYLo, Num4DyD ) + 1

   !-------------------------------------------------------------------------------------------------
   ! get values of Z for interpolation.  Linear interpolation; Nearest-neighbor extrapolation.
   !-------------------------------------------------------------------------------------------------
   Znorm = MIN( MAX( ( Zt + InputPosition(3) - ZRef )/Zmax, 0.0 ), 1.0 ) !bjj: define ZRef

   Zgrid = MIN( MAX( MOD( Znorm, DelZgrid ), 0.0 ), 1.0 )
   IZLo  = MAX( INT( Znorm*Num4DzD1 ) + 1, 1 )

      ! If we are located at the upper end of the z dimension, decrement the index by one and set the grid coordinate to 1.

   IF ( IZLo == Num4DzD )  THEN
      IZLo  = Num4DzD1
      Zgrid = 1.0
   ENDIF
   IZHi = IZLo + 1

      !..............................................................................................
      ! Find the equivalent Znorm (Znorm_w) for the w-component, which may be shifted vertically
      ! by half the original grid spacing.
      !..............................................................................................

   IF ( VertShft ) THEN
      Znorm_w = MAX( Znorm - 0.5*DelZgrid/FD_DF_Z, 0.0 )
   ELSE
      Znorm_w = Znorm
   ENDIF

   Zgrid_w = MIN( MAX( MOD( Znorm_w, DelZgrid ), 0.0 ), 1.0 )
   IZLo_w  = MAX( INT( Znorm_w*Num4DzD1 ) + 1, 1 )

   IF ( IZLo_w == Num4DzD )  THEN
      IZLo_w  = Num4DzD1
      Zgrid_w = 1.0
   ENDIF

   IZHi_w = IZLo_w + 1


   !-------------------------------------------------------------------------------------------------
   ! Interpolate for u component of wind within the grid.
   !-------------------------------------------------------------------------------------------------

   Iylz  = ( FDu(IXLo,IYLo,IZHi,Ind4Dold) - FDu(IXLo,IYLo,IZLo,Ind4Dold) )*Zgrid + FDu(IXLo,IYLo,IZLo,Ind4Dold)
   Iyhz  = ( FDu(IXLo,IYHi,IZHi,Ind4Dold) - FDu(IXLo,IYHi,IZLo,Ind4Dold) )*Zgrid + FDu(IXLo,IYHi,IZLo,Ind4Dold)
   Ixlyz = ( Iyhz - Iylz )*Ygrid + Iylz

   Iylz  = ( FDu(IXHi,IYLo,IZHi,Ind4Dold) - FDu(IXHi,IYLo,IZLo,Ind4Dold) )*Zgrid + FDu(IXHi,IYLo,IZLo,Ind4Dold)
   Iyhz  = ( FDu(IXHi,IYHi,IZHi,Ind4Dold) - FDu(IXHi,IYHi,IZLo,Ind4Dold) )*Zgrid + FDu(IXHi,IYHi,IZLo,Ind4Dold)
   Ixhyz = ( Iyhz - Iylz )*Ygrid + Iylz

   Ixyzo = ( Ixhyz - Ixlyz )*Xgrid + Ixlyz

   Iylz  = ( FDu(IXLo,IYLo,IZHi,Ind4Dnew) - FDu(IXLo,IYLo,IZLo,Ind4Dnew) )*Zgrid + FDu(IXLo,IYLo,IZLo,Ind4Dnew)
   Iyhz  = ( FDu(IXLo,IYHi,IZHi,Ind4Dnew) - FDu(IXLo,IYHi,IZLo,Ind4Dnew) )*Zgrid + FDu(IXLo,IYHi,IZLo,Ind4Dnew)
   Ixlyz = ( Iyhz - Iylz )*Ygrid + Iylz

   Iylz  = ( FDu(IXHi,IYLo,IZHi,Ind4Dnew) - FDu(IXHi,IYLo,IZLo,Ind4Dnew) )*Zgrid + FDu(IXHi,IYLo,IZLo,Ind4Dnew)
   Iyhz  = ( FDu(IXHi,IYHi,IZHi,Ind4Dnew) - FDu(IXHi,IYHi,IZLo,Ind4Dnew) )*Zgrid + FDu(IXHi,IYHi,IZLo,Ind4Dnew)
   Ixhyz = ( Iyhz - Iylz )*Ygrid + Iylz

   Ixyzn = ( Ixhyz - Ixlyz )*Xgrid + Ixlyz

!   FD_GetWindSpeed%Velocity(1) = ( Ixyzn - Ixyzo )*Tgrid + Ixyzo
   TempWindSpeed(1) = ( Ixyzn - Ixyzo )*Tgrid + Ixyzo

   !-------------------------------------------------------------------------------------------------
   ! Interpolate for v component of wind within the grid.
   !-------------------------------------------------------------------------------------------------

   Iylz  = ( FDv(IXLo,IYLo,IZHi,Ind4Dold) - FDv(IXLo,IYLo,IZLo,Ind4Dold) )*Zgrid + FDv(IXLo,IYLo,IZLo,Ind4Dold)
   Iyhz  = ( FDv(IXLo,IYHi,IZHi,Ind4Dold) - FDv(IXLo,IYHi,IZLo,Ind4Dold) )*Zgrid + FDv(IXLo,IYHi,IZLo,Ind4Dold)
   Ixlyz = ( Iyhz - Iylz )*Ygrid + Iylz

   Iylz  = ( FDv(IXHi,IYLo,IZHi,Ind4Dold) - FDv(IXHi,IYLo,IZLo,Ind4Dold) )*Zgrid + FDv(IXHi,IYLo,IZLo,Ind4Dold)
   Iyhz  = ( FDv(IXHi,IYHi,IZHi,Ind4Dold) - FDv(IXHi,IYHi,IZLo,Ind4Dold) )*Zgrid + FDv(IXHi,IYHi,IZLo,Ind4Dold)
   Ixhyz = ( Iyhz - Iylz )*Ygrid + Iylz

   Ixyzo = ( Ixhyz - Ixlyz )*Xgrid + Ixlyz

   Iylz  = ( FDv(IXLo,IYLo,IZHi,Ind4Dnew) - FDv(IXLo,IYLo,IZLo,Ind4Dnew) )*Zgrid + FDv(IXLo,IYLo,IZLo,Ind4Dnew)
   Iyhz  = ( FDv(IXLo,IYHi,IZHi,Ind4Dnew) - FDv(IXLo,IYHi,IZLo,Ind4Dnew) )*Zgrid + FDv(IXLo,IYHi,IZLo,Ind4Dnew)
   Ixlyz = ( Iyhz - Iylz )*Ygrid + Iylz

   Iylz  = ( FDv(IXHi,IYLo,IZHi,Ind4Dnew) - FDv(IXHi,IYLo,IZLo,Ind4Dnew) )*Zgrid + FDv(IXHi,IYLo,IZLo,Ind4Dnew)
   Iyhz  = ( FDv(IXHi,IYHi,IZHi,Ind4Dnew) - FDv(IXHi,IYHi,IZLo,Ind4Dnew) )*Zgrid + FDv(IXHi,IYHi,IZLo,Ind4Dnew)
   Ixhyz = ( Iyhz - Iylz )*Ygrid + Iylz

   Ixyzn = ( Ixhyz - Ixlyz )*Xgrid + Ixlyz

!FIXME:delete
!   FD_GetWindSpeed%Velocity(2) = ( Ixyzn - Ixyzo )*Tgrid + Ixyzo
   TempWindSpeed(2) = ( Ixyzn - Ixyzo )*Tgrid + Ixyzo

   !-------------------------------------------------------------------------------------------------
   ! Interpolate for w component of wind within the grid.
   !-------------------------------------------------------------------------------------------------
   !bjj: should Zgrid actually be Zgrid_w here?  I changed it so that it's consistent

   Iylz  = ( FDw(IXLo,IYLo,IZHi_w,Ind4Dold) - FDw(IXLo,IYLo,IZLo_w,Ind4Dold) )*Zgrid_w + FDw(IXLo,IYLo,IZLo_w,Ind4Dold)
   Iyhz  = ( FDw(IXLo,IYHi,IZHi_w,Ind4Dold) - FDw(IXLo,IYHi,IZLo_w,Ind4Dold) )*Zgrid_w + FDw(IXLo,IYHi,IZLo_w,Ind4Dold)
   Ixlyz = ( Iyhz - Iylz )*Ygrid + Iylz

   Iylz  = ( FDw(IXHi,IYLo,IZHi_w,Ind4Dold) - FDw(IXHi,IYLo,IZLo_w,Ind4Dold) )*Zgrid_w + FDw(IXHi,IYLo,IZLo_w,Ind4Dold)
   Iyhz  = ( FDw(IXHi,IYHi,IZHi_w,Ind4Dold) - FDw(IXHi,IYHi,IZLo_w,Ind4Dold) )*Zgrid_w + FDw(IXHi,IYHi,IZLo_w,Ind4Dold)
   Ixhyz = ( Iyhz - Iylz )*Ygrid + Iylz

   Ixyzo = ( Ixhyz - Ixlyz )*Xgrid + Ixlyz

   Iylz  = ( FDw(IXLo,IYLo,IZHi_w,Ind4Dnew) - FDw(IXLo,IYLo,IZLo_w,Ind4Dnew) )*Zgrid_w + FDw(IXLo,IYLo,IZLo_w,Ind4Dnew)
   Iyhz  = ( FDw(IXLo,IYHi,IZHi_w,Ind4Dnew) - FDw(IXLo,IYHi,IZLo_w,Ind4Dnew) )*Zgrid_w + FDw(IXLo,IYHi,IZLo_w,Ind4Dnew)
   Ixlyz = ( Iyhz - Iylz )*Ygrid + Iylz

   Iylz  = ( FDw(IXHi,IYLo,IZHi_w,Ind4Dnew) - FDw(IXHi,IYLo,IZLo_w,Ind4Dnew) )*Zgrid_w + FDw(IXHi,IYLo,IZLo_w,Ind4Dnew)
   Iyhz  = ( FDw(IXHi,IYHi,IZHi_w,Ind4Dnew) - FDw(IXHi,IYHi,IZLo_w,Ind4Dnew) )*Zgrid_w + FDw(IXHi,IYHi,IZLo_w,Ind4Dnew)
   Ixhyz = ( Iyhz - Iylz )*Ygrid + Iylz

   Ixyzn = ( Ixhyz - Ixlyz )*Xgrid + Ixlyz

!FIXME:delete
!   FD_GetWindSpeed%Velocity(3) = ( Ixyzn - Ixyzo )*Tgrid + Ixyzo
   TempWindSpeed(3) = ( Ixyzn - Ixyzo )*Tgrid + Ixyzo

   ! Copy the windspeed info to the output
   FD_GetWindSpeed = TempWindSpeed

   !-------------------------------------------------------------------------------------------------
   ! Set the previous time here to compare with later...
   !-------------------------------------------------------------------------------------------------
   PrevTime = Time

   RETURN

END FUNCTION FD_GetWindSpeed
!====================================================================================================
SUBROUTINE FD_Terminate( ErrStat )
! This subroutine deallocates arrays, closes files, and un-sets the initialization flag.
!----------------------------------------------------------------------------------------------------

   INTEGER,    INTENT(OUT)    :: ErrStat           ! return 0 if no errors; non-zero otherwise


   CLOSE( FDunit )

   ErrStat = 0

   IF ( ALLOCATED( FDu       ) )   DEALLOCATE( FDu,       STAT=ErrStat )
   IF ( ALLOCATED( FDv       ) )   DEALLOCATE( FDv,       STAT=ErrStat )
   IF ( ALLOCATED( FDw       ) )   DEALLOCATE( FDw,       STAT=ErrStat )
   IF ( ALLOCATED( FDuData   ) )   DEALLOCATE( FDuData,   STAT=ErrStat )
   IF ( ALLOCATED( FDvData   ) )   DEALLOCATE( FDvData,   STAT=ErrStat )
   IF ( ALLOCATED( FDwData   ) )   DEALLOCATE( FDwData,   STAT=ErrStat )
   IF ( ALLOCATED( Times4D   ) )   DEALLOCATE( Times4D,   STAT=ErrStat )
   IF ( ALLOCATED( Times4DIx ) )   DEALLOCATE( Times4DIx, STAT=ErrStat )
   IF ( ALLOCATED( AdvFiles  ) )   DEALLOCATE( AdvFiles,  STAT=ErrStat )

   Initialized = .FALSE.

END SUBROUTINE FD_Terminate
!====================================================================================================
END MODULE FDWind
