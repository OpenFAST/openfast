MODULE CTWind
! This module uses reads coherent turbulence parameter (CTP) files and processes the data in them
! to get coherent turbulence which is later superimposed on a background wind field (the super-
! positioning occurs elsewhere).  The turbulence in this module is part of the KH billow, which
! can be read using FDWind.  As a result, the scaling here should be similiar to FDWind.
!
! This module assumes that the origin, (0,0,0), is located at the tower centerline at ground level,
! and that all units are specified in the metric system (using meters and seconds).
! Data is shifted by half the grid width when used with FFWind.
!
!  Created 25-Sept-2009 by B. Jonkman, National Renewable Energy Laboratory
!     using subroutines and modules from AeroDyn v12.58
!
!  Modified Jan-2013 by A. Platt, National Renewable Energy Laboratory
!     to fit the modularization framework used by FAST
!
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


!fixme: this invokes SAVE
   INTEGER, PARAMETER           :: NumComps  = 3                              ! number of components

            ! CT_Wind
   REAL(ReKi)                   :: DelYCTgrid                                 ! The nondimensional distance between grid points in the y direction.
   REAL(ReKi)                   :: DelZCTgrid                                 ! The nondimensional distance between grid points in the z direction.
   REAL(ReKi)                   :: CTDistSc                                   ! Disturbance scale (ratio of wave height to rotor diameter).
   REAL(ReKi)                   :: CTOffset (NumComps)                        ! Offsets to convert integer data to actual wind speeds.
   REAL(ReKi)                   :: CTScale  (NumComps)                        ! Scaling factors to convert integer data to actual wind speeds.


   REAL(ReKi), ALLOCATABLE      :: CTvelU   (:,:,:)                           ! The y-z grid velocity data (U components) for the lower- and upper-bound time slices
   REAL(ReKi), ALLOCATABLE      :: CTvelV   (:,:,:)                           ! The y-z grid velocity data (V components) for the lower- and upper-bound time slices
   REAL(ReKi), ALLOCATABLE      :: CTvelW   (:,:,:)                           ! The y-z grid velocity data (W components) for the lower- and upper-bound time slices
   REAL(ReKi)                   :: CTLy                                       ! Fractional location of tower centerline from right (looking downwind) to left side of the dataset.
   REAL(ReKi)                   :: CTLz                                       ! Fractional location of hub height from bottom to top of dataset.
   REAL(ReKi)                   :: CTScaleVel                                 ! Scaling velocity, U0.  2*U0 is the difference in wind speed between the top and bottom of the wave.
   REAL(ReKi), ALLOCATABLE      :: Tdata    (:)                               ! The list of times for the CT-wind input files.

   REAL(ReKi)                   :: CT_Zref                                    ! The reference height for the CT file (the bottom of the billow)
   REAL(ReKi)                   :: CTYHWid                                    ! The half the width of the background dataset, used to compute the CTwind time offset
   REAL(ReKi)                   :: CTYmax                                     ! The dimensional lateral width of the dataset.
   REAL(ReKi)                   :: CTYt                                       ! Distance of the tower from the right side of the dataset (looking downwind).
   REAL(ReKi)                   :: CTZmax                                     ! The dimensional vertical height of the dataset.
   REAL(ReKi)                   :: InvMCTWS                                   ! The multiplicative inverse of the mean hub height wind speed for the CT wind data

   INTEGER                      :: CT_DF_Y                                    ! The decimation factor for the CT wind data in the y direction.
   INTEGER                      :: CT_DF_Z                                    ! The decimation factor for the CT wind data in the z direction.
   INTEGER                      :: CTvel_files(2)                             ! Times for the CT wind files stored in CTvel arrays.

   INTEGER                      :: IndCT_hi                                   ! An index into the 3rd dimension of the CTvel arrays, indicating the upper time slice (allows us to avoid copying array)
   INTEGER                      :: IndCT_lo                                   ! An index into the 3rd dimension of the CTvel arrays, indicating the lower time slice (allows us to avoid copying array)

   INTEGER                      :: NumCTt                                     ! The number of CT wind grids, no more than one grid per time step.
   INTEGER                      :: NumCTy                                     ! The number of CT wind grid points in the y direction.
   INTEGER                      :: NumCTyD                                    ! The decimated number of CT wind grid points in the y direction.
   INTEGER                      :: NumCTyD1                                   ! The decimated number of CT wind grid points in the y direction minus 1.
   INTEGER                      :: NumCTz                                     ! The number of CT wind grid points in the z direction.
   INTEGER                      :: NumCTzD                                    ! The decimated number of CT wind grid points in the z direction.
   INTEGER                      :: NumCTzD1                                   ! The decimated number of CT wind grid points in the z direction minus 1.
!FIXME: move to otherstate
   INTEGER, SAVE                :: TimeIndx  = 0                              ! Index into the time array
   INTEGER, ALLOCATABLE         :: TimeStpCT (:)                              ! The list of time steps from the original LE simulation, associated with the CT-wind times.

   INTEGER                      :: CTWindUnit                                 ! unit number used to read the wind files at each call to CT_GetWindSpeed()

   LOGICAL                      :: CTVertShft                                 ! Flag to indicate whether or not to shift the z values for the w component.

   CHARACTER(3)                 :: CText                                      ! The extension used for coherent turbulence data files. (usually "les" or "dns")
   CHARACTER(1024)              :: CTSpath                                    ! The path to the CT wind files.

!FIXME: move this to types -- parameters!
   TYPE :: CTWindFiles
      CHARACTER(1024)           :: CTfile                                   ! The name of the file containing the time-step history of the wind files.
      CHARACTER(1024)           :: CTbackgr                                   ! The name of the background wind data
   END TYPE CTWindFiles

!FIXME: this should be moved to parameters.
   TYPE, PUBLIC :: CT_Backgr
!   TYPE :: CT_Backgr
      CHARACTER(1024)           :: WindFile                                   ! The name of the background wind file
      INTEGER                   :: WindFileType                               ! The type of background wind file (currently only FF)
      LOGICAL                   :: CoherentStr                                ! If the coherent time step file is blank or doesn't exist, this is FALSE (use the background only)
   END TYPE CT_Backgr


   PUBLIC                       :: CT_Init
   PUBLIC                       :: CT_GetWindSpeed
   PUBLIC                       :: CT_SetRefVal
   PUBLIC                       :: CT_Terminate

CONTAINS
!====================================================================================================
SUBROUTINE CT_Init(UnWind, WindFile, BackGrndValues, ErrStat, ErrMsg)
!  This subroutine is called at the beginning of a simulation.  It reads the CTP file to obtain
!  the name of the CTS file, the path locating the binary KH files, and decimation factors.
!  It returns the background wind file and type; it also returns a flag that determines if CT wind
!  files are ACTUALLY to be used (e.g., if the CTS file is blank or there is one line of zero in the
!  CTS time array).
!----------------------------------------------------------------------------------------------------

      ! Passed Variables:

   INTEGER,          INTENT(IN)     :: UnWind                        ! unit number for reading wind files
   CHARACTER(*),     INTENT(IN)     :: WindFile                      ! Name of the CTP (.ctp) wind file
   TYPE(CT_Backgr),  INTENT(OUT)    :: BackGrndValues                ! output background values
   INTEGER,          INTENT(OUT)    :: ErrStat                       ! return ErrID_None if no errors
   CHARACTER(*),     INTENT(OUT)    :: ErrMsg                        ! message if there was an error

      ! Local Variables:

   TYPE(CTWindFiles)                :: CTP_files
   CHARACTER(3)                     :: CT_SC_ext                     ! extension of the scaling file
   LOGICAL                          :: EmptyFileStat                 ! temporary variable indicating the CT file was empty / non-existent

      ! Temporary error handling Variables

   INTEGER                          :: TmpErrStat                    ! Temporary Error Status
   CHARACTER(LEN(ErrMsg))           :: TmpErrMsg                     ! Temporary error message returned



   !-------------------------------------------------------------------------------------------------
   ! Initialize temporary variables
   !-------------------------------------------------------------------------------------------------

   TmpErrStat  = ErrID_None
   TmpErrMsg   = ''

   !-------------------------------------------------------------------------------------------------
   ! Check that the module hasn't already been initialized.
   !-------------------------------------------------------------------------------------------------

   IF ( TimeIndx /= 0 ) THEN
      ErrMsg   = ' CTWind has already been initialized.'
      ErrStat  = ErrID_Fatal
      RETURN
   ELSE
      ErrMsg   = ''
      ErrStat  = ErrID_None
   END IF


   !-------------------------------------------------------------------------------------------------
   ! Read the CTP file and set the background data info to be returned later
   !-------------------------------------------------------------------------------------------------

   CALL ReadCTP( UnWind, WindFile, CTP_files, TmpErrStat, TmpErrMsg )

      ! Errors occured?
   ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
   ErrStat  = MAX( ErrStat, TmpErrStat )
   IF ( ErrStat >= AbortErrLev )    RETURN


   BackGrndValues%WindFile     = CTP_files%CTbackgr
   BackGrndValues%WindFileType = FF_Wind             !bjj: perhaps we should check the wind type here


   !-------------------------------------------------------------------------------------------------
   ! Read the CT file to get the time step and file number arrays
   !-------------------------------------------------------------------------------------------------

   CALL ReadCT( UnWind, CTP_files%CTfile, CT_SC_ext, EmptyFileStat, TmpErrStat, TmpErrMsg )

      ! Errors check
   ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
   ErrStat  = MAX(ErrStat, TmpErrStat)
   IF ( ErrStat >= AbortErrLev ) RETURN


         ! If the CTP_files%CTfile was empty or non-existent, we continue on without it.

   IF ( EmptyFileStat ) THEN

         ! The file is missing, blank (or possibly incomplete), or has only 1 time step line (which
         ! is zero); Go on without the CT file, using just the background

      ErrMsg = TRIM(ErrMsg)//' '//'Coherent turbulence wind file will be turned off.'

      BackGrndValues%CoherentStr  = .FALSE.
      CALL CT_Terminate( TmpErrStat, TmpErrMsg )

         ! Error check
      ErrStat  = MAX( ErrStat, TmpErrStat )
      ErrMsg = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
      IF ( ErrStat >= AbortErrLev ) RETURN


      RETURN

   ELSEIF ( (ErrStat < AbortErrLev) .AND. NumCTt > 1) THEN
      BackGrndValues%CoherentStr  = .TRUE.

      !-------------------------------------------------------------------------------------------------
      ! Read file containing scaling for the binary large-eddy files
      !-------------------------------------------------------------------------------------------------
      CALL ReadCTScales( UnWind, TRIM( CTSpath )//'\Scales.'//TRIM( CT_SC_ext ), TmpErrStat, TmpErrMsg )

         ! Errors occured?
      ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
      ErrStat  = MAX(ErrStat, TmpErrStat)
      IF ( ErrStat >= AbortErrLev ) RETURN


      CTScale(:)  = CTScaleVel*CTScale(:)
      CTOffset(:) = CTScaleVel*CTOffset(:)

   END IF


   !-------------------------------------------------------------------------------------------------
   ! Set some values that don't change during the run
   !-------------------------------------------------------------------------------------------------

   CTYHWid        = 0.0                                                    ! This value is used to perform a time shift (the equivalent distance of FFYHWid [approx. rotor radius])
   CT_Zref        = -1.0                                                   ! This value needs to be set after the corresponding background turbulence has been read (or the CTS file should be changed)

   NumCTyD        = ( NumCTy + CT_DF_Y - 1 )/CT_DF_Y                       ! The decimated number of CT wind grid points in the y direction.
   NumCTzD        = ( NumCTz + CT_DF_Z - 1 )/CT_DF_Z                       ! The decimated number of CT wind grid points in the z direction.
   NumCTyD1       = NumCTyD - 1                                            ! The decimated number of CT wind grid points in the y direction minus 1.
   NumCTzD1       = NumCTzD - 1                                            ! The decimated number of CT wind grid points in the z direction minus 1.

   CTYt           = CTYmax*CTLy                                            ! Distance of the tower from the right side of the dataset (looking downwind).
!   CTZt           = CTZmax*CTLz                                            ! Distance of the hub from the bottom of the dataset.
   DelYCTgrid     = 1.0/NumCTyD1                                           ! The nondimensional distance between grid points in the y direction.
   DelZCTgrid     = 1.0/NumCTzD1                                           ! The nondimensional distance between grid points in the z direction.



   !-------------------------------------------------------------------------------------------------
   ! Allocate the wind array and initialize it
   !-------------------------------------------------------------------------------------------------

   IF (.NOT. ALLOCATED(CTvelU) ) THEN
      ALLOCATE ( CTvelU(NumCTyD,NumCTzD,2), STAT=TmpErrStat )

      IF ( TmpErrStat /= 0 )  THEN
         ErrMsg   = TRIM(ErrMsg)//' Error allocating memory for the CTvelU array.'
         ErrStat  = ErrID_Fatal
         RETURN
      END IF
   END IF

   IF (.NOT. ALLOCATED(CTvelV) ) THEN
!      CALL AllocAry( CTvelV, NumCTyD, NumCTzD, 2, 'CTvelV', ErrStat ) !AllRAry3 AllocAry
      ALLOCATE ( CTvelV(NumCTyD,NumCTzD,2), STAT=TmpErrStat )

      IF ( TmpErrStat /= 0 )  THEN
         ErrMsg   = TRIM(ErrMsg)//' Error allocating memory for the CTvelV array.'
         ErrStat  = ErrID_Fatal
         RETURN
      END IF
   END IF

   IF (.NOT. ALLOCATED(CTvelW) ) THEN
!      CALL AllocAry( CTvelW, NumCTyD, NumCTzD, 2, 'CTvelW', ErrStat ) !AllRAry3 AllocAry
      ALLOCATE ( CTvelW(NumCTyD,NumCTzD,2), STAT=TmpErrStat )

      IF ( TmpErrStat /= 0 )  THEN
         ErrMsg   = TRIM(ErrMsg)//' Error allocating memory for the CTvelW array.'
         ErrStat  = ErrID_Fatal
         RETURN
      END IF
   END IF

   CTvelU(:,:,:) = 0.0                                                    ! the original velocity data
   CTvelV(:,:,:) = 0.0                                                    ! the original velocity data
   CTvelW(:,:,:) = 0.0                                                    ! the original velocity data

   !-------------------------------------------------------------------------------------------------
   ! Initialize the arrays and set the initialization flag
   !-------------------------------------------------------------------------------------------------
   CTvel_files(:) = 0                                                      ! the name of the files currently in the CTvel array
   CTWindUnit     = UnWind                                                 ! This unit is needed to open the binary files at each step
   TimeIndx       = 1

   RETURN

END SUBROUTINE CT_Init
!====================================================================================================
SUBROUTINE CT_SetRefVal(Height, HWidth, ErrStat, ErrMsg)

   REAL(ReKi),    INTENT(IN)           :: Height                           ! a reference height (should be hub height)
   REAL(ReKi),    INTENT(IN), OPTIONAL :: HWidth                           ! a reference offset (should be half grid width [~rotor radius])
   INTEGER,       INTENT(OUT)          :: ErrStat                          ! returns ErrID_None if no errors, nonzero otherwise
   CHARACTER(*),  INTENT(OUT)          :: ErrMsg                           ! Message to return about the error


   !-------------------------------------------------------------------------------------------------
   ! Check that we've initialized everything first
   !-------------------------------------------------------------------------------------------------

   IF ( TimeIndx == 0 ) THEN
      ErrMsg   = ' Initialialize the CTWind module before calling its subroutines.'
      ErrStat  = ErrID_Fatal
      RETURN
   ELSE IF ( CT_Zref >= 0 ) THEN
      ErrMsg   = ' Cannot reset the CTWind reference height in the middle of a simulation.'
      ErrStat  = ErrID_Fatal
      RETURN
   ELSE
      ErrStat = ErrID_None
   END IF


   !-------------------------------------------------------------------------------------------------
   ! Set the grid shift using the half-width
   !-------------------------------------------------------------------------------------------------
   IF ( PRESENT( HWidth ) ) THEN
      CTYHWid = HWidth

      IF ( CTYHWid < 0 ) THEN
         ErrMsg   = TRIM(ErrMsg)//' Reference width in CTWind cannot be negative.'
         CTYHWid  = 0
         ErrStat  = ErrID_Fatal
      END IF
   END IF


   !-------------------------------------------------------------------------------------------------
   ! Set the reference height (bottom of the KH billow) using the input hub-height
   !-------------------------------------------------------------------------------------------------
      ! CTZt = CTZmax*CTLz             ! the distance between the hub and the bottom of the dataset

   CT_Zref = Height - CTZmax*CTLz      ! the height of the bottom of the KH billow

   IF ( CT_Zref < 0 ) THEN
      ErrMsg   = TRIM(ErrMsg)//' Reference height in CTWind cannot be negative.'
      CT_Zref  = 0
      ErrStat  = ErrID_Fatal
   END IF


END SUBROUTINE CT_SetRefVal
!====================================================================================================
FUNCTION CT_GetWindSpeed(Time, InputPosition, ErrStat, ErrMsg)
! This function receives time and position (in InputInfo) where (undisturbed) velocities are are
! requested.  It returns the velocities at the specified time and space that are superimposed on
! a background wind flow.  This function interpolates into the full-field CT wind arrays, performing
! a time shift based on the average windspeed. The modified time is used to decide which pair of time
! slices to interpolate within and between. After finding the two time slices, it decides which four
! grid points bound the (Y,Z) pair. It does a bilinear interpolation for (Y,Z) on each bounding time
! slice, then linearly interpolates between the 2 time slices. This routine assumes that X is downwind,
! Y is to the left when looking downwind and Z is up.  In the time (X) and Z directions, steady winds
! are used when extrapolation is required.  The dataset is assumed to be periodic in the Y direction.
!----------------------------------------------------------------------------------------------------

      ! Passed variables:

   REAL(DbKi),          INTENT(IN) :: Time                                 ! the time
   REAL(ReKi),          INTENT(IN) :: InputPosition(3)                     ! the position (X,Y,Z)
   INTEGER,             INTENT(OUT):: ErrStat                              ! returns ErrID_None if no error; non-zero otherwise
   CHARACTER(*),        INTENT(OUT) :: ErrMsg                              ! Message about the error
   REAL(ReKi)                    :: CT_GetWindSpeed                        ! the resultant wind speed


      ! Local Variables:

   REAL(ReKi)                    :: Iyz_th                                 ! Temporary interpolated value. (time hi, all y, all z)
   REAL(ReKi)                    :: Iyz_tl                                 ! Temporary interpolated value. (time lo, all y, all z)
   REAL(ReKi)                    :: Iyhz                                   ! Temporary interpolated value. (y hi, all z)
   REAL(ReKi)                    :: Iylz                                   ! Temporary interpolated value. (y lo, all z)
   REAL(ReKi)                    :: TimeShifted                            ! Shifted time (necessary because we're keeping x constant)
   REAL(ReKi)                    :: Tgrid                                  ! Fractional distance between time grids.
   REAL(ReKi)                    :: Ygrid                                  ! Fractional distance between grids in the y direction.
   REAL(ReKi)                    :: Ynorm                                  ! Nondimensional lateral distance of the analysis point from right side of dataset (looking downwind).
   REAL(ReKi)                    :: Zgrid(3)                               ! Fractional distance between grids in the z direction.
   REAL(ReKi)                    :: Znorm                                  ! Nondimensional vertical distance of the analysis point from bottom of dataset.

   INTEGER                       :: I
   INTEGER                       :: IYHi
   INTEGER                       :: IYLo
   INTEGER                       :: IZHi(3)
   INTEGER                       :: IZLo(3)

      ! Temporary error handling
   INTEGER                       :: TmpErrStat
   CHARACTER(LEN(ErrMsg))        :: TmpErrMsg


   !-------------------------------------------------------------------------------------------------
   ! Check that we've initialized everything first
   !-------------------------------------------------------------------------------------------------

   IF ( TimeIndx == 0 ) THEN
      ErrMsg   = ' Initialialize the CTWind module before calling its subroutines.'
      ErrStat  = ErrID_Fatal
      RETURN
   ELSE IF ( CT_Zref < 0 ) THEN
      ErrMsg   = ' Set the reference height in the CTWind module before calling CT_GetWindSpeed.'
      ErrStat  = ErrID_Fatal
      RETURN
   ELSE
      ErrStat = ErrID_None
   END IF


   !-------------------------------------------------------------------------------------------------
   ! Perform the time shift. At time=0, a point half the grid width downstream will index into the zero
   ! time slice.  CTYHWid is used to shift the CT wind the same as FF wind is shifted.
   ! This assumes that the coherent turbulence events are moving at MCTWS
   !-------------------------------------------------------------------------------------------------

   TimeShifted = TIME + ( CTYHWid - InputPosition(1) )*InvMCTWS


   !-------------------------------------------------------------------------------------------------
   ! Find the bounding time slices:
   ! Linearly interpolate in time (or set to 0 before and/or after)
   ! (compare with NWTC_Num.f90\InterpStpReal)
   !-------------------------------------------------------------------------------------------------

      ! Let's check the limits first.

   IF ( TimeShifted <= Tdata(1) )  THEN

      TimeIndx = 1
      Tgrid    = 0.0

!      CT_GetWindSpeed%Velocity(:) = 0.0
!      RETURN

   ELSE IF ( TimeShifted >= Tdata(NumCTt) )  THEN

      TimeIndx = NumCTt - 1
      Tgrid    = 1.0

!      CT_GetWindSpeed%Velocity(:) = 0.0
!      RETURN

   ELSE

         ! Let's interpolate!

      TimeIndx = MAX( MIN( TimeIndx, NumCTt-1 ), 1 )


      DO

         IF ( TimeShifted < Tdata(TimeIndx) )  THEN

            TimeIndx = TimeIndx - 1

         ELSE IF ( TimeShifted >= Tdata(TimeIndx+1) )  THEN

            TimeIndx = TimeIndx + 1

         ELSE

            Tgrid = MIN( MAX( ( TimeShifted - Tdata(TimeIndx) )/( Tdata(TimeIndx+1) - Tdata(TimeIndx) ), 0.0 ), 1.0 )
            EXIT

         END IF

      END DO

   END IF


   !-------------------------------------------------------------------------------------------------
   ! Read the data at the two time steps, if necessary
   !-------------------------------------------------------------------------------------------------

   IF ( TimeStpCT(TimeIndx) == CTvel_files(2) ) THEN
      IndCT_lo = 2
      IndCT_hi = 1

   ELSE
      IndCT_lo = 1
      IndCT_hi = 2

      IF ( TimeStpCT(TimeIndx) /= CTvel_files(IndCT_lo) ) THEN
         CTvel_files(IndCT_lo) = TimeStpCT(TimeIndx)
         CALL ReadCTData ( CTWindUnit, CTvel_files(IndCT_lo), IndCT_lo, TmpErrStat, TmpErrMsg )

            ! Errors occured?
         ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
         ErrStat  = MAX(ErrStat, TmpErrStat)
         IF ( ErrStat >= AbortErrLev ) RETURN
      END IF

   END IF


   IF ( CTvel_files(IndCT_hi) /= TimeStpCT(TimeIndx+1) ) THEN

      CTvel_files(IndCT_hi) = TimeStpCT(TimeIndx+1)
      CALL ReadCTData ( CTWindUnit, CTvel_files(IndCT_hi), IndCT_hi, TmpErrStat, TmpErrMsg  )

         ! Errors occured?
      ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
      ErrStat  = MAX(ErrStat, TmpErrStat)
      IF ( ErrStat >= AbortErrLev ) RETURN
   END IF


   !-------------------------------------------------------------------------------------------------
   ! Calculate the y values;                   The lower-right corner is (1,1) when looking downwind.
   ! note that the KH data is periodic in this direction
   !-------------------------------------------------------------------------------------------------

   Ynorm = ( CTYt + InputPosition(2) )/CTYmax

      ! Ensure Ynorm is not negative.  The wave is periodic in y.

   IF ( Ynorm < 0.0 ) THEN
      Ynorm = 1.0 + MOD(Ynorm, 1.0)
   ENDIF

   Ygrid = MIN( MAX( MOD( Ynorm, DelYCTgrid ), 0.0 ), 1.0 )
   IYLo  = MAX( MOD( INT( Ynorm*NumCTyD1 ) + 1, NumCTyD1 ), 1 )
   IYHi =  MOD( IYLo, NumCTyD ) + 1


   !-------------------------------------------------------------------------------------------------
   ! Calculate the z values                   The lower-right corner is (1,1) when looking downwind.
   ! Note: the equivalent Znorm for the w-component may be shifted vertically by half the original
   ! grid spacing. (the K-H data staggers w differently than u & v).  We store IZLo, IZHi, and
   ! Zgrid in an array to account for this difference.
   !-------------------------------------------------------------------------------------------------

   Znorm = MIN( MAX( ( InputPosition(3) - CT_Zref )/CTZmax, 0.0 ), 1.0 ) ! non-dimensional height (CT_Zref is the bottom of the billow)

      ! Find out fractionally how far we are between grids in time and between grid points in each direction.
      ! Limit values to avoid extrapolation.  We need this for interpolation later on.

   Zgrid(1:2)   = MIN( MAX( MOD( Znorm, DelZCTgrid ), 0.0 ), 1.0 )
   IZLo(1:2)    = MAX( INT( Znorm*NumCTzD1 ) + 1, 1 )            ! Make sure the lowest possible value is 1.

      ! If we are located at the upper end of the z dimension, decrement the index by one and set the grid coordinate to 1.

   IF ( IZLo(1) == NumCTzD )  THEN
      IZLo(1:2)  = NumCTzD1
      Zgrid(1:2) = 1.0
   ENDIF


   !-------------------------------------------------------------------------------------------------
   ! Find the equivalent Znorm for the w-component, which may be shifted vertically by half
   ! the original grid spacing. (This is necessary due to the fact that the K-H data staggers w
   ! differently than u & v).  LES and DNS scale differently.
   !-------------------------------------------------------------------------------------------------

   IF ( CTVertShft )  THEN
      Znorm = MAX( Znorm - 0.5*DelZCTgrid/CT_DF_Z, 0.0 )

      Zgrid(3) = MIN( MAX( MOD( Znorm, DelZCTgrid ), 0.0 ), 1.0 )
      IZLo(3)  = MAX( INT( Znorm*NumCTzD1 ) + 1, 1 )       ! Make sure the lowest possible value is 1.


         ! If we are located at the upper end of the z dimension, decrement the index by one and set the grid coordinate to 1.

      IF ( IZLo(3) == NumCTzD )  THEN
         IZLo(3)  = NumCTzD1
         Zgrid(3) = 1.0
      ENDIF

   ELSE
      IZLo(3) = IZLo(1)
      Zgrid(3)= Zgrid(1)
   ENDIF

   IZHi(:) = IZLo(:) + 1

!bjj: old versions used Zgrid(3) = Zgrid(1) without regard to CTVertShft. It seemed wrong to me so I changed it.

   !-------------------------------------------------------------------------------------------------
   !  Interpolate for U component of wind within the grid.
   !-------------------------------------------------------------------------------------------------

      I = 1
         ! linearaly interpolate in the lower time slice
      Iylz   = ( CTvelU(IYLo,IZHi(I),IndCT_lo) - CTvelU(IYLo,IZLo(I),IndCT_lo) )*Zgrid(I) + CTvelU(IYLo,IZLo(I),IndCT_lo)
      Iyhz   = ( CTvelU(IYHi,IZHi(I),IndCT_lo) - CTvelU(IYHi,IZLo(I),IndCT_lo) )*Zgrid(I) + CTvelU(IYHi,IZLo(I),IndCT_lo)
      Iyz_tl = ( Iyhz - Iylz )*Ygrid + Iylz

         ! linearaly interpolate in the upper time slice
      Iylz   = ( CTvelU(IYLo,IZHi(I),IndCT_hi) - CTvelU(IYLo,IZLo(I),IndCT_hi) )*Zgrid(I) + CTvelU(IYLo,IZLo(I),IndCT_hi)
      Iyhz   = ( CTvelU(IYHi,IZHi(I),IndCT_hi) - CTvelU(IYHi,IZLo(I),IndCT_hi) )*Zgrid(I) + CTvelU(IYHi,IZLo(I),IndCT_hi)
      Iyz_th = ( Iyhz - Iylz )*Ygrid + Iylz

      CT_GetWindSpeed = ( Iyz_th - Iyz_tl )*Tgrid + Iyz_tl


   !-------------------------------------------------------------------------------------------------
   !  Interpolate for V component of wind within the grid.
   !-------------------------------------------------------------------------------------------------

      I = 2

         ! linearaly interpolate in the lower time slice
      Iylz   = ( CTvelV(IYLo,IZHi(I),IndCT_lo) - CTvelV(IYLo,IZLo(I),IndCT_lo) )*Zgrid(I) + CTvelV(IYLo,IZLo(I),IndCT_lo)
      Iyhz   = ( CTvelV(IYHi,IZHi(I),IndCT_lo) - CTvelV(IYHi,IZLo(I),IndCT_lo) )*Zgrid(I) + CTvelV(IYHi,IZLo(I),IndCT_lo)
      Iyz_tl = ( Iyhz - Iylz )*Ygrid + Iylz

         ! linearaly interpolate in the upper time slice
      Iylz   = ( CTvelV(IYLo,IZHi(I),IndCT_hi) - CTvelV(IYLo,IZLo(I),IndCT_hi) )*Zgrid(I) + CTvelV(IYLo,IZLo(I),IndCT_hi)
      Iyhz   = ( CTvelV(IYHi,IZHi(I),IndCT_hi) - CTvelV(IYHi,IZLo(I),IndCT_hi) )*Zgrid(I) + CTvelV(IYHi,IZLo(I),IndCT_hi)
      Iyz_th = ( Iyhz - Iylz )*Ygrid + Iylz

      CT_GetWindSpeed = ( Iyz_th - Iyz_tl )*Tgrid + Iyz_tl


   !-------------------------------------------------------------------------------------------------
   !  Interpolate for W component of wind within the grid.
   !-------------------------------------------------------------------------------------------------

      I = 3

         ! linearaly interpolate in the lower time slice
      Iylz   = ( CTvelW(IYLo,IZHi(I),IndCT_lo) - CTvelW(IYLo,IZLo(I),IndCT_lo) )*Zgrid(I) + CTvelW(IYLo,IZLo(I),IndCT_lo)
      Iyhz   = ( CTvelW(IYHi,IZHi(I),IndCT_lo) - CTvelW(IYHi,IZLo(I),IndCT_lo) )*Zgrid(I) + CTvelW(IYHi,IZLo(I),IndCT_lo)
      Iyz_tl = ( Iyhz - Iylz )*Ygrid + Iylz

         ! linearaly interpolate in the upper time slice
      Iylz   = ( CTvelW(IYLo,IZHi(I),IndCT_hi) - CTvelW(IYLo,IZLo(I),IndCT_hi) )*Zgrid(I) + CTvelW(IYLo,IZLo(I),IndCT_hi)
      Iyhz   = ( CTvelW(IYHi,IZHi(I),IndCT_hi) - CTvelW(IYHi,IZLo(I),IndCT_hi) )*Zgrid(I) + CTvelW(IYHi,IZLo(I),IndCT_hi)
      Iyz_th = ( Iyhz - Iylz )*Ygrid + Iylz

      CT_GetWindSpeed = ( Iyz_th - Iyz_tl )*Tgrid + Iyz_tl


   RETURN

END FUNCTION CT_GetWindSpeed
!====================================================================================================
SUBROUTINE ReadCTData ( UnWind, CTFileNo, Itime, ErrStat, ErrMsg )
!    This subroutine is used to read one time-step's worth of large-eddy
!    zero-mean wind data for each wind component from a file.
!----------------------------------------------------------------------------------------------------


      ! Passed variables.

   INTEGER,       INTENT(IN)     :: UnWind                                    ! The I/O unit of the input file
   INTEGER,       INTENT(IN)     :: CTFileNo                                  ! The number of the file to read
   INTEGER,       INTENT(IN)     :: Itime                                     ! The index of the time slice
   INTEGER,       INTENT(OUT)    :: ErrStat                                   ! returns ErrID_None if no error; non-zero otherwise
   CHARACTER(*),  INTENT(OUT)    :: ErrMsg                                    ! Message to return about the error

      ! Local variables.

!   CHARACTER(1),PARAMETER        :: Comp(NumComps) = (/'u', 'v', 'w' /)       ! the wind components
   CHARACTER(5)                  :: CTnum                                     ! string equivalent of input variable CTFileNo
   CHARACTER(1024)               :: FileName                                  ! The name of the input data file

      ! Temporary error handling Variables

   INTEGER                          :: TmpErrStat                             ! Temporary Error Status
   CHARACTER(LEN(ErrMsg))           :: TmpErrMsg                              ! Temporary error message returned


      !------------------------------------------------------------------------
      ! Initialize the error handling
      !------------------------------------------------------------------------
   ErrStat  = ErrID_None
   ErrMsg   = ''


   IF ( CTFileNo == 0 ) THEN

      CTvelU(:,:,Itime) = 0.0
      CTvelV(:,:,Itime) = 0.0
      CTvelW(:,:,Itime) = 0.0

   ELSE
         ! Loop through the components

      WRITE( CTnum, '(I5.5)' ) CTFileNo


      FileName = TRIM( CTSpath )//'\u\u_16i_'//CTnum//'.'//TRIM( CText )
      CALL LoadCTData( UnWind, TRIM(FileName), Itime, 1, CTvelU, TmpErrStat, TmpErrMsg )

         ! Errors occured?
      ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
      ErrStat  = MAX(ErrStat, TmpErrStat)
      IF ( ErrStat >= AbortErrLev ) RETURN


      FileName = TRIM( CTSpath )//'\v\v_16i_'//CTnum//'.'//TRIM( CText )
      CALL LoadCTData( UnWind, TRIM(FileName), Itime, 2, CTvelV, TmpErrStat, TmpErrMsg )

         ! Errors occured?
      ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
      ErrStat  = MAX(ErrStat, TmpErrStat)
      IF ( ErrStat >= AbortErrLev ) RETURN


      FileName = TRIM( CTSpath )//'\w\w_16i_'//CTnum//'.'//TRIM( CText )
      CALL LoadCTData( UnWind, TRIM(FileName), Itime, 3, CTvelW, TmpErrStat, TmpErrMsg )

         ! Errors occured?
      ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
      ErrStat  = MAX(ErrStat, TmpErrStat)
      IF ( ErrStat >= AbortErrLev ) RETURN


   END IF

   RETURN

END SUBROUTINE ReadCTData
!====================================================================================================
SUBROUTINE LoadCTData( UnWind, FileName, ITime, IComp, Vel, ErrStat, ErrMsg )
!  This function is used to read the input parameters for the coherent turbulence events,
!  based on the large-eddy simulation.
!----------------------------------------------------------------------------------------------------

      ! Passed variables.

   INTEGER,       INTENT(IN)     :: UnWind                                    ! The I/O unit of the input file
   CHARACTER(*),  INTENT(IN)     :: FileName                                  ! The name of the file to open
   INTEGER,       INTENT(IN)     :: Itime                                     ! The index of the time slice
   INTEGER,       INTENT(IN)     :: IComp                                     ! The index of the component
   REAL(ReKi),    INTENT(INOUT)  :: Vel    (NumCTyD,NumCTzD,2)                ! returns the velocity array (don't use INTENT OUT!)
   INTEGER,       INTENT(OUT)    :: ErrStat                                   ! returns ErrID_None if no error; non-zero otherwise
   CHARACTER(*),  INTENT(OUT)    :: ErrMsg                                    ! A message abouth the error


      ! Local Variables

   INTEGER(B2Ki)                 :: Com    (NumCTy)                           ! Temporary array to hold component's integer values for a given Z.
   INTEGER                       :: IY                                        ! A DO index for indexing the arrays in the y direction.
   INTEGER                       :: IYK                                       ! An index for the decimated arrays in the y direction.
   INTEGER                       :: IZ                                        ! A DO index for indexing the arrays in the z direction.
   INTEGER                       :: IZK                                       ! An index for the decimated arrays in the z direction.

      ! Temporary error handling
   INTEGER                       :: TmpErrStat
   CHARACTER(LEN(ErrMsg))        :: TmpErrMsg


   !-------------------------------------------------------------------------------------------------
   ! Set temporary Error info
   !-------------------------------------------------------------------------------------------------

   ErrStat  = ErrID_None
   ErrMsg   = ''
!FIXME: remove with add ErrMsg
   TmpErrMsg   = ''


   !-------------------------------------------------------------------------------------------------
   ! Open the input file
   !-------------------------------------------------------------------------------------------------

   CALL OpenUInBEFile( UnWind, TRIM(FileName), 2*NumCTy, TmpErrStat )   ! add ErrMsg
   IF ( ErrStat >= AbortErrLev ) THEN
      ErrMsg = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
      RETURN
   ELSE
      ErrMsg = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
   ENDIF


   !-------------------------------------------------------------------------------------------------
   ! Read the data and fill the arrays
   !-------------------------------------------------------------------------------------------------

   IZK = 0                          ! the Z index into the array (necessary b/c of decimation factor)
   DO IZ=1,NumCTz,CT_DF_Z

      READ (UnWind,REC=IZ,IOSTAT=TmpErrStat)  Com

      IF ( TmpErrStat /= 0 )  THEN

         ErrMsg = TRIM(ErrMsg)//' Error reading record '//TRIM( Num2LStr( IZ ) )//' of the binary CT wind file, "' &
                           //TRIM( FileName )//'."'
         ErrStat=ErrID_Fatal
         RETURN

      ENDIF

      IZK = IZK + 1
      IYK = 0                       ! the Y index into the array (necessary b/c of decimation factor)

      DO IY=1,NumCTy,CT_DF_Y
         IYK = IYK + 1
         Vel(IYK,IZK,ITime) = CTScale(IComp)*Com(IY) + CTOffset(IComp)
      ENDDO ! IY

   ENDDO ! IZ


   !-------------------------------------------------------------------------------------------------
   ! Close the file
   !-------------------------------------------------------------------------------------------------
      CLOSE ( UnWind )

   RETURN


END SUBROUTINE LoadCTData
!====================================================================================================
SUBROUTINE ReadCTP( UnWind, FileName, CTPscaling, ErrStat, ErrMsg )
!  This function is used to read the input parameters for the coherent turbulence events,
!  based on the large-eddy simulation.
!----------------------------------------------------------------------------------------------------


      ! Passed variables.

   INTEGER,             INTENT(IN)  :: UnWind                                 ! The I/O unit of the input file
   CHARACTER(*),        INTENT(IN)  :: FileName                               ! The name of the input data file
   TYPE(CTWindFiles),   INTENT(OUT) :: CTPscaling                             ! The file names contained in the CTP file
   INTEGER,             INTENT(OUT) :: ErrStat                                ! returns 0 if no error; non-zero otherwise'
   CHARACTER(*),        INTENT(OUT) :: ErrMsg                                 ! Error message to return


      ! Local variables.

   CHARACTER(1024)                  :: HeaderLine                             ! The header text in the file
   CHARACTER(1024)                  :: TmpPath

      ! Temporary error handling Variables

   INTEGER                          :: TmpErrStat                             ! Temporary Error Status
   CHARACTER(LEN(ErrMsg))           :: TmpErrMsg                              ! Temporary error message returned

   !-------------------------------------------------------------------------------------------------
   ! Set temporary Error info
   !-------------------------------------------------------------------------------------------------

   ErrStat  = ErrID_None
   ErrMsg   = ''


   !-------------------------------------------------------------------------------------------------
   ! Open the CTP input file
   !-------------------------------------------------------------------------------------------------

   CALL OpenFInpFile ( UnWind, TRIM( FileName ), TmpErrStat )  ! add ErrMsg


      ! Errors occured?
   ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
   ErrStat  = MAX(ErrStat, TmpErrStat)
   IF ( ErrStat >= AbortErrLev ) RETURN


   !-------------------------------------------------------------------------------------------------
   ! Read the CTP input file
   !-------------------------------------------------------------------------------------------------
   CALL ReadStr( UnWind, TRIM( FileName ), HeaderLine, 'Header line', 'The header line in the CTP file', TmpErrStat )  ! add ErrMsg

      ! Errors occured?
   ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
   ErrStat  = MAX(ErrStat, TmpErrStat)
   IF ( ErrStat > AbortErrLev ) THEN
      ErrMsg   = TRIM(ErrMsg)//' Error reading header line of '//TRIM(FileName)//'.'
      RETURN
   ENDIF
   CALL WrScr ( ' Heading of the CT-wind-parameter file: "'//TRIM(HeaderLine)//'"' )


   CALL ReadCom( UnWind, TRIM( FileName ), 'parameter header line', TmpErrStat )  ! add ErrMsg

      ! Errors occured?
   ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
   ErrStat  = MAX( ErrStat, TmpErrStat )
   IF ( ErrStat >= AbortErrLev ) THEN
      ErrMsg   = TRIM(ErrMsg)//' Error reading parameter header line of '//TRIM(FileName)//'.'
      RETURN
   ENDIF


   CALL ReadVar( UnWind, TRIM( FileName ), CTSpath,  'CTSpath',  &
                  'Location (path) of the binary coherent turbulence dataset', TmpErrStat )  ! add ErrMsg

      ! Errors occured?
   ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
   ErrStat  = MAX( ErrStat, TmpErrStat )
   IF ( ErrStat >= AbortErrLev ) THEN
      ErrMsg   = TRIM(ErrMsg)//' Error reading location (path) of the binary coherent turbulence dataset '//TRIM(FileName)//'.'
      RETURN
   ENDIF


   CALL ReadVar( UnWind, TRIM( FileName ), CTPscaling%CTfile, 'CTfile', &
                  'File containing the time steps for the coherent turbulence events (.cts)', TmpErrStat )  ! add ErrMsg

      ! Errors occured?
   ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
   ErrStat  = MAX( ErrStat, TmpErrStat )
   IF ( ErrStat >= AbortErrLev ) THEN
      ErrMsg   = TRIM(ErrMsg)//' Error reading name of file containing time steps for the coherent turbulent events from ' &
               //TRIM(FileName)//'.'
      RETURN
   ENDIF


   IF ( PathIsRelative( CTPscaling%CTfile ) ) THEN
      CALL GetPath( FileName, TmpPath )
      CTPscaling%CTfile = TRIM(TmpPath)//TRIM(CTPscaling%CTfile)
   END IF


   CALL ReadVar( UnWind, TRIM( FileName ), CTPscaling%CTbackgr, 'CTbackgr', 'File containing the background wind', TmpErrStat )  ! add ErrMsig

      ! Errors occured?
   ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
   ErrStat  = MAX( ErrStat, TmpErrStat )
   IF ( ErrStat >= AbortErrLev ) THEN
      ErrMsg   = TRIM(ErrMsg)//' Error reading name of file containing the background wind from '//TRIM(FileName)//'.'
      RETURN
   ENDIF

   IF ( PathIsRelative( CTPscaling%CTbackgr ) ) THEN
      CALL GetPath( FileName, TmpPath )
      CTPscaling%CTbackgr = TRIM(TmpPath)//TRIM(CTPscaling%CTbackgr)
   END IF


   CALL ReadVar( UnWind, TRIM( FileName ), CT_DF_Y, 'CT_DF_Y', 'Decimation factor for wind data in the Y direction', TmpErrStat )  ! add ErrMsg

      ! Errors occured?
   ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
   ErrStat  = MAX( ErrStat, TmpErrStat )
   IF ( ErrStat >= AbortErrLev ) THEN
      ErrMsg   = TRIM(ErrMsg)//' Error reading the decimation factor for wind data in the Y direction from ' &
               //TRIM(FileName)//'.'
      RETURN
   ENDIF


   CALL ReadVar( UnWind, TRIM( FileName ), CT_DF_Z, 'CT_DF_Z', 'Decimation factor for wind data in the Z direction', TmpErrStat )  ! add ErrMsg

      ! Errors occured?
   ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
   ErrStat  = MAX( ErrStat, TmpErrStat )
   IF ( ErrStat >= AbortErrLev ) THEN
      ErrMsg   = TRIM(ErrMsg)//' Error reading the decimation factor for wind data in the Z direction from ' &
               //TRIM(FileName)//'.'
      RETURN
   ENDIF


   !-------------------------------------------------------------------------------------------------
   ! Close the CTP input file
   !-------------------------------------------------------------------------------------------------

   CLOSE( UnWind )


END SUBROUTINE ReadCTP
!====================================================================================================
SUBROUTINE ReadCT ( UnWind, FileName, CT_SC_ext, EmptyFileStat, ErrStat, ErrMsg )
!  This subroutine is used to read the input parameters calculated in TurbSim for the scaling of
!  coherent turbulence events.  It reads the .cts file and saves the time step and file number arrays.
!----------------------------------------------------------------------------------------------------

! EmptyFileStat:   The CT file may not be exist or may be empty. TurbSim will do this in certain conditions.
!                 This is not a problem for program execution and used to be handled by Errstat<0. Set a
!                 warning in this case.


      ! Passed variables.

   INTEGER,             INTENT(IN)  :: UnWind                                 ! The I/O unit of the input file
   CHARACTER(*),        INTENT(IN)  :: FileName                               ! The name of the input data file
   CHARACTER(3),        INTENT(OUT) :: CT_SC_ext                              ! The extension used for coherent turbulence scale files.(usually "les", "dns", or "dat")
   LOGICAL,             INTENT(OUT) :: EmptyFileStat                           ! Special case for this file type. see note above
   INTEGER,             INTENT(OUT) :: ErrStat                                ! returns ErrID_Warn if can't open the file
   CHARACTER(*),        INTENT(OUT) :: ErrMsg                                 ! Message about what happened

      ! Local variables

   INTEGER                          :: IT                                     ! Loop counter

      ! Temporary error handling variables

   INTEGER                          :: TmpErrStat                             ! temporary ErrStat
   CHARACTER(LEN(ErrMsg))           :: TmpErrMsg                              ! temporary returned error message

   !-------------------------------------------------------------------------------------------------
   ! Initialize variables
   !-------------------------------------------------------------------------------------------------

   NumCTt = 0
   ErrMsg = ''

   !-------------------------------------------------------------------------------------------------
   ! Open the CTS input file -- can proceed if we can't open this file.
   !-------------------------------------------------------------------------------------------------

   CALL OpenFInpFile ( UnWind, TRIM( FileName ), TmpErrStat )  ! add ErrMsg when available
   IF (TmpErrStat /= 0) THEN
      EmptyFileStat = .TRUE.
      ErrMsg   = ' Error opening '//TRIM(FileName)//', ignoring it.'
      ErrStat  = ErrID_Warn
      RETURN
   ELSE
      EmptyFileStat = .FALSE.
   ENDIF

   !-------------------------------------------------------------------------------------------------
   ! Read the header of the CTS input file
   !-------------------------------------------------------------------------------------------------

      ! Check to see if the first value is numeric (old) or the file type (new) and start again

   READ  ( UnWind, *, IOSTAT=TmpErrStat ) CTScaleVel
   REWIND( UnWind )

   IF ( TmpErrStat /= 0 )  THEN   ! try again

      CALL ReadVar( UnWind, TRIM( FileName ), CText, 'CText', 'FileType ', TmpErrStat ) ! add ErrMsg

         ! Errors occured? Can proceeed if the file is empty.
      IF ( TmpErrStat < 0 ) THEN    ! end of record / end of file condition
         ErrMsg   = TRIM(ErrMsg)//' File '//TRIM(FileName)//' is empty; ignoring it'
         ErrStat  = ErrID_Warn
         EmptyFileStat = .TRUE.
         RETURN
      ELSE  ! positive, so something bad happened.
         ErrMsg   = TRIM(ErrMsg)//' Error reading from file '//TRIM(TmpErrMsg)
         ErrStat  = ErrID_Fatal
      ENDIF
      IF ( ErrStat >= AbortErrLev )    RETURN

      CT_SC_ext = CText

      CALL ReadVar( UnWind, TRIM( FileName ), CTScaleVel, 'CTScaleVel', ' ', TmpErrStat ) ! add ErrMsg

         ! Errors occured?
      ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
      ErrStat  = MAX( ErrStat, TmpErrStat )
      IF ( ErrStat >= AbortErrLev )    RETURN

   ELSE  ! assume LES files

      CALL ReadVar( UnWind, TRIM( FileName ), CTScaleVel, 'CTScaleVel', ' ', TmpErrStat ) ! add ErrMsg

         ! Errors occured?
      ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
      ErrStat  = MAX( ErrStat, TmpErrStat )
      IF ( ErrStat >= AbortErrLev )    RETURN

      CText     = 'les'
      CT_SC_ext = 'dat'
   END IF

   CALL ReadVar( UnWind, TRIM( FileName ), InvMCTWS, 'MeanCTWS', ' ', TmpErrStat )  ! add ErrMsg

      ! Errors occured?
   ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
   ErrStat  = MAX( ErrStat, TmpErrStat )
   IF ( ErrStat >= AbortErrLev )    RETURN

   InvMCTWS = 1.0 / InvMCTWS


   CALL ReadVar( UnWind, TRIM( FileName ), CTYmax, 'CTYmax', ' ', TmpErrStat )   ! add ErrMsg

      ! Errors occured?
   ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
   ErrStat  = MAX( ErrStat, TmpErrStat )
   IF ( ErrStat >= AbortErrLev )    RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), CTZmax, 'CTZmax', ' ', TmpErrStat )  ! add ErrMsg

      ! Errors occured?
   ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
   ErrStat  = MAX( ErrStat, TmpErrStat )
   IF ( ErrStat >= AbortErrLev )    RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), CTDistSc, 'CTDistSc', ' ', TmpErrStat )  ! add ErrMsg

      ! Errors occured?
   ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
   ErrStat  = MAX( ErrStat, TmpErrStat )
   IF ( ErrStat >= AbortErrLev )    RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), CTLy, 'CTLy', ' ', TmpErrStat )  ! add ErrMsg

      ! Errors occured?
   ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
   ErrStat  = MAX( ErrStat, TmpErrStat )
   IF ( ErrStat >= AbortErrLev )    RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), CTLz, 'CTLz', ' ', TmpErrStat )  ! add ErrMsg

      ! Errors occured?
   ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
   ErrStat  = MAX( ErrStat, TmpErrStat )
   IF ( ErrStat >= AbortErrLev )    RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), NumCTt, 'NumCTt', ' ', TmpErrStat )  ! add ErrMsg

      ! Errors occured?
   ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
   ErrStat  = MAX( ErrStat, TmpErrStat )
   IF ( ErrStat >= AbortErrLev )    RETURN


   !-------------------------------------------------------------------------------------------------
   ! Allocate space for the arrays
   !-------------------------------------------------------------------------------------------------

   IF (.NOT. ALLOCATED(Tdata) ) THEN
      ALLOCATE ( Tdata(NumCTt) , STAT=TmpErrStat )

         ! Errors occured?
      IF (TmpErrStat /= 0) THEN
         ErrMsg   = TRIM(ErrMsg)//' Error allocating memory for the Tdata array.'
         ErrStat  = ErrID_Fatal
         RETURN
      ENDIF

   END IF

   IF (.NOT. ALLOCATED(TimeStpCT) ) THEN
      ALLOCATE ( TimeStpCT(NumCTt) , STAT=TmpErrStat )

         ! Errors occured?
      IF (TmpErrStat /= 0) THEN
         ErrMsg   = TRIM(ErrMsg)//' Error allocating memory for the TimeStpCT array.'
         ErrStat  = ErrID_Fatal
         RETURN
      ENDIF

   END IF


   !-------------------------------------------------------------------------------------------------
   ! Read the arrays from the CTS input file
   !-------------------------------------------------------------------------------------------------

   DO IT=1,NumCTt

      READ (UnWind,*,IOSTAT=TmpErrStat)  Tdata(IT), TimeStpCT(IT)

         ! Errors occured?
      IF ( TmpErrStat /=0 ) THEN
         ErrStat=ErrID_Fatal
         ErrMsg   = TRIM(ErrMsg)//' Error reading record '//TRIM( Num2LStr( IT ) )//' of the CT-wind time-steps file, "' &
                         //TRIM( FileName )//'."'

         NumCTt = IT - 1
         RETURN
      ENDIF

   ENDDO ! IT


   !-------------------------------------------------------------------------------------------------
   ! Close the CTS input file
   !-------------------------------------------------------------------------------------------------
   CLOSE( UnWind )


   RETURN

END SUBROUTINE ReadCT
!====================================================================================================
SUBROUTINE ReadCTScales ( UnWind, FileName, ErrStat, ErrMsg )
!  This subroutine is used to read the input parameters for the coherent turbulence events, based
!  on the large-eddy simulation.
!----------------------------------------------------------------------------------------------------

      ! Passed variables

   INTEGER,             INTENT(IN)  :: UnWind                                 ! The I/O unit of the input file
   CHARACTER(*),        INTENT(IN)  :: FileName                               ! The name of the input data file
   INTEGER,             INTENT(OUT) :: ErrStat                                ! returns ErrID_None if no error; non-zero otherwise
   CHARACTER(*),        INTENT(OUT) :: ErrMsg                                 ! Message about the error


      ! Local variables

   INTEGER                          :: I                                      ! Array counter


      ! Temporary error handling Variables

   INTEGER                          :: TmpErrStat                             ! Temporary Error Status
   CHARACTER(LEN(ErrMsg))           :: TmpErrMsg                              ! Temporary error message returned

   !-------------------------------------------------------------------------------------------------
   ! Open the file with the scales (les or dns)
   !-------------------------------------------------------------------------------------------------

   CALL OpenFInpFile ( UnWind, TRIM( FileName ), ErrStat)   ! add ErrMsg
   IF (ErrStat /= 0) RETURN


   !-------------------------------------------------------------------------------------------------
   ! Read the file with the scales (les or dns)
   !-------------------------------------------------------------------------------------------------

   CALL ReadCom( UnWind, TRIM( FileName ), 'First line', TmpErrStat )   ! add ErrMsg

      ! Errors occured?
   ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
   ErrStat  = MAX(ErrStat, TmpErrStat)
   IF ( ErrStat >= AbortErrLev ) RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), CTVertShft, 'CTVertShft', ' ', TmpErrStat ) ! add ErrMsg

      ! Errors occured?
   ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
   ErrStat  = MAX(ErrStat, TmpErrStat)
   IF ( ErrStat >= AbortErrLev ) RETURN


   DO I = 1,3
      CALL ReadVar( UnWind, TRIM( FileName ), CTScale(I), 'CTScale('//TRIM(Num2LStr(I))//')', ' ', TmpErrStat ) ! add ErrMsg

         ! Errors occured?
      ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
      ErrStat  = MAX(ErrStat, TmpErrStat)
      IF ( ErrStat >= AbortErrLev ) RETURN


      CALL ReadVar( UnWind, TRIM( FileName ), CTOffset(I), 'CTOffset('//TRIM(Num2LStr(I))//')', ' ', TmpErrStat ) ! add ErrMsg

         ! Errors occured?
      ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
      ErrStat  = MAX(ErrStat, TmpErrStat)
      IF ( ErrStat >= AbortErrLev ) RETURN

   END DO !I


   CALL ReadVar( UnWind, TRIM( FileName ), NumCTy, 'NumCTy', ' ', TmpErrStat ) ! add ErrMsg

      ! Errors occured?
   ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
   ErrStat  = MAX(ErrStat, TmpErrStat)
   IF ( ErrStat >= AbortErrLev ) RETURN


   CALL ReadVar( UnWind, TRIM( FileName ), NumCTz, 'NumCTz', ' ', TmpErrStat ) ! add ErrMsg

      ! Errors occured?
   ErrMsg   = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
   ErrStat  = MAX(ErrStat, TmpErrStat)
   IF ( ErrStat >= AbortErrLev ) RETURN


   !-------------------------------------------------------------------------------------------------
   ! Close the file with the scales (les or dns)
   !-------------------------------------------------------------------------------------------------

   CLOSE( UnWind )


   RETURN

END SUBROUTINE ReadCTScales
!====================================================================================================
SUBROUTINE CT_Terminate( ErrStat, ErrMsg )
! This subroutine closes files, deallocates memory, and un-sets the initialization flag
!----------------------------------------------------------------------------------------------------

   INTEGER,       INTENT(OUT) :: ErrStat           ! return ErrID_None if no errors, ErrID level otherwise.
   CHARACTER(*),  INTENT(OUT) :: ErrMsg            ! Message about the error that occurred.


      ! Temporary error handling Variables

   INTEGER                          :: TmpErrStat                             ! Temporary Error Status
   CHARACTER(LEN(ErrMsg))           :: TmpErrMsg                              ! Temporary Error Message


      ! Initialize variables
   ErrStat  = ErrID_None
   ErrMsg   = ''
   TmpErrMsg   = ''
   TmpErrStat  = 0


   CLOSE( CTWindUnit )

   IF ( ALLOCATED( CTvelU    ) )  DEALLOCATE( CTvelU,    STAT=TmpErrStat,  ERRMSG=TmpErrMsg )
   IF (TmpErrStat /= 0) THEN
      ErrMsg   = TRIM(ErrMsg)//' InflowWind:CTWind: Error occured during de-allocation of CTvelU.'//TmpErrMsg
      ErrStat  = ErrID_Fatal
   ENDIF

   IF ( ALLOCATED( CTvelV    ) )  DEALLOCATE( CTvelV,    STAT=TmpErrStat,  ERRMSG=TmpErrMsg )
   IF (TmpErrStat /= 0) THEN
      ErrMsg   = TRIM(ErrMsg)//' InflowWind:CTWind: Error occured during de-allocation of CTvelV.'//TmpErrMsg
      ErrStat  = ErrID_Fatal
   ENDIF

   IF ( ALLOCATED( CTvelW    ) )  DEALLOCATE( CTvelW,    STAT=TmpErrStat,  ERRMSG=TmpErrMsg )
   IF (TmpErrStat /= 0) THEN
      ErrMsg   = TRIM(ErrMsg)//' InflowWind:CTWind: Error occured during de-allocation of CTvelW.'//TmpErrMsg
      ErrStat  = ErrID_Fatal
   ENDIF

   IF ( ALLOCATED( Tdata     ) )  DEALLOCATE( Tdata,     STAT=TmpErrStat,  ERRMSG=TmpErrMsg )
   IF (TmpErrStat /= 0) THEN
      ErrMsg   = TRIM(ErrMsg)//' InflowWind:CTWind: Error occured during de-allocation of Tdata.'//TmpErrMsg
      ErrStat  = ErrID_Fatal
   ENDIF

   IF ( ALLOCATED( TimeStpCT ) )  DEALLOCATE( TimeStpCT, STAT=TmpErrStat,  ERRMSG=TmpErrMsg )
   IF (TmpErrStat /= 0) THEN
      ErrMsg   = TRIM(ErrMsg)//' InflowWind:CTWind: Error occured during de-allocation of TimeStpCT.'//TmpErrMsg
      ErrStat  = ErrID_Fatal
   ENDIF

   TimeIndx = 0

END SUBROUTINE CT_Terminate
!====================================================================================================
END MODULE CTWind
