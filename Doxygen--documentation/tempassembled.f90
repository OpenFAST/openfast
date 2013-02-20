PROGRAM HHWind_Test
!-----------------------------------------------------------------------------------------------------------------------------------
!  This program is for testing the wind modules during development.
!
!
!  v1.00.00a-adp  -- HHWind tested
!
!
!
!
!
!
!
!-----------------------------------------------------------------------------------------------------------------------------------


   USE NWTC_Library
   USE IfW_HHWind



   IMPLICIT NONE

   TYPE( ProgDesc ), PARAMETER                        :: ProgInfo = ProgDesc("Wind_Test","v1.00.00a-adp","6-Feb-2013")



      ! Error handling
   CHARACTER(2048)                                    :: ErrMsg
   INTEGER(IntKi)                                     :: ErrStat


      ! Local variables
   CHARACTER(1024)                                    :: WindFileName
   INTEGER(IntKi)                                     :: UnWind         !FIXME: this should be removed when fully converted to the modular framework
   REAL(DbKi)                                         :: Time
   REAL(ReKi)                                         :: WindPosition(3)
   REAL(ReKi)                                         :: WindVelocity(3)


!FIXME: remove this one
   TYPE(HH_Info)                                      :: HHInitInfo



   !-=- Setup some things  -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

   CALL NWTC_Init
   CALL DispNVD( ProgInfo )

   Time = 2.0

      ! setup the file info
   WindFileName      = "../../Samples/Steady.wnd"           ! HHWind file
   HHInitInfo%ReferenceHeight = 80.                         ! meters
   HHInitInfo%Width           = 100.                        ! meters

   WindPosition(1)   = 0.0                                  ! longitudinal front/back of tower
   WindPosition(2)   = 0.0                                  ! lateral position left/right of tower
   WindPosition(3)   = HHInitInfo%ReferenceHeight           ! Height above ground

      ! find a unit number to use, then check the errors
   CALL GetNewUnit(UnWind,ErrStat,ErrMsg)

   IF ( ErrStat >= ErrID_Severe ) THEN
      CALL ProgAbort(ErrMsg)
   ELSEIF ( ErrStat /= ErrID_None ) THEN
      CALL ProgWarn(ErrMsg)
   ENDIF
   ErrMsg   = ""
   ErrStat  = ErrID_None



   !-=- Initialize the module  -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

   CALL WrScr(NewLine//" Initializing HHWind"//NewLine)

   CALL IfW_HHWind_Init( UnWind, WindFileName, HHInitInfo, ErrStat, ErrMsg )

   IF ( ErrStat >= ErrID_Severe ) THEN
      CALL ProgAbort(ErrMsg)
   ELSEIF ( ErrStat /= ErrID_None ) THEN
      CALL ProgWarn(ErrMsg)
   ENDIF
   ErrMsg   = ""
   ErrStat  = ErrID_None


   !-=- Simple call to get windspeed at just the hub -=-=-=-=-=-=-=-

   CALL WrScr(" Calculating wind velocity:")

   WindVelocity = IfW_HHWind_GetWindSpeed( Time, WindPosition, ErrStat, ErrMsg )

   IF ( ErrStat >= ErrID_Severe ) THEN
      CALL ProgAbort(ErrMsg)
   ELSEIF ( ErrStat /= ErrID_None ) THEN
      CALL ProgWarn(ErrMsg)
   ENDIF
   ErrMsg   = ""
   ErrStat  = ErrID_None


   !-=- Write out some info about what we just did -=-=-=-=-=-=-=-=-

   CALL WrScr("   Time: "//TRIM(Num2LStr(Time)))
   CALL WrScr("          (x, y, z)          (U, V, W)")

   CALL WrScr("          ("//TRIM(Num2LStr(WindPosition(1)))//", "//TRIM(Num2LStr(WindPosition(2)))//", " &
                           //TRIM(Num2LStr(WindPosition(3)))//")        (" &
                           //TRIM(Num2LStr(WindVelocity(1)))//", "//TRIM(Num2LStr(WindVelocity(2)))//", " &
                           //TRIM(Num2LStr(WindVelocity(3)))//")")

   !-=- Close everything -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

   CALL IfW_HHWind_Terminate( ErrStat, ErrMsg )

   IF ( ErrStat >= ErrID_Severe ) THEN
      CALL ProgAbort(ErrMsg)
   ELSEIF ( ErrStat /= ErrID_None ) THEN
      CALL ProgWarn(ErrMsg)
   ENDIF
   ErrMsg   = ""
   ErrStat  = ErrID_None


END PROGRAM


! TODO:
!  -- Fix the error handling so it is correct
!  -- Rename routines to match framework
!  -- Convert the Init routine to the framework
!  -- Convert the Calc routine to the framework -- unfunction it
!  -- Convert the End  routine to the framework
!  -- Figure out how to store the info we need (types)
MODULE IfW_HHWind
! This module contains all the data and procedures that define hub-height wind files. This could
! more accurately be called a point wind file since the wind speed at any point is calculated by
! shear applied to the point where wind is defined.  It is basically uniform wind over the rotor disk.
! The entire file is read on initialization, then the columns that make up the wind file are
! interpolated to the time requested, and wind is calculated based on the location in space.
!
! the file contains header information (rows that contain "!"), followed by numeric data stored in
! 8 columns:   (1) Time                                  [s]
!              (2) Horizontal wind speed       (V)       [m/s]
!              (3) Wind direction              (Delta)   [deg]
!              (4) Vertical wind speed         (VZ)      [m/s]
!              (5) Horizontal linear shear     (HLinShr) [-]
!              (6) Vertical power-law shear    (VShr)    [-]
!              (7) Vertical linear shear       (VLinShr) [-]
!              (8) Gust (horizontal) velocity  (VGust)   [m/s]
!
! The horizontal wind speed at (X, Y, Z) is then calculated using the interpolated columns by
!   Vh = V * ( Z/RefHt ) ** VShr                                        ! power-law wind shear
!      + V * HLinShr/RefWid * ( Y * COS(Delta) + X * SIN(Delta) )       ! horizontal linear shear
!      + V * VLinShr/RefWid * ( Z-RefHt )                               ! vertical linear shear
!      + VGust                                                          ! gust speed
!----------------------------------------------------------------------------------------------------

   USE                           NWTC_Library
   USE                           SharedInflowDefs
   USE                           InflowWind_Module_Types

   IMPLICIT                      NONE
   PRIVATE


   REAL(ReKi), ALLOCATABLE      :: Tdata  (:)                              ! Time array from the HH wind file
   REAL(ReKi), ALLOCATABLE      :: DELTA  (:)                              ! HH Wind direction (angle)
   REAL(ReKi), ALLOCATABLE      :: V      (:)                              ! HH horizontal wind speed
   REAL(ReKi), ALLOCATABLE      :: VZ     (:)                              ! wind, including tower shadow, along the Z axis
   REAL(ReKi), ALLOCATABLE      :: HSHR   (:)                              ! HH Horizontal linear shear
   REAL(ReKi), ALLOCATABLE      :: VSHR   (:)                              ! HH vertical shear exponent
   REAL(ReKi), ALLOCATABLE      :: VLINSHR(:)                              ! HH vertical linear shear
   REAL(ReKi), ALLOCATABLE      :: VGUST  (:)                              ! HH wind gust

   REAL(ReKi)                   :: LinearizeDels(7)                        ! The delta values for linearization -- perhaps at some point, this could be T/F and we determine the deltas by sqrt(eps) or something similar
   REAL(ReKi)                   :: RefHt                                   ! reference height; was HH (hub height); used to center the wind
   REAL(ReKi)                   :: RefWid                                  ! reference width; was 2*R (=rotor diameter); used to scale the linear shear

   INTEGER                      :: NumDataLines
!FIXME: move to otherstates
   INTEGER, SAVE                :: TimeIndx = 0                            ! An index into the Tdata array (to allow us faster searching, starting search from previous one)

!FIXME: move to parameters
   LOGICAL, SAVE                :: Linearize = .FALSE.                     ! If this is TRUE, we are linearizing

!FIXME: move this elsewhere
   TYPE, PUBLIC                 :: HH_Info
      REAL(ReKi)                :: ReferenceHeight
      REAL(ReKi)                :: Width
   END TYPE HH_Info

   PUBLIC                       :: IfW_HHWind_Init
   PUBLIC                       :: IfW_HHWind_GetWindSpeed
   PUBLIC                       :: IfW_HHWind_Terminate
   PUBLIC                       :: IfW_HHWind_SetLinearizeDels
!   PUBLIC                       :: HH_Get_ADhack_WindSpeed                  ! REMOVE THIS!!!!

CONTAINS
!====================================================================================================
SUBROUTINE IfW_HHWind_Init(UnWind, WindFile, WindInfo, ErrStat, ErrMsg)
! A subroutine to initialize the HHWind module.  It reads the HH file and stores the data in an
! array to use later.  It requires an initial reference height (hub height) and width (rotor diameter),
! both in meters, which are used to define the volume where wind velocities will be calculated.  This
! information is necessary because of the way the shears are defined.
!----------------------------------------------------------------------------------------------------

      ! Passed Variables:

   INTEGER,       INTENT(IN)        :: UnWind                        ! unit number for reading wind files
   TYPE(HH_Info), INTENT(IN)        :: WindInfo                      ! Additional information needed to initialize this wind type

   CHARACTER(*),  INTENT(IN)        :: WindFile                      ! Name of the text HH wind file

      ! Error handling
   INTEGER,       INTENT(OUT)       :: ErrStat                       ! determines if an error has been encountered
   CHARACTER(*),  INTENT(OUT)       :: ErrMsg                        ! A message about the error

      ! local variables

   INTEGER,       PARAMETER         :: NumCols = 8                   ! Number of columns in the HH file
   REAL(ReKi)                       :: TmpData(NumCols)              ! Temp variable for reading all columns from a line
   REAL(ReKi)                       :: DelDiff                       ! Temp variable for storing the direction difference

   INTEGER                          :: I
   INTEGER                          :: NumComments
   INTEGER                          :: ILine                         ! Counts the line number in the file
   INTEGER, PARAMETER               :: MaxTries = 100
   CHARACTER(1024)                  :: Line                          ! Temp variable for reading whole line from file

      ! Temporary variables for error handling
   INTEGER                          :: TmpErrStat                    ! Temp variable for the error status
   CHARACTER(1024)                  :: TmpErrMsg                     ! Temp variable for the error message


   !-------------------------------------------------------------------------------------------------
   ! Check that it's not already initialized
   !-------------------------------------------------------------------------------------------------

   IF ( TimeIndx /= 0 ) THEN
      CALL WrScr( ' HHWind has already been initialized.' )
      ErrStat = 1
      RETURN
   ELSE
      ErrStat = 0
!      CALL NWTC_Init()    ! Initialized in IfW_Init

      LinearizeDels(:) = 0.0
      Linearize        = .FALSE.
   END IF


   !-------------------------------------------------------------------------------------------------
   ! Open the file for reading
   !-------------------------------------------------------------------------------------------------
   CALL OpenFInpFile (UnWind, TRIM(WindFile), ErrStat)

   IF ( ErrStat /= 0 ) RETURN

   !-------------------------------------------------------------------------------------------------
   ! Find the number of comment lines
   !-------------------------------------------------------------------------------------------------
   LINE = '!'                          ! Initialize the line for the DO WHILE LOOP
   NumComments = -1

   DO WHILE (INDEX( LINE, '!' ) > 0 ) ! Lines containing "!" are treated as comment lines
      NumComments = NumComments + 1

      READ(UnWind,'( A )',IOSTAT=ErrStat) LINE

      IF ( ErrStat /=0 ) THEN
         CALL WrScr ( ' Error reading from HH wind file on line '//TRIM(Num2LStr(NumComments))//'.' )
         RETURN
      END IF

   END DO !WHILE

   !-------------------------------------------------------------------------------------------------
   ! Find the number of data lines
   !-------------------------------------------------------------------------------------------------
   NumDataLines = 0

   READ(LINE,*,IOSTAT=ErrStat) ( TmpData(I), I=1,NumCols )

   DO WHILE (ErrStat == 0)  ! read the rest of the file (until an error occurs)
      NumDataLines = NumDataLines + 1

      READ(UnWind,*,IOSTAT=ErrStat) ( TmpData(I), I=1,NumCols )

   END DO !WHILE


   IF (NumDataLines < 1) THEN
      CALL WrScr ( ' Error reading data from HH wind file on line '//TRIM(Num2LStr(NumDataLines+NumComments))//'.' )
      RETURN
   ELSE
      CALL WrScr ( ' Reading '//TRIM(Num2LStr(NumDataLines))//' lines of data from the HH wind file.' )
   END IF


   !-------------------------------------------------------------------------------------------------
   ! Allocate arrays for the HH data
   !-------------------------------------------------------------------------------------------------
   ! BJJ note: If the subroutine AllocAry() is called, the CVF compiler with A2AD does not work
   !   properly.  The arrays are not properly read even though they've been allocated.
   !-------------------------------------------------------------------------------------------------

   IF (.NOT. ALLOCATED(Tdata) ) THEN
      ALLOCATE ( Tdata(NumDataLines) , STAT=ErrStat )
      IF ( ErrStat /=0 ) THEN
         CALL WrScr( 'Error allocating memory for the HH time array.' )
         RETURN
      END IF
   END IF

   IF (.NOT. ALLOCATED(V) ) THEN
      ALLOCATE ( V(NumDataLines) , STAT=ErrStat )
      IF ( ErrStat /=0 ) THEN
         CALL WrScr( 'Error allocating memory for the HH horizontal wind speed array.' )
         RETURN
      END IF
   END IF

   IF (.NOT. ALLOCATED(Delta) ) THEN
      ALLOCATE ( Delta(NumDataLines) , STAT=ErrStat )
      IF ( ErrStat /=0 ) THEN
         CALL WrScr( 'Error allocating memory for the HH wind direction array.' )
         RETURN
      END IF
   END IF

   IF (.NOT. ALLOCATED(VZ) ) THEN
      ALLOCATE ( VZ(NumDataLines) , STAT=ErrStat )
      IF ( ErrStat /=0 ) THEN
         CALL WrScr( 'Error allocating memory for the HH vertical wind speed array.' )
         RETURN
      END IF
   END IF

   IF (.NOT. ALLOCATED(HShr) ) THEN
      ALLOCATE ( HShr(NumDataLines) , STAT=ErrStat )
      IF ( ErrStat /=0 ) THEN
         CALL WrScr( 'Error allocating memory for the HH horizontal linear shear array.' )
         RETURN
      END IF
   END IF

   IF (.NOT. ALLOCATED(VShr) ) THEN
      ALLOCATE ( VShr(NumDataLines) , STAT=ErrStat )
      IF ( ErrStat /=0 ) THEN
         CALL WrScr( 'Error allocating memory for the HH vertical power-law shear exponent array.' )
         RETURN
      END IF
   END IF

   IF (.NOT. ALLOCATED(VLinShr) ) THEN
      ALLOCATE ( VLinShr(NumDataLines) , STAT=ErrStat )
      IF ( ErrStat /=0 ) THEN
         CALL WrScr( 'Error allocating memory for the HH vertical linear shear array.' )
         RETURN
      END IF
   END IF

   IF (.NOT. ALLOCATED(VGust) ) THEN
      ALLOCATE ( VGust(NumDataLines) , STAT=ErrStat )
      IF ( ErrStat /=0 ) THEN
         CALL WrScr( 'Error allocating memory for the HH gust velocity array.' )
         RETURN
      END IF
   END IF


   !-------------------------------------------------------------------------------------------------
   ! Rewind the file (to the beginning) and skip the comment lines
   !-------------------------------------------------------------------------------------------------
   REWIND( UnWind )

   DO I=1,NumComments
      CALL ReadCom( UnWind, TRIM(WindFile), 'Header line #'//TRIM(Num2LStr(I)), ErrStat )
      IF ( ErrStat /= 0 ) RETURN
   END DO !I


   !-------------------------------------------------------------------------------------------------
   ! Read the data arrays
   !-------------------------------------------------------------------------------------------------

   DO I=1,NumDataLines

      CALL ReadAry( UnWind, TRIM(WindFile), TmpData(1:NumCols), NumCols, 'TmpData', &
                'Data from HH line '//TRIM(Num2LStr(NumComments+I)), ErrStat )
      IF (ErrStat /= 0) RETURN

      Tdata(  I) = TmpData(1)
      V(      I) = TmpData(2)
      Delta(  I) = TmpData(3)*D2R
      VZ(     I) = TmpData(4)
      HShr(   I) = TmpData(5)
      VShr(   I) = TmpData(6)
      VLinSHR(I) = TmpData(7)
      VGust(  I) = TmpData(8)

   END DO !I


   !-------------------------------------------------------------------------------------------------
   ! Make sure the wind direction isn't jumping more than 180 degrees between any 2 consecutive
   ! input times.  (Avoids interpolation errors with modular arithemetic.)
   !-------------------------------------------------------------------------------------------------

   DO I=2,NumDataLines

      ILine = 1

      DO WHILE ( ILine < MaxTries )

         DelDiff = ( Delta(I) - Delta(I-1) )

         IF ( ABS( DelDiff ) < Pi ) EXIT  ! exit inner loop

         Delta(I) = Delta(I) - SIGN( TwoPi, DelDiff )

         ILine = ILine + 1

      END DO

      IF ( ILine >= MaxTries ) THEN
         CALL WrScr( ' Error calculating wind direction from HH file. Delta(' &
               // TRIM(Num2LStr(I  )) // ') = ' // TRIM(Num2LStr(Delta(I))) // '; Delta(' &
               // TRIM(Num2LStr(I+1)) // ') = ' // TRIM(Num2LStr(Delta(I+1))) )
         ErrStat = 1
      END IF


   END DO !I


   !-------------------------------------------------------------------------------------------------
   ! Close the file
   !-------------------------------------------------------------------------------------------------

   CLOSE( UnWind )


   !-------------------------------------------------------------------------------------------------
   ! Print warnings and messages
   !-------------------------------------------------------------------------------------------------
!   CALL WrScr ( ' Processed '//TRIM( Num2LStr( NumDataLines ) )//' records of HH data' )


   IF ( Tdata(1) > 0.0 ) THEN
      CALL ProgWarn( 'The hub-height wind file : "'//TRIM(ADJUSTL(WindFile))//'" starts at a time '// &
                     'greater than zero. Interpolation errors may result.')
   ENDIF

   IF ( NumDataLines == 1 ) THEN
      CALL WrScr( ' Only 1 line in HH wind file. Steady, hub-height horizontal wind speed = '//TRIM(Num2LStr(V(1)))//' m/s.' )
   END IF


   !-------------------------------------------------------------------------------------------------
   ! Set the initial index into the time array (it indicates that we've initialized the module, too)
   ! and initialize the spatial scaling for the wind calculations
   !-------------------------------------------------------------------------------------------------
   TimeIndx = 1

   RefHt  = WindInfo%ReferenceHeight
   RefWid = WindInfo%Width


   RETURN

END SUBROUTINE IfW_HHWind_Init
!====================================================================================================
FUNCTION IfW_HHWind_GetWindSpeed(Time, InputPosition, ErrStat, ErrMsg)
! This subroutine linearly interpolates the columns in the HH input file to get the values for
! the requested time, then uses the interpolated values to calclate the wind speed at a point
! in space represented by InputPosition.
!----------------------------------------------------------------------------------------------------

   REAL(DbKi),          INTENT(IN)  :: Time                 ! time from the start of the simulation
   REAL(ReKi),          INTENT(IN)  :: InputPosition(3)     ! input information: positions X,Y,Z
   INTEGER,             INTENT(OUT) :: ErrStat              ! error status
   CHARACTER(*),        INTENT(OUT) :: ErrMsg               ! The error message
!FIXME: make this go away when convert to subroutine
   REAL(ReKi)                       :: IfW_HHwind_GetWindSpeed(3)   ! return velocities (U,V,W)

   REAL(ReKi)                       :: CosDelta             ! cosine of Delta_tmp
   REAL(ReKi)                       :: Delta_tmp            ! interpolated Delta   at input TIME
   REAL(ReKi)                       :: HShr_tmp             ! interpolated HShr    at input TIME
   REAL(ReKi)                       :: P                    ! temporary storage for slope (in time) used in linear interpolation
   REAL(ReKi)                       :: SinDelta             ! sine of Delta_tmp
   REAL(ReKi)                       :: V_tmp                ! interpolated V       at input TIME
   REAL(ReKi)                       :: VGust_tmp            ! interpolated VGust   at input TIME
   REAL(ReKi)                       :: VLinShr_tmp          ! interpolated VLinShr at input TIME
   REAL(ReKi)                       :: VShr_tmp             ! interpolated VShr    at input TIME
   REAL(ReKi)                       :: VZ_tmp               ! interpolated VZ      at input TIME
   REAL(ReKi)                       :: V1                   ! temporary storage for horizontal velocity



   !-------------------------------------------------------------------------------------------------
   ! verify the module was initialized first
   !-------------------------------------------------------------------------------------------------

   IF ( TimeIndx == 0 ) THEN
      ErrMsg   = ' Error: Call HH_Init() before getting wind speed.'
      ErrStat  = ErrID_Fatal         ! Fatal since no data returned
      RETURN
   ELSE
      ErrStat = 0
   END IF


   !-------------------------------------------------------------------------------------------------
   ! Linearly interpolate in time (or used nearest-neighbor to extrapolate)
   ! (compare with NWTC_Num.f90\InterpStpReal)
   !-------------------------------------------------------------------------------------------------

    IF ( Linearize ) THEN  !get the perturbed wind speed

      TimeIndx      = 1
      V_tmp         = V      (1) + LinearizeDels(1)
      Delta_tmp     = Delta  (1) + LinearizeDels(2)
      VZ_tmp        = VZ     (1) + LinearizeDels(3)
      HShr_tmp      = HShr   (1) + LinearizeDels(4)
      VShr_tmp      = VShr   (1) + LinearizeDels(5)
      VLinShr_tmp   = VLinShr(1) + LinearizeDels(6)
      VGust_tmp     = VGust  (1) + LinearizeDels(7)

      ! Let's check the limits.
   ELSE IF ( Time <= Tdata(1) .OR. NumDataLines == 1 )  THEN

      TimeIndx      = 1
      V_tmp         = V      (1)
      Delta_tmp     = Delta  (1)
      VZ_tmp        = VZ     (1)
      HShr_tmp      = HShr   (1)
      VShr_tmp      = VShr   (1)
      VLinShr_tmp   = VLinShr(1)
      VGust_tmp     = VGust  (1)

   ELSE IF ( Time >= Tdata(NumDataLines) )  THEN

      TimeIndx      = NumDataLines - 1
      V_tmp         = V      (NumDataLines)
      Delta_tmp     = Delta  (NumDataLines)
      VZ_tmp        = VZ     (NumDataLines)
      HShr_tmp      = HShr   (NumDataLines)
      VShr_tmp      = VShr   (NumDataLines)
      VLinShr_tmp   = VLinShr(NumDataLines)
      VGust_tmp     = VGust  (NumDataLines)

   ELSE

         ! Let's interpolate!

      TimeIndx = MAX( MIN( TimeIndx, NumDataLines-1 ), 1 )

      DO

         IF ( Time < Tdata(TimeIndx) )  THEN

            TimeIndx = TimeIndx - 1

         ELSE IF ( Time >= Tdata(TimeIndx+1) )  THEN

            TimeIndx = TimeIndx + 1

         ELSE
            P           = ( Time - Tdata(TimeIndx) )/( Tdata(TimeIndx+1) - Tdata(TimeIndx) )
            V_tmp       = ( V(      TimeIndx+1) - V(      TimeIndx) )*P + V(      TimeIndx)
            Delta_tmp   = ( Delta(  TimeIndx+1) - Delta(  TimeIndx) )*P + Delta(  TimeIndx)
            VZ_tmp      = ( VZ(     TimeIndx+1) - VZ(     TimeIndx) )*P + VZ(     TimeIndx)
            HShr_tmp    = ( HShr(   TimeIndx+1) - HShr(   TimeIndx) )*P + HShr(   TimeIndx)
            VShr_tmp    = ( VShr(   TimeIndx+1) - VShr(   TimeIndx) )*P + VShr(   TimeIndx)
            VLinShr_tmp = ( VLinShr(TimeIndx+1) - VLinShr(TimeIndx) )*P + VLinShr(TimeIndx)
            VGust_tmp   = ( VGust(  TimeIndx+1) - VGust(  TimeIndx) )*P + VGust(  TimeIndx)
            EXIT

         END IF

      END DO

   END IF


   !-------------------------------------------------------------------------------------------------
   ! calculate the wind speed at this time
   !-------------------------------------------------------------------------------------------------

   CosDelta = COS( Delta_tmp )
   SinDelta = SIN( Delta_tmp )

   V1 = V_tmp * ( ( InputPosition(3)/RefHt ) ** VShr_tmp &                                  ! power-law wind shear
        + ( HShr_tmp   * ( InputPosition(2) * CosDelta + InputPosition(1) * SinDelta ) &    ! horizontal linear shear
        +  VLinShr_tmp * ( InputPosition(3)-RefHt ) )/RefWid  ) &                           ! vertical linear shear
        + VGUST_tmp                                                                         ! gust speed
   IfW_HHWind_GetWindSpeed(1) =  V1 * CosDelta
   IfW_HHWind_GetWindSpeed(2) = -V1 * SinDelta
   IfW_HHWind_GetWindSpeed(3) =  VZ_tmp


   RETURN

END FUNCTION IfW_HHWind_GetWindSpeed
!!====================================================================================================
!FUNCTION HH_Get_ADHack_WindSpeed(Time, InputPosition, ErrStat, ErrMsg)
!! This subroutine linearly interpolates the columns in the HH input file to get the values for
!! the requested time, then uses the interpolated values to calclate the wind speed at a point
!! in space represented by InputPosition. THIS FUNCTION SHOULD BE REMOVED!!!!! (used for DISK VEL ONLY)
!!----------------------------------------------------------------------------------------------------
!
!   REAL(ReKi),          INTENT(IN)  :: Time                 ! time from the start of the simulation
!   REAL(ReKi),          INTENT(IN)  :: InputPosition(3)     ! input information: positions X,Y,Z   -   NOT USED HERE!!!
!   INTEGER,             INTENT(OUT) :: ErrStat              ! error status
!   CHARACTER(*),        INTENT(OUT) :: ErrMsg               ! The error message
!   REAL(ReKi)                       :: HH_Get_ADHack_WindSpeed(3)      ! return velocities (U,V,W)
!
!   REAL(ReKi)                       :: Delta_tmp            ! interpolated Delta   at input TIME
!   REAL(ReKi)                       :: P                    ! temporary storage for slope (in time) used in linear interpolation
!   REAL(ReKi)                       :: V_tmp                ! interpolated V       at input TIME
!   REAL(ReKi)                       :: VZ_tmp               ! interpolated VZ      at input TIME
!
!
!   !-------------------------------------------------------------------------------------------------
!   ! verify the module was initialized first
!   !-------------------------------------------------------------------------------------------------
!
!   IF ( TimeIndx == 0 ) THEN
!      ErrMsg   = ' Error: Call HH_Init() before getting wind speed.'
!      ErrStat  = ErrID_Fatal      ! Fatal since no data returned
!      RETURN
!   ELSE
!      ErrStat = 0
!   END IF
!
!
!   !-------------------------------------------------------------------------------------------------
!   ! Linearly interpolate in time (or use nearest-neighbor to extrapolate)
!   ! (compare with NWTC_Num.f90\InterpStpReal)
!   !-------------------------------------------------------------------------------------------------
!
!     ! Let's check the limits.
!
!   IF ( Time <= Tdata(1) .OR. NumDataLines == 1)  THEN
!
!      TimeIndx      = 1
!      V_tmp         = V      (1)
!      Delta_tmp     = Delta  (1)
!      VZ_tmp        = VZ     (1)
!
!   ELSE IF ( Time >= Tdata(NumDataLines) )  THEN
!
!      TimeIndx      = NumDataLines - 1
!      V_tmp         = V      (NumDataLines)
!      Delta_tmp     = Delta  (NumDataLines)
!      VZ_tmp        = VZ     (NumDataLines)
!
!   ELSE
!
!         ! Let's interpolate!
!
!      TimeIndx = MAX( MIN( TimeIndx, NumDataLines-1 ), 1 )
!
!      DO
!
!         IF ( Time < Tdata(TimeIndx) )  THEN
!
!            TimeIndx = TimeIndx - 1
!
!         ELSE IF ( Time >= Tdata(TimeIndx+1) )  THEN
!
!            TimeIndx = TimeIndx + 1
!
!         ELSE
!            P           = ( Time - Tdata(TimeIndx) )/( Tdata(TimeIndx+1) - Tdata(TimeIndx) )
!            V_tmp       = ( V(      TimeIndx+1) - V(      TimeIndx) )*P + V(      TimeIndx)
!            Delta_tmp   = ( Delta(  TimeIndx+1) - Delta(  TimeIndx) )*P + Delta(  TimeIndx)
!            VZ_tmp      = ( VZ(     TimeIndx+1) - VZ(     TimeIndx) )*P + VZ(     TimeIndx)
!            EXIT
!
!         END IF
!
!      END DO
!
!   END IF
!
!   !-------------------------------------------------------------------------------------------------
!   ! calculate the wind speed at this time
!   !-------------------------------------------------------------------------------------------------
!   HH_Get_ADHack_WindSpeed(1) =  V_tmp * COS( Delta_tmp )
!   HH_Get_ADHack_WindSpeed(2) = -V_tmp * SIN( Delta_tmp )
!   HH_Get_ADHack_WindSpeed(3) =  VZ_tmp
!
!
!   RETURN
!
!END FUNCTION HH_Get_ADHack_WindSpeed
!====================================================================================================

!FIXME: might need to move this into states???
SUBROUTINE IfW_HHWind_SetLinearizeDels( Perturbations, ErrStat, ErrMsg )
! This subroutine sets the perturbation values for the linearization scheme.
!----------------------------------------------------------------------------------------------------

   REAL(ReKi),          INTENT(IN)  :: Perturbations(7)     ! purturbations for each of the 7 input parameters
   INTEGER,             INTENT(OUT) :: ErrStat              ! time from the start of the simulation
   CHARACTER(*),        INTENT(OUT) :: ErrMsg               ! Error Message

   !-------------------------------------------------------------------------------------------------
   ! verify the module was initialized first
   !-------------------------------------------------------------------------------------------------

   IF ( TimeIndx == 0 ) THEN
      ErrMsg   = ' Error: Call HH_Init() before getting wind speed.'
      ErrStat  = ErrID_Fatal        ! Fatal since no data returned
      RETURN
   ELSE
      ErrStat = 0
   END IF

   Linearize = .TRUE.
   LinearizeDels(:) = Perturbations(:)

   RETURN

END SUBROUTINE IfW_HHWind_SetLinearizeDels
!====================================================================================================
SUBROUTINE IfW_HHWind_Terminate(ErrStat,ErrMsg)

      ! Error Handling

   INTEGER,          INTENT(OUT)    :: ErrStat                       ! determines if an error has been encountered
   CHARACTER(1024),  INTENT(OUT)    :: ErrMsg                        ! Message about errors

      ! Local Variables

   INTEGER                          :: SumErrs  !FIXME: this is depricated!!!!

      !-=- Initializ the routine -=-

   SumErrs = 0

   IF ( ALLOCATED(Tdata  ) ) DEALLOCATE( Tdata,   STAT=ErrStat )
   SumErrs = SumErrs + ABS(ErrStat)

   IF ( ALLOCATED(DELTA  ) ) DEALLOCATE( DELTA,   STAT=ErrStat )
   SumErrs = SumErrs + ABS(ErrStat)

   IF ( ALLOCATED(V      ) ) DEALLOCATE( V,       STAT=ErrStat )
   SumErrs = SumErrs + ABS(ErrStat)

   IF ( ALLOCATED(VZ     ) ) DEALLOCATE( VZ,      STAT=ErrStat )
   SumErrs = SumErrs + ABS(ErrStat)

   IF ( ALLOCATED(HSHR   ) ) DEALLOCATE( HSHR,    STAT=ErrStat )
   SumErrs = SumErrs + ABS(ErrStat)

   IF ( ALLOCATED(VSHR   ) ) DEALLOCATE( VSHR,    STAT=ErrStat )
   SumErrs = SumErrs + ABS(ErrStat)

   IF ( ALLOCATED(VGUST  ) ) DEALLOCATE( VGUST,   STAT=ErrStat )
   SumErrs = SumErrs + ABS(ErrStat)

   IF ( ALLOCATED(VLINSHR) ) DEALLOCATE( VLINSHR, STAT=ErrStat )
   SumErrs = SumErrs + ABS(ErrStat)

   ErrStat  = SumErrs
   TimeIndx = 0

END SUBROUTINE IfW_HHWind_Terminate
!====================================================================================================
END MODULE IfW_HHWind
!MODULARIZATION NOTES:
!  This file will eventually be replaced with an autogenerated one by the registry program
!  So, for now we will modify this as we develop, then setup the txt file for registry program to use
!     1) modify this to match new style
!     2) setup the InflowWind.txt file for the registry program
!     3) run the registry program to generate the InflowWind_Types.f90 file
!     4) verify that works correctly before removing this (compile with each and verify output)
!
!
!----------------------------------------------------------------------------------------------------
!FIXME: rename
MODULE SharedInflowDefs
! This module is used to define shared types and parameters that are used in the module InflowWind.
! 7 Oct 2009    B. Jonkman, NREL/NWTC
!----------------------------------------------------------------------------------------------------
!
!..................................................................................................................................
! LICENSING
! Copyright (C) 2012  National Renewable Energy Laboratory
!
!    This file is part of InflowWind.
!
!    InflowWind is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with InflowWind.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************

   USE NWTC_Library                                               ! Precision module
   IMPLICIT NONE
!   PRIVATE


      !---- Initialization data ---------------------------------------------------------------------
   TYPE, PUBLIC :: IfW_InitInputType
      ! Define inputs that the initialization routine may need here:
      ! e.g., the name of the input file, the file root name, etc.
         ! Filename and type
      CHARACTER(1024)               :: WindFileName
      INTEGER                       :: WindFileType
         ! Configuration Info
      REAL(ReKi)                    :: ReferenceHeight        ! reference height for HH and/or 4D winds (was hub height), in meters
      REAL(ReKi)                    :: Width                  ! width of the HH file (was 2*R), in meters
   END TYPE IfW_InitInputType


      ! ..... States ..............................................................................................................

   TYPE, PUBLIC :: IfW_ContinuousStateType
         ! Define continuous (differentiable) states here:
      REAL(ReKi) :: DummyContState                                            ! If you have continuous states, remove this variable
         ! If you have loose coupling with a variable-step integrator, store the actual time associated with the continuous states
         !    here (eg., REAL(DbKi) :: ContTime):
   END TYPE IfW_ContinuousStateType


   TYPE, PUBLIC :: IfW_DiscreteStateType
         ! Define discrete (nondifferentiable) states here:
      REAL(ReKi) :: DummyDiscState                                            ! If you have discrete states, remove this variable
   END TYPE IfW_DiscreteStateType


   TYPE, PUBLIC :: IfW_ConstraintStateType
         ! Define constraint states here:
      REAL(ReKi) :: DummyConstrState                                          !    If you have constraint states, remove this variable
   END TYPE IfW_ConstraintStateType


   TYPE, PUBLIC :: IfW_OtherStateType
         ! Define any data that are not considered actual "states" here:
         ! e.g. data used only for optimization purposes (indices for searching in an array, copies of previous calculations of output at a given time, etc.)
      INTEGER(IntKi) :: DummyOtherState                                       !    If you have other/optimzation states, remove this variable
   END TYPE IfW_OtherStateType


      ! ..... Parameters ..........................................................................................................

   TYPE, PUBLIC :: IfW_ParameterType

         ! Define parameters here that might need to be accessed from the outside world:

         ! Filename info
      CHARACTER(1024)               :: WindFileName
      CHARACTER(1024)               :: WindFileNameRoot
      CHARACTER(3)                  :: WindFileNameExt
      INTEGER                       :: WindFileType       = 0 ! This should be initially set to match the Undef_Wind parameter in the module

         ! Location
      REAL(ReKi)                    :: ReferenceHeight        ! reference height for HH and/or 4D winds (was hub height), in meters
      REAL(ReKi)                    :: Width                  ! width of the HH file, in meters
!NOTE: might be only for HH file
      REAL(ReKi)                    :: HalfWidth              ! half the width of the HH file (was 2*R), in meters

         ! Flags
      LOGICAL                       :: CT_Flag        = .FALSE.   ! determines if coherent turbulence is used
      LOGICAL                       :: Initialized    = .FALSE.   ! did we run the initialization?


   END TYPE IfW_ParameterType


      ! ..... Inputs ..............................................................................................................

   TYPE, PUBLIC :: IfW_InputType
         ! Define inputs that are contained on the mesh here:
!     TYPE(MeshType)                            ::
         ! Define inputs that are not on this mesh here:
      Real(ReKi),ALLOCATABLE        :: Position(:,:)
   END TYPE IfW_InputType


      ! ..... Outputs .............................................................................................................

   TYPE, PUBLIC :: IfW_OutputType
         ! Define outputs that are contained on the mesh here:
!     TYPE(MeshType)                            ::
         ! Define outputs that are not on this mesh here:
      REAL(ReKi),ALLOCATABLE        :: Velocity(:,:)
   END TYPE IfW_OutputType

!-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-
!     Original stuff below that hasn't yet been moved.
!-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-



!This was commented out before.
!   TYPE, PUBLIC :: InflLoc
!      REAL(ReKi)                    :: Position(3)                ! X, Y, Z
!   END TYPE InflLoc


!This was moved to the submodules types
!   TYPE, PUBLIC :: InflIntrpOut
!      REAL(ReKi)                    :: Velocity(3)                ! U, V, W
!   END TYPE InflIntrpOut

   !-------------------------------------------------------------------------------------------------
   ! Shared parameters, defining the wind types
   ! THEY MUST BE UNIQUE!
   !-------------------------------------------------------------------------------------------------




END MODULE SharedInflowDefs



!!FIXME: add these in at some point.
!      ! ..... Jacobians ...........................................................................................................
!
!   TYPE, PUBLIC :: IfW_PartialOutputPInputType
!
!         ! Define the Jacobian of the output equations (Y) with respect to the inputs (u), dY/du (or Partial Y / Partial u):
!
!      TYPE(IfW_InputType) :: DummyOutput                                    ! If you have output equations and input data, update this variable
!
!   END TYPE IfW_PartialOutputPInputType
!
!
!   TYPE, PUBLIC :: IfW_PartialContStatePInputType
!
!         ! Define the Jacobian of the continuous state equations (X) with respect to the inputs (u), dX/du (or Partial X / Partial u):
!
!      TYPE(IfW_InputType) :: DummyContState                                 ! If you have continuous state equations and input data, update this variable
!
!   END TYPE IfW_PartialContStatePInputType
!
!
!   TYPE, PUBLIC :: IfW_PartialDiscStatePInputType
!
!         ! Define the Jacobian of the discrete state equations (Xd) with respect to the inputs (u), dXd/du (or Partial Xd / Partial u):
!
!      TYPE(IfW_InputType) :: DummyDiscState                                 ! If you have discrete state equations and input data, update this variable
!
!   END TYPE IfW_PartialDiscStatePInputType
!
!
!   TYPE, PUBLIC :: IfW_PartialConstrStatePInputType
!
!         ! Define the Jacobian of the constraint state equations (Z) with respect to the inputs (u), dZ/du (or Partial Z / Partial u):
!
!      TYPE(IfW_InputType) :: DummyConstrState                                ! If you have constraint state equations and input data, update this variable
!
!   END TYPE IfW_PartialConstrStatePInputType
!
!
!   TYPE, PUBLIC :: IfW_PartialOutputPContStateType
!
!         ! Define the Jacobian of the output equations (Y) with respect to the continuous states (x), dY/dx (or Partial Y / Partial x):
!
!      TYPE(IfW_ContinuousStateType) :: DummyOutput                                    ! If you have output equations and continuous states, update this variable
!
!   END TYPE IfW_PartialOutputPContStateType
!
!
!   TYPE, PUBLIC :: IfW_PartialContStatePContStateType
!
!         ! Define the Jacobian of the continuous state equations (X) with respect to the continuous states (x), dX/dx (or Partial X / Partial x):
!
!      TYPE(IfW_ContinuousStateType) :: DummyContState                                 ! If you have continuous state equations and continuous states, update this variable
!
!   END TYPE IfW_PartialContStatePContStateType
!
!
!   TYPE, PUBLIC :: IfW_PartialDiscStatePContStateType
!
!         ! Define the Jacobian of the discrete state equations (Xd) with respect to the continuous states (x), dXd/dx (or Partial Xd / Partial x):
!
!      TYPE(IfW_ContinuousStateType) :: DummyDiscState                                 ! If you have discrete state equations and continuous states, update this variable
!
!   END TYPE IfW_PartialDiscStatePContStateType
!
!
!   TYPE, PUBLIC :: IfW_PartialConstrStatePContStateType
!
!         ! Define the Jacobian of the constraint state equations (Z) with respect to the continuous states (x), dZ/dx (or Partial Z / Partial x):
!
!      TYPE(IfW_ContinuousStateType) :: DummyConstrState                                ! If you have constraint state equations and continuous states, update this variable
!
!   END TYPE IfW_PartialConstrStatePContStateType
!
!
!   TYPE, PUBLIC :: IfW_PartialOutputPDiscStateType
!
!         ! Define the Jacobian of the output equations (Y) with respect to the discrete states (xd), dY/dxd (or Partial Y / Partial xd):
!
!      TYPE(IfW_DiscreteStateType) :: DummyOutput                                    ! If you have output equations and discrete states, update this variable
!
!   END TYPE IfW_PartialOutputPDiscStateType
!
!
!   TYPE, PUBLIC :: IfW_PartialContStatePDiscStateType
!
!         ! Define the Jacobian of the continuous state equations (X) with respect to the discrete states (xd), dX/dxd (or Partial X / Partial xd):
!
!      TYPE(IfW_DiscreteStateType) :: DummyContState                                 ! If you have continuous state equations and discrete states, update this variable
!
!   END TYPE IfW_PartialContStatePDiscStateType
!
!
!   TYPE, PUBLIC :: IfW_PartialDiscStatePDiscStateType
!
!         ! Define the Jacobian of the discrete state equations (Xd) with respect to the discrete states (xd), dXd/dxd (or Partial Xd / Partial xd):
!
!      TYPE(IfW_DiscreteStateType) :: DummyDiscState                                 ! If you have discrete state equations and discrete states, update this variable
!
!   END TYPE IfW_PartialDiscStatePDiscStateType
!
!
!   TYPE, PUBLIC :: IfW_PartialConstrStatePDiscStateType
!
!         ! Define the Jacobian of the constraint state equations (Z) with respect to the discrete states (xd), dZ/dxd (or Partial Z / Partial xd):
!
!      TYPE(IfW_DiscreteStateType) :: DummyConstrState                                ! If you have constraint state equations and discrete states, update this variable
!
!   END TYPE IfW_PartialConstrStatePDiscStateType
!
!
!   TYPE, PUBLIC :: IfW_PartialOutputPConstrStateType
!
!         ! Define the Jacobian of the output equations (Y) with respect to the constraint states (z), dY/dz (or Partial Y / Partial z):
!
!      TYPE(IfW_ConstraintStateType) :: DummyOutput                                    ! If you have output equations and constraint states, update this variable
!
!   END TYPE IfW_PartialOutputPConstrStateType
!
!
!   TYPE, PUBLIC :: IfW_PartialContStatePConstrStateType
!
!         ! Define the Jacobian of the continuous state equations (X) with respect to the constraint states (z), dX/dz (or Partial X / Partial z):
!
!      TYPE(IfW_ConstraintStateType) :: DummyContState                                 ! If you have continuous state equations and constraint states, update this variable
!
!   END TYPE IfW_PartialContStatePConstrStateType
!
!
!   TYPE, PUBLIC :: IfW_PartialDiscStatePConstrStateType
!
!         ! Define the Jacobian of the discrete state equations (Xd) with respect to the constraint states (z), dXd/dz (or Partial Xd / Partial z):
!
!      TYPE(IfW_ConstraintStateType) :: DummyDiscState                                 ! If you have discrete state equations and constraint states, update this variable
!
!   END TYPE IfW_PartialDiscStatePConstrStateType
!
!
!   TYPE, PUBLIC :: IfW_PartialConstrStatePConstrStateType
!
!         ! Define the Jacobian of the constraint state equations (Z) with respect to the constraint states (z), dZ/dz (or Partial Z / Partial z):
!
!      TYPE(IfW_ConstraintStateType) :: DummyConstrState                                ! If you have constraint state equations and constraint states, update this variable
!
!   END TYPE IfW_PartialConstrStatePConstrStateType
!
MODULE   InflowWind_Module_Types
!FIXME: check on the name of this module. Can I call It this, or will it interfere with the one from SharedDefs?
!..................................................................................................................................
! LICENSING
! Copyright (C) 2012  National Renewable Energy Laboratory
!
!    This file is part of InflowWind.
!
!    InflowWind is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with InflowWind.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************

   USE NWTC_Library

   ! The parameters here are not considered private, but are not accessable unless the module is called.
   IMPLICIT NONE


   INTEGER,PARAMETER          :: DEFAULT_Wind = -1    ! Undetermined wind type; calls internal routine to guess what type it is
   INTEGER,PARAMETER          :: Undef_Wind   =  0    ! This is the code for an undefined WindFileType
   INTEGER,PARAMETER          :: HH_Wind      =  1    ! Hub-Height wind file
   INTEGER,PARAMETER          :: FF_Wind      =  2    ! Binary full-field wind file
   INTEGER,PARAMETER          :: UD_Wind      =  3    ! User-defined wind
   INTEGER,PARAMETER          :: FD_Wind      =  4    ! 4-dimensional wind (LES)
   INTEGER,PARAMETER          :: CTP_Wind     =  5    ! Coherent turbulence wind field (superimpose KH billow on background wind)
   INTEGER,PARAMETER          :: HAWC_Wind    =  6    ! Binary full-field wind file in HAWC format


END MODULE   InflowWind_Module_Types
