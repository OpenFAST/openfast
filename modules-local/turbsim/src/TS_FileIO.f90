!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2014  National Renewable Energy Laboratory
!
!    This file is part of TurbSim.
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
MODULE TS_FileIO

   USE                     NWTC_Library   

   use TS_Profiles
   use TSSubs
   use ts_errors
   use TS_RandNum

   
   IMPLICIT                NONE



CONTAINS

!=======================================================================
SUBROUTINE GetInput(ErrStat, ErrMsg)


  ! This subroutine is used to read parameters from the input file.

USE                  TSMods

IMPLICIT             NONE

INTEGER(IntKi)  ,             INTENT(OUT)   :: ErrStat     ! allocation status
CHARACTER(*) ,                INTENT(OUT)   :: ErrMsg      ! error message

   ! Local variables

REAL(ReKi)        :: InCVar     (2)  ! Contains the coherence parameters (used for input)
REAL(ReKi)        :: RefHt           ! Height for reference wind speed.
REAL(ReKi)        :: tmp             ! variable for estimating Ustar
REAL(ReKi)        :: TmpUary (3)     !Temporary vector to store windSpeed(z) values
REAL(ReKi)        :: TmpUstar(3)     !Temporary vector to store ustar(z) values
REAL(ReKi)        :: TmpUstarD       !Temporary ustarD value
REAL(ReKi)        :: TmpZary (3)     !Temporary vector to store height(z) values
REAL(ReKi)        :: TmpZLary(3)     !Temporary vector to store zL(z) values
REAL(ReKi)        :: URef            ! Wind speed at the reference height.

INTEGER           :: IOS             ! I/O status
INTEGER           :: TmpIndex        ! Contains the index number when searching for substrings

LOGICAL           :: getPLExp        ! Whether the power law exponent needs to be calculated
LOGICAL           :: Randomize       ! Whether to randomize the coherent turbulence scaling
LOGICAL           :: UseDefault      ! Whether or not to use a default value


CHARACTER(99)     :: Line            ! An input line
CHARACTER(1)      :: Line1           ! The first character of an input line


!UnEc = US
!Echo = .FALSE.       ! Do not echo the input into a separate file

   !==========================================================
   ! Initialize temporary variables
   !==========================================================
TmpUary    = (/ 0.0, 0.0, 0.0 /)
TmpUstarD  = 0.0
TmpUstar   = (/ 0.0, 0.0, 0.0 /)
UstarOffset= 0.0
UstarSlope = 1.0
zlOffset   = 0.0

   ! Open input file.
CALL OpenFInpFile( UI, InFile )

CALL WrScr1(' Reading the input file "'//TRIM(InFile)//'".' )

   !==========================================================
   ! Read the runtime options.
   !==========================================================

WRITE (US,"( / 'Runtime Options:' / )")
CALL ReadCom( UI, InFile, "File Heading Line 1" )
CALL ReadCom( UI, InFile, "File Heading Line 2" )
CALL ReadCom( UI, InFile, "Runtime Options Heading" )
!READ (UI,'(//)')


   ! ---------- Read Random seed 1 -----------------------

CALL ReadVar( UI, InFile, p_RandNum%RandSeed(1), "RandSeed(1)", "Random seed #1")

WRITE (US,"( I10 , 2X , 'Random seed #1' )")  p_RandNum%RandSeed(1)


   ! ---------- Read Random seed 2 -----------------------

CALL ReadVar( UI, InFile, Line, "RandSeed(2)", "Random seed #2")

   ! Check if alternate random number generator is to be used

   READ (Line,*,IOSTAT=IOS) Line1

   CALL Conv2UC( Line1 )

   IF ( (Line1 == 'T') .OR.  (Line1 == 'F') ) THEN
      CALL TS_Abort (' Invalid RNG type. ')
   ENDIF

   READ (Line,*,IOSTAT=IOS)  p_RandNum%RandSeed(2)

   IF (IOS == 0) THEN

      FormStr = "( I10 , 2X , 'Random seed #2' )"
      WRITE (US,FormStr)  p_RandNum%RandSeed(2)
      RNG_type = "NORMAL"
      p_RandNum%pRNG = pRNG_INTRINSIC

   ELSE

      RNG_type = ADJUSTL( Line )
      CALL Conv2UC( RNG_type )

      IF ( RNG_type == "RANLUX") THEN
         p_RandNum%pRNG = pRNG_RANLUX
      ELSE IF ( RNG_type == "RNSNLW") THEN
         p_RandNum%pRNG = pRNG_SNLW3
      ELSE
         CALL TS_Abort( 'Invalid alternative random number generator.' )
      ENDIF

      FormStr  = "( 4X, A6, 2X, 'Type of random number generator' )"
      WRITE (US,FormStr)  RNG_type

   ENDIF

   ! Initialize the RNG (for computing "default" values that contain random variates)

CALL RandNum_Init(p_RandNum, OtherSt_RandNum, ErrStat, ErrMsg)


   ! --------- Read the flag for writing the binary HH (GenPro) turbulence parameters. -------------

CALL ReadVar( UI, InFile, WrBHHTP, "the flag for writing the binary HH turbulence parameters", &
                                    "Output binary HH turbulence parameters?")

FormStr = "( L10 , 2X , 'Output binary HH turbulence parameters?' )"
WRITE (US,FormStr)  WrBHHTP


   ! --------- Read the flag for writing the formatted turbulence parameters. ----------------------

CALL ReadVar( UI, InFile, WrFHHTP, "the flag for writing the formatted turbulence parameters", &
                                    "Output formatted turbulence parameters?")

FormStr = "( L10 , 2X , 'Output formatted turbulence parameters?' )"
WRITE (US,FormStr)  WrFHHTP


   ! ---------- Read the flag for writing the AeroDyn HH files. -------------------------------------

CALL ReadVar( UI, InFile, WrADHH, "the flag for writing AeroDyn HH files", "Output AeroDyn HH files?")

FormStr = "( L10 , 2X , 'Output AeroDyn HH files?' )"
WRITE (US,FormStr)  WrADHH

   ! ---------- Read the flag for writing the AeroDyn FF files. ---------------------------------------

CALL ReadVar( UI, InFile, WrADFF, "the flag for writing AeroDyn FF files", "Output AeroDyn FF files?")

FormStr = "( L10 , 2X , 'Output AeroDyn FF files?' )"
WRITE (US,FormStr)  WrADFF

   ! ---------- Read the flag for writing the BLADED FF files. -----------------------------------------

CALL ReadVar( UI, InFile, WrBLFF, "the flag for writing BLADED FF files", "Output BLADED FF files?")

FormStr = "( L10 , 2X , 'Output BLADED FF files?' )"
WRITE (US,FormStr)  WrBLFF


   ! ---------- Read the flag for writing the AeroDyn tower files. --------------------------------------

CALL ReadVar( UI, InFile, WrADTWR, "the flag for writing tower data", "Output tower data?")

FormStr = "( L10 , 2X , 'Output tower data?' )"
WRITE (US,FormStr)  WrADTWR


   ! ---------- Read the flag for writing the formatted FF files. ---------------------------------------

CALL ReadVar( UI, InFile, WrFmtFF, "the flag for writing formatted FF files", "Output formatted FF files?")

FormStr = "( L10 , 2X , 'Output formatted FF files?' )"
WRITE (US,FormStr)  WrFmtFF


   ! ---------- Read the flag for writing coherent time series files. --------------------------------------

CALL ReadVar( UI, InFile, WrACT, "the flag for writing coherent time series files", "Output coherent time series files?")

FormStr = "( L10 , 2X , 'Output coherent turbulence time step file?' )"
WRITE (US,FormStr)  WrACT


   ! ---------- Read the flag for turbine rotation. -----------------------------------------------------------

CALL ReadVar( UI, InFile, Clockwise, "the flag for direction of turbine rotation", "Clockwise rotation when looking downwind?")

FormStr = "( L10 , 2X , 'Clockwise rotation when looking downwind?' )"
WRITE (US,FormStr)  Clockwise

   ! ---------- Read the flag for determining IEC scaling -----------------------------------------------------
CALL ReadVar( UI, InFile, p_IEC%ScaleIEC, "ScaleIEC, the switch for scaling IEC turbulence", &
               "Scale IEC turbulence models to specified standard deviation?")

FormStr = "( I2, ' - ', A5, 2X , 'IEC turbulence models scaled to exact specified standard deviation' )"
   SELECT CASE ( p_IEC%ScaleIEC )
      CASE (0)
         WRITE (US,FormStr)  p_IEC%ScaleIEC, "NONE"
      CASE (1, -1)   ! included the -1 for reading t/f on other systems
         WRITE (US,FormStr)  p_IEC%ScaleIEC, "HUB"
         p_IEC%ScaleIEC = 1;
      CASE (2)
         WRITE (US,FormStr)  p_IEC%ScaleIEC, "ALL"
      CASE DEFAULT
         CALL TS_Abort ( 'The value for parameter ScaleIEC must be 0, 1, or 2.' )
   ENDSELECT


   !==================================================================================
   ! Read the turbine/model specifications.
   !===================================================================================

FormStr = "( // 'Turbine/Model Specifications:' / )"
WRITE (US,FormStr)
CALL ReadCom( UI, InFile, "Turbine/Model Specifications Heading Line 1" )
CALL ReadCom( UI, InFile, "Turbine/Model Specifications Heading Line 2" )
!READ (UI,'(/)')


   ! ------------ Read in the vertical matrix dimension. ---------------------------------------------

CALL ReadVar( UI, InFile, p_grid%NumGrid_Z, "the vertical matrix dimension", "Vertical grid-point matrix dimension")

   IF ( p_grid%NumGrid_Z < 2 )  THEN
      CALL TS_Abort ( 'The matrix must be >= 2x2.' )
   ENDIF

FormStr = "( I10 , 2X , 'Vertical grid-point matrix dimension' )"
WRITE (US,FormStr)  p_grid%NumGrid_Z


   ! ------------ Read in the lateral matrix dimension. ---------------------------------------------

CALL ReadVar( UI, InFile, p_grid%NumGrid_Y, "the horizontal matrix dimension", "Horizontal grid-point matrix dimension")

   IF ( p_grid%NumGrid_Y < 2 )  THEN
      CALL TS_Abort ( 'The matrix must be >= 2x2.' )
   ENDIF

FormStr = "( I10 , 2X , 'Horizontal grid-point matrix dimension' )"
WRITE (US,FormStr)  p_grid%NumGrid_Y


   ! ------------ Read in the time step. ---------------------------------------------

CALL ReadVar( UI, InFile, p_grid%TimeStep, "the time step", "Time step [seconds]")

   IF ( p_grid%TimeStep <=  0.0 )  THEN
      CALL TS_Abort ( 'The time step must be greater than zero.' )
   ENDIF

FormStr = "( F10.3 , 2X , 'Time step [seconds]' )"
WRITE (US,FormStr)  p_grid%TimeStep


   ! ------------ Read in the analysis time. ---------------------------------------------

CALL ReadVar( UI, InFile, p_grid%AnalysisTime, "the analysis time", "Analysis time [seconds]")

   IF ( p_grid%AnalysisTime <=  0.0 )  THEN
      CALL TS_Abort ( 'The analysis time must be greater than zero.' )
   ENDIF

FormStr = "( F10.3 , 2X , 'Analysis time [seconds]' )"
WRITE (US,FormStr)  p_grid%AnalysisTime


   ! ------------ Read in the usable time. ---------------------------------------------
CALL ReadVar( UI, InFile, Line, "the usable output time", "Usable output time [seconds]")

   READ( Line, *, IOSTAT=IOS) UsableTime
   
   IF ( IOS /= 0 ) THEN ! Line didn't contian a number
      CALL Conv2UC( Line )
      IF ( TRIM(Line) == 'ALL' ) THEN
         Periodic   = .TRUE.
         UsableTime = p_grid%AnalysisTime
      ELSE
         CALL TS_Abort ( 'The usable output time must be a number greater than zero (or the string "ALL").' )
      END IF
   
   ELSE
      IF ( UsableTime <=  0.0 )  THEN
         CALL TS_Abort ( 'The usable output time must be a number greater than zero (or the string "ALL").' )
      ENDIF
      Periodic = .FALSE.
   END IF
   

FormStr = "( F10.3 , 2X , 'Usable output time [seconds]' )"
WRITE (US,FormStr)  UsableTime


   ! ------------ Read in the hub height. ---------------------------------------------

CALL ReadVar( UI, InFile, p_grid%HubHt, "the hub height", "Hub height [m]")

   IF ( p_grid%HubHt <=  0.0 )  THEN
      CALL TS_Abort ( 'The hub height must be greater than zero.' )
   ENDIF

FormStr = "( F10.3 , 2X , 'Hub height [m]' )"
WRITE (US,FormStr)  p_grid%HubHt


   ! ------------ Read in the grid height. ---------------------------------------------

CALL ReadVar( UI, InFile, p_grid%GridHeight, "the grid height", "Grid height [m]")

   IF ( 0.5*p_grid%GridHeight > p_grid%HubHt  )THEN
      CALL TS_Abort( 'The hub must be higher than half of the grid height.')
   ENDIF

FormStr = "( F10.3 , 2X , 'Grid height [m]' )"
WRITE (US,FormStr)  p_grid%GridHeight


   ! ------------ Read in the grid width. ---------------------------------------------

CALL ReadVar( UI, InFile, p_grid%GridWidth, "the grid width", "Grid width [m]")

   IF ( p_grid%GridWidth <=  0.0 )  THEN
      CALL TS_Abort ( 'The grid width must be greater than zero.' )
   ENDIF

FormStr = "( F10.3 , 2X , 'Grid width [m]' )"
WRITE (US,FormStr)  p_grid%GridWidth


   ! ***** Calculate the diameter of the rotor disk *****

RotorDiameter = MIN( p_grid%GridWidth, p_grid%GridHeight )


   ! ------------ Read in the vertical flow angle. ---------------------------------------------

CALL ReadVar( UI, InFile, VFlowAng, "the vertical flow angle", "Vertical flow angle [degrees]")

   IF ( ABS( VFlowAng ) > 45.0 )  THEN
      CALL TS_Abort ( 'The vertical flow angle must not exceed +/- 45 degrees.' )
   ENDIF

FormStr = "( F10.3 , 2X , 'Vertical flow angle [degrees]' )"
WRITE (US,FormStr)  VFlowAng


   ! ------------ Read in the horizontal flow angle. ---------------------------------------------

CALL ReadVar( UI, InFile, HFlowAng, "the horizontal flow angle", "Horizontal flow angle [degrees]")

FormStr = "( F10.3 , 2X , 'Horizontal flow angle [degrees]' )"
WRITE (US,FormStr)  HFlowAng


   !==========================================================
   ! Read the meteorological boundary conditions.
   !==========================================================

FormStr = "( // 'Meteorological Boundary Conditions:' / )"
WRITE (US,FormStr)
!READ (UI,'(/)')
CALL ReadCom( UI, InFile, "Meteorological Boundary Conditions Heading Line 1" )
CALL ReadCom( UI, InFile, "Meteorological Boundary Conditions Heading Line 2" )


   ! ------------ Read in the turbulence model. ---------------------------------------------

CALL ReadVar( UI, InFile, TurbModel, "the spectral model", "spectral model")

   TurbModel = ADJUSTL( TurbModel )
   CALL Conv2UC( TurbModel )

   SELECT CASE ( TRIM(TurbModel) )
      CASE ( 'IECKAI' )
         TMName = 'IEC Kaimal'
         SpecModel = SpecModel_IECKAI
      CASE ( 'IECVKM' )
         TMName = 'IEC von Karman'
         SpecModel = SpecModel_IECVKM
      CASE ( 'TIDAL' )
         TMName = 'Tidal Channel Turbulence'
         SpecModel = SpecModel_TIDAL
      CASE ( 'RIVER' )
         TMName = 'River Turbulence'
         SpecModel = SpecModel_RIVER
      CASE ( 'SMOOTH' )
         TMName = 'RISO Smooth Terrain'
         SpecModel = SpecModel_SMOOTH
      CASE ( 'WF_UPW' )
         TMName = 'NREL Wind Farm Upwind'
         SpecModel = SpecModel_WF_UPW
      CASE ( 'WF_07D' )
         TMName = 'NREL 7D Spacing Wind Farm'
         SpecModel = SpecModel_WF_07D
      CASE ( 'WF_14D' )
         TMName = 'NREL 14D Spacing Wind Farm'
         SpecModel = SpecModel_WF_14D
      CASE ( 'NONE'   )
         TMName = 'No fluctuating wind components'
         SpecModel = SpecModel_NONE
      CASE ( 'MODVKM' )
         TMName = 'Modified von Karman'
         SpecModel = SpecModel_MODVKM
      CASE ( 'API' )
         TMName = 'API'
         SpecModel = SpecModel_API
      CASE ( 'NWTCUP' )
         TMName = 'NREL National Wind Technology Center'
         SpecModel = SpecModel_NWTCUP
      CASE ( 'GP_LLJ' )
         TMName = 'Great Plains Low-Level Jet'
         SpecModel = SpecModel_GP_LLJ
      CASE ( 'USRVKM' )
         TMName = 'von Karman model with user-defined specifications'
         SpecModel = SpecModel_USRVKM
      CASE ( 'USRINP' )
         TMName = 'User-input '
         CALL GetUSRspec("UsrSpec.inp", ErrStat, ErrMsg)      ! bjj: before documenting, please fix this hard-coded name!
         SpecModel = SpecModel_USER
      CASE DEFAULT
!BONNIE: add the UsrVKM model to this list when the model is complete as well as USRINP
         CALL TS_Abort ( 'The turbulence model must be "IECKAI", "IECVKM", "SMOOTH",' &
                    //' "WF_UPW", "WF_07D", "WF_14D", "NWTCUP", "GP_LLJ", "TIDAL", "API", or "NONE".' )

      END SELECT  ! TurbModel
      IF (ErrStat >= AbortErrLev) RETURN

FormStr = "( 4X , A6 , 2X , '"//TRIM( TMName )//" spectral model' )"
WRITE (US,FormStr)  TurbModel

   ! ------------ Read in the IEC standard and edition numbers. ---------------------------------------------

CALL ReadVar( UI, InFile, Line, "the number of the IEC standard", "Number of the IEC standard")

   IF ( SpecModel == SpecModel_IECKAI .or. SpecModel == SpecModel_IECVKM .OR. SpecModel == SpecModel_API ) THEN  !bjj: SpecModel==SpecModel_MODVKM is not in the IEC standard

      CALL Conv2UC( LINE )

      IF ( (Line(1:1) == 'T') .OR.  (Line(1:1) == 'F') ) THEN
         CALL TS_Abort ( 'The number of the IEC standard must be either "1", "2", or "3"' &
                          // ' with an optional IEC 61400-1 edition number ("1-ED2").' )
      ENDIF

      TmpIndex = INDEX(Line, "-ED")

      IF ( TmpIndex > 0 ) THEN
         READ ( Line(TmpIndex+3:),*,IOSTAT=IOS ) p_IEC%IECedition

         CALL CheckIOS( IOS, InFile, 'the IEC edition number', NumType )

         IF ( p_IEC%IECedition < 1 .OR. p_IEC%IECedition > SIZE(IECeditionSTR) ) THEN
            CALL TS_Abort( 'Invalid IEC edition number.' )
         ENDIF

         Line = Line(1:TmpIndex-1)
      ELSE
         p_IEC%IECedition = 0
      ENDIF

      READ ( Line,*,IOSTAT=IOS ) p_IEC%IECstandard

      SELECT CASE ( p_IEC%IECstandard )
         CASE ( 1 )
            IF (p_IEC%IECedition < 1 ) THEN  ! Set up the default
               IF ( SpecModel == SpecModel_IECVKM .OR. SpecModel == SpecModel_USRVKM ) THEN
                  p_IEC%IECedition = 2   ! The von Karman model is not specified in edition 3 of the -1 standard
               ELSE
                  p_IEC%IECedition = 3
               ENDIF
            ELSE
               IF ( p_IEC%IECedition < 2 ) THEN
                  CALL TS_Abort( 'The IEC edition number must be 2 or 3' )
               ENDIF
            ENDIF

         CASE ( 2 )
               ! The scaling is the same as 64100-1, Ed. 2 with "A" or user-specified turbulence
            IF (p_IEC%IECedition < 1 ) THEN  ! Set up the default
               p_IEC%IECedition = 2
            ELSE
               CALL TS_Abort( ' The edition number cannot be specified for the 61400-2 standard.')
            ENDIF
            IECeditionSTR(p_IEC%IECedition) = 'IEC 61400-2 Ed. 2: 2005'

         CASE ( 3 )
               ! The scaling for 61400-3 (Offshore) is the same as 61400-1 except it has a different power law exponent
            IF (p_IEC%IECedition < 1 ) THEN  ! Set up the default

               IF ( SpecModel /= SpecModel_IECKAI ) THEN
                  CALL TS_Abort( ' The von Karman model is not valid for the 61400-3 standard.')
               ENDIF
               p_IEC%IECedition = 3   ! This is the edition of the -1 standard

            ELSE
               CALL TS_Abort( ' The edition number cannot be specified for the 61400-3 standard.')
            ENDIF
            IECeditionSTR(p_IEC%IECedition) = 'IEC 61400-3 Ed. 1: 2006'

         CASE DEFAULT
            CALL TS_Abort( 'The number of the IEC 61400-x standard must be 1, 2, or 3.')

      END SELECT

      FormStr = "( 7X, I3, 2X, 'IEC standard: ', A )"
      WRITE (US,FormStr)  p_IEC%IECstandard, TRIM(IECeditionSTR(p_IEC%IECedition))

   ELSE ! NOT IEC
      p_IEC%IECstandard = 0
      p_IEC%IECedition  = 0

      FormStr = "( 7X, A3, 2X, 'IEC standard' )"
      WRITE (US,FormStr)  'N/A'

   ENDIF ! IEC


   ! ------------ Read in the IEC turbulence characteristic. ---------------------------------------------

CALL ReadVar( UI, InFile, Line, "the IEC turbulence characteristic", "IEC turbulence characteristic")
!!$! -- begin block --
!!$! This block reads turbulence intensity, and stores it in the variable TurbIntH20.  This variable is not currently used, but will be soon for user-specified turbulence intensity for the HYDRO spectral models.
!!$! This code is copied from TurbSim.f90
!!$READ (Line,*,IOSTAT=IOS)  PerTurbInt ! This is to read the
!!$IECTurbC = ADJUSTL( Line )
!!$IF ( IOS /= 0 ) THEN
!!$   CALL Conv2UC(IECTurbC)
!!$   IF ( IECTurbC == 'A' ) THEN
!!$      TurbIntH20  = 0.16
!!$   ELSEIF ( IECTurbC == 'B' ) THEN
!!$      TurbIntH20  = 0.14
!!$   ELSEIF ( IECTurbC == 'C' ) THEN
!!$      TurbIntH20  = 0.12
!!$   ELSEIF ( IECTurbC == 'K' ) THEN
!!$      TurbIntH20  = 0.16
!!$   ELSE   ! We should never get here, but just to be complete...
!!$      !print *, IECTurbC
!!$      CALL TS_Abort( ' Invalid IEC turbulence characteristic.' )
!!$   ENDIF
!!$ELSE
!!$   TurbIntH20 = 0.01*PerTurbInt
!!$ENDIF
!!$! -- end block --



   IF ( SpecModel == SpecModel_IECKAI .or. SpecModel == SpecModel_IECVKM .OR. SpecModel == SpecModel_MODVKM .OR. SpecModel == SpecModel_API ) THEN

      READ (Line,*,IOSTAT=IOS) Line1

      CALL Conv2UC( Line1 )

      IF ( (Line1 == 'T') .OR.  (Line1 == 'F') ) THEN
         CALL TS_Abort ( 'The IEC turbulence characteristic must be either "A", "B", "C", or a real number.')
      ENDIF

      ! Check to see if the entry was a number.

      READ (Line,*,IOSTAT=IOS)  p_IEC%PerTurbInt

      IF ( IOS == 0 )  THEN

         ! Let's use turbulence value.

         p_IEC%NumTurbInp = .TRUE.
         FormStr    = "( F10.3 , 2X , 'Percent turbulence intensity, ', A )"
         WRITE (US,FormStr)  p_IEC%PerTurbInt, TRIM(IECeditionSTR(p_IEC%IECedition))

      ELSE

         ! Let's use one of the standard turbulence values (A or B or C).

         p_IEC%NumTurbInp = .FALSE.
         p_IEC%IECTurbC   = ADJUSTL( Line )
         CALL Conv2UC(  p_IEC%IECTurbC )

         SELECT CASE ( p_IEC%IECTurbC )
            CASE ( 'A' )
            CASE ( 'B' )
               IF ( p_IEC%IECstandard == 2 ) THEN
                  CALL TS_Abort (' The IEC 61400-2 turbulence characteristic must be either "A" or a real number.' )
               ENDIF
            CASE ( 'C' )
               IF ( p_IEC%IECstandard == 2 ) THEN
                  CALL TS_Abort (' The IEC 61400-2 turbulence characteristic must be either "A" or a real number.' )
               ELSEIF ( p_IEC%IECedition < 3 ) THEN
                  CALL TS_Abort (' The turbulence characteristic for '//TRIM(IECeditionSTR(p_IEC%IECedition) )// &
                                  ' must be either "A", "B", or a real number.' )
               ENDIF
            CASE DEFAULT
               CALL TS_Abort ( 'The IEC turbulence characteristic must be either "A", "B", "C", or a real number.' )
         END SELECT  ! IECTurbC

         FormStr  = "( 9X , A1 , 2X , 'IEC turbulence characteristic' )"
         WRITE (US,FormStr)  p_IEC%IECTurbC

      ENDIF

   ELSE  !not IECKAI and not IECVKM and not MODVKM

      Line = ADJUSTL( Line )
      CALL Conv2UC( Line )

      IF ( Line(1:6) == 'KHTEST' ) THEN
         KHtest = .TRUE.
         FormStr = "( 4X, A6, 2X, 'Kelvin-Helmholtz billow test case' )"
         WRITE (US,FormStr)  'KHTEST'


         IF ( SpecModel /= SpecModel_NWTCUP ) THEN
            CALL TS_Abort( 'The KH test can be used with the "NWTCUP" spectral model only.' )
         ENDIF

         IF ( .NOT. WrACT ) THEN
            CALL WRScr( ' Coherent turbulence time step files must be generated when using the "KHTEST" option.' )
            WRACT  = .TRUE.
         ENDIF

      ELSE
         FormStr = "( 7X, A3, 2X, 'IEC turbulence characteristic' )"
         WRITE (US,FormStr)  'N/A'
      ENDIF


         ! These variables are not used for non-IEC turbulence

      p_IEC%IECedition = 0
      p_IEC%NumTurbInp = .TRUE.
      p_IEC%PerTurbInt = 0.0

   ENDIF

   ! ------------ Read in the IEC wind turbulence type ---------------------------------------------

CALL ReadVar( UI, InFile, Line, "the IEC turbulence type", "IEC turbulence type")

   IF ( SpecModel == SpecModel_IECKAI .or. SpecModel == SpecModel_IECVKM .OR. SpecModel == SpecModel_API ) THEN

      CALL Conv2UC( Line )

      p_IEC%IECTurbE   = Line(1:1)

      ! Let's see if the first character is a number (for the ETM case)
      SELECT CASE ( p_IEC%IECTurbE )
         CASE ('1')
            p_IEC%Vref = 50.0
            Line = Line(2:)
         CASE ('2')
            p_IEC%Vref = 42.5
            Line = Line(2:)
         CASE ('3')
            p_IEC%Vref = 37.5
            Line = Line(2:)
         CASE DEFAULT
               ! There's no number at the start of the string so let's move on (it's NTM).
            p_IEC%Vref = -999.9
            p_IEC%IECTurbE = ' '
      END SELECT

      SELECT CASE ( TRIM( Line ) )
         CASE ( 'NTM'   )
            p_IEC%IEC_WindType = IEC_NTM
            p_IEC%IEC_WindDesc = 'Normal Turbulence Model'
         CASE ( 'ETM'   )
            p_IEC%IEC_WindType = IEC_ETM
            p_IEC%IEC_WindDesc = 'Extreme Turbulence Model'
         CASE ( 'EWM1'  )
            p_IEC%IEC_WindType = IEC_EWM1
            p_IEC%IEC_WindDesc = 'Extreme 1-Year Wind Speed Model'
         CASE ( 'EWM50' )
            p_IEC%IEC_WindType = IEC_EWM50
            p_IEC%IEC_WindDesc = 'Extreme 50-Year Wind Speed Model'
         CASE DEFAULT
            CALL TS_Abort ( ' Valid entries for the IEC wind turbulence are "NTM", "xETM", "xEWM1", or "xEWM50", '// &
                             'where x is the wind turbine class (1, 2, or 3).' )
      END SELECT

      IF ( p_IEC%IEC_WindType /= IEC_NTM ) THEN

         IF (p_IEC%IECedition /= 3 .OR. p_IEC%IECstandard == 2) THEN
            CALL TS_Abort ( ' The extreme turbulence and extreme wind speed models are available with '// &
                         'the IEC 61400-1 Ed. 3 or 61400-3 scaling only.')
         ENDIF

         IF (p_IEC%Vref < 0. ) THEN
            CALL TS_Abort ( ' A wind turbine class (1, 2, or 3) must be specified with the '// &
                         'extreme turbulence and extreme wind types. (i.e. "1ETM")')
         ENDIF

         IF ( p_IEC%NumTurbInp ) THEN
            CALL TS_Abort ( ' When the turbulence intensity is entered as a percent, the IEC wind type must be "NTM".' )
         ENDIF

      ELSE

         p_IEC%IECTurbE = ' '

      ENDIF

      FormStr = "( 4X, A6 , 2X , 'IEC ', A )"
      WRITE (US,FormStr)  TRIM(p_IEC%IECTurbE)//TRIM(Line), TRIM(p_IEC%IEC_WindDesc)

   ELSE
      p_IEC%IEC_WindType = IEC_NTM

      FormStr = "( A10 , 2X , 'IEC turbulence type' )"
      WRITE (US,FormStr)  'N/A'
   ENDIF

   ! ------------ Read in the ETM c parameter (IEC 61400-1, Ed 3: Section 6.3.2.3, Eq. 19) ----------------------
UseDefault = .TRUE.
p_IEC%ETMc = 2;
CALL ReadRVarDefault( UI, InFile, p_IEC%ETMc, "the ETM c parameter", 'IEC Extreme Turbulence Model (ETM) "c" parameter [m/s]', US, &
                      UseDefault, IGNORE=(p_IEC%IEC_WindType /= IEC_ETM ))

   IF ( p_IEC%ETMc <= 0. ) THEN
      CALL TS_Abort('The ETM "c" parameter must be a positive number');
   ENDIF

   ! ------------ Read in the wind profile type -----------------------------------------------------------------

UseDefault = .TRUE.         ! Calculate the default value
SELECT CASE ( SpecModel )
   CASE ( SpecModel_GP_LLJ )
      WindProfileType = 'JET'
   CASE ( SpecModel_IECKAI,SpecModel_IECVKM,SpecModel_MODVKM )
      WindProfileType = 'IEC'
   CASE ( SpecModel_TIDAL )
      WindProfileType = 'H2L'
   CASE ( SpecModel_USRVKM )
      WindProfileType = 'USR'
   CASE ( SpecModel_API )
      WindProfileType = 'API'  ! ADDED BY YG
   CASE DEFAULT
      WindProfileType = 'IEC'
END SELECT

CALL ReadCVarDefault( UI, InFile, WindProfileType, "the wind profile type", "Wind profile type", US, UseDefault ) !converts WindProfileType to upper case

      ! Make sure the variable is valid for this turbulence model

   SELECT CASE ( TRIM(WindProfileType) )
      CASE ( 'JET','J' )
         IF ( SpecModel /= SpecModel_GP_LLJ ) THEN
            CALL TS_Abort( 'The jet wind profile is available with the GP_LLJ spectral model only.')
         ENDIF
      CASE ( 'LOG','L' )
         IF (p_IEC%IEC_WindType /= IEC_NTM ) THEN
            CALL TS_Abort( ' The IEC turbulence type must be NTM for the logarithmic wind profile.' )
!bjj check that IEC_WindType == IEC_NTM for non-IEC
         ENDIF
      CASE ( 'PL',  'P' )
      CASE ( 'H2L', 'H' )
         IF ( SpecModel /= SpecModel_TIDAL ) THEN
            CALL TS_Abort(  'The "H2L" mean profile type should be used only with the "TIDAL" spectral model.' )
         ENDIF
      CASE ( 'IEC', 'N/A' )
      CASE ( 'USR', 'U' )
      CASE ( 'API', 'A' )   ! ADDED BY Y.GUO
      CASE DEFAULT
         CALL TS_Abort( 'The wind profile type must be "JET", "LOG", "PL", "IEC", "USR", "H2L", or default.' )
   END SELECT

   IF ( SpecModel == SpecModel_TIDAL .AND. TRIM(WindProfileType) /= "H2L" ) THEN
      WindProfileType = 'H2L'
      CALL TS_Warn  ( 'Overwriting wind profile type to "H2L" for the "TIDAL" spectral model.', .TRUE.)
   ENDIF

   IF ( KHTest ) THEN
      IF ( TRIM(WindProfileType) /= 'IEC' .AND. TRIM(WindProfileType) /= 'PL' ) THEN
         WindProfileType = 'IEC'
         CALL TS_Warn  ( 'Overwriting wind profile type for the KH test.', .FALSE.)
      ENDIF
   ENDIF


   ! ------------ Read in the height for the reference wind speed. ---------------------------------------------

CALL ReadVar( UI, InFile, RefHt, "the reference height", "Reference height [m]")

   IF ( RefHt <=  0.0 .AND. WindProfileType(1:1) /= 'U' )  THEN
      CALL TS_Abort ( 'The reference height must be greater than zero.' )
   ENDIF

FormStr = "( F10.3 , 2X , 'Reference height [m]' )"
WRITE (US,FormStr)  RefHt
H_ref = RefHt ! Define the variable H_ref, for later use in HYDRO spectral models (RefHt gets modified later in the code)

   ! ------------ Read in the reference wind speed. -----------------------------------------------------

UseDefault = .FALSE.
URef       = -999.9

! If we specify a Ustar (i.e. if Ustar /= "default") then we can enter "default" here,
! otherwise, we get circular logic...

   CALL ReadRVarDefault( UI, InFile, URef, "the reference wind speed", "Reference wind speed [m/s]", US, UseDefault, &
                     IGNORE=(p_IEC%IEC_WindType == IEC_EWM1 .OR. p_IEC%IEC_WindType == IEC_EWM50 .OR. WindProfileType(1:1) == 'U') )

   NumUSRz = 0  ! initialize the number of points in a user-defined wind profile

   IF ( ( WindProfileType(1:1) /= 'J' .OR. .NOT. UseDefault) .AND. &
        (p_IEC%IEC_WindType /= IEC_EWM1 .AND. p_IEC%IEC_WindType /= IEC_EWM50 .AND. WindProfileType(1:1) /= 'U') ) THEN
      IF ( URef <=  0.0 )  THEN
         CALL TS_Abort ( 'The reference wind speed must be greater than zero.' )
      ENDIF

   ELSEIF ( WindProfileType(1:1) == 'U' ) THEN
      RefHt = p_grid%HubHt
      CALL GetUSR( UI, InFile, 37, ErrStat, ErrMsg ) !Read the last several lines of the file, then return to line 37
      IF (ErrStat >= AbortErrLev) RETURN
      URef = getWindSpeed(URef, RefHt, RefHt, RotorDiameter, PROFILE_TYPE=WindProfileType) !This is UHub
   ENDIF   ! Otherwise, we're using a Jet profile with default wind speed (for now it's -999.9)


      ! We need to save the reference wind speed for the API model.

   IF ( WindProfileType(1:3) == 'API' )  THEN
      U_Ref = URef
   END IF


   ! ------------ Read in the jet height -------------------------------------------------------------

UseDefault = .FALSE.
ZJetMax    = -999.9

CALL ReadRVarDefault( UI, InFile, ZJetMax, "the jet height", "Jet height [m]", US, UseDefault, IGNORE=WindProfileType(1:1) /= 'J')

   IF ( WindProfileType(1:1) == 'J' .AND. .NOT. UseDefault ) THEN
      IF ( ZJetMax <  70.0 .OR. ZJetMax > 490.0 )  THEN
         CALL TS_Abort ( 'The height of the maximum jet wind speed must be between 70 and 490 m.' )
      ENDIF
   ENDIF


   ! ------------ Read in the power law exponent, PLExp ---------------------------------------------

SELECT CASE ( SpecModel )
   CASE (SpecModel_WF_UPW, SpecModel_WF_07D, SpecModel_WF_14D, SpecModel_NWTCUP)
      IF ( KHtest ) THEN
         UseDefault = .TRUE.
         PLExp      = 0.3
      ELSE
         UseDefault = .FALSE.             ! This case needs the Richardson number to get a default
         PLExp      = 0.
      ENDIF

   CASE DEFAULT
      UseDefault = .TRUE.
      PLExp      = PowerLawExp( RICH_NO )  ! These cases do not use the Richardson number to get a default

END SELECT
getPLExp = .NOT. UseDefault

CALL ReadRVarDefault( UI, InFile, PLExp, "the power law exponent", "Power law exponent", US, UseDefault, &
               IGNORE=( ((INDEX( 'JLUH', WindProfileType(1:1) ) > 0.)) .OR. &
                         p_IEC%IEC_WindType == IEC_EWM1 .OR. p_IEC%IEC_WindType == IEC_EWM50 ) )

   IF ( .NOT. UseDefault ) THEN  ! We didn't use a default (we entered a number on the line)
      getPLExp = .FALSE.

      IF ( KHtest ) THEN
         IF ( PLExp /= 0.3 ) THEN
            PLExp = 0.3
            CALL TS_Warn  ( 'Overwriting the power law exponent for KH test.', .FALSE. )
         ENDIF
      ENDIF
   ENDIF


   ! ------------ Read in the surface roughness length, Z0 (that's z-zero) ---------------------------------------------

UseDefault = .TRUE.
SELECT CASE ( SpecModel )
   CASE (SpecModel_SMOOTH)
      Z0 = 0.010
   CASE (SpecModel_GP_LLJ )
      Z0 = 0.005
   CASE (SpecModel_WF_UPW )
      Z0 = 0.018
   CASE (SpecModel_NWTCUP )
      Z0 = 0.021
   CASE (SpecModel_WF_07D )
      Z0 = 0.233
   CASE (SpecModel_WF_14D )
      Z0 = 0.064
   CASE DEFAULT !IEC values
      Z0 = 0.030 ! Represents smooth, homogenous terrain
END SELECT

CALL ReadRVarDefault( UI, InFile, Z0, "the roughness length", "Surface roughness length [m]", US, UseDefault, &
                       IGNORE=SpecModel==SpecModel_TIDAL)

   IF ( Z0 <= 0.0 ) THEN
      CALL TS_Abort ( 'The surface roughness length must be a positive number or "default".')
   ENDIF


   !=================================================================================
   ! Read the meteorological boundary conditions for non-IEC models. !
   !=================================================================================

IF ( SpecModel /= SpecModel_IECKAI .AND. SpecModel /= SpecModel_IECVKM .AND. SpecModel /= SpecModel_API ) THEN  ! Modified by Y.Guo

   FormStr = "( // 'Non-IEC Meteorological Boundary Conditions:' / )"
   WRITE (US,FormStr)
   !READ (UI,'(/)')
   CALL ReadCom( UI, InFile, "Non-IEC Meteorological Boundary Conditions Heading Line 1" )
   CALL ReadCom( UI, InFile, "Non-IEC Meteorological Boundary Conditions Heading Line 2" )



      ! ------------ Read in the site latitude, LATITUDE. ---------------------------------------------

   UseDefault = .TRUE.
   Latitude   = 45.0

   CALL ReadRVarDefault( UI, InFile, Latitude, "the site latitude", "Site latitude [degrees]", US, UseDefault)

      IF ( ABS(Latitude) < 5.0 .OR. ABS(Latitude) > 90.0 ) THEN
         CALL TS_Abort( 'The latitude must be between -90 and 90 degrees but not between -5 and 5 degrees.' )
      ENDIF

   Fc = 2.0 * Omega * SIN( ABS(Latitude*D2R) )  ! Calculate Coriolis parameter from latitude

ELSE

   Latitude = 0.0                               !Not used in IEC specs
   Fc = 0.0

ENDIF    ! Not IECKAI and Not IECVKM

IF ( SpecModel /= SpecModel_IECKAI .AND. &
     SpecModel /= SpecModel_IECVKM .AND. &
     SpecModel /= SpecModel_MODVKM .AND. &
     SpecModel /= SpecModel_API    ) THEN


      ! ------------ Read in the gradient Richardson number, RICH_NO. ---------------------------------------------

   CALL ReadVar( UI, InFile, RICH_NO, "the gradient Richardson number", "Gradient Richardson number")

   IF ( KHtest ) THEN
      IF ( RICH_NO /= 0.02 ) THEN
         RICH_NO = 0.02
         CALL TS_Warn ( 'Overwriting the Richardson Number for KH test.', .FALSE. )
      ENDIF
   ENDIF

   IF ( SpecModel == SpecModel_USER .OR. SpecModel == SpecModel_USRVKM ) THEN
      IF ( RICH_NO /= 0.0 ) THEN
         RICH_NO = 0.0
         CALL TS_Warn ( 'Overwriting the Richardson Number for the '//TRIM(TurbModel)//' model.', .FALSE. )
      ENDIF
   ENDIF

   IF ( SpecModel == SpecModel_TIDAL ) THEN
      WRITE (US,"( A10, 2X, A )")  "N/A", 'Gradient Richardson number'
      RICH_NO = 0
   ELSE
      FormStr = "( F10.3 , 2X , 'Gradient Richardson number' )"
      WRITE (US,FormStr)  RICH_NO
   END IF

      ! ***** Calculate M-O z/L parameter   :  z/L is a number in (-inf, 1] *****

   IF ( SpecModel == SpecModel_NWTCUP ) THEN
         ! Calculate disk averaged Z/L from turbine layer Ri for NWTC/LIST experiment

      RICH_NO = MIN( MAX( RICH_NO, REAL(-1.0,ReKi) ), REAL(1.0,ReKi) )  ! Ensure that: -1 <= RICH_NO <= 1

      IF ( RICH_NO <= -0.1 ) THEN
         ZL = -0.254 + 1.047*RICH_NO
      ELSEIF ( RICH_NO < 0 ) THEN
         ZL = 10.369*RICH_NO/(1.0 - 19.393*RICH_NO)
      ELSE  !( RICH_NO < 0.155 ) THEN
         ZL = 2.535*MIN( RICH_NO, REAL(0.155, ReKi) ) / (1.0 - 6.252*MIN( RICH_NO, REAL(0.155,ReKi) ))
      ENDIF


   ELSEIF (SpecModel == SpecModel_GP_LLJ) THEN

      RICH_NO = MIN( MAX( RICH_NO, REAL(-1.0,ReKi) ), REAL(1.0,ReKi) )  ! Ensure that: -1 <= RICH_NO <= 1

      IF ( RICH_NO <= -0.1 ) THEN
         ZL = -0.047 + 1.054*RICH_NO
      ELSEIF ( RICH_NO < 0 ) THEN
         ZL = 2.213*RICH_NO/(1.0 - 4.698*RICH_NO)
      ELSE  !( RICH_NO < 0.1367 ) THEN
         ZL = 3.132*MIN( RICH_NO, REAL(0.1367,ReKi) ) / (1.0 - 6.762*MIN( RICH_NO, REAL(0.1367,ReKi) ))
      ENDIF

   ELSE ! see Businger, J.A.; Wyngaard, J.C.; Izumi, Y.; Bradley, E.F. (1971). "Flux-Profile Relationships in the Atmospheric Surface Layer." Journal of the Atmospheric Sciences (28); pp.181-189.

      IF ( RICH_NO <= 0.0 ) THEN
         ZL = RICH_NO
         !PhiM = (1.0 - 16.0*ZL)**-0.25
      ELSEIF ( RICH_NO < 0.16667 ) THEN
         ZL = MIN( RICH_NO / ( 1.0 - 5.0*RICH_NO ), REAL(1.0,ReKi) )  ! The MIN() will take care of rounding issues.
         !PhiM = (1.0 + 5.0*ZL)
      ELSE
         ZL = 1.0
      ENDIF

   ENDIF !SpecModels

   ZL = MIN( ZL, REAL(1.0,ReKi) )

      ! ***** Calculate M-O length scale, L [meters] *****
      ! L should be constant in the surface layer

   IF ( ZL /= 0.0 ) THEN
        L = p_grid%HubHt / ZL ! Since ZL is the average ZL over the rotor disk, we should use HubHt to estimate L instead
   ELSE
      L = HUGE( L )
   ENDIF

      ! ***** Calculate power law exponent, if needed *****

   IF ( getPLExp ) THEN
      PLExp = PowerLawExp( RICH_NO )
   ENDIF

      ! ------------ Read in the shear/friction velocity, Ustar, first calculating UstarDiab ------------------------

         ! Set up the heights for the zl- and ustar-profile averages across the rotor disk
      TmpZary  = (/ p_grid%HubHt-RotorDiameter/2., p_grid%HubHt, p_grid%HubHt+RotorDiameter/2. /)
      IF (TmpZary(3) .GE. profileZmin .AND. TmpZary(1) .LE. profileZmax ) THEN  !set height limits so we don't extrapolate too far
         DO TmpIndex = 1,3
            TmpZary(TmpIndex) = MAX( MIN(TmpZary(TmpIndex), profileZmax), profileZmin)
         ENDDO
      ENDIF


   UseDefault = .TRUE.

   UstarDiab  = getUstarDiab(URef, RefHt)
   Ustar      = UstarDiab

   SELECT CASE ( SpecModel )

      CASE (SpecModel_WF_UPW)

         IF ( ZL < 0.0 ) THEN
            Ustar = 1.162 * UstarDiab**( 2.0 / 3.0 )
         ELSE ! Include the neutral case to avoid strange discontinuities
            Ustar = 0.911 * UstarDiab**( 2.0 / 3.0 )
         ENDIF

      CASE ( SpecModel_WF_07D, SpecModel_WF_14D  )

         IF ( ZL < 0.0 ) THEN
            Ustar = 1.484 * UstarDiab**( 2.0 / 3.0 )
         ELSE ! Include the neutral case with the stable one to avoid strange discontinuities
            Ustar = 1.370 * UstarDiab**( 2.0 / 3.0 )
         ENDIF

      CASE (SpecModel_GP_LLJ )
         TmpUstarD  = -1.0
         IF ( URef < 0 ) THEN ! (1) We can't get a wind speed because Uref was default
            UseDefault = .FALSE. ! We'll calculate the default value later
         ELSE
            Ustar = 0.17454 + 0.72045*UstarDiab**1.36242
         ENDIF

      CASE ( SpecModel_NWTCUP  )
         UStar = 0.2716 + 0.7573*UstarDiab**1.2599

      CASE ( SpecModel_TIDAL , SpecModel_RIVER )
         ! Use a constant drag coefficient for the HYDRO spectral models.
         UStar = Uref*0.05 ! This corresponds to a drag coefficient of 0.0025.
         !UStar = Uref*0.04 ! This corresponds to a drag coefficient of 0.0016.

   END SELECT

   CALL ReadRVarDefault( UI, InFile, UStar, "the friction or shear velocity", "Friction or shear velocity [m/s]", US, UseDefault )

   IF ( Uref < 0.0 .AND. UseDefault ) THEN  ! This occurs if "default" was entered for both GP_LLJ wind speed and UStar
      CALL TS_Abort( 'The reference wind speed and friction velocity cannot both be "default."')
   ELSEIF (UStar <= 0) THEN
      CALL TS_Abort( 'The friction velocity must be a positive number.')
   ENDIF


         ! ***** Calculate wind speed at hub height *****

   IF ( WindProfileType(1:1) == 'J' ) THEN
      IF ( ZJetMax < 0 ) THEN ! Calculate a default value
         ZJetMax = -14.820561*Rich_No + 56.488123*zl + 166.499069*uStar + 188.253377
         ZJetMax = 1.9326*ZJetMax - 252.7267  ! Correct with the residual

         CALL RndJetHeight( p_RandNum, OtherSt_RandNum, tmp ) ! Add a random amount

         ZJetMax = MIN( MAX(ZJetMax + tmp, REAL(70.,ReKi) ), REAL(490.,ReKi) )
      ENDIF

      IF ( URef < 0 ) THEN ! Calculate a default value

         UJetMax = MAX( -21.5515 + 6.6827*LOG(ZJetMax), REAL(5.0,ReKi) ) !Jet max must be at least 5 m/s (occurs ~50 m); shouldn't happen, but just in case....

         CALL Rnd3ParmNorm( p_RandNum, OtherSt_RandNum, tmp, REAL(0.1076,ReKi), REAL(-0.1404,ReKi), REAL(3.6111,ReKi),  REAL(-15.,ReKi), REAL(20.,ReKi) )

         IF (UJetMax + tmp > 0 ) UJetMax = UJetMax + tmp

         CALL GetChebCoefs( UJetMax, ZJetMax ) ! These coefficients are a function of UJetMax, ZJetMax, RICH_NO, and uStar

         URef = getWindSpeed(UJetMax, ZJetMax, RefHt, RotorDiameter, PROFILE_TYPE=WindProfileType)

      ELSE
         CALL GetChebCoefs(URef, RefHt)
      ENDIF

   ENDIF !Jet wind profile

   UHub = getWindSpeed(URef, RefHt, p_grid%HubHt, RotorDiameter, PROFILE_TYPE=WindProfileType)

         ! ***** Get uStar- and zl-profile values, if required, and determine offsets *****
      IF ( TmpUstarD == 0.0 ) THEN
         TmpUstarD = Ustar
      ELSE
         IF ( TmpUstarD < 0.0 ) THEN  ! If > 0, we've already calculated these things...
            UstarDiab = getUstarDiab(URef, RefHt) !bjj: is this problematic for anything else?
            TmpUary   = getWindSpeed(URef, RefHt, TmpZary, RotorDiameter, PROFILE_TYPE=WindProfileType)
            TmpUstar  = getUstarARY( TmpUary,     TmpZary )
            TmpUstarD = SUM(TmpUstar) / SIZE(TmpUstar)    ! The average of those 3 points
         ENDIF
         UstarOffset = Ustar - TmpUstarD
         TmpUstar(:) = TmpUstar(:) + UstarOffset
      ENDIF

      TmpZLary = getZLARY(TmpUary, TmpZary)
      zlOffset = zL - SUM(TmpZLary) / SIZE(TmpUstar)


      ! ------------- Read in the mixing layer depth, ZI ---------------------------------------------

   UseDefault = .TRUE.
   IF ( ZL >= 0.0 .AND. SpecModel /= SpecModel_GP_LLJ ) THEN  !We must calculate ZI for stable GP_LLJ. z/L profile can change signs so ZI must be defined for spectra.
      ZI = 0.0
   ELSE
      IF ( Ustar < UstarDiab ) THEN
         ZI = ( 0.04 * Uref ) / ( 1.0E-4 * LOG10( RefHt / Z0 ) )  !for "very" windy days
      ELSE
         !Should Wind Farm models use the other definition since that was what was used in creating those models?
         ZI = Ustar / (6.0 * Fc)
      ENDIF
   ENDIF

   CALL ReadRVarDefault( UI, InFile, ZI, "the mixing layer depth", "Mixing layer depth [m]", US, UseDefault, &
                                                                  IGNORE=(ZL>=0. .AND. SpecModel /= SpecModel_GP_LLJ) )

      IF ( ( ZL < 0.0 ) .AND. ( ZI <= 0.0 ) ) THEN
         CALL TS_Abort ( 'The mixing layer depth must be a positive number for unstable flows.')
      ENDIF


      ! Get the default mean Reynolds stresses
   UWskip     = .FALSE.
   UVskip     = .FALSE.
   VWskip     = .FALSE.
   PC_UW      = TmpUstar(2)**2   ! Used only for GP-LLJ in GetDefaultRS()

   CALL GetDefaultRS(  PC_UW, PC_UV, PC_VW )


       ! ----------- Read in the mean hub u'w' Reynolds stress, PC_UW ---------------------------------------------
   UseDefault = .TRUE.
   CALL ReadRVarDefault( UI, InFile, PC_UW, "the mean hub u'w' Reynolds stress", &
                                            "Mean hub u'w' Reynolds stress", US, UseDefault, IGNORESTR = UWskip )
   IF (.NOT. UWskip) THEN
      TmpUstarD = ( TmpUstar(1)- 2.0*TmpUstar(2) + TmpUstar(3) )

      IF ( TmpUstarD /= 0.0 ) THEN
         UstarSlope  = 3.0*(Ustar -  SQRT( ABS(PC_UW) ) ) / TmpUstarD
         UstarOffset = SQRT( ABS(PC_UW) ) - UstarSlope*(TmpUstar(2) - UstarOffset)
      ELSE
         UstarSlope  = 0.0
         UstarOffset = SQRT( ABS(PC_UW) )
      ENDIF
   ENDIF

      ! ------------ Read in the mean hub u'v' Reynolds stress, PC_UV ---------------------------------------------
   UseDefault = .TRUE.
   CALL ReadRVarDefault( UI, InFile, PC_UV, "the mean hub u'v' Reynolds stress", &
                                            "Mean hub u'v' Reynolds stress", US, UseDefault, IGNORESTR = UVskip )

      ! ------------ Read in the mean hub v'w' Reynolds stress, PC_VW ---------------------------------------------
   UseDefault = .TRUE.
   CALL ReadRVarDefault( UI, InFile, PC_VW, "the mean hub v'w' Reynolds stress", &
                                            "Mean hub v'w' Reynolds stress", US, UseDefault, IGNORESTR = VWskip )

      ! ------------ Read in the u component coherence decrement (coh-squared def), InCDec(1) = InCDecU ------------
   CALL GetDefaultCoh( UHub, p_grid%HubHt )

   UseDefault = .TRUE.
   InCVar(1) = InCDec(1)
   InCVar(2) = InCohB(1)

   CALL ReadRAryDefault( UI, InFile, InCVar, "the u-component coherence parameters", &
                                             "u-component coherence parameters", US, UseDefault)
   InCDec(1) = InCVar(1)
   InCohB(1) = InCVar(2)

      IF ( InCDec(1) <= 0.0 ) THEN
         CALL TS_Abort ( 'The u-component coherence decrement must be a positive number.')
      ENDIF

      ! ------------ Read in the v component coherence decrement (coh-squared def), InCDec(2) = InCDecV ----------

   UseDefault = .TRUE.
   InCVar(1) = InCDec(2)
   InCVar(2) = InCohB(2)

   CALL ReadRAryDefault( UI, InFile, InCVar, "the v-component coherence parameters",  &
                                             "v-component coherence parameters", US, UseDefault)

   InCDec(2) = InCVar(1)
   InCohB(2) = InCVar(2)

      IF ( InCDec(2) <= 0.0 ) THEN
         CALL TS_Abort ( 'The v-component coherence decrement must be a positive number.')
      ENDIF

      ! ------------ Read in the w component coherence decrement (coh-squared def), InCDec(3) = InCDecW -------

   UseDefault = .TRUE.
   InCVar(1) = InCDec(3)
   InCVar(2) = InCohB(3)

   CALL ReadRAryDefault( UI, InFile, InCVar, "the w-component coherence parameters", &
                                             "w-component coherence parameters", US, UseDefault)

   InCDec(3) = InCVar(1)
   InCohB(3) = InCVar(2)

      IF ( InCDec(3) <= 0.0 ) THEN
         CALL TS_Abort ( 'The w-component coherence decrement must be a positive number.')
      ENDIF

         ! ------------ Read in the coherence exponent, COHEXP -----------------------------------

   UseDefault = .TRUE.
   CohExp     = 0.0    ! was 0.25
   CALL ReadRVarDefault( UI, InFile, CohExp, "the coherence exponent", "Coherence exponent", US, UseDefault)

      IF ( COHEXP < 0.0 ) THEN
         CALL TS_Abort ( 'The coherence exponent must be non-negative.')
      ENDIF


         !=================================================================================
         ! Read the Coherent Turbulence Scaling Parameters, if necessary.  !
         !=================================================================================

   IF ( WrACT ) THEN

      FormStr = "( // 'Coherent Turbulence Scaling Parameters:' / )"
      WRITE (US,FormStr)
      !READ (UI,'(/)')                        ! Read header line
      CALL ReadCom( UI, InFile, "Coherent Turbulence Scaling Parameters Heading Line 1" )
      CALL ReadCom( UI, InFile, "Coherent Turbulence Scaling Parameters Heading Line 2" )


         ! ------------ Read the name of the path containg event file definitions, CTEventPath --------------------------

      CALL ReadVar( UI, InFile, p_CohStr%CTEventPath, "the path of the coherent turbulence events", "Coherence events path")

      IF ( LEN( TRIM(p_CohStr%CTEventPath) ) <= 10 )  THEN
         FormStr = "( A10 , 2X , 'Name of the path containing the coherent turbulence data files' )"
      ELSE
         FormStr = "( A, /, 12X , 'Name of the path containing the coherent turbulence data files' )"
      ENDIF

      WRITE (US,FormStr)  TRIM(p_CohStr%CTEventPath)

      CALL ReadVar( UI, InFile, Line, "the event file type", "Event file type")


      IF ( KHtest ) THEN

         p_CohStr%CText = 'les'
         p_CohStr%CTEventFile = TRIM(p_CohStr%CTEventPath)//PathSep//'Events.xtm'

         CALL WrScr( ' LES events will be used for the KH test.' )

      ELSE

         p_CohStr%CText = Line  !This will preserve the case formatting, in case it matters.

         CALL Conv2UC( Line )

         IF (Line(1:6) == "RANDOM") THEN
             CALL RndUnif( p_RandNum, OtherSt_RandNum, tmp )

             IF ( tmp <= 0.5 ) THEN
                 p_CohStr%CText = 'les'
             ELSE
                 p_CohStr%CText = 'dns'
             ENDIF
         ENDIF

         p_CohStr%CTEventFile = TRIM(p_CohStr%CTEventPath)//PathSep//'Events.'//TRIM(p_CohStr%CText)

      ENDIF

      FormStr = "( 7X, A3, 2X, 'Type of coherent turbulence data files' )"

      WRITE (US,FormStr) TRIM(p_CohStr%CText)


         ! ------------ Read the Randomization Flag, Randomize -----------------------------------

      CALL ReadVar( UI, InFile, Randomize, "the randomization flag", "Randomize CT Scaling")


      IF ( KHtest ) THEN
         Randomize = .FALSE.
         CALL WrScr( ' Billow will cover rotor disk for KH test. ' )
      ENDIF

      FormStr = "( L10 , 2X , 'Randomize the disturbance scale and location?' )"
      WRITE (US,FormStr)  Randomize


         ! ------------ Read the Disturbance Scale, DistScl ---------------------------------------------

      CALL ReadVar( UI, InFile, p_CohStr%DistScl, "the disturbance scale", "Disturbance scale")


         ! ------------ Read the Lateral Fractional Location of tower centerline in wave, CTLy ----------

      CALL ReadVar( UI, InFile, p_CohStr%CTLy, "the fractional location of tower centerline from right", &
                           "Location of tower centerline")

         ! ------------ Read the Vertical Fraction Location of hub in wave, CTLz ------------------------

      CALL ReadVar( UI, InFile, p_CohStr%CTLz, "the fractional location of hub height from the bottom", &
                     "Location of hub height")

      IF ( KHtest ) THEN
            p_CohStr%DistScl = 1.0
            p_CohStr%CTLy    = 0.5
            p_CohStr%CTLz    = 0.5
      ELSEIF ( Randomize ) THEN

         CALL RndUnif( p_RandNum, OtherSt_RandNum, tmp )

            ! Assume a 75% chance of coherent turbulence being the size of the rotor
            ! If the rotor is too small, assume a 100% chance.
            ! If the turbulence is not the size of the rotor, assume it's half the size
            ! of the disk, with equal probablilty of being in the lower or upper half.

         IF ( tmp > 0.25 .OR. RotorDiameter <= 30.0 ) THEN

            p_CohStr%DistScl = 1.0
            p_CohStr%CTLy    = 0.5
            p_CohStr%CTLz    = 0.5

         ELSE

            p_CohStr%DistScl = 0.5
            p_CohStr%CTLy    = 0.5

            IF ( tmp < 0.125 ) THEN
               p_CohStr%CTLz = 0.0 ! The hub is on the bottom of the dataset (i.e. the billow is on the top of the disk)
            ELSE
               p_CohStr%CTLz = 1.0 ! The hub is on the top of the dataset
            ENDIF

         ENDIF

      ELSE  !Don't randomize:

         IF ( p_CohStr%DistScl < 0.0 ) THEN
            CALL TS_Abort ('The disturbance scale must be a positive.')
         ELSEIF ( RotorDiameter <= 30.0 .AND. p_CohStr%DistScl < 1.0 ) THEN
            CALL TS_Abort ('The disturbance scale must be at least 1.0 for rotor diameters less than 30.')
         ELSEIF ( RotorDiameter*p_CohStr%DistScl <= 15.0  ) THEN
            CALL TS_Abort ('The coherent turbulence must be greater than 15 meters in height.  '//&
                        'Increase the rotor diameter or the disturbance scale. ')
         ENDIF

      ENDIF


      FormStr = "( F10.3 , 2X , 'Disturbance scale (ratio of wave height to rotor diameter)' )"
      WRITE (US,FormStr)  p_CohStr%DistScl

      FormStr = "( F10.3 , 2X , 'Fractional location of tower centerline from right' )"
      WRITE (US,FormStr)  p_CohStr%CTLy

      FormStr = "( F10.3 , 2X , 'Fractional location of hub height from the bottom of the dataset' )"
      WRITE (US,FormStr)  p_CohStr%CTLz

         ! ---------- Read the Minimum event start time, CTStartTime --------------------------------------------

      CALL ReadVar( UI, InFile, p_CohStr%CTStartTime, "the minimum start time for coherent structures", "CTS Start Time")

      p_CohStr%CTStartTime = MAX( p_CohStr%CTStartTime, REAL(0.0,ReKi) ) ! A Negative start time doesn't really make sense...

      FormStr = "( F10.3 , 2X , 'Minimum start time for coherent structures [seconds]' )"
      WRITE (US,FormStr)  p_CohStr%CTStartTime

   ENDIF    ! WrACT


ELSE  ! IECVKM or IECKAI models

   RICH_NO = 0.0                       ! Richardson Number in neutral conditions
   COHEXP  = 0.0                       ! Coherence exponent

   IF ( p_IEC%IECedition == 3 ) THEN
      IncDec = (/ 12.00, 0.0, 0.0 /)   ! u-, v-, and w-component coherence decrement for IEC Ed. 3
   ELSE
      IncDec = (/  8.80, 0.0, 0.0 /)   ! u-, v-, and w-component coherence decrement
   ENDIF


      ! The following variables are not used in the IEC calculations

   ZL      = 0.0                       ! M-O z/L parameter
   L       = 0.0                       ! M-O length scale
   Ustar   = 0.0                       ! Shear or friction velocity
   ZI      = 0.0                       ! Mixing layer depth
   PC_UW   = 0.0                       ! u'w' x-axis correlation coefficient
   PC_UV   = 0.0                       ! u'v' x-axis correlation coefficient
   PC_VW   = 0.0                       ! v'w' x-axis correlation coefficient

   IF ( p_IEC%NumTurbInp .AND. p_IEC%PerTurbInt == 0.0 ) THEN    ! This will produce constant winds, instead of an error when the transfer matrix is singular
      TurbModel = 'NONE'
      SpecModel = SpecModel_NONE
   ENDIF

      ! Calculate wind speed at hub height

   UHub    = getWindSpeed(URef, RefHt, p_grid%HubHt, RotorDiameter, PROFILE_TYPE=WindProfileType)


ENDIF


   ! Done reading the input file.

CLOSE (UI)


RETURN
END SUBROUTINE GetInput
!=======================================================================
SUBROUTINE GetFiles

  ! This subroutine is used to open the summary output file.

USE              TSMods

IMPLICIT         NONE


CALL GetRoot( InFile, RootName )


   ! Open summary file.

CALL OpenFOutFile( US, TRIM( RootName )//'.sum' ) ! Formatted output file


   ! Write the program name and version, date and time into the summary file.

   ! Let's make sure the binary file and the full-field file have the same date and time.
DescStr = 'generated by '//TRIM( ProgName )//TRIM( ProgVer )//' on '//CurDate()//' at '//CurTime()//'.'

FormStr = "( / 'This summary file was ', A / )"
WRITE (US,FormStr)  TRIM(DescStr)

   ! Capitalize the first letter of the string.

DescStr = 'This full-field file was '//TRIM(DescStr)


RETURN
END SUBROUTINE GetFiles
!=======================================================================
SUBROUTINE GetUSR(U_in, FileName, NLines, ErrStat, ErrMsg)

   USE                                   TSMods, ONLY: p_grid
   USE                                   TSMods, ONLY: L_USR             ! von Karman length scale
   USE                                   TSMods, ONLY: NumUSRz           ! Number of user-defined heights for the profiles
   USE                                   TSMods, ONLY: Sigma_USR         ! standard deviation profile
   USE                                   TSMods, ONLY: StdScale          ! scaling for standard deviation profile
   USE                                   TSMods, ONLY: TurbModel, SpecModel         ! Type of turbulence model
   USE                                   TSMods, ONLY: U_USR             ! wind speed profile
   USE                                   TSMods, ONLY: WindDir_USR       ! wind direction profile
   USE                                   TSMods, ONLY: Z_USR             ! heights corresponding to the profiles
   IMPLICIT                              NONE

   INTEGER(IntKi),                  intent(  out) :: ErrStat                         ! Error level
   CHARACTER(*),                    intent(  out) :: ErrMsg                          ! Message describing error
   
   INTEGER(IntKi)                                 :: ErrStat2                        ! Error level (local)
   CHARACTER(MaxMsgLen)                           :: ErrMsg2                         ! Message describing error (local)
   
   
   CHARACTER(*), INTENT(IN)           :: FileName                       ! Name of the input file
   CHARACTER(200)                     :: LINE

   REAL(ReKi)                         :: L_Usr_Tmp
   REAL(ReKi)                         :: Sigma_USR_Tmp
   REAL(ReKi)                         :: U_USR_Tmp
   REAL(ReKi)                         :: WindDir_USR_Tmp
   REAL(ReKi)                         :: Z_USR_Tmp

   INTEGER                            :: I
   INTEGER                            :: Indx
   INTEGER                            :: J
   INTEGER                            :: IOAstat                        ! Input/Output/Allocate status
   INTEGER, INTENT(IN), OPTIONAL      :: NLines                         ! Number of lines to be skipped, if the file must be rewound
   INTEGER, INTENT(IN)                :: U_in                           ! Input unit.  This file is assumed to be already open

   LOGICAL                            :: ReadSigL                       ! Whether or not to read the last 2 columns

   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
      ! Find the end of the input file, where the "User-Defined Variables" are located
   READ ( U_in, '(A)', IOSTAT=IOAstat ) LINE

   IF ( IOAstat /= 0 ) THEN
      CALL TS_Abort( 'Could not read entire input file for user-defined variables.' )
   ENDIF

   CALL Conv2UC ( LINE )

   DO WHILE ( INDEX( LINE, 'USER-DEFINED' ) == 0 )

      READ ( U_in, '(A)', IOSTAT=IOAstat ) LINE

      IF ( IOAstat /= 0 ) THEN
         CALL TS_Abort( 'Could not read entire input file for user-defined variables.' )
      ENDIF

      CALL Conv2UC( LINE )

   ENDDO


      ! ---------- Read the size of the arrays --------------------------------------------
   CALL ReadVar( U_in, FileName, NumUSRz, "NumUSRz", "Number of heights in the user-defined profiles" )

   IF ( NumUSRz < 1 ) THEN
      CALL TS_Abort( 'The number of heights specified in the user-defined profiles must be at least 1.')
   ENDIF

   DO I=1,3
         ! ---------- Read the size of the arrays --------------------------------------------
      CALL ReadVar( U_in, FileName, StdScale(I), "StdScale", "Scaling value for user-defined standard deviation profile" )

      IF ( StdScale(I) <= 0. ) THEN
         CALL TS_Abort( 'The scaling value for the user-defined standard deviation profile must be positive.')
      ENDIF
   ENDDO

      ! Allocate the data arrays
   CALL AllocAry(Z_USR,       NumUSRz, 'Z_Usr (user-defined height)',               ErrStat2, ErrMsg2); CALL SetErrStat(ErrSTat2, ErrMsg2, ErrStat, ErrMsg, 'GetUSR')
   CALL AllocAry(U_USR,       NumUSRz, 'U_USR (user-defined wind speed)',           ErrStat2, ErrMsg2); CALL SetErrStat(ErrSTat2, ErrMsg2, ErrStat, ErrMsg, 'GetUSR')
   CALL AllocAry(WindDir_USR, NumUSRz, 'WindDir_USR (user-defined wind direction)', ErrStat2, ErrMsg2); CALL SetErrStat(ErrSTat2, ErrMsg2, ErrStat, ErrMsg, 'GetUSR')


   IF ( SpecModel == SpecModel_USRVKM ) THEN
      ReadSigL = .TRUE.

      CALL AllocAry(Sigma_USR, NumUSRz, 'Sigma_USR (user-defined sigma)',    ErrStat2, ErrMsg2); CALL SetErrStat(ErrSTat2, ErrMsg2, ErrStat, ErrMsg, 'GetUSR')
      CALL AllocAry(L_USR,     NumUSRz, 'L_USR (user-defined length scale)', ErrStat2, ErrMsg2); CALL SetErrStat(ErrSTat2, ErrMsg2, ErrStat, ErrMsg, 'GetUSR')
      
   ELSE
      ReadSigL = .FALSE.
   ENDIF

      ! ---------- Skip 4 lines --------------------------------------------
   DO I=1,4
      CALL ReadCom( U_in, FileName, "Headers for user-defined variables" )
   ENDDO

   DO I=1,NumUSRz

      IF ( ReadSigL ) THEN
         READ( U_in, *, IOSTAT=IOAstat ) Z_USR(I), U_USR(I), WindDir_USR(I), Sigma_USR(I), L_USR(I)
      ELSE
         READ( U_in, *, IOSTAT=IOAstat ) Z_USR(I), U_USR(I), WindDir_USR(I)
      ENDIF

      IF ( IOAstat /= 0 ) THEN
         CALL TS_Abort( 'Could not read entire user-defined variable list on line '//Int2LStr(I)//'.' )
      ENDIF

      IF ( ReadSigL ) THEN
         IF ( Sigma_USR(I) <= REAL( 0., ReKi ) ) THEN
            CALL TS_Abort( 'The standard deviation must be a positive number.' );
         ELSEIF ( L_USR(I) <= REAL( 0., ReKi ) ) THEN
            CALL TS_Abort( 'The length scale must be a positive number.' );
         ENDIF
      ENDIF

      IF ( WindDir_USR(I) > 360. ) THEN
         J = INT ( WindDir_USR(I) / 360. )
         WindDir_USR(I) = WindDir_USR(I) - J * 360.
      ELSEIF ( WindDir_USR(I) < 0. ) THEN
         J = INT ( -WindDir_USR(I) / 360. ) +1
         WindDir_USR(I) = WindDir_USR(I) + J * 360.
      ENDIF
   ENDDO

      ! Sort the arrays
   DO I=2,NumUSRz
      IF ( Z_USR(I) < Z_USR(I-1) ) THEN

         Indx = 1
         DO J=I-2,1,-1
            IF ( Z_USR(I) > Z_USR(J) ) THEN
               Indx = J+1
               EXIT
            ELSEIF ( Z_USR(I) == Z_USR(J) ) THEN
               CALL TS_Abort( 'Error: user-defined values must contain unique heights.' )
            ENDIF
         ENDDO

         Z_USR_Tmp       = Z_USR(I)
         U_USR_Tmp       = U_USR(I)
         WindDir_USR_Tmp = WindDir_USR(I)

         DO J=I,Indx+1,-1
            Z_USR(J)       = Z_USR(J-1)
            U_USR(J)       = U_USR(J-1)
            WindDir_USR(J) = WindDir_USR(J-1)
         ENDDO

         Z_USR(Indx)       = Z_USR_Tmp
         U_USR(Indx)       = U_USR_Tmp
         WindDir_USR(Indx) = WindDir_USR_Tmp

         IF ( ReadSigL ) THEN
            Sigma_USR_Tmp   = Sigma_USR(I)
            L_USR_Tmp       = L_USR(I)

            DO J=I,Indx+1,-1
               Sigma_USR(J)   = Sigma_USR(J-1)
               L_USR(J)       = L_USR(J-1)
            ENDDO

            Sigma_USR(Indx)   = Sigma_USR_Tmp
            L_USR(Indx)       = L_USR_Tmp
         ENDIF ! ReadSigL

      ENDIF
   ENDDO

      ! Rewind the file, if necessary.
   IF ( PRESENT(NLines) ) THEN
      REWIND( U_in , IOSTAT=IOAstat )

      IF ( IOAstat /= 0 ) THEN
         CALL TS_Abort( 'Error rewinding the file '//TRIM(FileName)//'.' )
      ENDIF

      DO I = 1,NLines
         CALL ReadCom( U_in, FileName, "Line "//Int2LStr(I) )
      ENDDO
   ENDIF

END SUBROUTINE GetUSR
!=======================================================================
SUBROUTINE GetUSRSpec(FileName,ErrStat,ErrMsg)

   USE                                   TSMods, ONLY: TurbModel, SpecModel         ! Type of turbulence model
   USE                                   TSMods, ONLY: NumUSRf           ! Length of user-defined spectra (# frequencies)
   USE                                   TSMods, ONLY: USpec             ! Unit number for the user-defined spectra file
   USE                                   TSMods, ONLY: Freq_USR          ! Frequencies for the user-defined spectra
   USE                                   TSMods, ONLY: Uspec_USR         ! User-defined u-component spectra
   USE                                   TSMods, ONLY: Vspec_USR         ! User-defined v-component spectra
   USE                                   TSMods, ONLY: Wspec_USR         ! User-defined w-component spectra

   IMPLICIT                              NONE

   INTEGER(IntKi),                  intent(  out) :: ErrStat                         ! Error level
   CHARACTER(*),                    intent(  out) :: ErrMsg                          ! Message describing error
   CHARACTER(*), INTENT(IN)           :: FileName                       ! Name of the input file
   CHARACTER(200)                     :: LINE

   REAL(ReKi)                         :: Freq_USR_Tmp
   REAL(ReKi)                         :: U_USR_Tmp
   REAL(ReKi)                         :: V_USR_Tmp
   REAL(ReKi)                         :: W_USR_Tmp
   REAL(ReKi)                         :: SpecScale (3)

   INTEGER                            :: I
   INTEGER                            :: Indx
   INTEGER                            :: J
   INTEGER                            :: IOAstat                        ! Input/Output/Allocate status

   INTEGER(IntKi)                                 :: ErrStat2                         ! Error level (local)
   CHARACTER(MaxMsgLen)                           :: ErrMsg2                          ! Message describing error (local)

   ErrStat = ErrID_None
   ErrMSg  = ""
   
      ! --------- Open the file ---------------

   CALL OpenFInpFile( USpec, FileName )

   CALL WrScr1(' Reading the user-defined spectra input file "'//TRIM(FileName)//'".' )


      ! --------- Read the comment lines at the beginning of the file ---------------
   DO I=1,3
      READ ( USpec, '(A)', IOSTAT=IOAstat ) LINE

      IF ( IOAstat /= 0 ) THEN
         CALL TS_Abort( 'Could not read entire input file for user-defined spectra.' )
      ENDIF
   ENDDO


      ! ---------- Read the size of the arrays --------------------------------------------
   CALL ReadVar( USpec, FileName, NumUSRf, "NumUSRf", "Number of frequencies in the user-defined spectra" )

   IF ( NumUSRf < 3 ) THEN
      CALL TS_Abort( 'The number of frequencies specified in the user-defined spectra must be at least 3.')
   ENDIF

   DO I=1,3
         ! ---------- Read the scaling for the arrays --------------------------------------------
      CALL ReadVar( USpec, FileName, SpecScale(I), "SpecScale", "Scaling value for user-defined standard deviation profile" )

      IF ( SpecScale(I) <= 0. ) THEN
         CALL TS_Abort( 'The scaling value for the user-defined spectra must be positive.')
      ENDIF
   ENDDO

      ! Allocate the data arrays
   CALL AllocAry( Freq_USR,  NumUSRf, 'Freq_USR (user-defined frequencies)', ErrStat2, ErrMsg2); CALL SetErrStat(ErrSTat2, ErrMsg2, ErrStat, ErrMsg, 'GetUSR')
   CALL AllocAry( Uspec_USR, NumUSRf, 'Uspec_USR (user-defined u spectra)' , ErrStat2, ErrMsg2); CALL SetErrStat(ErrSTat2, ErrMsg2, ErrStat, ErrMsg, 'GetUSR')
   CALL AllocAry( Vspec_USR, NumUSRf, 'Vspec_USR (user-defined v spectra)' , ErrStat2, ErrMsg2); CALL SetErrStat(ErrSTat2, ErrMsg2, ErrStat, ErrMsg, 'GetUSR')
   CALL AllocAry( Wspec_USR, NumUSRf, 'Wspec_USR (user-defined w spectra)' , ErrStat2, ErrMsg2); CALL SetErrStat(ErrSTat2, ErrMsg2, ErrStat, ErrMsg, 'GetUSR')



      ! ---------- Skip 4 lines --------------------------------------------
   DO I=1,4
      CALL ReadCom( USpec, FileName, "Headers for user-defined variables" )
   ENDDO

      ! ---------- Read the data lines --------------------------------------
   DO I=1,NumUSRf

      READ( USpec, *, IOSTAT=IOAstat ) Freq_USR(I), Uspec_USR(I), Vspec_USR(I), Wspec_USR(I)

      IF ( IOAstat /= 0 ) THEN
         CALL TS_Abort( 'Could not read entire user-defined spectra on line '//Int2LStr(I)//'.' )
      ENDIF

      IF ( ( Uspec_USR(I) <= REAL( 0., ReKi ) ) .OR. &
           ( Vspec_USR(I) <= REAL( 0., ReKi ) ) .OR. &
           ( Wspec_USR(I) <= REAL( 0., ReKi ) ) ) THEN
         CALL TS_Abort( 'The spectra must contain positive numbers.' );
!      ELSEIF ( Freq_USR(I) <= REAL( 0., ReKi ) ) THEN
!         CALL TS_Abort( 'The frequencies must be a positive number.' );
      ENDIF

         ! Scale by the factors earlier in the input file

      Uspec_USR(I) = Uspec_USR(I)*SpecScale(1)
      Vspec_USR(I) = Vspec_USR(I)*SpecScale(2)
      Wspec_USR(I) = Wspec_USR(I)*SpecScale(3)

   ENDDO

      ! ------- Sort the arrays by frequency -----------------------------------
   DO I=2,NumUSRf
      IF ( Freq_USR(I) < Freq_USR(I-1) ) THEN

         Indx = 1
         DO J=I-2,1,-1
            IF ( Freq_USR(I) > Freq_USR(J) ) THEN
               Indx = J+1
               EXIT
            ELSEIF ( Freq_USR(I) == Freq_USR(J) ) THEN
               CALL TS_Abort( 'Error: user-defined spectra must contain unique frequencies.' )
            ENDIF
         ENDDO

         Freq_USR_Tmp    = Freq_USR(I)
         U_USR_Tmp       = Uspec_USR(I)
         V_USR_Tmp       = Vspec_USR(I)
         W_USR_Tmp       = Wspec_USR(I)

         DO J=I,Indx+1,-1
            Freq_USR(J)    = Freq_USR(J-1)
            Uspec_USR(J)   = Uspec_USR(J-1)
            Vspec_USR(J)   = Vspec_USR(J-1)
            Wspec_USR(J)   = Wspec_USR(J-1)
         ENDDO

         Freq_USR(Indx)    = Freq_USR_Tmp
         Uspec_USR(I)      = U_USR_Tmp
         Vspec_USR(I)      = V_USR_Tmp
         Wspec_USR(I)      = W_USR_Tmp

      ENDIF
   ENDDO

      ! --------- Close the file ---------------------------------------

   CLOSE( USpec )

END SUBROUTINE GetUSRSpec

!=======================================================================
SUBROUTINE ReadCVarDefault ( UnIn, Fil, CharVar, VarName, VarDescr, US, Def, IGNORE )


      ! This routine reads a single character variable from the next line of the input file.
      ! The input is allowed to be "default"

      ! Argument declarations:


   INTEGER, INTENT(IN)          :: UnIn                                            ! I/O unit for input file.
   INTEGER, INTENT(IN)          :: US                                              ! I/O unit for echo/summary file.

   LOGICAL, INTENT(INOUT)       :: Def                                             ! - on input whether or not to use the default - on output, whether a default was used
   LOGICAL, INTENT(IN), OPTIONAL:: IGNORE                                          ! whether to ignore this input

   CHARACTER(250)               :: CharLine                                        ! Character string being read.
   CHARACTER(*), INTENT(INOUT)  :: CharVar                                         ! Character variable being read.
   CHARACTER( *), INTENT(IN)    :: Fil                                             ! Name of the input file.
   CHARACTER( *), INTENT(IN)    :: VarDescr                                        ! Text string describing the variable.
   CHARACTER( *), INTENT(IN)    :: VarName                                         ! Text string containing the variable name.

      ! Local declarations:

!  CHARACTER(38)                :: Frmt = "( 2X, ES11.4e2, 2X, A, T30, ' - ', A )" ! Output format for real parameters
   CHARACTER(38)                :: Frmt = "( A10, 2X, A )"                         ! Output format for real parameters


   CALL ReadVar( UnIn, Fil, CharLine, VarName, VarDescr )

   IF ( PRESENT(IGNORE) ) THEN
      IF ( IGNORE ) THEN
         WRITE (US,Frmt)  "N/A", VarDescr
         Def = .TRUE.
         RETURN
      ENDIF

   ENDIF

   CALL Conv2UC( CharLine )

   IF ( TRIM(CharLine) == 'DEFAULT' ) THEN

      CALL WrScr ( '    A default value will be used for '//TRIM(VarName)//'.' )

      IF ( Def ) THEN  ! use the value as a default
         WRITE (US,Frmt)  CharVar, VarDescr
      ELSE
         WRITE (US,Frmt)  "CALCULATED", VarDescr
      ENDIF

      Def = .TRUE.

   ELSE

      CharVar = CharLine

      WRITE (US,Frmt)  CharVar, VarDescr

      Def = .FALSE.

   ENDIF

   RETURN
END SUBROUTINE ReadCVarDefault ! ( UnIn, Fil, RealVar, VarName, VarDescr )
!=======================================================================
SUBROUTINE ReadRAryDefault ( UnIn, Fil, RealAry, VarName, VarDescr, US, Def, IGNORE )

      ! This routine reads a real array from the next line of the input file.
      ! The input is allowed to be "default"

      ! Argument declarations:

   REAL(ReKi), INTENT(INOUT)    :: RealAry (:)                                     ! Real variable being read.

   INTEGER, INTENT(IN)          :: UnIn                                            ! I/O unit for input file.
   INTEGER, INTENT(IN)          :: US                                              ! I/O unit for echo/summary file.

   LOGICAL, INTENT(INOUT)       :: Def                                             ! - on input whether or not to use the default - on output, whether a default was used
   LOGICAL, INTENT(IN), OPTIONAL:: IGNORE                                          ! whether or not to ignore this input

   CHARACTER(250)               :: CharLine                                        ! Character string being read.
   CHARACTER( *), INTENT(IN)    :: Fil                                             ! Name of the input file.
   CHARACTER( *), INTENT(IN)    :: VarDescr                                        ! Text string describing the variable.
   CHARACTER( *), INTENT(IN)    :: VarName                                         ! Text string containing the variable name.

      ! Local declarations:

   INTEGER                      :: IOS                                             ! I/O status returned from the read statement.

!  CHARACTER(38)                :: Frmt = "( 2X, ES11.4e2, 2X, A, T30, ' - ', A )" ! Output format for real parameters
   CHARACTER(38)                :: Frmt = "( F10.3, 2X, , A )"                     ! Output format for real parameters

   WRITE(Frmt,"(I2)") SIZE(RealAry)-1

   Frmt = "( '(',F9.3,"//TRIM(Frmt)//"(',',G10.3),')',2X , A )"

   CALL ReadVar( UnIn, Fil, CharLine, VarName, VarDescr ) !Maybe I should read this in explicitly...

   IF ( PRESENT(IGNORE) ) THEN
      IF ( IGNORE ) THEN
         WRITE (US,"( A10, 2X, A )")  "N/A", VarDescr
         Def = .TRUE.
         RETURN
      ENDIF

   ENDIF

   CALL Conv2UC( CharLine )

   IF ( TRIM(CharLine) == 'DEFAULT' ) THEN

      CALL WrScr ( '    A default value will be used for '//TRIM(VarName)//'.' )

      IF ( Def ) THEN  ! use the value as a default
         WRITE (US,Frmt)  RealAry, VarDescr
      ELSE
         WRITE (US,"( A10, 2X, A )")  "CALCULATED", VarDescr
      ENDIF

      Def = .TRUE.

   ELSE

      IF ( INDEX( CharLine(1:1), 'TF') > 0 ) THEN    ! We don't want 'T' or 'F' read as -1 or 0.
         CALL WrScr1 ( ' Invalid numerical input for "'//TRIM( VarName )//'".' )
      ENDIF

      READ (CharLine,*,IOSTAT=IOS)  RealAry

      IF (IOS /=0) THEN
         READ (CharLine,*,IOSTAT=IOS)  RealAry(1)  ! Try reading only the first element
      ENDIF

      CALL CheckIOS ( IOS, Fil, VarName, NumType )

      WRITE (US,Frmt)  RealAry, VarDescr

      Def = .FALSE.

   ENDIF


   RETURN

END SUBROUTINE ReadRAryDefault
!=======================================================================
SUBROUTINE ReadRVarDefault ( UnIn, Fil, RealVar, VarName, VarDescr, US, Def, IGNORE, IGNORESTR )

      ! This routine reads a single real variable from the next line of the input file.
      ! The input is allowed to be "default"

      ! Argument declarations:

   REAL(ReKi), INTENT(INOUT)    :: RealVar                                         ! Real variable being read.

   INTEGER, INTENT(IN)          :: UnIn                                            ! I/O unit for input file.
   INTEGER, INTENT(IN)          :: US                                              ! I/O unit for echo/summary file.

   LOGICAL, INTENT(INOUT)       :: Def                                             ! - on input whether or not to use the default - on output, whether a default was used
   LOGICAL, INTENT(IN), OPTIONAL:: IGNORE                                          ! whether or not to ignore this input
   LOGICAL, INTENT(OUT),OPTIONAL:: IGNORESTR                                       ! whether or not user requested to ignore this input

   CHARACTER(250)               :: CharLine                                        ! Character string being read.
   CHARACTER( *), INTENT(IN)    :: Fil                                             ! Name of the input file.
   CHARACTER( *), INTENT(IN)    :: VarDescr                                        ! Text string describing the variable.
   CHARACTER( *), INTENT(IN)    :: VarName                                         ! Text string containing the variable name.

      ! Local declarations:

   INTEGER                      :: IOS                                             ! I/O status returned from the read statement.

!  CHARACTER(38)                :: Frmt = "( 2X, ES11.4e2, 2X, A, T30, ' - ', A )" ! Output format for real parameters
   CHARACTER(38)                :: Frmt = "( F10.3, 2X, A )"                       ! Output format for real parameters


   CALL ReadVar( UnIn, Fil, CharLine, VarName, VarDescr )

   IF ( PRESENT(IGNORE) ) THEN

      IF ( IGNORE ) THEN
         WRITE (US,"( A10, 2X, A )")  "N/A", VarDescr
         Def = .TRUE.
         RETURN
      ENDIF

   ENDIF


   CALL Conv2UC( CharLine )

   IF ( PRESENT(IGNORESTR) ) THEN
      IF ( TRIM( CharLine ) == 'NONE' ) THEN
         IGNORESTR = .TRUE.
         WRITE (US,"( A10, 2X, A )")  "N/A", VarDescr
         Def = .TRUE.
         RealVar = 0.0  ! This is set for the Reynolds stress inputs, but if IGNORESTR is used for other inputs, it may need to be changed
         RETURN
      ENDIF
   ENDIF

   IF ( TRIM(CharLine) == 'DEFAULT' ) THEN

      CALL WrScr ( '    A default value will be used for '//TRIM(VarName)//'.' )

      IF ( PRESENT(IGNORESTR) ) THEN
         IF ( IGNORESTR ) THEN  !We've told it to ignore this input, by default
            WRITE (US,"( A10, 2X, A )")  "N/A", VarDescr
            RealVar = 0.0  ! This is set for the Reynolds stress inputs, but if IGNORESTR is used for other inputs, it may need to be changed
            RETURN
         ENDIF
      ENDIF

      IF ( Def ) THEN  ! use the value as a default
         WRITE (US,Frmt)  RealVar, VarDescr
      ELSE
         WRITE (US,"( A10, 2X, A )")  "CALCULATED", VarDescr
      ENDIF

      Def = .TRUE.

   ELSE

      IF ( INDEX( CharLine(1:1), 'TF') > 0 ) THEN    ! We don't want 'T' or 'F' read as -1 or 0.
         CALL WrScr1 ( ' Invalid numerical input for "'//TRIM( VarName )//'".' )
      ENDIF

      READ (CharLine,*,IOSTAT=IOS)  RealVar

      CALL CheckIOS ( IOS, Fil, VarName, NumType )

      WRITE (US,Frmt)  RealVar, VarDescr

      Def = .FALSE.

      IF ( PRESENT(IGNORESTR) ) THEN
         IGNORESTR = .FALSE.
      ENDIF

   ENDIF


   RETURN
END SUBROUTINE ReadRVarDefault
!=======================================================================

SUBROUTINE WrBinBLADED(USig, VSig, WSig)

USE                         TSMods


IMPLICIT                    NONE

REAL(ReKi)                  :: U_C1                          ! Scale for converting BLADED U data
REAL(ReKi)                  :: U_C2                          ! Offset for converting BLADED U data
REAL(ReKi),INTENT(INOUT)    :: USig                          ! Standard deviation of U
REAL(ReKi)                  :: V_C                           ! Scale for converting BLADED V data
REAL(ReKi),INTENT(INOUT)    :: VSig                          ! Standard deviation of V
REAL(ReKi)                  :: W_C                           ! Scale for converting BLADED W data
REAL(ReKi),INTENT(INOUT)    :: WSig                          ! Standard deviation of W
REAL(ReKi)                  :: TI(3)                         ! Turbulence intensity for scaling data
REAL(ReKi)                  :: ZTmp                          ! This is the vertical center of the grid

INTEGER(B4Ki)               :: CFirst
INTEGER(B4Ki)               :: CLast
INTEGER(B4Ki)               :: CStep
INTEGER(B4Ki)               :: I
INTEGER(B4Ki)               :: II
INTEGER(B4Ki)               :: IT
INTEGER(B4Ki)               :: IY
INTEGER(B4Ki)               :: IZ

INTEGER(B4Ki)               :: IP
INTEGER(B2Ki)               :: TmpVarray(3*p_grid%NumGrid_Y*p_grid%NumGrid_Z) ! This array holds the normalized velocities before being written to the binary file
INTEGER(B2Ki),ALLOCATABLE   :: TmpTWRarray(:)                   ! This array holds the normalized tower velocities

INTEGER                     :: AllocStat



      ! Put normalizing factors into the summary file.  The user can use them to
      ! tell a simulation program how to rescale the data.

   USig  = MAX(100.0*Tolerance, USig)
   VSig  = MAX(100.0*Tolerance, VSig)
   WSig  = MAX(100.0*Tolerance, WSig)

   TI(1) = USig / UHub
   TI(2) = VSig / UHub
   TI(3) = WSig / UHub

   FormStr = "(//,'Normalizing Parameters for Binary Data (approximate statistics):',/)"
   WRITE (US,FormStr)

   FormStr = "(3X,A,' =',F9.4,A)"
   WRITE (US,FormStr)  'UBar ', UHub, ' m/s'
   WRITE (US,FormStr)  'TI(u)', 100.0*TI(1), ' %'
   WRITE (US,FormStr)  'TI(v)', 100.0*TI(2), ' %'
   WRITE (US,FormStr)  'TI(w)', 100.0*TI(3), ' %'

   Ztmp  = ( p_grid%HubHt - p_grid%GridHeight / 2.0 - p_grid%Z(1) )  ! This is the grid offset

   WRITE (US,'()')
   WRITE (US,FormStr)  'Height Offset', Ztmp,        ' m'
   WRITE (US,FormStr)  'Grid Base    ', p_grid%Z(1), ' m'
   
   WRITE (US,'()'   )
   WRITE (US,'( A)' ) 'Creating a PERIODIC output file.'
 
      ! Calculate some numbers for normalizing the data.

   U_C1 = 1000.0/( UHub*TI(1) )
   U_C2 = 1000.0/TI(1)
   V_C  = 1000.0/( UHub*TI(2) )
   W_C  = 1000.0/( UHub*TI(3) )


   ZTmp     = p_grid%Z(1) + p_grid%GridHeight/2.0  !This is the vertical center of the grid

   IF ( WrBLFF )  THEN

      CALL WrScr ( ' Generating BLADED binary time-series file "'//TRIM( RootName )//'.wnd"' )

               ! Put header information into the binary data file.

      WRITE (UBFFW)   INT(  -99          , B2Ki )               ! -99 = New Bladed format
      WRITE (UBFFW)   INT(    4          , B2Ki )               ! 4 = improved von karman (but needed for next 7 inputs)
      WRITE (UBFFW)   INT(    3          , B4Ki )               ! 3 = number of wind components
      WRITE (UBFFW)  REAL( Latitude      , SiKi )               ! Latitude (degrees)
      WRITE (UBFFW)  REAL(   z0          , SiKi )               ! Roughness length (m)
      WRITE (UBFFW)  REAL( Ztmp          , SiKi )               ! Reference Height (m) ( Z(1) + GridHeight / 2.0 )
      WRITE (UBFFW)  REAL( 100.0*TI(1)   , SiKi )               ! Longitudinal turbulence intensity (%)
      WRITE (UBFFW)  REAL( 100.0*TI(2)   , SiKi )               ! Lateral turbulence intensity (%)
      WRITE (UBFFW)  REAL( 100.0*TI(3)   , SiKi )               ! Vertical turbulence intensity (%)

      WRITE (UBFFW)  REAL( p_grid%GridRes_Z     , SiKi )               ! grid spacing in vertical direction, in m
      WRITE (UBFFW)  REAL( p_grid%GridRes_Y     , SiKi )               ! grid spacing in lateral direction, in m
      WRITE (UBFFW)  REAL( p_grid%TimeStep*UHub , SiKi )               ! grid spacing in longitudinal direciton, in m
      WRITE (UBFFW)   INT( p_grid%NumOutSteps/2 , B4Ki )               ! half the number of points in alongwind direction
      WRITE (UBFFW)  REAL( UHub          , SiKi )               ! the mean wind speed in m/s
      WRITE (UBFFW)  REAL( 0             , SiKi )               ! the vertical length scale of the longitudinal component in m
      WRITE (UBFFW)  REAL( 0             , SiKi )               ! the lateral length scale of the longitudinal component in m
      WRITE (UBFFW)  REAL( 0             , SiKi )               ! the longitudinal length scale of the longitudinal component in m
      WRITE (UBFFW)   INT( 0             , B4Ki )               ! an unused integer
      WRITE (UBFFW)   INT( p_RandNum%RandSeed(1)   , B4Ki )               ! the random number seed
      WRITE (UBFFW)   INT( p_grid%NumGrid_Z     , B4Ki )               ! the number of grid points vertically
      WRITE (UBFFW)   INT( p_grid%NumGrid_Y     , B4Ki )               ! the number of grid points laterally
      WRITE (UBFFW)   INT( 0             , B4Ki )               ! the vertical length scale of the lateral component, not used
      WRITE (UBFFW)   INT( 0             , B4Ki )               ! the lateral length scale of the lateral component, not used
      WRITE (UBFFW)   INT( 0             , B4Ki )               ! the longitudinal length scale of the lateral component, not used
      WRITE (UBFFW)   INT( 0             , B4Ki )               ! the vertical length scale of the vertical component, not used
      WRITE (UBFFW)   INT( 0             , B4Ki )               ! the lateral length scale of the vertical component, not used
      WRITE (UBFFW)   INT( 0             , B4Ki )               ! the longitudinal length scale of the vertical component, not used


         ! Compute parameters for ordering output for FF AeroDyn files. (This is for BLADED compatibility.)

      IF ( Clockwise )  THEN
         CFirst = p_grid%NumGrid_Y
         CLast  = 1
         CStep  = -1
      ELSE
         CFirst = 1
         CLast  = p_grid%NumGrid_Y
         CStep  = 1
      ENDIF


         ! Loop through time.

      DO IT=1,p_grid%NumOutSteps  !Use only the number of timesteps requested originally

            ! Write out grid data in binary form.
         IP = 1
         DO IZ=1,p_grid%NumGrid_Z
            DO IY=CFirst,CLast,CStep

               II = ( IZ - 1 )*p_grid%NumGrid_Y + IY

               TmpVarray(IP)   = NINT( U_C1*V(IT,II,1) - U_C2, B2Ki )  ! Put the data into a temp array so that the WRITE() command works faster
               TmpVarray(IP+1) = NINT( V_C *V(IT,II,2)       , B2Ki )
               TmpVarray(IP+2) = NINT( W_C *V(IT,II,3)       , B2Ki )

               IP = IP + 3;
            ENDDO ! IY
         ENDDO ! IZ

         WRITE ( UBFFW )  TmpVarray(:) ! bjj: We cannot write the array including time because of stack overflow errors.. otherwise use compile option to put this on the heap instead of the stack

      ENDDO ! IT

      CLOSE ( UBFFW )


   ENDIF ! WrBLFF


   IF ( WrADTWR .AND. ( WrBLFF .OR. .NOT. WrADFF ) ) THEN
      IF ( p_grid%ExtraHubPT ) THEN
         IZ = p_grid%ZLim - p_grid%NumGrid_Z - 1
         I  = p_grid%NumGrid_Z*p_grid%NumGrid_Y + 2
      ELSE
         IZ = p_grid%ZLim - p_grid%NumGrid_Z
         I  = p_grid%NumGrid_Z*p_grid%NumGrid_Y + 1
      ENDIF

      IF ( p_grid%ExtraTwrPT ) THEN
         IY = I
         I  = I + 1
      ELSE
         IZ = IZ + 1
         IY = (p_grid%NumGrid_Y / 2) + 1      !The grid location of the top tower point
      ENDIF


      CALL WrScr ( ' Generating tower binary time-series file "'//TRIM( RootName )//'.twr"' )


      WRITE (UATWR)  REAL( p_grid%GridRes_Z ,SiKi )         ! grid spacing in vertical direction, in m
      WRITE (UATWR)  REAL( p_grid%TimeStep*UHub ,   SiKi )         ! grid spacing in longitudinal direciton, in m
      WRITE (UATWR)  REAL( p_grid%Z(1) ,     SiKi )         ! The vertical location of the highest tower grid point in m
      WRITE (UATWR)   INT( p_grid%NumOutSteps ,     B4Ki )         ! The number of points in alongwind direction
      WRITE (UATWR)   INT( IZ ,              B4Ki )         ! the number of grid points vertically
      WRITE (UATWR)  REAL( UHub ,            SiKi )         ! the mean wind speed in m/s
      WRITE (UATWR)  REAL( 100.0*TI(1),      SiKi )         ! Longitudinal turbulence intensity
      WRITE (UATWR)  REAL( 100.0*TI(2),      SiKi )         ! Lateral turbulence intensity
      WRITE (UATWR)  REAL( 100.0*TI(3),      SiKi )         ! Vertical turbulence intensity


      ALLOCATE ( TmpTWRarray( 3*(p_grid%NPoints-I+2) ) , STAT=AllocStat )

      IF ( AllocStat /= 0 )  THEN
         CALL TS_Abort ( 'Error allocating memory for temporary tower wind speed array.' )
      ENDIF


      DO IT=1,p_grid%NumOutSteps

         TmpTWRarray(1) = NINT( U_C1*V(IT,IY,1) - U_C2 , B2Ki )
         TmpTWRarray(2) = NINT( V_C *V(IT,IY,2)        , B2Ki )
         TmpTWRarray(3) = NINT( W_C *V(IT,IY,3)        , B2Ki )

         IP = 4
         DO II = I,p_grid%NPoints    ! This assumes we have points in a single line along the center of the hub
            TmpTWRarray(IP  ) = NINT( U_C1*V(IT,II,1) - U_C2 , B2Ki )
            TmpTWRarray(IP+1) = NINT( V_C *V(IT,II,2)        , B2Ki )
            TmpTWRarray(IP+2) = NINT( W_C *V(IT,II,3)        , B2Ki )

            IP = IP + 3
         ENDDO    ! II

         WRITE (UATWR) TmpTWRarray(:)

      ENDDO ! IT

      CLOSE ( UATWR )

      IF ( ALLOCATED( TmpTWRarray) ) DEALLOCATE( TmpTWRarray )

   ENDIF !WrADWTR


END SUBROUTINE WrBinBLADED
!=======================================================================
SUBROUTINE WrBinTURBSIM

USE                       TSMods


IMPLICIT                  NONE

REAL(SiKi), PARAMETER     :: IntMax   =  32767.0
REAL(SiKi), PARAMETER     :: IntMin   = -32768.0
REAL(SiKi), PARAMETER     :: IntRng   = IntMax - IntMin ! Max Range of 2-byte integer

REAL(SiKi)                :: UOff                       ! Offset for the U component
REAL(SiKi)                :: UScl                       ! Slope  for the U component
REAL(ReKi)                :: VMax(3)                    ! Maximum value of the 3 wind components
REAL(ReKi)                :: VMin(3)                    ! Minimum value of the 3 wind components
REAL(SiKi)                :: VOff                       ! Offset for the V component
REAL(SiKi)                :: VScl                       ! Slope  for the V component
REAL(SiKi)                :: WOff                       ! Offset for the W component
REAL(SiKi)                :: WScl                       ! Slope  for the W component

INTEGER, PARAMETER        :: DecRound  = 3              ! Number of decimal places to round to
INTEGER(B2Ki)             :: FileID                     ! File ID, determines specific output format (if periodic or not)
INTEGER                   :: IC                         ! counter for the velocity component of V
INTEGER                   :: II                         ! counter for the point on the grid/tower
INTEGER                   :: IT                         ! counter for the timestep
INTEGER(B4Ki)             :: LenDesc                    ! Length of the description string
INTEGER(B4Ki)             :: NumGrid                    ! Number of points on the grid
INTEGER(B4Ki)             :: NumTower                   ! Number of points on the tower
INTEGER(B4Ki)             :: TwrStart                   ! First index of a tower point
INTEGER(B4Ki)             :: TwrTop                     ! The index of top of the tower (it could be on the grid instead of at the end)

INTEGER(B4Ki)             :: IP
INTEGER(B2Ki),ALLOCATABLE :: TmpVarray(:)                ! This array holds the normalized velocities before being written to the binary file

INTEGER                   :: AllocStat



      ! Set the file format ID
      
   IF ( Periodic ) THEN 
      FileID = 8
   ELSE
      FileID = 7
   END IF      


      ! Find the range of our velocity

   DO IC=1,3

         ! Initialize the Min/Max values

      VMin(IC) = V(1,1,IC)
      VMax(IC) = V(1,1,IC)

      DO II=1,p_grid%NPoints   ! Let's check all of the points
         DO IT=1,p_grid%NumOutSteps  ! Use only the number of timesteps requested originally

            IF ( V(IT,II,IC) > VMax(IC) ) THEN

               VMax(IC) = V(IT,II,IC)

            ELSEIF ( V(IT,II,IC) < VMin(IC) ) THEN

               VMin(IC) = V(IT,II,IC)

            ENDIF

         ENDDO !IT
      ENDDO !II

   ENDDO !IC


      ! Calculate the scaling parameters for each component


   IF ( VMax(1) == VMin(1) ) THEN
      UScl = 1
   ELSE
      UScl = IntRng/REAL( VMax(1) - VMin(1) , SiKi )
   ENDIF

   IF ( VMax(2) == VMin(2) ) THEN
      VScl = 1
   ELSE
      VScl = IntRng/REAL( VMax(2) - VMin(2) , SiKi )
   ENDIF

   IF ( VMax(3) == VMin(3) ) THEN
      WScl = 1
   ELSE
      WScl = IntRng/REAL( VMax(3) - VMin(3) , SiKi )
   ENDIF


   UOff = IntMin - UScl*REAL( VMin(1)    , SiKi )
   VOff = IntMin - VScl*REAL( VMin(2)    , SiKi )
   WOff = IntMin - WScl*REAL( VMin(3)    , SiKi )


      ! Find the first tower point

   NumGrid  = p_grid%NumGrid_Y*p_grid%NumGrid_Z

   IF ( WrADTWR ) THEN

      TwrStart = NumGrid + 1

      IF ( p_grid%ExtraHubPT ) THEN
         TwrStart = TwrStart + 1
      ENDIF

      IF ( p_grid%ExtraTwrPT ) THEN
         TwrTop   = TwrStart
         TwrStart = TwrStart + 1
      ELSE
         TwrTop = INT(p_grid%NumGrid_Y / 2) + 1      ! The top tower point is on the grid where Z = 1
      ENDIF

      NumTower = p_grid%NPoints - TwrStart + 2

   ELSE

      NumTower = 0

   ENDIF


   LenDesc = LEN_TRIM( DescStr )             ! Length of the string that contains program name, version, date, and time

   CALL WrScr ( ' Generating AeroDyn binary time-series file "'//TRIM( RootName )//'.bts"' )



      ! Write the header

   WRITE (UAFFW)   INT( FileID             , B2Ki )          ! TurbSim format (7=not PERIODIC, 8=PERIODIC)

   WRITE (UAFFW)   INT( p_grid%NumGrid_Z          , B4Ki )          ! the number of grid points vertically
   WRITE (UAFFW)   INT( p_grid%NumGrid_Y          , B4Ki )          ! the number of grid points laterally
   WRITE (UAFFW)   INT( NumTower           , B4Ki )          ! the number of tower points
   WRITE (UAFFW)   INT( p_grid%NumOutSteps        , B4Ki )          ! the number of time steps

   WRITE (UAFFW)  REAL( p_grid%GridRes_Z          , SiKi )          ! grid spacing in vertical direction, in m
   WRITE (UAFFW)  REAL( p_grid%GridRes_Y          , SiKi )          ! grid spacing in lateral direction, in m
   WRITE (UAFFW)  REAL( p_grid%TimeStep           , SiKi )          ! grid spacing in delta time, in m/s
   WRITE (UAFFW)  REAL( UHub               , SiKi )          ! the mean wind speed in m/s at hub height
   WRITE (UAFFW)  REAL( p_grid%HubHt              , SiKi )          ! the hub height, in m
   WRITE (UAFFW)  REAL( p_grid%Z(1)               , SiKi )          ! the height of the grid bottom, in m

   WRITE (UAFFW)  REAL( UScl               , SiKi )          ! the U-component slope for scaling
   WRITE (UAFFW)  REAL( UOff               , SiKi )          ! the U-component offset for scaling
   WRITE (UAFFW)  REAL( VScl               , SiKi )          ! the V-component slope for scaling
   WRITE (UAFFW)  REAL( VOff               , SiKi )          ! the V-component offset for scaling
   WRITE (UAFFW)  REAL( WScl               , SiKi )          ! the W-component slope for scaling
   WRITE (UAFFW)  REAL( WOff               , SiKi )          ! the W-component offset for scaling

   WRITE (UAFFW)   INT( LenDesc            , B4Ki )          ! the number of characters in the string, max 200

   DO II=1,LenDesc

      WRITE (UAFFW)  INT( IACHAR( DescStr(II:II) ), B1Ki )   ! converted ASCII characters

   ENDDO

      ALLOCATE ( TmpVarray( 3*(p_grid%NumGrid_Z*p_grid%NumGrid_Y + NumTower) ) , STAT=AllocStat )

      IF ( AllocStat /= 0 )  THEN
         CALL TS_Abort ( 'Error allocating memory for temporary wind speed array.' )
      ENDIF

      ! Loop through time.

   DO IT=1,p_grid%NumOutSteps  !Use only the number of timesteps requested originally

         ! Write out grid data in binary form. II = (IZ - 1)*NumGrid_Y + IY, IY varies most rapidly

      IP = 1

      DO II=1,NumGrid

         TmpVarray(IP)   =  NINT( Max( Min( REAL(UScl*V(IT,II,1) + UOff, SiKi), IntMax ),IntMin) , B2Ki )
         TmpVarray(IP+1) =  NINT( Max( Min( REAL(VScl*V(IT,II,2) + VOff, SiKi), IntMax ),IntMin) , B2Ki )
         TmpVarray(IP+2) =  NINT( Max( Min( REAL(WScl*V(IT,II,3) + Woff, SiKi), IntMax ),IntMin) , B2Ki )

         IP = IP + 3
      ENDDO ! II


      IF ( WrADTWR ) THEN

            ! Write out the tower data in binary form

            ! Value at the top of the tower (bottom of grid)
         TmpVarray(IP)   =  NINT( Max( Min( REAL(UScl*V(IT,TwrTop,1) + UOff, SiKi), IntMax ),IntMin) , B2Ki )
         TmpVarray(IP+1) =  NINT( Max( Min( REAL(VScl*V(IT,TwrTop,2) + VOff, SiKi), IntMax ),IntMin) , B2Ki )
         TmpVarray(IP+2) =  NINT( Max( Min( REAL(WScl*V(IT,TwrTop,3) + Woff, SiKi), IntMax ),IntMin) , B2Ki )

         IP = IP + 3
         DO II=TwrStart,p_grid%NPoints
                ! Values of tower data
            TmpVarray(IP)   =  NINT( Max( Min( REAL(UScl*V(IT,II,1) + UOff, SiKi), IntMax ),IntMin) , B2Ki )
            TmpVarray(IP+1) =  NINT( Max( Min( REAL(VScl*V(IT,II,2) + VOff, SiKi), IntMax ),IntMin) , B2Ki )
            TmpVarray(IP+2) =  NINT( Max( Min( REAL(WScl*V(IT,II,3) + Woff, SiKi), IntMax ),IntMin) , B2Ki )

            IP = IP + 3
         ENDDO ! II

      ENDIF

      WRITE ( UAFFW ) TmpVarray(:)
   ENDDO ! IT

   CLOSE ( UAFFW )

   IF ( ALLOCATED( TmpVarray ) ) DEALLOCATE( TmpVarray )


END SUBROUTINE WrBinTURBSIM
!=======================================================================
SUBROUTINE WrFormattedFF(HubIndx)

USE                     TSMods

IMPLICIT                NONE

REAL(ReKi), ALLOCATABLE      :: ZRow       (:)                           ! The horizontal locations of the grid points (NumGrid_Y) at each height.

INTEGER                      :: AllocStat
INTEGER, INTENT(IN)          :: HubIndx
INTEGER                      :: II
INTEGER                      :: IT
INTEGER                      :: IVec
INTEGER                      :: IY
INTEGER                      :: IZ

CHARACTER(200)               :: FormStr3                                 ! String used to store format specifiers.
CHARACTER(200)               :: FormStr4                                 ! String used to store format specifiers.
CHARACTER(200)               :: FormStr5                                 ! String used to store format specifiers.
CHARACTER(200)               :: FormStr6                                 ! String used to store format specifiers.



   FormStr  = "( / 'This full-field turbulence file was generated by ' , A , A , ' on ' , A , ' at ' , A , '.' / )"
   FormStr1 = "( ' | ', A,'-comp |  Y  x  Z  | Grid Resolution (Y x Z) | Time-step | Hub Elev | Mean U |')"
   FormStr2 = "(I14,I6,F11.3,F11.3,F15.3,F11.2,F10.2)"
   FormStr3 = "(/,' Z Coordinates (m):')"
   FormStr4 = "(/,' Y Coordinates (m):')"
   FormStr5 = "(1X,98(F8.3),:)"
   FormStr6 = "(/,1X,2(F8.3))"

      ! Allocate the array of wind speeds.

   ALLOCATE ( ZRow(p_grid%NumGrid_Y) , STAT=AllocStat )

   IF ( AllocStat /= 0 )  THEN
      CALL TS_Abort ( 'Error allocating memory for array of wind speeds.' )
   ENDIF


   DO IVec=1,3

      CALL WrScr ( ' Generating full-field formatted file "'//TRIM(RootName)//'.'//Comp(IVec)//'".' )

      CALL OpenFOutFile ( UFFF, TRIM( RootName )//'.'//Comp(IVec) )


         ! Create file header.

      WRITE (UFFF,FormStr )  TRIM(ProgName), TRIM( ProgVer ), CurDate(), CurTime()

      WRITE (UFFF,FormStr1)  Comp(IVec)

      WRITE (UFFF,FormStr2)  p_grid%NumGrid_Y, p_grid%NumGrid_Z, p_grid%GridRes_Y, p_grid%GridRes_Z, p_grid%TimeStep, p_grid%HubHt, UHub
      WRITE (UFFF,FormStr3)
      WRITE (UFFF,FormStr5)  ( p_grid%Z(IZ)-p_grid%HubHt, IZ=1,p_grid%NumGrid_Z )
      WRITE (UFFF,FormStr4)
      WRITE (UFFF,FormStr5)  ( p_grid%Y(IY), IY=1,p_grid%NumGrid_Y )

         ! Write out elapsed time & hub-level value before component grid.

      DO IT=1,p_grid%NumOutSteps

         WRITE(UFFF,FormStr6)  p_grid%TimeStep*( IT - 1 ), V(IT,HubIndx,IVec)

         DO IZ=1,p_grid%NumGrid_Z  ! From the top to the bottom

            II = ( p_grid%NumGrid_Z - IZ )*p_grid%NumGrid_Y

            DO IY=1,p_grid%NumGrid_Y  ! From the left to the right
               ZRow(IY) = V(IT,II+IY,IVec)
            ENDDO ! IY

            WRITE (UFFF,FormStr5)  ( ZRow(IY), IY=1,p_grid%NumGrid_Y )

         ENDDO ! IZ

      ENDDO ! IT

      CLOSE ( UFFF )

   ENDDO ! IVec

      ! Deallocate the array of wind speeds.

   IF ( ALLOCATED( ZRow ) )  DEALLOCATE( ZRow )

END SUBROUTINE WrFormattedFF
!=======================================================================

SUBROUTINE WrSum_Heading1( US  )


USE TSMods, only: NumUSRz
USE TSMods, only: Z_USR
USE TSMods, only: U_USR
USE TSMods, only: WindDir_USR
USE TSMods, only: Sigma_USR
USE TSMods, only: StdScale
USE TSMods, only: L_USR

   INTEGER, INTENT(IN) :: US
   integer             :: i 


   IF ( NumUSRz > 0 ) THEN
      WRITE (US,"( // 'User-Defined profiles:' / )")
   
      IF ( ALLOCATED( L_USR ) ) THEN
         WRITE (US,"(A)") '  Height   Wind Speed   Horizontal Angle   u Std. Dev.   v Std. Dev.   w Std. Dev.   Length Scale'
         WRITE (US,"(A)") '   (m)        (m/s)          (deg)            (m/s)         (m/s)         (m/s)          (m)     '
         WRITE (US,"(A)") '  ------   ----------   ----------------   -----------   -----------   -----------   ------------'
      
         DO I=NumUSRz,1,-1
            WRITE (US,"( 1X,F7.2, 2X,F9.2,2X, 3X,F10.2,6X, 3(4X,F7.2,3X), 3X,F10.2 )")  &
                                Z_USR(I), U_USR(I), WindDir_USR(I), &
                                Sigma_USR(I)*StdScale(1), Sigma_USR(I)*StdScale(2), Sigma_USR(I)*StdScale(3), L_USR(I)      
         ENDDO   
      ELSE
         WRITE (US,"(A)") '  Height   Wind Speed   Horizontal Angle'
         WRITE (US,"(A)") '   (m)        (m/s)          (deg)      '
         WRITE (US,"(A)") '  ------   ----------   ----------------'
     
         DO I=NumUSRz,1,-1
            WRITE (US,"( 1X,F7.2, 2X,F9.2,2X, 3X,F10.2,6X)") Z_USR(I), U_USR(I), WindDir_USR(I)      
         ENDDO   
      ENDIF
   
   ENDIF

END SUBROUTINE WrSum_Heading1
!=======================================================================
SUBROUTINE WrSum_SpecModel(US, p_IEC, GridVertCenter )

!USE TSMods, only: FormStr
!USE TSMods, only: FormStr1
!USE TSMods, only: FormStr2
!USE TSMods, only: TMName
!
!USE TSMods, only: IEC_NTM
!USE TSMods, only: IECeditionSTR
use TSMods


   INTEGER,    INTENT(IN) :: US
   TYPE(IEC_ParameterType), INTENT(IN) :: p_IEC                       ! parameters for IEC models   
   REAL(ReKi), INTENT(IN)  ::  GridVertCenter                  ! Position of vertical "middle" grid point (to determine if it is the hub height)
   
   ! local variables:   

   REAL(ReKi)              ::  HalfRotDiam                     ! Half of the rotor diameter
   
   REAL(ReKi)              ::  CVFA                            ! Cosine of the vertical flow angle
   REAL(ReKi)              ::  SVFA                            ! Sine of the vertical flow angle
   REAL(ReKi)              ::  CHFA                            ! Cosine of the horizontal flow angle
   REAL(ReKi)              ::  SHFA                            ! Sine of the horizontal flow angle
   
   REAL(ReKi)              ::  UTmp                            ! The best fit of observed peak Uh at het height vs jet height
   REAL(ReKi)              ::  U_zb                            ! The velocity at the bottom of the rotor disk (for estimating log fit)
   REAL(DbKi)              ::  U_zt                            ! The velocity at the top of the rotor disk (for estimating log fit)
      
   INTEGER                 :: iz, jz                           ! loop counters
   LOGICAL                 ::  HubPr                           ! Flag to indicate if the hub height is to be printed separately in the summary file
   
   ! write to the summary file:
   
   
   WRITE (US,"( // 'Turbulence Simulation Scaling Parameter Summary:' / )")
   WRITE (US,    "('   Turbulence model used                            =  ' , A )")  TRIM(TMName)

   FormStr  = "('   ',A,' =' ,F9.3,A)"
   FormStr1 = "('   ',A,' =' ,I9  ,A)"
   FormStr2 = "('   ',A,' =  ',A)"
   
   
      ! Write out a parameter summary to the summary file.

IF ( ( SpecModel  == SpecModel_IECKAI ) .OR. &
     ( SpecModel  == SpecModel_IECVKM ) .OR. &
     ( SpecModel  == SpecModel_MODVKM ) .OR. &
     ( SpecModel  == SpecModel_API    ) )  THEN  ! ADDED BY YGUO on April 192013 snow day!!!
      
      
   IF ( p_IEC%NumTurbInp ) THEN
      WRITE (US,FormStr2)      "Turbulence characteristic                       ", "User-specified"
   ELSE
      WRITE (US,FormStr2)      "Turbulence characteristic                       ", TRIM(p_IEC%IECTurbE)//p_IEC%IECTurbC
      WRITE (US,FormStr2)      "IEC turbulence type                             ", TRIM(p_IEC%IEC_WindDesc)
      
      IF ( p_IEC%IEC_WindType /= IEC_NTM ) THEN       
         WRITE (US,FormStr)    "Reference wind speed average over 10 minutes    ", p_IEC%Vref,                " m/s"
         WRITE (US,FormStr)    "Annual wind speed average at hub height         ", p_IEC%Vave,                " m/s"
      ENDIF
   ENDIF      
   
   WRITE (US,FormStr2)         "IEC standard                                    ", IECeditionSTR(p_IEC%IECedition)
   
   IF ( SpecModel  /= SpecModel_MODVKM ) THEN
      ! Write out a parameter summary to the summary file.

      WRITE (US,FormStr)       "Mean wind speed at hub height                   ", UHub,                      " m/s"

      IF (.NOT. p_IEC%NumTurbInp) THEN ! "A", "B", or "C" turbulence:
         IF ( p_IEC%IECedition == 2 ) THEN
            WRITE (US,FormStr) "Char value of turbulence intensity at 15 m/s    ", 100.0*p_IEC%TurbInt15,     "%"
            WRITE (US,FormStr) "Standard deviation slope                        ", p_IEC%SigmaSlope,          ""
         ELSE                                                                                                 
               ! This is supposed to be the expected value of what is measured at a site.                     
               ! We actually calculate the 90th percentile value to use in the code as the                    
               ! "Characteristic Value".                                                                      
            WRITE (US,FormStr) "Expected value of turbulence intensity at 15 m/s", 100.0*p_IEC%TurbInt15,           "%"
         ENDIF                                                                                                
                                                                                                              
      ENDIF                                                                                                   
                                                                                                              
      WRITE (US,FormStr)       "Characteristic value of standard deviation      ", p_IEC%SigmaIEC(1),          " m/s"
      WRITE (US,FormStr)       "Turbulence scale                                ", p_IEC%Lambda(1),            " m"
                                                                                                              
      IF ( SpecModel  == SpecModel_IECKAI )  THEN                                                                     
         WRITE (US,FormStr)    "u-component integral scale                      ", p_IEC%IntegralScale(1),     " m"
         WRITE (US,FormStr)    "Coherency scale                                 ", p_IEC%LC,                   " m"
      ELSEIF ( SpecModel  == SpecModel_IECVKM )  THEN                                                                 
         WRITE (US,FormStr)    "Isotropic integral scale                        ", p_IEC%IntegralScale(1),     " m"
      ENDIF                                                                                                   
                                                                                                              
      WRITE (US,FormStr)       "Characteristic value of hub turbulence intensity", 100.0*p_IEC%TurbInt,              "%"
                                                                                                              
   ELSE   ! ModVKM                                                                                                    
      WRITE (US,FormStr1)      "Boundary layer depth                            ", NINT(h),                   " m"
      WRITE (US,FormStr)       "Site Latitude                                   ", Latitude,                  " degs"
      WRITE (US,FormStr)       "Hub mean streamwise velocity                    ", UHub,                      " m/s"
      WRITE (US,FormStr)       "Hub local u*                                    ", UStar,                     " m/s" !BONNIE: is this LOCAL? of Disk-avg
      WRITE (US,FormStr)       "Target IEC Turbulence Intensity                 ", 100.0*p_IEC%TurbInt,             "%"
      WRITE (US,FormStr)       "Target IEC u-component standard deviation       ", p_IEC%SigmaIEC(1),         " m/s"
      WRITE (US,FormStr)       "u-component integral scale                      ", p_IEC%Lambda(1),                 " m"
      WRITE (US,FormStr)       "v-component integral scale                      ", p_IEC%Lambda(2),                 " m"
      WRITE (US,FormStr)       "w-component integral scale                      ", p_IEC%Lambda(3),                 " m"
      WRITE (US,FormStr)       "Isotropic integral scale                        ", p_IEC%LC,                        " m"
   ENDIF                                                                                                      
   WRITE (US,FormStr)          "Gradient Richardson number                      ", 0.0,                       ""

! Ustar = SigmaIEC/2.15 ! Value based on equating original Kaimal spectrum with IEC formulation

ELSEIF ( SpecModel == SpecModel_TIDAL ) THEN
   WRITE (US,FormStr2)         "Gradient Richardson number                      ", "N/A"
   WRITE (US,FormStr)          "Mean velocity at hub height                     ", UHub,                      " m/s"     
   
ELSE   
 
   WRITE (US,FormStr)          "Gradient Richardson number                      ", RICH_NO,                   ""
   WRITE (US,FormStr)          "Monin-Obukhov (M-O) z/L parameter               ", ZL,                        ""
                                                                                                              
   IF ( ZL /= 0.0 ) THEN                                                                                      
      WRITE (US,FormStr)       "Monin-Obukhov (M-O) length scale                ", L,                         " m"
   ELSE                                                                                                       
      WRITE (US,FormStr2)      "Monin-Obukhov (M-O) length scale                ", "Infinite"                 
   ENDIF                                                                                                      
   WRITE (US,FormStr)          "Mean wind speed at hub height                   ", UHub,                      " m/s"     
    
ENDIF !  TurbModel  == 'IECKAI', 'IECVKM', or 'MODVKM'


HalfRotDiam = 0.5*RotorDiameter
U_zt        = getWindSpeed(UHub,p_grid%HubHt,p_grid%HubHt+HalfRotDiam,RotorDiameter,PROFILE_TYPE=WindProfileType)   !Velocity at the top of rotor
U_zb        = getWindSpeed(UHub,p_grid%HubHt,p_grid%HubHt-HalfRotDiam,RotorDiameter,PROFILE_TYPE=WindProfileType)   !Velocity at the bottom of the rotor
      
WRITE(US,'()')   ! A BLANK LINE

SELECT CASE ( TRIM(WindProfileType) )
   CASE ('JET','J')
      PLExp = LOG( U_zt/U_zb ) / LOG( (p_grid%HubHt+HalfRotDiam)/(p_grid%HubHt-HalfRotDiam) )  !TmpReal = RotorDiameter/2
      UTmp  = 0.0422*ZJetMax+10.1979 ! Best fit of observed peak Uh at jet height vs jet height
      
      WRITE (US,FormStr2)      "Wind profile type                               ", "Low-level jet"      
      WRITE (US,FormStr)       "Jet height                                      ",  ZJetMax,                  " m"
      WRITE (US,FormStr)       "Jet wind speed                                  ",  UJetMax,                  " m/s"
      WRITE (US,FormStr)       "Upper limit of observed jet wind speed          ",  UTmp,                     " m/s"
      WRITE (US,FormStr)       "Equivalent power law exponent across rotor disk ",  PLExp,                    ""
      
      IF ( UTmp < UJetMax ) THEN
         CALL TS_Warn( 'The computed jet wind speed is larger than the ' &
                     //'maximum observed jet wind speed at this height.', .FALSE. )
      ENDIF            
                    
   CASE ('LOG','L')
      PLExp = LOG( U_zt/U_zb ) / LOG( (p_grid%HubHt+HalfRotDiam)/(p_grid%HubHt-HalfRotDiam) )  
      
      WRITE (US,FormStr2)      "Wind profile type                               ", "Logarithmic"      
      WRITE (US,FormStr)       "Equivalent power law exponent across rotor disk ",  PLExp,                    ""

   CASE ('H2L','H')
      PLExp = LOG( U_zt/U_zb ) / LOG( (p_grid%HubHt+HalfRotDiam)/(p_grid%HubHt-HalfRotDiam) ) 
      
      WRITE (US,FormStr2)      "Velocity profile type                           ", "Logarithmic (H2L)"      
      WRITE (US,FormStr)       "Equivalent power law exponent across rotor disk ",  PLExp,                    ""

   CASE ('PL','P')
      WRITE (US,FormStr2)      "Wind profile type                               ", "Power law"      
      WRITE (US,FormStr)       "Power law exponent                              ",  PLExp,                    ""
      
   CASE ('USR','U')
      PLExp = LOG( U_zt/U_zb ) / LOG( (p_grid%HubHt+HalfRotDiam)/(p_grid%HubHt-HalfRotDiam) )  

      WRITE (US,FormStr2)      "Wind profile type                               ", "Linear interpolation of user-defined profile"
      WRITE (US,FormStr)       "Equivalent power law exponent across rotor disk ",  PLExp,                    ""
                               
   CASE DEFAULT                
      WRITE (US,FormStr2)      "Wind profile type                               ", "Power law on rotor disk, logarithmic elsewhere"
      WRITE (US,FormStr)       "Power law exponent                              ",  PLExp,                    ""
      
END SELECT

WRITE(US,FormStr)              "Mean shear across rotor disk                    ", (U_zt-U_zb)/RotorDiameter, " (m/s)/m"
WRITE(US,FormStr)              "Assumed rotor diameter                          ", RotorDiameter,             " m"      
WRITE(US,FormStr)              "Surface roughness length                        ", z0,                        " m"      
WRITE(US,'()')                                                                                                 ! A BLANK LINE
WRITE(US,FormStr1)             "Number of time steps in the FFT                 ", p_grid%NumSteps,                  ""       
WRITE(US,FormStr1)             "Number of time steps output                     ", p_grid%NumOutSteps,               ""          

IF (KHtest) THEN
   WRITE(US,"(/'KH Billow Test Parameters:' / )") ! HEADER
   WRITE(US,FormStr)           "Gradient Richardson number                      ", RICH_NO,                   ""
   WRITE(US,FormStr)           "Power law exponent                              ", PLexp,                     ""
   WRITE(US,FormStr)           "Length of coherent structures                   ", UsableTime / 2.0,          " s"
   WRITE(US,FormStr)           "Minimum coherent TKE                            ", 30.0,                      " (m/s)^2"
ENDIF


   ! Write mean flow angles and wind speed profile to the summary file.

FormStr = "(//,'Mean Flow Angles:',/)"
WRITE(US,FormStr)

FormStr = "(3X,A,F6.1,' degrees')"
WRITE(US,FormStr)  'Vertical   =', VFlowAng
WRITE(US,FormStr)  'Horizontal =', HFlowAng


FormStr = "(/'Mean Wind Speed Profile:')"
WRITE(US,FormStr)

IF ( ALLOCATED( ZL_profile ) .AND. ALLOCATED( Ustar_profile ) ) THEN
   FormStr = "(/,'   Height    Wind Speed   Horizontal Angle  U-comp (X)   V-comp (Y)   W-comp (Z)   z/L(z)    u*(z)')"
   WRITE(US,FormStr)
   FormStr = "(  '     (m)        (m/s)         (degrees)       (m/s)        (m/s)        (m/s)       (-)      (m/s)')"
   WRITE(US,FormStr)
   FormStr = "(  '   ------    ----------   ----------------  ----------   ----------   ----------   ------   ------')"
   WRITE(US,FormStr)

   FormStr = '(1X,F8.1,1X,F11.2,5x,F11.2,4x,3(2X,F8.2,3X),2(1X,F8.3))'
ELSE
   FormStr = "(/,'   Height    Wind Speed   Horizontal Angle  U-comp (X)   V-comp (Y)   W-comp (Z)')"
   WRITE(US,FormStr)
   FormStr = "(  '     (m)        (m/s)         (degrees)       (m/s)        (m/s)        (m/s)   ')"
   WRITE(US,FormStr)
   FormStr = "(  '   ------    ----------   ----------------  ----------   ----------   ----------')"
   WRITE(US,FormStr)

   FormStr = '(1X,F8.1,1X,F11.2,5x,F11.2,4x,3(2X,F8.2,3X))'
ENDIF



   ! Get the angles to rotate the wind components from streamwise orientation to the X-Y-Z grid at the Hub
            
CVFA = COS( VFlowAng*D2R )
SVFA = SIN( VFlowAng*D2R ) 
CHFA = COS( HFlowAng*D2R )
SHFA = SIN( HFlowAng*D2R )

HubPr = ( ABS( p_grid%HubHt - GridVertCenter ) > Tolerance )     !If the hub height is not on the z-grid, print it, too.


   ! Write out the grid points & the hub

DO IZ = p_grid%NumGrid_Z,1, -1
   
   IF ( HubPr  .AND. ( p_grid%Z(IZ) < p_grid%HubHt ) ) THEN
   
      JZ = p_grid%NumGrid_Z+1  ! This is the index of the Hub-height parameters if the hub height is not on the grid
      
      IF ( ALLOCATED( WindDir_profile ) ) THEN      
         CHFA = COS( WindDir_profile(JZ)*D2R )
         SHFA = SIN( WindDir_profile(JZ)*D2R )
         
         IF ( ALLOCATED( ZL_profile ) ) THEN
            
            WRITE(US,FormStr)  p_grid%Z(JZ), U(JZ), WindDir_profile(JZ), U(JZ)*CHFA*CVFA, U(JZ)*SHFA*CVFA, U(JZ)*SVFA, &
                              ZL_profile(JZ), UStar_profile(JZ)
         ELSE
            WRITE(US,FormStr)  p_grid%Z(JZ), U(JZ), WindDir_profile(JZ), U(JZ)*CHFA*CVFA, U(JZ)*SHFA*CVFA, U(JZ)*SVFA
         ENDIF
      ELSE
         IF ( ALLOCATED( ZL_profile ) ) THEN
            WRITE(US,FormStr)  p_grid%Z(JZ), U(JZ), HFlowAng, U(JZ)*CHFA*CVFA, U(JZ)*SHFA*CVFA, U(JZ)*SVFA, &
                              ZL_profile(JZ), UStar_profile(JZ)
         ELSE
            WRITE(US,FormStr)  p_grid%Z(JZ), U(JZ), HFlowAng, U(JZ)*CHFA*CVFA, U(JZ)*SHFA*CVFA, U(JZ)*SVFA
         ENDIF
      ENDIF
   
      HubPr = .FALSE.
   ENDIF
   
   IF ( ALLOCATED( WindDir_profile ) ) THEN
      CHFA = COS( WindDir_profile(IZ)*D2R )
      SHFA = SIN( WindDir_profile(IZ)*D2R )

      IF ( ALLOCATED( ZL_profile ) ) THEN
         WRITE(US,FormStr)  p_grid%Z(IZ), U(IZ), WindDir_profile(IZ), U(IZ)*CHFA*CVFA, U(IZ)*SHFA*CVFA, U(IZ)*SVFA, &
                            ZL_profile(IZ), UStar_profile(IZ)
      ELSE
         WRITE(US,FormStr)  p_grid%Z(IZ), U(IZ), WindDir_profile(IZ), U(IZ)*CHFA*CVFA, U(IZ)*SHFA*CVFA, U(IZ)*SVFA
      ENDIF
   ELSE
      IF ( ALLOCATED( ZL_profile ) ) THEN
         WRITE(US,FormStr)  p_grid%Z(IZ), U(IZ), HFlowAng, U(IZ)*CHFA*CVFA, U(IZ)*SHFA*CVFA, U(IZ)*SVFA, &
                            ZL_profile(IZ), UStar_profile(IZ)
      ELSE
         WRITE(US,FormStr)  p_grid%Z(IZ), U(IZ), HFlowAng, U(IZ)*CHFA*CVFA, U(IZ)*SHFA*CVFA, U(IZ)*SVFA
      ENDIF
   ENDIF                

ENDDO ! IZ
   
   ! Write out the tower points
   
DO IZ = p_grid%NumGrid_Z,p_grid%ZLim

   IF ( p_grid%Z(IZ) < p_grid%Z(1) ) THEN
      IF ( ALLOCATED( WindDir_profile ) ) THEN
         CHFA = COS( WindDir_profile(IZ)*D2R )
         SHFA = SIN( WindDir_profile(IZ)*D2R )

         IF ( ALLOCATED( ZL_profile ) ) THEN
            WRITE(US,FormStr)  p_grid%Z(IZ), U(IZ), WindDir_profile(IZ), U(IZ)*CHFA*CVFA, U(IZ)*SHFA*CVFA, U(IZ)*SVFA, &
                               ZL_profile(IZ), UStar_profile(IZ)
         ELSE
            WRITE(US,FormStr)  p_grid%Z(IZ), U(IZ), WindDir_profile(IZ), U(IZ)*CHFA*CVFA, U(IZ)*SHFA*CVFA, U(IZ)*SVFA
         ENDIF
      ELSE
         IF ( ALLOCATED( ZL_profile ) ) THEN
            WRITE(US,FormStr)  p_grid%Z(IZ), U(IZ), HFlowAng, U(IZ)*CHFA*CVFA, U(IZ)*SHFA*CVFA, U(IZ)*SVFA, &
                               ZL_profile(IZ), UStar_profile(IZ)
         ELSE
            WRITE(US,FormStr)  p_grid%Z(IZ), U(IZ), HFlowAng, U(IZ)*CHFA*CVFA, U(IZ)*SHFA*CVFA, U(IZ)*SVFA
         ENDIF
      ENDIF                
   ENDIF

ENDDO ! IZ


END SUBROUTINE WrSum_SpecModel
!=======================================================================
SUBROUTINE WrHH_ADtxtfile(TurbInt)

USE TSMods

   REAL(ReKi), INTENT(IN)  ::  TurbInt                         ! IEC target Turbulence Intensity 


   REAL(ReKi)              ::  CVFA                            ! Cosine of the vertical flow angle
   REAL(ReKi)              ::  SVFA                            ! Sine of the vertical flow angle
   REAL(ReKi)              ::  CHFA                            ! Cosine of the horizontal flow angle
   REAL(ReKi)              ::  SHFA                            ! Sine of the horizontal flow angle

   REAL(ReKi)              :: V_Inertial(3)                    ! U,V,W components (inertial)
   REAL(ReKi)              :: UH                               ! horizontal wind speed (U+V components)

   REAL(ReKi)              :: Time                             ! The instantaneous Time (s)
   INTEGER(IntKi)          :: IT                               ! loop counter (time step)

   
   
   CHFA = COS( HH_HFlowAng*D2R )
   SHFA = SIN( HH_HFlowAng*D2R )

   CVFA = COS( VFlowAng*D2R )
   SVFA = SIN( VFlowAng*D2R )


   CALL OpenFOutFile ( UAHH, TRIM( RootName)//'.hh' )

   CALL WrScr ( ' Hub-height AeroDyn data were written to "'//TRIM( RootName )//'.hh"' )

   FormStr = "( '! This hub-height wind-speed file was generated by ' , A , A , ' on ' , A , ' at ' , A , '.' )"
   WRITE (UAHH,FormStr)   TRIM(ProgName), TRIM( ProgVer ), CurDate(), CurTime()

   WRITE (UAHH,"( '!' )")
   WRITE (UAHH,"( '! The requested statistics for this data were:' )")
   WRITE (UAHH,"( '!    Mean Total Wind Speed = ' , F8.3 , ' m/s' )")  UHub
IF ( SpecModel == SpecModel_IECKAI .OR. SpecModel == SpecModel_IECVKM .OR. SpecModel == SpecModel_MODVKM ) THEN
   WRITE (UAHH,"( '!    Turbulence Intensity  = ' , F8.3 , '%' )")  100.0*TurbInt
ENDIF
   WRITE (UAHH,"( '!' )")
   WRITE (UAHH,"( '!   Time  HorSpd  WndDir  VerSpd  HorShr  VerShr  LnVShr  GstSpd' )")
   WRITE (UAHH,"( '!  (sec)   (m/s)   (deg)   (m/s)     (-)     (-)     (-)   (m/s)' )")

   DO IT = 1, p_grid%NumOutSteps
      
      Time = p_grid%TimeStep*( IT - 1 )
            
      CALL CalculateWindComponents(V(IT,p_grid%HubIndx,:), UHub, HH_HFlowAng, VFlowAng, V_Inertial, UH)      
      
      WRITE (UAHH,'(F8.3,3F8.2,3F8.3,F8.2)')  Time, UH, -1.0*R2D*ATAN2( V_Inertial(2) , V_Inertial(1) ), &
                                                    V_Inertial(3), 0.0, PLExp, 0.0, 0.0 
!bjj: Should we output instantaneous horizontal shear, instead of 0?  
!     Should the power law exponent be an instantaneous value, too?
!                                                   
   END DO

   CLOSE(UAHH)


END SUBROUTINE WrHH_ADtxtfile
!=======================================================================
SUBROUTINE WrHH_binary(ErrStat, ErrMsg)

   ! Output HH binary turbulence parameters for GenPro analysis.
   ! Output order:  Time,U,Uh,Ut,V,W,u',v',w',u'w',u'v',v'w',TKE,CTKE.


USE TSMods

   INTEGER(IntKi),                  intent(  out) :: ErrStat                         ! Error level
   CHARACTER(*),                    intent(  out) :: ErrMsg                          ! Message describing error


   ! local variables
   
   REAL(ReKi)              :: V_Inertial(3)                    ! U,V,W components (inertial)
   REAL(ReKi)              :: UH                               ! horizontal wind speed (U+V components)
   REAL(ReKi)              :: UT                               ! total wind speed (U+V+W components)
   REAL(ReKi)              :: uv                               ! The instantaneous u'v' Reynolds stress at the hub
   REAL(ReKi)              :: uw                               ! The instantaneous u'w' Reynolds stress at the hub
   REAL(ReKi)              :: vw                               ! The instantaneous v'w' Reynolds stress at the hub
   REAL(ReKi)              :: TKE                              ! The instantaneous TKE at the hub
   REAL(ReKi)              :: CTKE                             ! The instantaneous CTKE the hub
                                                               
   REAL(ReKi)              :: Time                             ! The instantaneous Time (s)
   INTEGER(IntKi)          :: IT                               ! loop counter (time step)



   CALL OpenUOutfile ( UGTP , TRIM( RootName)//'.bin', ErrStat, ErrMsg ) 
   IF (ErrStat >= AbortErrLev) RETURN

   CALL WrScr ( ' Hub-height binary turbulence parameters were written to "'//TRIM( RootName )//'.bin"' )

   DO IT = 1, p_grid%NumOutSteps
      
      Time = p_grid%TimeStep*( IT - 1 )
      
      CALL CalculateWindComponents(V(IT,p_grid%HubIndx,:), UHub, HH_HFlowAng, VFlowAng, V_Inertial, UH, UT)
      CALL CalculateStresses(      V(IT,p_grid%HubIndx,:), uv, uw, vw, TKE, CTKE )
      
      WRITE (UGTP)  REAL(Time,SiKi), REAL(V_Inertial(1),SiKi), REAL(UH,SiKi), REAL(UT,SiKi), &      
                                     REAL(V_Inertial(2),SiKi), REAL(V_Inertial(3),SiKi), &
                                     REAL(V(IT,p_grid%HubIndx,1),SiKi), &
                                     REAL(V(IT,p_grid%HubIndx,2),SiKi), &
                                     REAL(V(IT,p_grid%HubIndx,3),SiKi), &      
                     REAL(uw,SiKi), REAL(uv,SiKi), REAL(vw,SiKi), REAL(TKE,SiKi), REAL(CTKE,SiKi)                        
                                                   
   END DO

   CLOSE(UGTP)
!WrBHHTP

END SUBROUTINE WrHH_binary
!=======================================================================
SUBROUTINE WrHH_text()

   ! Output HH text turbulence parameters
   ! Output order:  Time,U,Uh,Ut,V,W,u',v',w',u'w',u'v',v'w',TKE,CTKE.


USE TSMods

   ! local variables
   
   REAL(ReKi)              :: V_Inertial(3)                    ! U,V,W components (inertial)
   REAL(ReKi)              :: UH                               ! horizontal wind speed (U+V components)
   REAL(ReKi)              :: UT                               ! total wind speed (U+V+W components)
   REAL(ReKi)              :: uv                               ! The instantaneous u'v' Reynolds stress at the hub
   REAL(ReKi)              :: uw                               ! The instantaneous u'w' Reynolds stress at the hub
   REAL(ReKi)              :: vw                               ! The instantaneous v'w' Reynolds stress at the hub
   REAL(ReKi)              :: TKE                              ! The instantaneous TKE at the hub
   REAL(ReKi)              :: CTKE                             ! The instantaneous CTKE the hub
                                                               
   REAL(ReKi)              :: Time                             ! The instantaneous Time (s)
   INTEGER(IntKi)          :: IT                               ! loop counter (time step)


   ! WrFHHTP

   CALL OpenFOutFile ( UFTP, TRIM( RootName)//'.dat' )

   CALL WrScr ( ' Hub-height formatted turbulence parameters were written to "'//TRIM( RootName )//'.dat"' )

   WRITE (UFTP,"( / 'This hub-height turbulence-parameter file was generated by ' , A , A , ' on ' , A , ' at ' , A , '.' / )") &
                TRIM(ProgName), TRIM( ProgVer ), CurDate(), CurTime()
 
   WRITE (UFTP,"('   Time',6X,'U',7X,'Uh',7X,'Ut',8X,'V',8X,'W',8X,'u''',7X,'v''',7X,'w'''," &
                        //"6X,'u''w''',5X,'u''v''',5X,'v''w''',5X,'TKE',6X,'CTKE')")   
   
   
   DO IT = 1, p_grid%NumOutSteps
      
      Time = p_grid%TimeStep*( IT - 1 )
      
      CALL CalculateWindComponents(V(IT,p_grid%HubIndx,:), UHub, HH_HFlowAng, VFlowAng, V_Inertial, UH, UT)
      CALL CalculateStresses(      V(IT,p_grid%HubIndx,:), uv, uw, vw, TKE, CTKE )
      
                             
      WRITE(UFTP,'(F7.2,13F9.3)') Time,V_Inertial(1),UH,UT,V_Inertial(2),V_Inertial(3), &
                                    V(IT,p_grid%HubIndx,1), V(IT,p_grid%HubIndx,2), V(IT,p_grid%HubIndx,3), &
                                    uw, uv, vw, TKE, CTKE
      
      
   END DO

   CLOSE(UFTP)

END SUBROUTINE WrHH_text
!=======================================================================
SUBROUTINE WrSum_Stats(USig, VSig, WSig, UXBar, UXSig, ErrStat, ErrMsg)


use TSMods

REAL(ReKi), INTENT(OUT) ::  USig                            ! Standard deviation of the u-component wind speed at the hub
REAL(ReKi), INTENT(OUT) ::  VSig                            ! Standard deviation of the v-component wind speed at the hub
REAL(ReKi), INTENT(OUT) ::  WSig                            ! Standard deviation of the w-component wind speed at the hub
REAL(ReKi), INTENT(OUT) ::  UXSig                           ! Standard deviation of the U-component wind speed at the hub
REAL(DbKi), INTENT(OUT) ::  UXBar                           ! The mean U-component (u rotated; x-direction) wind speed at the hub

   INTEGER(IntKi),                  intent(  out) :: ErrStat                         ! Error level
   CHARACTER(*),                    intent(  out) :: ErrMsg                          ! Message describing error


REAL(DbKi)              ::  SumS                            ! Sum of the velocity-squared, used for calculating standard deviations in the summary file
REAL(DbKi)              ::  UBar                            ! The mean u-component wind speed at the hub
REAL(DbKi)              ::  UHBar                           ! The mean horizontal wind speed at the hub
REAL(DbKi)              ::  UHSum2                          ! The sum of the squared horizontal wind speed at the hub
REAL(DbKi)              ::  UHTmp                           ! The instantaneous horizontal wind speed at the hub
REAL(DbKi)              ::  UHTmp2                          ! The instantaneous squared horizontal wind speed at the hub
REAL(DbKi)              ::  USum2                           ! The sum of the squared u-component wind speed at the hub
REAL(DbKi)              ::  UTBar                           ! The mean total wind speed at the hub
REAL(DbKi)              ::  UTmp                            ! The instantaneous u-component wind speed at the hub
REAL(DbKi)              ::  UTmp2                           ! The instantaneous squared u-component wind speed at the hub
REAL(DbKi)              ::  UTSum2                          ! The sum of the squared total wind speed at the hub
REAL(DbKi)              ::  UTTmp                           ! The instantaneous total wind speed at the hub
REAL(DbKi)              ::  UTTmp2                          ! The instantaneous squared total wind speed at the hub
REAL(DbKi)              ::  UXSum                           ! The sum of the U-component (u rotated) wind speed at the hub
REAL(DbKi)              ::  UXSum2                          ! The sum of the squared U-component (u rotated) wind speed at the hub
REAL(DbKi)              ::  UXTmp                           ! The instantaneous U-component (u rotated) wind speed at the hub
REAL(DbKi)              ::  UXTmp2                          ! The instantaneous squared U-component (u rotated) wind speed at the hub
REAL(DbKi)              ::  UYBar                           ! The mean V-component (v rotated; y-direction) wind speed at the hub
REAL(DbKi)              ::  UYSum                           ! The sum of the V-component (v rotated) wind speed at the hub
REAL(DbKi)              ::  UYSum2                          ! The sum of the squared V-component (v rotated) wind speed at the hub
REAL(DbKi)              ::  UYTmp                           ! The instantaneous V-component (v rotated) wind speed at the hub
REAL(DbKi)              ::  UYTmp2                          ! The instantaneous squared V-component (v rotated) wind speed at the hub
REAL(DbKi)              ::  UZBar                           ! The mean W-component (w rotated; z-direction) wind speed at the hub
REAL(DbKi)              ::  UZSum                           ! The sum of the W-component (w rotated) wind speed at the hub
REAL(DbKi)              ::  UZSum2                          ! The sum of the squared W-component (w rotated) wind speed at the hub
REAL(DbKi)              ::  UZTmp                           ! The instantaneous W-component (w rotated) wind speed at the hub
REAL(DbKi)              ::  UZTmp2                          ! The instantaneous squared W-component (w rotated) wind speed at the hub
REAL(DbKi)              ::  VBar                            ! The mean v-component wind speed at the hub
REAL(DbKi)              ::  VSum2                           ! The sum of the squared v-component wind speed at the hub
REAL(DbKi)              ::  VTmp                            ! The instantaneous v-component wind speed at the hub
REAL(DbKi)              ::  VTmp2                           ! The instantaneous squared v-component wind speed at the hub
REAL(DbKi)              ::  WBar                            ! The mean w-component wind speed at the hub
REAL(DbKi)              ::  WSum2                           ! The sum of the squared w-component wind speed at the hub
REAL(DbKi)              ::  WTmp                            ! The instantaneous w-component wind speed at the hub
REAL(DbKi)              ::  WTmp2                           ! The instantaneous squared w-component wind speed at the hub

REAL(ReKi)              ::  CHFA                            ! Cosine of the Horizontal Flow Angle
REAL(ReKi)              ::  CVFA                            ! Cosine of the Vertical Flow Angle
REAL(ReKi)              ::  SVFA                            ! Sine of the Vertical Flow Angle
REAL(ReKi)              ::  SHFA                            ! Sine of the Horizontal Flow Angle
REAL(ReKi)              ::  CTKEmax                         ! Maximum instantaneous Coherent Turbulent Kenetic Energy at the hub
REAL(ReKi)              ::  TKEmax                          ! Maximum instantaneous Turbulent Kenetic Energy at the hub


REAL(ReKi)              ::  UHmax                           ! Maximum horizontal wind speed at the hub
REAL(ReKi)              ::  UHmin                           ! Minimum horizontal wind speed at the hub
REAL(ReKi)              ::  Umax                            ! Maximum u-component wind speed at the hub
REAL(ReKi)              ::  Umin                            ! Minimum u-component wind speed at the hub
REAL(ReKi)              ::  UTSig                           ! Standard deviation of the total wind speed at the hub
REAL(ReKi)              ::  UT_TI                           ! Turbulent Intensity of the total wind speed at the hub
REAL(ReKi)              ::  UTmax                           ! Maximum total wind speed at the hub
REAL(ReKi)              ::  UTmin                           ! Minimum total wind speed at the hub
REAL(ReKi)              ::  UVMax                           ! Maximum u'v' Reynolds Stress at the hub
REAL(ReKi)              ::  UVMin                           ! Minimum u'v' Reynolds Stress at the hub
REAL(ReKi)              ::  UVTmp                           ! The instantaneous u'v' Reynolds stress at the hub
REAL(ReKi)              ::  UV_RS                           ! The average u'v' Reynolds stress at the hub
REAL(ReKi)              ::  UVcor                           ! The u-v cross component correlation coefficient at the hub
REAL(ReKi)              ::  UVsum                           ! The sum of the u'v' Reynolds stress component at the hub
REAL(ReKi)              ::  UWMax                           ! Maximum u'w' Reynolds Stress at the hub
REAL(ReKi)              ::  UWMin                           ! Minimum u'w' Reynolds Stress at the hub
REAL(ReKi)              ::  UWTmp                           ! The instantaneous u'w' Reynolds stress at the hub
REAL(ReKi)              ::  UW_RS                           ! The average u'w' Reynolds stress at the hub
REAL(ReKi)              ::  UWcor                           ! The u-w cross component correlation coefficient at the hub
REAL(ReKi)              ::  UWsum                           ! The sum of the u'w' Reynolds stress component at the hub
REAL(ReKi)              ::  UXmax                           ! Maximum U-component (X-direction) wind speed at the hub
REAL(ReKi)              ::  UXmin                           ! Minimum U-component wind speed at the hub
REAL(ReKi)              ::  UYmax                           ! Maximum V-component (Y-direction) wind speed at the hub
REAL(ReKi)              ::  UYmin                           ! Minimum V-component wind speed at the hub
REAL(ReKi)              ::  UYSig                           ! Standard deviation of the V-component wind speed at the hub
REAL(ReKi)              ::  UZmax                           ! Maximum W-component (Z-direction) wind speed at the hub
REAL(ReKi)              ::  UZmin                           ! Minimum W-component wind speed at the hub
REAL(ReKi)              ::  UZSig                           ! Standard deviation of the W-component wind speed at the hub
REAL(ReKi)              ::  U_TI                            ! The u-component turbulence intensity at the hub
REAL(ReKi)              ::  Vmax                            ! Maximum v-component wind speed at the hub
REAL(ReKi)              ::  Vmin                            ! Minimum v-component wind speed at the hub
REAL(ReKi)              ::  VWMax                           ! Maximum v'w' Reynolds Stress at the hub
REAL(ReKi)              ::  VWMin                           ! Minimum v'w' Reynolds Stress at the hub
REAL(ReKi)              ::  VWTmp                           ! The instantaneous v'w' Reynolds stress at the hub
REAL(ReKi)              ::  VW_RS                           ! The average v'w' Reynolds stress at the hub
REAL(ReKi)              ::  VWcor                           ! The v-w cross component correlation coefficient at the hub
REAL(ReKi)              ::  VWsum                           ! The sum of the v'w' Reynolds stress component at the hub
REAL(ReKi)              ::  V_TI                            ! The v-component turbulence intensity at the hub
REAL(ReKi)              ::  Wmax                            ! Maximum w-component wind speed at the hub
REAL(ReKi)              ::  Wmin                            ! Minimum w-component wind speed at the hub
REAL(ReKi)              ::  W_TI                            ! The w-component turbulence intensity at the hub

REAL(ReKi)              ::  CTKE                            ! Coherent Turbulent Kenetic Energy at the hub
REAL(ReKi)              ::  TKE                             ! Turbulent Kenetic Energy at the hub
REAL(ReKi)              ::  INumSteps                       ! Multiplicative Inverse of the Number of time Steps
REAL(ReKi)              ::  UHSig                           ! Approximate sigma of the horizontal wind speed at the hub point
REAL(ReKi)              ::  SUstar                          ! Simulated U-star at the hub

REAL(ReKi)              ::  UH_TI                           ! TI of the horizontal wind speed at the hub point


INTEGER(IntKi)          :: IT, IVec, IY, IZ, II

   ! Initialize statistical quantities for hub-height turbulence parameters.

   CALL WrScr ( ' Computing hub-height statistics' )

   INumSteps   = 1.0/p_grid%NumSteps

   CTKEmax = -HUGE( CTKEmax )
   TKEmax  = -HUGE(  TKEmax )
   UBar    =      0.0
   UHBar   =      0.0
   UHmax   = -HUGE( UHmax )
   UHmin   =  HUGE( UHmin )
   UHSum2  =      0.0
   Umax    = -HUGE( Umax )
   Umin    =  HUGE( Umin )
   USum2   =      0.0
   UTBar   =      0.0
   UTmax   = -HUGE( UTmax )
   UTmin   =  HUGE( UTmin )
   UTSum2  =      0.0
   UV_RS   =      0.0
   UVMax   =  V(1,p_grid%HubIndx,1)*V(1,p_grid%HubIndx,2) 
   UVMin   =  HUGE( UVMin )
   UVsum   =      0.0
   UW_RS   =      0.0
   UWMax   =  V(1,p_grid%HubIndx,1)*V(1,p_grid%HubIndx,3) 
   UWMin   =  HUGE( UWMin )
   UWsum   =      0.0
   VBar    =      0.0
   Vmax    = -HUGE( Vmax )
   Vmin    =  HUGE( Vmin )
   VSum2   =      0.0
   VW_RS   =      0.0
   VWMax   =  V(1,p_grid%HubIndx,2)*V(1,p_grid%HubIndx,3)
   VWMin   =  HUGE( VWMin )
   VWsum   =      0.0
   WBar    =      0.0
   Wmax    = -HUGE( Wmax )
   Wmin    =  HUGE( Wmin )
   WSum2   =      0.0
   UXBar   =      0.0
   UXmax   = -HUGE( UXmax )
   UXmin   =  HUGE( UXmin )
   UXSum   =      0.0
   UXSum2  =      0.0
   UXTmp   =      0.0
   UXTmp2  =      0.0
   UYBar   =      0.0
   UYmax   = -HUGE( UYmax )
   UYmin   =  HUGE( UYmin )
   UYSum   =      0.0
   UYSum2  =      0.0
   UYTmp   =      0.0
   UYTmp2  =      0.0
   UZBar   =      0.0
   UZmax   = -HUGE( UZmax )
   UZmin   =  HUGE( UZmin )
   UZSum   =      0.0
   UZSum2  =      0.0
   UZTmp   =      0.0
   UZTmp2  =      0.0

   CHFA = COS( HH_HFlowAng*D2R )
   SHFA = SIN( HH_HFlowAng*D2R )

   CVFA = COS( VFlowAng*D2R )
   SVFA = SIN( VFlowAng*D2R )

   DO IT=1,p_grid%NumSteps

         ! Calculate longitudinal (UTmp), lateral (VTmp), and upward (WTmp)
         ! values for hub station, as well as rotated (XTmp, YTmp, ZTmp) 
         ! components applying specified flow angles.

         ! Add mean wind speed to the streamwise component
      UTmp = V(IT,p_grid%HubIndx,1) + UHub
      VTmp = V(IT,p_grid%HubIndx,2)
      WTmp = V(IT,p_grid%HubIndx,3)
   
         ! Rotate the wind components from streamwise orientation to the X-Y-Z grid at the Hub      
      UXTmp = UTmp*CHFA*CVFA - VTmp*SHFA - WTmp*CHFA*SVFA
      UYTmp = UTmp*SHFA*CVFA + VTmp*CHFA - WTmp*SHFA*SVFA  
      UZTmp = UTmp*SVFA                  + WTmp*CVFA

         ! Calculate hub horizontal wind speed (UHTmp) and Total wind speed (UTTmp)
      UTmp2 = UTmp*UTmp          !flow coordinates
      VTmp2 = VTmp*VTmp
      WTmp2 = WTmp*WTmp
   
      UXTmp2 = UXTmp*UXTmp       !inertial frame coordinates
      UYTmp2 = UYTmp*UYTmp
      UZTmp2 = UZTmp*UZTmp
   
      UHTmp2 = UXTmp2 + UYTmp2   !inertial frame coordinates
      UTTmp2 = UHTmp2 + UZTmp2

      UHTmp = SQRT( UHTmp2 )     !inertial frame coordinates
      UTTmp = SQRT( UTTmp2 )

         ! Form running sums for hub standard deviations

      UBar   = UBar   + UTmp     !flow coordinates
      VBar   = VBar   + VTmp     !flow coordinates
      WBar   = WBar   + WTmp     !flow coordinates

      USum2  = USum2  + UTmp2    !flow coordinates
      VSum2  = VSum2  + VTmp2    !flow coordinates
      WSum2  = WSum2  + WTmp2    !flow coordinates

      UXBar   = UXBar + UXTmp
      UYBar   = UYBar + UYTmp
      UZBar   = UZBar + UZTmp

      UXSum2  = UXSum2 + UXTmp2
      UYSum2  = UYSum2 + UYTmp2
      UZSum2  = UZSum2 + UZTmp2

      UHBar  = UHBar  + UHTmp
      UTBar  = UTBar  + UTTmp

      UHSum2 = UHSum2 + UHTmp2
      UTSum2 = UTSum2 + UTTmp2


         ! Determine hub extremes.
      
      IF ( UTmp  > Umax  )  Umax  = UTmp     !flow coordinates,
      IF ( UTmp  < Umin  )  Umin  = UTmp     !flow coordinates,

      IF ( VTmp  > Vmax  )  Vmax  = VTmp     !flow coordinates,
      IF ( VTmp  < Vmin  )  Vmin  = VTmp     !flow coordinates,

      IF ( WTmp  > Wmax  )  Wmax  = WTmp     !flow coordinates,
      IF ( WTmp  < Wmin  )  Wmin  = WTmp     !flow coordinates,

      IF ( UXTmp > UXmax )  UXmax = UXTmp
      IF ( UXTmp < UXmin )  UXmin = UXTmp

      IF ( UYTmp > UYmax )  UYmax = UYTmp
      IF ( UYTmp < UYmin )  UYmin = UYTmp

      IF ( UZTmp > UZmax )  UZmax = UZTmp
      IF ( UZTmp < UZmin )  UZmin = UZTmp

      IF ( UHTmp > UHmax )  UHmax = UHTmp
      IF ( UHTmp < UHmin )  UHmin = UHTmp

      IF ( UTTmp > UTmax )  UTmax = UTTmp
      IF ( UTTmp < UTmin )  UTmin = UTTmp

         ! Find maxes and mins of instantaneous hub Reynolds stresses u'w', u'v', and v'w'

      UVTmp = V(IT,p_grid%HubIndx,1)*V(IT,p_grid%HubIndx,2)
      UWTmp = V(IT,p_grid%HubIndx,1)*V(IT,p_grid%HubIndx,3)
      VWTmp = V(IT,p_grid%HubIndx,2)*V(IT,p_grid%HubIndx,3)

      IF     ( UVTmp < UVMin )  THEN
         UVMin = UVTmp
      ELSEIF ( UVTmp > UVMax )  THEN
         UVMax = UVTmp
      ENDIF

      IF     ( UWTmp < UWMin )  THEN
         UWMin = UWTmp
      ELSEIF ( UWTmp > UWMax )  THEN
         UWMax = UWTmp
      ENDIF

      IF     ( VWTmp < VWMin )  THEN
         VWMin = VWTmp
      ELSEIF ( VWTmp > VWMax )  THEN
         VWMax = VWTmp
      ENDIF

         ! Find maximum of instantaneous TKE and CTKE.

      TKE  = 0.5*(V(IT,p_grid%HubIndx,1)*V(IT,p_grid%HubIndx,1) + V(IT,p_grid%HubIndx,2)*V(IT,p_grid%HubIndx,2) + V(IT,p_grid%HubIndx,3)*V(IT,p_grid%HubIndx,3))
      CTKE = 0.5*SQRT(UVTmp*UVTmp + UWTmp*UWTmp + VWTmp*VWTmp)

      IF (CTKE > CTKEmax) CTKEmax = CTKE
      IF ( TKE >  TKEmax)  TKEmax =  TKE

         ! Find sums for mean and square Reynolds stresses for hub-level simulation.
      UVsum = UVsum + UVTmp
      UWsum = UWsum + UWTmp
      VWsum = VWsum + VWTmp

   ENDDO ! IT



      ! Calculate mean hub-height Reynolds stresses.
   UW_RS =  UWsum*INumSteps
   UV_RS =  UVsum*INumSteps
   VW_RS =  VWsum*INumSteps

      ! Simulated Hub UStar.
   SUstar = SQRT( ABS( UW_RS ) )

      ! Calculate mean values for hub station.

   UBar = UBar*INumSteps
   VBar = VBar*INumSteps
   WBar = WBar*INumSteps

   UXBar = UXBar*INumSteps
   UYBar = UYBar*INumSteps
   UZBar = UZBar*INumSteps

   UHBar = UHBar*INumSteps
   UTBar = UTBar*INumSteps


      ! Calculate the standard deviations for hub station.
      ! (SNWind/SNLwind-3D) NOTE: This algorithm is the approximate algorithm.
      ! bjj: do the algebra and you'll find that it's std() using the 1/n definition   

   USig  = SQRT( MAX( USum2 *INumSteps-UBar *UBar , 0.0_DbKi ) )
   VSig  = SQRT( MAX( VSum2 *INumSteps-VBar *VBar , 0.0_DbKi ) )
   WSig  = SQRT( MAX( WSum2 *INumSteps-WBar *WBar , 0.0_DbKi ) )

   UXSig = SQRT( MAX( UXSum2*INumSteps-UXBar*UXBar, 0.0_DbKi ) )
   UYSig = SQRT( MAX( UYSum2*INumSteps-UYBar*UYBar, 0.0_DbKi ) )
   UZSig = SQRT( MAX( UZSum2*INumSteps-UZBar*UZBar, 0.0_DbKi ) )

   UHSig = SQRT( MAX( UHSum2*INumSteps-UHBar*UHBar, 0.0_DbKi ) )
   UTSig = SQRT( MAX( UTSum2*INumSteps-UTBar*UTBar, 0.0_DbKi ) )


      ! Calculate Cross-component correlation coefficients
   UWcor = ( UW_RS ) / (USig * WSig) ! this definition assumes u' and w' have zero mean
   UVcor = ( UV_RS ) / (USig * VSig)
   VWcor = ( VW_RS ) / (VSig * WSig)


      ! Calculate turbulence intensities.
   U_TI = USig/UBar    
   V_TI = VSig/UBar
   W_TI = WSig/UBar

   UH_TI = UHSig/UHBar
   UT_TI = UTSig/UTBar


      ! Write out the hub-level stats to the summary file.

   CALL WrScr ( ' Writing statistics to summary file' )

   FormStr = "(//,'Hub-Height Simulated Turbulence Statistical Summary:')"
   WRITE(US,FormStr)

   FormStr = "(/,3X,'Type of Wind        Min (m/s)   Mean (m/s)    Max (m/s)  Sigma (m/s)       TI (%)')"
   WRITE(US,FormStr)

   FormStr = "(  3X,'----------------    ---------   ----------    ---------  -----------       ------')"
   WRITE(US,FormStr)

   FormStr = "(3X,A,F13.2,2F13.2,2F13.3)"
   !bjj for analysis, extra precision: FormStr = "(3X,A,F13.2,2F13.2,2F13.6)"

   WRITE (US,FormStr)  'Longitudinal (u)',  Umin,  UBar,  Umax,  USig, 100.0* U_TI
   WRITE (US,FormStr)  'Lateral (v)     ',  Vmin,  VBar,  Vmax,  VSig, 100.0* V_TI
   WRITE (US,FormStr)  'Vertical (w)    ',  Wmin,  WBar,  Wmax,  WSig, 100.0* W_TI
   WRITE (US,FormStr)  'U component     ', UXmin, UXBar, UXmax, UXSig, 100.0*UXSig/UXBar
   WRITE (US,FormStr)  'V component     ', UYmin, UYBar, UYmax, UYSig, 100.0*UYSig/UXBar
   WRITE (US,FormStr)  'W component     ', UZmin, UZBar, UZmax, UZSig, 100.0*UZSig/UXBar
   WRITE (US,FormStr)  'Horizontal (U&V)', UHmin, UHBar, UHmax, UHSig, 100.0*UH_TI
   WRITE (US,FormStr)  'Total           ', UTmin, UTBar, UTmax, UTSig, 100.0*UT_TI

   WRITE(US,"(/,3X,'                    Min Reynolds     Mean Reynolds    Max Reynolds    Correlation')")
   WRITE(US,"(  3X,'Product             Stress (m/s)^2   Stress (m/s)^2   Stress (m/s)^2  Coefficient')")
   WRITE(US,"(  3X,'----------------    --------------   --------------   --------------  -----------')")

   FormStr = "(3X,A,3(3X,F12.3,3X),F11.3)"
   WRITE (US,FormStr)  "u'w'            ",  UWMin,  UW_RS, UWMax, UWcor
   WRITE (US,FormStr)  "u'v'            ",  UVMin,  UV_RS, UVMax, UVcor
   WRITE (US,FormStr)  "v'w'            ",  VWMin,  VW_RS, VWMax, VWcor

   FormStr = "(3X,A,' = ',F10.3,A)"
   WRITE(US,"(/)")   ! blank line
   WRITE(US,FormStr)  "Friction Velocity (Ustar) ", SUstar,  " m/s"
   WRITE(US,FormStr)  "Maximum Instantaneous TKE ", TKEmax,  " (m/s)^2"
   WRITE(US,FormStr)  "Maximum Instantaneous CTKE", CTKEmax, " (m/s)^2"

      !  Allocate the array of standard deviations.

   CALL AllocAry( SDary, p_grid%NumGrid_Y, 'SDary (standard deviations)', ErrStat, ErrMsg)
   IF (ErrStat >= AbortErrLev) RETURN


      ! Calculate standard deviations for each grid point.  Write them to summary file.

   WRITE(US,"(//,'Grid Point Variance Summary:',/)")
   WRITE(US,"(3X,'Y-coord',"//TRIM(Num2LStr(p_grid%NumGrid_Y))//"F8.2)")  p_grid%Y(1:p_grid%NumGrid_Y)

   FormStr1 = "(/,3X'Height   Standard deviation at grid points for the ',A,' component:')"
   FormStr2 = "(F9.2,1X,"//TRIM(Num2LStr(p_grid%NumGrid_Y))//"F8.3)"

   UTmp = 0
   VTmp = 0
   WTmp = 0

   DO IVec=1,3

      WRITE(US,FormStr1)  Comp(IVec)

      DO IZ=p_grid%NumGrid_Z,1,-1

         DO IY=1,p_grid%NumGrid_Y

            II   = (IZ-1)*p_grid%NumGrid_Y+IY
            SumS = 0.0

            DO IT=1,p_grid%NumSteps
               SumS = SumS + V(IT,II,IVec)**2            
            ENDDO ! IT         

            SDary(IY) = SQRT(SumS*INumSteps)   !  Was:  SDary(IZ,IY) = SQRT(SumS*INumSteps)/U(IZ,NumGrid/2)

         ENDDO ! IY

         WRITE(US,FormStr2) p_grid%Z(IZ), SDary(1:p_grid%NumGrid_Y)

         IF ( IVec == 1 ) THEN
            UTmp = UTmp + SUM( SDary )
         ELSEIF ( IVec == 2 ) THEN
            VTmp = VTmp + SUM( SDary )
         ELSE
            WTmp = WTmp + SUM( SDary )
         ENDIF
      ENDDO ! IZ

   ENDDO ! Ivec

   FormStr2 = "(6X,A,' component: ',F8.3,' m/s')"
   WRITE(US,"(/'   Mean standard deviation across all grid points:')")
   WRITE(US,FormStr2) Comp(1), UTmp / ( p_grid%NumGrid_Y*p_grid%NumGrid_Z )
   WRITE(US,FormStr2) Comp(2), VTmp / ( p_grid%NumGrid_Y*p_grid%NumGrid_Z ) 
   WRITE(US,FormStr2) Comp(3), WTmp / ( p_grid%NumGrid_Y*p_grid%NumGrid_Z ) 


      !  Deallocate the array of standard deviations.

   IF ( ALLOCATED( SDary ) )  DEALLOCATE( SDary )



END SUBROUTINE WrSum_Stats
!=======================================================================

!=======================================================================
END MODULE TS_FileIO
