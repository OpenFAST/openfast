!=======================================================================
PROGRAM TurbSim

   ! A turbulence simulator developed by contractors and staff at the
   ! National Renewable Energy Laboratory, Golden, Colorado
   !
   ! v1.0a-bjj  15-Mar-2004  B. Jonkman  
   ! v1.0        4-Nov-2005  B. Jonkman
   ! v1.01      24-Jan-2006  B. Jonkman  (NWTC subs v1.00a-mlb)
   ! v1.10      10-Apr-2006  B. Jonkman  (NWTC subs v1.12)
   ! v1.20      20-Oct-2006  B. Jonkman  (NWTC subs v1.12)
   ! v1.21       1-Feb-2007  B. Jonkman  (NWTC subs v1.12)
   ! v1.30       4-Apr-2008  B. Jonkman  (NWTC subs v1.01.09)
   ! v1.40      12-Sep-2008  B. Jonkman  (NWTC subs v1.01.09)
   ! v1.41h     11-Jun-2009  B. Jonkman  (NWTC subs v1.01.09)
   ! v1.50      25-Sep-2009  B. Jonkman  (NWTC subs v1.01.09)
   ! V1.06.00   21-Sep-2012  L. Kilcher & B. Jonkman (NWTC Library v1.04.01)
   !
   ! This program simulates a full field of turbulent winds at points in space
   ! in a rectangular Cartesian plane perpendicular to the mean wind direction.
   ! It emulates, as closely as possible, the specifications of the IEC
   ! Document 61400-1 with a choice of the Kaimal or von Karman spectral models.
   ! It also includes spectral models of smooth terrain, complex terrain, and 
   ! flow in and around a wind farm.
   !
   ! This program is based upon, and includes, code developed by Paul Veers
   ! of Sandia National Laboratories (SNLWIND), Neil Kelley of NREL
   ! (SNLWIND-3D and SNLWIND-IEC) and Marshall Buhl, also of NREL (SNWIND).
   !
   ! Example usage:  at a command prompt, type 
   !     turbsim turbsim.inp
   !     This will read the text input file named "turbsim.inp"

   ! -------------------------------------------------------------------------------------------------------
   ! This is for dimension 2 of V(:,:,:)....
   ! 
   ! The grid of points on the Cartesian plane is numbered in the following way (notice that the first 
   ! height starts at the bottom of the grid): 
   !
   !               Y(1)                        Y(2)             Y(3)             ... Y(NumGrid_Y)
   !              -------------------------------------------------------------------------------------------
   ! Z(NumGrid_Z):|V(NumGrid_Y*(NumGrid_Z-1)+1) ...                                  V(NumGrid_Y*NumGrid_Z) |
   ! ...          |...                                                                                      |
   ! Z(2)        :|V(NumGrid_Y + 1)            V(NumGrid_Y + 2) V(NumGrid_Y + 3) ... V(2*NumGrid_Y)         |
   ! Z(1)        :|V(1)                        V(2)             V(3)             ... V(NumGrid_Y)           |
   !              -------------------------------------------------------------------------------------------
   ! 
   ! Z(i) < Z(i+1) for all integers i, 1 <= i < NumGrid_Z
   ! Y(j) < Y(j+1) for all integers j, 1 <= j < NumGrid_Y
   !
   ! If an extra hub point is necessary because the point does not fall on the grid, 
   ! then it is added immediately following the regular grid points, i.e. 
   ! Hub point = NumGrid_Y * NumGrid_Z + 1.
   !
   ! If the tower wind file output is selected, those extra points (in a single vertical 
   ! line) are added at the end, after the hub point.
   ! --------------------------------------------------------------------------------------------------------


USE                        FFT_Module
USE                        Ran_Lux_Mod
USE                        TSMods
USE                        TSsubs

!BONNIE:*****************************
! USE    IFPORT, ONLY: TIMEF ! Wall Clock Time
!BONNIE:*****************************    


IMPLICIT                   NONE


   ! Declare local variables

REAL(DbKi)              ::  CGridSum                        ! The sums of the velocity components at the points surrounding the hub (or at the hub if it's on the grid)
REAL(DbKi)              ::  CGridSum2                       ! The sums of the squared velocity components at the points surrouding the hub 
REAL(DbKi), PARAMETER   ::  DblZero = 0.0                   ! Zero, used to compare 8-byte reals against zero in MAX() functions (to conform to F95 standards)
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
REAL(DbKi)              ::  UXBar                           ! The mean U-component (u rotated; x-direction) wind speed at the hub
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
REAL(ReKi)              ::  CPUtime                         ! Contains the number of seconds since the start of the program
REAL(ReKi)              ::  CVFA                            ! Cosine of the Vertical Flow Angle
REAL(ReKi)              ::  CTKE                            ! Coherent Turbulent Kenetic Energy at the hub
REAL(ReKi)              ::  CTKEmax                         ! Maximum instantaneous Coherent Turbulent Kenetic Energy at the hub
REAL(ReKi)              ::  DelF5                           ! half of the delta frequency, used to discretize the continuous PSD at each point
REAL(ReKi)              ::  DelF                            ! Delta frequency
REAL(ReKi)              ::  INumSteps                       ! Multiplicative Inverse of the Number of time Steps
REAL(ReKi), PARAMETER   ::  KHT_LES_dT = 0.036335           ! The average time step in the LES test file, used here for the KH test
REAL(ReKi), PARAMETER   ::  KHT_LES_Zm = 6.35475            ! The non-dimensional z dimension defined in LES test file, used here for the KH test
REAL(ReKi)              ::  LC                              ! IEC coherency scale parameter
REAL(ReKi)              ::  Lambda1                         ! IEC turbulence scale parameter
REAL(ReKi)              ::  MaxEvtCTKE                      ! Maximum non-dimensional CTKE in events (used for scaling coherent events to appropriate size)
REAL(ReKi)              ::  ScaleWid                        ! Scaling width for LE coherent turbulence (RotDiam in AeroDyn FD_Wind)
REAL(ReKi)              ::  ScaleVel                        ! Scaling velocity for LE coherent turbulence, U0.  2*U0 is the difference in wind speed between the top and bottom of the wave.
REAL(ReKi)              ::  SHFA                            ! Sine of the Horizontal Flow Angle
REAL(ReKi)              ::  SigmaSlope                      ! Slope used with IEC models to determine target sigma and turbulent intensity
REAL(ReKi)              ::  SumS                            ! Sum of the velocity-squared, used for calculating standard deviations in the summary file
REAL(ReKi)              ::  SUstar                          ! Simulated U-star at the hub
REAL(ReKi)              ::  SVFA                            ! Sine of the Vertical Flow Angle
REAL(ReKi)              ::  TargetSigma                     ! The target standard deviation for the IEC models, if scaling is requested.
REAL(ReKi)              ::  Time                            ! The time at each step (Used in the binary HH file for GenPro)
REAL(ReKi)              ::  TKE                             ! Turbulent Kenetic Energy at the hub
REAL(ReKi)              ::  TKEmax                          ! Maximum instantaneous Turbulent Kenetic Energy at the hub
!REAL(ReKi)             ::  TmpPLExp                        ! Temporarily holds the Power Law Exponent for KH scaling
REAL(ReKi)              ::  TmpU                            ! Temporarily holds the value of the u component
REAL(ReKi)              ::  TmpV                            ! Temporarily holds the value of the v component
REAL(ReKi)              ::  TmpW                            ! Temporarily holds the value of the w component
REAL(ReKi)              ::  TmpY                            ! Temp variable for interpolated hub point
REAL(ReKi)              ::  TmpZ                            ! Temp variable for interpolated hub point
REAL(ReKi)              ::  Tmp_YL_Z                        ! Temp variable for interpolated hub point
REAL(ReKi)              ::  Tmp_YH_Z                        ! Temp variable for interpolated hub point
REAL(ReKi)              ::  TurbInt                         ! IEC target Turbulence Intensity 
REAL(ReKi)              ::  TurbInt15                       ! Turbulence Intensity at hub height with a mean wind speed of 15 m/s
REAL(ReKi)              ::  UGridMean                       ! Average wind speed at the points surrounding the hub
REAL(ReKi)              ::  UGridSig                        ! Standard deviation of the wind speed at the points surrounding the hub
REAL(ReKi)              ::  UGridTI                         ! Turbulent Intensity of the points surrounding the hub
REAL(ReKi)              ::  UHSig                           ! Approximate sigma of the horizontal wind speed at the hub point
REAL(ReKi)              ::  UH_TI                           ! TI of the horizontal wind speed at the hub point
REAL(ReKi)              ::  UHmax                           ! Maximum horizontal wind speed at the hub
REAL(ReKi)              ::  UHmin                           ! Minimum horizontal wind speed at the hub
REAL(ReKi)              ::  Umax                            ! Maximum u-component wind speed at the hub
REAL(ReKi)              ::  Umin                            ! Minimum u-component wind speed at the hub
REAL(ReKi)              ::  USig                            ! Standard deviation of the u-component wind speed at the hub
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
REAL(ReKi)              ::  Uwave                           ! Wind speed at center of the k-h wave 
REAL(ReKi)              ::  UWMax                           ! Maximum u'w' Reynolds Stress at the hub
REAL(ReKi)              ::  UWMin                           ! Minimum u'w' Reynolds Stress at the hub
REAL(ReKi)              ::  UWTmp                           ! The instantaneous u'w' Reynolds stress at the hub
REAL(ReKi)              ::  UW_RS                           ! The average u'w' Reynolds stress at the hub
REAL(ReKi)              ::  UWcor                           ! The u-w cross component correlation coefficient at the hub
REAL(ReKi)              ::  UWsum                           ! The sum of the u'w' Reynolds stress component at the hub
REAL(ReKi)              ::  UXmax                           ! Maximum U-component (X-direction) wind speed at the hub
REAL(ReKi)              ::  UXmin                           ! Minimum U-component wind speed at the hub
REAL(ReKi)              ::  UXSig                           ! Standard deviation of the U-component wind speed at the hub
REAL(ReKi)              ::  UYmax                           ! Maximum V-component (Y-direction) wind speed at the hub
REAL(ReKi)              ::  UYmin                           ! Minimum V-component wind speed at the hub
REAL(ReKi)              ::  UYSig                           ! Standard deviation of the V-component wind speed at the hub
REAL(ReKi)              ::  UZmax                           ! Maximum W-component (Z-direction) wind speed at the hub
REAL(ReKi)              ::  UZmin                           ! Minimum W-component wind speed at the hub
REAL(ReKi)              ::  UZSig                           ! Standard deviation of the W-component wind speed at the hub
REAL(ReKi)              ::  U_TI                            ! The u-component turbulence intensity at the hub
REAL(ReKi)              ::  Vmax                            ! Maximum v-component wind speed at the hub
REAL(ReKi)              ::  Vmin                            ! Minimum v-component wind speed at the hub
REAL(ReKi)              ::  VSig                            ! Standard deviation of the v-component wind speed at the hub
REAL(ReKi)              ::  VWMax                           ! Maximum v'w' Reynolds Stress at the hub
REAL(ReKi)              ::  VWMin                           ! Minimum v'w' Reynolds Stress at the hub
REAL(ReKi)              ::  VWTmp                           ! The instantaneous v'w' Reynolds stress at the hub
REAL(ReKi)              ::  VW_RS                           ! The average v'w' Reynolds stress at the hub
REAL(ReKi)              ::  VWcor                           ! The v-w cross component correlation coefficient at the hub
REAL(ReKi)              ::  VWsum                           ! The sum of the v'w' Reynolds stress component at the hub
REAL(ReKi)              ::  V_TI                            ! The v-component turbulence intensity at the hub
REAL(ReKi)              ::  Wmax                            ! Maximum w-component wind speed at the hub
REAL(ReKi)              ::  Wmin                            ! Minimum w-component wind speed at the hub
REAL(ReKi)              ::  WSig                            ! Standard deviation of the w-component wind speed at the hub
REAL(ReKi)              ::  W_TI                            ! The w-component turbulence intensity at the hub
REAL(ReKi)              ::  Zbottom                         ! The height of the lowest point on the grid (before tower points are added), equal to Z(1)
REAL(ReKi)              ::  TmpReal                         ! A temporary variable holding a Z (height) value, and other misc. variables
REAL(ReKi)              ::  HorVar                          ! Variables used when DEBUG_OUT is set
REAL(ReKi)              ::  ROT                             ! Variables used when DEBUG_OUT is set
REAL(ReKi)              ::  Total                           ! Variables used when DEBUG_OUT is set
REAL(ReKi)              ::  TotalU                          ! Variables used when DEBUG_OUT is set             
REAL(ReKi)              ::  TotalV                          ! Variables used when DEBUG_OUT is set
REAL(ReKi)              ::  TotalW                          ! Variables used when DEBUG_OUT is set
INTEGER                 ::  AllocStat                       ! The status of allocating space for variables
INTEGER                 ::  HubIndx                         ! Index that tells where the hub point is in the V matrix
INTEGER                 ::  I                               ! A loop counter
INTEGER                 ::  IFreq                           ! Index for frequency
INTEGER                 ::  II                              ! An index for points in the velocity matrix
INTEGER                 ::  Indx                            ! An index for points in the velocity matrix
INTEGER                 ::  IT                              ! Index for time step
INTEGER                 ::  IVec                            ! Index for u-, v-, or w-wind component
INTEGER                 ::  IY                              ! An index for the Y position of a point I
INTEGER                 ::  IZ                              ! An index for the Z position of a point I
INTEGER                 ::  JZ                              ! An index for the Z position of a point J
INTEGER                 ::  NSize                           ! Size of the spectral matrix at each frequency

INTEGER                 ::  NumGrid_Y2                      ! Y Index of the hub (or the nearest point left the hub if hub does not fall on the grid)
INTEGER                 ::  NumGrid_Z2                      ! Z Index of the hub (or the nearest point below the hub if hub does not fall on the grid) 
INTEGER                 ::  NumSteps4                       ! one-fourth the number of steps
INTEGER                 ::  ZHi_YHi                         ! Index for interpolation of hub point, if necessary
INTEGER                 ::  ZHi_YLo                         ! Index for interpolation of hub point, if necessary
INTEGER                 ::  ZLo_YHi                         ! Index for interpolation of hub point, if necessary
INTEGER                 ::  ZLo_YLo                         ! Index for interpolation of hub point, if necessary

LOGICAL                 ::  HubPr                           ! Flag to indicate if the hub height is to be printed separately in the summary file

CHARACTER(1)            ::  Comp (3) = (/ 'u', 'v', 'w' /)  ! The names of the wind components

REAL(ReKi), EXTERNAL    ::  FindZ0                          ! An external function for the modified von Karman spectra

!BONNIE:*****************************
!    Time = TIMEF() ! Initialize the Wall Clock Time counter
!BONNIE:*****************************   

!Beep = .FALSE.
!print *, B1Ki, B2Ki, B4Ki, B8Ki
!print *, SiKi, DbKi, QuKi, ReKi 

   ! Initialize the Pi constants from NWTC_Library
   
CALL NWTC_Init

   ! Set the version number.

CALL SetVersion


   ! Open the console for output.

CALL OpenCon


   ! Print out program name, version, and date.

CALL DispNVD


   ! Check for command line arguments.

CALL CheckArgs( InFile )


   ! Open input file and summary file.

CALL GetFiles


   ! Get input parameters.

CALL GetInput

IF ( NumUSRz > 0 ) THEN
   FormStr = "( // 'User-Defined profiles:' / )"
   WRITE (US,FormStr)
   
   IF ( ALLOCATED( L_USR ) ) THEN
      FormStr = "(A97)"
      WRITE (US,FormStr) '  Height   Wind Speed   Horizontal Angle   u Std. Dev.   v Std. Dev.   w Std. Dev.   Length Scale'
      WRITE (US,FormStr) '   (m)        (m/s)          (deg)            (m/s)         (m/s)         (m/s)          (m)    '
      WRITE (US,FormStr) '  ------   ----------   ----------------   -----------   -----------   -----------   ------------'
      
      FormStr = "( 1X,F7.2, 2X,F9.2,2X, 3X,F10.2,6X, 3(4X,F7.2,3X), 3X,F10.2 )"
      DO I=NumUSRz,1,-1
         WRITE (US,FormStr)  Z_USR(I), U_USR(I), WindDir_USR(I), Sigma_USR(I)*StdScale(1), Sigma_USR(I)*StdScale(2), &
                             Sigma_USR(I)*StdScale(3), L_USR(I)      
      ENDDO   
   ELSE
      FormStr = "(A40)"
      WRITE (US,FormStr) '  Height   Wind Speed   Horizontal Angle'
      WRITE (US,FormStr) '   (m)        (m/s)          (deg)      '
      WRITE (US,FormStr) '  ------   ----------   ----------------'
     
      FormStr = "( 1X,F7.2, 2X,F9.2,2X, 3X,F10.2,6X)"
      DO I=NumUSRz,1,-1
         WRITE (US,FormStr) Z_USR(I), U_USR(I), WindDir_USR(I)      
      ENDDO   
   ENDIF
   
ENDIF


   ! Determine if coherent turbulence should be generated

!BONNIE: UPPER LIMIT ON RICH_NO?
IF ( WrACT ) THEN
   IF ( TurbModel(1:3) == 'IEC' .OR. TurbModel == 'MODVKM' .OR. TurbModel == 'USRINP' .OR. TurbModel(1:5) == 'TIDAL' ) THEN
      FormStr = "( // 'Coherent turbulence time step files are not available for IEC or TIDAL spectral models.')"

      WRITE (US, FormStr)
      CALL TS_Warn( ' A coherent turbulence time step file cannot be generated with the '//TRIM(TurbModel)//' model.', .TRUE. )
      WrACT = .FALSE.
      
   ELSEIF ( Rich_No <= -0.05 ) THEN
      FormStr = "( // 'Coherent turbulence time step files are not available when "// &
                      "the Richardson number is less than or equal to -0.05.')"
      
      WRITE (US, FormStr)
      CALL TS_Warn( ' A coherent turbulence time step file cannot be generated for RICH_NO <= -0.05.', .TRUE. )
      WrACT = .FALSE.
      
   ELSEIF ( .NOT. ( WrADFF .OR. WrBLFF ) ) THEN
      WrBLFF = .TRUE.
      CALL WrScr1( ' AeroDyn/BLADED Full-Field files will be generated along with the coherent turbulence file.' )
!BJJ fix in next release: We can't default to this file type at this point because AeroDyn can't read the .bts files yet.      
      !WrADFF = .TRUE.
      !CALL WrScr1( ' AeroDyn Full-Field files will be generated along with the coherent turbulence file.' )
   ENDIF
      
ENDIF !WrAct

      ! Warn if EWM is used with incompatible times
      
IF ( ( IEC_WindType == IEC_EWM1 .OR. IEC_WindType == IEC_EWM50 ) .AND. & 
      ABS( REAL(600.0,ReKi) - MAX(AnalysisTime,UsableTime) ) > REAL(90.0, ReKi) ) THEN
   CALL TS_Warn( ' The EWM parameters are valid for 10-min simulations only.', .TRUE. )
ENDIF        


   ! Open appropriate output files.  We will open formatted FF files later, if requested.
   ! Mention the files in the summary file.

IF ( WrBHHTP .OR. WrFHHTP .OR. WrADHH .OR. WrADFF .OR. WrFmtFF .OR. WrADTWR .OR. WrACT .OR. WrBLFF )  THEN

   FormStr = "( // 'You have requested that the following file(s) be generated:' / )"
   WRITE (US,FormStr)
   CALL WrScr1 ( ' You have requested that the following file(s) be generated:' )

   IF ( WrBHHTP )  THEN   

      CALL OpenUOutfile ( UGTP , TRIM( RootName)//'.bin' )
      FormStr = "( 3X , A , ' (hub-height binary turbulence-parameter file)' )"
      WRITE (US,FormStr)  TRIM( RootName)//'.bin'
      CALL WrScr ( '    '//TRIM( RootName)//'.bin (a binary hub-height turbulence-parameter file)' )

   ENDIF

   IF ( WrFHHTP )  THEN     

      CALL OpenFOutFile ( UFTP, TRIM( RootName)//'.dat' )
      FormStr = "( 3X , A , ' (formatted turbulence-parameter file)' )"
      WRITE (US,FormStr)  TRIM( RootName)//'.dat'
      CALL WrScr ( '    '//TRIM( RootName)//'.dat (a formatted turbulence-parameter file)' )

   ENDIF

   IF ( WrADHH )  THEN     

      CALL OpenFOutFile ( UAHH, TRIM( RootName)//'.hh' )
      FormStr = "( 3X , A , '  (AeroDyn hub-height file)' )"
      WRITE (US,FormStr)  TRIM( RootName)//'.hh'
      CALL WrScr ( '    '//TRIM( RootName)//'.hh  (an AeroDyn hub-height file)' )

   ENDIF

   IF ( WrADFF )  THEN

      CALL OpenBin ( UAFFW, TRIM(RootName)//'.bts', 2 )

      FormStr = "( 3X , A , ' (AeroDyn/TurbSim full-field wind file)' )"
      WRITE (US,FormStr)  TRIM( RootName)//'.bts'
      CALL WrScr ( '    '//TRIM( RootName)//'.bts (an AeroDyn/TurbSim full-field wind file)' )

   ENDIF

   IF ( WrBLFF )  THEN

      CALL OpenBin ( UBFFW, TRIM(RootName)//'.wnd', 2 )

      FormStr = "( 3X , A , ' (AeroDyn/BLADED full-field wnd file)' )"
      WRITE (US,FormStr)  TRIM( RootName)//'.wnd'
      CALL WrScr ( '    '//TRIM( RootName)//'.wnd (an AeroDyn/BLADED full-field wnd file)' )

   ENDIF

   IF ( WrADTWR .AND. (WrBLFF .OR. .NOT. WrADFF) )  THEN

      CALL OpenBin ( UATWR, TRIM( RootName )//'.twr', 2 )

      FormStr = "( 3X , A , ' (Binary tower twr file)' )"
      WRITE (US,FormStr)  TRIM( RootName)//'.twr'
      CALL WrScr ( '    '//TRIM( RootName)//'.twr (a binary tower file)' )

   ENDIF

   IF ( WrACT ) THEN

      CALL OpenFInpFile ( UACT,  CTEventFile  )
      CALL OpenFOutFile ( UACTTS, TRIM( RootName )//'.cts' )         
      
      FormStr = "( 3X , A , ' (AeroDyn coherent turbulence time step file)' )"
      WRITE (US,FormStr)  TRIM( RootName)//'.cts'
      CALL WrScr ( '    '//TRIM( RootName)//'.cts (a coherent turbulence time step file)' )

   ENDIF

   IF ( WrFmtFF )  THEN
      FormStr = "( 3X , A , ' (formatted full-field U-component file)' )"
      WRITE (US,FormStr)  TRIM( RootName)//'.u'
      CALL WrScr ( '    '//TRIM( RootName)//'.u (a formatted full-field U-component file)' )

      FormStr = "( 3X , A , ' (formatted full-field V-component file)' )"
      WRITE (US,FormStr)  TRIM( RootName)//'.v'
      CALL WrScr ( '    '//TRIM( RootName)//'.v (a formatted full-field V-component file)' )

      FormStr = "( 3X , A , ' (formatted full-field W-component file)' )"
      WRITE (US,FormStr)  TRIM( RootName)//'.w'
      CALL WrScr ( '    '//TRIM( RootName)//'.w (a formatted full-field W-component file)' )
   ENDIF

ELSE

   CALL TS_Abort ( 'You have requested no output.' )

ENDIF


   ! Calculate Total time and NumSteps.
   ! Find the product of small factors that is larger than NumSteps (prime #9 = 23).
   ! Make sure it is a multiple of 4 too.

NumOutSteps = CEILING( ( UsableTime + GridWidth/UHub )/TimeStep )
NumSteps    = MAX( CEILING( AnalysisTime / TimeStep ), NumOutSteps )
NumSteps4   = ( NumSteps - 1 )/4 + 1
NumSteps    = 4*PSF( NumSteps4 , 9 )  ! >= 4*NumSteps4 = NumOutSteps + 3 - MOD(NumOutSteps-1,4) >= NumOutSteps
INumSteps   = 1.0/NumSteps

   ! If we are going to generate debug output, open the debug file.

IF (DEBUG_OUT) THEN
   CALL OpenFOutFile ( UD, TRIM( RootName)//'.dbg' )
   
   WRITE (UD,*)  'UHub=',UHub
   WRITE (UD,*)  'NumOutSteps, NumSteps, INumSteps=', NumOutSteps, NumSteps, INumSteps
ENDIF


IF ( NumSteps < 8 )  THEN
   CALL TS_Abort ( 'Too few (less than 8) time steps.  Increase the usable length of the time series or decrease the time step.' )
ENDIF

   ! Define the other parameters for the time series.

GridRes_Y = GridWidth  / REAL( NumGrid_Y - 1, ReKi )
GridRes_Z = GridHeight / REAL( NumGrid_Z - 1, ReKi )      

Zbottom = HubHt + 0.5*RotorDiameter                         ! height of the highest grid points
Zbottom = Zbottom - GridRes_Z * REAL(NumGrid_Z - 1, ReKi)   ! height of the lowest grid points

IF ( Zbottom <= 0.0 ) THEN
   CALL TS_Abort ( 'The lowest grid point ('//TRIM(Flt2LStr(Zbottom))// ' m) must be above the ground. '//&
                   'Adjust the appropriate values in the input file.' )
ENDIF

NumGrid_Y2 = INT( ( NumGrid_Y + 1 ) / 2 )                    ! These are the hub indicies, unless the hub is an extra point
NumGrid_Z2 = INT( Tolerance + ( HubHt - Zbottom ) / GridRes_Z ) + 1 

ExtraTwrPt = .FALSE.

IF ( MOD(NumGrid_Y, 2) == 0 ) THEN
   ExtraTwrPt = .TRUE.
   ExtraHubPt = .TRUE.
ELSEIF ( ABS((NumGrid_Z2-1)*GridRes_Z + Zbottom - HubHt) > Tolerance ) THEN
   ExtraHubPt = .TRUE.
ELSE
   ExtraHubPt = .FALSE.
ENDIF

NumFreq = NumSteps/2
DelF    = 1.0/( NumSteps*TimeStep )
DelF5   = 0.5*DelF
NTot    = NumGrid_Y*NumGrid_Z                ! Number of points in the regular grid

IF (ExtraHubPT) THEN
   NTot = NTot + 1                           ! Add the hub point if necessary
   ZLim = NumGrid_Z+1
   YLim = NumGrid_Y+1
ELSE
   ZLim = NumGrid_Z
   YLim = NumGrid_Y
ENDIF

IF ( WrADTWR ) THEN

      ! Compute the number of points between the bottom of the grid and the ground 
      ! ( but we don't want to be on the ground, just more than "Tolerance" from it )
 
   IZ = INT( ( Zbottom - Tolerance ) / GridRes_Z )

   IF ( ExtraTwrPt ) THEN 
      IZ = IZ + 1  ! Let's add the point on the bottom of the grid so tower interpolation is easier in AeroDyn
   ENDIF

   IF ( IZ > 0 ) THEN

      ZLim = ZLim + IZ
      NTot = NTot + IZ                       ! Add the number of tower points

      IF ( .NOT. ExtraHubPT ) THEN
         YLim = YLim + 1
      ENDIF

   ELSE

      CALL TS_Warn ( ' There are no extra tower data points below the grid. Tower output will be turned off.', .TRUE.)  

      WrADTWR = .FALSE.
      CLOSE( UATWR )

   ENDIF

ENDIF

NSize   = NTot*( NTot + 1 )/2

IF (DEBUG_OUT) THEN
   ROT     = 1.0/( NumGrid_Y*TimeStep ) 
   WRITE (UD,*)  'TimeStep,DelF,NumFreq,ROT=',TimeStep,DelF,NumFreq,ROT
   WRITE (UD,*)  'PI,Deg2Rad=',Pi, D2R
ENDIF


   !  Allocate the array of frequencies.

ALLOCATE ( Freq(NumFreq) , STAT=AllocStat )

IF ( AllocStat /= 0 )  THEN
   CALL TS_Abort ( 'Error allocating '//TRIM( Int2LStr( ReKi*NumFreq/1024**2 ) )//' MB for the frequency array.' )
ENDIF


   ! Initialize the arrays.

DO IFreq=1,NumFreq
   Freq(IFreq) = IFreq*DelF
ENDDO 

IF (DEBUG_OUT) THEN
   WRITE (UD,*)  'DelF=',DelF
   WRITE (UD,*)  'Freq(1) = ',Freq(1),'  Freq(NumFreq) = ',Freq(NumFreq)
   WRITE (UD,*)  'GridRes=', GridRes_Y, ' x ', GridRes_Z
   WRITE (UD,*)  ' '
   WRITE (UD,*)  'U='
ENDIF


   !  Allocate the array of the steady, u-component winds.

ALLOCATE ( U(ZLim) , STAT=AllocStat )

IF ( AllocStat /= 0 )  THEN
   CALL TS_Abort ( 'Error allocating memory for the array of the steady, u-component winds.' )
ENDIF


   !  Allocate the array of lateral locations of the grid points.

ALLOCATE ( Y(YLim) , STAT=AllocStat )

IF ( AllocStat /= 0 )  THEN
   CALL TS_Abort ( 'Error allocating memory for the lateral locations of the grid points.' )
ENDIF


   !  Allocate the array of vertical locations of the grid points.

ALLOCATE ( Z(ZLim) , STAT=AllocStat )

IF ( AllocStat /= 0 )  THEN
   CALL TS_Abort ( 'Error allocating memory for the vertical locations of the grid points.' )
ENDIF

IF ( TurbModel == 'GP_LLJ') THEN

      ! Allocate the array for the z/l profile
      
   ALLOCATE ( ZL_profile(ZLim) , STAT=AllocStat )

   IF ( AllocStat /= 0 )  THEN
      CALL TS_Abort ( 'Error allocating memory for the z/l profile.' )
   ENDIF
      
         ! Allocate the array for the ustar profile      
   
   ALLOCATE ( Ustar_profile(ZLim) , STAT=AllocStat )

   IF ( AllocStat /= 0 )  THEN
      CALL TS_Abort ( 'Error allocating memory for the friction velocity profile.' )
   ENDIF
     
ENDIF

IF ( INDEX( 'JU', WindProfileType(1:1) ) > 0 ) THEN
   
         ! Allocate the array for the wind direction profile      
     
   ALLOCATE ( WindDir_profile(ZLim) , STAT=AllocStat )

   IF ( AllocStat /= 0 )  THEN
      CALL TS_Abort ( 'Error allocating memory for the wind direction profile.' )
   ENDIF
   
ENDIF   

   !  Allocate the array of vertical locations of the grid points.

ALLOCATE ( IYmax(ZLim) , STAT=AllocStat )

IF ( AllocStat /= 0 )  THEN
   CALL TS_Abort ( 'Error allocating memory for the number of horizontal locations at each height.' )
ENDIF


   ! Initialize cartesian Y,Z values of the grid.

DO IY = 1,NumGrid_Y
   Y(IY)    = -0.5*GridWidth  + GridRes_Y*( IY - 1 )
ENDDO

DO IZ = 1,NumGrid_Z
   Z(IZ)     = Zbottom + GridRes_Z*( IZ - 1 )
   IYmax(IZ) = NumGrid_Y           ! Number of lateral points at this height
ENDDO

IF (ExtraHubPT) THEN

   Y(NumGrid_Y+1)     = 0.0
   Z(NumGrid_Z+1)     = HubHt
   IYmax(NumGrid_Z+1) = 1

   HubIndx = NumGrid_Y*NumGrid_Z + 1

   JZ = NumGrid_Z + 2              ! The start of tower points, if they exist

ELSE

   HubIndx = NumGrid_Y*(NumGrid_Z2-1) + NumGrid_Y2

   JZ = NumGrid_Z + 1              ! The start of tower points, if they exist

ENDIF

IF ( WrADTWR ) THEN

   Y(YLim) = 0.0
   
   TmpReal  = Z(1)

   IF ( ExtraTwrPt ) THEN 
      TmpReal = TmpReal + GridRes_Z
   ENDIF

   DO IZ = JZ,ZLim
      Z(IZ)     = TmpReal - GridRes_Z
      TmpReal   = Z(IZ)

      IYmax(IZ) = 1                 ! The number of lateral points at this height
   ENDDO

ENDIF


IF ( INDEX( 'JU', WindProfileType(1:1) ) > 0 ) THEN
   U(1:ZLim) = getWindSpeed( UHub, HubHt, Z(1:ZLim), RotorDiameter, PROFILE=WindProfileType, UHANGLE=WindDir_profile) 
ELSE 
   U(1:ZLim) = getWindSpeed( UHub, HubHt, Z(1:ZLim), RotorDiameter, PROFILE=WindProfileType)
ENDIF
IF ( ALLOCATED( ZL_profile ) )        ZL_profile(:) = getZLARY(    U(1:ZLim), Z(1:ZLim) )
IF ( ALLOCATED( Ustar_profile ) )  Ustar_profile(:) = getUstarARY( U(1:ZLim), Z(1:ZLim) )


IF ( WrACT ) THEN
   ScaleWid = RotorDiameter * DistScl           !  This is the scaled height of the coherent event data set
   Zbottom  = HubHt - CTLz*ScaleWid             !  This is the height of the bottom of the wave in the scaled/shifted coherent event data set

   IF ( KHtest ) THEN      
         ! for LES test case....
      ScaleVel = ScaleWid * KHT_LES_dT /  KHT_LES_Zm    
      ScaleVel = 50 * ScaleVel                  ! We want 25 hz bandwidth so multiply by 50
   ELSE
   !   TmpPLExp = PLExp 
   !   PLExp    = MIN( 2.0, 1.35*PLExp )        ! Increase the shear of the background (?)

      ScaleVel =            getWindSpeed(UHub,HubHt,Zbottom+ScaleWid,RotorDiameter,PROFILE=WindProfileType)   ! Velocity at the top of the wave
      ScaleVel = ScaleVel - getWindSpeed(UHub,HubHt,Zbottom,         RotorDiameter,PROFILE=WindProfileType)   ! Shear across the wave
      ScaleVel = 0.5 * ScaleVel                                                                               ! U0 is half the difference between the top and bottom of the billow
      
   !   PLExp = TmpPLExp
   ENDIF

   Uwave = getWindSpeed(UHub,HubHt,Zbottom+0.5*ScaleWid,RotorDiameter,PROFILE=WindProfileType)                 ! WindSpeed at center of wave

!BONNIE: MAYBE WE SHOULDN'T OPEN THIS FILE UNTIL WE NEED TO WRITE TO IT
   IF (ScaleVel < 0. ) THEN
      CALL TS_Warn( ' A coherent turbulence time step file cannot be generated with negative shear.', .TRUE. )
      WrACT = .FALSE.
   ENDIF
ENDIF


IF (DEBUG_OUT) THEN
   DO IZ = 1,ZLim
      WRITE (UD,*)  Z(IZ), U(IZ)
   ENDDO
ENDIF

FormStr = "( // 'Turbulence Simulation Scaling Parameter Summary:' / )"
WRITE (US,FormStr)
FormStr = "('   Turbulence model used                            =  ' , A )"
WRITE (US,FormStr)  TRIM(TMName)

FormStr  = "('   ',A,' =' ,F9.3,A)"
FormStr1 = "('   ',A,' =' ,I9  ,A)"
FormStr2 = "('   ',A,' =  ',A)"

IF ( ( TurbModel  == 'IECKAI' ) .OR. ( TurbModel  == 'IECVKM' ) .OR. ( TurbModel  == 'MODVKM' ) )  THEN

      ! If IECKAI or IECVKM spectral models are specified, determine turb intensity 
      ! and slope of Sigma wrt wind speed from IEC turbulence characteristic, 
      ! IECTurbC = A, B, or C or from user specified quantity.
      

   IF ( NumTurbInp )  THEN
   
      TurbInt  = 0.01*PerTurbInt
      SigmaIEC = TurbInt*UHub

   ELSE

      SELECT CASE (IECedition)

         CASE ( 2 )

            IF ( IECTurbC  == 'A' ) THEN
               TurbInt15  = 0.18
               SigmaSlope = 2.0
            ELSEIF ( IECTurbC  == 'B' ) THEN
               TurbInt15  = 0.16
               SigmaSlope = 3.0
            ELSE   ! We should never get here, but just to be complete...
               CALL TS_Abort( ' Invalid IEC turbulence characteristic.' )
            ENDIF

            SigmaIEC = TurbInt15*( ( 15.0 + SigmaSlope*UHub ) / ( SigmaSlope + 1.0 ) )
            TurbInt  = SigmaIEC/UHub
         
         CASE ( 3 )

            IF ( IECTurbC == 'A' ) THEN
               TurbInt15  = 0.16
            ELSEIF ( IECTurbC == 'B' ) THEN
               TurbInt15  = 0.14
            ELSEIF ( IECTurbC == 'C' ) THEN
               TurbInt15  = 0.12
            ELSE   ! We should never get here, but just to be complete...
               CALL TS_Abort( ' Invalid IEC turbulence characteristic.' )
            ENDIF  

                   
            SELECT CASE ( IEC_WindType )
               CASE ( IEC_NTM )
                  SigmaIEC = TurbInt15*( 0.75*UHub + 5.6 )                                      ! [IEC-1 Ed3 6.3.1.3 (11)]
               CASE ( IEC_ETM )
                  Vave     = 0.2*Vref                                                           ! [IEC-1 Ed3 6.3.1.1 ( 9)]
                  SigmaIEC = ETMc*TurbInt15*( 0.072*(Vave/ETMc + 3.0)*(Uhub/ETMc - 4.0)+10.0 )  ! [IEC-1 Ed3 6.3.2.3 (19)]
               CASE ( IEC_EWM1, IEC_EWM50 )
                  Vave     = 0.2*Vref                                                           ! [IEC-1 Ed3 6.3.1.1 ( 9)]
                  SigmaIEC = 0.11*Uhub                                                          ! [IEC-1 Ed3 6.3.2.1 (16)]
               CASE DEFAULT 
                  CALL TS_Abort( 'Invalid IEC wind type.')
            END SELECT           
            TurbInt  = SigmaIEC/UHub     
            
         CASE DEFAULT ! Likewise, this should never happen...

            CALL TS_Abort( 'Invalid IEC 61400-1 edition number.' )
            
         END SELECT                            
      

   ENDIF

   IF (DEBUG_OUT) THEN
      WRITE (UD,*)  'SigmaIEC,TurbInt=',SigmaIEC,TurbInt
   ENDIF

      ! IEC turbulence scale parameter, Lambda1, and IEC coherency scale parameter, LC

   IF ( IECedition == 2 ) THEN
      IF ( HubHt < 30.0 )  THEN
         Lambda1 = 0.7*HubHt
      ELSE
         Lambda1 = 21.0
      ENDIF

      LC = 3.5*Lambda1

   ELSE !IF (IECedition == 3 )
      IF ( HubHt < 60.0 )  THEN
         Lambda1 = 0.7*HubHt
      ELSE
         Lambda1 = 42.0
      ENDIF

      LC = 8.1*Lambda1

   ENDIF

   IF ( MVK .AND. TurbModel  == 'MODVKM' ) THEN
      z0 = FindZ0(HubHt, SigmaIEC, UHub, Fc)
      CALL ScaleMODVKM(HubHt, UHub, TmpU, TmpV, TmpW)
   ENDIF


      ! Write out a parameter summary to the summary file.

   IF ( NumTurbInp ) THEN
      WRITE (US,FormStr2)      "Turbulence characteristic                       ", "User-specified"
   ELSE
      WRITE (US,FormStr2)      "Turbulence characteristic                       ", TRIM(IECTurbE)//IECTurbC
      WRITE (US,FormStr2)      "IEC turbulence type                             ", TRIM(IEC_WindDesc)
      
      IF ( IEC_WindType /= IEC_NTM ) THEN       
         WRITE (US,FormStr)    "Reference wind speed average over 10 minutes    ", Vref,                      " m/s"
         WRITE (US,FormStr)    "Annual wind speed average at hub height         ", Vave,                      " m/s"
      ENDIF
   ENDIF      
   
   WRITE (US,FormStr2)         "IEC standard                                    ", IECeditionSTR(IECedition)
   
   IF ( TurbModel  /= 'MODVKM' ) THEN
      ! Write out a parameter summary to the summary file.

      WRITE (US,FormStr)       "Mean wind speed at hub height                   ", UHub,                      " m/s"

      IF (.NOT. NumTurbInp) THEN
         IF ( IECedition == 2 ) THEN
            WRITE (US,FormStr) "Char value of turbulence intensity at 15 m/s    ", 100.0*TurbInt15,           "%"
            WRITE (US,FormStr) "Standard deviation slope                        ", SigmaSlope,                ""
         ELSE                                                                                                 
               ! This is supposed to be the expected value of what is measured at a site.                     
               ! We actually calculate the 90th percentile value to use in the code as the                    
               ! "Characteristic Value".                                                                      
            WRITE (US,FormStr) "Expected value of turbulence intensity at 15 m/s", 100.0*TurbInt15,           "%"
         ENDIF                                                                                                
                                                                                                              
      ENDIF                                                                                                   
                                                                                                              
      WRITE (US,FormStr)       "Characteristic value of standard deviation      ", SigmaIEC,                  " m/s"
      WRITE (US,FormStr)       "Turbulence scale                                ", Lambda1,                   " m"
                                                                                                              
      IF ( TurbModel  == 'IECKAI' )  THEN                                                                     
         WRITE (US,FormStr)    "u-component integral scale                      ", Lambda1*8.1,               " m"
         WRITE (US,FormStr)    "Coherency scale                                 ", LC,                        " m"
      ELSEIF ( TurbModel  == 'IECVKM' )  THEN                                                                 
         WRITE (US,FormStr)    "Isotropic integral scale                        ", LC,                        " m"
      ENDIF                                                                                                   
                                                                                                              
      WRITE (US,FormStr)       "Characteristic value of hub turbulence intensity", 100.0*TurbInt,             "%"
                                                                                                              
   ELSE                                                                                                       
      WRITE (US,FormStr1)      "Boundary layer depth                            ", NINT(h),                   " m"
      WRITE (US,FormStr)       "Site Latitude                                   ", Latitude,                  " degs"
      WRITE (US,FormStr)       "Hub mean streamwise velocity                    ", UHub,                      " m/s"
      WRITE (US,FormStr)       "Hub local u*                                    ", UStar,                     " m/s" !BONNIE: is this LOCAL? of Disk-avg
      WRITE (US,FormStr)       "Target IEC Turbulence Intensity                 ", 100.0*TurbInt,             "%"
      WRITE (US,FormStr)       "Target IEC u-component standard deviation       ", SigmaIEC,                  " m/s"
      WRITE (US,FormStr)       "u-component integral scale                      ", TmpU,                      " m"
      WRITE (US,FormStr)       "v-component integral scale                      ", TmpV,                      " m"
      WRITE (US,FormStr)       "w-component integral scale                      ", TmpW,                      " m"
      WRITE (US,FormStr)       "Isotropic integral scale                        ", LC,                        " m"
   ENDIF                                                                                                      
   WRITE (US,FormStr)          "Gradient Richardson number                      ", 0.0,                       ""

! Ustar = SigmaIEC/2.15 ! Value based on equating original Kaimal spectrum with IEC formulation

ELSEIF ( TRIM(TurbModel) == 'TIDAL' ) THEN
   WRITE (US,FormStr2)         "Gradient Richardson number                      ", "N/A"
   WRITE (US,FormStr)          "Mean velocity at hub height                     ", UHub,                      " m/s"     
   
ELSE   
   LC = 0.0    ! The length scale is not defined for the non-IEC models
 
   WRITE (US,FormStr)          "Gradient Richardson number                      ", RICH_NO,                   ""
   WRITE (US,FormStr)          "Monin-Obukhov (M-O) z/L parameter               ", ZL,                        ""
                                                                                                              
   IF ( ZL /= 0.0 ) THEN                                                                                      
      WRITE (US,FormStr)       "Monin-Obukhov (M-O) length scale                ", L,                         " m"
   ELSE                                                                                                       
      WRITE (US,FormStr2)      "Monin-Obukhov (M-O) length scale                ", "Infinite"                 
   ENDIF                                                                                                      
   WRITE (US,FormStr)          "Mean wind speed at hub height                   ", UHub,                      " m/s"     
    
ENDIF !  TurbModel  == 'IECKAI', 'IECVKM', or 'MODVKM'

TmpReal = 0.5*RotorDiameter
WTmp    = getWindSpeed(UHub,HubHt,HubHt+TmpReal,RotorDiameter,PROFILE=WindProfileType)   !Velocity at the top of rotor
VTmp    = getWindSpeed(UHub,HubHt,HubHt-TmpReal,RotorDiameter,PROFILE=WindProfileType)   !Velocity at the bottom of the rotor
      
WRITE(US,'()')   ! A BLANK LINE

SELECT CASE ( TRIM(WindProfileType) )
   CASE ('JET','J')
      PLExp = LOG( WTmp/VTmp ) / LOG( (HubHt+TmpReal)/(HubHt-TmpReal) )  !TmpReal = RotorDiameter/2
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
      PLExp = LOG( WTmp/VTmp ) / LOG( (HubHt+TmpReal)/(HubHt-TmpReal) )  !TmpReal = RotorDiameter/2
      
      WRITE (US,FormStr2)      "Wind profile type                               ", "Logarithmic"      
      WRITE (US,FormStr)       "Equivalent power law exponent across rotor disk ",  PLExp,                    ""

   CASE ('H2L','H')
      PLExp = LOG( WTmp/VTmp ) / LOG( (HubHt+TmpReal)/(HubHt-TmpReal) )  !TmpReal = RotorDiameter/2
      
      WRITE (US,FormStr2)      "Velocity profile type                           ", "Logarithmic (H2L)"      
      WRITE (US,FormStr)       "Equivalent power law exponent across rotor disk ",  PLExp,                    ""

   CASE ('PL','P')
      WRITE (US,FormStr2)      "Wind profile type                               ", "Power law"      
      WRITE (US,FormStr)       "Power law exponent                              ",  PLExp,                    ""
      
   CASE ('USR','U')
      PLExp = LOG( WTmp/VTmp ) / LOG( (HubHt+TmpReal)/(HubHt-TmpReal) )  !TmpReal = RotorDiameter/2    

      WRITE (US,FormStr2)      "Wind profile type                               ", "Linear interpolation of user-defined profile"
      WRITE (US,FormStr)       "Equivalent power law exponent across rotor disk ",  PLExp,                    ""
                               
   CASE DEFAULT                
      WRITE (US,FormStr2)      "Wind profile type                               ", "Power law on rotor disk, logarithmic elsewhere"
      WRITE (US,FormStr)       "Power law exponent                              ",  PLExp,                    ""
      
END SELECT

WRITE(US,FormStr)              "Mean shear across rotor disk                    ", (WTmp-VTmp)/RotorDiameter, " (m/s)/m"
WRITE(US,FormStr)              "Assumed rotor diameter                          ", RotorDiameter,             " m"      
WRITE(US,FormStr)              "Surface roughness length                        ", z0,                        " m"      
WRITE(US,'()')                                                                                                 ! A BLANK LINE
WRITE(US,FormStr1)             "Number of time steps in the FFT                 ", NumSteps,                  ""       
WRITE(US,FormStr1)             "Number of time steps output                     ", NumOutSteps,               ""          

IF (KHtest) THEN
   WRITE(US,"(/'KH Billow Test Parameters:' / )") ! HEADER
   WRITE(US,FormStr)           "Gradient Richardson number                      ", RICH_NO,                   ""
   WRITE(US,FormStr)           "Power law exponent                              ", PLexp,                     ""
   WRITE(US,FormStr)           "Length of coherent structures                   ", UsableTime / 2.0,          " s"
   WRITE(US,FormStr)           "Minimum coherent TKE                            ", 30.0,                      " (m/s)^2"
ENDIF

!bjj: This doesn't need to be part of the summary file, so put it in the debug file
IF ( DEBUG_OUT .AND. WindProfileType(1:1) == 'J' ) THEN
   FormStr = "(//'Jet wind profile Chebyshev coefficients:' / )"
   WRITE (UD,FormStr) 
   
   FormStr = "( 3X,'Order: ',11(1X,I10) )"
   WRITE(US,FormStr)  ( IY, IY=0,10 )
   FormStr = "( 3X,'------ ',11(1X,'----------') )"
   WRITE(UD,FormStr)
    
   FormStr = "( 3X,'Speed: ',11(1X,E10.3) )"
   WRITE(UD,FormStr)  ( ChebyCoef_WS(IY), IY=1,11 )

   FormStr = "( 3X,'Angle: ',11(1X,E10.3) )"
   WRITE(UD,FormStr)  ( ChebyCoef_WD(IY), IY=1,11 )   
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
HubPr = ( ABS( HubHt - Z(NumGrid_Z2) ) > Tolerance )     !If the hub height is not on the z-grid, print it, too.

   ! Get the angles to rotate the wind components from streamwise orientation to the X-Y-Z grid at the Hub
            
CVFA = COS( VFlowAng*D2R )
SVFA = SIN( VFlowAng*D2R ) 
CHFA = COS( HFlowAng*D2R )
SHFA = SIN( HFlowAng*D2R )

   ! Write out the grid points & the hub

DO IZ = NumGrid_Z,1, -1
   
   IF ( HubPr  .AND. ( Z(IZ) < HubHt ) ) THEN
   
      JZ = NumGrid_Z+1  ! This is the index of the Hub-height parameters if the hub height is not on the grid
      
      IF ( ALLOCATED( WindDir_profile ) ) THEN      
         CHFA = COS( WindDir_profile(JZ)*D2R )
         SHFA = SIN( WindDir_profile(JZ)*D2R )
         
         IF ( ALLOCATED( ZL_profile ) ) THEN
            
            WRITE(US,FormStr)  Z(JZ), U(JZ), WindDir_profile(JZ), U(JZ)*CHFA*CVFA, U(JZ)*SHFA*CVFA, U(JZ)*SVFA, &
                              ZL_profile(JZ), UStar_profile(JZ)
         ELSE
            WRITE(US,FormStr)  Z(JZ), U(JZ), WindDir_profile(JZ), U(JZ)*CHFA*CVFA, U(JZ)*SHFA*CVFA, U(JZ)*SVFA
         ENDIF
      ELSE
         IF ( ALLOCATED( ZL_profile ) ) THEN
            WRITE(US,FormStr)  Z(JZ), U(JZ), HFlowAng, U(JZ)*CHFA*CVFA, U(JZ)*SHFA*CVFA, U(JZ)*SVFA, &
                              ZL_profile(JZ), UStar_profile(JZ)
         ELSE
            WRITE(US,FormStr)  Z(JZ), U(JZ), HFlowAng, U(JZ)*CHFA*CVFA, U(JZ)*SHFA*CVFA, U(JZ)*SVFA
         ENDIF
      ENDIF
   
      HubPr = .FALSE.
   ENDIF
   
   IF ( ALLOCATED( WindDir_profile ) ) THEN
      CHFA = COS( WindDir_profile(IZ)*D2R )
      SHFA = SIN( WindDir_profile(IZ)*D2R )

      IF ( ALLOCATED( ZL_profile ) ) THEN
         WRITE(US,FormStr)  Z(IZ), U(IZ), WindDir_profile(IZ), U(IZ)*CHFA*CVFA, U(IZ)*SHFA*CVFA, U(IZ)*SVFA, &
                            ZL_profile(IZ), UStar_profile(IZ)
      ELSE
         WRITE(US,FormStr)  Z(IZ), U(IZ), WindDir_profile(IZ), U(IZ)*CHFA*CVFA, U(IZ)*SHFA*CVFA, U(IZ)*SVFA
      ENDIF
   ELSE
      IF ( ALLOCATED( ZL_profile ) ) THEN
         WRITE(US,FormStr)  Z(IZ), U(IZ), HFlowAng, U(IZ)*CHFA*CVFA, U(IZ)*SHFA*CVFA, U(IZ)*SVFA, &
                            ZL_profile(IZ), UStar_profile(IZ)
      ELSE
         WRITE(US,FormStr)  Z(IZ), U(IZ), HFlowAng, U(IZ)*CHFA*CVFA, U(IZ)*SHFA*CVFA, U(IZ)*SVFA
      ENDIF
   ENDIF                

ENDDO ! IZ
   
   ! Write out the tower points
   
DO IZ = NumGrid_Z,ZLim

   IF ( Z(IZ) < Z(1) ) THEN
      IF ( ALLOCATED( WindDir_profile ) ) THEN
         CHFA = COS( WindDir_profile(IZ)*D2R )
         SHFA = SIN( WindDir_profile(IZ)*D2R )

         IF ( ALLOCATED( ZL_profile ) ) THEN
            WRITE(US,FormStr)  Z(IZ), U(IZ), WindDir_profile(IZ), U(IZ)*CHFA*CVFA, U(IZ)*SHFA*CVFA, U(IZ)*SVFA, &
                               ZL_profile(IZ), UStar_profile(IZ)
         ELSE
            WRITE(US,FormStr)  Z(IZ), U(IZ), WindDir_profile(IZ), U(IZ)*CHFA*CVFA, U(IZ)*SHFA*CVFA, U(IZ)*SVFA
         ENDIF
      ELSE
         IF ( ALLOCATED( ZL_profile ) ) THEN
            WRITE(US,FormStr)  Z(IZ), U(IZ), HFlowAng, U(IZ)*CHFA*CVFA, U(IZ)*SHFA*CVFA, U(IZ)*SVFA, &
                               ZL_profile(IZ), UStar_profile(IZ)
         ELSE
            WRITE(US,FormStr)  Z(IZ), U(IZ), HFlowAng, U(IZ)*CHFA*CVFA, U(IZ)*SHFA*CVFA, U(IZ)*SVFA
         ENDIF
      ENDIF                
   ENDIF

ENDDO ! IZ


   !  Allocate the turbulence PSD array.

ALLOCATE ( S(NumFreq,NTot,3) , STAT=AllocStat )

IF ( AllocStat /= 0 )  THEN
   CALL TS_Abort ( 'Error allocating '//TRIM( Int2LStr( ReKi*NumFreq*NTot*3/1024**2 ) )//' MB for the turbulence PSD array.' )
ENDIF


   !  Allocate the work array.

ALLOCATE ( Work(NumFreq,3) , STAT=AllocStat ) !NumSteps2

IF ( AllocStat /= 0 )  THEN
   CALL TS_Abort ( 'Error allocating '//TRIM( Int2LStr( ReKi*NumFreq*3/1024**2 ) )//' MB for the work array.' )
ENDIF

IF (PSD_OUT) THEN
   CALL OpenFOutFile ( UP, TRIM( RootName )//'.psd')
   WRITE (UP,"(A)")  'PSDs '
   FormStr = "( A4,'"//TAB//"',A4,"//TRIM( Int2LStr( NumFreq ) )//"('"//TAB//"',G10.4) )"
   WRITE (UP, FormStr)  'Comp','Ht', Freq(:)
   FormStr = "( I4,"//TRIM( Int2LStr( NumFreq+1 ) )//"('"//TAB//"',G10.4) )"
ENDIF


   ! Allocate and initialize the DUDZ array for MHK models (TIDAL and RIVER)

IF ( TurbModel(1:5) == 'TIDAL' .OR. TurbModel(1:5) == 'RIVER' ) THEN
      ! Calculate the shear, DUDZ, for all heights.
   ALLOCATE ( DUDZ(ZLim) , STAT=AllocStat ) ! Shear
   DUDZ(1)=(U(2)-U(1))/(Z(2)-Z(1))
   DUDZ(ZLim)=(U(ZLim)-U(ZLim-1))/(Z(ZLim)-Z(ZLim-1))
   DO I = 2,ZLim-1
      DUDZ(I)=(U(I+1)-U(I-1))/(Z(I+1)-Z(I-1))
   ENDDO
ENDIF

   ! Calculate the single point Power Spectral Densities. 

JZ = 0   ! The index for numbering the points on the grid
DO IZ=1,ZLim

         ! The continuous, one-sided PSDs are evaluated at discrete
         ! points and the results are stored in the "Work" matrix.


!bonnie: fix this so the the IEC runs differently?? It doesn't need as large an array here anymore...
      IF ( TurbModel(1:3) == 'IEC' ) THEN
         CALL PSDcal( HubHt, UHub  )
      ELSEIF ( TurbModel(1:5) == 'TIDAL' .OR. TurbModel(1:5) == 'RIVER' ) THEN ! HydroTurbSim specific.
         !print *, Ustar
         !Sigma_U2=(TurbIntH20*U(IZ))**2 ! A fixed value of the turbulence intensity.  Do we want to implement this?
         Sigma_U2=4.5*Ustar*Ustar*EXP(-2*Z(IZ)/H_ref)
         Sigma_V2=0.5*Sigma_U2
         Sigma_W2=0.2*Sigma_U2
         CALL PSDCal( Z(IZ) , DUDZ(IZ) )
      ELSEIF ( ALLOCATED( ZL_profile ) ) THEN
         CALL PSDcal( Z(IZ), U(IZ), ZL_profile(IZ), Ustar_profile(IZ) )
      ELSE
         CALL PSDcal( Z(IZ), U(IZ) )
      ENDIF               

      IF (DEBUG_OUT) THEN
            TotalU = 0.0
            TotalV = 0.0
            TotalW = 0.0

            DO IFreq=1,NumFreq
               TotalU = TotalU + Work(IFreq,1)
               TotalV = TotalV + Work(IFreq,2)
               TotalW = TotalW + Work(IFreq,3)
            ENDDO ! IFreq

            Total  = SQRT( TotalU*TotalU + TotalV*TotalV )
            TotalU = TotalU*DelF
            TotalV = TotalV*DelF
            TotalW = TotalW*DelF
            HorVar = Total *DelF

            WRITE(UD,*)  'At H=', Z(IZ)
            WRITE(UD,*)  '   TotalU=', TotalU, '  TotalV=', TotalV, '  TotalW=', TotalW, ' (m/s^2)'
            WRITE(UD,*)  '   HorVar=',HorVar,' (m/s^2)', '   HorSigma=',SQRT(HorVar),' (m/s)'
      ENDIF


         ! Discretize the continuous PSD and store it in matrix "S"
           
      DO IVec=1,3
         
         Indx = JZ

         DO IY=1,IYmax(IZ)   
          
            Indx = Indx + 1
 
            DO IFreq=1,NumFreq
               S(IFreq,Indx,IVec) = Work(IFreq,IVec)*DelF5
            ENDDO ! IFreq

         ENDDO !IY

      ENDDO ! IVec

      JZ = JZ + IYmax(IZ)     ! The next starting index at height IZ + 1

ENDDO ! IZ

IF ( PSD_OUT ) THEN
   CLOSE( UP )
ENDIF

   ! Deallocate memory for work array

IF ( ALLOCATED( Work  ) )  DEALLOCATE( Work  )

IF ( ALLOCATED( DUDZ            ) )  DEALLOCATE( DUDZ            )
IF ( ALLOCATED( Z_USR           ) )  DEALLOCATE( Z_USR           )
IF ( ALLOCATED( U_USR           ) )  DEALLOCATE( U_USR           )
IF ( ALLOCATED( WindDir_USR     ) )  DEALLOCATE( WindDir_USR     )
IF ( ALLOCATED( Sigma_USR       ) )  DEALLOCATE( Sigma_USR       )
IF ( ALLOCATED( L_USR           ) )  DEALLOCATE( L_USR           )

IF ( ALLOCATED( ZL_profile      ) )  DEALLOCATE( ZL_profile      )
IF ( ALLOCATED( Ustar_profile   ) )  DEALLOCATE( Ustar_profile   )

IF ( ALLOCATED( Freq_USR   ) )  DEALLOCATE( Freq_USR   )
IF ( ALLOCATED( Uspec_USR  ) )  DEALLOCATE( Uspec_USR  )
IF ( ALLOCATED( Vspec_USR  ) )  DEALLOCATE( Vspec_USR  )
IF ( ALLOCATED( Wspec_USR  ) )  DEALLOCATE( Wspec_USR  )


   ! Allocate memory for random number array

ALLOCATE ( RandNum(NTot*NumFreq*3) , STAT=AllocStat )

IF ( AllocStat /= 0 )  THEN
   CALL TS_Abort ( 'Error allocating '//TRIM( Int2LStr( ReKi*NumFreq*3/1024**2 ) )//' MB for the random-number array.' )
ENDIF

   ! Reinitialize the random number generator ( it was initialized when the
   ! seeds were read in ) so that the same seed always generates the same
   ! random phases, regardless of previous randomizations in this code.

RandSeed(1) = RandSeedTmp
CALL RndInit()

   ! Let's go ahead and get all the random numbers we will need for the entire
   ! run.  This (hopefully) will be faster than getting them one at a time,
   ! but it will use more memory.
   ! These pRNGs have been initialized in the GetInput() subroutine

IF (RNG_type == 'NORMAL') THEN  

      !The first two real numbers in the RandSeed array are used as seeds
      !The number of seeds needed are compiler specific, thus we can't assume only 2 seeds anymore

   CALL RANDOM_NUMBER ( RandNum )

   ! Let's harvest the random seeds so that they can be used for the next run if desired.
   ! Write them to the summary file.

   CALL RANDOM_SEED ( GET=RandSeedAry )

   FormStr = "(//,'Harvested Random Seeds after Generation of the Random Numbers:',/)"
   WRITE(US,FormStr)

   DO I = 1,SIZE( RandSeedAry )
      FormStr = "(I13,' Harvested seed #',I2)"
      WRITE(US,FormStr)  RandSeedAry(I), I
   END DO

   DEALLOCATE(RandSeedAry, STAT=AllocStat)


ELSEIF (RNG_type == 'RANLUX') THEN

   CALL RanLux ( RandNum )

   CALL RLuxAT ( LuxLevel, RandSeed(1), I, II )

   FormStr = "(//,'Harvested Random Seeds after Generation of the Random Numbers:')"
   WRITE(US,FormStr)

   FormStr = "(/,I13,' K1')"
   WRITE(US,FormStr)  I

   FormStr = "(I13,' K2')"
   WRITE(US,FormStr)  II

ELSE
 
   I = NTot*NumFreq

   Indx = 1
   DO IVec = 1,3
      CALL ARand( RandSeed(IVec), RandNum, I,  Indx)
      Indx = Indx + I
   ENDDO

ENDIF


   !  Allocate the transfer-function matrix. 
   
ALLOCATE ( TRH( MAX(NumSteps,NSize) ) , STAT=AllocStat )

IF ( AllocStat /= 0 )  THEN
   CALL TS_Abort ( 'Error allocating '//TRIM( Int2LStr( ReKi*MAX(NumSteps,NSize)/1024**2 ) )// &
                    ' MB for the temporary transfer-function matrix.' )   
ENDIF


   !  Allocate the array that contains the velocities.

ALLOCATE ( V(NumSteps,NTot,3), STAT=AllocStat )

IF ( AllocStat /= 0 )  THEN
   CALL TS_Abort ( 'Error allocating '//TRIM( Int2LStr( ReKi*NumSteps*NTot*3/1024**2 ) )//' MB for velocity array.' )
ENDIF


   ! Calculate the transfer function matrices from the spectral matrix (the fourier coefficients).

CALL WrScr ( ' Calculating the spectral and transfer function matrices:' )

V(:,:,:) = 0.0    ! initialize the velocity matrix

IF (TurbModel(1:4) /= 'NONE') THEN
   CALL CohSpecVMat(LC, NSize, Comp)    
ENDIF


   ! Deallocate the Freq, S, and RandNum arrays and the spectral matrix

IF ( ALLOCATED( Freq    ) )  DEALLOCATE( Freq    )
IF ( ALLOCATED( S       ) )  DEALLOCATE( S       )
IF ( ALLOCATED( RandNum ) )  DEALLOCATE( RandNum )

   !  Allocate the FFT working storage and initialize its variables

CALL InitFFT( NumSteps )


   ! Get the stationary-point time series.

CALL WrScr ( ' Generating time series for all points:' )

DO IVec=1,3

   CALL WrScr ( '    '//Comp(IVec)//'-component' )

   DO Indx=1,NTot    !NTotB

         ! Overwrite the first point with zero.  This sets the real (and 
         ! imaginary) part of the steady-state value to zero so that we 
         ! can add in the mean value later.

      TRH(1)=0.0

      DO IT=2,NumSteps-1
         TRH(IT) = V(IT-1,Indx,IVec)
      ENDDO ! IT

         ! Now, let's add a complex zero to the end to set the power in the Nyquist
         ! frequency to zero.

      TRH(NumSteps) = 0.0


         ! perform FFT

      CALL ApplyFFT( TRH )
        
      V(1:NumSteps,Indx,IVec) = TRH(1:NumSteps)

   ENDDO ! Indx

ENDDO ! IVec

   ! Deallocate the TRH array and the FFT working storage.
IF ( ALLOCATED( TRH  ) )  DEALLOCATE( TRH  )

CALL ExitFFT

   ! Crossfeed cross-axis components to u', v', w' components and scale IEC models if necessary

SELECT CASE ( TRIM(TurbModel) )
      
   CASE ('GP_LLJ','NWTCUP', 'SMOOTH', 'WF_UPW', 'WF_07D', 'WF_14D', 'USRVKM', 'TIDAL', 'RIVER') ! Do reynolds stress for HYDRO also.
               
            ! Calculate coefficients for obtaining "correct" Reynold's stresses at the hub
         UWsum = 0.0
         UVsum = 0.0
         VWsum = 0.0
         USum2 = 0.0
         VSum2 = 0.0
         WSum2 = 0.0
         
         DO IT = 1,NumSteps
            UWsum = UWsum + V(IT,HubIndx,1) * V(IT,HubIndx,3)
            UVsum = UVsum + V(IT,HubIndx,1) * V(IT,HubIndx,2)
            VWsum = VWsum + V(IT,HubIndx,2) * V(IT,HubIndx,3)
            USum2 = USum2 + V(IT,HubIndx,1) * V(IT,HubIndx,1)
            VSum2 = VSum2 + V(IT,HubIndx,2) * V(IT,HubIndx,2)
            WSum2 = WSum2 + V(IT,HubIndx,3) * V(IT,HubIndx,3)
         ENDDO
         UWsum = UWsum * INumSteps  !These "sums" are now "means"
         UVsum = UVsum * INumSteps  !These "sums" are now "means"
         VWsum = VWsum * INumSteps  !These "sums" are now "means"
         USum2 = USum2 * INumSteps  !These "sums" are now "means"
         VSum2 = VSum2 * INumSteps  !These "sums" are now "means"
         WSum2 = WSum2 * INumSteps  !These "sums" are now "means"
            
            !BJJ: this is for v=alpha1, w=alpha2, u=alpha3 using derivation equations
         UWMax = ( PC_UW - UWsum ) / USum2                                             !alpha23
         VWMax = ( USum2*(PC_VW - VWsum - UWMax*UVsum) - PC_UW*(PC_UV - UVsum) ) / &   !alpha12
                 ( USum2*(WSum2 + UWMax*UWsum) - UWsum*PC_UW )
         UVMax = ( PC_UV - UVsum - VWMax*UWsum) / USum2                                !alpha13         
         
                  
            ! if we enter "none" for any of the Reynolds-stress terms, don't scale that component:
         IF (UWskip) UWMax = 0.0
         IF (UVskip) UVMax = 0.0
         IF (VWskip) VWMax = 0.0
         
            !bjj: I'm implementing limits on the range of values here so that the spectra don't get too
            !     out of whack.  We'll display a warning in this case.
            
         IF ( ABS(UWMax) > 1.0 .OR. ABS(UVMax) > 1.0 .OR. ABS(VWMax) > 1.0 ) THEN
            CALL TS_Warn( "Scaling terms exceed 1.0.  Reynolds stresses may be affected.", .FALSE.)
         ENDIF
         
         UWMax = MAX( MIN( UWMax, 1.0 ), -1.0 )
         UVMax = MAX( MIN( UVMax, 1.0 ), -1.0 )
         VWMax = MAX( MIN( VWMax, 1.0 ), -1.0 )
                  
         DO Indx = 1,NTot
            DO IT = 1, NumSteps
               TmpU = V(IT,Indx,1)
               TmpV = V(IT,Indx,2)
               TmpW = V(IT,Indx,3)
                  
                  !BJJ: this is for v=alpha1, w=alpha2, u=alpha3
               V(IT,Indx,2) = UVMax*TmpU + TmpV + VWmax*TmpW 
               V(IT,Indx,3) = UWmax*TmpU +              TmpW
            ENDDO          
         ENDDO

         FormStr = "(//,'Scaling statistics from the hub grid point:',/)"
         WRITE( US, FormStr )
                  
         FormStr = "(3X,'Cross-Component  Scaling Factor')"
         WRITE( US, FormStr )
         FormStr = "(3X,'---------------  --------------')"
         WRITE( US, FormStr )
         FormStr = "(3X,A,2X,E14.5)"
         WRITE( US, FormStr ) "u'w'           ", UWmax
         WRITE( US, FormStr ) "u'v'           ", UVmax
         WRITE( US, FormStr ) "v'w'           ", VWmax

   CASE ( 'IECKAI', 'IECVKM' )

      IF (ScaleIEC > 0) THEN

         FormStr = "(//,'Scaling statistics from the hub grid point:',/)"
         WRITE( US, FormStr )
                  
         FormStr = "(2X,'Component  Target Sigma (m/s)  Simulated Sigma (m/s)  Scaling Factor')"
         WRITE( US, FormStr )
         FormStr = "(2X,'---------  ------------------  ---------------------  --------------')"
         WRITE( US, FormStr )
         FormStr = "(2X,3x,A,7x,f11.3,9x,f12.3,11x,f10.3)"                        

         DO IVec = 1,3
            CGridSum  = 0.0
            CGridSum2 = 0.0
                           
            DO IT=1,NumSteps !BJJ: NumOutSteps  -- scale to the output value?
               CGridSum  = CGridSum  + V( IT, HubIndx, IVec )
               CGridSum2 = CGridSum2 + V( IT, HubIndx, IVec )* V( IT, HubIndx, IVec )
            ENDDO ! IT
               
            UGridMean = CGridSum/NumSteps !BJJ: NumOutSteps  -- scale to the output value?
            UGridSig  = SQRT( ABS( (CGridSum2/NumSteps) - UGridMean*UGridMean ) )
            
            IF ( IVec == 1 .OR. TurbModel == 'IECVKM' ) THEN
               TargetSigma = SigmaIEC
            ELSEIF (IVec == 2) THEN
               TargetSigma = 0.8*SigmaIEC             !IEC 61400-1, Appendix B
            ELSE  
               TargetSigma = 0.5*SigmaIEC             !IEC 61400-1, Appendix B
            ENDIF

            WRITE( US, FormStr) Comp(IVec)//"'", TargetSigma, UGridSig, TargetSigma/UGridSig
            
            IF (ScaleIEC == 1 .OR. IVec > 1) THEN ! v and w have no coherence, thus all points have same std, so we'll save some calculations
               V(:,:,IVec) =     (TargetSigma / UGridSig) * V(:,:,IVec)
            ELSE  ! Scale each point individually
               DO Indx = 1,NTot             
                  CGridSum  = 0.0
                  CGridSum2 = 0.0
                                    
                  DO IT=1,NumSteps !BJJ: NumOutSteps  -- scale to the output value?
                     CGridSum  = CGridSum  + V( IT, Indx, IVec )
                     CGridSum2 = CGridSum2 + V( IT, Indx, IVec )* V( IT, Indx, IVec )
                  ENDDO ! IT
                  
                  UGridMean = CGridSum/NumSteps !BJJ: NumOutSteps  -- scale to the output value?
                  UGridSig  = SQRT( ABS( (CGridSum2/NumSteps) - UGridMean*UGridMean ) )
            
                  V(:,Indx,IVec) = (TargetSigma / UGridSig) * V(:,Indx,IVec)   
               ENDDO ! Indx
            ENDIF            

         ENDDO !IVec
           
      ENDIF !ScaleIEC

END SELECT

      ! Write Headers for Hub-Height files

  ! Generate header for HH AeroDyn file.

IF ( WrADHH )  THEN

   CALL WrScr ( ' Hub-height AeroDyn data were written to "'//TRIM( RootName )//'.hh"' )

   FormStr = "( '! This hub-height wind-speed file was generated by ' , A , A , ' on ' , A , ' at ' , A , '.' )"
   WRITE (UAHH,FormStr)   TRIM(ProgName), TRIM( ProgVer ), CurDate(), CurTime()

   FormStr = "( '!' )"
   WRITE (UAHH,FormStr)

   FormStr = "( '! The requested statistics for this data were:' )"
   WRITE (UAHH,FormStr)

   FormStr = "( '!    Mean Total Wind Speed = ' , F8.3 , ' m/s' )"
   WRITE (UAHH,FormStr)  UHub

   IF ( (TurbModel == 'IECKAI') .OR. (TurbModel == 'IECVKM') .OR. (TurbModel == 'MODVKM') ) THEN
      FormStr = "( '!    Turbulence Intensity  = ' , F8.3 , '%' )"
      WRITE (UAHH,FormStr)  100.0*TurbInt
   ENDIF

   FormStr = "( '!' )"
   WRITE (UAHH,FormStr)

   FormStr = "( '!   Time  HorSpd  WndDir  VerSpd  HorShr  VerShr  LnVShr  GstSpd' )"
   WRITE (UAHH,FormStr)

   FormStr = "( '!  (sec)   (m/s)   (deg)   (m/s)     (-)     (-)     (-)   (m/s)' )"
   WRITE (UAHH,FormStr)

ENDIF ! ( WrADHH )


   ! Generate hub-height formatted turbulence parameters file header.

IF ( WrFHHTP )  THEN

   CALL WrScr ( ' Hub-height formatted turbulence parameters were written to "'//TRIM( RootName )//'.dat"' )

   FormStr = "( / 'This hub-height turbulence-parameter file was generated by ' , A , A , ' on ' , A , ' at ' , A , '.' / )"
   WRITE (UFTP,FormStr)  TRIM(ProgName), TRIM( ProgVer ), CurDate(), CurTime()
 
   FormStr = "('   Time',6X,'U',7X,'Uh',7X,'Ut',8X,'V',8X,'W',8X,'u''',7X,'v''',7X,'w'''," &
           //"6X,'u''w''',5X,'u''v''',5X,'v''w''',5X,'TKE',6X,'CTKE')"
   WRITE (UFTP,FormStr)

ENDIF ! WrFHHTP


   ! Inform user if HH binary turbulence parameters are being written.

IF ( WrBHHTP )  THEN
   CALL WrScr ( ' Hub-height binary turbulence parameters were written to "'//TRIM( RootName )//'.bin"' )
ENDIF


   ! Initialize statistical quantities for hub-height turbulence parameters.

CALL WrScr ( ' Computing hub-height statistics' )

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
UVMax   =  V(1,HubIndx,1)*V(1,HubIndx,2) 
UVMin   =  HUGE( UVMin )
UVsum   =      0.0
UW_RS   =      0.0
UWMax   =  V(1,HubIndx,1)*V(1,HubIndx,3) 
UWMin   =  HUGE( UWMin )
UWsum   =      0.0
VBar    =      0.0
Vmax    = -HUGE( Vmax )
Vmin    =  HUGE( Vmin )
VSum2   =      0.0
VW_RS   =      0.0
VWMax   =  V(1,HubIndx,2)*V(1,HubIndx,3)
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

IF ( ALLOCATED( WindDir_profile ) ) THEN
   IF (ExtraHubPT) THEN    
      JZ = NumGrid_Z+1  ! This is the index of the Hub-height parameters if the hub height is not on the grid
   ELSE
      JZ = NumGrid_Z2
   ENDIF
   CHFA = COS( WindDir_profile(JZ)*D2R )
   SHFA = SIN( WindDir_profile(JZ)*D2R )
ELSE
   CHFA = COS( HFlowAng*D2R )
   SHFA = SIN( HFlowAng*D2R )
ENDIF
         
CVFA = COS( VFlowAng*D2R )
SVFA = SIN( VFlowAng*D2R )

DO IT=1,NumSteps

   Time = TimeStep*( IT - 1 )

      ! Calculate longitudinal (UTmp), lateral (VTmp), and upward (WTmp)
      ! values for hub station, as well as rotated (XTmp, YTmp, ZTmp) 
      ! components applying specified flow angles.

      ! Add mean wind speed to the streamwise component
   UTmp = V(IT,HubIndx,1) + UHub
   VTmp = V(IT,HubIndx,2)
   WTmp = V(IT,HubIndx,3)
   
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

   UVTmp = V(IT,HubIndx,1)*V(IT,HubIndx,2)
   UWTmp = V(IT,HubIndx,1)*V(IT,HubIndx,3)
   VWTmp = V(IT,HubIndx,2)*V(IT,HubIndx,3)

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

   TKE  = 0.5*(V(IT,HubIndx,1)*V(IT,HubIndx,1) + V(IT,HubIndx,2)*V(IT,HubIndx,2) + V(IT,HubIndx,3)*V(IT,HubIndx,3))
   CTKE = 0.5*SQRT(UVTmp*UVTmp + UWTmp*UWTmp + VWTmp*VWTmp)

   IF (CTKE > CTKEmax) CTKEmax = CTKE
   IF ( TKE >  TKEmax)  TKEmax =  TKE

      ! Find sums for mean and square Reynolds stresses for hub-level simulation.
   UVsum = UVsum + UVTmp
   UWsum = UWsum + UWTmp
   VWsum = VWsum + VWTmp

         ! Are we generating HH AeroDyn files?

   IF ( IT <= NumOutSteps ) THEN
      IF ( WrADHH )  THEN
         WRITE (UAHH,'(F8.3,3F8.2,3F8.3,F8.2)')  Time, UHTmp, -1.0*R2D*ATAN2( UYTmp , UXTmp ), &
                                                   UZTmp, 0.0, PLExp, 0.0, 0.0 
!bjj: Should we output instantaneous horizontal shear, instead of 0?  
!     Should the power law exponent be an instantaneous value, too?
!                                                   
      ENDIF ! ( WrADHH )


         ! Output HH binary turbulence parameters for GenPro analysis.
         ! Output order:  Time,U,Uh,Ut,V,W,u',v',w',u'w',u'v',v'w',TKE,CTKE.
   
      IF ( WrBHHTP )  THEN

            !bjj: Vtmp = V(IT,HubIndx,2); WTmp = V(IT,HubIndx,3) so it's redundant to output them twice:

         WRITE (UGTP)  REAL(Time,SiKi), REAL(UXTmp,SiKi), REAL(UHTmp,SiKi), REAL(UTTmp,SiKi), &
                                        REAL(UYTmp,SiKi), REAL(UZTmp,SiKi), &
                       REAL(V(IT,HubIndx,1),SiKi), REAL(V(IT,HubIndx,2),SiKi), REAL(V(IT,HubIndx,3),SiKi), &
                       REAL(UWTmp,SiKi), REAL(UVTmp,SiKi), REAL(VWTmp,SiKi), REAL(TKE,SiKi), REAL(CTKE,SiKi)                        

      ENDIF


         ! Output HH formatted turbulence parameters for plotting and analysis.
         ! Output order:  ET,U,Uh,Ut,V,W,u',v',w',u'w',u'v',v'w', TKE, CTKE.

      IF ( WrFHHTP )  THEN
         WRITE(UFTP,'(F7.2,13F9.3)') Time,UXTmp,UHTmp,UTTmp,UYTmp,UZTmp, &
                                     V(IT,HubIndx,1), V(IT,HubIndx,2), V(IT,HubIndx,3), &
                                     UWTmp, UVTmp, VWTmp, TKE, CTKE
      ENDIF

   ENDIF !IT <= NumOutSteps

ENDDO ! IT

   ! Close hub-height output files
   
IF ( WrADHH  ) CLOSE(UAHH)
IF ( WrFHHTP ) CLOSE(UFTP)
IF ( WrBHHTP ) CLOSE(UGTP)


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

USig  = SQRT( MAX( USum2 *INumSteps-UBar *UBar , DblZero ) )
VSig  = SQRT( MAX( VSum2 *INumSteps-VBar *VBar , DblZero ) )
WSig  = SQRT( MAX( WSum2 *INumSteps-WBar *WBar , DblZero ) )

UXSig = SQRT( MAX( UXSum2*INumSteps-UXBar*UXBar, DblZero ) )
UYSig = SQRT( MAX( UYSum2*INumSteps-UYBar*UYBar, DblZero ) )
UZSig = SQRT( MAX( UZSum2*INumSteps-UZBar*UZBar, DblZero ) )

UHSig = SQRT( MAX( UHSum2*INumSteps-UHBar*UHBar, DblZero ) )
UTSig = SQRT( MAX( UTSum2*INumSteps-UTBar*UTBar, DblZero ) )


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

FormStr = "(/,3X,'                    Min Reynolds     Mean Reynolds    Max Reynolds    Correlation')"
WRITE(US,FormStr)
FormStr = "(  3X,'Product             Stress (m/s)^2   Stress (m/s)^2   Stress (m/s)^2  Coefficient')"
WRITE(US,FormStr)
FormStr = "(  3X,'----------------    --------------   --------------   --------------  -----------')"
WRITE(US,FormStr)

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

ALLOCATE ( SDary(NumGrid_Y) , STAT=AllocStat )

IF ( AllocStat /= 0 )  THEN
   CALL TS_Abort ( 'Error allocating memory for the array of standard deviations.' )
ENDIF


   ! Calculate standard deviations for each grid point.  Write them to summary file.

FormStr = "(//,'Grid Point Variance Summary:',/)"
WRITE(US,FormStr)

FormStr = "(3X,'Y-coord',"//TRIM(FLT2LSTR(REAL(NumGrid_Y,ReKi)))//"F8.2)"

WRITE(US,FormStr)  Y(1:NumGrid_Y)

FormStr1 = "(/,3X'Height   Standard deviation at grid points for the ',A,' component:')"
FormStr2 = "(F9.2,1X,"//TRIM(FLT2LSTR(REAL(NumGrid_Y,ReKi)))//"F8.3)"


UTmp = 0
VTmp = 0
WTmp = 0

DO IVec=1,3

   WRITE(US,FormStr1)  Comp(IVec)

   DO IZ=NumGrid_Z,1,-1

      DO IY=1,NumGrid_Y

         II   = (IZ-1)*NumGrid_Y+IY
         SumS = 0.0

         DO IT=1,NumSteps
            SumS = SumS + V(IT,II,IVec)**2            
         ENDDO ! IT         

         SDary(IY) = SQRT(SumS*INumSteps)   !  Was:  SDary(IZ,IY) = SQRT(SumS*INumSteps)/U(IZ,NumGrid/2)

      ENDDO ! IY

      WRITE(US,FormStr2) Z(IZ), SDary(1:NumGrid_Y)

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
WRITE(US,FormStr2) Comp(1), UTmp / ( NumGrid_Y*NumGrid_Z )
WRITE(US,FormStr2) Comp(2), VTmp / ( NumGrid_Y*NumGrid_Z ) 
WRITE(US,FormStr2) Comp(3), WTmp / ( NumGrid_Y*NumGrid_Z ) 


   !  Deallocate the array of standard deviations.

IF ( ALLOCATED( SDary ) )  DEALLOCATE( SDary )

   ! Add mean wind to u' components.

II = 0
DO IZ=1,ZLim   

   IF ( ALLOCATED( WindDir_profile ) ) THEN  ! The horizontal flow angle changes with height
      CHFA = COS( WindDir_profile(IZ)*D2R )
      SHFA = SIN( WindDir_profile(IZ)*D2R )
   ENDIF      

  DO IY=1,IYmax(IZ)  

      II = II + 1

      DO IT=1,NumSteps

            ! Add mean wind speed to the streamwise component
            
         TmpU = V(IT,II,1) + U(IZ)
         TmpV = V(IT,II,2)
         TmpW = V(IT,II,3)
         
            ! Rotate the wind to the X-Y-Z (inertial) reference frame coordinates
                     
         V(IT,II,1) = TmpU*CHFA*CVFA - TmpV*SHFA - TmpW*CHFA*SVFA
         V(IT,II,2) = TmpU*SHFA*CVFA + TmpV*CHFA - TmpW*SHFA*SVFA  
         V(IT,II,3) = TmpU*SVFA                  + TmpW*CVFA                  

      ENDDO ! IT

  ENDDO ! IY

ENDDO ! IZ

TmpU = MAX( ABS(MAXVAL(U)-UHub), ABS(MINVAL(U)-UHub) )  !Get the range of wind speed values for scaling in BLADED-format .wnd files

   ! Deallocate memory for the matrix of the steady, u-component winds.

IF ( ALLOCATED( U               ) )  DEALLOCATE( U               )
IF ( ALLOCATED( WindDir_profile ) )  DEALLOCATE( WindDir_profile )

   ! Are we generating a coherent turbulent timestep file?
   
IF ( WrACT ) THEN

   CALL WrScr ( ' Generating coherent turbulent time step file "'//TRIM( RootName )//'.cts"' )

   IF ( .NOT. KHtest ) THEN

         ! If the coherent structures do not cover the whole disk, increase the shear

      IF ( DistScl < 1.0 ) THEN ! Increase the shear by up to two when the wave is half the size of the disk...
         CALL RndUnif( TmpReal )
         ScaleVel = ScaleVel * ( 1.0 + TmpReal * (1 - DistScl) / DistScl )
      ENDIF

         !Apply a scaling factor to account for short inter-arrival times getting wiped out due to long events

      ScaleVel =  ScaleVel*( 1.0 + 323.1429 * EXP( -MAX(Uwave,10.0) / 2.16617 ) )

      !TSclFact = ScaleWid / (ScaleVel * Zm_max)

   ENDIF


         ! Determine the maximum predicted CTKE
         
         SELECT CASE ( TRIM(TurbModel) )
         
            CASE ( 'NWTCUP',  'NONE', 'USRVKM' )
            
               IF (KHtest) THEN
                  CTKE = 30.0 !Scale for large coherence
                  CALL RndNWTCpkCTKE( CTKE )
               ELSE    
               
                     ! Increase the Scaling Velocity for computing U,V,W in AeroDyn
                     ! These numbers are based on LIST/ART data (58m-level sonic anemometer)

                  CTKE =  0.616055*Rich_No - 0.242143*Uwave + 23.921801*WSig - 11.082978
            
                     ! Add up to +/- 10% or +/- 6 m^2/s^2 (uniform distribution)
                  CALL RndUnif( TmpReal )
                  CTKE = MAX( CTKE + (2.0 * TmpReal - 1.0) * 6.0, 0.0 )

                  IF ( CTKE > 0.0 ) THEN
                     IF ( CTKE > 20.0)  THEN    ! Correct with residual
                        CTKE = CTKE + ( 0.11749127 * (CTKE**1.369025) - 7.5976449 )
                     ENDIF

                     IF ( CTKE >= 30.0 .AND. Rich_No >= 0.0 .AND. Rich_No <= 0.05 ) THEN
                        CALL RndNWTCpkCTKE( CTKE )
                     ENDIF
                  ENDIF
                  
               ENDIF !KHTest
               
            CASE ( 'GP_LLJ', 'SMOOTH' , 'TIDAL', 'RIVER' )          

               CTKE = pkCTKE_LLJ(Zbottom+0.5*ScaleWid)
               
            CASE ( 'WF_UPW' )
               CTKE = -2.964523*Rich_No - 0.207382*Uwave + 25.640037*WSig - 10.832925
               
            CASE ( 'WF_07D' )
               CTKE = 9.276618*Rich_No + 6.557176*Ustar + 3.779539*WSig - 0.106633

               IF ( (Rich_No > -0.025) .AND. (Rich_No < 0.05) .AND. (UStar > 1.0) .AND. (UStar < 1.56) ) THEN
                  CALL RndpkCTKE_WFTA( TmpReal )  ! Add a random residual
                  CTKE = CTKE + TmpReal
               ENDIF
               
               
            CASE ( 'WF_14D' )
               CTKE = 1.667367*Rich_No - 0.003063*Uwave + 19.653682*WSig - 11.808237
               
            CASE DEFAULT   ! This case should not happen            
               CALL TS_Abort( 'Invalid turbulence model in coherent structure analysis.' )            
               
         END SELECT                    

         CTKE   = MAX( CTKE, 1.0 )     ! make sure CTKE is not negative and, so that we don't divide by zero in ReadEventFile, set it to some arbitrary low number


      ! Read and allocate coherent event start times and lengths, calculate TSclFact

   CALL ReadEventFile( UACT, ScaleWid, ScaleVel, CTKE )

   CLOSE ( UACT )  
   
   CALL CalcEvents( REAL( Uwave, ReKi ), MaxEvtCTKE,  Zbottom+0.5*ScaleWid) 

   WRITE(US, '(A,F8.3," (m/s)^2")') 'Maximum predicted event CTKE         = ', CTKE

      ! Write event data to the time step output file (opened at the beginnig)

   WRITE (UACTTS, "( A14,   ' = FileType')") CTExt
   WRITE (UACTTS, "( G14.7, ' = ScaleVel')") ScaleVel
   WRITE (UACTTS, "( G14.7, ' = MHHWindSpeed')") UHub
   WRITE (UACTTS, "( G14.7, ' = Ymax')") ScaleWid*Ym_max/Zm_max
   WRITE (UACTTS, "( G14.7, ' = Zmax')") ScaleWid
   WRITE (UACTTS, "( G14.7, ' = DistScl')") DistScl
   WRITE (UACTTS, "( G14.7, ' = CTLy')") CTLy
   WRITE (UACTTS, "( G14.7, ' = CTLz')") CTLz
   WRITE (UACTTS, "( G14.7, ' = NumCTt')") NumCTt
   
   CALL WriteEvents ( UACTTS, UACT, TSclFact )

   CLOSE ( UACTTS )

      ! Deallocate the coherent event times arrays.

   IF ( ALLOCATED( EventName ) )  DEALLOCATE( EventName )
   IF ( ALLOCATED( EventTS   ) )  DEALLOCATE( EventTS   )
   IF ( ALLOCATED( EventLen  ) )  DEALLOCATE( EventLen  )
   IF ( ALLOCATED( pkCTKE    ) )  DEALLOCATE( pkCTKE    )

ENDIF !WrACT


   ! Are we generating FF AeroDyn/BLADED files OR Tower Files?

IF ( WrBLFF .OR. WrADTWR .OR. WrADFF )  THEN

      ! Calculate mean value & turb intensity of U-component of the interpolated hub point (for comparison w/ AeroDyn output)

   IF (ExtraHubPT) THEN

         ! Get points for bi-linear interpolation
      ZLo_YLo   = ( NumGrid_Z2 - 1 )*NumGrid_Y + NumGrid_Y2
      ZHi_YLo   = ( NumGrid_Z2     )*NumGrid_Y + NumGrid_Y2
      ZLo_YHi   = ( NumGrid_Z2 - 1 )*NumGrid_Y + NumGrid_Y2 + 1
      ZHi_YHi   = ( NumGrid_Z2     )*NumGrid_Y + NumGrid_Y2 + 1
    
      TmpZ      = (HubHt - Z(NumGrid_Z2))/GridRes_Z
      TmpY      = ( 0.0  - Y(NumGrid_Y2))/GridRes_Y
      CGridSum  = 0.0
      CGridSum2 = 0.0

      DO IT=1,NumSteps
         
      ! Interpolate within the grid for this time step.

         Tmp_YL_Z  = ( V( IT, ZHi_YLo, 1 ) - V( IT, ZLo_YLo, 1 ) )*TmpZ + V( IT, ZLo_YLo, 1 )
         Tmp_YH_Z  = ( V( IT, ZHi_YHi, 1 ) - V( IT, ZLo_YHi, 1 ) )*TmpZ + V( IT, ZLo_YHi, 1 )
         TmpV      = ( Tmp_YH_Z - Tmp_YL_Z )*TmpY + Tmp_YL_Z

         CGridSum  = CGridSum  + TmpV
         CGridSum2 = CGridSum2 + TmpV*TmpV
      ENDDO ! IT

      UGridMean = CGridSum/NumSteps
      UGridSig  = SQRT( ABS( (CGridSum2/NumSteps) - UGridMean*UGridMean ) )
      UGridTI   = 100.0*UGridSig/UGridMean

      FormStr = "(//,'U-component (X) statistics from the interpolated hub point:',/)"

   ELSE

      UGridMean = UXBar
      UGridTI = 100.0*UXSig/UXBar      
      FormStr = "(//,'U-component (X) statistics from the hub grid point:',/)"

   ENDIF

      ! Put the average statistics of the four center points in the summary file.

   WRITE (US,FormStr)

   FormStr = "(3X,A,' =',F9.4,A)"
   WRITE(US,FormStr)  'Mean' , UGridMean, ' m/s'
   WRITE(US,FormStr)  'TI  ' , UGridTI  , ' %'

   IF ( WrADFF ) THEN
      CALL WrBinTURBSIM
   END IF   
   
   IF ( WrBLFF .OR. (WrADTWR .AND. .NOT. WrADFF) ) THEN

         ! We need to take into account the shear across the grid in the sigma calculations for scaling the data, 
         ! and ensure that 32.767*Usig >= |V-UHub| so that we don't get values out of the range of our scaling values
         ! in this BLADED-style binary output.  TmpU is |V-UHub|
      USig = MAX(USig,0.05*TmpU)
      CALL WrBinBLADED(USig, VSig, WSig)

   ENDIF
   
ENDIF


   ! Are we generating FF formatted files?

IF ( WrFmtFF )  THEN
   CALL WrFormattedFF(HubIndx)
ENDIF ! ( WrFmtFF )


   ! Deallocate the temporary V, Y, Z, and IYmax arrays.

IF ( ALLOCATED( V               ) )  DEALLOCATE( V               )
IF ( ALLOCATED( Y               ) )  DEALLOCATE( Y               )
IF ( ALLOCATED( Z               ) )  DEALLOCATE( Z               )
IF ( ALLOCATED( IYmax           ) )  DEALLOCATE( IYmax           )


IF (DEBUG_OUT) THEN
   CLOSE( UD )       ! Close the debugging file
ENDIF

WRITE ( US, '(/"Nyquist frequency of turbulent wind field =      ",F8.3, " Hz")' ) 1.0 / (2.0 * TimeStep)
IF ( WrACT .AND. EventTimeStep > 0.0 ) THEN
   WRITE ( US, '( "Nyquist frequency of coherent turbulent events = ",F8.3, " Hz")' ) 1.0 / (2.0 * EventTimeStep)
ENDIF


   ! Request CPU-time used.

CALL CPU_TIME ( CPUtime )

FormStr = "(//,'Processing complete.  ',A,' CPU seconds used.')"
WRITE (US,FormStr)  TRIM( Flt2LStr( CPUtime ) )

!BONNIE: ****************************
!Time = TIMEF()
!PRINT *, Time
!FormStr = "(//,'BONNIE TEST.  ',A,' seconds for completion.')"
!WRITE (US,FormStr)  TRIM( Flt2LStr( CPUtime ) )
!END BONNIE ************************************

CLOSE ( US )


CALL WrScr1  ( ' Processing complete.  '//TRIM( Flt2LStr( CPUtime ) )//' CPU seconds used.' )
CALL NormStop


END PROGRAM
