!=======================================================================
PROGRAM TurbSim

   ! A turbulence simulator developed by contractors and staff at the
   ! National Renewable Energy Laboratory, Golden, Colorado
   !
   ! v1.0a-bjj       15-Mar-2004  B. Jonkman  
   ! v1.0             4-Nov-2005  B. Jonkman
   ! v1.01           24-Jan-2006  B. Jonkman  (NWTC subs v1.00a-mlb)
   ! v1.10           10-Apr-2006  B. Jonkman  (NWTC subs v1.12)
   ! v1.20           20-Oct-2006  B. Jonkman  (NWTC subs v1.12)
   ! v1.21            1-Feb-2007  B. Jonkman  (NWTC subs v1.12)
   ! v1.30            4-Apr-2008  B. Jonkman  (NWTC subs v1.01.09)
   ! v1.40           12-Sep-2008  B. Jonkman  (NWTC subs v1.01.09)
   ! v1.41h          11-Jun-2009  B. Jonkman  (NWTC subs v1.01.09)
   ! v1.50           25-Sep-2009  B. Jonkman  (NWTC subs v1.01.09)
   ! V1.06.00        21-Sep-2012  L. Kilcher & B. Jonkman (NWTC Library v1.04.01)
   ! v1.07.00a-bjj    9-Jul-2014  B. Jonkman  (NWTC Library v1.04.02a-bjj)
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
USE                        ModifiedvKrm_mod
USE                        TSMods
USE                        TSsubs
USE TS_FileIO
   USE TS_RandNum
   USE TS_Profiles
   USE TS_VelocitySpectra
use TS_CohStructures

!BONNIE:*****************************
! USE    IFPORT, ONLY: TIMEF ! Wall Clock Time
!BONNIE:*****************************    


IMPLICIT                   NONE


   ! Declare local variables

REAL(DbKi)              ::  CGridSum                        ! The sums of the velocity components at the points surrounding the hub (or at the hub if it's on the grid)
REAL(DbKi)              ::  CGridSum2                       ! The sums of the squared velocity components at the points surrouding the hub 

REAL(ReKi)              ::  UVMax                           ! Maximum u'v' Reynolds Stress at the hub
REAL(ReKi)              ::  UWMax                           ! Maximum u'w' Reynolds Stress at the hub
REAL(ReKi)              ::  VWMax                           ! Maximum v'w' Reynolds Stress at the hub

REAL(ReKi)              ::  USig                            ! Standard deviation of the u-component wind speed at the hub
REAL(ReKi)              ::  VSig                            ! Standard deviation of the v-component wind speed at the hub
REAL(ReKi)              ::  WSig                            ! Standard deviation of the w-component wind speed at the hub
REAL(ReKi)              ::  UXSig                           ! Standard deviation of the U-component wind speed at the hub



REAL(ReKi)              ::  UVsum                           ! The sum of the u'v' Reynolds stress component at the hub
REAL(ReKi)              ::  UWsum                           ! The sum of the u'w' Reynolds stress component at the hub
REAL(ReKi)              ::  VWsum                           ! The sum of the v'w' Reynolds stress component at the hub

REAL(DbKi)              ::  USum2                           ! The sum of the squared u-component wind speed at the hub
REAL(DbKi)              ::  VSum2                           ! The sum of the squared v-component wind speed at the hub
REAL(DbKi)              ::  WSum2                           ! The sum of the squared w-component wind speed at the hub
REAL(DbKi)              ::  UXBar                           ! The mean U-component (u rotated; x-direction) wind speed at the hub

REAL(ReKi)              ::  UGridMean                       ! Average wind speed at the points surrounding the hub
REAL(ReKi)              ::  UGridSig                        ! Standard deviation of the wind speed at the points surrounding the hub
REAL(ReKi)              ::  UGridTI                         ! Turbulent Intensity of the points surrounding the hub


REAL(ReKi)              ::  CPUtime                         ! Contains the number of seconds since the start of the program
REAL(ReKi)              ::  DelF5                           ! half of the delta frequency, used to discretize the continuous PSD at each point
REAL(ReKi)              ::  DelF                            ! Delta frequency
REAL(ReKi)              ::  INumSteps                       ! Multiplicative Inverse of the Number of time Steps
REAL(ReKi)              ::  TmpU                            ! Temporarily holds the value of the u component
REAL(ReKi)              ::  TmpV                            ! Temporarily holds the value of the v component
REAL(ReKi)              ::  TmpW                            ! Temporarily holds the value of the w component
REAL(ReKi)              ::  TmpY                            ! Temp variable for interpolated hub point
REAL(ReKi)              ::  TmpZ                            ! Temp variable for interpolated hub point
REAL(ReKi)              ::  Tmp_YL_Z                        ! Temp variable for interpolated hub point
REAL(ReKi)              ::  Tmp_YH_Z                        ! Temp variable for interpolated hub point

REAL(ReKi)              ::  HorVar                          ! Variables used when DEBUG_OUT is set
REAL(ReKi)              ::  ROT                             ! Variables used when DEBUG_OUT is set
REAL(ReKi)              ::  Total                           ! Variables used when DEBUG_OUT is set
REAL(ReKi)              ::  TotalU                          ! Variables used when DEBUG_OUT is set             
REAL(ReKi)              ::  TotalV                          ! Variables used when DEBUG_OUT is set
REAL(ReKi)              ::  TotalW                          ! Variables used when DEBUG_OUT is set

REAL(ReKi)              :: this_HFlowAng                    ! Horizontal flow angle.
REAL(ReKi)              :: v3(3)                            ! temporary 3-component array containing velocity
REAL(ReKi)              ::  ActualSigma(3)                  ! actual standard deviations at the hub
REAL(ReKi)              ::  HubFactor(3)                    ! factor used to scale standard deviations at the hub point



INTEGER                 ::  AllocStat                       ! The status of allocating space for variables
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


INTEGER(IntKi)          :: ErrStat     ! allocation status
CHARACTER(MaxMsgLen)    :: ErrMsg      ! error message


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


   ! Print out program name, version, and date.

CALL DispNVD


   ! Check for command line arguments.

CALL CheckArgs( InFile )


   ! Open input file and summary file.

CALL GetFiles


   ! Get input parameters.

CALL GetInput(ErrStat, ErrMsg)
CALL CheckError()


call WrSum_Heading1(US)


   ! Determine if coherent turbulence should be generated

!BONNIE: UPPER LIMIT ON RICH_NO?
IF ( WrACT ) THEN
   IF ( SpecModel == SpecModel_IECKAI .OR. &
        SpecModel == SpecModel_IECVKM .OR. &
        SpecModel == SpecModel_MODVKM .OR. &
        SpecModel == SpecModel_USER   .OR. &
        SpecModel == SpecModel_TIDAL   ) THEN
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
      WrADFF = .TRUE.
      CALL WrScr1( ' AeroDyn Full-Field files(.bts) will be generated along with the coherent turbulence file.' )
   ENDIF
      
ENDIF !WrAct

      ! Warn if EWM is used with incompatible times
      
IF ( ( p_IEC%IEC_WindType == IEC_EWM1 .OR. p_IEC%IEC_WindType == IEC_EWM50 ) .AND. & 
      ABS( REAL(600.0,ReKi) - MAX(p_grid%AnalysisTime,UsableTime) ) > REAL(90.0, ReKi) ) THEN
   CALL TS_Warn( ' The EWM parameters are valid for 10-min simulations only.', .TRUE. )
ENDIF        


      ! Warn if Periodic is used with incompatible settings
      
IF ( Periodic .AND. .NOT. EqualRealNos(p_grid%AnalysisTime, UsableTime) ) THEN
   CALL TS_Warn( ' Periodic output files will not be generated when AnalysisTime /= UsableTime. Setting Periodic = .FALSE.', &
                 WrSum=.TRUE.)
   Periodic = .FALSE.
END IF
      

   ! Open appropriate output files.  We will open formatted FF files later, if requested.
   ! Mention the files in the summary file.

IF ( WrBHHTP .OR. WrFHHTP .OR. WrADHH .OR. WrADFF .OR. WrFmtFF .OR. WrADTWR .OR. WrACT .OR. WrBLFF )  THEN

   FormStr = "( // 'You have requested that the following file(s) be generated:' / )"
   WRITE (US,FormStr)
   CALL WrScr1 ( ' You have requested that the following file(s) be generated:' )

   IF ( WrBHHTP )  THEN   

      CALL OpenUOutfile ( UGTP , TRIM( RootName)//'.bin', ErrStat, ErrMsg )  ! just making sure it can be opened (not locked elsewhere)
      CLOSE(UGTP)
      CALL CheckError()
!      CALL OpenBOutFile ( UGTP, TRIM( RootName)//'.bin', ErrStat, ErrMsg )
      FormStr = "( 3X , A , ' (hub-height binary turbulence-parameter file)' )"
      WRITE (US,FormStr)  TRIM( RootName)//'.bin'
      CALL WrScr ( '    '//TRIM( RootName)//'.bin (a binary hub-height turbulence-parameter file)' )

   ENDIF

   IF ( WrFHHTP )  THEN     

      CALL OpenFOutFile ( UFTP, TRIM( RootName)//'.dat', ErrStat, ErrMsg ) ! just making sure it can be opened (not locked elsewhere)
      CLOSE( UFTP )
      CALL CheckError()
      FormStr = "( 3X , A , ' (formatted turbulence-parameter file)' )"
      WRITE (US,FormStr)  TRIM( RootName)//'.dat'
      CALL WrScr ( '    '//TRIM( RootName)//'.dat (a formatted turbulence-parameter file)' )

   ENDIF

   IF ( WrADHH )  THEN     

      CALL OpenFOutFile ( UAHH, TRIM( RootName)//'.hh' ) ! just making sure it can be opened (not locked elsewhere)
      CLOSE( UAHH )
      FormStr = "( 3X , A , '  (AeroDyn hub-height file)' )"
      WRITE (US,FormStr)  TRIM( RootName)//'.hh'
      CALL WrScr ( '    '//TRIM( RootName)//'.hh  (an AeroDyn hub-height file)' )

   ENDIF

   IF ( WrADFF )  THEN

      CALL OpenBOutFile ( UAFFW, TRIM(RootName)//'.bts', ErrStat, ErrMsg )
      CALL CheckError()

      FormStr = "( 3X , A , ' (AeroDyn/TurbSim full-field wind file)' )"
      WRITE (US,FormStr)  TRIM( RootName)//'.bts'
      CALL WrScr ( '    '//TRIM( RootName)//'.bts (an AeroDyn/TurbSim full-field wind file)' )

   ENDIF

   IF ( WrBLFF )  THEN

      CALL OpenBOutFile ( UBFFW, TRIM(RootName)//'.wnd', ErrStat, ErrMsg )
      CALL CheckError()

      FormStr = "( 3X , A , ' (AeroDyn/BLADED full-field wnd file)' )"
      WRITE (US,FormStr)  TRIM( RootName)//'.wnd'
      CALL WrScr ( '    '//TRIM( RootName)//'.wnd (an AeroDyn/BLADED full-field wnd file)' )

   ENDIF

   IF ( WrADTWR .AND. (WrBLFF .OR. .NOT. WrADFF) )  THEN

      CALL OpenBOutFile ( UATWR, TRIM( RootName )//'.twr', ErrStat, ErrMsg )
      CALL CheckError()

      FormStr = "( 3X , A , ' (Binary tower twr file)' )"
      WRITE (US,FormStr)  TRIM( RootName)//'.twr'
      CALL WrScr ( '    '//TRIM( RootName)//'.twr (a binary tower file)' )

   ENDIF

   IF ( WrACT ) THEN

      CALL OpenFOutFile ( UACTTS, TRIM( RootName )//'.cts', ErrStat, ErrMsg  )         ! just making sure it can be opened (not locked elsewhere)
      close(UACTTS)
      CALL CheckError()
      
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

IF ( Periodic ) THEN
   p_grid%NumSteps    = CEILING( p_grid%AnalysisTime / p_grid%TimeStep )
   NumSteps4   = ( p_grid%NumSteps - 1 )/4 + 1
   p_grid%NumSteps    = 4*PSF( NumSteps4 , 9 )  ! >= 4*NumSteps4 = NumOutSteps + 3 - MOD(NumOutSteps-1,4) >= NumOutSteps
   p_grid%NumOutSteps = p_grid%NumSteps
ELSE
   p_grid%NumOutSteps = CEILING( ( UsableTime + p_grid%GridWidth/UHub )/p_grid%TimeStep )
   p_grid%NumSteps    = MAX( CEILING( p_grid%AnalysisTime / p_grid%TimeStep ), p_grid%NumOutSteps )
   NumSteps4   = ( p_grid%NumSteps - 1 )/4 + 1
   p_grid%NumSteps    = 4*PSF( NumSteps4 , 9 )  ! >= 4*NumSteps4 = NumOutSteps + 3 - MOD(NumOutSteps-1,4) >= NumOutSteps
END IF

INumSteps   = 1.0/p_grid%NumSteps

p_grid%NumFreq     = p_grid%NumSteps/2
DelF        = 1.0/( p_grid%NumSteps*p_grid%TimeStep )
DelF5       = 0.5*DelF


   !  Allocate/Initialize the array of frequencies.

call AllocAry( p_grid%Freq, p_grid%NumFreq, 'Freq (frequency array)', ErrStat, ErrMsg)
CALL CheckError()

DO IFreq=1,p_grid%NumFreq
   p_grid%Freq(IFreq) = IFreq*DelF
ENDDO 


   ! If we are going to generate debug output, open the debug file.

IF (DEBUG_OUT) THEN
   CALL OpenFOutFile ( UD, TRIM( RootName)//'.dbg' )
   
   WRITE (UD,*)  'UHub=',UHub
   WRITE (UD,*)  'NumOutSteps, NumSteps, INumSteps=', p_grid%NumOutSteps, p_grid%NumSteps, INumSteps
   
   ROT     = 1.0/( p_grid%NumGrid_Y*p_grid%TimeStep ) 
   WRITE (UD,*)  'TimeStep,DelF,NumFreq,ROT=',p_grid%TimeStep,DelF,p_grid%NumFreq,ROT
   WRITE (UD,*)  'PI,Deg2Rad=',Pi, D2R
   
   WRITE (UD,*)  'DelF=',DelF
   WRITE (UD,*)  'Freq(1) = ',p_grid%Freq(1),'  Freq(NumFreq) = ',p_grid%Freq(p_grid%NumFreq)
   
ENDIF


IF ( p_grid%NumSteps < 8 )  THEN
   CALL TS_Abort ( 'Too few (less than 8) time steps.  Increase the usable length of the time series or decrease the time step.' )
ENDIF



   ! Define the other parameters for the time series.
CALL CreateGrid( p_grid, NumGrid_Y2, NumGrid_Z2, NSize, ErrStat, ErrMsg )
CALL CheckError()
      


IF (DEBUG_OUT) THEN
   WRITE (UD,*)  'GridRes=', p_grid%GridRes_Y, ' x ', p_grid%GridRes_Z
   WRITE (UD,*)  ' '
   WRITE (UD,*)  'U='   
ENDIF


   !  Wind speed:
CALL AllocAry(U,     p_grid%ZLim, 'u (steady, u-component winds)', ErrStat, ErrMsg )
CALL CheckError()

IF ( WindProfileType(1:3) == 'API' )  THEN
   U = getVelocityProfile( U_Ref, H_Ref, p_grid%Z, RotorDiameter, PROFILE_TYPE=WindProfileType)
ELSE 
   U = getVelocityProfile( UHub,  p_grid%HubHt, p_grid%Z, RotorDiameter, PROFILE_TYPE=WindProfileType) 
ENDIF

   ! Wind Direction:
IF ( INDEX( 'JU', WindProfileType(1:1) ) > 0 ) THEN
   CALL AllocAry(WindDir_profile, p_grid%ZLim, 'WindDir_profile (wind direction profile)', ErrStat, ErrMsg )                  ! Allocate the array for the wind direction profile      
   CALL CheckError()
   
   WindDir_profile = getDirectionProfile(p_grid%Z, WindProfileType)
END IF

IF ( SpecModel == SpecModel_GP_LLJ) THEN

      ! Allocate the arrays for the z/l and ustar profile
      
   CALL AllocAry(ZL_profile,    p_grid%ZLim, 'ZL_profile (z/l profile)', ErrStat, ErrMsg )
   CALL CheckError()
   CALL AllocAry(Ustar_profile, p_grid%ZLim, 'Ustar_profile (friction velocity profile)', ErrStat, ErrMsg )         
   CALL CheckError()

   ZL_profile(:)    = getZLARY(    U, p_grid%Z )
   Ustar_profile(:) = getUstarARY( U, p_grid%Z )
   
END IF


IF ( INDEX( 'JU', WindProfileType(1:1) ) > 0 ) THEN
   IF (p_grid%ExtraHubPT) THEN    
      JZ = p_grid%NumGrid_Z+1  ! This is the index of the Hub-height parameters if the hub height is not on the grid
   ELSE
      JZ = NumGrid_Z2
   ENDIF
   HH_HFlowAng = WindDir_profile(JZ)   
ELSE
   HH_HFlowAng = HFlowAng
END IF



IF (DEBUG_OUT) THEN
   DO IZ = 1,p_grid%ZLim
      WRITE (UD,*)  p_grid%Z(IZ), U(IZ)
   ENDDO
ENDIF


  
IF ( ( SpecModel  == SpecModel_IECKAI ) .OR. &
     ( SpecModel  == SpecModel_IECVKM ) .OR. &
     ( SpecModel  == SpecModel_MODVKM ) .OR. &
     ( SpecModel  == SpecModel_API    ) )  THEN  ! ADDED BY YGUO on April 192013 snow day!!!

      ! If IECKAI or IECVKM spectral models are specified, determine turb intensity 
      ! and slope of Sigma wrt wind speed from IEC turbulence characteristic, 
      ! IECTurbC = A, B, or C or from user specified quantity.
      
   CALL CalcIECScalingParams(p_IEC)
                  
   IF (DEBUG_OUT) THEN
      WRITE (UD,*)  'SigmaIEC,TurbInt=',p_IEC%SigmaIEC(1),p_IEC%TurbInt
   ENDIF

ELSE
    p_IEC%SigmaIEC = 0
    p_IEC%Lambda = 0
    p_IEC%IntegralScale = 0
    p_IEC%LC = 0.0    ! The length scale is not defined for the non-IEC models
ENDIF !  TurbModel  == 'IECKAI', 'IECVKM', 'API', or 'MODVKM'
   
CALL WrSum_SpecModel( US, p_IEC, p_grid%Z(NumGrid_Z2) )




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






   !  Allocate the turbulence PSD array.

CALL AllocAry( S,    p_grid%NumFreq,p_grid%NPoints,3, 'S (turbulence PSD)',ErrStat, ErrMsg )
CALL CheckError()

   !  Allocate the work array.

CALL AllocAry( Work, p_grid%NumFreq,3,             'Work',              ErrStat, ErrMsg )
CALL CheckError()

IF (PSD_OUT) THEN
   CALL OpenFOutFile ( UP, TRIM( RootName )//'.psd')
   WRITE (UP,"(A)")  'PSDs '
   FormStr = "( A4,'"//TAB//"',A4,"//TRIM( Int2LStr( p_grid%NumFreq ) )//"('"//TAB//"',G10.4) )"
   WRITE (UP, FormStr)  'Comp','Ht', p_grid%Freq(:)
   FormStr = "( I4,"//TRIM( Int2LStr( p_grid%NumFreq+1 ) )//"('"//TAB//"',G10.4) )"
ENDIF


   ! Allocate and initialize the DUDZ array for MHK models (TIDAL and RIVER)

IF ( SpecModel == SpecModel_TIDAL .OR. SpecModel == SpecModel_RIVER ) THEN
      ! Calculate the shear, DUDZ, for all heights.
   ALLOCATE ( DUDZ(p_grid%ZLim) , STAT=AllocStat ) ! Shear
   DUDZ(1)=(U(2)-U(1))/(p_grid%Z(2)-p_grid%Z(1))
   DUDZ(p_grid%ZLim)=(U(p_grid%ZLim)-U(p_grid%ZLim-1))/(p_grid%Z(p_grid%ZLim)-p_grid%Z(p_grid%ZLim-1))
   DO I = 2,p_grid%ZLim-1
      DUDZ(I)=(U(I+1)-U(I-1))/(p_grid%Z(I+1)-p_grid%Z(I-1))
   ENDDO
ENDIF

   ! Calculate the single point Power Spectral Densities. 

JZ = 0   ! The index for numbering the points on the grid
DO IZ=1,p_grid%ZLim

         ! The continuous, one-sided PSDs are evaluated at discrete
         ! points and the results are stored in the "Work" matrix.


!bonnie: fix this so the the IEC runs differently?? It doesn't need as large an array here anymore...
      IF ( SpecModel == SpecModel_IECKAI .or. SpecModel == SpecModel_IECVKM ) THEN
         CALL PSDcal( p_grid%HubHt, UHub)   ! Added by Y.G.  !!!! only hub height is acceptable !!!!
      ELSEIF ( SpecModel == SpecModel_TIDAL .OR. SpecModel == SpecModel_RIVER ) THEN ! HydroTurbSim specific.
         !print *, Ustar
         !Sigma_U2=(TurbIntH20*U(IZ))**2 ! A fixed value of the turbulence intensity.  Do we want to implement this?
         Sigma_U2=4.5*Ustar*Ustar*EXP(-2*p_grid%Z(IZ)/H_ref)
         Sigma_V2=0.5*Sigma_U2
         Sigma_W2=0.2*Sigma_U2
         CALL PSDCal( p_grid%Z(IZ) , DUDZ(IZ) )
      ELSEIF ( ALLOCATED( ZL_profile ) ) THEN
         CALL PSDcal( p_grid%Z(IZ), U(IZ), ZL_profile(IZ), Ustar_profile(IZ) )
      ELSE
         CALL PSDcal( p_grid%Z(IZ), U(IZ) )
      ENDIF               

      IF (DEBUG_OUT) THEN
            TotalU = 0.0
            TotalV = 0.0
            TotalW = 0.0

            DO IFreq=1,p_grid%NumFreq
               TotalU = TotalU + Work(IFreq,1)
               TotalV = TotalV + Work(IFreq,2)
               TotalW = TotalW + Work(IFreq,3)
            ENDDO ! IFreq

            Total  = SQRT( TotalU*TotalU + TotalV*TotalV )
            TotalU = TotalU*DelF
            TotalV = TotalV*DelF
            TotalW = TotalW*DelF
            HorVar = Total *DelF

            WRITE(UD,*)  'At H=', p_grid%Z(IZ)
            WRITE(UD,*)  '   TotalU=', TotalU, '  TotalV=', TotalV, '  TotalW=', TotalW, ' (m/s^2)'
            WRITE(UD,*)  '   HorVar=',HorVar,' (m/s^2)', '   HorSigma=',SQRT(HorVar),' (m/s)'
      ENDIF


         ! Discretize the continuous PSD and store it in matrix "S"
           
      DO IVec=1,3
         
         Indx = JZ

         DO IY=1,p_grid%IYmax(IZ)   
          
            Indx = Indx + 1
 
            DO IFreq=1,p_grid%NumFreq
               S(IFreq,Indx,IVec) = Work(IFreq,IVec)*DelF5
            ENDDO ! IFreq

         ENDDO !IY

      ENDDO ! IVec

      JZ = JZ + p_grid%IYmax(IZ)     ! The next starting index at height IZ + 1

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

CALL AllocAry( PhaseAngles, p_grid%NPoints, p_grid%NumFreq, 3, 'Random Phases', ErrStat, ErrMsg )
CALL CheckError()

CALL RndPhases(p_RandNum, OtherSt_RandNum, PhaseAngles, p_grid%NPoints, p_grid%NumFreq, ErrStat, ErrMsg)
CALL CheckError()

WRITE(US,"(//,'Harvested Random Seeds after Generation of the Random Numbers:',/)")

DO I = 1,SIZE( OtherSt_RandNum%nextSeed  )
   WRITE(US,"(I13,' Harvested seed #',I2)")  OtherSt_RandNum%nextSeed (I), I
END DO

IF (ALLOCATED(OtherSt_RandNum%nextSeed) ) DEALLOCATE(OtherSt_RandNum%nextSeed)


   
CALL AllocAry( TRH, MAX(p_grid%NumSteps,NSize), 'TRH (transfer-function matrix)', ErrStat, ErrMsg) !  Allocate the transfer-function matrix.
CALL CheckError()

CALL AllocAry( V, p_grid%NumSteps, p_grid%NPoints, 3, 'V (velocity)', ErrStat, ErrMsg) !  Allocate the array that contains the velocities.
CALL CheckError()


   ! Calculate the transfer function matrices from the spectral matrix (the fourier coefficients).

CALL WrScr ( ' Calculating the spectral and transfer function matrices:' )


IF (SpecModel /= SpecModel_NONE) THEN                         ! MODIFIED BY Y GUO
    IF (SpecModel == SpecModel_API) THEN
        CALL CohSpecVMat_API(p_IEC%LC, NSize) 
    ELSE
        CALL CohSpecVMat(p_IEC%LC, NSize)
    ENDIF
ENDIF


   ! Deallocate the Freq, S, and RandPhases arrays and the spectral matrix

IF ( ALLOCATED( p_grid%Freq        ) )  DEALLOCATE( p_grid%Freq        )
IF ( ALLOCATED( S           ) )  DEALLOCATE( S           )
IF ( ALLOCATED( PhaseAngles ) )  DEALLOCATE( PhaseAngles )

   !  Allocate the FFT working storage and initialize its variables

CALL InitFFT( p_grid%NumSteps, FFT_Data, ErrStat=ErrStat )
CALL CheckError()


   ! Get the stationary-point time series.

CALL WrScr ( ' Generating time series for all points:' )

DO IVec=1,3

   CALL WrScr ( '    '//Comp(IVec)//'-component' )

   DO Indx=1,p_grid%NPoints    !NTotB

         ! Overwrite the first point with zero.  This sets the real (and 
         ! imaginary) part of the steady-state value to zero so that we 
         ! can add in the mean value later.

      TRH(1)=0.0

      DO IT=2,p_grid%NumSteps-1
         TRH(IT) = V(IT-1,Indx,IVec)
      ENDDO ! IT

         ! Now, let's add a complex zero to the end to set the power in the Nyquist
         ! frequency to zero.

      TRH(p_grid%NumSteps) = 0.0


         ! perform FFT

      CALL ApplyFFT( TRH, FFT_Data, ErrStat )
      CALL CheckError()
        
      V(1:p_grid%NumSteps,Indx,IVec) = TRH(1:p_grid%NumSteps)

   ENDDO ! Indx

ENDDO ! IVec

   ! Deallocate the TRH array and the FFT working storage.
IF ( ALLOCATED( TRH  ) )  DEALLOCATE( TRH  )

CALL ExitFFT(FFT_Data, ErrStat )
CALL CheckError()


   ! Crossfeed cross-axis components to u', v', w' components and scale IEC models if necessary

SELECT CASE ( SpecModel )
!MLB: There does not seem to be a CASE for TurbModel=="API".
      
CASE (SpecModel_GP_LLJ, &
      SpecModel_NWTCUP, &
      SpecModel_SMOOTH, &
      SpecModel_WF_UPW, &
      SpecModel_WF_07D, &
      SpecModel_WF_14D, &
      SpecModel_USRVKM, &
      SpecModel_TIDAL,  &
      SpecModel_RIVER   ) ! Do reynolds stress for HYDRO also.
               
            ! Calculate coefficients for obtaining "correct" Reynold's stresses at the hub
         UWsum = 0.0
         UVsum = 0.0
         VWsum = 0.0
         USum2 = 0.0
         VSum2 = 0.0
         WSum2 = 0.0
         
         DO IT = 1,p_grid%NumSteps
            UWsum = UWsum + V(IT,p_grid%HubIndx,1) * V(IT,p_grid%HubIndx,3)
            UVsum = UVsum + V(IT,p_grid%HubIndx,1) * V(IT,p_grid%HubIndx,2)
            VWsum = VWsum + V(IT,p_grid%HubIndx,2) * V(IT,p_grid%HubIndx,3)
            USum2 = USum2 + V(IT,p_grid%HubIndx,1) * V(IT,p_grid%HubIndx,1)
            VSum2 = VSum2 + V(IT,p_grid%HubIndx,2) * V(IT,p_grid%HubIndx,2)
            WSum2 = WSum2 + V(IT,p_grid%HubIndx,3) * V(IT,p_grid%HubIndx,3)
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
                  
         DO Indx = 1,p_grid%NPoints
            DO IT = 1, p_grid%NumSteps
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

   CASE ( SpecModel_IECKAI , SpecModel_IECVKM )

      IF (p_IEC%ScaleIEC > 0) THEN

         CALL TimeSeriesScaling_IEC(p_grid, V, p_IEC%ScaleIEC, p_IEC%SigmaIEC, ActualSigma, HubFactor)
         
         WRITE( US, "(//,'Scaling statistics from the hub grid point:',/)" )                  
         WRITE( US, "(2X,'Component  Target Sigma (m/s)  Simulated Sigma (m/s)  Scaling Factor')" )
         WRITE( US, "(2X,'---------  ------------------  ---------------------  --------------')" )
         
         FormStr = "(2X,3x,A,7x,f11.3,9x,f12.3,11x,f10.3)"                        
         DO IVec = 1,3
            WRITE( US, FormStr) Comp(IVec)//"'", p_IEC%SigmaIEC(IVec), ActualSigma(IVec), HubFactor(IVec)
         END DO
           
      ENDIF !ScaleIEC

END SELECT



!..............................................................................
! Write hub-height output files (before adding mean and rotating final results)
!..............................................................................

IF ( WrADHH )  THEN
   CALL WrHH_ADtxtfile(p_IEC%TurbInt)   
END IF

IF ( WrBHHTP )  THEN
   CALL WrHH_binary(ErrStat, ErrMsg)
   CALL CheckError()
END IF

IF ( WrFHHTP )  THEN
   CALL WrHH_text()   
END IF

   ! Write statistics of the run to the summary file:
CALL WrSum_Stats(USig, VSig, WSig, UXBar, UXSig, ErrStat, ErrMsg)
CALL CheckError()
p_CohStr%WSig=WSig

!..............................................................................
! Add mean wind to u' components and rotate to inertial reference  
!  frame coordinate system
!..............................................................................
II = 0
DO IZ=1,p_grid%ZLim   

   IF ( ALLOCATED( WindDir_profile ) ) THEN  ! The horizontal flow angle changes with height
      this_HFlowAng = WindDir_profile(IZ)
   ELSE
      this_HFlowAng = HFlowAng
   ENDIF      

  DO IY=1,p_grid%IYmax(IZ)  

      II = II + 1

      DO IT=1,p_grid%NumSteps

            ! Add mean wind speed to the streamwise component and
            ! Rotate the wind to the X-Y-Z (inertial) reference frame coordinates:
            
         v3 = V(IT,II,:)
         CALL CalculateWindComponents( v3, U(IZ), this_HFlowAng, VFlowAng, V(IT,II,:) )
                        
      ENDDO ! IT
  ENDDO ! IY
ENDDO ! IZ

TmpU = MAX( ABS(MAXVAL(U)-UHub), ABS(MINVAL(U)-UHub) )  !Get the range of wind speed values for scaling in BLADED-format .wnd files

   ! Deallocate memory for the matrix of the steady, u-component winds.

IF ( ALLOCATED( U               ) )  DEALLOCATE( U               )
IF ( ALLOCATED( WindDir_profile ) )  DEALLOCATE( WindDir_profile )

!..............................................................................
! Are we generating a coherent turbulent timestep file?
!..............................................................................
   
IF ( WrACT ) THEN
   
   CALL CohStr_WriteCTS(ErrStat, ErrMsg)
   CALL CheckError()
   
         
      ! Write the number of separate events to the summary file

   IF (KHtest) THEN
      WRITE ( US,'(/)' )
   ELSE
      WRITE ( US,'(//A,F8.3," seconds")' ) 'Average expected time between events = ',y_CohStr%lambda
   ENDIF

   WRITE ( US, '(A,I8)'   )            'Number of coherent events            = ', y_CohStr%NumCTEvents_separate
   WRITE ( US, '(A,F8.3," seconds")')  'Predicted length of coherent events  = ', y_CohStr%ExpectedTime
   WRITE ( US, '(A,F8.3," seconds")')  'Length of coherent events            = ', y_CohStr%EventTimeSum
   WRITE ( US, '(A,F8.3," (m/s)^2")')  'Maximum predicted event CTKE         = ', y_CohStr%CTKE

   

ENDIF !WrACT

!..............................................................................
! Are we generating FF AeroDyn/BLADED files OR Tower Files?
!..............................................................................

IF ( WrBLFF .OR. WrADFF .OR. WrADTWR )  THEN !bjj: isn't WrADTWR redundant here?

      ! Calculate mean value & turb intensity of U-component of the interpolated hub point (for comparison w/ AeroDyn output)

   IF (p_grid%ExtraHubPT) THEN

         ! Get points for bi-linear interpolation
      ZLo_YLo   = ( NumGrid_Z2 - 1 )*p_grid%NumGrid_Y + NumGrid_Y2
      ZHi_YLo   = ( NumGrid_Z2     )*p_grid%NumGrid_Y + NumGrid_Y2
      ZLo_YHi   = ( NumGrid_Z2 - 1 )*p_grid%NumGrid_Y + NumGrid_Y2 + 1
      ZHi_YHi   = ( NumGrid_Z2     )*p_grid%NumGrid_Y + NumGrid_Y2 + 1
    
      TmpZ      = (p_grid%HubHt - p_grid%Z(NumGrid_Z2))/p_grid%GridRes_Z
      TmpY      = ( 0.0  - p_grid%Y(NumGrid_Y2))/p_grid%GridRes_Y
      CGridSum  = 0.0
      CGridSum2 = 0.0

      DO IT=1,p_grid%NumSteps
         
      ! Interpolate within the grid for this time step.

         Tmp_YL_Z  = ( V( IT, ZHi_YLo, 1 ) - V( IT, ZLo_YLo, 1 ) )*TmpZ + V( IT, ZLo_YLo, 1 )
         Tmp_YH_Z  = ( V( IT, ZHi_YHi, 1 ) - V( IT, ZLo_YHi, 1 ) )*TmpZ + V( IT, ZLo_YHi, 1 )
         TmpV      = ( Tmp_YH_Z - Tmp_YL_Z )*TmpY + Tmp_YL_Z

         CGridSum  = CGridSum  + TmpV
         CGridSum2 = CGridSum2 + TmpV*TmpV
      ENDDO ! IT

      UGridMean = CGridSum/p_grid%NumSteps
      UGridSig  = SQRT( ABS( (CGridSum2/p_grid%NumSteps) - UGridMean*UGridMean ) )
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
   CALL WrFormattedFF(p_grid%HubIndx)
ENDIF ! ( WrFmtFF )


   ! Deallocate the temporary V, Y, Z, and IYmax arrays.

IF ( ALLOCATED( V               ) )  DEALLOCATE( V               )
IF ( ALLOCATED( p_grid%Y        ) )  DEALLOCATE( p_grid%Y        )
IF ( ALLOCATED( p_grid%Z        ) )  DEALLOCATE( p_grid%Z        )
IF ( ALLOCATED( p_grid%IYmax    ) )  DEALLOCATE( p_grid%IYmax    )


IF (DEBUG_OUT) THEN
   CLOSE( UD )       ! Close the debugging file
ENDIF

WRITE ( US, '(/"Nyquist frequency of turbulent wind field =      ",F8.3, " Hz")' ) 1.0 / (2.0 * p_grid%TimeStep)
IF ( WrACT .AND. y_CohStr%EventTimeStep > 0.0 ) THEN
   WRITE ( US, '( "Nyquist frequency of coherent turbulent events = ",F8.3, " Hz")' ) 1.0 / (2.0 * y_CohStr%EventTimeStep)
ENDIF


   ! Request CPU-time used.

CALL CPU_TIME ( CPUtime )

FormStr = "(//,'Processing complete.  ',A,' CPU seconds used.')"
WRITE (US,FormStr)  TRIM( Num2LStr( CPUtime ) )

!BONNIE: ****************************
!Time = TIMEF()
!PRINT *, Time
!FormStr = "(//,'BONNIE TEST.  ',A,' seconds for completion.')"
!WRITE (US,FormStr)  TRIM( Flt2LStr( CPUtime ) )
!END BONNIE ************************************

CLOSE ( US )


CALL WrScr1  ( ' Processing complete.  '//TRIM( Num2LStr( CPUtime ) )//' CPU seconds used.' )
CALL NormStop


CONTAINS
!...........................................................
SUBROUTINE CheckError()
   IF (ErrStat /= ErrID_None) THEN
   
      IF (ErrStat >= AbortErrLev) THEN
         CALL TS_end()
         CALL TS_Abort( TRIM(ErrMSg) )
      ELSE
         CALL WrScr(TRIM(ErrMsg))
      END IF
   END IF
END SUBROUTINE CheckError

END PROGRAM
