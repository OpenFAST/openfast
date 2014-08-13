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


REAL(ReKi)              ::  USig                            ! Standard deviation of the u-component wind speed at the hub
REAL(ReKi)              ::  VSig                            ! Standard deviation of the v-component wind speed at the hub
REAL(ReKi)              ::  WSig                            ! Standard deviation of the w-component wind speed at the hub
REAL(ReKi)              ::  UXSig                           ! Standard deviation of the U-component wind speed at the hub




REAL(DbKi)              ::  UXBar                           ! The mean U-component (u rotated; x-direction) wind speed at the hub

REAL(ReKi)              ::  UGridMean                       ! Average wind speed at the points surrounding the hub
REAL(ReKi)              ::  UGridSig                        ! Standard deviation of the wind speed at the points surrounding the hub
REAL(ReKi)              ::  UGridTI                         ! Turbulent Intensity of the points surrounding the hub


REAL(ReKi)              ::  CPUtime                         ! Contains the number of seconds since the start of the program
REAL(ReKi)              ::  TmpU                            ! Temporarily holds the value of the u component
REAL(ReKi)              ::  TmpV                            ! Temporarily holds the value of the v component
REAL(ReKi)              ::  TmpY                            ! Temp variable for interpolated hub point
REAL(ReKi)              ::  TmpZ                            ! Temp variable for interpolated hub point
REAL(ReKi)              ::  Tmp_YL_Z                        ! Temp variable for interpolated hub point
REAL(ReKi)              ::  Tmp_YH_Z                        ! Temp variable for interpolated hub point




INTEGER                 ::  IT                              ! Index for time step
INTEGER                 ::  IY                              ! An index for the Y position of a point I
INTEGER                 ::  JZ                              ! An index for the Z position of a point J
INTEGER                 ::  NSize                           ! Size of the spectral matrix at each frequency

INTEGER                 ::  NumGrid_Y2                      ! Y Index of the hub (or the nearest point left the hub if hub does not fall on the grid)
INTEGER                 ::  NumGrid_Z2                      ! Z Index of the hub (or the nearest point below the hub if hub does not fall on the grid) 
INTEGER                 ::  ZHi_YHi                         ! Index for interpolation of hub point, if necessary
INTEGER                 ::  ZHi_YLo                         ! Index for interpolation of hub point, if necessary
INTEGER                 ::  ZLo_YHi                         ! Index for interpolation of hub point, if necessary
INTEGER                 ::  ZLo_YLo                         ! Index for interpolation of hub point, if necessary
INTEGER                 ::  UnOut                           ! unit for output files


CHARACTER(200)          :: InFile                           ! Name of the TurbSim input file.
CHARACTER(200)          :: FormStr                          ! String used to store format specifiers.


REAL(ReKi), ALLOCATABLE          :: PhaseAngles (:,:,:)                      ! The array that holds the random phases [number of points, number of frequencies, number of wind components=3].
REAL(ReKi), ALLOCATABLE          :: S           (:,:,:)                      ! The turbulence PSD array (NumFreq,NPoints,3).
REAL(ReKi), ALLOCATABLE          :: V           (:,:,:)                      ! An array containing the summations of the rows of H (NumSteps,NPoints,3).
REAL(ReKi), ALLOCATABLE          :: U           (:)                          ! The steady u-component wind speeds for the grid (ZLim).
REAL(ReKi), ALLOCATABLE          :: HWindDir    (:)                          ! A profile of horizontal wind angle (measure of wind direction with height)

!TYPE( TurbSim_ParameterType )    :: p

INTEGER(IntKi)          :: ErrStat     ! allocation status
CHARACTER(MaxMsgLen)    :: ErrMsg      ! error message


!BONNIE:*****************************
!    Time = TIMEF() ! Initialize the Wall Clock Time counter
!BONNIE:*****************************   



   ! ... Initialize NWTC Library (open console, set pi constants) ...
CALL NWTC_Init( ProgNameIN=TurbSim_Ver%Name, EchoLibVer=.FALSE. )       


   ! Print out program name, version, and date.

CALL DispNVD(TurbSim_Ver)


   ! Check for command line arguments.
InFile = 'TurbSim.inp'  ! default name for input file
CALL CheckArgs( InFile )


   ! Open input file and summary file.

CALL GetFiles( InFile, p%RootName, p%DescStr )


   ! Get input parameters.

CALL GetInput(InFile, ErrStat, ErrMsg)
CALL CheckError()

CALL WrSum_EchoInputs() 
call WrSum_UserInput(p%met,US)

call TS_ValidateInput(ErrStat, ErrMsg)
CALL CheckError()


   ! Open appropriate output files.  We will open formatted FF files later, if requested.
   ! Mention the files in the summary file.

IF ( ANY (p%WrFile) )  THEN
   CALL GetNewUnit( UnOut )

   WRITE (US,"( // 'You have requested that the following file(s) be generated:' / )")
   CALL WrScr1  ( ' You have requested that the following file(s) be generated:' )

   IF ( p%WrFile(FileExt_BIN) )  THEN   

!      CALL OpenBOutFile ( UnOut, TRIM( p%RootName)//'.bin', ErrStat, ErrMsg )
      CALL OpenUOutfile ( UnOut , TRIM( p%RootName)//'.bin', ErrStat, ErrMsg )  ! just making sure it can be opened (not locked elsewhere)
      CLOSE(UnOut)
      CALL CheckError()
      
      WRITE (US,"( 3X ,'"//TRIM( p%RootName)//".bin (a binary hub-height turbulence-parameter file)' )")  
      CALL WrScr ( '    '//TRIM( p%RootName)//'.bin (a binary hub-height turbulence-parameter file)' )

   ENDIF

   IF ( p%WrFile(FileExt_DAT) )  THEN     

      CALL OpenFOutFile ( UnOut, TRIM( p%RootName)//'.dat', ErrStat, ErrMsg ) ! just making sure it can be opened (not locked elsewhere)
      CLOSE( UnOut )
      CALL CheckError()
      
      WRITE (US, "( 3X ,'"//TRIM( p%RootName)//".dat (a formatted turbulence-parameter file)' )")  
      CALL WrScr ( '     '//TRIM( p%RootName)//'.dat (a formatted turbulence-parameter file)' )

   ENDIF

   IF ( p%WrFile(FileExt_HH) )  THEN     

      CALL OpenFOutFile ( UnOut, TRIM( p%RootName)//'.hh', ErrStat, ErrMsg ) ! just making sure it can be opened (not locked elsewhere)
      CLOSE( UnOut )
      CALL CheckError()
      
      WRITE (US,"( 3X ,'"//TRIM( p%RootName)//".hh  (an AeroDyn hub-height file)' )")
      CALL WrScr ( '    '//TRIM( p%RootName)//'.hh  (an AeroDyn hub-height file)' )

   ENDIF

   IF ( p%WrFile(FileExt_BTS) )  THEN

      CALL OpenBOutFile ( UnOut, TRIM(p%RootName)//'.bts', ErrStat, ErrMsg )
      CLOSE( UnOut )
      CALL CheckError()

      WRITE (US,"( 3X ,'"//TRIM( p%RootName)//".bts (an AeroDyn/TurbSim full-field wind file)' )")  
      CALL WrScr ( '    '//TRIM( p%RootName)//'.bts (an AeroDyn/TurbSim full-field wind file)' )

   ENDIF

   IF ( p%WrFile(FileExt_WND) )  THEN

      CALL OpenBOutFile ( UnOut, TRIM(p%RootName)//'.wnd', ErrStat, ErrMsg )
      CLOSE(UnOut)
      CALL CheckError()

      WRITE (US,"( 3X ,'"//TRIM( p%RootName)//".wnd (an AeroDyn/BLADED full-field wnd file)' )")  
      CALL WrScr ( '    '//TRIM( p%RootName)//'.wnd (an AeroDyn/BLADED full-field wnd file)' )

   ENDIF
   
   IF ( p%WrFile(FileExt_TWR) .AND. p%WrFile(FileExt_WND) )  THEN

      CALL OpenBOutFile ( UnOut, TRIM( p%RootName )//'.twr', ErrStat, ErrMsg )
      CLOSE(UnOut)       
      CALL CheckError()

      WRITE (US,"( 3X ,'"//TRIM( p%RootName)//".twr (a binary tower file)' )")  
      CALL WrScr ( '    '//TRIM( p%RootName)//'.twr (a binary tower file)' )

   ENDIF

   IF ( p%WrFile(FileExt_CTS) ) THEN      
      WRITE (US,"( 3X ,'"//TRIM( p%RootName)//".cts (a coherent turbulence time step file)' )")  
      CALL WrScr ( '    '//TRIM( p%RootName)//'.cts (a coherent turbulence time step file)' )
   ENDIF

   IF ( p%WrFile(FileExt_UVW) )  THEN
      WRITE (US,"( 3X ,'"//TRIM( p%RootName)//".u (a formatted full-field U-component file)' )")  
      CALL WrScr ( '    '//TRIM( p%RootName)//'.u (a formatted full-field U-component file)' )

      WRITE (US,"( 3X ,'"//TRIM( p%RootName)//".v (a formatted full-field V-component file)' )")  
      CALL WrScr ( '    '//TRIM( p%RootName)//'.v (a formatted full-field V-component file)' )

      WRITE (US,"( 3X ,'"//TRIM( p%RootName)//".w (a formatted full-field W-component file)' )")  
      CALL WrScr ( '    '//TRIM( p%RootName)//'.w (a formatted full-field W-component file)' )
   ENDIF

ELSE
   ErrStat = ErrID_Fatal
   ErrMsg  = 'You have requested no output.'
   CALL CheckError()
ENDIF



   ! Define the other parameters for the time series.
CALL CreateGrid( p%grid, NumGrid_Y2, NumGrid_Z2, NSize, ErrStat, ErrMsg )
CALL CheckError()
      

   !  Wind speed:
CALL AllocAry(U,     p%grid%ZLim, 'u (steady, u-component winds)', ErrStat, ErrMsg )
CALL CheckError()

IF ( p%met%WindProfileType(1:3) == 'API' )  THEN
   U = getVelocityProfile( p%met%URef, p%met%RefHt, p%grid%Z, p%grid%RotorDiameter)
ELSE 
   U = getVelocityProfile(   p%UHub,  p%grid%HubHt, p%grid%Z, p%grid%RotorDiameter) 
ENDIF

   ! Wind Direction:
CALL AllocAry(HWindDir, p%grid%ZLim, 'HWindDir (wind direction profile)', ErrStat, ErrMsg )                  ! Allocate the array for the wind direction profile      
CALL CheckError()
   
HWindDir = getDirectionProfile(p%grid%Z)
   
IF (p%grid%ExtraHubPT) THEN    
   JZ = p%grid%NumGrid_Z+1  ! This is the index of the Hub-height parameters if the hub height is not on the grid
ELSE
   JZ = NumGrid_Z2
ENDIF
p%met%HH_HFlowAng = HWindDir(JZ)



IF ( p%met%TurbModel_ID == SpecModel_GP_LLJ) THEN

      ! Allocate the arrays for the z/l and ustar profile
      
   CALL AllocAry(p%met%ZL_profile,    p%grid%ZLim, 'ZL_profile (z/l profile)', ErrStat, ErrMsg )
   CALL CheckError()
   CALL AllocAry(p%met%Ustar_profile, p%grid%ZLim, 'Ustar_profile (friction velocity profile)', ErrStat, ErrMsg )         
   CALL CheckError()

   p%met%ZL_profile(:)    = getZLARY(    U, p%grid%Z, p%met%Rich_No, p%met%ZL, p%met%L, p%met%ZLOffset, p%met%WindProfileType )
   p%met%Ustar_profile(:) = getUstarARY( U, p%grid%Z, p%met%UStarOffset, p%met%UStarSlope )
   
END IF

  
IF ( p%met%IsIECModel )  THEN  

      ! If IECKAI or IECVKM spectral models are specified, determine turb intensity 
      ! and slope of Sigma wrt wind speed from IEC turbulence characteristic, 
      ! IECTurbC = A, B, or C or from user specified quantity.
      
   CALL CalcIECScalingParams(p%IEC, p%grid%HubHt, p%UHub, p%met%InCDec, p%met%InCohB, p%met%TurbModel_ID)
                  
ELSE
    p%IEC%SigmaIEC = 0
    p%IEC%Lambda = 0
    p%IEC%IntegralScale = 0
    p%IEC%LC = 0.0    ! The length scale is not defined for the non-IEC models
ENDIF !  TurbModel  == 'IECKAI', 'IECVKM', 'API', or 'MODVKM'
   
CALL WrSum_SpecModel( US, U, HWindDir, p%IEC, p%grid%Z(NumGrid_Z2) )




IF (DEBUG_OUT) THEN
   ! If we are going to generate debug output, open the debug file.
   !UD = -1
   !CALL GetNewUnit( UD, ErrStat, ErrMsg )
   CALL OpenFOutFile ( UD, TRIM( p%RootName)//'.dbg', ErrStat, ErrMsg )
      CALL CheckError()      
   
   !bjj: This doesn't need to be part of the summary file, so put it in the debug file
   IF (p%met%WindProfileType(1:1) == 'J' ) THEN
      WRITE(UD,"(//'Jet wind profile Chebyshev coefficients:' / )") 
      WRITE(UD,"( 3X,'Order: ',11(1X,I10) )")  ( IY, IY=0,10 )
      WRITE(UD,"( 3X,'------ ',11(1X,'----------') )")
      WRITE(UD,"( 3X,'Speed: ',11(1X,E10.3) )")  ( p%met%ChebyCoef_WS(IY), IY=1,11 )
      WRITE(UD,"( 3X,'Angle: ',11(1X,E10.3) )")  ( p%met%ChebyCoef_WD(IY), IY=1,11 )   
   ENDIF
ENDIF



   !  Allocate the turbulence PSD array.

CALL AllocAry( S,    p%grid%NumFreq,p%grid%NPoints,3, 'S (turbulence PSD)',ErrStat, ErrMsg )
CALL CheckError()




! Calculate the single-point power spectral densities

CALL CalcTargetPSD(p, S, U, ErrStat, ErrMsg)


IF ( ALLOCATED( p%met%USR_Z         ) )  DEALLOCATE( p%met%USR_Z           )
IF ( ALLOCATED( p%met%USR_U         ) )  DEALLOCATE( p%met%USR_U           )
IF ( ALLOCATED( p%met%USR_WindDir   ) )  DEALLOCATE( p%met%USR_WindDir     )
IF ( ALLOCATED( p%met%USR_Sigma     ) )  DEALLOCATE( p%met%USR_Sigma       )
IF ( ALLOCATED( p%met%USR_L         ) )  DEALLOCATE( p%met%USR_L           )

IF ( ALLOCATED( p%met%ZL_profile    ) )  DEALLOCATE( p%met%ZL_profile      )
IF ( ALLOCATED( p%met%Ustar_profile ) )  DEALLOCATE( p%met%Ustar_profile   )

IF ( ALLOCATED( p%met%USR_Freq      ) )  DEALLOCATE( p%met%USR_Freq   )
IF ( ALLOCATED( p%met%USR_Uspec     ) )  DEALLOCATE( p%met%USR_Uspec  )
IF ( ALLOCATED( p%met%USR_Vspec     ) )  DEALLOCATE( p%met%USR_Vspec  )
IF ( ALLOCATED( p%met%USR_Wspec     ) )  DEALLOCATE( p%met%USR_Wspec  )


   ! Allocate memory for random number array

CALL AllocAry( PhaseAngles, p%grid%NPoints, p%grid%NumFreq, 3, 'Random Phases', ErrStat, ErrMsg )
CALL CheckError()


   ! Get the phase angles
CALL RndPhases(p%RNG, OtherSt_RandNum, PhaseAngles, p%grid%NPoints, p%grid%NumFreq, US, ErrStat, ErrMsg)
CALL CheckError()

IF (ALLOCATED(OtherSt_RandNum%nextSeed) ) DEALLOCATE(OtherSt_RandNum%nextSeed)  

CALL AllocAry( V, p%grid%NumSteps, p%grid%NPoints, 3, 'V (velocity)', ErrStat, ErrMsg) !  Allocate the array that contains the velocities.
CALL CheckError()


   ! Calculate the transfer function matrices from the spectral matrix (the fourier coefficients).

CALL WrScr ( ' Calculating the spectral and transfer function matrices:' )


IF (p%met%TurbModel_ID /= SpecModel_NONE) THEN                         ! MODIFIED BY Y GUO
    IF (p%met%TurbModel_ID == SpecModel_API) THEN
        CALL CalcFourierCoeffs_API( NSize, U, PhaseAngles, S, V, ErrStat, ErrMsg) 
    ELSE
        CALL CalcFourierCoeffs( NSize, U, PhaseAngles, S, V, ErrStat, ErrMsg)
    ENDIF
    CALL CheckError()
ENDIF


   ! Deallocate the Freq, S, and RandPhases arrays and the spectral matrix

IF ( ALLOCATED( p%grid%Freq ) )  DEALLOCATE( p%grid%Freq )
IF ( ALLOCATED( S           ) )  DEALLOCATE( S           )
IF ( ALLOCATED( PhaseAngles ) )  DEALLOCATE( PhaseAngles )


   ! Create the time series:
   
CALL Coeffs2TimeSeries( V, p%grid%NumSteps, p%grid%NPoints, ErrStat, ErrMsg)
CALL CheckError()


   ! Scale time series (if desired) for cross-component correlation or IEC statistics:
CALL TimeSeriesScaling(p, V, US)
CALL CheckError()



!..............................................................................
! Write hub-height output files (before adding mean and rotating final results)
!..............................................................................

IF ( p%WrFile(FileExt_HH) )  THEN
   CALL WrHH_ADtxtfile(V, p%RootName, p%IEC%TurbInt, ErrStat, ErrMsg)   
   CALL CheckError()
END IF

IF ( p%WrFile(FileExt_BIN) )  THEN
   CALL WrHH_binary(V, p%RootName, ErrStat, ErrMsg)
   CALL CheckError()
END IF

IF ( p%WrFile(FileExt_DAT) )  THEN
   CALL WrHH_text(V,  p%RootName, ErrStat, ErrMsg )   
   CALL CheckError()
END IF

   ! Write statistics of the run to the summary file:
CALL WrSum_Stats(V, USig, VSig, WSig, UXBar, UXSig, ErrStat, ErrMsg)
CALL CheckError()

!..............................................................................
! Add mean wind to u' components and rotate to inertial reference  
!  frame coordinate system
!..............................................................................
CALL AddMeanAndRotate(p, V, U, HWindDir)

TmpU = MAX( ABS(MAXVAL(U)-p%UHub), ABS(MINVAL(U)-p%UHub) )  !Get the range of wind speed values for scaling in BLADED-format .wnd files

   ! Deallocate memory for the matrix of the steady, u-component winds.

IF ( ALLOCATED( U        ) )  DEALLOCATE( U        )
IF ( ALLOCATED( HWindDir ) )  DEALLOCATE( HWindDir )

!..............................................................................
! Are we generating a coherent turbulent timestep file?
!..............................................................................
   
IF ( p%WrFile(FileExt_CTS) ) THEN
   p_CohStr%WSig=WSig

   CALL CohStr_WriteCTS(p%met, p%RNG, p%grid, p_CohStr, OtherSt_RandNum, y_CohStr, p%RootName, ErrStat, ErrMsg)
   CALL CheckError()
   
         
      ! Write the number of separate events to the summary file

   IF (p%met%KHtest) THEN
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

IF ( p%WrFile(FileExt_WND) .OR. p%WrFile(FileExt_BTS)  )  THEN 

      ! Calculate mean value & turb intensity of U-component of the interpolated hub point (for comparison w/ AeroDyn output)

   IF (p%grid%ExtraHubPT) THEN

         ! Get points for bi-linear interpolation
      ZLo_YLo   = ( NumGrid_Z2 - 1 )*p%grid%NumGrid_Y + NumGrid_Y2
      ZHi_YLo   = ( NumGrid_Z2     )*p%grid%NumGrid_Y + NumGrid_Y2
      ZLo_YHi   = ( NumGrid_Z2 - 1 )*p%grid%NumGrid_Y + NumGrid_Y2 + 1
      ZHi_YHi   = ( NumGrid_Z2     )*p%grid%NumGrid_Y + NumGrid_Y2 + 1
    
      TmpZ      = (p%grid%HubHt - p%grid%Z(NumGrid_Z2))/p%grid%GridRes_Z
      TmpY      = ( 0.0  - p%grid%Y(NumGrid_Y2))/p%grid%GridRes_Y
      CGridSum  = 0.0
      CGridSum2 = 0.0

      DO IT=1,p%grid%NumSteps
         
      ! Interpolate within the grid for this time step.

         Tmp_YL_Z  = ( V( IT, ZHi_YLo, 1 ) - V( IT, ZLo_YLo, 1 ) )*TmpZ + V( IT, ZLo_YLo, 1 )
         Tmp_YH_Z  = ( V( IT, ZHi_YHi, 1 ) - V( IT, ZLo_YHi, 1 ) )*TmpZ + V( IT, ZLo_YHi, 1 )
         TmpV      = ( Tmp_YH_Z - Tmp_YL_Z )*TmpY + Tmp_YL_Z

         CGridSum  = CGridSum  + TmpV
         CGridSum2 = CGridSum2 + TmpV*TmpV
      ENDDO ! IT

      UGridMean = CGridSum/p%grid%NumSteps
      UGridSig  = SQRT( ABS( (CGridSum2/p%grid%NumSteps) - UGridMean*UGridMean ) )
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

   IF ( p%WrFile(FileExt_BTS) ) THEN
      CALL WrBinTURBSIM(V, p%RootName, ErrStat, ErrMsg)
      CALL CheckError()
   END IF   
   
   IF ( p%WrFile(FileExt_WND) ) THEN

         ! We need to take into account the shear across the grid in the sigma calculations for scaling the data, 
         ! and ensure that 32.767*Usig >= |V-UHub| so that we don't get values out of the range of our scaling values
         ! in this BLADED-style binary output.  TmpU is |V-UHub|
      USig = MAX(USig,0.05*TmpU)
      CALL WrBinBLADED(V, USig, VSig, WSig, ErrStat, ErrMsg)
      CALL CheckError()
   ENDIF
   
ENDIF


   ! Are we generating FF formatted files?

IF ( p%WrFile(FileExt_UVW) )  THEN
   CALL WrFormattedFF(p%RootName, p%grid, p%UHub, V)
ENDIF ! ( WrFile(FileExt_UVW) )


   ! Deallocate the temporary V, Y, Z, and IYmax arrays.

IF ( ALLOCATED( V               ) )  DEALLOCATE( V               )
IF ( ALLOCATED( p%grid%Y        ) )  DEALLOCATE( p%grid%Y        )
IF ( ALLOCATED( p%grid%Z        ) )  DEALLOCATE( p%grid%Z        )
IF ( ALLOCATED( p%grid%IYmax    ) )  DEALLOCATE( p%grid%IYmax    )


IF (DEBUG_OUT) CLOSE( UD )       ! Close the debugging file

WRITE ( US, '(/"Nyquist frequency of turbulent wind field =      ",F8.3, " Hz")' ) 1.0 / (2.0 * p%grid%TimeStep)
IF ( p%WrFile(FileExt_CTS) .AND. y_CohStr%EventTimeStep > 0.0 ) THEN
   WRITE ( US, '( "Nyquist frequency of coherent turbulent events = ",F8.3, " Hz")' ) 1.0 / (2.0 * y_CohStr%EventTimeStep)
ENDIF


   ! Request CPU-time used.

CALL CPU_TIME ( CPUtime )

WRITE (US,"(//,'Processing complete.  ',A,' CPU seconds used.')")  TRIM( Num2LStr( CPUtime ) )

!BONNIE: ****************************
!Time = TIMEF()
!PRINT *, Time
!WRITE (US,"(//,'BONNIE TEST.  ',A,' seconds for completion.')")  TRIM( Flt2LStr( CPUtime ) )
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
