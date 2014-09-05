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

TYPE(CohStr_ParameterType)       :: p_CohStr
TYPE(RandNum_OtherStateType)     :: OtherSt_RandNum             ! other states for random numbers (next seed, etc)
TYPE(CohStr_OutputType)          :: y_CohStr
TYPE( TurbSim_ParameterType )    :: p


REAL(ReKi)              ::  USig                            ! Standard deviation of the u-component wind speed at the hub
REAL(ReKi)              ::  VSig                            ! Standard deviation of the v-component wind speed at the hub
REAL(ReKi)              ::  WSig                            ! Standard deviation of the w-component wind speed at the hub
REAL(ReKi)              ::  UXSig                           ! Standard deviation of the U-component wind speed at the hub

REAL(DbKi)              ::  UXBar                           ! The mean U-component (u rotated; x-direction) wind speed at the hub
REAL(ReKi)              ::  CPUtime                         ! Contains the number of seconds since the start of the program

#ifdef DEBUG_TS
INTEGER                 ::  IFreq ! for debugging                              
#endif

INTEGER                 ::  IY  , I                         ! An index for the Y position of a point I

INTEGER                 ::  UnOut                           ! unit for output files


CHARACTER(200)          :: InFile                           ! Name of the TurbSim input file.


REAL(ReKi), ALLOCATABLE          :: PhaseAngles (:,:,:)                      ! The array that holds the random phases [number of points, number of frequencies, number of wind components=3].
REAL(ReKi), ALLOCATABLE          :: S           (:,:,:)                      ! The turbulence PSD array (NumFreq,NPoints,3).
REAL(ReKi), ALLOCATABLE          :: V           (:,:,:)                      ! An array containing the summations of the rows of H (NumSteps,NPoints,3).
REAL(ReKi), ALLOCATABLE          :: U           (:)                          ! The steady u-component wind speeds for the grid (NPOints).
REAL(ReKi), ALLOCATABLE          :: HWindDir    (:)                          ! A profile of horizontal wind angle (measure of wind direction with height)

!TYPE( TurbSim_ParameterType )    :: p

INTEGER(IntKi)          :: ErrStat     ! allocation status
CHARACTER(MaxMsgLen)    :: ErrMsg      ! error message


!BONNIE:*****************************
!    Time = TIMEF() ! Initialize the Wall Clock Time counter
!BONNIE:*****************************   

p%US = -1

   ! ... Initialize NWTC Library (open console, set pi constants) ...
CALL NWTC_Init( ProgNameIN=TurbSim_Ver%Name, EchoLibVer=.FALSE. )       


   ! Print out program name, version, and date.

CALL DispNVD(TurbSim_Ver)


   ! Check for command line arguments.
InFile = 'TurbSim.inp'  ! default name for input file
CALL CheckArgs( InFile )

CALL GetRoot( InFile, p%RootName )

   ! Open input file and summary file.

CALL OpenSummaryFile( p%RootName, p%US, p%DescStr, ErrStat, ErrMsg )
CALL CheckError()

   ! Get input parameters.

CALL ReadInputFile(InFile, p, p_cohStr, OtherSt_RandNum, ErrStat, ErrMsg)
CALL CheckError()

CALL WrSum_EchoInputs(p, p_CohStr) 
call WrSum_UserInput(p%met,p%US)

call TS_ValidateInput(p, ErrStat, ErrMsg)
CALL CheckError()


   ! Open appropriate output files.  We will open formatted FF files later, if requested.
   ! Mention the files in the summary file.

IF ( ANY (p%WrFile) )  THEN
   CALL GetNewUnit( UnOut )

   WRITE (p%US,"( // 'You have requested that the following file(s) be generated:' / )")
   CALL WrScr1  ( ' You have requested that the following file(s) be generated:' )

   IF ( p%WrFile(FileExt_BIN) )  THEN   

!      CALL OpenBOutFile ( UnOut, TRIM( p%RootName)//'.bin', ErrStat, ErrMsg )
      CALL OpenUOutfile ( UnOut , TRIM( p%RootName)//'.bin', ErrStat, ErrMsg )  ! just making sure it can be opened (not locked elsewhere)
      CLOSE(UnOut)
      CALL CheckError()
      
      WRITE (p%US,"( 3X ,'"//TRIM( p%RootName)//".bin (a binary hub-height turbulence-parameter file)' )")  
      CALL WrScr ( '    '//TRIM( p%RootName)//'.bin (a binary hub-height turbulence-parameter file)' )

   ENDIF

   IF ( p%WrFile(FileExt_DAT) )  THEN     

      CALL OpenFOutFile ( UnOut, TRIM( p%RootName)//'.dat', ErrStat, ErrMsg ) ! just making sure it can be opened (not locked elsewhere)
      CLOSE( UnOut )
      CALL CheckError()
      
      WRITE (p%US, "( 3X ,'"//TRIM( p%RootName)//".dat (a formatted turbulence-parameter file)' )")  
      CALL WrScr ( '     '//TRIM( p%RootName)//'.dat (a formatted turbulence-parameter file)' )

   ENDIF

   IF ( p%WrFile(FileExt_HH) )  THEN     

      CALL OpenFOutFile ( UnOut, TRIM( p%RootName)//'.hh', ErrStat, ErrMsg ) ! just making sure it can be opened (not locked elsewhere)
      CLOSE( UnOut )
      CALL CheckError()
      
      WRITE (p%US,"( 3X ,'"//TRIM( p%RootName)//".hh  (an AeroDyn hub-height file)' )")
      CALL WrScr ( '    '//TRIM( p%RootName)//'.hh  (an AeroDyn hub-height file)' )

   ENDIF

   IF ( p%WrFile(FileExt_BTS) )  THEN

      CALL OpenBOutFile ( UnOut, TRIM(p%RootName)//'.bts', ErrStat, ErrMsg )
      CLOSE( UnOut )
      CALL CheckError()

      WRITE (p%US,"( 3X ,'"//TRIM( p%RootName)//".bts (an AeroDyn/TurbSim full-field wind file)' )")  
      CALL WrScr ( '    '//TRIM( p%RootName)//'.bts (an AeroDyn/TurbSim full-field wind file)' )

   ENDIF

   IF ( p%WrFile(FileExt_WND) )  THEN

      CALL OpenBOutFile ( UnOut, TRIM(p%RootName)//'.wnd', ErrStat, ErrMsg )
      CLOSE(UnOut)
      CALL CheckError()

      WRITE (p%US,"( 3X ,'"//TRIM( p%RootName)//".wnd (an AeroDyn/BLADED full-field wnd file)' )")  
      CALL WrScr ( '    '//TRIM( p%RootName)//'.wnd (an AeroDyn/BLADED full-field wnd file)' )

   ENDIF
   
   IF ( p%WrFile(FileExt_TWR) .AND. p%WrFile(FileExt_WND) )  THEN

      CALL OpenBOutFile ( UnOut, TRIM( p%RootName )//'.twr', ErrStat, ErrMsg )
      CLOSE(UnOut)       
      CALL CheckError()

      WRITE (p%US,"( 3X ,'"//TRIM( p%RootName)//".twr (a binary tower file)' )")  
      CALL WrScr ( '    '//TRIM( p%RootName)//'.twr (a binary tower file)' )

   ENDIF

   IF ( p%WrFile(FileExt_CTS) ) THEN      
      WRITE (p%US,"( 3X ,'"//TRIM( p%RootName)//".cts (a coherent turbulence time step file)' )")  
      CALL WrScr ( '    '//TRIM( p%RootName)//'.cts (a coherent turbulence time step file)' )
   ENDIF

   IF ( p%WrFile(FileExt_UVW) )  THEN
      WRITE (p%US,"( 3X ,'"//TRIM( p%RootName)//".u (a formatted full-field U-component file)' )")  
      CALL WrScr ( '    '//TRIM( p%RootName)//'.u (a formatted full-field U-component file)' )

      WRITE (p%US,"( 3X ,'"//TRIM( p%RootName)//".v (a formatted full-field V-component file)' )")  
      CALL WrScr ( '    '//TRIM( p%RootName)//'.v (a formatted full-field V-component file)' )

      WRITE (p%US,"( 3X ,'"//TRIM( p%RootName)//".w (a formatted full-field W-component file)' )")  
      CALL WrScr ( '    '//TRIM( p%RootName)//'.w (a formatted full-field W-component file)' )
   ENDIF

ELSE
   ErrStat = ErrID_Fatal
   ErrMsg  = 'You have requested no output.'
   CALL CheckError()
ENDIF


!..................................................................................................................................
! Define the spatial grid
!..................................................................................................................................

   ! Define the other parameters for the time series.
CALL CreateGrid( p%grid, p%usr, p%UHub, p%WrFile(FileExt_TWR), ErrStat, ErrMsg )
CALL CheckError()
      
CALL CalcIECScalingParams(p%IEC, p%grid%HubHt, p%UHub, p%met%InCDec, p%met%InCohB, p%met%TurbModel_ID, p%met%IsIECModel)                  


!..................................................................................................................................
! Calculate mean velocity and direction profiles:
!..................................................................................................................................

   !  Wind speed:
CALL AllocAry(U,     SIZE(p%grid%Z), 'u (steady, u-component winds)', ErrStat, ErrMsg )
CALL CheckError()

IF ( p%met%WindProfileType(1:3) == 'API' )  THEN
   CALL getVelocityProfile( p, p%met%URef, p%met%RefHt,  p%grid%Z, U, ErrStat, ErrMsg)
ELSE 
   CALL getVelocityProfile( p, p%UHub,     p%grid%HubHt, p%grid%Z, U, ErrStat, ErrMsg) 
ENDIF
CALL CheckError()

   ! Wind Direction:
CALL AllocAry(HWindDir, SIZE(p%grid%Z), 'HWindDir (wind direction profile)', ErrStat, ErrMsg )                  ! Allocate the array for the wind direction profile      
CALL CheckError()
   
CALL getDirectionProfile(p, p%grid%Z, HWindDir, ErrStat, ErrMsg)
CALL CheckError()
   
p%met%HH_HFlowAng = HWindDir( p%grid%HubIndx )

!..................................................................................................................................
! Calculate remaining parameters required for simulation:
!..................................................................................................................................


IF ( p%met%TurbModel_ID == SpecModel_GP_LLJ) THEN

      ! Allocate the arrays for the z/l and ustar profile
      
   CALL AllocAry(p%met%ZL_profile,    SIZE(p%grid%Z), 'ZL_profile (z/l profile)', ErrStat, ErrMsg )
   CALL CheckError()
   CALL AllocAry(p%met%Ustar_profile, SIZE(p%grid%Z), 'Ustar_profile (friction velocity profile)', ErrStat, ErrMsg )         
   CALL CheckError()

   p%met%ZL_profile(:)    = getZLARY(    U, p%grid%Z, p%met%Rich_No, p%met%ZL, p%met%L, p%met%ZLOffset, p%met%WindProfileType )
   p%met%Ustar_profile(:) = getUstarARY( p, U, p%grid%Z, p%met%UStarOffset, p%met%UStarSlope )
   
END IF

           
CALL WrSum_SpecModel( p, U, HWindDir, ErrStat, ErrMsg )
CALL CheckError()

IF (DEBUG_OUT) THEN
   ! If we are going to generate debug output, open the debug file.
   !UD = -1
   !CALL GetNewUnit( UD, ErrStat, ErrMsg )
   CALL OpenFOutFile ( UD, TRIM( p%RootName)//'.dbg', ErrStat, ErrMsg )
      CALL CheckError()      
   
   !bjj: This doesn't need to be part of the summary file, so put it in the debug file
   IF (p%met%WindProfileType(1:1) == 'J' ) THEN
      WRITE(UD,"(//'Jet wind profile Chebyshev coefficients:' / )") 
      WRITE(UD,"( 3X,'Order: ',11(1X,I10) )")  ( I, I=0,10 )
      WRITE(UD,"( 3X,'------ ',11(1X,'----------') )")
      WRITE(UD,"( 3X,'Speed: ',11(1X,E10.3) )")  ( p%met%ChebyCoef_WS(I), I=1,11 )
      WRITE(UD,"( 3X,'Angle: ',11(1X,E10.3) )")  ( p%met%ChebyCoef_WD(I), I=1,11 )   
   ENDIF
ENDIF

!..................................................................................................................................
! Get the single-point power spectral densities
!..................................................................................................................................

CALL AllocAry( S,    p%grid%NumFreq,p%grid%NPoints,3, 'S (turbulence PSD)',ErrStat, ErrMsg )
CALL CheckError()

CALL CalcTargetPSD(p, S, U, ErrStat, ErrMsg)
CALL CheckError()


IF ( ALLOCATED( p%met%USR_Z         ) )  DEALLOCATE( p%met%USR_Z           )
IF ( ALLOCATED( p%met%USR_U         ) )  DEALLOCATE( p%met%USR_U           )
IF ( ALLOCATED( p%met%USR_WindDir   ) )  DEALLOCATE( p%met%USR_WindDir     )
IF ( ALLOCATED( p%met%USR_Sigma     ) )  DEALLOCATE( p%met%USR_Sigma       )
IF ( ALLOCATED( p%met%USR_L         ) )  DEALLOCATE( p%met%USR_L           )

IF ( ALLOCATED( p%met%ZL_profile    ) )  DEALLOCATE( p%met%ZL_profile      )
IF ( ALLOCATED( p%met%Ustar_profile ) )  DEALLOCATE( p%met%Ustar_profile   )

!IF ( ALLOCATED( p%usr%f            ) )  DEALLOCATE( p%usr%f               ) bjj: do we need to keep these for phase angles?
IF ( ALLOCATED( p%usr%S             ) )  DEALLOCATE( p%usr%S               )

!..................................................................................................................................
! Get the phase angles
!..................................................................................................................................
   
CALL AllocAry( PhaseAngles, p%grid%NPoints, p%grid%NumFreq, 3, 'Random Phases', ErrStat, ErrMsg )
CALL CheckError()
        
CALL SetPhaseAngles( p, OtherSt_RandNum, PhaseAngles, ErrStat, ErrMsg )
CALL CheckError()


#ifdef DEBUG_TS
DO iFreq=1, p%grid%NumFreq
   WRITE( 73, '(7(F15.6," "))') iFreq*(1.0/p%grid%AnalysisTime), ( S(iFreq,1,iY), phaseAngles(1,iFreq,iY), iY=1,3 )
END DO
#endif


IF ( ALLOCATED(OtherSt_RandNum%nextSeed ) ) DEALLOCATE( OtherSt_RandNum%nextSeed )  
IF ( ALLOCATED(p%usr%PhaseAngles        ) ) DEALLOCATE( p%usr%PhaseAngles        )  
IF ( ALLOCATED(p%usr%f                  ) ) DEALLOCATE( p%usr%f                  ) ! bjj: do we need to keep these for phase angles or should we destroy earlier?

!..................................................................................................................................
! Get the Fourier Coefficients
!..................................................................................................................................
CALL AllocAry( V, p%grid%NumSteps, p%grid%NPoints, 3, 'V (velocity)', ErrStat, ErrMsg) !  Allocate the array that contains the velocities.
CALL CheckError()


   ! Calculate the transfer function matrices from the spectral matrix (the fourier coefficients).

CALL WrScr ( ' Calculating the spectral and transfer function matrices:' )


IF (p%met%TurbModel_ID /= SpecModel_NONE) THEN                         ! MODIFIED BY Y GUO
    IF (p%met%TurbModel_ID == SpecModel_API) THEN
        CALL CalcFourierCoeffs_API( p, U, PhaseAngles, S, V, ErrStat, ErrMsg) 
    ELSE
        CALL CalcFourierCoeffs( p, U, PhaseAngles, S, V, ErrStat, ErrMsg)
    ENDIF
    CALL CheckError()
ELSE
   V = 0.0_ReKi
ENDIF


   ! Deallocate the Freq, S, and RandPhases arrays and the spectral matrix

IF ( ALLOCATED( p%grid%Freq ) )  DEALLOCATE( p%grid%Freq )
IF ( ALLOCATED( S           ) )  DEALLOCATE( S           )
IF ( ALLOCATED( PhaseAngles ) )  DEALLOCATE( PhaseAngles )

!..................................................................................................................................
! Create the time series
!..................................................................................................................................  
CALL Coeffs2TimeSeries( V, p%grid%NumSteps, p%grid%NPoints, ErrStat, ErrMsg)
CALL CheckError()

!..................................................................................................................................
! Scale time series (if desired) for cross-component correlation or IEC statistics:
!..................................................................................................................................  
CALL ScaleTimeSeries(p, V, ErrStat, ErrMsg)
CALL CheckError()


!..................................................................................................................................
! Write statistics of the run to the summary file:
!..................................................................................................................................
CALL WrSum_Stats(p, V, USig, VSig, WSig, UXBar, UXSig, ErrStat, ErrMsg)
CALL CheckError()

!..................................................................................................................................
! Write hub-height output files (before adding mean and rotating final results)
!..................................................................................................................................

IF ( p%WrFile(FileExt_HH) )  THEN
   CALL WrHH_ADtxtfile(p, V, p%IEC%TurbInt, ErrStat, ErrMsg)   
   CALL CheckError()
END IF

IF ( p%WrFile(FileExt_BIN) )  THEN
   CALL WrHH_binary(p, V, ErrStat, ErrMsg)
   CALL CheckError()
END IF

IF ( p%WrFile(FileExt_DAT) )  THEN
   CALL WrHH_text(p, V, ErrStat, ErrMsg )   
   CALL CheckError()
END IF

!..................................................................................................................................
! Add mean wind to u' components and rotate to inertial reference  
!  frame coordinate system
!..................................................................................................................................
CALL AddMeanAndRotate(p, V, U, HWindDir)


   ! Deallocate memory for the matrix of the steady, u-component winds.

IF ( ALLOCATED( U        ) )  DEALLOCATE( U        )
IF ( ALLOCATED( HWindDir ) )  DEALLOCATE( HWindDir )

!..................................................................................................................................
! Generate coherent turbulence if desired:
!..................................................................................................................................   
IF ( p%WrFile(FileExt_CTS) ) THEN
   p_CohStr%WSig=WSig

   CALL CohStr_WriteCTS(p, p_CohStr, OtherSt_RandNum, y_CohStr, ErrStat, ErrMsg)
   CALL CheckError()
   
         
      ! Write the number of separate events to the summary file

   IF (p%met%KHtest) THEN
      WRITE ( p%US,'(/)' )
   ELSE
      WRITE ( p%US,'(//A,F8.3," seconds")' ) 'Average expected time between events = ',y_CohStr%lambda
   ENDIF

   WRITE ( p%US, '(A,I8)'   )            'Number of coherent events            = ', y_CohStr%NumCTEvents_separate
   WRITE ( p%US, '(A,F8.3," seconds")')  'Predicted length of coherent events  = ', y_CohStr%ExpectedTime
   WRITE ( p%US, '(A,F8.3," seconds")')  'Length of coherent events            = ', y_CohStr%EventTimeSum
   WRITE ( p%US, '(A,F8.3," (m/s)^2")')  'Maximum predicted event CTKE         = ', y_CohStr%CTKE
   
ENDIF !WrACT

!..................................................................................................................................
! Generate output files:
!..................................................................................................................................

   ! Are we generating FF files?
IF ( p%WrFile(FileExt_BTS) .OR. p%WrFile(FileExt_WND) ) THEN
   CALL WrSum_InterpolatedHubStats(p, V)

      
   IF ( p%WrFile(FileExt_BTS) ) THEN
      CALL WrBinTURBSIM(p, V, ErrStat, ErrMsg)
      CALL CheckError()
   END IF   
   
   IF ( p%WrFile(FileExt_WND) ) THEN      
      CALL WrBinBLADED(p, V, USig, VSig, WSig, ErrStat, ErrMsg)
      CALL CheckError()
   END IF
END IF


   ! Are we generating FF formatted files?

IF ( p%WrFile(FileExt_UVW) )  THEN
   CALL WrFormattedFF(p%RootName, p%grid, p%UHub, V)
ENDIF ! ( WrFile(FileExt_UVW) )



!..................................................................................................................................
! End:
!..................................................................................................................................

   ! Deallocate the V, Y, Z, and IYmax arrays.

IF ( ALLOCATED( V               ) )  DEALLOCATE( V               )
IF ( ALLOCATED( p%grid%Y        ) )  DEALLOCATE( p%grid%Y        )
IF ( ALLOCATED( p%grid%Z        ) )  DEALLOCATE( p%grid%Z        )


IF (DEBUG_OUT) CLOSE( UD )       ! Close the debugging file

WRITE ( p%US, '(/"Nyquist frequency of turbulent wind field =      ",F8.3, " Hz")' ) 1.0 / (2.0 * p%grid%TimeStep)
IF ( p%WrFile(FileExt_CTS) .AND. y_CohStr%EventTimeStep > 0.0 ) THEN
   WRITE ( p%US, '( "Nyquist frequency of coherent turbulent events = ",F8.3, " Hz")' ) 1.0 / (2.0 * y_CohStr%EventTimeStep)
ENDIF


   ! Request CPU-time used.

CALL CPU_TIME ( CPUtime )

WRITE (p%US,"(//,'Processing complete.  ',A,' CPU seconds used.')")  TRIM( Num2LStr( CPUtime ) )

!BONNIE: ****************************
!Time = TIMEF()
!PRINT *, Time
!WRITE (US,"(//,'BONNIE TEST.  ',A,' seconds for completion.')")  TRIM( Flt2LStr( CPUtime ) )
!END BONNIE ************************************

CLOSE ( p%US )
p%US = -1

CALL WrScr1  ( ' Processing complete.  '//TRIM( Num2LStr( CPUtime ) )//' CPU seconds used.' )
CALL NormStop


CONTAINS
!..................................................................................................................................
SUBROUTINE CheckError()
   IF (ErrStat /= ErrID_None) THEN
   
      IF (ErrStat >= AbortErrLev) THEN
         CALL TS_end(p)
         
         WRITE (p%US, "(/'ERROR:  ', A / )") TRIM(ErrMSg)
         WRITE (p%US, "('ABORTING PROGRAM.')" )
         
         CALL ProgAbort ( TRIM(ErrMSg), .FALSE., 5.0_ReKi )
         
      ELSE
         CALL WrScr(TRIM(ErrMsg))
      END IF
   END IF
END SUBROUTINE CheckError

END PROGRAM
