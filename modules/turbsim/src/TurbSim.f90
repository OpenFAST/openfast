!=======================================================================
PROGRAM TurbSim

   ! A turbulence simulator developed at the
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
   ! v2.00.00a-bjj      Oct-2014  B. Jonkman  (NWTC Library 2.0)
   !
   ! This program simulates a full field of turbulent winds at points in space
   ! in a rectangular Cartesian plane perpendicular to the mean wind direction.
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
   !               Yb(1)                                    Yb(2)                        Yb(3)                         ... Yb(NumGrid_Y)
   !              --------------------------------------------------------------------------------------------------------------------------------------------
   ! Zb(NumGrid_Z):|V(GridPtIndx(NumGrid_Y*(NumGrid_Z-1)+1)) ...                                                           V(GridPtIndx(NumGrid_Z*NumGrid_Y)) |
   ! ...           |...                                                                                                                                       |
   ! Zb(2)        :|V(GridPtIndx(NumGrid_Y +             1)) V(GridPtIndx(NumGrid_Y + 2)) V(GridPtIndx(NumGrid_Y + 3)) ... V(GridPtIndx(        2*NumGrid_Y)) |
   ! Zb(1)        :|V(GridPtIndx(                        1)) V(GridPtIndx(            2)) V(GridPtIndx(            3)) ... V(GridPtIndx(          NumGrid_Y)) |
   !              --------------------------------------------------------------------------------------------------------------------------------------------
   ! 
   ! Zb(i) < Zb(i+1) for all integers i, 1 <= i < NumGrid_Z
   ! Yb(j) < Yb(j+1) for all integers j, 1 <= j < NumGrid_Y
   ! note that the Y and Z arrays used in the code are NOT necessarially in the order of Yb and Zb described here.
   !
   ! If an extra hub point is necessary because the point does not fall on the grid, 
   ! then it is added immediately following the regular grid points.
   !
   ! If the tower wind file output is selected, those extra points (in a single vertical 
   ! line) are added at the end, at the end (after the grid and hub point).
   !
   ! Any user-defined time-series points are stored at the BEGINNING, before the grid points.
   ! --------------------------------------------------------------------------------------------------------


USE TSsubs
USE TS_FileIO
USE TS_Profiles
use TS_CohStructures
use VersionInfo

IMPLICIT                   NONE


   ! Declare local variables

TYPE( TurbSim_ParameterType )    :: p                       ! TurbSim parameters
TYPE(RandNum_OtherStateType)     :: OtherSt_RandNum         ! other states for random numbers (next seed, etc)

REAL(ReKi), ALLOCATABLE          :: PhaseAngles (:,:,:)     ! The array that holds the random phases [number of points, number of frequencies, number of wind components=3].
REAL(ReKi), ALLOCATABLE          :: S           (:,:,:)     ! The turbulence PSD array (NumFreq,NPoints,3).
REAL(ReKi), ALLOCATABLE          :: V           (:,:,:)     ! An array containing the summations of the rows of H (NumSteps,NPoints,3).
REAL(ReKi), ALLOCATABLE          :: U           (:)         ! The steady u-component wind speeds for the grid (NPoints).
REAL(ReKi), ALLOCATABLE          :: HWindDir    (:)         ! A profile of horizontal wind angle (NPoints) (measure of wind direction with height)
REAL(ReKi), ALLOCATABLE          :: VWindDir    (:)         ! A profile of vretical wind angle (NPoints) (measure of wind vertical angle with height)

REAL(ReKi)                       :: CPUtime                 ! Contains the number of seconds since the start of the program

REAL(ReKi)                       ::  USig                   ! Standard deviation of the u-component wind speed at the hub (used for scaling WND files)
REAL(ReKi)                       ::  VSig                   ! Standard deviation of the v-component wind speed at the hub (used for scaling WND files)
REAL(ReKi)                       ::  WSig                   ! Standard deviation of the w-component wind speed at the hub (used for scaling WND files and CTS files)

INTEGER(IntKi)                   :: ErrStat                 ! allocation status
CHARACTER(MaxMsgLen)             :: ErrMsg                  ! error message
CHARACTER(200)                   :: InFile                  ! Name of the TurbSim input file.
CHARACTER(20)                    :: FlagArg                 ! flag argument from command line


!BONNIE:*****************************
!    Time = TIMEF() ! Initialize the Wall Clock Time counter
!BONNIE:*****************************   

p%US = -1

   ! ... Initialize NWTC Library (open console, set pi constants) ...
CALL NWTC_Init( ProgNameIN=TurbSim_Ver%Name, EchoLibVer=.FALSE. )       

   ! Check for command line arguments.
InFile = 'TurbSim.inp'  ! default name for input file
CALL CheckArgs( InFile, Flag=FlagArg )
IF ( LEN( TRIM(FlagArg) ) > 0 ) CALL NormStop()

   ! Print out program name, version, and date.

   ! Display the copyright notice
   CALL DispCopyrightLicense( TurbSim_Ver%Name )
      ! Obtain OpenFAST git commit hash
   TurbSim_Ver%Ver = 'from OpenFAST-'//TRIM(QueryGitVersion())
      ! Tell our users what they're running
   CALL WrScr( ' Running '//TRIM( GetNVD(TurbSim_Ver) )//NewLine )
   
CALL GetRoot( InFile, p%RootName )

   ! Open input file and summary file.

CALL OpenSummaryFile( p%RootName, p%US, p%DescStr, ErrStat, ErrMsg )
CALL CheckError(ErrStat, ErrMsg)

   ! Get input parameters.

CALL ReadInputFile(InFile, p, OtherSt_RandNum, ErrStat, ErrMsg)
CALL CheckError(ErrStat, ErrMsg)

CALL WrSum_EchoInputs(p) 
call WrSum_UserInput(p%met,p%usr, p%US)

CALL TS_ValidateInput(p, ErrStat, ErrMsg)
CALL CheckError(ErrStat, ErrMsg)


!..................................................................................................................................
! Define the spatial grid
!..................................................................................................................................

   ! Define the other parameters for the time series.
CALL CreateGrid( p%grid, p%usr, p%UHub, p%WrFile(FileExt_TWR), ErrStat, ErrMsg )
CALL CheckError(ErrStat, ErrMsg)
      
!..................................................................................................................................
! Calculate mean velocity and direction profiles:
!..................................................................................................................................

   !  Wind speed:
CALL AllocAry(U,     SIZE(p%grid%Z), 'u (steady, u-component winds)', ErrStat, ErrMsg )
CALL CheckError(ErrStat, ErrMsg)

CALL getVelocityProfile( p, p%UHub,     p%grid%HubHt, p%grid%Z, U, ErrStat, ErrMsg) 
CALL CheckError(ErrStat, ErrMsg)

   ! Wind Direction:
CALL AllocAry(HWindDir, SIZE(p%grid%Z), 'HWindDir (wind direction profile)', ErrStat, ErrMsg )                  ! Allocate the array for the wind direction profile      
CALL CheckError(ErrStat, ErrMsg)

CALL AllocAry(VWindDir, SIZE(p%grid%Z), 'VWindDir (vertical wind angle profile)', ErrStat, ErrMsg )             ! Allocate the array for the vertical wind profile      
CALL CheckError(ErrStat, ErrMsg)

CALL getDirectionProfile(p, p%grid%Z, HWindDir, VWindDir, ErrStat, ErrMsg)
CALL CheckError(ErrStat, ErrMsg)
   
p%met%HH_HFlowAng = HWindDir( p%grid%HubIndx )
p%met%HH_VFlowAng = VWindDir( p%grid%HubIndx )

!..................................................................................................................................
! Calculate remaining parameters required for simulation:
!..................................................................................................................................


IF ( p%met%TurbModel_ID == SpecModel_GP_LLJ) THEN

      ! Allocate the arrays for the z/l and ustar profile
      
   CALL AllocAry(p%met%ZL_profile,    SIZE(p%grid%Z), 'ZL_profile (z/l profile)', ErrStat, ErrMsg )
   CALL CheckError(ErrStat, ErrMsg)
   CALL AllocAry(p%met%Ustar_profile, SIZE(p%grid%Z), 'Ustar_profile (friction velocity profile)', ErrStat, ErrMsg )         
   CALL CheckError(ErrStat, ErrMsg)

   p%met%ZL_profile(:)    = getZLProfile(       U, p%grid%Z, p%met%Rich_No, p%met%ZL, p%met%L, p%met%ZLOffset, p%met%WindProfileType )
   p%met%Ustar_profile(:) = getUStarProfile( p, U, p%grid%Z, p%met%UStarOffset, p%met%UStarSlope )
   
END IF

           
CALL WrSum_SpecModel( p, U, HWindDir, VWindDir, ErrStat, ErrMsg )
CALL CheckError(ErrStat, ErrMsg)


!..................................................................................................................................
! Get the single-point power spectral densities
!..................................................................................................................................

CALL AllocAry( S,    p%grid%NumFreq,p%grid%NPoints,3, 'S (turbulence PSD)',ErrStat, ErrMsg )
CALL CheckError(ErrStat, ErrMsg)

CALL CalcTargetPSD(p, S, U, ErrStat, ErrMsg)
CALL CheckError(ErrStat, ErrMsg)

   ! we don't need these arrays any more, so deallocate to save some space
IF ( ALLOCATED( p%met%USR_Z         ) )  DEALLOCATE( p%met%USR_Z           )
IF ( ALLOCATED( p%met%USR_U         ) )  DEALLOCATE( p%met%USR_U           )
IF ( ALLOCATED( p%met%USR_WindDir   ) )  DEALLOCATE( p%met%USR_WindDir     )
IF ( ALLOCATED( p%met%USR_Sigma     ) )  DEALLOCATE( p%met%USR_Sigma       )
IF ( ALLOCATED( p%met%USR_L         ) )  DEALLOCATE( p%met%USR_L           )

IF ( ALLOCATED( p%met%ZL_profile    ) )  DEALLOCATE( p%met%ZL_profile      )
IF ( ALLOCATED( p%met%Ustar_profile ) )  DEALLOCATE( p%met%Ustar_profile   )

!IF ( ALLOCATED( p%usr%f            ) )  DEALLOCATE( p%usr%f               ) bjj: do we need to keep these for phase angles?
IF ( ALLOCATED( p%usr%S             ) )  DEALLOCATE( p%usr%S               )
IF ( ALLOCATED( p%usr%meanU         ) )  DEALLOCATE( p%usr%meanU           )
IF ( ALLOCATED( p%usr%meanDir       ) )  DEALLOCATE( p%usr%meanDir         )
IF ( ALLOCATED( p%usr%meanVAng      ) )  DEALLOCATE( p%usr%meanVAng        )

!..................................................................................................................................
! Get the phase angles
!..................................................................................................................................
   
CALL AllocAry( PhaseAngles, p%grid%NPoints, p%grid%NumFreq, 3, 'Random Phases', ErrStat, ErrMsg )
CALL CheckError(ErrStat, ErrMsg)
        
CALL SetPhaseAngles( p, OtherSt_RandNum, PhaseAngles, ErrStat, ErrMsg )
CALL CheckError(ErrStat, ErrMsg)

   ! we don't need these arrays any more, so deallocate to save some space
IF ( ALLOCATED(OtherSt_RandNum%nextSeed ) ) DEALLOCATE( OtherSt_RandNum%nextSeed )  
IF ( ALLOCATED(p%usr%PhaseAngles        ) ) DEALLOCATE( p%usr%PhaseAngles        )  
IF ( ALLOCATED(p%usr%f                  ) ) DEALLOCATE( p%usr%f                  ) ! bjj: do we need to keep these for phase angles or should we destroy earlier?

!..................................................................................................................................
! Get the Fourier Coefficients
!..................................................................................................................................
CALL AllocAry( V, p%grid%NumSteps, p%grid%NPoints, 3, 'V (velocity)', ErrStat, ErrMsg) !  Allocate the array that contains the velocities.
CALL CheckError(ErrStat, ErrMsg)


   ! Calculate the transfer function matrices from the spectral matrix (the fourier coefficients).

CALL WrScr ( ' Calculating the spectral and transfer function matrices:' )

CALL CalcFourierCoeffs(  p, U, PhaseAngles, S, V, ErrStat, ErrMsg )
CALL CheckError(ErrStat, ErrMsg)


   ! we don't need these arrays any more, so deallocate to save some space
IF ( ALLOCATED( p%grid%Freq ) )  DEALLOCATE( p%grid%Freq )
IF ( ALLOCATED( S           ) )  DEALLOCATE( S           )
IF ( ALLOCATED( PhaseAngles ) )  DEALLOCATE( PhaseAngles )

!..................................................................................................................................
! Create the time series
!..................................................................................................................................  
CALL Coeffs2TimeSeries( V, p%grid%NumSteps, p%grid%NPoints, p%usr%NPoints, ErrStat, ErrMsg)
CALL CheckError(ErrStat, ErrMsg)

!..................................................................................................................................
! Scale time series (if desired) for cross-component correlation or IEC statistics:
!..................................................................................................................................  
CALL ScaleTimeSeries(p, V, ErrStat, ErrMsg)
CALL CheckError(ErrStat, ErrMsg)


!..................................................................................................................................
! Write statistics of the run to the summary file:
!..................................................................................................................................
CALL WrSum_Stats(p, V, USig, VSig, WSig, ErrStat, ErrMsg)
CALL CheckError(ErrStat, ErrMsg)

!..................................................................................................................................
! Write hub-height output files (before adding mean and rotating final results)
!..................................................................................................................................

IF ( p%WrFile(FileExt_HH) )  THEN
   CALL WrHH_ADtxtfile(p, V, p%IEC%TurbInt, ErrStat, ErrMsg)   
   CALL CheckError(ErrStat, ErrMsg)
END IF

IF ( p%WrFile(FileExt_BIN) )  THEN
   CALL WrHH_binary(p, V, ErrStat, ErrMsg)
   CALL CheckError(ErrStat, ErrMsg)
END IF

IF ( p%WrFile(FileExt_DAT) )  THEN
   CALL WrHH_text(p, V, ErrStat, ErrMsg )   
   CALL CheckError(ErrStat, ErrMsg)
END IF


!..................................................................................................................................
! Add mean wind to u' components and rotate to inertial reference frame coordinate system
!..................................................................................................................................
CALL AddMeanAndRotate(p, V, U, HWindDir, VWindDir)

   ! Deallocate memory for the matrix of the steady, u-component winds.

IF ( ALLOCATED( U              ) )  DEALLOCATE( U              )
IF ( ALLOCATED( HWindDir       ) )  DEALLOCATE( HWindDir       )
IF ( ALLOCATED( VWindDir       ) )  DEALLOCATE( VWindDir       )

!..................................................................................................................................
! Generate coherent turbulence if desired:
!..................................................................................................................................   
IF ( p%WrFile(FileExt_CTS) ) THEN
   
   CALL CohStr_WriteCTS(p, WSig, OtherSt_RandNum, ErrStat, ErrMsg)
   CALL CheckError(ErrStat, ErrMsg)
               
ENDIF !WrACT

!..................................................................................................................................
! Generate full-field output files:
!..................................................................................................................................

   ! Are we generating binary FF files?
IF ( p%WrFile(FileExt_BTS) .OR. p%WrFile(FileExt_WND) .OR. p%WrFile(FileExt_HAWC) ) THEN
   CALL WrSum_InterpolatedHubStats(p, V)
      
   IF ( p%WrFile(FileExt_BTS) ) THEN
      CALL WrBinTURBSIM(p, V, ErrStat, ErrMsg)
      CALL CheckError(ErrStat, ErrMsg)
   END IF   
   
   IF ( p%WrFile(FileExt_WND) ) THEN      
      CALL WrBinBLADED(p, V, USig, VSig, WSig, ErrStat, ErrMsg)
      CALL CheckError(ErrStat, ErrMsg)
   END IF
   
   IF ( p%WrFile(FileExt_HAWC) ) THEN
      CALL WrBinHAWC( p, V, USig, VSig, WSig, ErrStat, ErrMsg )
      CALL CheckError(ErrStat, ErrMsg)
   END IF
   
END IF


   ! Are we generating formatted (text) FF files?
IF ( p%WrFile(FileExt_UVW) )  THEN
   CALL WrFormattedFF(p%RootName, p%grid, p%UHub, V)
ENDIF ! ( WrFile(FileExt_UVW) )

!..................................................................................................................................
! End:
!..................................................................................................................................

IF ( ALLOCATED( V  ) )  DEALLOCATE( V )


   ! Request CPU-time used.

CALL CPU_TIME ( CPUtime )
WRITE (p%US,"(//,'Processing complete.  ',A,' CPU seconds used.')")  TRIM( Num2LStr( CPUtime ) )
CALL WrScr1  (  ' Processing complete.  '//TRIM( Num2LStr( CPUtime ) )//' CPU seconds used.' )

CALL TS_End( p, OtherSt_RandNum )
CALL NormStop


CONTAINS
!..................................................................................................................................
SUBROUTINE CheckError(ErrID,Msg)
   INTEGER(IntKi), INTENT(IN)           :: ErrID       ! The error identifier (ErrStat)
   CHARACTER(*),   INTENT(IN)           :: Msg         ! The error message (ErrMsg)


   IF (ErrID /= ErrID_None) THEN
   
      IF (ErrID >= AbortErrLev) THEN

         IF (ALLOCATED(PhaseAngles)) DEALLOCATE(PhaseAngles)
         IF (ALLOCATED(S          )) DEALLOCATE(S          )
         IF (ALLOCATED(V          )) DEALLOCATE(V          )
         IF (ALLOCATED(U          )) DEALLOCATE(U          )
         IF (ALLOCATED(HWindDir   )) DEALLOCATE(HWindDir   )
         IF (ALLOCATED(VWindDir   )) DEALLOCATE(VWindDir   )
         
         
         if (p%US > 0) then
            WRITE (p%US, "(/'ERROR:  ', A / )") TRIM(Msg)
            WRITE (p%US, "('ABORTING PROGRAM.')" )
         end if
         
         CALL TS_end(p, OtherSt_RandNum)
         
         CALL ProgAbort ( TRIM(Msg), .FALSE., 5.0_ReKi )
         
      ELSE
         CALL WrScr(TRIM(Msg))
      END IF
      
   END IF
END SUBROUTINE CheckError

END PROGRAM
