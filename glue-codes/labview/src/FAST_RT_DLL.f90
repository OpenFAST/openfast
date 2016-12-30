!  FAST_RT_DLL.f90 
! (c) 2009, 2012 National Renewable Energy Laboratory
!  Paul Fleming, National Wind Technology Center, September 2009, 2012
!  Bonnie Jonkman, National Wind Technology Center, October 2012
! 
!  Modification of FAST for Labview RT
!  Also includes code from FAST_Simulink Adaptation
!====================================================================================

subroutine FAST_RT_DLL_INIT (FileName_RT_Byte, FLen)

  ! Expose subroutine FAST_RT_DLL_INIT to users of this DLL
  !
  !DEC$ ATTRIBUTES DLLEXPORT::FAST_RT_DLL_INIT

USE                     NWTC_Library
USE							General, ONLY : PriFile, Cmpl4LV

USE                     FAST_IO_Subs  ! FAST_Input(), FAST_Begin()
USE                     FASTSubs      ! FAST_Initialize()

				! This sub-routine is called by RT to initialize all internal variables

IMPLICIT				NONE

INTEGER, PARAMETER         :: MaxFileNameLen = 100    
INTEGER(B1Ki)              :: FileName_RT_Byte(MaxFileNameLen)   ! FileName_RT_Byte

CHARACTER(MaxFileNameLen)  :: FileName_RT_Char        ! FileName_RT_Byte converted to ASCII characters
INTEGER						   :: FLen                    ! trim length of FileName_RT_Byte
INTEGER                    :: I                       ! temporary loop counter


IF ( FLen > MaxFileNameLen ) CALL ProgAbort('File name is too long in FAST_RT_DLL_INIT.')
DO I=1,FLen
   FileName_RT_Char(I:I) = ACHAR(FileName_RT_Byte(I))
END DO
!EQUIVALENCE(FileName_RT_Byte2,FileName_RT_Char)  !Make the character filename equivalent to incoming filename byte array

!FileName_RT_Byte2(:) = FileName_RT_Byte(:)



!Assign PriFile based on passed in string
PriFile = FileName_RT_Char(1:FLen)


     ! Open and read input files, initialize global parameters.
CALL FAST_Begin( PriFile, RootName, DirRoot )


!Set compiler flag for Simulink
Cmpl4LV   = .TRUE.

CALL FAST_Input()

   ! Set up initial values for all degrees of freedom.
CALL FAST_Initialize(p,x,y,OtherState)


end subroutine FAST_RT_DLL_INIT




!====================================================================================
subroutine FAST_RT_DLL_SIM (BlPitchCom_RT, YawPosCom_RT, YawRateCom_RT, ElecPwr_RT, GenTrq_RT, OutData_RT, Time_RT, HSSBrFrac_RT)


  ! Expose subroutine FAST_RT_DLL_SIM to users of this DLL
  !
  !DEC$ ATTRIBUTES DLLEXPORT::FAST_RT_DLL_SIM


USE                             SimCont      !ZTime

! These are needed for FirstTime = .FALSE.
USE                             DriveTrain   ! GenTrq and now also HSSBrFrac
USE                             TurbCont     ! BlPitch
USE                             TurbConf     ! NumBl
USE                             Blades       ! TipNode
USE                             Precision    ! ReKi
USE                             Features     ! CompAero
USE                             Output       ! for WrOutHdr

USE                           FASTSubs       ! TimeMarch()

IMPLICIT						NONE

				!  This sub-routine implements n-iterations of time step and returns outputs to Labview RT

  ! Variables
REAL(ReKi), INTENT(IN)       :: GenTrq_RT                          ! Mechanical generator torque.
REAL(ReKi), INTENT(IN)       :: ElecPwr_RT                         ! Electrical power
REAL(ReKi), INTENT(IN)       :: YawPosCom_RT                       ! Yaw position
REAL(ReKi), INTENT(IN)       :: YawRateCom_RT                      ! Yaw rate
REAL(ReKi), INTENT(IN)       :: BlPitchCom_RT  (*) 
REAL(ReKi), INTENT(OUT)      :: OutData_RT  (*) 
REAL(ReKi), INTENT(OUT)      :: Time_RT
REAL(ReKi), INTENT(IN)       :: HSSBrFrac_RT                       ! Brake Fraction

  !Copy in inputs from RT
   BlPitchCom = BlPitchCom_RT(1:NumBl)
   YawPosCom = YawPosCom_RT
   YawRateCom = YawRateCom_RT
   ElecPwr = ElecPwr_RT
   GenTrq= GenTrq_RT
   HSSBrFrac = HSSBrFrac_RT

   ! Set the command pitch angles to the actual pitch angles since we have no
   ! built-in pitch actuator:
    BlPitch = BlPitchCom


!Run simulation
Call TimeMarch( p_StrD, x_StrD, OtherSt_StrD, y_StrD, ErrStat, ErrMsg )


!Copy outputs
OutData_RT(1:p_StrD%NumOuts) = OutData(1:p_StrD%NumOuts)
Time_RT = ZTime;
OutData_RT(p_StrD%NumOuts+1) = TMax;
OutData_RT(p_StrD%NumOuts+2) = Time_RT;

end subroutine FAST_RT_DLL_SIM
