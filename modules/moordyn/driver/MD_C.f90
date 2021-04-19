!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2021 Nicole Mendoza
!
! This file is part of MoorDyn.
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
MODULE MoorDynAPI

    USE ISO_C_BINDING
    USE MoorDyn
    USE MoorDyn_Types
    USE NWTC_Library

IMPLICIT NONE

PUBLIC :: MD_INIT_C
PUBLIC :: MD_UPDATESTATES_C
PUBLIC :: MD_CALCOUTPUT_C
PUBLIC :: MD_END_C

! Global Variables
TYPE(MD_InitInputType)                  :: InitInp     !< Input data for initialization routine
TYPE(MD_InputType)                      :: u           !< An initial guess for the input; input mesh must be defined
TYPE(MD_ParameterType)                  :: p           !< Parameters
TYPE(MD_ContinuousStateType)            :: x           !< Initial continuous states
TYPE(MD_DiscreteStateType)              :: xd          !< Initial discrete states
TYPE(MD_ConstraintStateType)            :: z           !< Initial guess of the constraint states
TYPE(MD_OtherStateType)                 :: other       !< Initial other states
TYPE(MD_OutputType)                     :: y           !< Initial system outputs (outputs are not calculated; only the output mesh is initialized)
TYPE(MD_MiscVarType)                    :: m           !< Initial misc/optimization variables
TYPE(MD_InitOutputType)                 :: InitOutData !< Output for initialization routine

CONTAINS

!===============================================================================================================
!---------------------------------------------- MD INIT --------------------------------------------------------
!===============================================================================================================
SUBROUTINE MD_INIT_C(MD_InputFileName_C, InputFileNameLength_C, DT_C, NumChannels_C, OutputChannelNames_C, OutputChannelUnits_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='MD_INIT_C')

    INTEGER(C_INT)                                 , INTENT(IN   )   :: InputFileNameLength_C
    CHARACTER(KIND=C_CHAR)                         , INTENT(IN   )   :: MD_InputFileName_C(InputFileNameLength_C)
    REAL(C_DOUBLE)                                 , INTENT(IN   )   :: DT_C
    INTEGER(C_INT)                                 , INTENT(  OUT)   :: NumChannels_C
    TYPE(C_PTR)                                    , INTENT(  OUT)   :: OutputChannelNames_C
    TYPE(C_PTR)                                    , INTENT(  OUT)   :: OutputChannelUnits_C
    INTEGER(C_INT)                                 , INTENT(  OUT)   :: ErrStat_C
    CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: ErrMsg_C(1025) 

    ! Local Variables
    CHARACTER(KIND=C_CHAR, LEN=1), DIMENSION(:), POINTER             :: character_pointer
    CHARACTER(InputFileNameLength_C)                                 :: MD_InputFileName
    CHARACTER(CHANLEN+1), ALLOCATABLE, TARGET                        :: tmp_OutputChannelNames_C(:)
    CHARACTER(CHANLEN+1), ALLOCATABLE, TARGET                        :: tmp_OutputChannelUnits_C(:)
    REAL(DbKi)                                                       :: DTcoupling
    INTEGER                                                          :: ErrStat
    CHARACTER(ErrMsgLen)                                             :: ErrMsg
    INTEGER                                                          :: I

    ! Additional Inputs
    REAL(ReKi)                                                       :: g = -9.806       !< gravity constant [[m/s^2]]
    REAL(ReKi)                                                       :: rhoW = -1000.    !< sea density [[kg/m^3]]
    REAL(ReKi)                                                       :: WtrDepth = -100. !< depth of water [[m]]
    REAL(ReKi) , DIMENSION(:), ALLOCATABLE                           :: PtfmInit         !< initial position of platform(s) originally size 6 [-]
    REAL(ReKi) , DIMENSION(:,:,:), ALLOCATABLE                       :: WaveVel          !<  [-]
    REAL(ReKi) , DIMENSION(:,:,:), ALLOCATABLE                       :: WaveAcc          !<  [-]
    REAL(ReKi) , DIMENSION(:,:), ALLOCATABLE                         :: WavePDyn         !<  [-]
    REAL(ReKi) , DIMENSION(:,:), ALLOCATABLE                         :: WaveElev         !<  [-]
    REAL(DbKi) , DIMENSION(:), ALLOCATABLE                           :: WaveTime         !<  [-]
    
    ! Convert the MD input filename from C to Fortran
    MD_InputFileName = TRANSFER(MD_InputFileName_C, MD_InputFileName)
    I = INDEX(MD_InputFileName,C_NULL_CHAR) - 1            ! if this has a c null character at the end...
    IF ( I > 0 ) MD_InputFileName = MD_InputFileName(1:I)  ! remove it
    PRINT *, 'Inside MD_INIT_C: the passed input filename is ', MD_InputFileName

    ! Set other inputs for calling MD_Init
    ALLOCATE(PtfmInit(6))
    ALLOCATE(WaveVel())
    ALLOCATE(WaveAcc())
    ALLOCATE(WavePDyn())
    ALLOCATE(WaveElev())
    ALLOCATE(WaveTime())

    DTcoupling               = REAL(DT_C, DbKi)
    InitInp%FileName         = MD_InputFileName
    InitInp%RootName         = 'MDroot'
    InitInp%g                = g
    InitInp%rhoW             = rhoW
    InitInp%WtrDepth         = WtrDepth
    InitInp%PtfmInit         = PtfmInit
    InitInp%WaveVel          = WaveVel
    InitInp%WaveAcc          = WaveAcc
    InitInp%WavePDyn         = WavePDyn
    InitInp%WaveElev         = WaveElev
    InitInp%WaveTime         = WaveTime

    ! Call the main subroutine MD_Init
    CALL MD_Init(InitInp, u, p, x, xd, z, other, y, m, DTcoupling, InitOutData, ErrStat, ErrMsg)
    IF (ErrStat .NE. 0) THEN 
        PRINT *, "MD_INIT_C: Main MD_Init subroutine failed!"
        PRINT *, ErrMsg
     ELSE
        PRINT*, "MD_INIT_C: Successfully called MD_Init ....."
     END IF

    ! Convert the outputs of MD_Init from Fortran to C
     ALLOCATE(tmp_OutputChannelNames_C(size(InitOutData%writeOutputHdr)))
     ALLOCATE(tmp_OutputChannelUnits_C(size(InitOutData%writeOutputUnt)))
     NumChannels_C = size(InitOutData%writeOutputHdr)
     PRINT *, 'MD_INIT_C: The number of output channels is ', NumChannels_C

     DO I = 1,NumChannels_C
        tmp_OutputChannelNames_C(I) = TRANSFER(InitOutData%writeOutputHdr(I)//C_NULL_CHAR, tmp_OutputChannelNames_C(I))
        tmp_OutputChannelUnits_C(I) = TRANSFER(InitOutData%writeOutputUnt(I)//C_NULL_CHAR, tmp_OutputChannelUnits_C(I))
     END DO
     OutputChannelNames_C = C_LOC(tmp_OutputChannelNames_C)
     OutputChannelUnits_C = C_LOC(tmp_OutputChannelUnits_C)

    if (ErrStat /= 0) then
        ErrStat_C = ErrID_Fatal
     else
        ErrStat_C = ErrID_None
     end if
     ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )

    PRINT*, "DONE WITH MD_INIT_C!"

END SUBROUTINE MD_INIT_C

!===============================================================================================================
!---------------------------------------------- MD UPDATE STATES -----------------------------------------------
!===============================================================================================================
SUBROUTINE MD_UPDATESTATES_C() BIND (C, NAME='MD_UPDATESTATES_C')

! SUBROUTINE MD_UpdateStates( t, n, u, utimes, p, x, xd, z, other, m, ErrStat, ErrMsg)

    PRINT*, "DONE WITH MD_UPDATESTATES_C!"

END SUBROUTINE MD_UPDATESTATES_C

!===============================================================================================================
!---------------------------------------------- MD CALC OUTPUT -------------------------------------------------
!===============================================================================================================
SUBROUTINE MD_CALCOUTPUT_C() BIND (C, NAME='MD_CALCOUTPUT_C')

! SUBROUTINE MD_CalcOutput( t, u, p, x, xd, z, other, y, m, ErrStat, ErrMsg )

    PRINT*, "DONE WITH MD_CALCOUTPUT_C!"

END SUBROUTINE MD_CALCOUTPUT_C

!===============================================================================================================
!----------------------------------------------- MD END --------------------------------------------------------
!===============================================================================================================
SUBROUTINE MD_END_C() BIND (C, NAME='MD_END_C')

! SUBROUTINE MD_End(u, p, x, xd, z, other, y, m, ErrStat , ErrMsg)

    PRINT*, "DONE WITH MD_END_C!"

END SUBROUTINE MD_END_C

END MODULE