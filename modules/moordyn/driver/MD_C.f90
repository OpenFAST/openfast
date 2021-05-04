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

INTEGER, PARAMETER :: InputStringLength = 179          !< Fixed length for all lines of the string-based input file

CONTAINS

!===============================================================================================================
!---------------------------------------------- MD INIT --------------------------------------------------------
!===============================================================================================================
SUBROUTINE MD_INIT_C(MD_InputFileString_C, InputFileStringLength_C, DT_C, G_C, RHO_C, DEPTH_C, PtfmInit_C, NumChannels_C, OutputChannelNames_C, OutputChannelUnits_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='MD_INIT_C')

    !TEMPORARY hack until Waves handling is finalized
    USE WAVES, only: WaveGrid_n, WaveGrid_x0, WaveGrid_y0, WaveGrid_dx, WaveGrid_dy, WaveGrid_nx, WaveGrid_ny, WaveGrid_nz

    TYPE(C_PTR)                                    , INTENT(IN   )   :: MD_InputFileString_C
    INTEGER(C_INT)                                 , INTENT(IN   )   :: InputFileStringLength_C
    REAL(C_DOUBLE)                                 , INTENT(IN   )   :: DT_C
    REAL(C_DOUBLE)                                 , INTENT(IN   )   :: G_C
    REAL(C_DOUBLE)                                 , INTENT(IN   )   :: RHO_C
    REAL(C_DOUBLE)                                 , INTENT(IN   )   :: DEPTH_C
    REAL(C_FLOAT)                                  , INTENT(IN   )   :: PtfmInit_C(6)
    INTEGER(C_INT)                                 , INTENT(  OUT)   :: NumChannels_C
    TYPE(C_PTR)                                    , INTENT(  OUT)   :: OutputChannelNames_C
    TYPE(C_PTR)                                    , INTENT(  OUT)   :: OutputChannelUnits_C
    INTEGER(C_INT)                                 , INTENT(  OUT)   :: ErrStat_C
    CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: ErrMsg_C(1025) 

    ! Local Variables
    CHARACTER(InputStringLength), DIMENSION(InputFileStringLength_C) :: InputFileStrings
    CHARACTER(KIND=C_CHAR, LEN=1), DIMENSION(:), POINTER             :: character_pointer
    CHARACTER, DIMENSION(InputStringLength)                          :: single_line_character_array
    CHARACTER(InputStringLength)                                     :: single_line_chars
    CHARACTER(CHANLEN+1), ALLOCATABLE, TARGET                        :: tmp_OutputChannelNames_C(:)
    CHARACTER(CHANLEN+1), ALLOCATABLE, TARGET                        :: tmp_OutputChannelUnits_C(:)
    REAL(DbKi)                                                       :: DTcoupling
    INTEGER(IntKi)                                                   :: ErrStat, ErrStat2
    CHARACTER(ErrMsgLen)                                             :: ErrMsg, ErrMsg2
    INTEGER                                                          :: I, J

    ! NOTE: Wave info will be handled differently in the future.  So the following is a temporary hack until that is finalized
    ! Hard coded for 10 wave steps.  Doesn't actually matter since it will get zeroed
    INTEGER(IntKi)                                                   :: NStepWave = 10

    ! Initialize ErrStat
    ErrStat = ErrID_None
    ErrMsg  = ""

    ! Convert the string-input from C-style character arrays (char**) to a Fortran-style array of characters.
    ! TODO: Add error checking
    CALL C_F_pointer(MD_InputFileString_C, character_pointer, [InputFileStringLength_C * InputStringLength])
    DO i = 0, InputFileStringLength_C - 1
       single_line_character_array = character_pointer(i * InputStringLength + 1 : i * InputStringLength + InputStringLength)
       DO j = 1, InputStringLength
          single_line_chars(j:j) = single_line_character_array(j)
       END DO
       InputFileStrings(i + 1) = single_line_chars
    END DO
    PRINT *, 'InputFileStrings = ', InputFileStrings

    ! Store string-inputs as type FileInfoType
    CALL InitFileInfo(InputFileStrings, InitInp%PassedPrimaryInputData, ErrStat, ErrMsg)           
    IF (ErrStat .NE. 0) THEN 
       PRINT *, "MD_INIT_C: Failed to convert main input file to FileInfoType"
       PRINT *, ErrMsg
    END IF

    ! Set other inputs for calling MD_Init
    DTcoupling               = REAL(DT_C, DbKi)
    InitInp%FileName         = 'notUsed'
    InitInp%RootName         = 'MDroot'
    InitInp%UsePrimaryInputFile = .FALSE.

    ! Environment variables -- These should be passed in from C.
    InitInp%g                = REAL(G_C, DbKi)
    InitInp%rhoW             = REAL(RHO_C, DbKi)
    InitInp%WtrDepth         = REAL(DEPTH_C, DbKi)

    ! Platform position (x,y,z,Rx,Ry,Rz) -- where rotations are small angle assumption in radians.
    ! This data is used to set the CoupledKinematics mesh that will be used at each timestep call
    CALL AllocAry (InitInp%PtfmInit, 6, 'InitInp%PtfmInit', ErrStat2, ErrMsg2 ); IF (Failed()) RETURN
    DO I = 1,6
        InitInp%PtfmInit(I)  = REAL(PtfmInit_C(I),ReKi)
    END DO
    PRINT *, 'InitInp%PtfmInit = ', InitInp%PtfmInit

    ! Wave INformation - THIS IS A SHORT TERM HACK
    ! Fake wave info -- completely still, with no dynamic pressure terms
    ! Set wave info to zeros -- assume 10 timesteps for now (doesn't really matter since it isn't getting used)
    CALL AllocAry ( InitInp%WaveVel  ,NStepWave, WaveGrid_n, 3, 'InitInp%WaveVel' , ErrStat2, ErrMsg2 );    IF (Failed()) RETURN
    CALL AllocAry ( InitInp%WaveAcc  ,NStepWave, WaveGrid_n, 3, 'InitInp%WaveAcc' , ErrStat2, ErrMsg2 );    IF (Failed()) RETURN
    CALL AllocAry ( InitInp%WavePDyn ,NStepWave, WaveGrid_n,    'InitInp%WavePDyn', ErrStat2, ErrMsg2 );    IF (Failed()) RETURN
    CALL AllocAry ( InitInp%WaveElev ,NStepWave, WaveGrid_n,    'InitInp%WaveElev', ErrStat2, ErrMsg2 );    IF (Failed()) RETURN
    CALL AllocAry ( InitInp%WaveTime ,NStepWave,                'InitInp%WaveTime', ErrStat2, ErrMsg2 );    IF (Failed()) RETURN
    DO i=1,NStepWave
       InitInp%WaveTime(i) = DTcoupling * REAL(i-1, DbKi)
    END DO
    InitInp%WaveVel          = 0.0_ReKi
    InitInp%WaveAcc          = 0.0_ReKi
    InitInp%WavePDyn         = 0.0_ReKi
    InitInp%WaveElev         = 0.0_ReKi

    ! Call the main subroutine MD_Init
    CALL MD_Init(InitInp, u, p, x, xd, z, other, y, m, DTcoupling, InitOutData, ErrStat, ErrMsg)
    !FIXME: this may catch messages labelled as Info as fatal errors.  You probably don't want that.
    IF (ErrStat /= ErrID_None) THEN
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

    ! Clean up variables and set up for IFW_CALCOUTPUT_C
    CALL MD_DestroyInitInput( InitInp, ErrStat, ErrMsg )
    IF (ErrStat /= 0) THEN
        PRINT *, ErrMsg
    END IF

    IF (ErrStat /= 0) THEN
        ErrStat_C = ErrID_Fatal
    ELSE
        ErrStat_C = ErrID_None
    END IF
    ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )

    PRINT*, "DONE WITH MD_INIT_C!"

CONTAINS

    SUBROUTINE Cleanup()
        ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )
        ErrStat_C = ErrStat
    END SUBROUTINE Cleanup

    logical FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MD_Init_C')
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
    END FUNCTION Failed

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
SUBROUTINE MD_END_C(ErrStat_C,ErrMsg_C) BIND (C, NAME='MD_END_C')

    INTEGER(C_INT)                , INTENT(  OUT)      :: ErrStat_C
    CHARACTER(KIND=C_CHAR)        , INTENT(  OUT)      :: ErrMsg_C

    ! Local variables
    INTEGER                                            :: ErrStat
    CHARACTER(ErrMsgLen)                               :: ErrMsg

    ! Call the main subroutine MD_End
    CALL MD_End(u, p, x, xd, z, other, y, m, ErrStat , ErrMsg)
    IF (ErrStat .NE. 0) THEN
        PRINT *, "MD_END_C: MD_End failed"
        PRINT *, ErrMsg
    ELSE
        PRINT*, "MD_END_C: Successfully called MD_END ....."
    END IF

    ! Convert the outputs of MD_End from Fortran to C
    IF (ErrStat /= 0) THEN
        ErrStat_C = ErrID_Fatal
    ELSE
        ErrStat_C = ErrID_None
    END IF
    ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )

    PRINT*, "DONE WITH MD_END_C!"

END SUBROUTINE MD_END_C

END MODULE
