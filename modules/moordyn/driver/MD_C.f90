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
!FIXME: add g_C, rhoW_C, WtrDepth_C, and PtfmInit_C to the interface  -- make PtfmInit_C an array size 6
SUBROUTINE MD_INIT_C(MD_InputFileName_C, InputFileNameLength_C, DT_C, NumChannels_C, OutputChannelNames_C, OutputChannelUnits_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='MD_INIT_C')

   !TEMPORARY hack until Waves handling is finalized
   USE WAVES, only: WaveGrid_n, WaveGrid_x0, WaveGrid_dx, WaveGrid_nx, WaveGrid_y0, WaveGrid_dy, WaveGrid_ny, WaveGrid_nz

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
    INTEGER(IntKi)                                                   :: ErrStat
    CHARACTER(ErrMsgLen)                                             :: ErrMsg
    INTEGER                                                          :: I

    integer(IntKi)                                                   :: ErrStat2       ! temporary error status of the operation
    character(ErrMsgLen)                                             :: ErrMsg2        ! temporary error message

    ! NOTE: Wave info will be handled differently in the future.  So the following is a temporary hack until that is finalized
    ! Hard coded for 10 wave steps.  Doesn't actually matter since it will get zeroed
    integer(IntKi)                                                   :: NStepWave = 10

    ! Initialize ErrStat
    ErrStat = ErrID_None
    ErrMsg  = ""

    ! Convert the MD input filename from C to Fortran
    MD_InputFileName = TRANSFER(MD_InputFileName_C, MD_InputFileName)
    I = INDEX(MD_InputFileName,C_NULL_CHAR) - 1            ! if this has a c null character at the end...
    IF ( I > 0 ) MD_InputFileName = MD_InputFileName(1:I)  ! remove it
    PRINT *, 'Inside MD_INIT_C: the passed input filename is ', MD_InputFileName

    ! Set other inputs for calling MD_Init
    DTcoupling               = REAL(DT_C, DbKi)
    InitInp%FileName         = MD_InputFileName
    InitInp%RootName         = 'MDroot'

    ! Environment variables -- These should be passed in from C.
    InitInp%g                = -9.806        ! Set this from a value passed in from C   -- check with Matt Hall on sign
    InitInp%rhoW             = 1000.0_ReKi   ! Set this from a value passed in from C
    InitInp%WtrDepth         = -100.0_ReKi   ! Set this from a value passed in from C   -- check with Matt Hall on sign

    ! Platform position (x,y,z,Rx,Ry,Rz) where rotations are small angle assumption in radians.
    ! This data is used to set the CoupledKinematics mesh that will be used at each timestep call
    call AllocAry (InitInp%PtfmInit, 6, 'InitInp%PtfmInit', ErrStat2, ErrMsg2 ); if (Failed()) return
    InitInp%PtfmInit         = (/ 0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 0.0_ReKi /)

    ! THIS IS A SHORT TERM HACK
    ! Fake wave info -- completely still, with no dynamic pressure terms
    ! Set wave info to zeros -- assume 10 timesteps for now (doesn't really matter since it isn't getting used)
    call AllocAry ( InitInp%WaveVel  ,NStepWave, WaveGrid_n, 3, 'InitInp%WaveVel' , ErrStat2, ErrMsg2 );    if (Failed()) return
    call AllocAry ( InitInp%WaveAcc  ,NStepWave, WaveGrid_n, 3, 'InitInp%WaveAcc' , ErrStat2, ErrMsg2 );    if (Failed()) return
    call AllocAry ( InitInp%WavePDyn ,NStepWave, WaveGrid_n,    'InitInp%WavePDyn', ErrStat2, ErrMsg2 );    if (Failed()) return
    call AllocAry ( InitInp%WaveElev ,NStepWave, WaveGrid_n,    'InitInp%WaveElev', ErrStat2, ErrMsg2 );    if (Failed()) return
    call AllocAry ( InitInp%WaveTime ,NStepWave,                'InitInp%WaveTime', ErrStat2, ErrMsg2 );    if (Failed()) return
    do i=1,NStepWave
       InitInp%WaveTime(i) = DTcoupling * REAL(i-1, DbKi)
    enddo
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

    if (ErrStat /= 0) then
        ErrStat_C = ErrID_Fatal
     else
        ErrStat_C = ErrID_None
     end if
     ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )

    PRINT*, "DONE WITH MD_INIT_C!"

contains

   subroutine Cleanup()
      ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )
      ErrStat_C = ErrStat
   end subroutine Cleanup

   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MD_Init_C')
      Failed =  ErrStat >= AbortErrLev
      if (Failed) call CleanUp()
   end function Failed

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
