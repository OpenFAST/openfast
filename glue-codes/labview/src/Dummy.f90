MODULE WaveTankTesting

IMPLICIT NONE
SAVE

contains

subroutine NoOpNoArgs() bind (C, NAME="NoOpNoArgs")
implicit none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: NoOpNoArgs
!GCC$ ATTRIBUTES DLLEXPORT :: NoOpNoArgs
#endif

print *, "Hello from NoOpNoArgs"

end subroutine


! SUBROUTINE NoOp(ErrStat_C, ErrMsg_C) bind (C, NAME="WaveTank_NoOp")
! IMPLICIT NONE
! #ifndef IMPLICIT_DLLEXPORT
! !DEC$ ATTRIBUTES DLLEXPORT :: WaveTank_NoOp
! !GCC$ ATTRIBUTES DLLEXPORT :: WaveTank_NoOp
! #endif

!     INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_C
!     CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_C(ErrMsgLen_C)

!     ! Local variables
!     INTEGER(C_INT)                          :: ErrStat_C2
!     CHARACTER(KIND=C_CHAR, LEN=ErrMsgLen_C) :: ErrMsg_C2

!     ! Initialize error handling
!     ErrStat_C = ErrID_None
!     ErrMsg_C  = " "//C_NULL_CHAR

!     ! No operation
!     ErrStat_C2 = ErrID_Info
!     ErrMsg_C2 = "Hi Stephen - No op here."

!     CALL SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'WaveTank_NoOp')

! END SUBROUTINE

END MODULE
