! Simple example of use of HSL_MA57
Module HSL_M57_INTERFACE

CONTAINS

!-----------------------------------------------------------------------------------------

subroutine HSL_MA57_2007(MATRIX, B, X, N, NE)

   USE HSL_MA57_DOUBLE
   IMPLICIT NONE
   INTEGER I,INFO, N, NE
   TYPE(ZD11_TYPE) MATRIX
   TYPE(MA57_CONTROL) CONTROL
   TYPE(MA57_AINFO) AINFO
   TYPE(MA57_FINFO) FINFO
   TYPE(MA57_SINFO) SINFO
   TYPE(MA57_FACTORS) FACTORS

!   DOUBLE PRECISION, ALLOCATABLE :: B(:),X(:)
   real(8):: B(N), X(N)

! Initialize the structures
      CALL MA57_INITIALIZE(FACTORS,CONTROL)

! Analyse
      CALL MA57_ANALYSE(MATRIX,FACTORS,CONTROL,AINFO)
      IF(AINFO%FLAG<0) THEN
         WRITE(6,'(A,I2)') &
            ' Failure of MA57_ANALYSE with AINFO%FLAG=', AINFO%FLAG
         STOP
      END IF

! Factorize
      CALL MA57_FACTORIZE(MATRIX,FACTORS,CONTROL,FINFO)
      IF(FINFO%FLAG<0) THEN
         WRITE(6,'(A,I2)') &
            ' Failure of MA57_FACTORIZE with FINFO%FLAG=', FINFO%FLAG
         STOP
      END IF

! Solve without refinement
      X = B
      CALL MA57_SOLVE(MATRIX,FACTORS,X,CONTROL,SINFO)
!      IF(SINFO%FLAG==0)WRITE(6,'(A,/,(3F20.16))')  &
!         ' Solution without refinement is',X

! Perform one refinement
      CALL MA57_SOLVE(MATRIX,FACTORS,X,CONTROL,SINFO,B)
!      IF(SINFO%FLAG==0)WRITE(6,'(A,/,(3F20.16))') &
!          ' Solution after one refinement is',X

! Clean up
!      DEALLOCATE(MATRIX%VAL, MATRIX%ROW, MATRIX%COL)
      CALL MA57_FINALIZE(FACTORS,CONTROL,INFO)

END SUBROUTINE HSL_MA57_2007

END Module HSL_M57_INTERFACE


