      Module EA16_INTERFACE

	CONTAINS

	subroutine EA16_double(N, NZA, MatrixA, IIA, &
     &               NZM, MatrixM, IIM,          &
     &               WHICH, BLK, NV,             &
     &               NWANT, OMEGA, PHI)

      USE HSL_MA57_DOUBLE
      IMPLICIT NONE

      INTEGER N, BLK, NWANT, NV, WHICH, LDV, LDBV
      INTEGER LIW, LW
      PARAMETER (LIW=500, LW=2500)
      INTEGER IWORK(LIW)
      DOUBLE PRECISION WORK(LW)
      real(8), allocatable:: V(:,:), BV(:,:)

      INTEGER IDO, MODE, LIWORK, LWORK, NEINEG, IPOS(10), INFO(20)
      DOUBLE PRECISION SIGMA, RANGE(2)

      INTEGER ICNTL(20)
      DOUBLE PRECISION CNTL(15)

! Variables for MA57

      INTEGER I,INFO57, NZA, NZM
      TYPE(ZD11_TYPE) MatrixA, MatrixM
      TYPE(MA57_CONTROL) CONTROL, CONTROLF
      TYPE(MA57_AINFO) AINFO, AINFOF
      TYPE(MA57_FINFO) FINFO, FINFOF
      TYPE(MA57_SINFO) SINFO, SINFOF
      TYPE(MA57_FACTORS) FACTORS, FACTORSF

! VARIABLES FOR M AND A
	INTEGER:: IIA(N+1), IIM(N+1)
	REAL(8):: A(NZA), M(NZM), B(N,NV)

! EIGENVALUE AND EIGENVECTORS
    REAL(8) :: OMEGA(NWANT), PHI(N, NWANT)

	INTEGER, PARAMETER:: KMAX = 10000
	INTEGER:: K_INTERATION
    REAL(8):: TEMPA, TEMPM
	REAL(8), allocatable:: tA(:), tM(:)

! MATRIX F
    TYPE(ZD11_TYPE) MatrixF
    INTEGER:: NZF, j
!    INTEGER, allocatable:: iRF(:), jCF(:)
!    REAL(8), allocatable:: F(:)

      LDV = N
      LDBV = N
      ALLOCATE(V(N,NV), BV(N,BLK))

! Set the default values of the control parameters for MA57.
      CALL MA57_INITIALIZE(FACTORS, CONTROL)
      CALL MA57_INITIALIZE(FACTORSF, CONTROLF)

! Set the default values of the control parameters for EA16.
      CALL EA16ID(ICNTL,CNTL)

! By selecting ICNTL(5)=2, ICNTL(18)=1, the pole is changed at
!   every restart.
      ICNTL(5) = 2
      ICNTL(18) = 1
	  ICNTL(6) = 1000

! Set range
!     RANGE(1) = 0
!	RANGE(2) = 1.5

! Allow a change of MODE
      ICNTL(7) = 1

! Compute the amount of storage required for this routine.
      CALL EA16AD(N,BLK,NWANT,NV,LIWORK,LWORK,ICNTL,INFO)


! Check the size of the workspace.
      IF (LIWORK.GT.LIW) THEN
         WRITE (6,'(A,I4)') 'Increase LIW to at least ',LIWORK
         STOP
      ELSE IF (LWORK.GT.LW) THEN
         WRITE (6,'(A,I4)') 'Increase LW to at least ',LWORK
         STOP
      END IF

! We start with MODE=3  (standard mode).
      MODE = 3

! Factorise M using MA57.
! Store the factorization in FACTORS.
      CALL MA57_ANALYSE(MATRIXM,FACTORS,CONTROL,AINFO)
      IF(AINFO%FLAG<0) THEN
         WRITE(6,'(A,I2)') &
            ' Failure of MA57_ANALYSE with AINFO%FLAG=', AINFO%FLAG
         STOP
      END IF

      CALL MA57_FACTORIZE(MATRIXM,FACTORS,CONTROL,FINFO)
      IF(FINFO%FLAG<0) THEN
         WRITE(6,'(A,I2)') &
            ' Failure of MA57_FACTORIZE with FINFO%FLAG=', FINFO%FLAG
         STOP
      END IF

! PREPARE FOR F = A -sigma*M
! Initialize MatrixF
! construct tA(nzf), tM(nzf), so that F(i) = tA(i) - sigma*tM(i)
        allocate(tA(nza+nzm), tM(nza+nzm))
        if (.not.allocated(MatrixF%VAL)) then
           allocate(MatrixF%VAL(nza+nzm))
           allocate(MatrixF%COL(nza+nzm))
           allocate(MatrixF%ROW(nza+nzm))
        end if
        NZF = 0
        DO 20 i=1,N
        DO 30 j = i, N
            tempA = GETAIJ(MatrixA%VAL, iiA, MatrixA%COL, N, NZA, i, j)
            tempM = GETAIJ(MatrixM%VAL, iiM, MatrixM%COL, N, NZM, i, j)
            
            if ((ABS(tempA).GT.1.0d-25) .or. (ABS(tempM).GT.1.0d-25)) then
                NZF = NZF + 1
                tA(nzf) = tempA
                tM(nzf) = tempM
                MatrixF%VAL(NZF) = 0
                MatrixF%ROW(NZF) = i
                MatrixF%COL(NZF) = j
            end if

30           CONTINUE
20         CONTINUE
        MatrixF%N = N
        MatrixF%NE = NZF

!        DO i = 1, nzf
!            !write(*, '(I6, 3x, I6, 3x, I6, 3x, f10.4)') i, MatrixF%ROW(i), MatrixF%COL(i), MatrixF%VAL(i)
!            write(*, '(I6, 3x, f10.1, f10.1)') i, tA(i), tM(i)
!        END DO

!************************************************************************
!  Iteration
!************************************************************************
! Initialise reverse communication parameter IDO
      IDO  = 0
	
	K_INTERATION = 0
 1    CONTINUE
         CALL EA16BD(N, BLK, NWANT, NV, MODE, WHICH, IDO, IPOS, &
     &               V, LDV, BV, LDBV, RANGE, SIGMA, NEINEG,    &
     &               IWORK, LIWORK, WORK, LWORK, ICNTL, CNTL,   &
     &               INFO)

! Reverse communication action
         IF (IDO.EQ.100) THEN
! Finished
            GO TO 2
         ELSE IF (IDO.EQ.1) THEN

            IF (MODE.EQ.3) THEN

! Compute  V(:,IPOS(3):IPOS(4)) = M^-1 * A * V(:,IPOS(1):IPOS(2))
! We first compute  A * V(:,IPOS(1):IPOS(2)) and then
! solve a linear system with M, using the factorization in F.

               CALL MATVEC(MatrixA%VAL, N, NZA, MATRIXA%COL,  &
     &                     MATRIXA%ROW, V(:,IPOS(1):IPOS(2)), &
     &                     V(:,IPOS(3):IPOS(4)) )

! Compute  V(:,IPOS(3)) = inv(F) * V(:,IPOS(1)) using MA57
               CALL MA57_SOLVE(MATRIXM,FACTORS,V(:,IPOS(3):IPOS(4)),CONTROL,SINFO)

            ELSE IF (MODE.EQ.4) THEN

! Compute V(:,IPOS(3):IPOS(4)) = (A - SIGMA M)^-1 * M *
!                                V(:,IPOS(1):IPOS(2))
! with BV(1:N, IPOS(3):IPOS(4)) = M * V(1:N,IPOS(1):IPOS(2))
! Compute  V(:,IPOS(3):IPOS(4)) = inv(F) * BV using MA57
               ! write(*, *) " mode = 4"
               CALL DLACPY('All', N, BLK, BV, &
               &                     LDBV, V(1:N, IPOS(3):IPOS(4)), LDV)
               CALL MA57_SOLVE(MATRIXF,FACTORSF,V(1:N, IPOS(3):IPOS(4)),CONTROLF,SINFOF)

            END IF

         ELSE IF (IDO.EQ.2) THEN

! Compute  BV(:,IPOS(3):IPOS(4)) = M * V(:,IPOS(1):IPOS(2))

             CALL MATVEC(MatrixM%VAL, N, NZM, MATRIXM%COL, MATRIXM%ROW, &
     &                   V(1:N,IPOS(1):IPOS(2)), BV(1:N,IPOS(3):IPOS(4)) )

         ELSE IF (IDO.EQ.4) THEN

! Form F = A - SIGMA * M
! Factorize F using MA57

          DO i = 1, nzf
              MatrixF%VAL(i) = tA(i) - sigma*tM(i)
          END DO

          CALL MA57_ANALYSE(MATRIXF,FACTORSF,CONTROLF,AINFOF)
          IF(AINFOF%FLAG<0) THEN
             WRITE(6,'(A,I2)') &
                ' Failure of MA57_ANALYSE with AINFOF%FLAG=', AINFOF%FLAG
             STOP
          END IF
          
          CALL MA57_FACTORIZE(MATRIXF,FACTORSF,CONTROLF,FINFOF)
          IF(FINFOF%FLAG<0) THEN
             WRITE(6,'(A,I2)') &
                ' Failure of MA57_FACTORIZE with FINFOF%FLAG=', FINFOF%FLAG
                ! If failure, flag SIGMA as unusable
                IDO = -4
          END IF

          ! Set the number of negative eigenvalues
            NEINEG = finfof%neig

         ELSE IF (IDO.EQ.5) THEN

! Change the mode from 3 into 4.
            MODE = 4

         END IF

	K_INTERATION = K_INTERATION +1
   ! write(*, '(A, I6, A, e10.4)') "K_interation = ",  K_INTERATION, " sigma = ", sigma

	IF (K_INTERATION >= KMAX) THEN
		WRITE(*, *) "ITERATION =  KMAX ! "
		go to 2
	END IF

         GO TO 1

 2    CONTINUE
!************************************************************************
!  end Iteration
!************************************************************************

 
	write(*,*) "K_INTERATION = ", K_INTERATION

! Check for errors
      IF (INFO(1).LT.0) THEN
         WRITE (6,*) 'There was an error : INFO(1) = ', INFO(1)
         STOP
      END IF

! SET VALUES FOR OMEGA AND PHI
!     write(*, *) " "
!	write(*, *) "IPOS(5)", work(ipos(5))
!	write(*, *) "IPOS(6)", work(ipos(6))
!	write(*, *) " "
	DO I = 1, NWANT
		OMEGA(i) = WORK(ipos(5)+i-1)
	END DO

	PHI = V(:, 1:NWANT)

! Print some information
      WRITE (6,'(A,I5)') 'Number of iterations      : ', INFO(3)
      WRITE (6,'(A,I5)') 'Number of products with M : ', INFO(4)
      WRITE (6,'(A,I5)') 'Number of factorizations  : ', INFO(6)

! Print the computed eigenvalues
      WRITE (6, '(A/8(e15.7))') 'Eigenvalues : ', &
     &      (WORK(I),I=IPOS(5),IPOS(6))
      WRITE (6,'(A,e15.7,A)') 'All eigenvalues larger than ',range(1), &
     &               ' are computed.'

          if (allocated(MatrixF%VAL)) deallocate(MatrixF%VAL)
          if (allocated(MatrixF%ROW)) deallocate(MatrixF%ROW)
          if (allocated(MatrixF%COL)) deallocate(MatrixF%COL)

          if (allocated(MatrixA%VAL)) deallocate(MatrixA%VAL)
          if (allocated(MatrixA%ROW)) deallocate(MatrixA%ROW)
          if (allocated(MatrixA%COL)) deallocate(MatrixA%COL)

          if (allocated(MatrixM%VAL)) deallocate(MatrixM%VAL)
          if (allocated(MatrixM%ROW)) deallocate(MatrixM%ROW)
          if (allocated(MatrixM%COL)) deallocate(MatrixM%COL)

          deallocate(V)
          deallocate(BV)
          deallocate(tA)
          deallocate(tM)
! Clean up
!      DEALLOCATE(MATRIX%VAL, MATRIX%ROW, MATRIX%COL)
      CALL MA57_FINALIZE(FACTORS,CONTROL,INFO57)
      CALL MA57_FINALIZE(FACTORSf,CONTROLf, INFO57)
      END subroutine EA16_double

!*******************************************************
      SUBROUTINE MATVEC(MAT, N, NZ, ICN, IRN, X, Y)

      INTEGER N, NZ
      INTEGER ICN(NZ), IRN(NZ)
      DOUBLE PRECISION MAT(NZ)
      real(8), intent(in):: X(N)
      real(8), intent(out):: Y(N)

      INTEGER I
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.D0)

      DO 10 I=1,N
         Y(I) = ZERO
 10   CONTINUE

      DO 20 I=1,NZ
         Y(IRN(I)) = Y(IRN(I)) + MAT(I) * X(ICN(I))
         IF (IRN(I).NE.ICN(I)) &
     &      Y(ICN(I)) = Y(ICN(I)) + MAT(I) * X(IRN(I))
 20   CONTINUE
      RETURN
      END subroutine Matvec

!==================================================================
! get A(i, j)
!
!------------------------------------------------------------------
	Function GETAIJ(A, iA, jA, N, nnz, i, j)

	IMPLICIT NONE
	integer:: N, nnz
	real(8):: A(nnz), GetAij
	integer:: jA(nnz), iA(N+1)
	integer:: i, j, k1, k2, jj

	k1 = iA(i)
	k2 = iA(i+1)-1

	DO jj = k1, k2
		IF (j == jA(jj)) THEN
			GetAij = A(jj)
			return;
		END IF
	END DO
		
	GetAij = 0


	END function
	
	END	Module EA16_INTERFACE
