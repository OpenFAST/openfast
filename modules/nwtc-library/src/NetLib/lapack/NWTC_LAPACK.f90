!**********************************************************************************************************************************
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
!
!> This code provides a wrapper for the LAPACK routines currently used at the NWTC (mainly codes in the FAST framework). This 
!! enables us to call generic routines (not single- or double-precision specific ones) so that we don't have to change source
!! code to compile in double vs. single precision.   
!
!**********************************************************************************************************************************
MODULE NWTC_LAPACK

   USE NWTC_Base        ! we only need the precision and error level constants
   
   
!   USE, INTRINSIC               :: ISO_C_Binding, only: C_FLOAT, C_DOUBLE          ! this is included in NWTC_Library

      ! Notes:

         ! Your project must include the following files:
         ! From the NWTC Subroutine Library:
         !     SingPrec.f90          [from NWTC Library]
         !     Sys*.f90              [from NWTC Library]
         !     NWTC_Base.f90         [from NWTC Library]
         ! lapack library (preferably a binary, but available in source form from http://www.netlib.org/, too)
         ! This wrapper file:
         !     NWTC_LAPACK.f90

   !INTEGER, PARAMETER  :: Lib_ReKi = SiKi   !
   !INTEGER, PARAMETER  :: Lib_DbKi = R8Ki   ! DbKi
   !
   ! bjj: when using the built-in (or dynamic) lapack libraries, S=Real(SiKi); D=Real(R8Ki).
   !      if people are compiling the lapack source, S=real; D=double precision. (default real and doubles)
   !      we need to check this somehow to make sure the right routines are called.
   ! (or define a directive)

   ! http://www.netlib.org/lapack/explore-html/ 
   
   
   IMPLICIT  NONE

   !> Computes the solution to system of linear equations A * X = B for GB matrices.
   INTERFACE LAPACK_gbsv 
      MODULE PROCEDURE LAPACK_dgbsv
      MODULE PROCEDURE LAPACK_sgbsv
   END INTERFACE

   !> Computes scalar1*op( A )*op( B ) + scalar2*C where op(x) = x or op(x) = x**T for matrices A, B, and C.
   INTERFACE LAPACK_gemm   
      MODULE PROCEDURE LAPACK_dgemm
      MODULE PROCEDURE LAPACK_sgemm
   END INTERFACE
   
   !> Computes the solution to system of linear equations A * X = B for GE matrices.
   INTERFACE LAPACK_gesv 
      MODULE PROCEDURE LAPACK_dgesv
      MODULE PROCEDURE LAPACK_sgesv
   END INTERFACE

   !> Factor matrix into A=PLU.
   INTERFACE LAPACK_getrf 
      MODULE PROCEDURE LAPACK_dgetrf
      MODULE PROCEDURE LAPACK_sgetrf
   END INTERFACE

   !> Compute the inverse of a matrix using the LU factorization.
   INTERFACE LAPACK_getri 
      MODULE PROCEDURE LAPACK_dgetri
      MODULE PROCEDURE LAPACK_sgetri
   END INTERFACE

   !> Solve system(s) of linear equations Ax=PLUx=b.
   INTERFACE LAPACK_getrs 
      MODULE PROCEDURE LAPACK_dgetrs
      MODULE PROCEDURE LAPACK_sgetrs
      MODULE PROCEDURE LAPACK_dgetrs1
      MODULE PROCEDURE LAPACK_sgetrs1
   END INTERFACE

   !> Compute generalized eigenvalues and/or eigenvectors for a pair of N-by-N real nonsymmetric matrices (A,B).
   INTERFACE LAPACK_ggev 
      MODULE PROCEDURE LAPACK_dggev
      MODULE PROCEDURE LAPACK_sggev
   END INTERFACE

   !> Compute the solution to system of linear equations A * X = B for PO matrices.
   INTERFACE LAPACK_posv 
      MODULE PROCEDURE LAPACK_dposv
      MODULE PROCEDURE LAPACK_sposv
   END INTERFACE

   !> Compute the Cholesky factorization of a real symmetric positive definite matrix A stored in packed format.
   INTERFACE LAPACK_pptrf 
      MODULE PROCEDURE LAPACK_dpptrf
      MODULE PROCEDURE LAPACK_spptrf
   END INTERFACE

   !> Compute the SVD for a general matrix A = USV^T.
   INTERFACE LAPACK_gesvd
      MODULE PROCEDURE LAPACK_dgesvd
      MODULE PROCEDURE LAPACK_sgesvd
   END INTERFACE

   !> Unpack  packed (1D) to regular matrix format (2D)
   INTERFACE LAPACK_TPTTR  
      MODULE PROCEDURE LAPACK_STPTTR
      MODULE PROCEDURE LAPACK_DTPTTR
   END INTERFACE   
   
   !> straight-up lapack routines (from ExtPtfm_MCKF):
   INTERFACE LAPACK_COPY
       SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
           USE Precision, only: R8Ki
           INTEGER    :: INCX,INCY,N
           real(R8Ki) :: DX(*),DY(*)
       ENDSUBROUTINE
       SUBROUTINE SCOPY(N,X,INCX,Y,INCY)
           USE Precision, only: SiKi
           INTEGER    :: INCX,INCY,N
           real(SiKi) :: X(*),Y(*)
       ENDSUBROUTINE
   END INTERFACE
   
   INTERFACE LAPACK_GEMV
       SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
           USE Precision, only: R8Ki
           real(R8Ki) :: ALPHA,BETA
           integer    :: INCX,INCY,LDA,M,N
           character  :: TRANS
           real(R8Ki) :: A(LDA,*),X(*),Y(*)
       ENDSUBROUTINE
       SUBROUTINE SGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
           USE Precision, only: SiKi
           real(SiKi) :: ALPHA,BETA
           integer    :: INCX,INCY,LDA,M,N
           character  :: TRANS
           real(SiKi) :: A(LDA,*),X(*),Y(*)
       ENDSUBROUTINE
   END INTERFACE LAPACK_GEMV
   
   INTERFACE LAPACK_AXPY
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
           USE Precision, only: R8Ki
           real(R8Ki) :: DA
           integer    :: INCX,INCY,N
           real(R8Ki) :: DX(*),DY(*)
       ENDSUBROUTINE
      SUBROUTINE SAXPY(N,A,X,INCX,Y,INCY)
           USE Precision, only: SiKi
           real(SiKi) :: A
           integer    :: INCX,INCY,N
           real(SiKi) :: X(*),Y(*)
       ENDSUBROUTINE
   END INTERFACE

   CONTAINS

!=======================================================================
!> general banded solve: Computes the solution to system of linear equations A * X = B for GB (general, banded) matrices.
!! use LAPACK_GBSV (nwtc_lapack::lapack_gbsv) instead of this specific function.
   SUBROUTINE LAPACK_DGBSV( N, KL, KU, NRHS, AB, IPIV, B, ErrStat, ErrMsg )


      ! passed parameters

      INTEGER,         intent(in   ) :: KL                !< The number of subdiagonals within the band of A.  KL >= 0.
      INTEGER,         intent(in   ) :: KU                !< The number of superdiagonals within the band of A.  KU >= 0.
      INTEGER,         intent(in   ) :: N                 !< The number of linear equations, i.e., the order of the matrix A.  N >= 0.
      INTEGER,         intent(in   ) :: NRHS              !< The number of right hand sides, i.e., the number of columns of the matrix B.  NRHS >= 0.

      !     .. Array Arguments ..
      REAL(R8Ki)      ,intent(inout) :: AB( :, : )        !< On entry, the matrix A in band storage, in rows KL+1 to 2*KL+KU+1; rows 1 to KL of the array need not be set.
                                                          !! The j-th column of A is stored in the j-th column of the array AB as follows:
                                                          !!    AB(KL+KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+KL)
                                                          !! On exit, details of the factorization: U is stored as an upper triangular band matrix with KL+KU superdiagonals in
                                                          !! rows 1 to KL+KU+1, and the multipliers used during the factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
      REAL(R8Ki)      ,intent(inout) :: B( :, : )         !< On entry, the N-by-NRHS right hand side matrix B.  On exit, if INFO = 0, the N-by-NRHS solution matrix X.
      INTEGER,         intent(  out) :: IPIV( : )         !< The pivot indices that define the permutation matrix P; row i of the matrix was interchanged with row IPIV(i).

      INTEGER(IntKi),  intent(  out) :: ErrStat           !< Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            !< Message describing error

         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value; > 0: if INFO = i, U(i,i) is exactly zero.
                                                          ! > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed, but the factor U is exactly singular, so the solution could not be computed.
      INTEGER                        :: LDAB              ! The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
      INTEGER                        :: LDB               ! The leading dimension of the array B.   LDB  >= max(1,N).


      LDAB  = SIZE(AB,1)
      LDB   = SIZE(B, 1)



      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL dgbsv (N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO)

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_DGBSV: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_DGBSV: U('//TRIM(ErrMsg)//','//TRIM(ErrMsg)//')=0. Factor U is exactly singular.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_DGBSV
!=======================================================================
!> general banded solve: Computes the solution to system of linear equations A * X = B for GB (general, banded) matrices.
!! use LAPACK_GBSV (nwtc_lapack::lapack_gbsv) instead of this specific function.
   SUBROUTINE LAPACK_SGBSV( N, KL, KU, NRHS, AB, IPIV, B, ErrStat, ErrMsg )


      ! passed parameters

      INTEGER,         intent(in   ) :: KL                !< The number of subdiagonals within the band of A.  KL >= 0.
      INTEGER,         intent(in   ) :: KU                !< The number of superdiagonals within the band of A.  KU >= 0.
      INTEGER,         intent(in   ) :: N                 !< The number of linear equations, i.e., the order of the matrix A.  N >= 0.
      INTEGER,         intent(in   ) :: NRHS              !< The number of right hand sides, i.e., the number of columns of the matrix B.  NRHS >= 0.

      !     .. Array Arguments ..
      REAL(SiKi)      ,intent(inout) :: AB( :, : )        !< On entry, the matrix A in band storage, in rows KL+1 to 2*KL+KU+1; rows 1 to KL of the array need not be set.
                                                          !! The j-th column of A is stored in the j-th column of the array AB as follows:
                                                          !!    AB(KL+KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+KL)
                                                          !! On exit, details of the factorization: U is stored as an upper triangular band matrix with KL+KU superdiagonals in
                                                          !! rows 1 to KL+KU+1, and the multipliers used during the factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
      REAL(SiKi)      ,intent(inout) :: B( :, : )         !< On entry, the N-by-NRHS right hand side matrix B.  On exit, if INFO = 0, the N-by-NRHS solution matrix X.
      INTEGER,         intent(  out) :: IPIV( : )         !< The pivot indices that define the permutation matrix P; row i of the matrix was interchanged with row IPIV(i).

      INTEGER(IntKi),  intent(  out) :: ErrStat           !< Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            !< Message describing error

         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value; > 0: if INFO = i, U(i,i) is exactly zero.
                                                          ! > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed, but the factor U is exactly singular, so the solution could not be computed.

      INTEGER                        :: LDAB              ! The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
      INTEGER                        :: LDB               ! The leading dimension of the array B.   LDB  >= max(1,N).


      LDAB  = SIZE(AB,1)
      LDB   = SIZE(B, 1)


      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL SGBSV (N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO)

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_SGBSV: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_SGBSV: U('//TRIM(ErrMsg)//','//TRIM(ErrMsg)//')=0. Factor U is exactly singular.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_SGBSV
!=======================================================================
!> general matrix multiply: computes C = alpha*op( A )*op( B ) + beta*C where op(x) = x or op(x) = x**T for matrices A, B, and C
!! use LAPACK_GEMM (nwtc_lapack::lapack_gemm) instead of this specific function.
   SUBROUTINE LAPACK_DGEMM( TRANSA, TRANSB, ALPHA, A, B, BETA, C, ErrStat, ErrMsg )

         ! passed parameters

      CHARACTER(1),    intent(in   ) :: TRANSA            !< On entry, TRANSA specifies the form of op( A ) to be used in the matrix multiplication as follows:
                                                          !!     TRANSA = 'N' or 'n', op( A ) = A.
                                                          !!     TRANSA = 'T' or 't', op( A ) = A**T.
      CHARACTER(1),    intent(in   ) :: TRANSB            !< On entry, TRANSB specifies the form of op( A ) to be used in the matrix multiplication as follows:
                                                          !!     TRANSB = 'N' or 'n', op( B ) = B.
                                                          !!     TRANSB = 'T' or 't', op( B ) = B**T.

      REAL(R8Ki)      ,intent(in   ) :: ALPHA             !< On entry, ALPHA specifies the scalar alpha.
      REAL(R8Ki)      ,intent(in   ) :: BETA              !< On entry, BETA specifies the scalar beta. When BETA is supplied as zero then C need not be set on input.
      REAL(R8Ki)      ,intent(in   ) :: A( :, : )         !< Matrix A
      REAL(R8Ki)      ,intent(in   ) :: B( :, : )         !< Matrix B
      REAL(R8Ki)      ,intent(inout) :: C( :, : )         !< Matrix C: Before entry, C must contain the matrix C, except when beta is zero, in which case C need not
                                                          !< be set on entry. On exit, the array C is overwritten by the m by n matrix ( alpha*op( A )*op( B ) + beta*C ).

      INTEGER(IntKi),  intent(  out) :: ErrStat           !< Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            !< Message describing error
                                                          
         ! local variables                                                                   
                                                                                                                                                                                                                  
      INTEGER                        :: M                 ! M specifies the number of rows of the matrix op(A)
      INTEGER                        :: K                 ! K specifies the number of columns of the matrix op(A)
      INTEGER                        :: N                 ! N specifies the number of columns of the matrix op(B)
      INTEGER                        :: KB                ! KB specifies the number of rows of the matrix op(B)
                
      INTEGER                        :: LDA               ! LDA specifies the first dimension of A as declared in the calling (sub) program. When TRANSA = 'N' or 'n' then
                                                          ! LDA must be at least max( 1, m ), otherwise LDA must be at least max( 1, k ).

      INTEGER                        :: LDB               ! LDB specifies the first dimension of B as declared in the calling (sub) program. When TRANSB = 'N' or 'n' then
                                                          ! LDB must be at least max( 1, k ), otherwise LDB must be at least max( 1, n ).
                                                  
      CHARACTER(*), PARAMETER        :: RoutineName = 'LAPACK_DGEMM'
                                                                                                                                                                              
      LDA = SIZE(A,1)
      LDB = SIZE(B,1)
      
      IF (INDEX('Nn',TransA) > 0) THEN
         M = SIZE(A,1)
         K = SIZE(A,2)
      ELSE
         M = SIZE(A,2)
         K = SIZE(A,1)
      END IF

      IF (INDEX('Nn',TransB) > 0) THEN
         N = SIZE(B,2)
         Kb = SIZE(B,1)
      ELSE
         N = SIZE(B,1)
         Kb = SIZE(B,2)
      END IF

      ErrStat = ErrID_None
      ErrMsg  = ""

      
      IF ( K /= Kb ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg = RoutineName//":Size of Matrix A is incompatible with size of Matrix B."
         RETURN
      END IF
      
      IF ( M /= SIZE(C,1) .OR. N /= SIZE(C,2) ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg = RoutineName//":Size of Matrix C is incompatible with Matrix A and Matrix B."
         RETURN
      END IF
      
      IF ( M == 0 .or. N == 0 ) THEN
         ! this is a null case...
         RETURN
      END IF
      
      IF (K == 0) THEN
         C = C * beta ! A*B is null
         RETURN
      END IF
      
      
      CALL DGEMM (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, M)


   RETURN
   END SUBROUTINE LAPACK_DGEMM
!=======================================================================
!> general matrix multiply: computes C = alpha*op( A )*op( B ) + beta*C where op(x) = x or op(x) = x**T for matrices A, B, and C
!! use LAPACK_GEMM (nwtc_lapack::lapack_gemm) instead of this specific function.
   SUBROUTINE LAPACK_SGEMM( TRANSA, TRANSB, ALPHA, A, B, BETA, C, ErrStat, ErrMsg )

         ! passed parameters

      CHARACTER(1),    intent(in   ) :: TRANSA            !< On entry, TRANSA specifies the form of op( A ) to be used in the matrix multiplication as follows:
                                                          !!     TRANSA = 'N' or 'n', op( A ) = A.
                                                          !!     TRANSA = 'T' or 't', op( A ) = A**T.
      CHARACTER(1),    intent(in   ) :: TRANSB            !< On entry, TRANSB specifies the form of op( A ) to be used in the matrix multiplication as follows:
                                                          !!     TRANSB = 'N' or 'n', op( B ) = B.
                                                          !!     TRANSB = 'T' or 't', op( B ) = B**T.

      REAL(SiKi)      ,intent(in   ) :: ALPHA             !< On entry, ALPHA specifies the scalar alpha.
      REAL(SiKi)      ,intent(in   ) :: BETA              !< On entry, BETA specifies the scalar beta. When BETA is supplied as zero then C need not be set on input.
      REAL(SiKi)      ,intent(in   ) :: A( :, : )         !< Matrix A
      REAL(SiKi)      ,intent(in   ) :: B( :, : )         !< Matrix B
      REAL(SiKi)      ,intent(inout) :: C( :, : )         !< Matrix C: Before entry, C must contain the matrix C, except when beta is zero, in which case C need not
                                                          !! be set on entry. On exit, the array C is overwritten by the m by n matrix ( alpha*op( A )*op( B ) + beta*C ).

      INTEGER(IntKi),  intent(  out) :: ErrStat           !< Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            !< Message describing error
                                                          
         ! local variables                                                                   
                                                                                                                                                                                                                  
      INTEGER                        :: M                 ! M specifies the number of rows of the matrix op(A)
      INTEGER                        :: K                 ! K specifies the number of columns of the matrix op(A)
      INTEGER                        :: N                 ! N specifies the number of columns of the matrix op(B)
      INTEGER                        :: KB                ! KB specifies the number of rows of the matrix op(B)
                
      INTEGER                        :: LDA               ! LDA specifies the first dimension of A as declared in the calling (sub) program. When TRANSA = 'N' or 'n' then
                                                          ! LDA must be at least max( 1, m ), otherwise LDA must be at least max( 1, k ).

      INTEGER                        :: LDB               ! LDB specifies the first dimension of B as declared in the calling (sub) program. When TRANSB = 'N' or 'n' then
                                                          ! LDB must be at least max( 1, k ), otherwise LDB must be at least max( 1, n ).
                                                  
      CHARACTER(*), PARAMETER        :: RoutineName = 'LAPACK_SGEMM'
                                                                                                                                                                              
      LDA = SIZE(A,1)
      LDB = SIZE(B,1)
      
      IF (INDEX('Nn',TransA) > 0) THEN
         M = SIZE(A,1)
         K = SIZE(A,2)
      ELSE
         M = SIZE(A,2)
         K = SIZE(A,1)
      END IF

      IF (INDEX('Nn',TransB) > 0) THEN
         N = SIZE(B,2)
         Kb = SIZE(B,1)
      ELSE
         N = SIZE(B,1)
         Kb = SIZE(B,2)
      END IF

      ErrStat = ErrID_None
      ErrMsg  = ""

      
      IF ( K /= Kb ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg = RoutineName//":Size of Matrix A is incompatible with size of Matrix B."
         RETURN
      END IF
      
      IF ( M /= SIZE(C,1) .OR. N /= SIZE(C,2) ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg = RoutineName//":Size of Matrix C is incompatible with Matrix A and Matrix B."
         RETURN
      END IF
      
      IF ( M == 0 .or. N == 0 ) THEN
         ! this is a null case...
         RETURN
      END IF
      
      IF (K == 0) THEN
         C = C * beta ! A*B is null
         RETURN
      END IF
         
      
      CALL SGEMM (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, M)


   RETURN
   END SUBROUTINE LAPACK_SGEMM   
!=======================================================================
!> general solve: Computes the solution to system of linear equations A * X = B for GE matrices.
!! use LAPACK_GESV (nwtc_lapack::lapack_gesv) instead of this specific function.
   SUBROUTINE LAPACK_DGESV ( N, A, IPIV, B, ErrStat, ErrMsg )

      ! passed parameters

      INTEGER,         intent(in   ) :: N                 !< The number of linear equations, i.e., the order of the matrix A.  N >= 0.

      !     .. Array Arguments ..
      REAL(R8Ki)      ,intent(inout) :: A( :, : )         !< On entry, the N-by-N coefficient matrix A.  On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored.
      REAL(R8Ki)      ,intent(inout) :: B( :, : )         !< On entry, the N-by-NRHS matrix of right hand side matrix B.  On exit, if INFO = 0, the N-by-NRHS solution matrix X.
      INTEGER,         intent(  out) :: IPIV( : )         !< The pivot indices that define the permutation matrix P; row i of the matrix was interchanged with row IPIV(i).

      INTEGER(IntKi),  intent(  out) :: ErrStat           !< Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            !< Message describing error

         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value; > 0: if INFO = i, U(i,i) is exactly zero.
                                                          ! > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed, but the factor U is exactly singular, so the solution could not be computed.
      INTEGER                        :: NRHS              ! The number of right hand sides, i.e., the number of columns of the matrix B.  NRHS >= 0.
      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,N).
      INTEGER                        :: LDB               ! The leading dimension of the array B.  LDB >= max(1,N).

      
      LDA  = SIZE(A,1)
      LDB  = SIZE(B,1)
      NRHS = SIZE(B,2)

      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_DGESV: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_DGESV: U('//TRIM(ErrMsg)//','//TRIM(ErrMsg)//')=0. Factor U is exactly singular.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_DGESV
!=======================================================================
!> general solve: Computes the solution to system of linear equations A * X = B for GE matrices.
!! use LAPACK_GESV (nwtc_lapack::lapack_gesv) instead of this specific function.
   SUBROUTINE LAPACK_SGESV ( N, A, IPIV, B, ErrStat, ErrMsg )

      ! passed parameters

      INTEGER,         intent(in   ) :: N                 !< The number of linear equations, i.e., the order of the matrix A.  N >= 0.

      !     .. Array Arguments ..
      REAL(SiKi)      ,intent(inout) :: A( :, : )         !< On entry, the N-by-N coefficient matrix A.  On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored.
      REAL(SiKi)      ,intent(inout) :: B( :, : )         !< On entry, the N-by-NRHS matrix of right hand side matrix B.  On exit, if INFO = 0, the N-by-NRHS solution matrix X.
      INTEGER,         intent(  out) :: IPIV( : )         !< The pivot indices that define the permutation matrix P; row i of the matrix was interchanged with row IPIV(i).

      INTEGER(IntKi),  intent(  out) :: ErrStat           !< Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            !< Message describing error

         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value; > 0: if INFO = i, U(i,i) is exactly zero.
                                                          ! > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed, but the factor U is exactly singular, so the solution could not be computed.
      INTEGER                        :: NRHS              ! The number of right hand sides, i.e., the number of columns of the matrix B.  NRHS >= 0.
      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,N).
      INTEGER                        :: LDB               ! The leading dimension of the array B.  LDB >= max(1,N).
                        
      
      LDA  = SIZE(A,1)
      LDB  = SIZE(B,1)
      NRHS = SIZE(B,2)


      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL SGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_SGESV: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_SGESV: U('//TRIM(ErrMsg)//','//TRIM(ErrMsg)//')=0. Factor U is exactly singular.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_SGESV
!=======================================================================
!> general matrix factorization: Factor matrix into A=PLU.
!! use LAPACK_GETRF (nwtc_lapack::lapack_getrf) instead of this specific function.
   SUBROUTINE LAPACK_DGETRF( M, N, A, IPIV, ErrStat, ErrMsg )

      ! passed parameters

      INTEGER,         intent(in   ) :: M                 !< The number of rows of the matrix A.  M >= 0.
      INTEGER,         intent(in   ) :: N                 !< The number of columns of the matrix A.  N >= 0.

      !     .. Array Arguments ..
      REAL(R8Ki)      ,intent(inout) :: A( :, : )         !< On entry, the M-by-N matrix to be factored. On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored.
      INTEGER,         intent(  out) :: IPIV( : )         !< The pivot indices; for 1 <= i <= min(M,N), row i of the matrix was interchanged with row IPIV(i).

      INTEGER(IntKi),  intent(  out) :: ErrStat           !< Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            !< Message describing error

         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value; > 0: if INFO = i, U(i,i) is exactly zero. The factor U is exactly singular.
      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,M).

      LDA  = SIZE(A,1)


      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL DGETRF( M, N, A, LDA, IPIV, INFO )

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_DGETRF: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_DGETRF: U('//TRIM(ErrMsg)//','//TRIM(ErrMsg)//')=0. Factor U is exactly singular.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_DGETRF
!=======================================================================
!> general matrix factorization: Factor matrix into A=PLU.
!! use LAPACK_GETRF (nwtc_lapack::lapack_getrf) instead of this specific function.
   SUBROUTINE LAPACK_SGETRF( M, N, A, IPIV, ErrStat, ErrMsg )


      ! passed parameters

      INTEGER,        intent(in   ) :: M                 !< The number of rows of the matrix A.  M >= 0.
      INTEGER,        intent(in   ) :: N                 !< The number of columns of the matrix A.  N >= 0.

      !     .. Array Arguments ..
      REAL(SiKi)     ,intent(inout) :: A( :, : )         !< On entry, the M-by-N matrix to be factored. On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored.
      INTEGER,        intent(  out) :: IPIV( : )         !< The pivot indices; for 1 <= i <= min(M,N), row i of the matrix was interchanged with row IPIV(i).

      INTEGER(IntKi), intent(  out) :: ErrStat           !< Error level
      CHARACTER(*),   intent(  out) :: ErrMsg            !< Message describing error

         ! local variables
      INTEGER                       :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value; > 0: if INFO = i, U(i,i) is exactly zero. The factor U is exactly singular.
      INTEGER                       :: LDA               ! The leading dimension of the array A.  LDA >= max(1,M).

      LDA  = SIZE(A,1)



      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL SGETRF( M, N, A, LDA, IPIV, INFO )

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_SGETRF: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_SGETRF: U('//TRIM(ErrMsg)//','//TRIM(ErrMsg)//')=0. Factor U is exactly singular.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_SGETRF
!=======================================================================
!> general solve of factorized matrix: Solve system of linear equations Ax=PLUx=b.
!! use LAPACK_GETRS (nwtc_lapack::lapack_getrs) instead of this specific function.
   SUBROUTINE LAPACK_DGETRS( TRANS, N, A, IPIV, B, ErrStat, ErrMsg )

      ! passed parameters

      CHARACTER(1),    intent(in   ) :: TRANS             !< Specifies the form of the system of equations: = 'N':  A * X = B  (No transpose)
                                                          !!                                                = 'T':  A**T* X = B  (Transpose)
                                                          !!                                                = 'C':  A**T* X = B  (Conjugate transpose = Transpose)
      INTEGER,         intent(in   ) :: N                 !< The order of the matrix A.  N >= 0.

      !     .. Array Arguments ..
      INTEGER,         intent(in   ) :: IPIV( : )         !< The pivot indices from DGETRF (nwtc_lapack::lapack_getrf); for 1<=i<=N, row i of the matrix was interchanged with row IPIV(i).
      REAL(R8Ki)      ,intent(in   ) :: A( :, : )         !< The factors L and U from the factorization A = P*L*U as computed by DGETRF.
      REAL(R8Ki)      ,intent(inout) :: B( :, : )         !< On entry, the right hand side matrix B. On exit, the solution matrix X.

      INTEGER(IntKi),  intent(  out) :: ErrStat           !< Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            !< Message describing error

         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value;
      INTEGER                        :: NRHS              ! The number of right hand sides, i.e., the number of columns of the matrix B.  NRHS >= 0.
      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,N).
      INTEGER                        :: LDB               ! The leading dimension of the array B.  LDB >= max(1,N).
                        
      
      LDA  = SIZE(A,1)
      LDB  = SIZE(B,1)
      NRHS = SIZE(B,2)


      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) -INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_DGETRS: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_DGETRS: unknown error -'//TRIM(ErrMsg)//'.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_DGETRS
!=======================================================================
!> general solve of factorized matrix: Solve system of linear equations Ax=PLUx=b.
!! use LAPACK_GETRS (nwtc_lapack::lapack_getrs) instead of this specific function.
   SUBROUTINE LAPACK_DGETRS1( TRANS, N, A, IPIV, B, ErrStat, ErrMsg )

      ! passed parameters

      CHARACTER(1),    intent(in   ) :: TRANS             !< Specifies the form of the system of equations: = 'N':  A * X = B  (No transpose)
                                                          !!                                                = 'T':  A**T* X = B  (Transpose)
                                                          !!                                                = 'C':  A**T* X = B  (Conjugate transpose = Transpose)
      INTEGER,         intent(in   ) :: N                 !< The order of the matrix A.  N >= 0.

      !     .. Array Arguments ..
      INTEGER,         intent(in   ) :: IPIV( : )         !< The pivot indices from DGETRF (nwtc_lapack::lapack_getrf); for 1<=i<=N, row i of the matrix was interchanged with row IPIV(i).
      REAL(R8Ki)      ,intent(in   ) :: A( :, : )         !< The factors L and U from the factorization A = P*L*U as computed by DGETRF.
      REAL(R8Ki)      ,intent(inout) :: B( :    )         !< On entry, the right hand side matrix B. On exit, the solution matrix X.

      INTEGER(IntKi),  intent(  out) :: ErrStat           !< Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            !< Message describing error

         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value;
      INTEGER                        :: NRHS              ! The number of right hand sides, i.e., the number of columns of the matrix B.  NRHS >= 0.
      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,N).
      INTEGER                        :: LDB               ! The leading dimension of the array B.  LDB >= max(1,N).
                        
      
      LDA  = SIZE(A,1)
      LDB  = SIZE(B,1)
      NRHS = 1


      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) -INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_DGETRS1: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_DGETRS1: unknown error -'//TRIM(ErrMsg)//'.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_DGETRS1
!=======================================================================
!> general solve of factorized matrix: Solve system of linear equations Ax=PLUx=b.
!! use LAPACK_GETRS (nwtc_lapack::lapack_getrs) instead of this specific function.
   SUBROUTINE LAPACK_SGETRS( TRANS, N, A, IPIV, B, ErrStat, ErrMsg )

      ! passed parameters

      CHARACTER(1),   intent(in   )  :: TRANS             !< Specifies the form of the system of equations: = 'N':  A * X = B  (No transpose)
                                                          !!                                                = 'T':  A**T* X = B  (Transpose)
                                                          !!                                                = 'C':  A**T* X = B  (Conjugate transpose = Transpose)
      INTEGER,        intent(in   )  :: N                 !< The order of the matrix A.  N >= 0.
                                     
      !     .. Array Arguments ..    
      INTEGER,        intent(in   )  :: IPIV( : )         !< The pivot indices from DGETRF; for 1<=i<=N, row i of the matrix was interchanged with row IPIV(i).
      REAL(SiKi),     intent(in   )  :: A( :, : )         !< The factors L and U from the factorization A = P*L*U as computed by SGETRF.
      REAL(SiKi),     intent(inout)  :: B( :, : )         !< On entry, the right hand side matrix B. On exit, the solution matrix X.
                                     
      INTEGER(IntKi), intent(  out)  :: ErrStat           !< Error level
      CHARACTER(*),   intent(  out)  :: ErrMsg            !< Message describing error

         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value;
      INTEGER                        :: NRHS              ! The number of right hand sides, i.e., the number of columns of the matrix B.  NRHS >= 0.
      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,N).
      INTEGER                        :: LDB               ! The leading dimension of the array B.  LDB >= max(1,N).
                        
      
      LDA  = SIZE(A,1)
      LDB  = SIZE(B,1)
      NRHS = SIZE(B,2)


      ErrStat = ErrID_None
      ErrMsg  = ""

      !IF (ReKi == C_FLOAT) THEN
         CALL SGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      !ELSEIF (ReKi == C_DOUBLE) THEN
      !   CALL DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      !ELSE
      !   ErrStat = ErrID_FATAL
      !   ErrMsg  = "LAPACK_SGETRS: Matrix A is an invalid type."
      !   RETURN
      !END IF

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) -INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_SGETRS: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_SGETRS: unknown error -'//TRIM(ErrMsg)//'.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_SGETRS
!=======================================================================
!> general solve of factorized matrix: Solve system of linear equations Ax=PLUx=b.
!! use LAPACK_GETRS (nwtc_lapack::lapack_getrs) instead of this specific function.
   SUBROUTINE LAPACK_SGETRS1( TRANS, N, A, IPIV, B, ErrStat, ErrMsg )

      ! passed parameters

      CHARACTER(1),   intent(in   )  :: TRANS             !< Specifies the form of the system of equations: = 'N':  A * X = B  (No transpose)
                                                          !!                                                = 'T':  A**T* X = B  (Transpose)
                                                          !!                                                = 'C':  A**T* X = B  (Conjugate transpose = Transpose)
      INTEGER,        intent(in   )  :: N                 !< The order of the matrix A.  N >= 0.
                                     
      !     .. Array Arguments ..    
      INTEGER,        intent(in   )  :: IPIV( : )         !< The pivot indices from DGETRF; for 1<=i<=N, row i of the matrix was interchanged with row IPIV(i).
      REAL(SiKi),     intent(in   )  :: A( :, : )         !< The factors L and U from the factorization A = P*L*U as computed by SGETRF.
      REAL(SiKi),     intent(inout)  :: B( :    )         !< On entry, the right hand side matrix B. On exit, the solution matrix X.
                                     
      INTEGER(IntKi), intent(  out)  :: ErrStat           !< Error level
      CHARACTER(*),   intent(  out)  :: ErrMsg            !< Message describing error
                                     
         ! local variables           
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value;
      INTEGER                        :: NRHS              ! The number of right hand sides, i.e., the number of columns of the matrix B.  NRHS >= 0.
      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,N).
      INTEGER                        :: LDB               ! The leading dimension of the array B.  LDB >= max(1,N).
                        
      
      LDA  = SIZE(A,1)
      LDB  = SIZE(B,1)
      NRHS = 1


      ErrStat = ErrID_None
      ErrMsg  = ""

      
      !IF (ReKi == C_FLOAT) THEN
         CALL SGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      !ELSEIF (ReKi == C_DOUBLE) THEN
      !   CALL DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      !ELSE
      !   ErrStat = ErrID_FATAL
      !   ErrMsg  = "LAPACK_SGETRS: Matrix A is an invalid type."
      !   RETURN
      !END IF

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) -INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_SGETRS1: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_SGETRS1: unknown error -'//TRIM(ErrMsg)//'.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_SGETRS1
!=======================================================================
!> Compute the inverse of a general matrix using the LU factorization.
!! use LAPACK_GETRI (nwtc_lapack::lapack_getri) instead of this specific function.
   SUBROUTINE LAPACK_DGETRI( N, A, IPIV, WORK, LWORK, ErrStat, ErrMsg )

      ! passed parameters

      INTEGER,         intent(in   ) :: N                 !< The order of the matrix A.  N >= 0.
      INTEGER,         intent(in   ) :: LWORK             !< The dimension of the array WORK. LWORK >= max(1,N). For optimal performance LWORK >= N*NB, where NB is the optimal blocksize returned by ILAENV.
                                                          !! If LWORK = -1, then a workspace query is assumed; the routine only calculates the optimal size of the WORK array, returns this value as the first
                                                          !! entry of the WORK array, and no error message related to LWORK is issued by XERBLA.

      !     .. Array Arguments ..
      INTEGER,         intent(in   ) :: IPIV( : )         !< dimension (N). The pivot indices from DGETRF; for 1<=i<=N, row i of the matrix was interchanged with row IPIV(i).
      REAL(R8Ki)      ,intent(inout) :: A( :, : )         !< On entry, the factors L and U from the factorization A = P*L*U as computed by DGETRF. On exit, if INFO = 0, the inverse of the original matrix A.
      REAL(R8Ki)      ,intent(inout) :: WORK( : )         !< On exit, if INFO=0, then WORK(1) returns the optimal LWORK.

      INTEGER(IntKi),  intent(  out) :: ErrStat           !< Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            !< Message describing error

         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value; > 0: if INFO = i, U(i,i) is exactly zero. The matrix is singular and its inverse could not be computed.
      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,M).

      LDA  = SIZE(A,1)

      ErrStat = ErrID_None
      ErrMsg  = ""

      !IF (DbKi == C_DOUBLE) THEN
         CALL DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
      !ELSEIF (DbKi == C_FLOAT) THEN
      !   CALL DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
      !ELSE
      !   ErrStat = ErrID_FATAL
      !   ErrMsg  = "LAPACK_DGETRI: Matrix A is an invalid type."
      !   RETURN
      !END IF

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_DGETRI: illegal value in argument "//TRIM(ErrMsg)//"."
        ELSE
            ErrMsg = 'LAPACK_DGETRI: U('//TRIM(ErrMsg)//','//TRIM(ErrMsg)//')=0. Matrix is singular and its inverse cannot be computed.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_DGETRI
!=======================================================================
!> Compute the inverse of a general matrix using the LU factorization.
!! use LAPACK_GETRI (nwtc_lapack::lapack_getri) instead of this specific function.
   SUBROUTINE LAPACK_SGETRI( N, A, IPIV, WORK, LWORK, ErrStat, ErrMsg )

      ! passed parameters

      INTEGER,        intent(in   )  :: N                 !< The order of the matrix A.  N >= 0.
      INTEGER,        intent(in   )  :: LWORK             !< The dimension of the array WORK. LWORK >= max(1,N). For optimal performance LWORK >= N*NB, where NB is the optimal blocksize returned by ILAENV.
                                                          !! If LWORK = -1, then a workspace query is assumed; the routine only calculates the optimal size of the WORK array, returns this value as the first
                                                          !! entry of the WORK array, and no error message related to LWORK is issued by XERBLA.
                                     
      !     .. Array Arguments ..    
      INTEGER,        intent(in   )  :: IPIV( : )         !< dimension (N). The pivot indices from DGETRF; for 1<=i<=N, row i of the matrix was interchanged with row IPIV(i).
      REAL(SiKi),     intent(inout)  :: A( :, : )         !< On entry, the factors L and U from the factorization A = P*L*U as computed by SGETRF. On exit, if INFO = 0, the inverse of the original matrix A.
      REAL(SiKi),     intent(inout)  :: WORK( : )         !< On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
                                     
      INTEGER(IntKi), intent(  out)  :: ErrStat           !< Error level
      CHARACTER(*),   intent(  out)  :: ErrMsg            !< Message describing error
                                     
         ! local variables           
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value; > 0: if INFO = i, U(i,i) is exactly zero. The matrix is singular and its inverse could not be computed.
      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,M).

      LDA  = SIZE(A,1)

      ErrStat = ErrID_None
      ErrMsg  = ""


      CALL SGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_SGETRI: illegal value in argument "//TRIM(ErrMsg)//"."
        ELSE
            ErrMsg = 'LAPACK_SGETRI: U('//TRIM(ErrMsg)//','//TRIM(ErrMsg)//')=0. Matrix is singular and its inverse cannot be computed.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_SGETRI
!=======================================================================
!> Compute generalized eigenvalues and/or eigenvectors for a pair of N-by-N real nonsymmetric matrices (A,B).
!! use LAPACK_GGEV (nwtc_lapack::lapack_ggev) instead of this specific function.
   SUBROUTINE LAPACK_DGGEV(JOBVL, JOBVR, N, A, B, ALPHAR, ALPHAI, BETA, VL, VR, WORK, LWORK, ErrStat, ErrMsg)

      ! passed variables/parameters:

      CHARACTER(1),    intent(in   ) :: JOBVL             !< = 'N':  do not compute the left generalized eigenvectors; = 'V':  compute the left generalized eigenvectors.
      CHARACTER(1),    intent(in   ) :: JOBVR             !< = 'N':  do not compute the right generalized eigenvectors; = 'V':  compute the right generalized eigenvectors.

      INTEGER,         intent(in   ) :: N                 !< The order of the matrices A, B, VL, and VR.  N >= 0.

      INTEGER,         intent(in   ) :: LWORK             !< The dimension of the array WORK.  LWORK >= max(1,8*N). For good performance, LWORK must generally be larger.
                                                          !!   If LWORK = -1, then a workspace query is assumed; the routine only calculates the optimal size of the WORK array, returns
                                                          !!   this value as the first entry of the WORK array, and no error message related to LWORK is issued by XERBLA.


      REAL(R8Ki)      ,intent(inout) :: A( :, : )         !< dimension (LDA, N). On entry, the matrix A in the pair (A,B). On exit, A has been overwritten.
      REAL(R8Ki)      ,intent(inout) :: B( :, : )         !< dimension (LDB, N). On entry, the matrix B in the pair (A,B). On exit, B has been overwritten.

      REAL(R8Ki)      ,intent(  out) :: ALPHAR( : )       !< dimension (N). See comments for variable "Beta"
      REAL(R8Ki)      ,intent(  out) :: ALPHAI( : )       !< dimension (N). See comments for variable "Beta".
      REAL(R8Ki)      ,intent(  out) :: BETA( : )         !< On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will be the generalized eigenvalues.  If ALPHAI(j) is zero, then
                                                          !!   the j-th eigenvalue is real; if positive, then the j-th and (j+1)-st eigenvalues are a complex conjugate pair, with
                                                          !!   ALPHAI(j+1) negative.
                                                          !!
                                                          !!   Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j) may easily over- or underflow, and BETA(j) may even be zero.
                                                          !!   Thus, the user should avoid naively computing the ratio alpha/beta.  However, ALPHAR and ALPHAI will be always less
                                                          !!   than and usually comparable with norm(A) in magnitude, and BETA always less than and usually comparable with norm(B).


      REAL(R8Ki)      ,intent(  out) :: VL( :, : )        !< dimension (LDVL,N). If JOBVL = 'V', the left eigenvectors u(j) are stored one after another in the columns of VL, in the same
                                                          !!   order as their eigenvalues. If the j-th eigenvalue is real, then u(j) = VL(:,j), the j-th column of VL. If the j-th and
                                                          !!   (j+1)-th eigenvalues form a complex conjugate pair, then u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).
                                                          !!   Each eigenvector is scaled so the largest component has abs(real part)+abs(imag. part)=1. Not referenced if JOBVL = 'N'.
      REAL(R8Ki)      ,intent(  out) :: VR( :, : )        !< dimension (LDVR,N). If JOBVR = 'V', the right eigenvectors v(j) are stored one after another in the columns of VR, in the same
                                                          !!   order as their eigenvalues. If the j-th eigenvalue is real, then v(j) = VR(:,j), the j-th column of VR. If the j-th and
                                                          !!   (j+1)-th eigenvalues form a complex conjugate pair, then v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).
                                                          !!   Each eigenvector is scaled so the largest component has abs(real part)+abs(imag. part)=1. Not referenced if JOBVR = 'N'.

      REAL(R8Ki)      ,intent(inout) :: WORK( : )         !< dimension (MAX(1,LWORK)). On exit, if INFO = 0, WORK(1) returns the optimal LWORK.



      INTEGER(IntKi),  intent(  out) :: ErrStat           !< Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            !< Message describing error


         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit;
                                                          ! < 0:
                                                          !   = -i, the i-th argument had an illegal value;
                                                          ! > 0:
                                                          !   = 1,...,N: The QZ iteration failed.  No eigenvectors have been calculated, but ALPHAR(j), ALPHAI(j), and BETA(j) should be correct for j=INFO+1,...,N.
                                                          !   = N+1: other than QZ iteration failed in DHGEQZ.
                                                          !   = N+2: error return from DTGEVC.
      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,N).
      INTEGER                        :: LDB               ! The leading dimension of the array B.  LDB >= max(1,N).
      INTEGER                        :: LDVL              ! The leading dimension of the matrix VL. LDVL >= 1, and if JOBVL = 'V', LDVL >= N
      INTEGER                        :: LDVR              ! The leading dimension of the matrix VR. LDVR >= 1, and if JOBVR = 'V', LDVR >= N.
      CHARACTER(20)                  :: n_str
                        
      
      LDA  = SIZE(A,1)
      LDB  = SIZE(B,1)

      LDVL  = SIZE(VL,1)
      LDVR  = SIZE(VR,1)
      
      

      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_DGGEV: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSEIF (INFO <N) THEN
            !ErrStat = ErrID_Severe
            WRITE( ErrMsg, * ) INFO + 1
            WRITE( n_str, * ) n
            ErrMsg  = "LAPACK_DGGEV: The QZ iteration failed. No eigenvectors have been calculated, but ALPHAR(j), ALPHAI(j), and BETA(j) should be correct for j="&
                       //TRIM(ErrMsg)//",...,"//TRIM(n_str)//"."
         ELSEIF (INFO == N ) THEN
            ErrMsg  = "LAPACK_DGGEV: The QZ iteration failed. No eigenvectors have been calculated."
         ELSEIF (INFO == N+1) THEN
            ErrMsg  = "LAPACK_DGGEV: other than QZ iteration failed in DHGEQZ."
         ELSEIF (INFO == N+2) THEN
            ErrMsg  = "LAPACK_DGGEV: error return from DTGEVC."
         ELSE
            ErrMsg = 'LAPACK_DGGEV: unknown error '//TRIM(ErrMsg)//'.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_DGGEV
!=======================================================================
!> Compute generalized eigenvalues and/or eigenvectors for a pair of N-by-N real nonsymmetric matrices (A,B).
!! use LAPACK_GGEV (nwtc_lapack::lapack_ggev) instead of this specific function.
   SUBROUTINE LAPACK_SGGEV(JOBVL, JOBVR, N, A, B, ALPHAR, ALPHAI, BETA, VL, VR, WORK, LWORK, ErrStat, ErrMsg)

      ! subroutine arguments

      CHARACTER(1),   intent(in   )  :: JOBVL             !< = 'N':  do not compute the left generalized eigenvectors; = 'V':  compute the left generalized eigenvectors.
      CHARACTER(1),   intent(in   )  :: JOBVR             !< = 'N':  do not compute the right generalized eigenvectors; = 'V':  compute the right generalized eigenvectors.
                                    
      INTEGER,        intent(in   )  :: N                 !< The order of the matrices A, B, VL, and VR.  N >= 0.
      INTEGER,        intent(in   )  :: LWORK             !< The dimension of the array WORK.  LWORK >= max(1,8*N). For good performance, LWORK must generally be larger.
                                                          !!   If LWORK = -1, then a workspace query is assumed; the routine only calculates the optimal size of the WORK array, returns
                                                          !!   this value as the first entry of the WORK array, and no error message related to LWORK is issued by XERBLA.
                                    
                                    
      REAL(SiKi),     intent(inout)  :: A( :, : )         !< dimension (LDA, N). On entry, the matrix A in the pair (A,B). On exit, A has been overwritten.
      REAL(SiKi),     intent(inout)  :: B( :, : )         !< dimension (LDB, N). On entry, the matrix B in the pair (A,B). On exit, B has been overwritten.
                                    
                                    
      REAL(SiKi),     intent(  out)  :: ALPHAR( : )       !< dimension (N). See comments for variable "Beta"
      REAL(SiKi),     intent(  out)  :: ALPHAI( : )       !< dimension (N). See comments for variable "Beta".
      REAL(SiKi),     intent(  out)  :: BETA( : )         !< On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will be the generalized eigenvalues.  If ALPHAI(j) is zero, then
                                                          !!   the j-th eigenvalue is real; if positive, then the j-th and (j+1)-st eigenvalues are a complex conjugate pair, with
                                                          !!   ALPHAI(j+1) negative.
                                                          !!
                                                          !!   Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j) may easily over- or underflow, and BETA(j) may even be zero.
                                                          !!   Thus, the user should avoid naively computing the ratio alpha/beta.  However, ALPHAR and ALPHAI will be always less
                                                          !!   than and usually comparable with norm(A) in magnitude, and BETA always less than and usually comparable with norm(B).


      REAL(SiKi),     intent(  out)  :: VL( :, : )        !< dimension (LDVL,N). If JOBVL = 'V', the left eigenvectors u(j) are stored one after another in the columns of VL, in the same
                                                          !!   order as their eigenvalues. If the j-th eigenvalue is real, then u(j) = VL(:,j), the j-th column of VL. If the j-th and
                                                          !!   (j+1)-th eigenvalues form a complex conjugate pair, then u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).
                                                          !!   Each eigenvector is scaled so the largest component has abs(real part)+abs(imag. part)=1. Not referenced if JOBVL = 'N'.
      REAL(SiKi),     intent(  out)  :: VR( :, : )        !< dimension (LDVR,N). If JOBVR = 'V', the right eigenvectors v(j) are stored one after another in the columns of VR, in the same
                                                          !!   order as their eigenvalues. If the j-th eigenvalue is real, then v(j) = VR(:,j), the j-th column of VR. If the j-th and
                                                          !!   (j+1)-th eigenvalues form a complex conjugate pair, then v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).
                                                          !!   Each eigenvector is scaled so the largest component has abs(real part)+abs(imag. part)=1. Not referenced if JOBVR = 'N'.
                                     
      REAL(SiKi),     intent(inout)  :: WORK( : )         !< dimension (MAX(1,LWORK)). On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
                                     
                                     
                                     
      INTEGER(IntKi), intent(  out)  :: ErrStat           !< Error level
      CHARACTER(*),   intent(  out)  :: ErrMsg            !< Message describing error
                                     
         ! local variables           
      INTEGER                        :: INFO              ! = 0:  successful exit;
                                                          ! < 0:
                                                          !   = -i, the i-th argument had an illegal value;
                                                          ! > 0:
                                                          !   = 1,...,N: The QZ iteration failed.  No eigenvectors have been calculated, but ALPHAR(j), ALPHAI(j), and BETA(j) should be correct for j=INFO+1,...,N.
                                                          !   = N+1: other than QZ iteration failed in SHGEQZ.
                                                          !   = N+2: error return from STGEVC.
      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,N).
      INTEGER                        :: LDB               ! The leading dimension of the array B.  LDB >= max(1,N).
      INTEGER                        :: LDVL              ! The leading dimension of the matrix VL. LDVL >= 1, and if JOBVL = 'V', LDVL >= N
      INTEGER                        :: LDVR              ! The leading dimension of the matrix VR. LDVR >= 1, and if JOBVR = 'V', LDVR >= N.
      CHARACTER(20)                  :: n_str
                        
      
      LDA  = SIZE(A,1)
      LDB  = SIZE(B,1)

      LDVL  = SIZE(VL,1)
      LDVR  = SIZE(VR,1)


      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL SGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_SGGEV: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSEIF (INFO <=N) THEN
            !ErrStat = ErrID_Severe
            WRITE( ErrMsg, * ) INFO + 1
            WRITE( n_str, * ) n
            ErrMsg  = "LAPACK_SGGEV: The QZ iteration failed. No eigenvectors have been calculated, but ALPHAR(j), ALPHAI(j), and BETA(j) should be correct for j="&
                       //TRIM(ErrMsg)//",...,"//TRIM(n_str)//"."
         ELSEIF (INFO == N ) THEN
            ErrMsg  = "LAPACK_SGGEV: The QZ iteration failed. No eigenvectors have been calculated."
         ELSEIF (INFO == N+1) THEN
            ErrMsg  = "LAPACK_SGGEV: other than QZ iteration failed in SHGEQZ."
         ELSEIF (INFO == N+2) THEN
            ErrMsg  = "LAPACK_SGGEV: error return from STGEVC."
         ELSE
            ErrMsg = 'LAPACK_SGGEV: unknown error '//TRIM(ErrMsg)//'.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_SGGEV
!=======================================================================
!> Compute the solution to system of linear equations A * X = B for PO (positive-definite) matrices.
!! use LAPACK_POSV (nwtc_lapack::lapack_posv) instead of this specific function.
   SUBROUTINE LAPACK_DPOSV (UPLO, N, NRHS, A, B, ErrStat, ErrMsg)


      ! passed parameters

      INTEGER,         intent(in   ) :: N                !< The number of linear equations, i.e., the order of the matrix A.  N >= 0.
      INTEGER,         intent(in   ) :: NRHS             !< The number of right hand sides, i.e., the number of columns of the matrix B.  NRHS >= 0.

      !     .. Array Arguments ..
      REAL(R8Ki)      ,intent(inout) :: A( :, : )        !< On entry, the symmetric matrix A.  If UPLO = 'U', the leading N-by-N upper triangular part of A contains the upper
                                                         !! triangular part of the matrix A, and the strictly lower triangular part of A is not referenced.  If UPLO = 'L', the
                                                         !! leading N-by-N lower triangular part of A contains the lower triangular part of the matrix A, and the strictly upper
                                                         !! triangular part of A is not referenced.
                                                         !! On exit, if INFO = 0, the factor U or L from the Cholesky factorization A = U**T*U or A = L*L**T.
      REAL(R8Ki)      ,intent(inout) :: B( :, : )        !< On entry, the N-by-NRHS matrix of right hand side matrix B.  On exit, if INFO = 0, the N-by-NRHS solution matrix X.

      INTEGER(IntKi),  intent(  out) :: ErrStat          !< Error level
      CHARACTER(*),    intent(  out) :: ErrMsg           !< Message describing error
      CHARACTER(1),    intent(in   ) :: UPLO             !< 'U':  Upper triangle of A is stored; 'L':  Lower triangle of A is stored.

         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value; > 0: if INFO = i, U(i,i) is exactly zero.
                                                          ! > 0: if INFO = i, the leading minor of order i of A is not positive definite, so the factorization could not be
                                                          ! completed, and the solution has not been computed.

      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,N).
      INTEGER                        :: LDB               ! The leading dimension of the array B.  LDB >= max(1,N).


      LDA  = SIZE(A,1)
      LDB  = SIZE(B,1)
                                                          
                                                          
      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL DPOSV (UPLO, N, NRHS, A, LDA, B, LDB, INFO)

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_DPOSV: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_DPOSV: Leading minor order '//TRIM(ErrMsg)//' of A is not positive definite, so factorization could not be completed, and the solution has not been computed.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_DPOSV
!=======================================================================
!> Compute the solution to system of linear equations A * X = B for PO (positive-definite) matrices.
!! use LAPACK_POSV (nwtc_lapack::lapack_posv) instead of this specific function.
   SUBROUTINE LAPACK_SPOSV (UPLO, N, NRHS, A, B, ErrStat, ErrMsg)


      ! passed parameters

      INTEGER,         intent(in   ) :: N                 !< The number of linear equations, i.e., the order of the matrix A.  N >= 0.
      INTEGER,         intent(in   ) :: NRHS              !< The number of right hand sides, i.e., the number of columns of the matrix B.  NRHS >= 0.

      !     .. Array Arguments ..
      REAL(SiKi)      ,intent(inout) :: A( :, : )         !< On entry, the symmetric matrix A.  If UPLO = 'U', the leading N-by-N upper triangular part of A contains the upper
                                                          !! triangular part of the matrix A, and the strictly lower triangular part of A is not referenced.  If UPLO = 'L', the
                                                          !! leading N-by-N lower triangular part of A contains the lower triangular part of the matrix A, and the strictly upper
                                                          !! triangular part of A is not referenced.
                                                          !! On exit, if INFO = 0, the factor U or L from the Cholesky factorization A = U**T*U or A = L*L**T.
      REAL(SiKi)      ,intent(inout) :: B( :, : )         !< On entry, the N-by-NRHS matrix of right hand side matrix B.  On exit, if INFO = 0, the N-by-NRHS solution matrix X.

      INTEGER(IntKi),  intent(  out) :: ErrStat           !< Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            !< Message describing error
      CHARACTER(1),    intent(in   ) :: UPLO              !< 'U':  Upper triangle of A is stored; 'L':  Lower triangle of A is stored.

         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value; > 0: if INFO = i, U(i,i) is exactly zero.
                                                          ! > 0: if INFO = i, the leading minor of order i of A is not positive definite, so the factorization could not be
                                                          ! completed, and the solution has not been computed.
      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,N).
      INTEGER                        :: LDB               ! The leading dimension of the array B.  LDB >= max(1,N).


      LDA  = SIZE(A,1)
      LDB  = SIZE(B,1)



      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL SPOSV (UPLO, N, NRHS, A, LDA, B, LDB, INFO)

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_SPOSV: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_SPOSV: Leading minor order '//TRIM(ErrMsg)//' of A is not positive definite, so factorization could not be completed, and the solution has not been computed.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_SPOSV
!=======================================================================
!> Compute the Cholesky factorization of a real symmetric positive definite matrix A stored in packed format.
!! use LAPACK_PPTRF (nwtc_lapack::lapack_pptrf) instead of this specific function.
   SUBROUTINE LAPACK_DPPTRF (UPLO, N, AP, ErrStat, ErrMsg)

   ! DPPTRF computes the Cholesky factorization of a real symmetric
   ! positive definite matrix A stored in packed format.
   !
   ! The factorization has the form
   !      A = U**T * U,  if UPLO = 'U', or
   !      A = L  * L**T,  if UPLO = 'L',
   ! where U is an upper triangular matrix and L is lower triangular.
   

      ! passed parameters

      INTEGER,         intent(in   ) :: N                 !< The order of the matrix A.  N >= 0.

      !     .. Array Arguments ..
      REAL(R8Ki)      ,intent(inout) :: AP( : )           !< AP is REAL array, dimension (N*(N+1)/2)
                                                          !! On entry, the upper or lower triangle of the symmetric matrix A, packed columnwise in a linear array.  The j-th column of A
                                                          !! is stored in the array AP as follows:
                                                          !!    if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
                                                          !!    if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
                                                          !! See below for further details.
                                                          !! On exit, if INFO = 0, the triangular factor U or L from the Cholesky factorization A = U**T*U or A = L*L**T, in the same storage format as A.
                                                          !!
                                                          !! Further details:      
                                                          !!   The packed storage scheme is illustrated by the following example
                                                          !!   when N = 4, UPLO = 'U':
                                                          !! 
                                                          !!   Two-dimensional storage of the symmetric matrix A:
                                                          !! 
                                                          !!      a11 a12 a13 a14
                                                          !!          a22 a23 a24
                                                          !!              a33 a34     (aij = aji)
                                                          !!                  a44
                                                          !! 
                                                          !!   Packed storage of the upper triangle of A:
                                                          !! 
                                                          !!   AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]
                                                                                                                    
                                                          
      INTEGER(IntKi),  intent(  out) :: ErrStat           !< Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            !< Message describing error
      CHARACTER(1),    intent(in   ) :: UPLO              !< 'U':  Upper triangle of A is stored; 'L':  Lower triangle of A is stored.
      
      

         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value; 
                                                          ! > 0:  if INFO = i, the leading minor of order i is not positive definite, and the factorization could not be completed.

                                                          
      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL DPPTRF (UPLO, N, AP, INFO)

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_DPPTRF: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_DPPTRF: Leading minor order '//TRIM(ErrMsg)//' of A is not positive definite, so Cholesky factorization could not be completed.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_DPPTRF   
!=======================================================================
!> Compute the Cholesky factorization of a real symmetric positive definite matrix A stored in packed format.
!! use LAPACK_PPTRF (nwtc_lapack::lapack_pptrf) instead of this specific function.
   SUBROUTINE LAPACK_SPPTRF (UPLO, N, AP, ErrStat, ErrMsg)

   ! SPPTRF computes the Cholesky factorization of a real symmetric
   ! positive definite matrix A stored in packed format.
   !
   ! The factorization has the form
   !      A = U**T * U,  if UPLO = 'U', or
   !      A = L  * L**T,  if UPLO = 'L',
   ! where U is an upper triangular matrix and L is lower triangular.
   

      ! passed parameters

      INTEGER,         intent(in   ) :: N                 !< The order of the matrix A.  N >= 0.

      !     .. Array Arguments ..
      REAL(SiKi)      ,intent(inout) :: AP( : )           !< AP is REAL array, dimension (N*(N+1)/2)
                                                          !! On entry, the upper or lower triangle of the symmetric matrix A, packed columnwise in a linear array.  The j-th column of A
                                                          !! is stored in the array AP as follows:
                                                          !!    if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
                                                          !!    if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
                                                          !! See LAPACK_DPPTRF for further details.
                                                          !! On exit, if INFO = 0, the triangular factor U or L from the Cholesky factorization A = U**T*U or A = L*L**T, in the same storage format as A.

      INTEGER(IntKi),  intent(  out) :: ErrStat           !< Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            !< Message describing error
      CHARACTER(1),    intent(in   ) :: UPLO              !< 'U':  Upper triangle of A is stored; 'L':  Lower triangle of A is stored.

         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value; 
                                                          ! > 0:  if INFO = i, the leading minor of order i is not positive definite, and the factorization could not be completed


      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL SPPTRF (UPLO, N, AP, INFO)

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_SPPTRF: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_SPPTRF: Leading minor order '//TRIM(ErrMsg)//' of A is not positive definite, so Cholesky factorization could not be completed.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_SPPTRF
!=======================================================================
!> Compute singular value decomposition (SVD) for a general matrix, A.
!! use LAPACK_DGESVD (nwtc_lapack::lapack_dgesvd) instead of this specific function.
   SUBROUTINE LAPACK_DGESVD(JOBU, JOBVT, M, N, A, S, U, VT, WORK, LWORK, ErrStat, ErrMsg)

      ! passed variables/parameters:

      CHARACTER(1),    intent(in   ) :: JOBU              !<  'A':  all M columns of U are returned in array U;
                                                          !!  'S':  the first min(m,n) columns of U (the left singular
                                                          !!        vectors) are returned in the array U;
                                                          !!  'O':  the first min(m,n) columns of U (the left singular
                                                          !!        vectors) are overwritten on the array A;
                                                          !!  'N':  no columns of U (no left singular vectors) are
                                                          !!        computed.
      CHARACTER(1),    intent(in   ) :: JOBVT             !<  'A':  all N rows of V^T are returned in the array VT;
                                                          !!  'S':  the first min(m,n) rows of V^T (the right singular
                                                          !!        vectors) are returned in the array VT;
                                                          !!  'O':  the first min(m,n) rows of V^T (the right singular
                                                          !!        vectors) are overwritten on the array A;
                                                          !!  'N':  no rows of V**T (no right singular vectors) are
                                                          !!        computed.

      INTEGER,         intent(in   ) :: M                 !< The number of rows of the input matrix A.  M >= 0.
      INTEGER,         intent(in   ) :: N                 !< The number of columns of the input matrix A.  N >= 0.

      REAL(R8Ki),      intent(inout) :: A( :, : )         !< A is DOUBLE PRECISION array, dimension (LDA,N)
                                                          !! On entry, the M-by-N matrix A.
                                                          !! On exit,
                                                          !! if JOBU = 'O',  A is overwritten with the first min(m,n)
                                                          !!                 columns of U (the left singular vectors,
                                                          !!                 stored columnwise);
                                                          !! if JOBVT = 'O', A is overwritten with the first min(m,n)
                                                          !!                 rows of V**T (the right singular vectors,
                                                          !!                 stored rowwise);
                                                          !! if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
                                                          !!                 are destroyed.

      REAL(R8Ki),      intent(  out) :: S( : )            !< S is DOUBLE PRECISION array, dimension (min(M,N))
                                                          !! The singular values of A, sorted so that S(i) >= S(i+1).
      REAL(R8Ki),      intent(  out) :: U( :, : )         !< U is DOUBLE PRECISION array, dimension (LDU,UCOL)
                                                          !! (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
                                                          !! If JOBU = 'A', U contains the M-by-M orthogonal matrix U;
                                                          !! if JOBU = 'S', U contains the first min(m,n) columns of U
                                                          !! (the left singular vectors, stored columnwise);
                                                          !! if JOBU = 'N' or 'O', U is not referenced.
      REAL(R8Ki),      intent(  out) :: VT( :, : )        !< VT is DOUBLE PRECISION array, dimension (LDVT,N)
                                                          !! If JOBVT = 'A', VT contains the N-by-N orthogonal matrix
                                                          !! V**T;
                                                          !! if JOBVT = 'S', VT contains the first min(m,n) rows of
                                                          !! V**T (the right singular vectors, stored rowwise);
                                                          !! if JOBVT = 'N' or 'O', VT is not referenced.

      REAL(R8Ki),      intent(inout) :: WORK( : )         !< dimension (MAX(1,LWORK)). On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

      INTEGER,         intent(in   ) :: LWORK             !< The dimension of the array WORK.
                                                          !! LWORK >= MAX(1,5*MIN(M,N)) for the paths (see comments inside code):
                                                          !!    - PATH 1  (M much larger than N, JOBU='N')
                                                          !!    - PATH 1t (N much larger than M, JOBVT='N')
                                                          !! LWORK >= MAX(1,3*MIN(M,N) + MAX(M,N),5*MIN(M,N)) for the other paths
                                                          !! For good performance, LWORK should generally be larger.
                                                          !! If LWORK = -1, then a workspace query is assumed; the routine
                                                          !! only calculates the optimal size of the WORK array, returns
                                                          !! this value as the first entry of the WORK array, and no error
                                                          !! message related to LWORK is issued by XERBLA.
      INTEGER(IntKi),  intent(  out) :: ErrStat           !< Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            !< Message describing error

      ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit.
                                                          ! < 0:  if INFO = -i, the i-th argument had an illegal value.
                                                          ! > 0:  if DBDSQR did not converge, INFO specifies how many
                                                          !       superdiagonals of an intermediate bidiagonal form B
                                                          !       did not converge to zero. See the description of WORK
                                                          !       above for details.
      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,M).
      INTEGER                        :: LDU               ! The leading dimension of the array U.  LDU >= 1; if JOBU = 'S' or 'A', LDU >= M.
      INTEGER                        :: LDVT              ! The leading dimension of the array VT.  LDVT >= 1; if
                                                          !! JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).

      LDA  = SIZE(A,1)
      LDU  = SIZE(U,1)
      LDVT  = SIZE(VT,1)

      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL DGESVD(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO)

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_DGESVD: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSEIF (INFO > 0) THEN
            !ErrStat = ErrID_Severe
            ErrMsg  = "DBDSQR did not converge, INFO specifies how many superdiagonals of an intermediate bidiagonal form B did not converge to zero"&
                       //TRIM(ErrMsg)//",...,"//"."
         ELSE
            ErrMsg = 'LAPACK_DGESVD: unknown error '//TRIM(ErrMsg)//'.'
         END IF
      END IF

   RETURN
   END SUBROUTINE LAPACK_DGESVD
!=======================================================================
!> Compute singular value decomposition (SVD) for a general matrix, A.
!! use LAPACK_SGESVD (nwtc_lapack::lapack_sgesvd) instead of this specific function.
   SUBROUTINE LAPACK_SGESVD(JOBU, JOBVT, M, N, A, S, U, VT, WORK, LWORK, ErrStat, ErrMsg)

      ! passed variables/parameters:

      CHARACTER(1),    intent(in   ) :: JOBU              !<  'A':  all M columns of U are returned in array U;
                                                          !!  'S':  the first min(m,n) columns of U (the left singular
                                                          !!        vectors) are returned in the array U;
                                                          !!  'O':  the first min(m,n) columns of U (the left singular
                                                          !!        vectors) are overwritten on the array A;
                                                          !!  'N':  no columns of U (no left singular vectors) are
                                                          !!        computed.
      CHARACTER(1),    intent(in   ) :: JOBVT             !<  'A':  all N rows of V^T are returned in the array VT;
                                                          !!  'S':  the first min(m,n) rows of V^T (the right singular
                                                          !!        vectors) are returned in the array VT;
                                                          !!  'O':  the first min(m,n) rows of V^T (the right singular
                                                          !!        vectors) are overwritten on the array A;
                                                          !!  'N':  no rows of V**T (no right singular vectors) are
                                                          !!        computed.

      INTEGER,         intent(in   ) :: M                 !< The number of rows of the input matrix A.  M >= 0.
      INTEGER,         intent(in   ) :: N                 !< The number of columns of the input matrix A.  N >= 0.

      REAL(SiKi),      intent(inout) :: A( :, : )         !< A is SINGLE PRECISION array, dimension (LDA,N)
                                                          !! On entry, the M-by-N matrix A.
                                                          !! On exit,
                                                          !! if JOBU = 'O',  A is overwritten with the first min(m,n)
                                                          !!                 columns of U (the left singular vectors,
                                                          !!                 stored columnwise);
                                                          !! if JOBVT = 'O', A is overwritten with the first min(m,n)
                                                          !!                 rows of V**T (the right singular vectors,
                                                          !!                 stored rowwise);
                                                          !! if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
                                                          !!                 are destroyed.

      REAL(SiKi),      intent(  out) :: S( : )            !< S is SINGLE PRECISION array, dimension (min(M,N))
                                                          !! The singular values of A, sorted so that S(i) >= S(i+1).
      REAL(SiKi),      intent(  out) :: U( :, : )         !< U is SINGLE PRECISION array, dimension (LDU,UCOL)
                                                          !! (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
                                                          !! If JOBU = 'A', U contains the M-by-M orthogonal matrix U;
                                                          !! if JOBU = 'S', U contains the first min(m,n) columns of U
                                                          !! (the left singular vectors, stored columnwise);
                                                          !! if JOBU = 'N' or 'O', U is not referenced.
      REAL(SiKi),      intent(  out) :: VT( :, : )        !< VT is SINGLE PRECISION array, dimension (LDVT,N)
                                                          !! If JOBVT = 'A', VT contains the N-by-N orthogonal matrix
                                                          !! V**T;
                                                          !! if JOBVT = 'S', VT contains the first min(m,n) rows of
                                                          !! V**T (the right singular vectors, stored rowwise);
                                                          !! if JOBVT = 'N' or 'O', VT is not referenced.

      REAL(SiKi),      intent(inout) :: WORK( : )         !< dimension (MAX(1,LWORK)). On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

      INTEGER,         intent(in   ) :: LWORK             !< The dimension of the array WORK.
                                                          !! LWORK >= MAX(1,5*MIN(M,N)) for the paths (see comments inside code):
                                                          !!    - PATH 1  (M much larger than N, JOBU='N')
                                                          !!    - PATH 1t (N much larger than M, JOBVT='N')
                                                          !! LWORK >= MAX(1,3*MIN(M,N) + MAX(M,N),5*MIN(M,N)) for the other paths
                                                          !! For good performance, LWORK should generally be larger.
                                                          !! If LWORK = -1, then a workspace query is assumed; the routine
                                                          !! only calculates the optimal size of the WORK array, returns
                                                          !! this value as the first entry of the WORK array, and no error
                                                          !! message related to LWORK is issued by XERBLA.
      INTEGER(IntKi),  intent(  out) :: ErrStat           !< Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            !< Message describing error

      ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit.
                                                          ! < 0:  if INFO = -i, the i-th argument had an illegal value.
                                                          ! > 0:  if SBDSQR did not converge, INFO specifies how many
                                                          !       superdiagonals of an intermediate bidiagonal form B
                                                          !       did not converge to zero. See the description of WORK
                                                          !       above for details.
      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,M).
      INTEGER                        :: LDU               ! The leading dimension of the array U.  LDU >= 1; if JOBU = 'S' or 'A', LDU >= M.
      INTEGER                        :: LDVT              ! The leading dimension of the array VT.  LDVT >= 1; if
                                                          !! JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).

      LDA  = SIZE(A,1)
      LDU  = SIZE(U,1)
      LDVT  = SIZE(VT,1)

      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL SGESVD(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO)

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_SGESVD: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSEIF (INFO > 0) THEN
            !ErrStat = ErrID_Severe
            ErrMsg  = "DBDSQR did not converge, INFO specifies how many superdiagonals of an intermediate bidiagonal form B did not converge to zero"&
                       //TRIM(ErrMsg)//",...,"//"."
         ELSE
            ErrMsg = 'LAPACK_SGESVD: unknown error '//TRIM(ErrMsg)//'.'
         END IF
      END IF

   RETURN
   END SUBROUTINE LAPACK_SGESVD
   !=======================================================================
   !INTERFACE LAPACK_TPTTR:
   !>  Unpack a by-column-packed array into a 2D matrix format
   !!  See documentation in  DTPTTR/STPTTR source code.
   SUBROUTINE LAPACK_DTPTTR( UPLO, N, AP, A, LDA, ErrStat, ErrMsg )
      CHARACTER(1),   intent(in   ) :: UPLO     !< = 'U': A is an upper triangular matrix; 'L': A is a lower triangular matrix
      INTEGER,        intent(in   ) :: N        !< The order of matrix A and AP.
      INTEGER,        intent(in)    :: LDA      !< The leading dimension of the matrix A. LDA ? max(1,N)
      INTEGER(IntKi), intent(out)   :: ErrStat  !< Error level
      CHARACTER(*),   intent(out)   :: ErrMsg   !< Message describing error
      REAL(R8Ki),     intent(in)    :: AP( : )  !< Packed array
      REAL(R8Ki),     intent(out)   :: A( :,: ) !< Unpacked array : Note AP(1)=A(1,1); AP(2)=A(1,2); AP(3)=A(2,2); AP(4)=A(1,3) etc. by column, upper triang
      INTEGER :: INFO ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value 
      ErrStat = ErrID_None
      ErrMsg  = ""
      CALL DTPTTR( UPLO, N, AP, A, LDA, INFO )                
      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_DTPTTR: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_DTPTTR: Unknown error '//TRIM(ErrMsg)//'.'
         END IF
      END IF      
      RETURN
   END SUBROUTINE LAPACK_DTPTTR
   !=======================================================================
   !>  Unpack a by-column-packed array into a 2D matrix format
   SUBROUTINE LAPACK_STPTTR( UPLO, N, AP, A, LDA, ErrStat, ErrMsg )
      CHARACTER(1),   intent(in   ) :: UPLO     !< = 'U': A is an upper triangular matrix; 'L': A is a lower triangular matrix
      INTEGER,        intent(in   ) :: N        !< The order of matrix A and AP.
      INTEGER,        intent(in)    :: LDA      !< The leading dimension of the matrix A. LDA ? max(1,N)
      INTEGER(IntKi), intent(out)   :: ErrStat  !< Error level
      CHARACTER(*),   intent(out)   :: ErrMsg   !< Message describing error
      REAL(SiKi),     intent(in)    :: AP( : )  !< Packed array
      REAL(SiKi),     intent(out)   :: A( :,: ) !< Unpacked array : Note AP(1)=A(1,1); AP(2)=A(1,2); AP(3)=A(2,2); AP(4)=A(1,3) etc. by column, upper triang
      INTEGER :: INFO ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value 
      ErrStat = ErrID_None
      ErrMsg  = ""
      CALL STPTTR( UPLO, N, AP, A, LDA, INFO )                
      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_STPTTR: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_STPTTR: Unknown error '//TRIM(ErrMsg)//'.'
         END IF
      END IF      
      RETURN
   END SUBROUTINE LAPACK_STPTTR
   !=======================================================================
END MODULE NWTC_LAPACK
