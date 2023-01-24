!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
!
!    This file is part of SubDyn.
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
!**********************************************************************************************************************************
!> Standalone tools for beam-based finite element method (FEM)
!! No dependency with SubDyn types and representation
MODULE FEM
  USE NWTC_Library
  IMPLICIT NONE

  INTEGER, PARAMETER  :: FEKi = R8Ki  ! Define the kind to be used for FEM
  INTEGER, PARAMETER  :: LaKi = R8Ki  ! Define the kind to be used for LaPack

  INTERFACE FINDLOCI ! In the future, use FINDLOC from intrinsic
     MODULE PROCEDURE FINDLOCI_R8Ki
     MODULE PROCEDURE FINDLOCI_IntKi
     MODULE PROCEDURE FINDLOCI_SiKi
  END INTERFACE

 
CONTAINS
!------------------------------------------------------------------------------------------------------
!> Return eigenvalues, Omega, and eigenvectors

SUBROUTINE EigenSolve(K, M, N, bCheckSingularity, EigVect, Omega, ErrStat, ErrMsg )
   USE NWTC_LAPACK, only: LAPACK_ggev
   INTEGER       ,          INTENT(IN   )    :: N             !< Number of degrees of freedom, size of M and K
   REAL(LaKi),              INTENT(INOUT)    :: K(N, N)       !< Stiffness matrix 
   REAL(LaKi),              INTENT(INOUT)    :: M(N, N)       !< Mass matrix 
   LOGICAL,                 INTENT(IN   )    :: bCheckSingularity                  ! If True, the solver will fail if rigid modes are present 
   REAL(LaKi),              INTENT(INOUT)    :: EigVect(N, N) !< Returned Eigenvectors
   REAL(LaKi),              INTENT(INOUT)    :: Omega(N)      !< Returned Eigenvalues
   INTEGER(IntKi),          INTENT(  OUT)    :: ErrStat       !< Error status of the operation
   CHARACTER(*),            INTENT(  OUT)    :: ErrMsg        !< Error message if ErrStat /= ErrID_None
   ! LOCALS         
   REAL(LaKi), ALLOCATABLE                   :: WORK (:),  VL(:,:), AlphaR(:), AlphaI(:), BETA(:) ! eigensolver variables
   INTEGER                                   :: i  
   INTEGER                                   :: LWORK                          !variables for the eigensolver
   INTEGER,    ALLOCATABLE                   :: KEY(:)
   INTEGER(IntKi)                            :: ErrStat2
   CHARACTER(ErrMsgLen)                      :: ErrMsg2
   REAL(LaKi) :: normA
   REAL(LaKi) :: Omega2(N)  !< Squared eigenvalues
   REAL(LaKi), parameter :: MAX_EIGENVALUE = HUGE(1.0_ReKi) ! To avoid overflow when switching to ReKi
      
   ErrStat = ErrID_None
   ErrMsg  = ''
         
   ! allocate working arrays and return arrays for the eigensolver
   LWORK=8*N + 16  !this is what the eigensolver wants  >> bjj: +16 because of MKL ?ggev documenation ( "lwork >= max(1, 8n+16) for real flavors"), though LAPACK documenation says 8n is fine
   !bjj: there seems to be a memory problem in *GGEV, so I'm making the WORK array larger to see if I can figure it out
   CALL AllocAry( Work,    LWORK, 'Work',   ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'EigenSolve') 
   CALL AllocAry( AlphaR,  N,     'AlphaR', ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'EigenSolve')
   CALL AllocAry( AlphaI,  N,     'AlphaI', ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'EigenSolve')
   CALL AllocAry( Beta,    N,     'Beta',   ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'EigenSolve')
   CALL AllocAry( VL,      N,  N, 'VL',     ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'EigenSolve')
   CALL AllocAry( KEY,     N,     'KEY',    ErrStat2, ErrMsg2 ); if(Failed()) return
    
   ! --- Eigenvalue  analysis
   ! note: SGGEV seems to have memory issues in certain cases. The eigenvalues seem to be okay, but the eigenvectors vary wildly with different compiling options.
   !       DGGEV seems to work better, so I'm making these variables LaKi (which is set to R8Ki for now)   - bjj 4/25/2014
   ! bjj: This comes from the LAPACK documentation:
   !   Note: the quotients AlphaR(j)/BETA(j) and AlphaI(j)/BETA(j) may easily over- or underflow, and BETA(j) may even be zero.
   !   Thus, the user should avoid naively computing the ratio Alpha/beta.  However, AlphaR and AlphaI will be always less
   !   than and usually comparable with norm(A) in magnitude, and BETA always less than and usually comparable with norm(B).    
   ! Omega2=AlphaR/BETA  !Note this may not be correct if AlphaI<>0 and/or BETA=0 TO INCLUDE ERROR CHECK, also they need to be sorted
   CALL  LAPACK_ggev('N','V',N ,K, M, AlphaR, AlphaI, Beta, VL, EigVect, WORK, LWORK, ErrStat2, ErrMsg2)
   if(Failed()) return

   ! --- Determinign and sorting eigen frequencies 
   Omega2(:) =0.0_LaKi
   DO I=1,N !Initialize the key and calculate Omega
      KEY(I)=I
      !Omega2(I) = AlphaR(I)/Beta(I)
      if ( EqualRealNos(real(Beta(I),ReKi),0.0_ReKi) ) then
         ! --- Beta =0 
         if (bCheckSingularity) call WrScr('[WARN] Large eigenvalue found, system may be ill-conditioned')
         Omega2(I) = MAX_EIGENVALUE
      elseif ( EqualRealNos(real(AlphaI(I),ReKi),0.0_ReKi) ) THEN
         ! --- Real Eigenvalues
         IF ( AlphaR(I)<0.0_LaKi ) THEN
            if ( (AlphaR(I)/Beta(I))<1e-6_LaKi ) then
               ! Tolerating very small negative eigenvalues
               if (bCheckSingularity) call WrScr('[INFO] Negative eigenvalue found with small norm (system may contain rigid body mode)')
               Omega2(I)=0.0_LaKi
            else
               if (bCheckSingularity) call WrScr('[WARN] Negative eigenvalue found, system may be ill-conditioned.')
               Omega2(I)=AlphaR(I)/Beta(I)
            endif
         else
            Omega2(I) = AlphaR(I)/Beta(I)
         endif
      else
         ! --- Complex Eigenvalues
         normA = sqrt(AlphaR(I)**2 + AlphaI(I)**2)
         if ( (normA/Beta(I))<1e-6_LaKi ) then
            ! Tolerating very small eigenvalues with imaginary part
            if (bCheckSingularity) call WrScr('[WARN] Complex eigenvalue found with small norm, approximating as 0')
            Omega2(I) = 0.0_LaKi
         elseif ( abs(AlphaR(I))>1e3_LaKi*abs(AlphaI(I)) ) then
            ! Tolerating very small imaginary part compared to real part... (not pretty)
            if (bCheckSingularity) call WrScr('[WARN] Complex eigenvalue found with small Im compare to Re')
            Omega2(I) = AlphaR(I)/Beta(I)
         else
            if (bCheckSingularity) call WrScr('[WARN] Complex eigenvalue found with large imaginary value)')
            Omega2(I) = MAX_EIGENVALUE
         endif
         !call Fatal('Complex eigenvalue found, system may be ill-conditioned'); return
      endif
      ! Capping to avoid overflow
      if (Omega2(I)> MAX_EIGENVALUE) then
         Omega2(I) = MAX_EIGENVALUE
      endif
   enddo  

   ! Sorting. LASRT has issues for double precision 64 bit on windows
   !CALL ScaLAPACK_LASRT('I',N,Omega2,KEY,ErrStat2,ErrMsg2); if(Failed()) return 
   CALL sort_in_place(Omega2,KEY)
    
   ! --- Sorting eigen vectors
   ! KEEP ME: scaling of the eigenvectors using generalized mass =identity criterion
   ! ALLOCATE(normcoeff(N,N), STAT = ErrStat )
   ! result1 = matmul(M,EigVect)
   ! result2 = matmul(transpose(EigVect),result1)
   ! normcoeff=sqrt(result2)  !This should be a diagonal matrix which contains the normalization factors
   ! normcoeff=sqrt(matmul(transpose(EigVect),matmul(M,EigVect)))  !This should be a diagonal matrix which contains the normalization factors
   VL=EigVect  !temporary storage for sorting EigVect
   DO I=1,N 
      !EigVect(:,I)=VL(:,KEY(I))/normcoeff(KEY(I),KEY(I))  !reordered and normalized
      EigVect(:,I)=VL(:,KEY(I))  !just reordered as Huimin had a normalization outside of this one
   ENDDO
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

   ! --- Return Omega (capped by huge(ReKi)) and check for singularity
   Omega(:) = 0.0_LaKi
   do I=1,N 
      if (EqualRealNos(real(Omega2(I),ReKi), 0.0_ReKi)) then  ! NOTE: may be necessary for some corner numerics
         Omega(i)=0.0_LaKi
         if (bCheckSingularity) then
            call Fatal('Zero eigenvalue found, system may contain rigid body mode'); return
         endif
      elseif (Omega2(I)>0) then 
         Omega(i)=sqrt(Omega2(I))
      else
         ! Negative eigenfrequency
         print*,'>>> Wrong eigenfrequency, Omega^2=',Omega2(I) ! <<< This should never happen
         Omega(i)= 0.0_LaKi 
         call Fatal('Negative eigenvalue found, system may be ill-conditioned'); return
      endif
   enddo

   CALL CleanupEigen()
   RETURN

CONTAINS
   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'EigenSolve') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUpEigen()
   END FUNCTION Failed

   SUBROUTINE Fatal(ErrMsg_in)
      character(len=*), intent(in) :: ErrMsg_in
      CALL SetErrStat(ErrID_Fatal, ErrMsg_in, ErrStat, ErrMsg, 'EigenSolve');
      CALL CleanUpEigen()
   END SUBROUTINE Fatal

   SUBROUTINE CleanupEigen()
      IF (ALLOCATED(Work)  ) DEALLOCATE(Work)
      IF (ALLOCATED(AlphaR)) DEALLOCATE(AlphaR)
      IF (ALLOCATED(AlphaI)) DEALLOCATE(AlphaI)
      IF (ALLOCATED(Beta)  ) DEALLOCATE(Beta)
      IF (ALLOCATED(VL)    ) DEALLOCATE(VL)
      IF (ALLOCATED(KEY)   ) DEALLOCATE(KEY)
   END SUBROUTINE CleanupEigen
  
END SUBROUTINE EigenSolve

pure subroutine sort_in_place(a,key)
   real(LaKi), intent(inout), dimension(:) :: a
   integer(IntKi), intent(inout), dimension(:) :: key
   integer(IntKi) :: tempI
   real(LaKi) :: temp
   integer(IntKi) :: i, j
   do i = 2, size(a)
      j = i - 1
      temp  = a(i)
      tempI = key(i)
      do while (j>=1 .and. a(j)>temp)
         a(j+1) = a(j)
         key(j+1) = key(j)
         j = j - 1
         if (j<1) then
            exit
         endif
      end do
      a(j+1)   = temp
      key(j+1) = tempI
   end do
end subroutine sort_in_place

!> Compute the determinant of a real matrix using an LU factorization
FUNCTION Determinant(A, ErrStat, ErrMsg) result(det)
   use NWTC_LAPACK, only: LAPACK_GETRF
   REAL(FEKi),      INTENT(IN   ) :: A(:, :) !< Input matrix, no side effect
   INTEGER(IntKi),  INTENT(  OUT) :: ErrStat !< Error status of the operation
   CHARACTER(*),    INTENT(  OUT) :: ErrMsg  !< Error message if ErrStat /= ErrID_None
   real(FEKi)              :: det !< May easily overflow
   integer(IntKi)          :: i
   integer                 :: n
   integer, allocatable    :: ipiv(:)
   real(FEKi), allocatable :: PLU(:,:)
   real(FEKi) :: ScaleVal

   n = size(A(1,:))
   allocate(PLU(n,n))
   allocate(ipiv(n))
   ScaleVal= 1.0_FEKi
   PLU = A/ScaleVal
   ! general matrix factorization: Factor matrix into A=PLU.
   call LAPACK_GETRF( n, n, PLU, ipiv, ErrStat, ErrMsg ) !call dgetrf(n, n, PLU, n, ipiv, info)
   if (ErrStat==ErrID_Fatal) then
      print*,'Error in getrf'
      det = 0
      deallocate(PLU)
      deallocate(ipiv)
      return
   endif
   ! PLU now contains the LU of the factorization A = PLU
   ! As L has unit diagonal entries, the determinant can be computed
   ! from the product of U's diagonal entries. Additional sign changes
   ! stemming from the permutations P have to be taken into account as well.
   det = 1.0_FEKi
   do i = 1,n
      if(ipiv(i) /= i) then  ! additional sign change
         det = -det*PLU(i,i)
      else
         det =  det*PLU(i,i)
      endif
   end do
   deallocate(PLU)
   deallocate(ipiv)
   IF ( EqualRealNos(real(det, ReKi), 0.0_ReKi) ) THEN
      print*,'Det is zero'
      return 
   else
      det = det*(ScaleVal**n)
   endif
END FUNCTION Determinant
!------------------------------------------------------------------------------------------------------
!> Create a chessboard-like matrix with `valBlack` on the "black" cases, starting with black at (1,1)
!! As a generalization, "black" values may be spaced every `nSpace` squares
!! For instance, blackVal=9, whiteVal=0, nSpace=2
!!  [9 0 0 9 0 0 9]
!!  [0 9 0 0 9 0 0]
!!  [0 0 9 0 0 9 0]
!! Diagonal values may be overriden by `diagVal`
!! Matrix M does not need to be square
subroutine ChessBoard(M, blackVal, whiteVal, nSpace, diagVal)
   real(ReKi), dimension(:,:), intent(  out) :: M        !< Output matrix
   real(ReKi),                 intent(in   ) :: blackVal !< value for black squares
   real(ReKi),                 intent(in   ) :: whiteVal !< value for white squre
   integer(IntKi), optional,   intent(in   ) :: nSpace   !< spacing between black values, default 1
   real(ReKi), optional,       intent(in   ) :: diagVal  !< Value to override diagonal
   integer(IntKi) :: i, j, jFake, n
   ! Default value for spacing is 1 if not provided
   if (present(nSpace)) then; n=nSpace+1; else; n=2; endif
   ! Default values are white values
   M(:,:) = whiteVal
   ! Setting black values everyother n values
   do i=1,size(M,2)
      do jFake=1,size(M,2),n ! everyother n values
         j = mod(jFake+i-2, size(M,2)) +1
         !print*,'i,j',i,jFake,j
         M(i,j) = blackVal
      enddo
   enddo
   ! Forcing diagonal values
   if (present(diagVal)) then
      do i=1,size(M,1)
         do j=1,size(M,2) ! Matrix not necessarily square
            if (i==j) M(i,i) = diagVal
         enddo
      enddo
   endif
end subroutine ChessBoard
!------------------------------------------------------------------------------------------------------
!> Partition matrices and vectors into Boundary (R) and internal (L) nodes
!!  M = [ MRR, MRL ]
!!      [ sym, MLL ]
!! MRR = M(IDR, IDR),  KRR = M(IDR, IDR), FR = F(IDR)
!! MLL = M(IDL, IDL),  KRR = K(IDL, IDL), FL = F(IDL)
!! MRL = M(IDR, IDL),  KRR = K(IDR, IDL)
!! NOTE: generic code
SUBROUTINE BreakSysMtrx(MM, KK, IDR, IDL, nR, nL, MRR, MLL, MRL, KRR, KLL, KRL, FG, FGR, FGL, CC, CRR, CLL, CRL)
   REAL(FEKi),             INTENT(IN   )  :: MM(:,:)   !< Mass Matrix
   REAL(FEKi),             INTENT(IN   )  :: KK(:,:)   !< Stiffness matrix
   INTEGER(IntKi),         INTENT(IN   )  :: nR
   INTEGER(IntKi),         INTENT(IN   )  :: nL
   INTEGER(IntKi),         INTENT(IN   )  :: IDR(nR)   !< Indices of leader DOFs
   INTEGER(IntKi),         INTENT(IN   )  :: IDL(nL)   !< Indices of interior DOFs
   REAL(FEKi),             INTENT(  OUT)  :: MRR(nR, nR)
   REAL(FEKi),             INTENT(  OUT)  :: MLL(nL, nL) 
   REAL(FEKi),             INTENT(  OUT)  :: MRL(nR, nL)
   REAL(FEKi),             INTENT(  OUT)  :: KRR(nR, nR)
   REAL(FEKi),             INTENT(  OUT)  :: KLL(nL, nL)
   REAL(FEKi),             INTENT(  OUT)  :: KRL(nR, nL)
   REAL(FEKi), OPTIONAL,   INTENT(IN   )  :: FG(:)     !< Force vector
   REAL(FEKi), OPTIONAL,   INTENT(  OUT)  :: FGR(nR)
   REAL(FEKi), OPTIONAL,   INTENT(  OUT)  :: FGL(nL)
   REAL(FEKi), OPTIONAL,   INTENT(IN   )  :: CC(:,:)   !< Stiffness matrix
   REAL(FEKi), OPTIONAL,   INTENT(  OUT)  :: CRR(nR, nR)
   REAL(FEKi), OPTIONAL,   INTENT(  OUT)  :: CLL(nL, nL)
   REAL(FEKi), OPTIONAL,   INTENT(  OUT)  :: CRL(nR, nL)
   INTEGER(IntKi) :: I, J, II, JJ

   ! RR: Leader/Boundary DOFs
   DO I = 1, nR 
      II = IDR(I)
      DO J = 1, nR
         JJ = IDR(J)
         MRR(I, J) = MM(II, JJ)
         KRR(I, J) = KK(II, JJ)
      ENDDO
   ENDDO
   ! LL: Interior/follower DOFs
   DO I = 1, nL
      II = IDL(I)
      DO J = 1, nL
         JJ = IDL(J)
         MLL(I, J) = MM(II, JJ)
         KLL(I, J) = KK(II, JJ)
      ENDDO
   ENDDO
   ! RL: cross terms
   DO I = 1, nR 
      II = IDR(I)
      DO J = 1, nL
         JJ = IDL(J)
         MRL(I, J) = MM(II, JJ)
         KRL(I, J) = KK(II, JJ) 
      ENDDO 
   ENDDO
   ! Forces
   if (present(FG)) then
      if (present(FGR)) then
         do I = 1, nR 
            II = IDR(I)
            FGR(I) = FG(II)
         enddo
      endif
      if (present(FGL)) then
         do I = 1, nL
            II = IDL(I)
            FGL(I) = FG(II)
         enddo
      endif
   endif
   if (present(CC)) then
      ! RR: Leader/Boundary DOFs
      DO I = 1, nR 
         II = IDR(I)
         DO J = 1, nR
            JJ = IDR(J)
            CRR(I, J) = CC(II, JJ)
         ENDDO
      ENDDO
      ! LL: Interior/follower DOFs
      DO I = 1, nL
         II = IDL(I)
         DO J = 1, nL
            JJ = IDL(J)
            CLL(I, J) = CC(II, JJ)
         ENDDO
      ENDDO
      ! RL: cross terms
      DO I = 1, nR 
         II = IDR(I)
         DO J = 1, nL
            JJ = IDL(J)
            CRL(I, J) = CC(II, JJ) 
         ENDDO 
      ENDDO
   endif
END SUBROUTINE BreakSysMtrx

!------------------------------------------------------------------------------------------------------
!> Performs Craig-Bampton reduction of M and K matrices and optional Force vector
!! TODO: (Damping  optional)
!! Convention is: 
!!    "R": leader DOF     ->    "B": reduced leader DOF
!!    "L": interior DOF   ->    "M": reduced interior DOF (CB-modes)
!! NOTE: 
!!    - M_MM = Identity and K_MM = Omega*2 hence these matrices are not returned
!!    - Possibility to get more CB modes using the input nM_Out>nM
!!
!! NOTE: generic code
SUBROUTINE CraigBamptonReduction(MM, KK, IDR, nR, IDL, nL, nM, nM_Out, MBB, MBM, KBB, PhiL, PhiR, OmegaL, ErrStat, ErrMsg, FG, FGR, FGL, FGB, FGM, CC, CBB, CBM, CMM) 
   use NWTC_LAPACK, only: LAPACK_GEMV
   REAL(FEKi),             INTENT(IN   ) :: MM(:, :) !< Mass matrix
   REAL(FEKi),             INTENT(IN   ) :: KK(:, :) !< Stiffness matrix
   INTEGER(IntKi),         INTENT(IN   ) :: nR
   INTEGER(IntKi),         INTENT(IN   ) :: IDR(nR)   !< Indices of leader DOFs
   INTEGER(IntKi),         INTENT(IN   ) :: nL
   INTEGER(IntKi),         INTENT(IN   ) :: IDL(nL)   !< Indices of interior DOFs
   INTEGER(IntKi),         INTENT(IN   ) :: nM        !< Number of CB modes
   INTEGER(IntKi),         INTENT(IN   ) :: nM_Out    !< Number of modes returned for PhiL & OmegaL
   REAL(FEKi),             INTENT(  OUT) :: MBB( nR, nR)     !< Reduced Guyan Mass Matrix
   REAL(FEKi),             INTENT(  OUT) :: KBB( nR, nR)     !< Reduced Guyan Stiffness matrix
   REAL(FEKi),             INTENT(  OUT) :: MBM( nR, nM)     !< Cross term
   REAL(FEKi),             INTENT(  OUT) :: PhiR(nL, nR)     !< Guyan Modes   
   REAL(FEKi),             INTENT(  OUT) :: PhiL(nL, nM_out) !< Craig-Bampton modes
   REAL(FEKi),             INTENT(  OUT) :: OmegaL(nM_out)   !< Eigenvalues 
   REAL(FEKi), OPTIONAL,   INTENT(IN   ) :: FG(:)        !< Force vector (typically a constant force, like gravity)
   REAL(FEKi), OPTIONAL,   INTENT(  OUT) :: FGR(nR)      !< Force vector partitioned for R DOFs (TODO remove me)
   REAL(FEKi), OPTIONAL,   INTENT(  OUT) :: FGL(nL)      !< Force vector partitioned for L DOFs (TODO somehow for Static improvment..)
   REAL(FEKi), OPTIONAL,   INTENT(  OUT) :: FGB(nR)      !< Force vector in Guyan modes = FR+PhiR^t FL
   REAL(FEKi), OPTIONAL,   INTENT(  OUT) :: FGM(nM)      !< Force vector in CB modes    =    PhiM^t FL
   REAL(FEKi), OPTIONAL,   INTENT(IN   ) :: CC(:, :)     !< Damping matrix
   REAL(FEKi), OPTIONAL,   INTENT(  OUT) :: CBB(nR, nR)  !< Guyan Damping matrix
   REAL(FEKi), OPTIONAL,   INTENT(  OUT) :: CBM(nR, nM)  !< Coupling Damping matrix
   REAL(FEKi), OPTIONAL,   INTENT(  OUT) :: CMM(nM, nM)  !< Craig-Bampton Damping matrix
   INTEGER(IntKi),         INTENT(  OUT) :: ErrStat     ! Error status of the operation
   CHARACTER(*),           INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   INTEGER(IntKi)                        :: ErrStat2                                                                    
   CHARACTER(ErrMsgLen)                  :: ErrMsg2
   CHARACTER(*), PARAMETER               :: RoutineName = 'CraigBamptonReduction_FromPartition'
   ! Partitioned variables 
   real(FEKi), allocatable :: MRR(:, :)
   real(FEKi), allocatable :: MLL(:, :)
   real(FEKi), allocatable :: MRL(:, :)
   real(FEKi), allocatable :: KRR(:, :)
   real(FEKi), allocatable :: KLL(:, :)
   real(FEKi), allocatable :: KRL(:, :)
   real(FEKi), allocatable :: CRR(:, :)
   real(FEKi), allocatable :: CRL(:, :)
   real(FEKi), allocatable :: CLL(:, :)
   ! --- Break system
   CALL AllocAry(MRR, nR, nR, 'matrix MRR', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)  
   CALL AllocAry(MLL, nL, nL, 'matrix MLL', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)  
   CALL AllocAry(MRL, nR, nL, 'matrix MRL', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)  
   CALL AllocAry(KRR, nR, nR, 'matrix KRR', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)  
   CALL AllocAry(KLL, nL, nL, 'matrix KLL', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)  
   CALL AllocAry(KRL, nR, nL, 'matrix KRL', ErrStat2, ErrMsg2 ); if(Failed()) return
   if (present(CC)) then
      CALL AllocAry(CRR, nR, nR, 'matrix CRR', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)  
      CALL AllocAry(CLL, nL, nL, 'matrix CLL', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)  
      CALL AllocAry(CRL, nR, nL, 'matrix CRL', ErrStat2, ErrMsg2 ); if(Failed()) return
   endif
   call BreakSysMtrx(MM, KK, IDR, IDL, nR, nL, MRR, MLL, MRL, KRR, KLL, KRL, FG=FG, FGR=FGR, FGL=FGL, CC=CC, CRR=CRR, CLL=CLL, CRL=CRL)
   ! --- CB reduction
   call CraigBamptonReduction_FromPartition( MRR, MLL, MRL, KRR, KLL, KRL, nR, nL, nM, nM_Out,& !< Inputs 
                              MBB, MBM, KBB, PhiL, PhiR, OmegaL, ErrStat2, ErrMsg2, & !< Outputs
                              CRR=CRR, CLL=CLL, CRL=CRL,& !< Optional inputs
                              CBB=CBB, CBM=CBM, CMM=CMM)  !< Optional Outputs
   if(Failed()) return

   ! --- Reduction of force if provided
   if (present(FG).and.present(FGR).and.present(FGL)) then
      if (present(FGB)) then
         !FGB = FGR + matmul( transpose(PhiR), FGL)
         if (nL>0) then
            CALL LAPACK_GEMV('t', nL  , nR,  1.0_FeKi, PhiR, nL, FGL, 1, 0.0_FeKi, FGB, 1 )
            FGB = FGR + FGB
         else
            FGB = FGR 
         endif
      endif
      if (present(FGM)) then
         !FGM = matmul( FGL, PhiL(:,1:nM) ) != matmul( transpose(PhiM), FGL ) because FGL is 1-D
         if (nM>0) then
            CALL LAPACK_GEMV('t', nL  , nM,  1.0_FeKi, PhiL(:,1:nM), nL, FGL, 1, 0.0_FeKi, FGM, 1 )
         endif
      endif
   endif
   call CleanUp()

contains
   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CraigBamptonReduction') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   end function Failed
   subroutine CleanUp()
      IF(ALLOCATED(MRR)  ) DEALLOCATE(MRR) 
      IF(ALLOCATED(MLL)  ) DEALLOCATE(MLL) 
      IF(ALLOCATED(MRL)  ) DEALLOCATE(MRL) 
      IF(ALLOCATED(KRR)  ) DEALLOCATE(KRR) 
      IF(ALLOCATED(KLL)  ) DEALLOCATE(KLL) 
      IF(ALLOCATED(KRL)  ) DEALLOCATE(KRL) 
      IF(ALLOCATED(CRR)  ) DEALLOCATE(CRR) 
      IF(ALLOCATED(CLL)  ) DEALLOCATE(CLL) 
      IF(ALLOCATED(CRL)  ) DEALLOCATE(CRL) 
   end subroutine
END SUBROUTINE CraigBamptonReduction

!------------------------------------------------------------------------------------------------------
!> Performs Craig-Bampton reduction based on partitioned matrices M and K
!! Convention is: 
!!    "R": leader DOF     ->    "B": reduced leader DOF
!!    "L": interior DOF   ->    "M": reduced interior DOF (CB-modes)
!! NOTE: 
!!    - M_MM = Identity and K_MM = Omega*2 hence these matrices are not returned
!!    - Possibility to get more CB modes using the input nM_Out>nM (e.g. for static improvement)
!!
!! NOTE: generic code
SUBROUTINE CraigBamptonReduction_FromPartition( MRR, MLL, MRL, KRR, KLL, KRL, nR, nL, nM, nM_Out,&
                     MBB, MBM, KBB, PhiL, PhiR, OmegaL, ErrStat, ErrMsg,&
                     CRR, CLL, CRL, CBB, CBM, CMM)
   USE NWTC_LAPACK, only: LAPACK_getrs, LAPACK_getrf, LAPACK_gemm
   INTEGER(IntKi),         INTENT(  in)  :: nR
   INTEGER(IntKi),         INTENT(  in)  :: nL
   INTEGER(IntKi),         INTENT(  in)  :: nM_Out
   INTEGER(IntKi),         INTENT(  in)  :: nM
   REAL(FEKi),             INTENT(  IN)  :: MRR( nR, nR) !< Partitioned mass and stiffness matrices
   REAL(FEKi),             INTENT(  IN)  :: MLL( nL, nL) 
   REAL(FEKi),             INTENT(  IN)  :: MRL( nR, nL)
   REAL(FEKi),             INTENT(  IN)  :: KRR( nR, nR)
   REAL(FEKi),             INTENT(INOUT) :: KLL( nL, nL)  ! on exit, it has been factored (otherwise not changed)
   REAL(FEKi),             INTENT(  IN)  :: KRL( nR, nL)
   REAL(FEKi),             INTENT(  OUT) :: MBB( nR, nR)
   REAL(FEKi),             INTENT(  OUT) :: MBM( nR, nM)
   REAL(FEKi),             INTENT(  OUT) :: KBB( nR, nR)
   REAL(FEKi),             INTENT(  OUT) :: PhiR(nL, nR)     !< Guyan Modes   
   REAL(FEKi),             INTENT(  OUT) :: PhiL(nL, nM_Out) !< Craig-Bampton modes
   REAL(FEKi),             INTENT(  OUT) :: OmegaL(nM_Out)   !< Eigenvalues
   INTEGER(IntKi),         INTENT(  OUT) :: ErrStat          !< Error status of the operation
   CHARACTER(*),           INTENT(  OUT) :: ErrMsg           !< Error message if ErrStat /= ErrID_None
   REAL(FEKi), OPTIONAL,   INTENT(  IN)  :: CRR( nR, nR) !< Partitioned damping matrices
   REAL(FEKi), OPTIONAL,   INTENT(  IN)  :: CLL( nL, nL) 
   REAL(FEKi), OPTIONAL,   INTENT(  IN)  :: CRL( nR, nL)
   REAL(FEKi), OPTIONAL,   INTENT(  OUT) :: CBB( nR, nR)  !< Guyan damping matrix
   REAL(FEKi), OPTIONAL,   INTENT(  OUT) :: CBM( nR, nM)  !< Coupling damping matrix
   REAL(FEKi), OPTIONAL,   INTENT(  OUT) :: CMM( nM, nM)  !< CB damping matrix
   ! LOCAL VARIABLES
   REAL(FEKi) , allocatable :: Mu(:, :)          ! matrix for normalization Mu(p%nDOFL, p%nDOFL) [bjj: made allocatable to try to avoid stack issues]
   REAL(FEKi) , allocatable :: Temp(:, :)        ! temp matrix for intermediate steps [bjj: made allocatable to try to avoid stack issues]
   REAL(FEKi) , allocatable :: PhiR_T_MLL(:,:)   ! PhiR_T_MLL(nR,nL) = transpose of PhiR * MLL (temporary storage)
   INTEGER                  :: I        !counter
   INTEGER                  :: ipiv(nL) ! length min(m,n) (See LAPACK documentation)
   INTEGER(IntKi)           :: ErrStat2
   CHARACTER(ErrMsgLen)     :: ErrMsg2
   CHARACTER(*), PARAMETER  :: RoutineName = 'CraigBamptonReduction_FromPartition'
   ErrStat = ErrID_None 
   ErrMsg  = ''
   
   if (nM_out>nL) then
      ErrMsg2='Cannot request more modes than internal degrees of Freedom'; ErrStat2=ErrID_Fatal; 
      if(Failed()) return;
   endif
   if (nM_out<nM) then
      ErrMsg2='Cannot request more output modes than modes.'; ErrStat2=ErrID_Fatal; 
      if(Failed()) return;
   endif
   
   ! --- Compute CB modes (PhiL) and eigenvalues (OmegaL)
   if ( nM_out > 0 ) then 
      ! bCheckSingularity = True
      CALL EigenSolveWrap(KLL, MLL, nL, nM_out, .True., PhiL(:,1:nM_out), OmegaL(1:nM_out),  ErrStat2, ErrMsg2); if(Failed()) return
      ! --- Normalize PhiL
      ! MU = MATMUL ( MATMUL( TRANSPOSE(PhiL), MLL ), PhiL )
      CALL AllocAry( Temp , nM_out, nL     , 'Temp' , ErrStat2 , ErrMsg2); if(Failed()) return
      CALL AllocAry( MU   , nM_out, nM_out , 'Mu'   , ErrStat2 , ErrMsg2); if(Failed()) return
      CALL LAPACK_gemm( 'T', 'N', 1.0_FeKi, PhiL, MLL, 0.0_FeKi, Temp  , ErrStat2, ErrMsg2); if(Failed()) return
      CALL LAPACK_gemm( 'N', 'N', 1.0_FeKi, Temp, PhiL, 0.0_FeKi, MU  , ErrStat2, ErrMsg2); if(Failed()) return
      DEALLOCATE(Temp)
      ! PhiL = MATMUL( PhiL, MU ) ! this is the normalization (MU is diagonal)   
      DO I = 1, nM_out
         PhiL(:,I) = PhiL(:,I) / SQRT( MU(I, I) )
      ENDDO    
      DEALLOCATE(MU)
      if (present(CRR)) then
         ! CB damping CMM = PhiL^T  CLL PhiL
         CALL AllocAry( Temp , nM, nL    , 'Temp' , ErrStat2 , ErrMsg2); if(Failed()) return
         CALL LAPACK_gemm( 'T', 'N', 1.0_FeKi, PhiL(1:nL, 1:nM), CLL,  0.0_FeKi, Temp , ErrStat2, ErrMsg2); if(Failed()) return
         CALL LAPACK_gemm( 'N', 'N', 1.0_FeKi, Temp, PhiL(1:nL, 1:nM), 0.0_FeKi, CMM  , ErrStat2, ErrMsg2); if(Failed()) return
         DEALLOCATE(Temp)
      endif
   else
      PhiL   = 0.0_FEKi
      OmegaL = 0.0_FEKi
      if (present(CRR)) CMM  = 0.0_FEKi
   end if

   if (nL>0) then
      ! --- Compute Guyan Modes (PhiR)
      ! factor KLL to compute PhiR: KLL*PhiR=-TRANSPOSE(KRL)
      ! ** note this must be done after EigenSolveWrap() because it modifies KLL **
      CALL LAPACK_getrf( nL, nL, KLL, ipiv, ErrStat2, ErrMsg2); if(Failed()) return
      
      PhiR = -1.0_FEKi * TRANSPOSE(KRL) !set "b" in Ax=b  (solve KLL * PhiR = - TRANSPOSE( KRL ) for PhiR)
      CALL LAPACK_getrs( TRANS='N', N=nL, A=KLL, IPIV=ipiv, B=PhiR, ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
      
      ! --- Set MBB, MBM, and KBB from Eq. 4:
      CALL AllocAry( PhiR_T_MLL,  nR, nL, 'PhiR_T_MLL', ErrStat2, ErrMsg2); if(Failed()) return
      CALL AllocAry( Temp , nR, nR, 'Temp' , ErrStat2 , ErrMsg2); if(Failed()) return
         
      ! PhiR_T_MLL = TRANSPOSE(PhiR) * MLL
      CALL LAPACK_gemm( 'T', 'N', 1.0_FeKi, PhiR, MLL, 0.0_FeKi, PhiR_T_MLL  , ErrStat2, ErrMsg2); if(Failed()) return
      ! MBB1 = MATMUL(MRL, PhiR)
      CALL LAPACK_gemm( 'N', 'N', 1.0_FeKi, MRL, PhiR, 0.0_FeKi, MBB  , ErrStat2, ErrMsg2); if(Failed()) return
      ! MBB2 = MATMUL( PhiR_T_MLL, PhiR )
      CALL LAPACK_gemm( 'N', 'N', 1.0_FeKi, PhiR_T_MLL, PhiR, 0.0_FeKi, Temp  , ErrStat2, ErrMsg2); if(Failed()) return
      MBB = MRR + MBB + TRANSPOSE( MBB ) + Temp
      DEALLOCATE(Temp)
         
      IF ( nM == 0) THEN
         MBM = 0.0_FEKi
      ELSE
         CALL AllocAry( Temp , nR, nM, 'Temp' , ErrStat2 , ErrMsg2); if(Failed()) return
         !MBM = MATMUL( PhiR_T_MLL, PhiL(:,1:nM))  ! last half of operation
         CALL LAPACK_gemm( 'N', 'N', 1.0_FeKi, PhiR_T_MLL, PhiL(:,1:nM), 0.0_FeKi, MBM  , ErrStat2, ErrMsg2); if(Failed()) return
         ! Temp = MATMUL( MRL, PhiL(:,1:nM) )
         CALL LAPACK_gemm( 'N', 'N', 1.0_FeKi, MRL, PhiL(:,1:nM), 0.0_FeKi, Temp  , ErrStat2, ErrMsg2); if(Failed()) return
         MBM = Temp + MBM    !This had PhiM      
         DEALLOCATE(Temp)
      ENDIF
      
      !KBB = MATMUL(KRL, PhiR)   
      CALL LAPACK_gemm( 'N', 'N', 1.0_FeKi, KRL, PhiR, 0.0_FeKi, KBB  , ErrStat2, ErrMsg2); if(Failed()) return
      KBB = KBB + KRR

      if (present(CRR)) then
         ! Guyan damping CBB = CRR + (CRL*PhiR) + (CRL*PhiR)^T + PhiR^T*CLL*PhiR
         ! PhiR_T_CLL = TRANSPOSE(PhiR) * CLL
         CALL AllocAry( Temp , nR, nR, 'Temp' , ErrStat2 , ErrMsg2); if(Failed()) return
         CALL LAPACK_gemm( 'T', 'N', 1.0_FeKi, PhiR, MLL, 0.0_FeKi, PhiR_T_MLL  , ErrStat2, ErrMsg2); if(Failed()) return
         ! CBB = MATMUL(CRL, PhiR)
         CALL LAPACK_gemm( 'N', 'N', 1.0_FeKi, CRL, PhiR, 0.0_FeKi, CBB  , ErrStat2, ErrMsg2); if(Failed()) return
         ! CBB2 = MATMUL( PhiR_T_CLL, PhiR )
         CALL LAPACK_gemm( 'N', 'N', 1.0_FeKi, PhiR_T_MLL, PhiR, 0.0_FeKi, Temp  , ErrStat2, ErrMsg2); if(Failed()) return
         CBB = CRR + CBB + TRANSPOSE( CBB ) + Temp
         DEALLOCATE(Temp)
         ! Cross coupling CMB = PhiM^T*CLR + PhiM^T CLL PhiR
         !                CBM = CRL*PhiM + PhiR^T CLL^T PhiM (NOTE: assuming CLL symmetric)
         IF ( nM == 0) THEN
            CBM = 0.0_FEKi
            CMM = 0.0_FEKi
         ELSE
            CBM = MATMUL( PhiR_T_MLL, PhiL(:,1:nM))  ! last half of operation
            CBM = MATMUL( CRL, PhiL(:,1:nM) ) + CBM    !This had PhiM      
         ENDIF
      endif
   else
      PhiR(1:nL,1:nR) = 0.0_FEKi   ! Empty
      MBM (1:nR,1:nM) = 0.0_FEKi ! Empty
      MBB = MRR
      KBB = KRR
      if (present(CRR)) then
         CBB=CRR
         CBM=0.0_FEKi
      endif
   endif
        
   call CleanUp()
CONTAINS

   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CraigBamptonReduction_FromPartition') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   end function Failed
   
   subroutine CleanUp()
      if (allocated(Mu        )) DEALLOCATE(Mu        )
      if (allocated(Temp      )) DEALLOCATE(Temp      )
      if (allocated(PhiR_T_MLL)) DEALLOCATE(PhiR_T_MLL)
   end subroutine
END SUBROUTINE CraigBamptonReduction_FromPartition

!------------------------------------------------------------------------------------------------------
!> Wrapper function for eigen value analyses, for two cases:
!! Case1: K and M are taken "as is", this is used for the "LL" part of the matrix
!! Case2: K and M contain some constraints lines, and they need to be removed from the Mass/Stiffness matrix. Used for full system
SUBROUTINE EigenSolveWrap(K, M, nDOF, NOmega,  bCheckSingularity, EigVect, Omega, ErrStat, ErrMsg, bDOF )
   INTEGER,                INTENT(IN   )    :: nDOF                               ! Total degrees of freedom of the incoming system
   REAL(FEKi),             INTENT(IN   )    :: K(nDOF, nDOF)                      ! stiffness matrix 
   REAL(FEKi),             INTENT(IN   )    :: M(nDOF, nDOF)                      ! mass matrix 
   INTEGER,                INTENT(IN   )    :: NOmega                             ! No. of requested eigenvalues
   LOGICAL,                INTENT(IN   )    :: bCheckSingularity                  ! If True, the solver will fail if rigid modes are present 
   REAL(FEKi),             INTENT(  OUT)    :: EigVect(nDOF, NOmega)                  ! Returned Eigenvectors
   REAL(FEKi),             INTENT(  OUT)    :: Omega(NOmega)                      ! Returned Eigenvalues
   INTEGER(IntKi),         INTENT(  OUT)    :: ErrStat                            ! Error status of the operation
   CHARACTER(*),           INTENT(  OUT)    :: ErrMsg                             ! Error message if ErrStat /= ErrID_None
   LOGICAL,   OPTIONAL,    INTENT(IN   )    :: bDOF(nDOF)                         ! Optinal Mask for DOF to keep (True), or reduce (False)
   
   ! LOCALS         
   REAL(LaKi), ALLOCATABLE                   :: K_LaKi(:,:), M_LaKi(:,:) 
   REAL(LaKi), ALLOCATABLE                   :: EigVect_LaKi(:,:), Omega_LaKi(:) 
   INTEGER(IntKi)                            :: N
   INTEGER(IntKi)                            :: ErrStat2
   CHARACTER(ErrMsgLen)                      :: ErrMsg2
   ErrStat = ErrID_None
   ErrMsg  = ''
   EigVect=0.0_FeKi
   Omega=0.0_FeKi

   ! --- Unfortunate conversion to FEKi... TODO TODO consider storing M and K in FEKi
   if (present(bDOF)) then
      ! Remove unwanted DOFs
      call RemoveDOF(M, bDOF, M_LaKi, ErrStat2, ErrMsg2); if(Failed()) return
      call RemoveDOF(K, bDOF, K_LaKi, ErrStat2, ErrMsg2); if(Failed()) return
   else
      N=size(K,1)
      CALL AllocAry(K_LaKi      , N, N, 'K_FEKi',    ErrStat2, ErrMsg2); if(Failed()) return
      CALL AllocAry(M_LaKi      , N, N, 'M_FEKi',    ErrStat2, ErrMsg2); if(Failed()) return
      K_LaKi = real( K, LaKi )
      M_LaKi = real( M, LaKi )
   endif
   N=size(K_LaKi,1)

   ! Note:  NOmega must be <= N, which is the length of Omega2, Phi!
   if ( NOmega > nDOF ) then
      CALL SetErrStat(ErrID_Fatal,"NOmega must be less than or equal to N",ErrStat,ErrMsg,'EigenSolveWrap')
      CALL CleanupEigen()
      return
   end if

   ! --- Eigenvalue analysis
   CALL AllocAry(EigVect_LAKi, N, N, 'EigVect', ErrStat2, ErrMsg2); if(Failed()) return;
   CALL AllocAry(Omega_LaKi,     N , 'Omega', ErrStat2, ErrMsg2); if(Failed()) return; ! <<< NOTE: Needed due to dimension of Omega
   CALL EigenSolve(K_LaKi, M_LaKi, N, bCheckSingularity, EigVect_LaKi, Omega_LaKi, ErrStat2, ErrMsg2 ); if (Failed()) return;

   Omega(:)        = huge(1.0_ReKi)
   Omega(1:nOmega) = real(Omega_LaKi(1:nOmega), FEKi) !<<< nOmega<N

   ! --- Setting up Phi, and type conversion
   if (present(bDOF)) then
      ! Insert 0s where bDOF was false
      CALL InsertDOFRows(EigVect_LaKi(:,1:nOmega), bDOF, 0.0_FEKi, EigVect, ErrStat2, ErrMsg2 ); if(Failed()) return
   else
      EigVect=REAL( EigVect_LaKi(:,1:NOmega), LaKi )   ! eigenvectors
   endif
   CALL CleanupEigen()
   return
CONTAINS
   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'EigenSolveWrap') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUpEigen()
   END FUNCTION Failed

   SUBROUTINE CleanupEigen()
      IF (ALLOCATED(Omega_LaKi)  ) DEALLOCATE(Omega_LaKi) 
      IF (ALLOCATED(EigVect_LaKi)) DEALLOCATE(EigVect_LaKi)
      IF (ALLOCATED(K_LaKi)      ) DEALLOCATE(K_LaKi)
      IF (ALLOCATED(M_LaKi)      ) DEALLOCATE(M_LaKi)
   END SUBROUTINE CleanupEigen
  
END SUBROUTINE EigenSolveWrap
!------------------------------------------------------------------------------------------------------
!> Remove degrees of freedom from a matrix (lines and rows)
SUBROUTINE RemoveDOF(A, bDOF, Ared, ErrStat, ErrMsg )
   REAL(FEKi),             INTENT(IN   ) :: A(:, :)        ! full matrix
   logical,                INTENT(IN   ) :: bDOF(:)        ! Array of logical specifying whether a DOF is to be kept(True), or removed (False)
   REAL(LaKi),ALLOCATABLE, INTENT(  OUT) :: Ared(:,:)      ! reduced matrix
   INTEGER(IntKi),         INTENT(  OUT) :: ErrStat        ! Error status of the operation
   CHARACTER(*),           INTENT(  OUT) :: ErrMsg         ! Error message if ErrStat /= ErrID_None
   !locals
   INTEGER                               :: I, J           ! counters into full matrix
   INTEGER                               :: Ir, Jr         ! counters into reduced matrix
   INTEGER                               :: nr             ! number of reduced DOF
   ErrStat = ErrID_None
   ErrMsg  = ''    

   nr= count(bDOF)
   CALL AllocAry(Ared, nr, nr, 'Ared', ErrStat, ErrMsg ); if (ErrStat >= AbortErrLev) return

   ! Remove rows and columns from A when bDOF is 
   Jr=0
   do J = 1, size(A,1)
      if (bDOF(J)) then
         Jr=Jr+1
         Ir=0
         do I = 1, size(A,1)
            if (bDOF(I)) then
               Ir=Ir+1
               Ared(Ir, Jr) = REAL( A(I, J), FEKi )
            end if
         end do
      endif
   end do
END SUBROUTINE RemoveDOF

!> Expand a matrix to includes rows where bDOF is False (inverse behavior as RemoveDOF)
SUBROUTINE InsertDOFrows(Ared, bDOF, DefaultVal, A, ErrStat, ErrMsg )
   REAL(LaKi),             INTENT(IN   ) :: Ared(:, :)     ! Reduced matrix
   logical,                INTENT(IN   ) :: bDOF(:)        ! Array of logical specifying whether a DOF is to be kept(True), or removed (False)
   REAL(FEKi),             INTENT(IN   ) :: DefaultVal     ! Default value to fill the 
   REAL(FEKi)            , INTENT(INOUT) :: A(:,:)         ! Full matrix
   INTEGER(IntKi),         INTENT(  OUT) :: ErrStat        ! Error status of the operation
   CHARACTER(*),           INTENT(  OUT) :: ErrMsg         ! Error message if ErrStat /= ErrID_None
   !locals
   INTEGER                               :: I         ! counter into full matrix
   INTEGER                               :: Ir        ! counter into reduced matrix
   INTEGER                               :: n         ! number of DOF (fullsystem)
   ErrStat = ErrID_None
   ErrMsg  = ''    
   n= size(bDOF)
   IF ( size(Ared,1) > n) THEN
      ErrStat = ErrID_Fatal
      ErrMsg = 'InsertDOFrows: Number of reduced rows needs to be lower than full system rows'
      RETURN
   END IF
   IF ( size(Ared,2) /= size(A,2) ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg = 'InsertDOFrows: Inconsistent number of columns between A and Ared'
      RETURN
   END IF
   !CALL AllocAry(A, n, size(Ared,2), 'A', ErrStat, ErrMsg ); if (ErrStat >= AbortErrLev) return

   ! Use rows from Ared when bDOF is true, use default value otherwise
   ir=0 ! initialize 
   do i=1,n
      if (bDOF(i)) then
         ir =ir +1
         A(i,:)=Ared(ir,:)
      else
         A(i,:)=DefaultVal
      endif   
   enddo
END SUBROUTINE InsertDOFrows
!------------------------------------------------------------------------------------------------------
!> Returns index of val in Array (val is an integer!)
! NOTE: in the future use intrinsinc function findloc
FUNCTION FINDLOCI_R8Ki(Array, Val) result(i)
   real(R8Ki)    , dimension(:), intent(in) :: Array !< Array to search in
   integer(IntKi), intent(in)               :: val   !< Val
   integer(IntKi)                           :: i     !< Index of joint in joint table
   i = 1
   do while ( i <= size(Array) )
      if ( Val == NINT(Array(i)) ) THEN
         return ! Exit when found
      else
         i = i + 1
      endif
   enddo
   i=-1
END FUNCTION
!> Returns index of val in Array (val is an integer!)
! NOTE: in the future use intrinsinc function findloc
FUNCTION FINDLOCI_IntKi(Array, Val) result(i)
   integer(IntKi), dimension(:), intent(in) :: Array !< Array to search in
   integer(IntKi), intent(in)               :: val   !< Val
   integer(IntKi)                           :: i     !< Index of joint in joint table
   i = 1
   do while ( i <= size(Array) )
      if ( Val == Array(i) ) THEN
         return ! Exit when found
      else
         i = i + 1
      endif
   enddo
   i=-1
END FUNCTION

FUNCTION FINDLOCI_SiKi(Array, Val) result(i)
   real(SiKi), dimension(:), intent(in) :: Array !< Array to search in
   integer(IntKi), intent(in)               :: val   !< Val
   integer(IntKi)                           :: i     !< Index of joint in joint table
   i = 1
   do while ( i <= size(Array) )
      if ( Val == Array(i) ) THEN
         return ! Exit when found
      else
         i = i + 1
      endif
   enddo
   i=-1
END FUNCTION
!------------------------------------------------------------------------------------------------------
SUBROUTINE RigidTransformationLine(dx,dy,dz,iLine,Line)
   real(ReKi),               INTENT(IN)  :: dx,dy,dz
   integer(IntKi)     ,      INTENT(IN)  :: iLine 
   Real(ReKi), dimension(6), INTENT(OUT) :: Line
   SELECT CASE (iLine)
      CASE (1); Line = (/1.0_ReKi, 0.0_ReKi, 0.0_ReKi, 0.0_ReKi,       dz,      -dy/)
      CASE (2); Line = (/0.0_ReKi, 1.0_ReKi, 0.0_ReKi,      -dz, 0.0_ReKi,       dx/)
      CASE (3); Line = (/0.0_ReKi, 0.0_ReKi, 1.0_ReKi,       dy,      -dx, 0.0_ReKi/)
      CASE (4); Line = (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 1.0_ReKi, 0.0_ReKi, 0.0_ReKi/)
      CASE (5); Line = (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 1.0_ReKi, 0.0_ReKi/)
      CASE (6); Line = (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 1.0_ReKi/)
      CASE DEFAULT
         Line=-99999999_ReKi
         print*,'Error in RigidTransformationLine'
         STOP
!          ErrStat = ErrID_Fatal
!          ErrMsg  = 'Error calculating transformation matrix TI '
!          return
   END SELECT
END SUBROUTINE
!------------------------------------------------------------------------------------------------------
!> Rigid transformation matrix between DOFs of node j and k where node j is the leader node.
SUBROUTINE GetRigidTransformation(Pj, Pk, TRigid, ErrStat, ErrMsg)
   REAL(ReKi),       INTENT(IN   )  :: Pj(3)         ! (x,y,z) positions of leader node
   REAL(ReKi),       INTENT(IN   )  :: Pk(3)         ! (x,y,z) positions of follower node
   REAL(ReKi),       INTENT(  OUT)  :: TRigid(6,6)   ! Transformation matrix such that xk = T.xj
   INTEGER(IntKi),   INTENT(  OUT)  :: ErrStat       ! Error status of the operation
   CHARACTER(*),     INTENT(  OUT)  :: ErrMsg        ! Error message if ErrStat /= ErrID_None
   ! Local
   !REAL(ReKi) :: L             ! length of element
   !REAL(ReKi) :: DirCos(3, 3)  ! direction cosine matrix
   !REAL(ReKi) :: R0(3,3) 
   integer(IntKi) :: I
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! --- Formulation using Delta of Global coordinates
   Trigid=0; do I = 1,6; Trigid(I,I) = 1; enddo
   Trigid ( 1, 5 ) =  (Pk(3) - Pj(3))
   Trigid ( 1, 6 ) = -(Pk(2) - Pj(2))
   Trigid ( 2, 4 ) = -(Pk(3) - Pj(3))
   Trigid ( 2, 6 ) =  (Pk(1) - Pj(1))
   Trigid ( 3, 4 ) =  (Pk(2) - Pj(2))
   Trigid ( 3, 5 ) = -(Pk(1) - Pj(1))

   ! --- Formulation bty transforming the "local" matrix into a global one
   !call GetDirCos(Pj, Pk, R0, L, ErrStat, ErrMsg)
   !TRigid = 0 ; do I = 1,6; TRigid(I,I) = 1; enddo
   !TRigid (1, 5) =  L
   !TRigid (2, 4) = -L
   !TRigid(1:3,4:6) =  matmul( R0 , matmul(TRigid(1:3,4:6), transpose(R0)) )

   ! --- Formulation using L and Rotation matrix
   !TRigid = 0; do I = 1,6; TRigid(I,I) = 1; enddo
   !TRigid ( 1, 5 ) =  L*R0(3,3)
   !TRigid ( 1, 6 ) = -L*R0(2,3)
   !TRigid ( 2, 4 ) = -L*R0(3,3)
   !TRigid ( 2, 6 ) =  L*R0(1,3)
   !TRigid ( 3, 4 ) =  L*R0(2,3)
   !TRigid ( 3, 5 ) = -L*R0(1,3)
END SUBROUTINE GetRigidTransformation
!------------------------------------------------------------------------------------------------------
!> Computes directional cosine matrix DirCos
!! Transforms from element to global coordinates:  xg = DC.xe,  Kg = DC.Ke.DC^t
!! Assumes that the element main direction is along ze.
!!
!! bjj: note that this is the transpose of what is normally considered the Direction Cosine Matrix  
!!      in the FAST framework.
SUBROUTINE GetDirCos(P1, P2, DirCos, L_out, ErrStat, ErrMsg)
   REAL(ReKi) ,      INTENT(IN   )  :: P1(3), P2(3)      ! (x,y,z) global positions of two nodes making up an element
   REAL(FEKi) ,      INTENT(  OUT)  :: DirCos(3, 3)      ! calculated direction cosine matrix
   REAL(ReKi) ,      INTENT(  OUT)  :: L_out             ! length of element
   INTEGER(IntKi),   INTENT(  OUT)  :: ErrStat           ! Error status of the operation
   CHARACTER(*),     INTENT(  OUT)  :: ErrMsg            ! Error message if ErrStat /= ErrID_None
   REAL(FEKi)                       :: Dx, Dy, Dz, Dxy,L! distances between nodes
   ErrMsg  = ""
   ErrStat = ErrID_None
   
   Dx=P2(1)-P1(1)
   Dy=P2(2)-P1(2)
   Dz=P2(3)-P1(3)
   Dxy = sqrt( Dx**2 + Dy**2 )
   L   = sqrt( Dx**2 + Dy**2 + Dz**2)
   
   IF ( EqualRealNos(L, 0.0_FEKi) ) THEN
      ErrMsg = ' Same starting and ending location in the element.'
      ErrStat = ErrID_Fatal
      RETURN
   ENDIF
   
   IF ( EqualRealNos(Dxy, 0.0_FEKi) ) THEN 
      DirCos=0.0_FEKi    ! whole matrix set to 0
      IF ( Dz < 0) THEN  !x is kept along global x
         DirCos(1, 1) =  1.0_FEKi
         DirCos(2, 2) = -1.0_FEKi
         DirCos(3, 3) = -1.0_FEKi
      ELSE
         DirCos(1, 1) = 1.0_ReKi
         DirCos(2, 2) = 1.0_ReKi
         DirCos(3, 3) = 1.0_ReKi
      ENDIF 
   ELSE
      DirCos(1, 1) =  Dy/Dxy
      DirCos(1, 2) = +Dx*Dz/(L*Dxy)
      DirCos(1, 3) =  Dx/L
      
      DirCos(2, 1) = -Dx/Dxy
      DirCos(2, 2) = +Dz*Dy/(L*Dxy)
      DirCos(2, 3) =  Dy/L
     
      DirCos(3, 1) = 0.0_FEKi
      DirCos(3, 2) = -Dxy/L
      DirCos(3, 3) = +Dz/L
   ENDIF
   L_out= real(L, ReKi)

END SUBROUTINE GetDirCos
!------------------------------------------------------------------------------------------------------
!> Returns two vectors orthonormal to the input vector
SUBROUTINE GetOrthVectors(e1, e2, e3, ErrStat, ErrMsg)
   real(ReKi) ,      intent(in   )  :: e1(3) !<
   real(ReKi) ,      intent(  out)  :: e2(3) !<
   real(ReKi) ,      intent(  out)  :: e3(3) !<
   integer(IntKi),   intent(  out)  :: ErrStat           ! error status of the operation
   character(*),     intent(  out)  :: ErrMsg            ! error message if errstat /= errid_none
   real(ReKi) :: min_norm
   real(ReKi) :: e2_norm
   real(ReKi) :: e1b(3)
   ErrMsg  = ""
   ErrStat = ErrID_None

   min_norm = min( abs(e1(1)),  abs(e1(2)), abs(e1(3)) )
   ! Finding a good candidate for orthogonality
   if      (min_norm == abs(e1(1))) then; e2 = (/ 0._ReKi, -e1(3),  e1(2) /)
   else if (min_norm == abs(e1(2))) then; e2 = (/ e1(3)  , 0._ReKi, -e1(1) /)
   else if (min_norm == abs(e1(3))) then; e2 = (/-e1(2)  ,  e1(1), 0._ReKi /)
   endif
   ! Normalizing
   e2_norm=sqrt(e2(1)**2 + e2(2)**2 + e2(3)**2)
   if (abs(e2_norm)<1e-8) then
      ErrStat=ErrID_Fatal
      ErrMsg='Failed to determine orthogonal vector'
      e2=-99999._ReKi
      e3=-99999._ReKi
      return
   endif
   e2 = e2/e2_norm
   e1b= e1/sqrt(e1(1)**2 + e1(2)**2 + e1(3)**2)
   ! Third 
   e3 =  cross_product(e1b,e2)
END SUBROUTINE GetOrthVectors
!------------------------------------------------------------------------------------------------------
!> Element stiffness matrix for classical beam elements
!! shear is true  -- non-tapered Timoshenko beam 
!! shear is false -- non-tapered Euler-Bernoulli beam 
SUBROUTINE ElemK_Beam(A, L, Ixx, Iyy, Jzz, Shear, kappa_x, kappa_y, E, G, DirCos, K)
   REAL(ReKi), INTENT( IN) :: A, L, Ixx, Iyy, Jzz, E, G, kappa_x, kappa_y
   REAL(FEKi), INTENT( IN) :: DirCos(3,3) !< From element to global: xg = DC.xe,  Kg = DC.Ke.DC^t
   LOGICAL   , INTENT( IN) :: Shear
   REAL(FEKi), INTENT(OUT) :: K(12, 12) 
   ! Local variables
   REAL(FEKi)                            :: Ax, Ay, Kx, Ky
   REAL(FEKi)                            :: DC(12, 12)
   
   Ax = kappa_x*A
   Ay = kappa_y*A
   
   K(1:12,1:12) = 0.0_FEKi
   
   IF (Shear) THEN
      Kx = 12.0_FEKi*E*Iyy / (G*Ax*L*L)
      Ky = 12.0_FEKi*E*Ixx / (G*Ay*L*L)
   ELSE
      Kx = 0.0_FEKi
      Ky = 0.0_FEKi
   ENDIF
      
   K( 9,  9) = E*A/L
   K( 7,  7) = 12.0_FEKi*E*Iyy/( L*L*L*(1.0_FEKi + Kx) )
   K( 8,  8) = 12.0_FEKi*E*Ixx/( L*L*L*(1.0_FEKi + Ky) )
   K(12, 12) = G*Jzz/L
   K(10, 10) = (4.0_FEKi + Ky)*E*Ixx / ( L*(1.0_FEKi+Ky) )  
   K(11, 11) = (4.0_FEKi + Kx)*E*Iyy / ( L*(1.0_FEKi+Kx) )
   K( 2,  4) = -6._FEKi*E*Ixx / ( L*L*(1.0_FEKi+Ky) )
   K( 1,  5) =  6._FEKi*E*Iyy / ( L*L*(1.0_FEKi+Kx) )
   K( 4, 10) = (2.0_FEKi-Ky)*E*Ixx / ( L*(1.0_FEKi+Ky) )
   K( 5, 11) = (2.0_FEKi-Kx)*E*Iyy / ( L*(1.0_FEKi+Kx) )
   
   K( 3,  3)  = K(9,9)
   K( 1,  1)  = K(7,7)
   K( 2,  2)  = K(8,8)
   K( 6,  6)  = K(12,12)
   K( 4,  4)  = K(10,10)
   K(5,5)  = K(11,11)
   K(4,2)  = K(2,4)
   K(5,1)  = K(1,5)
   K(10,4) = K(4,10)
   K(11,5) = K(5,11)
   K(12,6)= -K(6,6)
   K(10,2)=  K(4,2)
   K(11,1)=  K(5,1)
   K(9,3) = -K(3,3)
   K(7,1) = -K(1,1)
   K(8,2) = -K(2,2)
   K(6, 12) = -K(6,6)
   K(2, 10) =  K(4,2)
   K(1, 11) =  K(5,1)
   K(3, 9)  = -K(3,3)
   K(1, 7)  = -K(1,1)
   K(2, 8)  = -K(2,2)
   K(11,7) = -K(5,1)
   K(10,8) = -K(4,2)
   K(7,11) = -K(5,1)
   K(8,10) = -K(4,2)
   K(7,5) = -K(5,1)
   K(5,7) = -K(5,1)
   K(8,4) = -K(4,2)
   K(4,8) = -K(4,2)
   
   DC = 0.0_FEKi
   DC( 1: 3,  1: 3) = DirCos
   DC( 4: 6,  4: 6) = DirCos
   DC( 7: 9,  7: 9) = DirCos
   DC(10:12, 10:12) = DirCos
   
   K = MATMUL( MATMUL(DC, K), TRANSPOSE(DC) ) ! TODO: change me if DirCos convention is  transposed
   
END SUBROUTINE ElemK_Beam
!------------------------------------------------------------------------------------------------------
!> Element stiffness matrix for pretension cable
!! Element coordinate system:  z along the cable!
SUBROUTINE ElemK_Cable(A, L, E, T0, DirCos, K)
   REAL(ReKi), INTENT( IN) :: A, L, E
   REAL(ReKi), INTENT( IN) :: T0 ! Pretension [N]
   REAL(FEKi), INTENT( IN) :: DirCos(3,3) !< From element to global: xg = DC.xe,  Kg = DC.Ke.DC^t
   REAL(FEKi), INTENT(OUT) :: K(12, 12) 
   ! Local variables
   REAL(FEKi) :: L0, Eps0, EAL0, EE
   REAL(FEKi) :: DC(12, 12)

   Eps0 = T0/(E*A)
   L0   = L/(1+Eps0)  ! "rest length" for which pretension would be 0
   EAL0 = E*A/L0
   EE   = EAL0* Eps0/(1+Eps0)

   K(1:12,1:12)=0.0_FEKi

   ! Note: only translational DOF involved (1-3, 7-9)
   K(1,1)= EE
   K(2,2)= EE
   K(3,3)= EAL0

   K(1,7)= -EE
   K(2,8)= -EE
   K(3,9)= -EAL0

   K(7,1)= -EE
   K(8,2)= -EE
   K(9,3)= -EAL0

   K(7,7)= EE
   K(8,8)= EE
   K(9,9)= EAL0


   DC = 0.0_FEKi
   DC( 1: 3,  1: 3) = DirCos
   DC( 4: 6,  4: 6) = DirCos
   DC( 7: 9,  7: 9) = DirCos
   DC(10:12, 10:12) = DirCos
   
   K = MATMUL( MATMUL(DC, K), TRANSPOSE(DC) ) ! TODO: change me if DirCos convention is  transposed
END SUBROUTINE ElemK_Cable
!------------------------------------------------------------------------------------------------------
!> Element mass matrix for classical beam elements
SUBROUTINE ElemM_Beam(A, L, Ixx, Iyy, Jzz, rho, DirCos, M)
   REAL(ReKi), INTENT( IN) :: A, L, Ixx, Iyy, Jzz, rho
   REAL(FEKi), INTENT( IN) :: DirCos(3,3) !< From element to global: xg = DC.xe,  Kg = DC.Ke.DC^t
   REAL(FEKi), INTENT(OUT) :: M(12, 12)

   REAL(FEKi) :: t, rx, ry, po
   REAL(FEKi) :: DC(12, 12)
   
   t = rho*A*L;
   rx = rho*Ixx;
   ry = rho*Iyy;
   po = rho*Jzz*L;   

   M(1:12,1:12) = 0.0_FEKi
      
   M( 9,  9) = t/3.0_FEKi
   M( 7,  7) = 13.0_FEKi*t/35.0_FEKi + 6.0_FEKi*ry/(5.0_FEKi*L)
   M( 8,  8) = 13.0_FEKi*t/35.0_FEKi + 6.0_FEKi*rx/(5.0_FEKi*L)
   M(12, 12) = po/3.0_FEKi
   M(10, 10) = t*L*L/105.0_FEKi + 2.0_FEKi*L*rx/15.0_FEKi
   M(11, 11) = t*L*L/105.0_FEKi + 2.0_FEKi*L*ry/15.0_FEKi
   M( 2,  4) = -11.0_FEKi*t*L/210.0_FEKi - rx/10.0_FEKi
   M( 1,  5) =  11.0_FEKi*t*L/210.0_FEKi + ry/10.0_FEKi
   M( 3,  9) = t/6.0_FEKi
   M( 5,  7) =  13._FEKi*t*L/420._FEKi - ry/10._FEKi
   M( 4,  8) = -13._FEKi*t*L/420._FEKi + rx/10._FEKi
   M( 6, 12) = po/6._FEKi
   M( 2, 10) =  13._FEKi*t*L/420._FEKi - rx/10._FEKi
   M( 1, 11) = -13._FEKi*t*L/420._FEKi + ry/10._FEKi
   M( 8, 10) =  11._FEKi*t*L/210._FEKi + rx/10._FEKi
   M( 7, 11) = -11._FEKi*t*L/210._FEKi - ry/10._FEKi
   M( 1,  7) =  9._FEKi*t/70._FEKi - 6._FEKi*ry/(5._FEKi*L)
   M( 2,  8) =  9._FEKi*t/70._FEKi - 6._FEKi*rx/(5._FEKi*L)
   M( 4, 10) = -L*L*t/140._FEKi - rx*L/30._FEKi 
   M( 5, 11) = -L*L*t/140._FEKi - ry*L/30._FEKi
   
   M( 3,  3) = M( 9,  9)
   M( 1,  1) = M( 7,  7)
   M( 2,  2) = M( 8,  8)
   M( 6,  6) = M(12, 12)
   M( 4,  4) = M(10, 10)
   M( 5,  5) = M(11, 11)
   M( 4,  2) = M( 2,  4)
   M( 5,  1) = M( 1,  5)
   M( 9,  3) = M( 3,  9)
   M( 7,  5) = M( 5,  7)
   M( 8,  4) = M( 4,  8)
   M(12,  6) = M( 6, 12)
   M(10,  2) = M( 2, 10)
   M(11,  1) = M( 1, 11)
   M(10,  8) = M( 8, 10)
   M(11,  7) = M( 7, 11)
   M( 7,  1) = M( 1,  7)
   M( 8,  2) = M( 2,  8)
   M(10,  4) = M( 4, 10)
   M(11,  5) = M( 5, 11)
   
   DC = 0.0_FEKi
   DC( 1: 3,  1: 3) = DirCos
   DC( 4: 6,  4: 6) = DirCos
   DC( 7: 9,  7: 9) = DirCos
   DC(10:12, 10:12) = DirCos
   
   M = MATMUL( MATMUL(DC, M), TRANSPOSE(DC) ) ! TODO change me if direction cosine is transposed

END SUBROUTINE ElemM_Beam
!------------------------------------------------------------------------------------------------------
!> Element stiffness matrix for pretension cable
SUBROUTINE ElemM_Cable(A, L, rho, DirCos, M)
   REAL(ReKi), INTENT( IN) :: A,rho
   REAL(FEKi), INTENT( IN) :: L
   REAL(FEKi), INTENT( IN) :: DirCos(3,3) !< From element to global: xg = DC.xe,  Kg = DC.Ke.DC^t
   REAL(FEKi), INTENT(OUT) :: M(12, 12) 
   ! Local variables
   REAL(FEKi) :: DC(12, 12)
   REAL(FEKi) :: t

   t = rho*A*L;

   M(1:12,1:12) = 0.0_FEKi

   M( 1,  1) = 13._FEKi/35._FEKi * t
   M( 2,  2) = 13._FEKi/35._FEKi * t
   M( 3,  3) = t/3.0_FEKi

   M( 7,  7) = 13._FEKi/35._FEKi * t
   M( 8,  8) = 13._FEKi/35._FEKi * t
   M( 9,  9) = t/3.0_FEKi

   M( 1,  7) =  9._FEKi/70._FEKi * t
   M( 2,  8) =  9._FEKi/70._FEKi * t
   M( 3,  9) = t/6.0_FEKi

   M( 7,  1) =  9._FEKi/70._FEKi * t 
   M( 8,  2) =  9._FEKi/70._FEKi * t
   M( 9,  3) = t/6.0_FEKi
   
   DC = 0.0_FEKi
   DC( 1: 3,  1: 3) = DirCos
   DC( 4: 6,  4: 6) = DirCos
   DC( 7: 9,  7: 9) = DirCos
   DC(10:12, 10:12) = DirCos
   
   M = MATMUL( MATMUL(DC, M), TRANSPOSE(DC) ) ! TODO: change me if DirCos convention is  transposed
END SUBROUTINE ElemM_Cable
!------------------------------------------------------------------------------------------------------
!> calculates the lumped forces and moments due to gravity on a given element:
!! the element has two nodes, with the loads for both elements stored in array F. Indexing of F is:
!!    Fx_n1=1,Fy_n1=2,Fz_n1=3,Mx_n1= 4,My_n1= 5,Mz_n1= 6,
!!    Fx_n2=7,Fy_n2=8,Fz_n2=9,Mx_n2=10,My_n2=11,Mz_n2=12
SUBROUTINE ElemG(A, L, rho, DirCos, F, g)
   REAL(ReKi), INTENT( IN ) :: A     !< area
   REAL(ReKi), INTENT( IN ) :: L     !< element length
   REAL(ReKi), INTENT( IN ) :: rho   !< density
   REAL(FEKi), INTENT( IN)  :: DirCos(3,3)      !< From element to global: xg = DC.xe,  Kg = DC.Ke.DC^t
   REAL(ReKi), INTENT( IN ) :: g     !< gravity
   REAL(FEKi), INTENT( OUT) :: F(12) !< returned loads. positions 1-6 are the loads for node 1 ; 7-12 are loads for node 2.
   REAL(FEKi) :: TempCoeff
   REAL(FEKi) :: w            ! weight per unit length
   
   F = 0.0_FEKi      ! initialize whole array to zero, then set the non-zero portions
   w = rho*A*g       ! weight per unit length
   
   ! lumped forces on both nodes (z component only):
   F(3) = -0.5_FEKi*L*w 
   F(9) = F(3)
          
   ! lumped moments on node 1 (x and y components only):
   ! bjj: note that RRD wants factor of 1/12 because of boundary conditions. Our MeshMapping routines use factor of 1/6 (assuming generic/different boundary  
   !      conditions), so we may have some inconsistent behavior. JMJ suggests using line2 elements for SubDyn's input/output meshes to improve the situation.
   TempCoeff = L*L*w/12.0_FEKi ! let's not calculate this twice  
   F(4) = -TempCoeff * DirCos(2,3) ! = -L*w*Dy/12._FEKi   !bjj: DirCos(2,3) = Dy/L
   F(5) =  TempCoeff * DirCos(1,3) ! =  L*w*Dx/12._FEKi   !bjj: DirCos(1,3) = Dx/L

      ! lumped moments on node 2: (note the opposite sign of node 1 moment)
   F(10) = -F(4)
   F(11) = -F(5)
   !F(12) is 0 for g along z alone
   
END SUBROUTINE ElemG
!------------------------------------------------------------------------------------------------------
!> 
SUBROUTINE ElemF_Cable(T0, DirCos, F)
   REAL(ReKi), INTENT( IN ) :: T0          !< Pretension load [N]
   REAL(FEKi), INTENT( IN)  :: DirCos(3,3) !< From element to global: xg = DC.xe,  Kg = DC.Ke.DC^t
   REAL(FEKi), INTENT( OUT) :: F(12)       !< returned loads. 1-6 for node 1; 7-12 for node 2.
   ! Local variables
   REAL(FEKi) :: DC(12, 12)

   F(1:12) = 0.0_FEKi  ! init 
   F(3) = +T0  
   F(9) = -T0 

   DC = 0.0_FEKi
   DC( 1: 3,  1: 3) = DirCos
   DC( 4: 6,  4: 6) = DirCos
   DC( 7: 9,  7: 9) = DirCos
   DC(10:12, 10:12) = DirCos

   F = MATMUL(DC, F)! TODO: change me if DirCos convention is  transposed

END SUBROUTINE ElemF_Cable
!------------------------------------------------------------------------------------------------------
!> Calculates the lumped gravity forces at the nodes given the element geometry
!! It assumes a linear variation of the dimensions from node 1 to node 2, thus the area may be quadratically varying if crat<>1
!! bjj: note this routine is a work in progress, intended for future version of SubDyn. Compare with ElemG.
SUBROUTINE LumpForces(Area1,Area2,crat,L,rho, g, DirCos, F)
   REAL(ReKi), INTENT( IN ) :: Area1,Area2,crat !< X-sectional areas at node 1 and node 2, t2/t1 thickness ratio
   REAL(ReKi), INTENT( IN ) :: g                !< gravity
   REAL(ReKi), INTENT( IN ) :: L                !< Length of element
   REAL(ReKi), INTENT( IN ) :: rho              !< density
   REAL(ReKi), INTENT( IN)  :: DirCos(3,3)      !< From element to global: xg = DC.xe,  Kg = DC.Ke.DC^t
   REAL(ReKi), INTENT( OUT) :: F(12)            !< Lumped forces
   !LOCALS
   REAL(ReKi)                         :: TempCoeff,a0,a1,a2  !coefficients of the gravity quadratically distributed force
   
   !Calculate quadratic polynomial coefficients
   a0 = a1
   print*,'Error: the function lumpforces is not ready to use'
   STOP

   !Calculate quadratic polynomial coefficients
   a0 = -99999 ! TODO: this is wrong
   a2 = ( (Area1+A2) - (Area1*crat+Area2/crat) )/L**2. ! *x**2
   a1 = (Area2-Area1)/L -a2*L                          ! *x
   
   !Now calculate the Lumped Forces
   F = 0
   F(3) = -(a0*L/2. +a1*L**2/6. +a2*L**3/12. )*rho*g  !Forces along z (must be negative on earth)
   F(9) = -(a0*L/2. +a1*L**2/3. +a2*L**3/4.  )*rho*g  !Forces along z (must be negative on earth)

   !Now calculate the Lumped Moments
   !HERE TO BE COMPLETED FOR THE BELOW
   TempCoeff = 1.0/12.0*g*L*L*rho*Area2  !RRD : I am changing this to >0 sign 6/10/13
      
   !F(4) = TempCoeff*( DirCos(1, 3)*DirCos(2, 1) - DirCos(1, 1)*DirCos(2, 3) ) !These do not work if convnetion on z2>z1, x2>x1, y2>y1 are not followed as I have discovered 7/23
   !F(5) = TempCoeff*( DirCos(1, 3)*DirCos(2, 2) - DirCos(1, 2)*DirCos(2, 3) ) 
   
   !RRD attempt at new dircos which keeps x in the X-Y plane
   F(4) = -TempCoeff * SQRT(1-DirCos(3,3)**2) * DirCos(1,1) !bjj: compare with ElemG() and verify this lumping is consistent
   F(5) = -TempCoeff * SQRT(1-DirCos(3,3)**2) * DirCos(2,1) !bjj: compare with ElemG() and verify this lumping is consistent
   !RRD ends
   F(10) = -F(4)
   F(11) = -F(5)
   !F(12) is 0 for g along z alone
END SUBROUTINE LumpForces

!------------------------------------------------------------------------------------------------------
!>
!! Method 1:  pinv_A = A \ eye(m) (matlab)
!!    call _GELSS to solve A.X=B 
!!    pinv(A) = B(1:n,1:m)
!! Method 2: [U,S,V] = svd(A); pinv_A = ( V / S ) * U'; (matlab) 
!     perform lapack GESVD and then  pinv(A) = V*(inv(S))*U'
SUBROUTINE PseudoInverse(A, Ainv, ErrStat, ErrMsg)
   use NWTC_LAPACK, only: LAPACK_GESVD, LAPACK_GEMM
   real(FEKi), dimension(:,:), intent(in)  :: A
   real(FEKi), dimension(:,:), allocatable :: Ainv
   INTEGER(IntKi),  INTENT(  OUT) :: ErrStat       ! < Error status of the operation
   CHARACTER(*),    INTENT(  OUT) :: ErrMsg        ! < Error message if ErrStat /    = ErrID_None
   !
   real(FEKi), dimension(:),   allocatable :: S
   real(FEKi), dimension(:,:), allocatable :: U
   real(FEKi), dimension(:,:), allocatable :: Vt
   real(FEKi), dimension(:),   allocatable :: WORK
   real(FEKi), dimension(:,:), allocatable :: Acopy
   integer :: j ! Loop indices
   integer :: M !< The number of rows of the input matrix A
   integer :: N !< The number of columns of the input matrix A
   integer :: K !< 
   integer :: L !< 
   integer :: LWORK !< 
   M = size(A,1)
   N = size(A,2)
   K = min(M,N)
   L = max(M,N)
   LWORK = MAX(1,3*K +L,5*K)
   allocate(S(K)); S = 0;
   !! LWORK >= MAX(1,3*MIN(M,N) + MAX(M,N),5*MIN(M,N)) for the other paths
   allocate(Work(LWORK)); Work=0
   allocate(U (M,K) ); U=0;
   allocate(Vt(K,N) ); Vt=0;
   allocate(Ainv(N,M)); Ainv=0;
   allocate(Acopy(M,N)); Acopy=A;

   ! --- Compute the SVD of A
   ! [U,S,V] = svd(A)
   !call DGESVD       ('S', 'S', M, N, A, M,  S, U, M  , Vt  , K,   WORK, LWORK, INFO)
   call LAPACK_GESVD('S', 'S', M, N, Acopy, S, U, Vt, WORK, LWORK, ErrStat, ErrMsg)

   !--- Compute PINV = V**T * SIGMA * U**T in two steps
   !  SIGMA = S^(-1)=1/S(j), S is diagonal
   do j = 1, K
      U(:,j) = U(:,j)/S(j)
   end do
   ! Compute Ainv = 1.0*V^t * U^t + 0.0*Ainv     V*(inv(S))*U' 
   !call DGEMM( 'T', 'T', N, M, K, 1.0, V, K, U, M, 0.0, Ainv, N)
   print*,'8'
   call LAPACK_GEMM( 'T', 'T', 1.0_FEKi, Vt, U, 0.0_FEKi, Ainv, ErrStat, ErrMsg)
   ! --- Compute rank
   !tol=maxval(shape(A))*epsilon(maxval(S))
   !rank=0
   !do i=1,K
   !   if(S(i) .gt. tol)then
   !      rank=rank+1
   !   end if
   !end do
   !print*,'Rank',rank
   !   Ainv=transpose(matmul(matmul(U(:,1:r),S_inv(1:r,1:r)),Vt(1:r,:)))
   END SUBROUTINE PseudoInverse

END MODULE FEM
