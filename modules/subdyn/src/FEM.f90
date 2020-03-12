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

  INTEGER,          PARAMETER  :: LAKi            = R8Ki                  ! Define the kind to be used for LAPACK routines for getting eigenvalues/vectors. Apparently there is a problem with SGGEV's eigenvectors
 
CONTAINS
!------------------------------------------------------------------------------------------------------
!> Return eigenvalues, Omega, and eigenvectors
SUBROUTINE EigenSolve(K, M, N, bCheckSingularity, EigVect, Omega2, ErrStat, ErrMsg )
   USE NWTC_LAPACK, only: LAPACK_ggev
   USE NWTC_ScaLAPACK, only : ScaLAPACK_LASRT
   INTEGER       ,          INTENT(IN   )    :: N             !< Number of degrees of freedom, size of M and K
   REAL(LaKi),              INTENT(INOUT)    :: K(N, N)       !< Stiffness matrix 
   REAL(LaKi),              INTENT(INOUT)    :: M(N, N)       !< Mass matrix 
   LOGICAL,                 INTENT(IN   )    :: bCheckSingularity !< If True, check for singularities
   REAL(LaKi),              INTENT(INOUT)    :: EigVect(N, N) !< Returned Eigenvectors
   REAL(LaKi),              INTENT(INOUT)    :: Omega2(N)      !< Returned Eigenvalues
   INTEGER(IntKi),          INTENT(  OUT)    :: ErrStat       !< Error status of the operation
   CHARACTER(*),            INTENT(  OUT)    :: ErrMsg        !< Error message if ErrStat /= ErrID_None
   ! LOCALS         
   REAL(LAKi), ALLOCATABLE                   :: WORK (:),  VL(:,:), ALPHAR(:), ALPHAI(:), BETA(:) ! eigensolver variables
   INTEGER                                   :: i  
   INTEGER                                   :: LWORK                          !variables for the eigensolver
   INTEGER,    ALLOCATABLE                   :: KEY(:)
   INTEGER(IntKi)                            :: ErrStat2
   CHARACTER(ErrMsgLen)                      :: ErrMsg2
      
   ErrStat = ErrID_None
   ErrMsg  = ''
         
   ! allocate working arrays and return arrays for the eigensolver
   LWORK=8*N + 16  !this is what the eigensolver wants  >> bjj: +16 because of MKL ?ggev documenation ( "lwork >= max(1, 8n+16) for real flavors"), though LAPACK documenation says 8n is fine
   !bjj: there seems to be a memory problem in *GGEV, so I'm making the WORK array larger to see if I can figure it out
   CALL AllocAry( Work,    LWORK, 'Work',   ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'EigenSolve') 
   CALL AllocAry( ALPHAR,  N,     'ALPHAR', ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'EigenSolve')
   CALL AllocAry( ALPHAI,  N,     'ALPHAI', ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'EigenSolve')
   CALL AllocAry( BETA,    N,     'BETA',   ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'EigenSolve')
   CALL AllocAry( VL,      N,  N, 'VL',     ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'EigenSolve')
   CALL AllocAry( KEY,     N,     'KEY',    ErrStat2, ErrMsg2 ); if(Failed()) return
    
   ! --- Eigenvalue  analysis
   ! note: SGGEV seems to have memory issues in certain cases. The eigenvalues seem to be okay, but the eigenvectors vary wildly with different compiling options.
   !       DGGEV seems to work better, so I'm making these variables LAKi (which is set to R8Ki for now)   - bjj 4/25/2014
   ! bjj: This comes from the LAPACK documentation:
   !   Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j) may easily over- or underflow, and BETA(j) may even be zero.
   !   Thus, the user should avoid naively computing the ratio alpha/beta.  However, ALPHAR and ALPHAI will be always less
   !   than and usually comparable with norm(A) in magnitude, and BETA always less than and usually comparable with norm(B).    
   ! Omega2=ALPHAR/BETA  !Note this may not be correct if ALPHAI<>0 and/or BETA=0 TO INCLUDE ERROR CHECK, also they need to be sorted
   CALL  LAPACK_ggev('N','V',N ,K, M, ALPHAR, ALPHAI, BETA, VL, EigVect, WORK, LWORK, ErrStat2, ErrMsg2)
   if(Failed()) return

   ! --- Determinign and sorting eigen frequencies 
   DO I=1,N !Initialize the key and calculate Omega
      KEY(I)=I
      IF ( AlphaR(I)<0.0_LAKi ) THEN
         !print*,'I', I, AlphaR(I), AlphaI(I), Beta(I), real(AlphaR(I)/Beta(I))
         if (bCheckSingularity) then
            ErrStat2=ErrID_Fatal
            ErrMsg2= 'Negative eigenvalue found, system may be singular (may contain rigid body modes)'
            if(Failed()) return
         endif
      !ELSE IF ( .not.EqualRealNos(AlphaI(I),0.0_LAKi )) THEN
      !   print*,'I', I, AlphaR(I), AlphaI(I), Beta(I), real(AlphaR(I)/Beta(I))
      !   ErrStat2=ErrID_Fatal
      !   ErrMsg2= 'Complex eigenvalue found, system may be singular (may contain rigid body modes)'
      !   if(Failed()) return
      ELSE IF ( EqualRealNos(real(Beta(I),ReKi),0.0_ReKi) ) THEN
         Omega2(I) = HUGE(Omega2)  ! bjj: should this be an error?
         if (bCheckSingularity) then
            ErrStat2=ErrID_Fatal
            ErrMsg2= 'Large eigenvalue found, system may be singular (may contain rigid body modes)'
            if(Failed()) return
         endif
      ELSE IF ( EqualRealNos(real(ALPHAR(I),ReKi),0.0_ReKi) ) THEN
         if (bCheckSingularity) then
            ErrStat2=ErrID_Fatal
            ErrMsg2= 'Eigenvalue zero found, system is singular (may contain rigid body modes)'
            if(Failed()) return
         endif
      ELSE
         Omega2(I) = REAL( ALPHAR(I)/BETA(I), ReKi )
      END IF           
   ENDDO  
   ! Sorting
   CALL ScaLAPACK_LASRT('I',N,Omega2,KEY,ErrStat2,ErrMsg2); if(Failed()) return 
    
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
   CALL CleanupEigen()
   RETURN

CONTAINS
   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'EigenSolve') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUpEigen()
   END FUNCTION Failed

   SUBROUTINE CleanupEigen()
      IF (ALLOCATED(Work)  ) DEALLOCATE(Work)
      IF (ALLOCATED(ALPHAR)) DEALLOCATE(ALPHAR)
      IF (ALLOCATED(ALPHAI)) DEALLOCATE(ALPHAI)
      IF (ALLOCATED(BETA)  ) DEALLOCATE(BETA)
      IF (ALLOCATED(VL)    ) DEALLOCATE(VL)
      IF (ALLOCATED(KEY)   ) DEALLOCATE(KEY)
   END SUBROUTINE CleanupEigen
  
END SUBROUTINE EigenSolve

!> Compute the determinant of a real matrix using an LU factorization
FUNCTION Determinant(A, ErrStat, ErrMsg) result(det)
   use NWTC_LAPACK, only: LAPACK_GETRF
   REAL(LaKi),      INTENT(IN   ) :: A(:, :) !< Input matrix, no side effect
   INTEGER(IntKi),  INTENT(  OUT) :: ErrStat !< Error status of the operation
   CHARACTER(*),    INTENT(  OUT) :: ErrMsg  !< Error message if ErrStat /= ErrID_None
   real(DbKi)              :: det !< May easily overflow
   integer(IntKi)          :: i
   integer                 :: n
   integer, allocatable    :: ipiv(:)
   real(LaKi), allocatable :: PLU(:,:)
   real(LaKi) :: ScaleVal

   n = size(A(1,:))
   allocate(PLU(n,n))
   allocate(ipiv(n))
   ScaleVal= 1.0_LaKi
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
   det = 1.0_LaKi
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
!> Remove degrees of freedom from a matrix (lines and rows)
SUBROUTINE RemoveDOF(A, bDOF, Ared, ErrStat, ErrMsg )
   REAL(ReKi),             INTENT(IN   ) :: A(:, :)        ! full matrix
   logical,                INTENT(IN   ) :: bDOF(:)        ! Array of logical specifying whether a DOF is to be kept(True), or removed (False)
   REAL(LAKi),ALLOCATABLE, INTENT(  OUT) :: Ared(:,:)      ! reduced matrix
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
               Ared(Ir, Jr) = REAL( A(I, J), LAKi )
            end if
         end do
      endif
   end do
END SUBROUTINE RemoveDOF

!> Expand a matrix to includes rows where bDOF is False (inverse behavior as RemoveDOF)
SUBROUTINE InsertDOFrows(Ared, bDOF, DefaultVal, A, ErrStat, ErrMsg )
   REAL(LAKi),             INTENT(IN   ) :: Ared(:, :)     ! Reduced matrix
   logical,                INTENT(IN   ) :: bDOF(:)        ! Array of logical specifying whether a DOF is to be kept(True), or removed (False)
   REAL(ReKi),             INTENT(IN   ) :: DefaultVal     ! Default value to fill the 
   REAL(ReKi)            , INTENT(INOUT) :: A(:,:)         ! Full matrix
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
         A(i,:)=real( Ared(ir,:), ReKi ) 
      else
         A(i,:)=DefaultVal
      endif   
   enddo
END SUBROUTINE InsertDOFrows
!------------------------------------------------------------------------------------------------------
!> Returns index of val in Array (val is an integer!)
! NOTE: in the future use intrinsinc function findloc
FUNCTION FINDLOCI_ReKi(Array, Val) result(i)
   real(ReKi)    , dimension(:), intent(in) :: Array !< Array to search in
   integer(IntKi), intent(in)               :: val   !< Val
   integer(IntKi)                           :: i     !< Index of joint in joint table
   logical :: found
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
   logical :: found
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
SUBROUTINE GetDirCos(P1, P2, DirCos, L, ErrStat, ErrMsg)
   REAL(ReKi) ,      INTENT(IN   )  :: P1(3), P2(3)      ! (x,y,z) global positions of two nodes making up an element
   REAL(ReKi) ,      INTENT(  OUT)  :: DirCos(3, 3)      ! calculated direction cosine matrix
   REAL(ReKi) ,      INTENT(  OUT)  :: L                 ! length of element
   INTEGER(IntKi),   INTENT(  OUT)  :: ErrStat           ! Error status of the operation
   CHARACTER(*),     INTENT(  OUT)  :: ErrMsg            ! Error message if ErrStat /= ErrID_None
   REAL(ReKi)                       :: Dx, Dy, Dz, Dxy   ! distances between nodes
   ErrMsg  = ""
   ErrStat = ErrID_None
   
   Dx=P2(1)-P1(1)
   Dy=P2(2)-P1(2)
   Dz=P2(3)-P1(3)
   Dxy = sqrt( Dx**2 + Dy**2 )
   L   = sqrt( Dx**2 + Dy**2 + Dz**2)
   
   IF ( EqualRealNos(L, 0.0_ReKi) ) THEN
      ErrMsg = ' Same starting and ending location in the element.'
      ErrStat = ErrID_Fatal
      RETURN
   ENDIF
   
   IF ( EqualRealNos(Dxy, 0.0_ReKi) ) THEN 
      DirCos=0.0_ReKi    ! whole matrix set to 0
      IF ( Dz < 0) THEN  !x is kept along global x
         DirCos(1, 1) =  1.0_ReKi
         DirCos(2, 2) = -1.0_ReKi
         DirCos(3, 3) = -1.0_ReKi
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
     
      DirCos(3, 1) = 0.0_ReKi
      DirCos(3, 2) = -Dxy/L
      DirCos(3, 3) = +Dz/L
   ENDIF

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
SUBROUTINE ElemK_Beam(A, L, Ixx, Iyy, Jzz, Shear, kappa, E, G, DirCos, K)
   REAL(ReKi), INTENT( IN) :: A, L, Ixx, Iyy, Jzz, E, G, kappa
   REAL(ReKi), INTENT( IN) :: DirCos(3,3) !< From element to global: xg = DC.xe,  Kg = DC.Ke.DC^t
   LOGICAL   , INTENT( IN) :: Shear
   REAL(ReKi), INTENT(OUT) :: K(12, 12) 
   ! Local variables
   REAL(ReKi)                            :: Ax, Ay, Kx, Ky
   REAL(ReKi)                            :: DC(12, 12)
   
   Ax = kappa*A
   Ay = kappa*A
   
   K(1:12,1:12) = 0
   
   IF (Shear) THEN
      Kx = 12.0*E*Iyy / (G*Ax*L*L)
      Ky = 12.0*E*Ixx / (G*Ay*L*L)
   ELSE
      Kx = 0.0
      Ky = 0.0
   ENDIF
      
   K( 9,  9) = E*A/L
   K( 7,  7) = 12.0*E*Iyy/( L*L*L*(1.0 + Kx) )
   K( 8,  8) = 12.0*E*Ixx/( L*L*L*(1.0 + Ky) )
   K(12, 12) = G*Jzz/L
   K(10, 10) = (4.0 + Ky)*E*Ixx / ( L*(1.0+Ky) )  
   K(11, 11) = (4.0 + Kx)*E*Iyy / ( L*(1.0+Kx) )
   K( 2,  4) = -6.*E*Ixx / ( L*L*(1.0+Ky) )
   K( 1,  5) =  6.*E*Iyy / ( L*L*(1.0+Kx) )
   K( 4, 10) = (2.0-Ky)*E*Ixx / ( L*(1.0+Ky) )
   K( 5, 11) = (2.0-Kx)*E*Iyy / ( L*(1.0+Kx) )
   
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
   
   DC = 0
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
   REAL(ReKi), INTENT( IN) :: DirCos(3,3) !< From element to global: xg = DC.xe,  Kg = DC.Ke.DC^t
   REAL(ReKi), INTENT(OUT) :: K(12, 12) 
   ! Local variables
   REAL(ReKi) :: L0, Eps0, EAL0, EE
   REAL(ReKi) :: DC(12, 12)

   Eps0 = T0/(E*A)
   L0   = L/(1+Eps0)  ! "rest length" for which pretension would be 0
   EAL0 = E*A*L0
   EE   = EAL0* Eps0/(1+Eps0)

   K(1:12,1:12)=0

   ! Note: only translatoinal DOF involved (1-3, 7-9)
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

   DC = 0
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
   REAL(ReKi), INTENT( IN) :: DirCos(3,3) !< From element to global: xg = DC.xe,  Kg = DC.Ke.DC^t
   REAL(ReKi), INTENT(OUT) :: M(12, 12)

   REAL(ReKi) :: t, rx, ry, po
   REAL(ReKi) :: DC(12, 12)
   
   t = rho*A*L;
   rx = rho*Ixx;
   ry = rho*Iyy;
   po = rho*Jzz*L;   

   M(1:12,1:12) = 0
      
   M( 9,  9) = t/3.0
   M( 7,  7) = 13.0*t/35.0 + 6.0*ry/(5.0*L)
   M( 8,  8) = 13.0*t/35.0 + 6.0*rx/(5.0*L)
   M(12, 12) = po/3.0
   M(10, 10) = t*L*L/105.0 + 2.0*L*rx/15.0
   M(11, 11) = t*L*L/105.0 + 2.0*L*ry/15.0
   M( 2,  4) = -11.0*t*L/210.0 - rx/10.0
   M( 1,  5) =  11.0*t*L/210.0 + ry/10.0
   M( 3,  9) = t/6.0
   M( 5,  7) =  13.*t*L/420. - ry/10.
   M( 4,  8) = -13.*t*L/420. + rx/10. 
   M( 6, 12) = po/6.
   M( 2, 10) =  13.*t*L/420. - rx/10. 
   M( 1, 11) = -13.*t*L/420. + ry/10.
   M( 8, 10) =  11.*t*L/210. + rx/10.
   M( 7, 11) = -11.*t*L/210. - ry/10. 
   M( 1,  7) =  9.*t/70. - 6.*ry/(5.*L)
   M( 2,  8) =  9.*t/70. - 6.*rx/(5.*L)
   M( 4, 10) = -L*L*t/140. - rx*L/30. 
   M( 5, 11) = -L*L*t/140. - ry*L/30.
   
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
   
   DC = 0
   DC( 1: 3,  1: 3) = DirCos
   DC( 4: 6,  4: 6) = DirCos
   DC( 7: 9,  7: 9) = DirCos
   DC(10:12, 10:12) = DirCos
   
   M = MATMUL( MATMUL(DC, M), TRANSPOSE(DC) ) ! TODO change me if direction cosine is transposed

END SUBROUTINE ElemM_Beam
!------------------------------------------------------------------------------------------------------
!> Element stiffness matrix for pretension cable
SUBROUTINE ElemM_Cable(A, L, rho, DirCos, M)
   REAL(ReKi), INTENT( IN) :: A, L, rho
   REAL(ReKi), INTENT( IN) :: DirCos(3,3) !< From element to global: xg = DC.xe,  Kg = DC.Ke.DC^t
   REAL(ReKi), INTENT(OUT) :: M(12, 12) 
   ! Local variables
   REAL(ReKi) :: DC(12, 12)
   REAL(ReKi) :: t

   t = rho*A*L;

   M(1:12,1:12) = 0

   M( 1,  1) = 13._ReKi/35._ReKi * t
   M( 2,  2) = 13._ReKi/35._ReKi * t
   M( 3,  3) = t/3.0_ReKi

   M( 7,  7) = 13._ReKi/35._ReKi * t
   M( 8,  8) = 13._ReKi/35._ReKi * t
   M( 9,  9) = t/3.0_ReKi

   M( 1,  7) =  9._ReKi/70._ReKi * t
   M( 2,  8) =  9._ReKi/70._ReKi * t
   M( 3,  9) = t/6.0_ReKi

   M( 7,  1) =  9._ReKi/70._ReKi * t 
   M( 8,  2) =  9._ReKi/70._ReKi * t
   M( 9,  3) = t/6.0_ReKi
   
   DC = 0
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
   REAL(ReKi), INTENT( IN)  :: DirCos(3,3)      !< From element to global: xg = DC.xe,  Kg = DC.Ke.DC^t
   REAL(ReKi), INTENT( IN ) :: g     !< gravity
   REAL(ReKi), INTENT( OUT) :: F(12) !< returned loads. positions 1-6 are the loads for node 1 ; 7-12 are loads for node 2.
   REAL(ReKi) :: TempCoeff
   REAL(ReKi) :: w            ! weight per unit length
   
   F = 0             ! initialize whole array to zero, then set the non-zero portions
   w = rho*A*g       ! weight per unit length
   
   ! lumped forces on both nodes (z component only):
   F(3) = -0.5*L*w 
   F(9) = F(3)
          
   ! lumped moments on node 1 (x and y components only):
   ! bjj: note that RRD wants factor of 1/12 because of boundary conditions. Our MeshMapping routines use factor of 1/6 (assuming generic/different boundary  
   !      conditions), so we may have some inconsistent behavior. JMJ suggests using line2 elements for SubDyn's input/output meshes to improve the situation.
   TempCoeff = L*L*w/12.0_ReKi ! let's not calculate this twice  
   F(4) = -TempCoeff * DirCos(2,3) ! = -L*w*Dy/12.   !bjj: DirCos(2,3) = Dy/L
   F(5) =  TempCoeff * DirCos(1,3) ! =  L*w*Dx/12.   !bjj: DirCos(1,3) = Dx/L

      ! lumped moments on node 2: (note the opposite sign of node 1 moment)
   F(10) = -F(4)
   F(11) = -F(5)
   !F(12) is 0 for g along z alone
   
END SUBROUTINE ElemG
!------------------------------------------------------------------------------------------------------
!> 
SUBROUTINE ElemF_Cable(T0, DirCos, F)
   REAL(ReKi), INTENT( IN ) :: T0          !< Pretension load [N]
   REAL(ReKi), INTENT( IN)  :: DirCos(3,3) !< From element to global: xg = DC.xe,  Kg = DC.Ke.DC^t
   REAL(ReKi), INTENT( OUT) :: F(12)       !< returned loads. 1-6 for node 1; 7-12 for node 2.
   ! Local variables
   REAL(ReKi) :: DC(12, 12)

   F(1:12) = 0  ! init 
   F(3) = +T0  
   F(9) = -T0 

   DC = 0
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
   use NWTC_LAPACK, only: LAPACK_DGESVD, LAPACK_GEMM
   real(LaKi), dimension(:,:), intent(in)  :: A
   real(LaKi), dimension(:,:), allocatable :: Ainv
   INTEGER(IntKi),  INTENT(  OUT) :: ErrStat       ! < Error status of the operation
   CHARACTER(*),    INTENT(  OUT) :: ErrMsg        ! < Error message if ErrStat /    = ErrID_None
   !
   real(LaKi), dimension(:),   allocatable :: S
   real(LaKi), dimension(:,:), allocatable :: U
   real(LaKi), dimension(:,:), allocatable :: Vt
   real(LaKi), dimension(:),   allocatable :: WORK
   real(LaKi), dimension(:,:), allocatable :: Acopy
   integer :: i, j ! Loop indices
   integer :: M !< The number of rows of the input matrix A
   integer :: N !< The number of columns of the input matrix A
   integer :: K !< 
   integer :: L !< 
   integer :: LWORK !< 
   integer :: INFO 
   real(ReKi) :: tol
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
   call LAPACK_DGESVD('S', 'S', M, N, Acopy, S, U, Vt, WORK, LWORK, ErrStat, ErrMsg)

   !--- Compute PINV = V**T * SIGMA * U**T in two steps
   !  SIGMA = S^(-1)=1/S(j), S is diagonal
   do j = 1, K
      U(:,j) = U(:,j)/S(j)
   end do
   ! Compute Ainv = 1.0*V^t * U^t + 0.0*Ainv     V*(inv(S))*U' 
   !call DGEMM( 'T', 'T', N, M, K, 1.0, V, K, U, M, 0.0, Ainv, N)
   call LAPACK_GEMM( 'T', 'T', 1.0_Laki, Vt, U, 0.0_LaKi, Ainv, ErrStat, ErrMsg)
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
