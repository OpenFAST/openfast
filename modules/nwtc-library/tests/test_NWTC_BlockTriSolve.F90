module test_NWTC_BlockTriSolve

use testdrive, only: new_unittest, unittest_type, error_type, check
use NWTC_Base
use NWTC_Library_Types
use NWTC_LAPACK
use ModReg
use ModMesh_Types, only: MESH_NEWCOPY
use BlockTriSolve

implicit none
private
public :: test_NWTC_BlockTriSolve_suite

real(R8Ki), parameter :: RelTol = 1.0e-11_R8Ki

contains

subroutine test_NWTC_BlockTriSolve_suite(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
               new_unittest("dense path matches LAPACK", test_dense_path), &
               new_unittest("block structure detected", test_block_structure), &
               new_unittest("block solve matches dense solve", test_block_solve), &
               new_unittest("automatic method selection", test_auto_selection), &
               new_unittest("diagonal and identity rows", test_diagonal), &
               new_unittest("small pivot falls back to dense", test_small_pivot_fallback), &
               new_unittest("singular block falls back to dense", test_singular_block_fallback), &
               new_unittest("singular matrix returns error", test_singular_error), &
               new_unittest("refactorization switches methods", test_refactor), &
               new_unittest("copy and pack round trip", test_copy_pack), &
               new_unittest("empty and scalar systems", test_trivial_sizes), &
               new_unittest("invalid arguments", test_invalid_args) &
               ]
end subroutine

!-------------------------------------------------------------------------------
! Tests
!-------------------------------------------------------------------------------

subroutine test_dense_path(error)
   type(error_type), allocatable, intent(out) :: error
   type(BlockTriSolveType)                    :: S
   real(R8Ki), allocatable                    :: A(:, :), b(:), x(:), xRef(:)
   integer(IntKi)                             :: ErrStat, seed
   character(ErrMsgLen)                       :: ErrMsg

   seed = 11
   call RandomMatrix(seed, 20, A)
   call RandomVector(seed, 20, b)
   call DenseReference(A, b, xRef)

   x = b
   call BTS_Factor(S, A, ErrStat, ErrMsg)
   call check(error, ErrStat, ErrID_None, message=trim(ErrMsg)); if (allocated(error)) return
   call check(error, S%Method, BTS_Method_Dense, message="small system should use dense method"); if (allocated(error)) return
   call BTS_Solve(S, A, x, ErrStat, ErrMsg)
   call check(error, ErrStat, ErrID_None, message=trim(ErrMsg)); if (allocated(error)) return

   ! Same LAPACK calls as the reference, so results must be identical
   call check(error, all(x == xRef), message="dense path differs from LAPACK reference")
end subroutine

subroutine test_block_structure(error)
   type(error_type), allocatable, intent(out) :: error
   type(BlockTriSolveType)                    :: S
   real(R8Ki), allocatable                    :: A(:, :)
   integer(IntKi), allocatable                :: Sizes(:), Owner(:), BlkOfPos(:), iPerm(:), GotSizes(:)
   integer(IntKi), parameter                  :: BlockSizes(10) = [1, 1, 3, 1, 5, 2, 1, 1, 4, 1]
   integer(IntKi)                             :: ErrStat, seed, n, b, i, j, k
   character(ErrMsgLen)                       :: ErrMsg
   logical                                    :: Valid

   seed = 21
   call BlockTriMatrix(seed, BlockSizes, 0.4_R8Ki, A, Owner)
   n = size(A, 1)

   call BTS_Factor(S, A, ErrStat, ErrMsg, Method=BTS_Method_BlockTri)
   call check(error, ErrStat, ErrID_None, message=trim(ErrMsg)); if (allocated(error)) return
   call check(error, S%Method, BTS_Method_BlockTri); if (allocated(error)) return
   call check(error, S%NumBlocks, size(BlockSizes), message="wrong number of blocks"); if (allocated(error)) return

   ! Block sizes match (as multisets)
   GotSizes = S%BlkStart(2:) - S%BlkStart(:S%NumBlocks)
   Sizes = BlockSizes
   call SortInt(GotSizes)
   call SortInt(Sizes)
   call check(error, all(GotSizes == Sizes), message="block sizes do not match"); if (allocated(error)) return

   ! Each detected block contains unknowns from exactly one constructed block
   Valid = .true.
   do b = 1, S%NumBlocks
      k = Owner(S%Perm(S%BlkStart(b)))
      Valid = Valid .and. all(Owner(S%Perm(S%BlkStart(b):S%BlkStart(b + 1) - 1)) == k)
   end do
   call check(error, Valid, message="detected block mixes constructed blocks"); if (allocated(error)) return

   ! Permutation is valid and the permuted matrix is block lower triangular
   allocate (iPerm(n), BlkOfPos(n))
   iPerm = 0
   iPerm(S%Perm) = [(k, k=1, n)]
   call check(error, all(iPerm > 0), message="permutation is not a bijection"); if (allocated(error)) return
   do b = 1, S%NumBlocks
      BlkOfPos(S%BlkStart(b):S%BlkStart(b + 1) - 1) = b
   end do
   Valid = .true.
   do j = 1, n
      do i = 1, n
         if (A(i, j) /= 0.0_R8Ki) Valid = Valid .and. (BlkOfPos(iPerm(j)) <= BlkOfPos(iPerm(i)))
      end do
   end do
   call check(error, Valid, message="permuted matrix is not block lower triangular")
end subroutine

subroutine test_block_solve(error)
   type(error_type), allocatable, intent(out) :: error
   type(BlockTriSolveType)                    :: S
   real(R8Ki), allocatable                    :: A(:, :), A0(:, :), b(:), x(:), xRef(:)
   integer(IntKi), allocatable                :: Owner(:)
   integer(IntKi), parameter                  :: BlockSizes(12) = [2, 1, 1, 7, 1, 3, 1, 1, 12, 1, 4, 1]
   integer(IntKi)                             :: ErrStat, seed, iRhs, n
   character(ErrMsgLen)                       :: ErrMsg

   seed = 31
   call BlockTriMatrix(seed, BlockSizes, 0.3_R8Ki, A, Owner)
   A0 = A
   n = size(A, 1)

   call BTS_Factor(S, A, ErrStat, ErrMsg, Method=BTS_Method_BlockTri)
   call check(error, ErrStat, ErrID_None, message=trim(ErrMsg)); if (allocated(error)) return
   call check(error, S%Method, BTS_Method_BlockTri); if (allocated(error)) return
   call check(error, all(A == A0), message="block method must not modify the matrix"); if (allocated(error)) return

   ! Several right-hand sides with the same factorization
   do iRhs = 1, 4
      call RandomVector(seed, n, b)
      call DenseReference(A0, b, xRef)
      x = b
      call BTS_Solve(S, A, x, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None, message=trim(ErrMsg)); if (allocated(error)) return
      call check(error, maxval(abs(x - xRef)) <= RelTol*maxval(abs(xRef)), &
                 message="block solution differs from dense solution"); if (allocated(error)) return
      call check(error, ResidualNorm(A0, x, b) <= RelTol, message="block solution residual too large"); if (allocated(error)) return
   end do
end subroutine

subroutine test_auto_selection(error)
   type(error_type), allocatable, intent(out) :: error
   type(BlockTriSolveType)                    :: S
   real(R8Ki), allocatable                    :: A(:, :), A0(:, :), b(:), x(:), xRef(:)
   integer(IntKi), allocatable                :: Owner(:), BlockSizes(:)
   integer(IntKi)                             :: ErrStat, seed, n
   character(ErrMsgLen)                       :: ErrMsg

   ! Mostly single-unknown blocks with a coupled core, similar to the OpenFAST solver Jacobian
   seed = 41
   BlockSizes = [spread(1_IntKi, 1, 60), 40_IntKi, spread(1_IntKi, 1, 80)]
   call BlockTriMatrix(seed, BlockSizes, 0.05_R8Ki, A, Owner)
   A0 = A
   n = size(A, 1)
   call BTS_Factor(S, A, ErrStat, ErrMsg)
   call check(error, ErrStat, ErrID_None, message=trim(ErrMsg)); if (allocated(error)) return
   call check(error, S%Method, BTS_Method_BlockTri, message="sparse block structure should use block method"); if (allocated(error)) return
   call RandomVector(seed, n, b)
   call DenseReference(A0, b, xRef)
   x = b
   call BTS_Solve(S, A, x, ErrStat, ErrMsg)
   call check(error, ErrStat, ErrID_None, message=trim(ErrMsg)); if (allocated(error)) return
   call check(error, maxval(abs(x - xRef)) <= RelTol*maxval(abs(xRef)), message="auto block solution incorrect"); if (allocated(error)) return
   call check(error, ResidualNorm(A0, x, b) <= RelTol, message="auto block solution residual too large"); if (allocated(error)) return

   ! Fully coupled dense matrix
   call RandomMatrix(seed, 100, A)
   call BTS_Factor(S, A, ErrStat, ErrMsg)
   call check(error, ErrStat, ErrID_None, message=trim(ErrMsg)); if (allocated(error)) return
   call check(error, S%Method, BTS_Method_Dense, message="dense matrix should use dense method"); if (allocated(error)) return

   ! Sparse but strongly connected matrix (single block): no benefit from the block form
   call SparseCoupledMatrix(seed, 120, A)
   call BTS_Factor(S, A, ErrStat, ErrMsg)
   call check(error, ErrStat, ErrID_None, message=trim(ErrMsg)); if (allocated(error)) return
   call check(error, S%Method, BTS_Method_Dense, message="single coupled block should use dense method")
end subroutine

subroutine test_diagonal(error)
   type(error_type), allocatable, intent(out) :: error
   type(BlockTriSolveType)                    :: S
   real(R8Ki), allocatable                    :: A(:, :), b(:), x(:)
   real(R8Ki)                                 :: d(7)
   integer(IntKi)                             :: ErrStat, i, seed
   character(ErrMsgLen)                       :: ErrMsg

   d = [1.0_R8Ki, -2.5_R8Ki, 3.0e6_R8Ki, 1.0_R8Ki, 7.0e-4_R8Ki, 1.0_R8Ki, -9.0_R8Ki]
   allocate (A(7, 7))
   A = 0.0_R8Ki
   do i = 1, 7
      A(i, i) = d(i)
   end do
   seed = 51
   call RandomVector(seed, 7, b)
   x = b

   call BTS_Factor(S, A, ErrStat, ErrMsg, Method=BTS_Method_BlockTri)
   call check(error, ErrStat, ErrID_None, message=trim(ErrMsg)); if (allocated(error)) return
   call check(error, S%Method, BTS_Method_BlockTri); if (allocated(error)) return
   call check(error, S%NumBlocks, 7, message="diagonal matrix should have one block per unknown"); if (allocated(error)) return
   call check(error, size(S%Blk), 0, message="diagonal matrix should have no multi-unknown blocks"); if (allocated(error)) return
   call BTS_Solve(S, A, x, ErrStat, ErrMsg)
   call check(error, ErrStat, ErrID_None, message=trim(ErrMsg)); if (allocated(error)) return
   call check(error, all(x == b/d), message="diagonal solve should be exact division")
end subroutine

subroutine test_small_pivot_fallback(error)
   type(error_type), allocatable, intent(out) :: error
   type(BlockTriSolveType)                    :: S
   real(R8Ki), allocatable                    :: A(:, :), A0(:, :), b(:), x(:), xRef(:)
   integer(IntKi), allocatable                :: Owner(:)
   integer(IntKi)                             :: ErrStat, seed, n, i
   character(ErrMsgLen)                       :: ErrMsg

   seed = 61
   call BlockTriMatrix(seed, [1, 4, 1, 1, 3, 1], 0.5_R8Ki, A, Owner)
   n = size(A, 1)

   ! Make one single-unknown pivot numerically zero relative to its row
   do i = 1, n
      if (count(Owner == Owner(i)) == 1 .and. count(A(i, :) /= 0.0_R8Ki) > 1) exit
   end do
   call check(error, i <= n, message="test matrix has no coupled single-unknown row"); if (allocated(error)) return
   A(i, i) = 1.0e-20_R8Ki*maxval(abs(A(i, :)))
   A0 = A

   call RandomVector(seed, n, b)
   call DenseReference(A0, b, xRef)
   call BTS_Factor(S, A, ErrStat, ErrMsg, Method=BTS_Method_BlockTri)
   call check(error, ErrStat, ErrID_None, message=trim(ErrMsg)); if (allocated(error)) return
   call check(error, S%Method, BTS_Method_Dense, message="tiny pivot should fall back to dense"); if (allocated(error)) return
   x = b
   call BTS_Solve(S, A, x, ErrStat, ErrMsg)
   call check(error, ErrStat, ErrID_None, message=trim(ErrMsg)); if (allocated(error)) return
   call check(error, all(x == xRef), message="fallback solution differs from LAPACK reference")
end subroutine

subroutine test_singular_block_fallback(error)
   type(error_type), allocatable, intent(out) :: error
   type(BlockTriSolveType)                    :: S
   real(R8Ki), allocatable                    :: A(:, :)
   integer(IntKi)                             :: ErrStat
   character(ErrMsgLen)                       :: ErrMsg

   ! 2x2 coupled block that is exactly singular, preceded by a single unknown
   allocate (A(3, 3))
   A = reshape([1.0_R8Ki, 0.5_R8Ki, 0.5_R8Ki, &
                0.0_R8Ki, 1.0_R8Ki, 2.0_R8Ki, &
                0.0_R8Ki, 2.0_R8Ki, 4.0_R8Ki], [3, 3])
   call BTS_Factor(S, A, ErrStat, ErrMsg, Method=BTS_Method_BlockTri)

   ! Falls back to dense LU, which reports the singular matrix like LAPACK does
   call check(error, ErrStat, ErrID_Fatal, message="singular matrix should return a fatal error"); if (allocated(error)) return
   call check(error, S%Method /= BTS_Method_BlockTri, message="singular block must not use block method")
end subroutine

subroutine test_singular_error(error)
   type(error_type), allocatable, intent(out) :: error
   type(BlockTriSolveType)                    :: S
   real(R8Ki), allocatable                    :: A(:, :), x(:)
   integer(IntKi)                             :: ErrStat, seed
   character(ErrMsgLen)                       :: ErrMsg

   seed = 71
   call RandomMatrix(seed, 70, A)
   A(13, :) = 0.0_R8Ki
   call BTS_Factor(S, A, ErrStat, ErrMsg)
   call check(error, ErrStat, ErrID_Fatal, message="zero row should return a fatal error"); if (allocated(error)) return

   ! Solve must not be attempted on an unfactored system
   allocate (x(70))
   x = 1.0_R8Ki
   call BTS_Solve(S, A, x, ErrStat, ErrMsg)
   call check(error, ErrStat, ErrID_Fatal, message="solve after failed factorization should fail")
end subroutine

subroutine test_refactor(error)
   type(error_type), allocatable, intent(out) :: error
   type(BlockTriSolveType)                    :: S
   real(R8Ki), allocatable                    :: A1(:, :), A2(:, :), A0(:, :), b(:), x(:), xRef(:)
   integer(IntKi), allocatable                :: Owner(:)
   integer(IntKi)                             :: ErrStat, seed, iPass
   character(ErrMsgLen)                       :: ErrMsg

   seed = 81
   do iPass = 1, 2

      ! Block triangular system
      call BlockTriMatrix(seed, [3, 1, 1, 6, 1, 2], 0.3_R8Ki, A1, Owner)
      A0 = A1
      call BTS_Factor(S, A1, ErrStat, ErrMsg, Method=BTS_Method_BlockTri)
      call check(error, S%Method, BTS_Method_BlockTri); if (allocated(error)) return
      call RandomVector(seed, size(A0, 1), b)
      call DenseReference(A0, b, xRef)
      x = b
      call BTS_Solve(S, A1, x, ErrStat, ErrMsg)
      call check(error, maxval(abs(x - xRef)) <= RelTol*maxval(abs(xRef)), message="block refactor solution incorrect"); if (allocated(error)) return

      ! Different size dense system with the same solver data
      call RandomMatrix(seed, 9, A2)
      A0 = A2
      call BTS_Factor(S, A2, ErrStat, ErrMsg)
      call check(error, S%Method, BTS_Method_Dense); if (allocated(error)) return
      call check(error, S%N, 9); if (allocated(error)) return
      call check(error, .not. allocated(S%Blk), message="block data should be cleared"); if (allocated(error)) return
      call RandomVector(seed, 9, b)
      call DenseReference(A0, b, xRef)
      x = b
      call BTS_Solve(S, A2, x, ErrStat, ErrMsg)
      call check(error, all(x == xRef), message="dense refactor solution incorrect"); if (allocated(error)) return
   end do
end subroutine

subroutine test_copy_pack(error)
   type(error_type), allocatable, intent(out) :: error
   type(BlockTriSolveType)                    :: S, S2, S3
   type(RegFile)                              :: RF
   real(R8Ki), allocatable                    :: A(:, :), b(:), x1(:), x2(:), x3(:)
   integer(IntKi), allocatable                :: Owner(:)
   integer(IntKi)                             :: ErrStat, seed, Un
   character(ErrMsgLen)                       :: ErrMsg
   character(*), parameter                    :: FileName = 'test_NWTC_BlockTriSolve_pack.bin'

   seed = 91
   call BlockTriMatrix(seed, [1, 5, 1, 1, 3, 1, 2], 0.4_R8Ki, A, Owner)
   call BTS_Factor(S, A, ErrStat, ErrMsg, Method=BTS_Method_BlockTri)
   call check(error, S%Method, BTS_Method_BlockTri); if (allocated(error)) return
   call RandomVector(seed, size(A, 1), b)
   x1 = b
   call BTS_Solve(S, A, x1, ErrStat, ErrMsg)

   ! Copy (used when duplicating turbine data)
   call NWTC_Library_CopyBlockTriSolveType(S, S2, MESH_NEWCOPY, ErrStat, ErrMsg)
   call check(error, ErrStat, ErrID_None, message=trim(ErrMsg)); if (allocated(error)) return
   x2 = b
   call BTS_Solve(S2, A, x2, ErrStat, ErrMsg)
   call check(error, all(x2 == x1), message="copied solver gives different result"); if (allocated(error)) return

   ! Pack and unpack (used for checkpoint/restart)
   open (newunit=Un, file=FileName, access='stream', form='unformatted', status='replace', action='write')
   call InitRegFile(RF, Un, ErrStat, ErrMsg)
   call check(error, ErrStat, ErrID_None, message=trim(ErrMsg)); if (allocated(error)) return
   call NWTC_Library_PackBlockTriSolveType(RF, S)
   call CloseRegFile(RF, ErrStat, ErrMsg)
   call check(error, ErrStat, ErrID_None, message=trim(ErrMsg)); if (allocated(error)) return

   open (newunit=Un, file=FileName, access='stream', form='unformatted', status='old', action='read')
   call OpenRegFile(RF, Un, ErrStat, ErrMsg)
   call check(error, ErrStat, ErrID_None, message=trim(ErrMsg)); if (allocated(error)) return
   call NWTC_Library_UnPackBlockTriSolveType(RF, S3)
   close (Un, status='delete')
   call check(error, RF%ErrStat, ErrID_None, message=trim(RF%ErrMsg)); if (allocated(error)) return

   x3 = b
   call BTS_Solve(S3, A, x3, ErrStat, ErrMsg)
   call check(error, ErrStat, ErrID_None, message=trim(ErrMsg)); if (allocated(error)) return
   call check(error, all(x3 == x1), message="unpacked solver gives different result")
end subroutine

subroutine test_trivial_sizes(error)
   type(error_type), allocatable, intent(out) :: error
   type(BlockTriSolveType)                    :: S
   real(R8Ki), allocatable                    :: A(:, :), x(:)
   integer(IntKi)                             :: ErrStat, iMethod
   character(ErrMsgLen)                       :: ErrMsg

   do iMethod = BTS_Method_Auto, BTS_Method_BlockTri
      allocate (A(0, 0), x(0))
      call BTS_Factor(S, A, ErrStat, ErrMsg, Method=iMethod)
      call check(error, ErrStat, ErrID_None, message=trim(ErrMsg)); if (allocated(error)) return
      call BTS_Solve(S, A, x, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None, message=trim(ErrMsg)); if (allocated(error)) return
      deallocate (A, x)

      allocate (A(1, 1), x(1))
      A = 4.0_R8Ki
      x = 2.0_R8Ki
      call BTS_Factor(S, A, ErrStat, ErrMsg, Method=iMethod)
      call check(error, ErrStat, ErrID_None, message=trim(ErrMsg)); if (allocated(error)) return
      call BTS_Solve(S, A, x, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None, message=trim(ErrMsg)); if (allocated(error)) return
      call check(error, x(1), 0.5_R8Ki, message="scalar solve incorrect"); if (allocated(error)) return
      deallocate (A, x)
   end do
end subroutine

subroutine test_invalid_args(error)
   type(error_type), allocatable, intent(out) :: error
   type(BlockTriSolveType)                    :: S
   real(R8Ki), allocatable                    :: A(:, :), x(:)
   integer(IntKi)                             :: ErrStat, seed
   character(ErrMsgLen)                       :: ErrMsg

   allocate (A(3, 4))
   A = 1.0_R8Ki
   call BTS_Factor(S, A, ErrStat, ErrMsg)
   call check(error, ErrStat, ErrID_Fatal, message="non-square matrix should fail"); if (allocated(error)) return

   seed = 101
   call RandomMatrix(seed, 5, A)
   call BTS_Factor(S, A, ErrStat, ErrMsg, Method=99)
   call check(error, ErrStat, ErrID_Fatal, message="invalid method should fail"); if (allocated(error)) return

   call BTS_Factor(S, A, ErrStat, ErrMsg)
   allocate (x(4))
   x = 1.0_R8Ki
   call BTS_Solve(S, A, x, ErrStat, ErrMsg)
   call check(error, ErrStat, ErrID_Fatal, message="wrong right-hand side size should fail")
end subroutine

!-------------------------------------------------------------------------------
! Helpers
!-------------------------------------------------------------------------------

!> Deterministic uniform random number in [-1, 1) (portable across compilers)
function RandUniform(seed) result(r)
   integer(IntKi), intent(inout) :: seed
   real(R8Ki)                    :: r
   integer(B8Ki)                 :: s
   s = modulo(1103515245_B8Ki*int(seed, B8Ki) + 12345_B8Ki, 2147483648_B8Ki)
   seed = int(s, IntKi)
   r = 2.0_R8Ki*real(s, R8Ki)/2147483648.0_R8Ki - 1.0_R8Ki
end function

subroutine RandomVector(seed, n, v)
   integer(IntKi), intent(inout)          :: seed
   integer(IntKi), intent(in)             :: n
   real(R8Ki), allocatable, intent(out)   :: v(:)
   integer(IntKi)                         :: i
   allocate (v(n))
   do i = 1, n
      v(i) = RandUniform(seed)
   end do
end subroutine

!> Dense, well-conditioned random matrix
subroutine RandomMatrix(seed, n, A)
   integer(IntKi), intent(inout)          :: seed
   integer(IntKi), intent(in)             :: n
   real(R8Ki), allocatable, intent(out)   :: A(:, :)
   integer(IntKi)                         :: i, j
   allocate (A(n, n))
   do j = 1, n
      do i = 1, n
         A(i, j) = RandUniform(seed)
      end do
      A(j, j) = A(j, j) + real(n, R8Ki)
   end do
end subroutine

!> Sparse matrix whose dependency graph is a single cycle plus random entries (one strongly connected block)
subroutine SparseCoupledMatrix(seed, n, A)
   integer(IntKi), intent(inout)          :: seed
   integer(IntKi), intent(in)             :: n
   real(R8Ki), allocatable, intent(out)   :: A(:, :)
   integer(IntKi)                         :: i, j
   allocate (A(n, n))
   A = 0.0_R8Ki
   do i = 1, n
      A(i, i) = 4.0_R8Ki + RandUniform(seed)
      A(i, modulo(i, n) + 1) = RandUniform(seed)
      j = 1 + int(0.5_R8Ki*(RandUniform(seed) + 1.0_R8Ki)*real(n - 1, R8Ki))
      A(i, j) = A(i, j) + 0.5_R8Ki*RandUniform(seed)
   end do
end subroutine

!> Randomly permuted block lower triangular matrix. Each multi-unknown diagonal block is strongly
!! connected (contains a cycle) and has small diagonals in places so that LU pivoting is exercised.
!! Owner(i) is the constructed block containing unknown i of the permuted matrix.
subroutine BlockTriMatrix(seed, BlockSizes, Density, A, Owner)
   integer(IntKi), intent(inout)             :: seed
   integer(IntKi), intent(in)                :: BlockSizes(:)
   real(R8Ki), intent(in)                    :: Density
   real(R8Ki), allocatable, intent(out)      :: A(:, :)
   integer(IntKi), allocatable, intent(out)  :: Owner(:)
   real(R8Ki), allocatable                   :: T(:, :)
   integer(IntKi), allocatable               :: Start(:), Blk(:), p(:)
   integer(IntKi)                            :: n, nb, b, i, j, m, k, tmp

   nb = size(BlockSizes)
   n = sum(BlockSizes)
   allocate (T(n, n), Start(nb + 1), Blk(n), p(n), Owner(n))
   T = 0.0_R8Ki
   Start(1) = 1
   do b = 1, nb
      Start(b + 1) = Start(b) + BlockSizes(b)
      Blk(Start(b):Start(b + 1) - 1) = b
   end do

   do b = 1, nb
      m = BlockSizes(b)
      if (m == 1) then
         T(Start(b), Start(b)) = sign(1.0_R8Ki + abs(RandUniform(seed)), RandUniform(seed))
      else
         do j = Start(b), Start(b + 1) - 1
            do i = Start(b), Start(b + 1) - 1
               T(i, j) = RandUniform(seed)
            end do
            T(j, j) = T(j, j) + real(m, R8Ki)
         end do
         ! First diagonal entry small to force a row interchange
         T(Start(b), Start(b)) = 1.0e-3_R8Ki*T(Start(b), Start(b))
      end if
   end do

   ! Coupling to earlier blocks, with a wide range of magnitudes
   do i = 1, n
      do j = 1, Start(Blk(i)) - 1
         if (0.5_R8Ki*(RandUniform(seed) + 1.0_R8Ki) < Density) T(i, j) = RandUniform(seed)*10.0_R8Ki**int(3.0_R8Ki*RandUniform(seed))
      end do
   end do

   ! Make diagonals dominate the coupling to earlier blocks so long chains stay well conditioned
   ! (except the first row of multi-unknown blocks, which keeps a small pivot)
   do i = 1, n
      if ((BlockSizes(Blk(i)) == 1) .or. (i /= Start(Blk(i)))) then
         T(i, i) = T(i, i) + sign(sum(abs(T(i, 1:Start(Blk(i)) - 1))), T(i, i))
      end if
   end do

   ! Random symmetric permutation (Fisher-Yates)
   p = [(i, i=1, n)]
   do i = n, 2, -1
      k = 1 + int(0.5_R8Ki*(RandUniform(seed) + 1.0_R8Ki)*real(i, R8Ki))
      k = min(k, i)
      tmp = p(i)
      p(i) = p(k)
      p(k) = tmp
   end do
   A = T(p, p)
   Owner = Blk(p)
end subroutine

!> Reference solution with dense LAPACK LU (same calls as the dense method)
subroutine DenseReference(A, b, x)
   real(R8Ki), intent(in)                 :: A(:, :), b(:)
   real(R8Ki), allocatable, intent(out)   :: x(:)
   real(R8Ki), allocatable                :: LU(:, :)
   integer, allocatable                   :: Piv(:)
   integer(IntKi)                         :: ErrStat
   character(ErrMsgLen)                   :: ErrMsg
   LU = A
   x = b
   allocate (Piv(size(A, 1)))
   call LAPACK_getrf(size(A, 1), size(A, 1), LU, Piv, ErrStat, ErrMsg)
   call LAPACK_getrs('N', size(A, 1), LU, Piv, x, ErrStat, ErrMsg)
end subroutine

!> Relative residual ||A*x - b||_inf / (||A||_inf*||x||_inf + ||b||_inf)
function ResidualNorm(A, x, b) result(r)
   real(R8Ki), intent(in) :: A(:, :), x(:), b(:)
   real(R8Ki)             :: r
   r = maxval(abs(matmul(A, x) - b))/(maxval(sum(abs(A), dim=2))*maxval(abs(x)) + maxval(abs(b)))
end function

subroutine SortInt(v)
   integer(IntKi), intent(inout) :: v(:)
   integer(IntKi)                :: i, j, t
   do i = 2, size(v)
      t = v(i)
      j = i - 1
      do while (j >= 1)
         if (v(j) <= t) exit
         v(j + 1) = v(j)
         j = j - 1
      end do
      v(j + 1) = t
   end do
end subroutine

end module
