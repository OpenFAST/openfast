!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2026  National Laboratory of the Rockies
!
!    This file is part of the NWTC Subroutine Library.
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

!> Direct solver for square linear systems A*x = b that exploits a block triangular form (BTF).
!!
!! BTS_Factor uses the exact nonzero pattern of A to find a symmetric permutation P such that
!! P*A*P^T is block lower triangular. The diagonal blocks are the strongly connected components
!! of the dependency graph (row i depends on unknown j if A(i,j) /= 0 for i /= j), found with
!! Tarjan's algorithm, which returns the blocks in a valid forward-substitution order. Diagonal
!! blocks with more than one unknown are factored with LAPACK getrf (partial pivoting within the
!! block); single-unknown blocks are stored as their diagonal value; entries that couple a block
!! to earlier blocks are stored as exact nonzeros in compressed sparse row (CSR) form. BTS_Solve
!! then performs block forward substitution.
!!
!! Numerical properties:
!! - No entries are dropped and no Schur complements are formed: the blocks above the diagonal of
!!   P*A*P^T are exactly zero by construction, so pivoting across blocks is unnecessary. Block
!!   forward substitution is backward stable when each diagonal-block solve is, and results differ
!!   from a dense LU solve only by round-off.
!! - Diagonal entries do not affect the block structure. Zero diagonal entries inside a
!!   multi-unknown block are handled by pivoting within the block, and an exactly zero
!!   single-unknown block means A is singular. No maximum-transversal (zero-free diagonal)
!!   permutation is applied, so the block structure is the finest possible only when the diagonal
!!   of A has no zeros, as for the OpenFAST tight-coupling Jacobian. Matrices with zero diagonals
!!   may yield coarser blocks, which is correct but slower.
!! - A is not scaled. Scaling the rows and columns of each block by powers of two scales the
!!   computed solution exactly; callers may equilibrate A beforehand if needed.
!!
!! BTS_Factor falls back to a dense LU factorization of the full matrix (LAPACK getrf in place on
!! A) when, in automatic mode, the block form is not expected to be faster, when a single-unknown
!! pivot is numerically zero, or when a multi-unknown diagonal block is exactly singular.
!!
!! BTS_Solve does not modify the solver data or A, so solves with the same factorization may run
!! concurrently if each uses its own right-hand side and workspace. BTS_Factor modifies the solver
!! data and must not run concurrently with other calls that use the same solver data.
module BlockTriSolve

use NWTC_Base
use NWTC_Library_Types
use NWTC_LAPACK

implicit none
private

public :: BTS_Factor
public :: BTS_Solve
public :: BTS_Method_Auto, BTS_Method_Dense, BTS_Method_BlockTri

integer(IntKi), parameter :: BTS_Method_Auto     = 0 !< Select the method automatically (BTS_Factor argument only)
integer(IntKi), parameter :: BTS_Method_Dense    = 1 !< Dense LU factorization of the full matrix
integer(IntKi), parameter :: BTS_Method_BlockTri = 2 !< Block triangular form with factored diagonal blocks

integer(IntKi), parameter :: MinSizeAuto  = 64        !< Smaller systems use dense LU in automatic mode
real(R8Ki), parameter     :: MaxCostRatio = 0.75_R8Ki !< Use BTF only if its estimated solve cost is below this fraction of dense
real(R8Ki), parameter     :: SparseCost   = 2.0_R8Ki  !< Relative per-entry cost of indexed (CSR) storage vs dense storage

!> Guard against dividing by a numerically zero single-unknown pivot: minimum |pivot| relative to the
!! largest entry in its row, below which the dense LU fallback is used. This is not a stability
!! threshold; substitution through single-unknown blocks is backward stable for any nonzero pivot.
real(R8Ki), parameter     :: PivTol       = 100.0_R8Ki*epsilon(1.0_R8Ki)

!> LAPACK solve using LU factors, called directly so that the factors of each block can be read
!! from the packed storage array without copying
interface
   subroutine DGETRS(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)
      use Precision, only: R8Ki
      character, intent(in)      :: TRANS
      integer, intent(in)        :: N, NRHS, LDA, LDB
      real(R8Ki), intent(in)     :: A(LDA, *)
      integer, intent(in)        :: IPIV(*)
      real(R8Ki), intent(inout)  :: B(LDB, *)
      integer, intent(out)       :: INFO
   end subroutine
end interface

contains

!> Factor the square matrix A for subsequent calls to BTS_Solve.
!! A must be passed unchanged to BTS_Solve: if the dense method is used, A is overwritten by its LU factors.
subroutine BTS_Factor(S, A, ErrStat, ErrMsg, Method, Fallback)
   type(BlockTriSolveType), intent(inout)    :: S        !< Solver data
   real(R8Ki), intent(inout)                 :: A(:, :)  !< Matrix to factor (overwritten by LU factors if the dense method is used)
   integer(IntKi), intent(out)               :: ErrStat  !< Error status
   character(*), intent(out)                 :: ErrMsg   !< Error message
   integer(IntKi), optional, intent(in)      :: Method   !< Requested method (default BTS_Method_Auto)
   logical, optional, intent(out)            :: Fallback !< True if the block triangular form was selected but a numerically
                                                         !! zero single-unknown pivot or a singular diagonal block required dense LU

   character(*), parameter                   :: RoutineName = 'BTS_Factor'
   integer(IntKi)                            :: ErrStat2
   character(ErrMsgLen)                      :: ErrMsg2
   integer(IntKi)                            :: ReqMethod, n, nnzOff, b, i0, i1, m, k, r, e, l
   integer(B8Ki)                             :: nnzA
   integer(IntKi), allocatable               :: rp(:), ci(:), iPerm(:), BlkOf(:)
   real(R8Ki), allocatable                   :: RowMax(:)
   real(R8Ki)                                :: CostBTF, CostDense
   logical                                   :: Success

   ErrStat = ErrID_None
   ErrMsg = ''
   if (present(Fallback)) Fallback = .false.

   ReqMethod = BTS_Method_Auto
   if (present(Method)) ReqMethod = Method
   if (all(ReqMethod /= [BTS_Method_Auto, BTS_Method_Dense, BTS_Method_BlockTri])) then
      call SetErrStat(ErrID_Fatal, 'Invalid factorization method requested.', ErrStat, ErrMsg, RoutineName)
      return
   end if

   n = size(A, 1)
   if (size(A, 2) /= n) then
      call SetErrStat(ErrID_Fatal, 'Matrix must be square.', ErrStat, ErrMsg, RoutineName)
      return
   end if

   ! Clear data from any previous factorization
   call ResetData(); if (Failed()) return

   ! Dense method if requested or if the system is small
   if ((ReqMethod == BTS_Method_Dense) .or. ((ReqMethod == BTS_Method_Auto) .and. (n < MinSizeAuto))) then
      call FactorDense()
      return
   end if

   ! Count nonzeros before building the pattern so that a nearly full matrix, which cannot benefit
   ! (every nonzero costs at least one dense entry in the block form), skips the pattern storage
   call CountPattern(A, rp, RowMax, nnzA)
   CostDense = real(n, R8Ki)**2
   if ((ReqMethod == BTS_Method_Auto) .and. (real(nnzA, R8Ki) >= MaxCostRatio*CostDense)) then
      call FactorDense()
      return
   end if

   ! Off-diagonal nonzero pattern in CSR form (original ordering)
   call FillPattern(A, rp, ci)

   ! Find the diagonal blocks and their solve order
   call FindBlocks(n, rp, ci, S%Perm, S%BlkStart, S%NumBlocks)

   ! Position of each original unknown in solve order, and block of each position
   allocate (iPerm(n), BlkOf(n))
   iPerm(S%Perm) = [(k, k=1, n)]
   do b = 1, S%NumBlocks
      BlkOf(S%BlkStart(b):S%BlkStart(b + 1) - 1) = b
   end do

   ! Count entries coupling blocks; entries may only reference the same or earlier blocks
   nnzOff = 0
   do r = 1, n
      do e = rp(r), rp(r + 1) - 1
         l = iPerm(ci(e))
         if (BlkOf(l) > BlkOf(iPerm(r))) then
            call SetErrStat(ErrID_Fatal, 'Internal error: block ordering is not lower triangular.', ErrStat, ErrMsg, RoutineName)
            return
         else if (BlkOf(l) < BlkOf(iPerm(r))) then
            nnzOff = nnzOff + 1
         end if
      end do
   end do

   ! Estimated cost per solve relative to dense LU
   CostBTF = real(n, R8Ki) + SparseCost*real(nnzOff, R8Ki)
   do b = 1, S%NumBlocks
      m = S%BlkStart(b + 1) - S%BlkStart(b)
      if (m > 1) CostBTF = CostBTF + real(m, R8Ki)**2
   end do
   if ((ReqMethod == BTS_Method_Auto) .and. (CostBTF >= MaxCostRatio*CostDense)) then
      call FactorDense()
      return
   end if

   ! Factor diagonal blocks; use dense LU if a pivot is numerically zero or a block is singular
   call FactorBlocks(Success)
   if (ErrStat >= AbortErrLev) return
   if (.not. Success) then
      if (present(Fallback)) Fallback = .true.
      call FactorDense()
      return
   end if

   ! Store entries coupling each row to earlier blocks (CSR, solve ordering)
   allocate (S%RowPtr(n + 1), S%ColInd(nnzOff), S%Val(nnzOff), stat=ErrStat2)
   if (ErrStat2 /= 0) then
      call SetErrStat(ErrID_Fatal, 'Error allocating off-diagonal block storage.', ErrStat, ErrMsg, RoutineName)
      return
   end if
   S%RowPtr(1) = 1
   do k = 1, n
      r = S%Perm(k)
      S%RowPtr(k + 1) = S%RowPtr(k)
      do e = rp(r), rp(r + 1) - 1
         l = iPerm(ci(e))
         if (BlkOf(l) < BlkOf(k)) then
            S%ColInd(S%RowPtr(k + 1)) = l
            S%Val(S%RowPtr(k + 1)) = A(r, ci(e))
            S%RowPtr(k + 1) = S%RowPtr(k + 1) + 1
         end if
      end do
   end do

   S%Method = BTS_Method_BlockTri

contains

   subroutine ResetData()
      call NWTC_Library_DestroyBlockTriSolveType(S, ErrStat2, ErrMsg2)
      S%N = n
      S%Method = 0
      S%NumBlocks = 0
   end subroutine

   subroutine FactorDense()
      call ResetData(); if (Failed()) return
      allocate (S%Piv(n), stat=ErrStat2)
      if (ErrStat2 /= 0) then
         call SetErrStat(ErrID_Fatal, 'Error allocating pivot array.', ErrStat, ErrMsg, RoutineName)
         return
      end if
      if (n > 0) then
         call LAPACK_getrf(n, n, A, S%Piv, ErrStat2, ErrMsg2); if (Failed()) return
      end if
      S%Method = BTS_Method_Dense
   end subroutine

   !> Factor the diagonal blocks. Success is false (with no error set) if a single-unknown pivot is
   !! numerically zero or a multi-unknown block is exactly singular; errors are only set for
   !! allocation failures.
   subroutine FactorBlocks(Success)
      logical, intent(out)    :: Success
      real(R8Ki), allocatable :: BlkA(:, :)
      integer(IntKi)          :: LUSize, p0

      Success = .false.

      ! Packed storage for the LU factors of all multi-unknown blocks
      LUSize = 0
      do b = 1, S%NumBlocks
         m = S%BlkStart(b + 1) - S%BlkStart(b)
         if (m > 1) LUSize = LUSize + m*m
      end do
      allocate (S%BlkLUStart(S%NumBlocks), S%BlkLU(LUSize), S%Diag(n), S%Piv(n), stat=ErrStat2)
      if (ErrStat2 /= 0) then
         call SetErrStat(ErrID_Fatal, 'Error allocating diagonal block storage.', ErrStat, ErrMsg, RoutineName)
         return
      end if
      S%BlkLUStart = 0
      S%Diag = 0.0_R8Ki
      S%Piv = 0

      p0 = 1
      do b = 1, S%NumBlocks
         i0 = S%BlkStart(b)
         i1 = S%BlkStart(b + 1) - 1
         m = i1 - i0 + 1
         if (m == 1) then
            r = S%Perm(i0)
            if (abs(A(r, r)) <= PivTol*RowMax(r)) return
            S%Diag(i0) = A(r, r)
         else
            BlkA = A(S%Perm(i0:i1), S%Perm(i0:i1))
            call LAPACK_getrf(m, m, BlkA, S%Piv(i0:i1), ErrStat2, ErrMsg2)
            if (ErrStat2 /= ErrID_None) then
               ! Singular block: the dense LU fallback reports the error for the full system
               ErrStat2 = ErrID_None
               ErrMsg2 = ''
               return
            end if
            S%BlkLUStart(b) = p0
            S%BlkLU(p0:p0 + m*m - 1) = reshape(BlkA, [m*m])
            p0 = p0 + m*m
         end if
      end do

      Success = .true.
   end subroutine

   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

!> Solve A*x = b using the factorization from BTS_Factor. On entry x contains b; on exit, the solution.
!! A must be the matrix that was passed to BTS_Factor (it holds the LU factors if the dense method was used).
!! Neither S nor A is modified.
subroutine BTS_Solve(S, A, x, ErrStat, ErrMsg, Work)
   type(BlockTriSolveType), intent(in)       :: S        !< Solver data from BTS_Factor
   real(R8Ki), intent(in)                    :: A(:, :)  !< Matrix passed to BTS_Factor
   real(R8Ki), intent(inout)                 :: x(:)     !< Right-hand side on entry, solution on exit
   integer(IntKi), intent(out)               :: ErrStat  !< Error status
   character(*), intent(out)                 :: ErrMsg   !< Error message
   real(R8Ki), optional, intent(inout)       :: Work(:)  !< Workspace of at least S%N values for the block method
                                                         !! (allocated internally if absent or too small)

   character(*), parameter                   :: RoutineName = 'BTS_Solve'
   integer(IntKi)                            :: ErrStat2
   character(ErrMsgLen)                      :: ErrMsg2
   real(R8Ki), allocatable                   :: y(:)

   ErrStat = ErrID_None
   ErrMsg = ''

   if (size(x) /= S%N) then
      call SetErrStat(ErrID_Fatal, 'Size of right-hand side does not match the factored system.', ErrStat, ErrMsg, RoutineName)
      return
   end if
   if (S%N == 0) return

   select case (S%Method)

   case (BTS_Method_Dense)
      call LAPACK_getrs('N', S%N, A, S%Piv, x, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   case (BTS_Method_BlockTri)
      if (present(Work)) then
         if (size(Work) >= S%N) then
            call SolveBlocks(S, S%N, x, Work, ErrStat2, ErrMsg2)
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            return
         end if
      end if
      allocate (y(S%N), stat=ErrStat2)
      if (ErrStat2 /= 0) then
         call SetErrStat(ErrID_Fatal, 'Error allocating solve workspace.', ErrStat, ErrMsg, RoutineName)
         return
      end if
      call SolveBlocks(S, S%N, x, y, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   case default
      call SetErrStat(ErrID_Fatal, 'Matrix has not been factored.', ErrStat, ErrMsg, RoutineName)
   end select
end subroutine

!> Block forward substitution for the block triangular method, using workspace y in solve order
subroutine SolveBlocks(S, n, x, y, ErrStat, ErrMsg)
   type(BlockTriSolveType), intent(in)       :: S
   integer(IntKi), intent(in)                :: n
   real(R8Ki), intent(inout)                 :: x(:)
   real(R8Ki), intent(inout)                 :: y(n)
   integer(IntKi), intent(out)               :: ErrStat
   character(*), intent(out)                 :: ErrMsg
   integer(IntKi)                            :: b, i0, i1, k, q, m, p0
   integer                                   :: Info
   real(R8Ki)                                :: Sum

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Reorder right-hand side into solve order
   do k = 1, n
      y(k) = x(S%Perm(k))
   end do

   do b = 1, S%NumBlocks
      i0 = S%BlkStart(b)
      i1 = S%BlkStart(b + 1) - 1

      ! Subtract contributions of unknowns in earlier blocks
      do k = i0, i1
         Sum = y(k)
         do q = S%RowPtr(k), S%RowPtr(k + 1) - 1
            Sum = Sum - S%Val(q)*y(S%ColInd(q))
         end do
         y(k) = Sum
      end do

      ! Solve the diagonal block
      if (S%BlkLUStart(b) == 0) then
         y(i0) = y(i0)/S%Diag(i0)
      else
         m = i1 - i0 + 1
         p0 = S%BlkLUStart(b)
         call DGETRS('N', m, 1, S%BlkLU(p0:p0 + m*m - 1), m, S%Piv(i0:i1), y(i0:i1), m, Info)
         if (Info /= 0) then
            ErrStat = ErrID_Fatal
            write (ErrMsg, '(A,I0,A)') 'SolveBlocks: LAPACK DGETRS illegal value in argument ', -Info, '.'
            return
         end if
      end if
   end do

   ! Return solution in original ordering
   do k = 1, n
      x(S%Perm(k)) = y(k)
   end do
end subroutine

!> Count the off-diagonal nonzeros of each row (returned as CSR row pointers), the maximum absolute
!! value in each row, and the total number of nonzeros (including the diagonal).
subroutine CountPattern(A, rp, RowMax, nnzA)
   real(R8Ki), intent(in)                    :: A(:, :)
   integer(IntKi), allocatable, intent(out)  :: rp(:)
   real(R8Ki), allocatable, intent(out)      :: RowMax(:)
   integer(B8Ki), intent(out)                :: nnzA
   integer(IntKi)                            :: n, i, j

   n = size(A, 1)
   allocate (rp(n + 1), RowMax(n))
   rp = 0
   RowMax = 0.0_R8Ki
   nnzA = 0

   ! Column-major traversal; rp(i+1) temporarily holds the count for row i
   do j = 1, n
      do i = 1, n
         if (A(i, j) /= 0.0_R8Ki) then
            nnzA = nnzA + 1
            RowMax(i) = max(RowMax(i), abs(A(i, j)))
            if (i /= j) rp(i + 1) = rp(i + 1) + 1
         end if
      end do
   end do

   rp(1) = 1
   do i = 1, n
      rp(i + 1) = rp(i + 1) + rp(i)
   end do
end subroutine

!> Fill the CSR column indices of the off-diagonal nonzeros of A (ascending within each row)
subroutine FillPattern(A, rp, ci)
   real(R8Ki), intent(in)                    :: A(:, :)
   integer(IntKi), intent(in)                :: rp(:)
   integer(IntKi), allocatable, intent(out)  :: ci(:)
   integer(IntKi), allocatable               :: Next(:)
   integer(IntKi)                            :: n, i, j

   n = size(A, 1)
   allocate (ci(rp(n + 1) - 1), Next(n))
   Next = rp(1:n)
   do j = 1, n
      do i = 1, n
         if ((i /= j) .and. (A(i, j) /= 0.0_R8Ki)) then
            ci(Next(i)) = j
            Next(i) = Next(i) + 1
         end if
      end do
   end do
end subroutine

!> Find the strongly connected components of the graph with edges i -> j for each
!! off-diagonal nonzero A(i,j), using an iterative form of Tarjan's algorithm. Components
!! are emitted only after all components they depend on, which is the forward-substitution
!! order of the block lower triangular form. Unknowns within a block keep their original order.
subroutine FindBlocks(n, rp, ci, Perm, BlkStart, NumBlocks)
   integer(IntKi), intent(in)                :: n, rp(:), ci(:)
   integer(IntKi), allocatable, intent(out)  :: Perm(:), BlkStart(:)
   integer(IntKi), intent(out)               :: NumBlocks
   integer(IntKi), allocatable               :: Idx(:), Low(:), Stk(:), CallNode(:), CallEdge(:), Comp(:), Next(:)
   logical, allocatable                      :: OnStk(:)
   integer(IntKi)                            :: Counter, nStk, nCall, v0, v, w, e, i, b

   allocate (Idx(n), Low(n), Stk(n), CallNode(n), CallEdge(n), Comp(n), OnStk(n))
   Idx = 0
   Low = 0
   OnStk = .false.
   Counter = 0
   nStk = 0
   NumBlocks = 0

   do v0 = 1, n
      if (Idx(v0) /= 0) cycle

      ! Visit starting node
      nCall = 1
      CallNode(1) = v0
      CallEdge(1) = rp(v0)
      call Visit(v0)

      do while (nCall > 0)
         v = CallNode(nCall)
         e = CallEdge(nCall)
         if (e < rp(v + 1)) then
            ! Process next edge of v
            CallEdge(nCall) = e + 1
            w = ci(e)
            if (Idx(w) == 0) then
               call Visit(w)
               nCall = nCall + 1
               CallNode(nCall) = w
               CallEdge(nCall) = rp(w)
            else if (OnStk(w)) then
               Low(v) = min(Low(v), Idx(w))
            end if
         else
            ! All edges of v processed; emit component if v is its root
            if (Low(v) == Idx(v)) then
               NumBlocks = NumBlocks + 1
               do
                  w = Stk(nStk)
                  nStk = nStk - 1
                  OnStk(w) = .false.
                  Comp(w) = NumBlocks
                  if (w == v) exit
               end do
            end if
            nCall = nCall - 1
            if (nCall > 0) Low(CallNode(nCall)) = min(Low(CallNode(nCall)), Low(v))
         end if
      end do
   end do

   ! Block start positions from block sizes
   allocate (BlkStart(NumBlocks + 1), Next(NumBlocks))
   Next = 0
   do i = 1, n
      Next(Comp(i)) = Next(Comp(i)) + 1
   end do
   BlkStart(1) = 1
   do b = 1, NumBlocks
      BlkStart(b + 1) = BlkStart(b) + Next(b)
   end do

   ! Solve order: blocks in emission order, unknowns ascending within each block
   allocate (Perm(n))
   Next = BlkStart(1:NumBlocks)
   do i = 1, n
      Perm(Next(Comp(i))) = i
      Next(Comp(i)) = Next(Comp(i)) + 1
   end do

contains

   subroutine Visit(u)
      integer(IntKi), intent(in) :: u
      Counter = Counter + 1
      Idx(u) = Counter
      Low(u) = Counter
      nStk = nStk + 1
      Stk(nStk) = u
      OnStk(u) = .true.
   end subroutine
end subroutine

end module BlockTriSolve
