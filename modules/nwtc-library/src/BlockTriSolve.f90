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
!! of the dependency graph (row i depends on unknown j if A(i,j) /= 0), found with Tarjan's
!! algorithm, which returns the blocks in a valid forward-substitution order. Diagonal blocks
!! with more than one unknown are factored with LAPACK getrf; single-unknown blocks are stored
!! as their diagonal value; entries that couple a block to earlier blocks are stored as exact
!! nonzeros in compressed sparse row (CSR) form. BTS_Solve then performs block forward
!! substitution. No entries are dropped, so results differ from a dense LU solve only by
!! round-off.
!!
!! BTS_Factor falls back to a dense LU factorization of the full matrix (LAPACK getrf in place
!! on A, as used previously by callers) when the block form is not expected to be faster, when
!! a single-unknown pivot is numerically zero, or when a diagonal block is singular.
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

integer(IntKi), parameter :: MinSizeAuto  = 64                        !< Smaller systems use dense LU in automatic mode
real(R8Ki), parameter     :: MaxCostRatio = 0.75_R8Ki                 !< Use BTF only if its estimated solve cost is below this fraction of dense
real(R8Ki), parameter     :: SparseCost   = 2.0_R8Ki                  !< Relative per-entry cost of indexed (CSR) storage vs dense storage
real(R8Ki), parameter     :: PivTol       = 100.0_R8Ki*epsilon(1.0_R8Ki) !< Minimum |single-unknown pivot| relative to the largest entry in its row

contains

!> Factor the square matrix A for subsequent calls to BTS_Solve.
!! A must be passed unchanged to BTS_Solve: if the dense method is used, A is overwritten by its LU factors.
subroutine BTS_Factor(S, A, ErrStat, ErrMsg, Method)
   type(BlockTriSolveType), intent(inout)    :: S        !< Solver data
   real(R8Ki), intent(inout)                 :: A(:, :)  !< Matrix to factor (overwritten by LU factors if the dense method is used)
   integer(IntKi), intent(out)               :: ErrStat  !< Error status
   character(*), intent(out)                 :: ErrMsg   !< Error message
   integer(IntKi), optional, intent(in)      :: Method   !< Requested method (default BTS_Method_Auto)

   character(*), parameter                   :: RoutineName = 'BTS_Factor'
   integer(IntKi)                            :: ErrStat2
   character(ErrMsgLen)                      :: ErrMsg2
   integer(IntKi)                            :: ReqMethod, n, nnzA, nnzOff, b, i0, i1, m, k, r, e, l
   integer(IntKi), allocatable               :: rp(:), ci(:), iPerm(:), BlkOf(:)
   real(R8Ki), allocatable                   :: RowMax(:)
   real(R8Ki)                                :: CostBTF, CostDense
   logical                                   :: Success

   ErrStat = ErrID_None
   ErrMsg = ''

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

   ! Nonzero pattern of the off-diagonal entries in CSR form (original ordering)
   call BuildPattern(A, rp, ci, RowMax, nnzA)

   ! Every nonzero costs at least one dense entry in the block form, so a matrix
   ! that is nearly full cannot benefit
   CostDense = real(n, R8Ki)**2
   if ((ReqMethod == BTS_Method_Auto) .and. (real(nnzA, R8Ki) >= MaxCostRatio*CostDense)) then
      call FactorDense()
      return
   end if

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

   ! Factor diagonal blocks; fall back to dense LU if any block is unsuitable
   call FactorBlocks(Success); if (Failed()) return
   if (.not. Success) then
      call FactorDense()
      return
   end if

   ! Store entries coupling each row to earlier blocks (CSR, solve ordering)
   allocate (S%RowPtr(n + 1), S%ColInd(nnzOff), S%Val(nnzOff), S%Work(n), stat=ErrStat2)
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
   S%Work = 0.0_R8Ki

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

   subroutine FactorBlocks(Success)
      logical, intent(out) :: Success
      integer(IntKi)       :: NumMulti, iBlk, j

      Success = .false.

      NumMulti = count((S%BlkStart(2:) - S%BlkStart(:S%NumBlocks)) > 1)
      allocate (S%BlkIdx(S%NumBlocks), S%Blk(NumMulti), S%Diag(n), stat=ErrStat2)
      if (ErrStat2 /= 0) then
         call SetErrStat(ErrID_Fatal, 'Error allocating diagonal block storage.', ErrStat, ErrMsg, RoutineName)
         return
      end if
      S%BlkIdx = 0
      S%Diag = 0.0_R8Ki

      iBlk = 0
      do b = 1, S%NumBlocks
         i0 = S%BlkStart(b)
         i1 = S%BlkStart(b + 1) - 1
         m = i1 - i0 + 1
         if (m == 1) then
            r = S%Perm(i0)
            if (abs(A(r, r)) <= PivTol*RowMax(r)) return
            S%Diag(i0) = A(r, r)
         else
            iBlk = iBlk + 1
            S%BlkIdx(b) = iBlk
            allocate (S%Blk(iBlk)%LU(m, m), S%Blk(iBlk)%Piv(m), stat=ErrStat2)
            if (ErrStat2 /= 0) then
               call SetErrStat(ErrID_Fatal, 'Error allocating diagonal block.', ErrStat, ErrMsg, RoutineName)
               return
            end if
            do j = 1, m
               S%Blk(iBlk)%LU(:, j) = A(S%Perm(i0:i1), S%Perm(i0 + j - 1))
            end do
            call LAPACK_getrf(m, m, S%Blk(iBlk)%LU, S%Blk(iBlk)%Piv, ErrStat2, ErrMsg2)
            if (ErrStat2 /= ErrID_None) return
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
subroutine BTS_Solve(S, A, x, ErrStat, ErrMsg)
   type(BlockTriSolveType), intent(inout)    :: S        !< Solver data (workspace is updated)
   real(R8Ki), intent(in)                    :: A(:, :)  !< Matrix passed to BTS_Factor
   real(R8Ki), intent(inout)                 :: x(:)     !< Right-hand side on entry, solution on exit
   integer(IntKi), intent(out)               :: ErrStat  !< Error status
   character(*), intent(out)                 :: ErrMsg   !< Error message

   character(*), parameter                   :: RoutineName = 'BTS_Solve'
   integer(IntKi)                            :: ErrStat2
   character(ErrMsgLen)                      :: ErrMsg2
   integer(IntKi)                            :: b, i0, i1, k, q
   real(R8Ki)                                :: Sum

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
      associate (y => S%Work, Perm => S%Perm, RowPtr => S%RowPtr, ColInd => S%ColInd, Val => S%Val)

         ! Reorder right-hand side into solve order
         do k = 1, S%N
            y(k) = x(Perm(k))
         end do

         ! Block forward substitution
         do b = 1, S%NumBlocks
            i0 = S%BlkStart(b)
            i1 = S%BlkStart(b + 1) - 1

            ! Subtract contributions of unknowns in earlier blocks
            do k = i0, i1
               Sum = y(k)
               do q = RowPtr(k), RowPtr(k + 1) - 1
                  Sum = Sum - Val(q)*y(ColInd(q))
               end do
               y(k) = Sum
            end do

            ! Solve the diagonal block
            if (S%BlkIdx(b) == 0) then
               y(i0) = y(i0)/S%Diag(i0)
            else
               associate (Blk => S%Blk(S%BlkIdx(b)))
                  call LAPACK_getrs('N', i1 - i0 + 1, Blk%LU, Blk%Piv, y(i0:i1), ErrStat2, ErrMsg2)
               end associate
               if (ErrStat2 /= ErrID_None) then
                  call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                  if (ErrStat >= AbortErrLev) return
               end if
            end if
         end do

         ! Return solution in original ordering
         do k = 1, S%N
            x(Perm(k)) = y(k)
         end do
      end associate

   case default
      call SetErrStat(ErrID_Fatal, 'Matrix has not been factored.', ErrStat, ErrMsg, RoutineName)
   end select
end subroutine

!> Build the CSR nonzero pattern of the off-diagonal entries of A, the maximum absolute
!! value in each row, and the total number of nonzeros (including the diagonal).
subroutine BuildPattern(A, rp, ci, RowMax, nnzA)
   real(R8Ki), intent(in)                    :: A(:, :)
   integer(IntKi), allocatable, intent(out)  :: rp(:), ci(:)
   real(R8Ki), allocatable, intent(out)      :: RowMax(:)
   integer(IntKi), intent(out)               :: nnzA
   integer(IntKi), allocatable               :: Next(:)
   integer(IntKi)                            :: n, i, j

   n = size(A, 1)
   allocate (rp(n + 1), RowMax(n), Next(n))
   Next = 0
   RowMax = 0.0_R8Ki
   nnzA = 0

   ! Count nonzeros per row (column-major traversal)
   do j = 1, n
      do i = 1, n
         if (A(i, j) /= 0.0_R8Ki) then
            nnzA = nnzA + 1
            RowMax(i) = max(RowMax(i), abs(A(i, j)))
            if (i /= j) Next(i) = Next(i) + 1
         end if
      end do
   end do

   rp(1) = 1
   do i = 1, n
      rp(i + 1) = rp(i) + Next(i)
   end do

   ! Fill column indices; ascending within each row
   allocate (ci(rp(n + 1) - 1))
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
