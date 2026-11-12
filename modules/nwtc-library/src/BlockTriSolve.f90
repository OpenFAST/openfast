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
!! Low-rank coupling: the rows of each multi-unknown block, and each run of consecutive
!! single-unknown blocks whose rows do not reference each other, form a panel whose coupling to
!! earlier unknowns is applied as a single matrix-vector product. When a panel is large, its
!! coupling matrix C is compressed with a truncated QR factorization with column pivoting,
!! C ~ U*W, if ||C - U*W||_F <= LowRankTol*||C||_F (verified after compression), the rank is at
!! most a small fraction of the panel dimensions, and the solve cost is at least halved. Otherwise
!! the exact nonzeros are kept. LowRankTol is a fixed, conservative tolerance below the round-off
!! introduced by a dense LU factorization of the same system.
!!
!! Numerical properties:
!! - No entries are dropped and no Schur complements are formed: the blocks above the diagonal of
!!   P*A*P^T are exactly zero by construction, so pivoting across blocks is unnecessary. Block
!!   forward substitution is backward stable when each diagonal-block solve is, and results differ
!!   from a dense LU solve only by round-off. Low-rank coupling perturbs the compressed panels by
!!   at most LowRankTol relative to their Frobenius norm.
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

!> Low-rank panel compression: maximum relative Frobenius-norm error of the compressed coupling
!! (below the round-off of a dense LU factorization of systems where compression applies)
real(R8Ki), parameter     :: LowRankTol      = 1.0e-13_R8Ki
integer(IntKi), parameter :: MinPanelNnz     = 4096      !< Minimum number of coupling nonzeros to attempt compression
integer(IntKi), parameter :: MinPanelDim     = 16        !< Minimum panel rows and columns to attempt compression
integer(IntKi), parameter :: MaxRankFraction = 8         !< Rank must not exceed min(rows, columns)/MaxRankFraction
real(R8Ki), parameter     :: MaxLRCostRatio  = 0.5_R8Ki  !< Low-rank solve cost must not exceed this fraction of CSR cost

!> Indices into LRInfo(:, panel)
integer(IntKi), parameter :: LR_Row0 = 1, LR_NRow = 2, LR_NCol = 3, LR_Rank = 4, LR_ColStart = 5, &
                             LR_UStart = 6, LR_WStart = 7, LR_NInfo = 7

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

   ! Compress large low-rank coupling panels
   call CompressPanels(S, ErrStat2, ErrMsg2); if (Failed()) return

   S%Method = BTS_Method_BlockTri

contains

   subroutine ResetData()
      call NWTC_Library_DestroyBlockTriSolveType(S, ErrStat2, ErrMsg2)
      S%N = n
      S%Method = 0
      S%NumBlocks = 0
      S%NumLR = 0
      S%NWork = 0
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
   real(R8Ki), optional, intent(inout)       :: Work(:)  !< Workspace of at least S%NWork values for the block method
                                                         !! (allocated internally if absent or too small)

   character(*), parameter                   :: RoutineName = 'BTS_Solve'
   integer(IntKi)                            :: ErrStat2
   character(ErrMsgLen)                      :: ErrMsg2
   real(R8Ki), allocatable                   :: w(:)

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
         if (size(Work) >= S%NWork) then
            call SolveBlocks(S, S%N, S%NWork, x, Work, ErrStat2, ErrMsg2)
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            return
         end if
      end if
      allocate (w(S%NWork), stat=ErrStat2)
      if (ErrStat2 /= 0) then
         call SetErrStat(ErrID_Fatal, 'Error allocating solve workspace.', ErrStat, ErrMsg, RoutineName)
         return
      end if
      call SolveBlocks(S, S%N, S%NWork, x, w, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   case default
      call SetErrStat(ErrID_Fatal, 'Matrix has not been factored.', ErrStat, ErrMsg, RoutineName)
   end select
end subroutine

!> Block forward substitution for the block triangular method. Workspace w holds the solution in
!! solve order (1:n), followed by the gathered panel inputs and the low-rank intermediate vector.
subroutine SolveBlocks(S, n, nw, x, w, ErrStat, ErrMsg)
   type(BlockTriSolveType), intent(in)       :: S
   integer(IntKi), intent(in)                :: n, nw
   real(R8Ki), intent(inout)                 :: x(:)
   real(R8Ki), intent(inout)                 :: w(nw)
   integer(IntKi), intent(out)               :: ErrStat
   character(*), intent(out)                 :: ErrMsg
   integer(IntKi)                            :: b, i0, i1, k, q, m, p0, ip
   integer(IntKi)                            :: r0, mR, nC, rk, cs, us, ws
   integer                                   :: Info
   real(R8Ki)                                :: Sum

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Reorder right-hand side into solve order
   do k = 1, n
      w(k) = x(S%Perm(k))
   end do

   ip = 1
   do b = 1, S%NumBlocks
      i0 = S%BlkStart(b)
      i1 = S%BlkStart(b + 1) - 1

      ! Apply a low-rank panel starting at this block: rows -= U*(W*y(cols))
      if (ip <= S%NumLR) then
         if (S%LRInfo(LR_Row0, ip) == i0) then
            r0 = S%LRInfo(LR_Row0, ip)
            mR = S%LRInfo(LR_NRow, ip)
            nC = S%LRInfo(LR_NCol, ip)
            rk = S%LRInfo(LR_Rank, ip)
            cs = S%LRInfo(LR_ColStart, ip)
            us = S%LRInfo(LR_UStart, ip)
            ws = S%LRInfo(LR_WStart, ip)
            do k = 1, nC
               w(n + k) = w(S%LRCols(cs + k - 1))
            end do
            call DGEMV('N', rk, nC, 1.0_R8Ki, S%LRW(ws:ws + rk*nC - 1), rk, w(n + 1:n + nC), 1, &
                             0.0_R8Ki, w(n + nC + 1:n + nC + rk), 1)
            call DGEMV('N', mR, rk, -1.0_R8Ki, S%LRU(us:us + mR*rk - 1), mR, w(n + nC + 1:n + nC + rk), 1, &
                             1.0_R8Ki, w(r0:r0 + mR - 1), 1)
            ip = ip + 1
         end if
      end if

      ! Subtract contributions of unknowns in earlier blocks
      do k = i0, i1
         Sum = w(k)
         do q = S%RowPtr(k), S%RowPtr(k + 1) - 1
            Sum = Sum - S%Val(q)*w(S%ColInd(q))
         end do
         w(k) = Sum
      end do

      ! Solve the diagonal block
      if (S%BlkLUStart(b) == 0) then
         w(i0) = w(i0)/S%Diag(i0)
      else
         m = i1 - i0 + 1
         p0 = S%BlkLUStart(b)
         call DGETRS('N', m, 1, S%BlkLU(p0:p0 + m*m - 1), m, S%Piv(i0:i1), w(i0:i1), m, Info)
         if (Info /= 0) then
            ErrStat = ErrID_Fatal
            write (ErrMsg, '(A,I0,A)') 'SolveBlocks: LAPACK DGETRS illegal value in argument ', -Info, '.'
            return
         end if
      end if
   end do

   ! Return solution in original ordering
   do k = 1, n
      x(S%Perm(k)) = w(k)
   end do
end subroutine

!> Find coupling panels and store large low-rank ones as U*W, removing their rows from the CSR
!! storage. A panel is a multi-unknown block, or a run of consecutive single-unknown blocks whose
!! rows only reference unknowns before the run. Also sets the solve workspace size S%NWork.
subroutine CompressPanels(S, ErrStat, ErrMsg)
   type(BlockTriSolveType), intent(inout)    :: S
   integer(IntKi), intent(out)               :: ErrStat
   character(*), intent(out)                 :: ErrMsg

   character(*), parameter                   :: RoutineName = 'CompressPanels'
   integer(IntKi)                            :: n, b, i0, i1, k, q, NumLR, MaxNC, MaxRank, nnzNew
   integer(IntKi), allocatable               :: ColMap(:), Cols(:), Info(:, :), LRCols(:), NewRowPtr(:), NewColInd(:)
   real(R8Ki), allocatable                   :: LRU(:), LRW(:), NewVal(:)
   logical, allocatable                      :: InPanel(:)

   ErrStat = ErrID_None
   ErrMsg = ''

   n = S%N
   S%NumLR = 0
   S%NWork = n
   if (S%NumBlocks == 0) return

   allocate (ColMap(n), InPanel(n), Info(LR_NInfo, S%NumBlocks), LRCols(0), LRU(0), LRW(0))
   ColMap = 0
   InPanel = .false.
   NumLR = 0
   MaxNC = 0
   MaxRank = 0

   b = 1
   do while (b <= S%NumBlocks)
      i0 = S%BlkStart(b)
      i1 = S%BlkStart(b + 1) - 1
      if (i1 == i0) then
         ! Extend the run while the next block is a single unknown that does not reference the run
         do while (b < S%NumBlocks)
            k = S%BlkStart(b + 1)
            if (S%BlkStart(b + 2) - k /= 1) exit
            if (any(S%ColInd(S%RowPtr(k):S%RowPtr(k + 1) - 1) >= i0)) exit
            b = b + 1
            i1 = k
         end do
      end if
      call TryPanel(i0, i1)
      if (ErrStat >= AbortErrLev) return
      b = b + 1
   end do

   if (NumLR == 0) return

   S%NumLR = NumLR
   S%LRInfo = Info(:, 1:NumLR)
   call move_alloc(LRCols, S%LRCols)
   call move_alloc(LRU, S%LRU)
   call move_alloc(LRW, S%LRW)
   S%NWork = n + MaxNC + MaxRank

   ! Remove the entries of compressed rows from the CSR storage
   nnzNew = 0
   do k = 1, n
      if (.not. InPanel(k)) nnzNew = nnzNew + S%RowPtr(k + 1) - S%RowPtr(k)
   end do
   allocate (NewRowPtr(n + 1), NewColInd(nnzNew), NewVal(nnzNew), stat=k)
   if (k /= 0) then
      call SetErrStat(ErrID_Fatal, 'Error allocating coupling storage.', ErrStat, ErrMsg, RoutineName)
      return
   end if
   NewRowPtr(1) = 1
   do k = 1, n
      NewRowPtr(k + 1) = NewRowPtr(k)
      if (InPanel(k)) cycle
      do q = S%RowPtr(k), S%RowPtr(k + 1) - 1
         NewColInd(NewRowPtr(k + 1)) = S%ColInd(q)
         NewVal(NewRowPtr(k + 1)) = S%Val(q)
         NewRowPtr(k + 1) = NewRowPtr(k + 1) + 1
      end do
   end do
   call move_alloc(NewRowPtr, S%RowPtr)
   call move_alloc(NewColInd, S%ColInd)
   call move_alloc(NewVal, S%Val)

contains

   !> Compress the coupling of rows r0:r1 (solve order) if it is large and low rank
   subroutine TryPanel(r0, r1)
      integer(IntKi), intent(in) :: r0, r1
      real(R8Ki), allocatable    :: C(:, :), C0(:, :), U(:, :), W(:, :)
      integer(IntKi)             :: mR, nC, nnzP, RMax, Rank, j, stat
      real(R8Ki)                 :: NormF, ErrF
      logical                    :: Success

      mR = r1 - r0 + 1
      nnzP = S%RowPtr(r1 + 1) - S%RowPtr(r0)
      if ((mR < MinPanelDim) .or. (nnzP < MinPanelNnz)) return

      ! Columns referenced by the panel
      allocate (Cols(nnzP))
      nC = 0
      do q = S%RowPtr(r0), S%RowPtr(r1 + 1) - 1
         j = S%ColInd(q)
         if (ColMap(j) == 0) then
            nC = nC + 1
            ColMap(j) = nC
            Cols(nC) = j
         end if
      end do

      ! Largest useful rank: well below the panel dimensions and at least halving the solve cost
      RMax = min(min(mR, nC)/MaxRankFraction, &
                 int((MaxLRCostRatio*SparseCost*real(nnzP, R8Ki) - real(nC, R8Ki))/real(mR + nC, R8Ki), IntKi))

      if ((nC >= MinPanelDim) .and. (RMax >= 1)) then

         ! Dense coupling matrix
         allocate (C(mR, nC), stat=stat)
         if (stat /= 0) then
            call SetErrStat(ErrID_Fatal, 'Error allocating coupling panel.', ErrStat, ErrMsg, RoutineName)
            return
         end if
         C = 0.0_R8Ki
         do k = r0, r1
            do q = S%RowPtr(k), S%RowPtr(k + 1) - 1
               C(k - r0 + 1, ColMap(S%ColInd(q))) = S%Val(q)
            end do
         end do
         C0 = C
         NormF = norm2(C0)

         call TruncatedQRCP(C, LowRankTol, RMax, Rank, U, W, Success)

         ! Verify the compressed panel against the exact coupling
         if (Success) Success = Rank >= 1
         if (Success) then
            ErrF = norm2(C0 - matmul(U, W))
            Success = ErrF <= LowRankTol*NormF
         end if

         if (Success) then
            NumLR = NumLR + 1
            Info(LR_Row0, NumLR) = r0
            Info(LR_NRow, NumLR) = mR
            Info(LR_NCol, NumLR) = nC
            Info(LR_Rank, NumLR) = Rank
            Info(LR_ColStart, NumLR) = size(LRCols) + 1
            Info(LR_UStart, NumLR) = size(LRU) + 1
            Info(LR_WStart, NumLR) = size(LRW) + 1
            LRCols = [LRCols, Cols(1:nC)]
            LRU = [LRU, reshape(U, [mR*Rank])]
            LRW = [LRW, reshape(W, [Rank*nC])]
            InPanel(r0:r1) = .true.
            MaxNC = max(MaxNC, nC)
            MaxRank = max(MaxRank, Rank)
         end if
      end if

      ColMap(Cols(1:nC)) = 0
      deallocate (Cols)
   end subroutine
end subroutine

!> Truncated QR factorization with column pivoting of C (m x nc), stopped at the smallest rank r
!! with ||C*P - Q(:,1:r)*R(1:r,:)||_F <= Tol*||C||_F, where the residual norm is computed exactly
!! from the trailing block at each step. On success, C ~ U*W with U = Q(:,1:r) (m x r, orthonormal
!! columns) and W = R(1:r,:)*P^T (r x nc). Success is false if the tolerance is not met with rank <= RMax.
!! C is overwritten.
subroutine TruncatedQRCP(C, Tol, RMax, Rank, U, W, Success)
   real(R8Ki), intent(inout)                 :: C(:, :)
   real(R8Ki), intent(in)                    :: Tol
   integer(IntKi), intent(in)                :: RMax
   integer(IntKi), intent(out)               :: Rank
   real(R8Ki), allocatable, intent(out)      :: U(:, :), W(:, :)
   logical, intent(out)                      :: Success

   real(R8Ki), allocatable                   :: ColNorm2(:), Tau(:), Tmp(:)
   integer(IntKi), allocatable               :: Jpvt(:)
   integer(IntKi)                            :: m, nc, k, j, ic, KMax, iTmp
   real(R8Ki)                                :: Tol2, Alpha, XNorm, Beta, s, dTmp

   m = size(C, 1)
   nc = size(C, 2)
   KMax = min(RMax, m, nc)
   Rank = 0
   Success = .false.

   allocate (ColNorm2(nc), Tau(KMax), Jpvt(nc))
   do ic = 1, nc
      ColNorm2(ic) = dot_product(C(:, ic), C(:, ic))
      Jpvt(ic) = ic
   end do
   Tol2 = Tol**2*sum(ColNorm2)

   do k = 1, KMax + 1

      ! Residual of the rank k-1 approximation is the trailing block C(k:m, k:nc)
      if (sum(ColNorm2(k:nc)) <= Tol2) then
         Success = .true.
         exit
      end if
      if (k > KMax) exit

      ! Move the column with the largest remaining norm to position k
      j = k - 1 + maxloc(ColNorm2(k:nc), 1)
      if (j /= k) then
         Tmp = C(:, k)
         C(:, k) = C(:, j)
         C(:, j) = Tmp
         dTmp = ColNorm2(k); ColNorm2(k) = ColNorm2(j); ColNorm2(j) = dTmp
         iTmp = Jpvt(k); Jpvt(k) = Jpvt(j); Jpvt(j) = iTmp
      end if

      ! Householder reflector H = I - Tau*v*v^T with v = [1; C(k+1:m, k)] annihilating C(k+1:m, k)
      Alpha = C(k, k)
      XNorm = norm2(C(k + 1:m, k))
      if (XNorm == 0.0_R8Ki) then
         Tau(k) = 0.0_R8Ki
      else
         Beta = -sign(hypot(Alpha, XNorm), Alpha)
         Tau(k) = (Beta - Alpha)/Beta
         C(k + 1:m, k) = C(k + 1:m, k)/(Alpha - Beta)
         C(k, k) = Beta
      end if

      ! Apply the reflector to the trailing columns and recompute their remaining norms
      do ic = k + 1, nc
         if (Tau(k) /= 0.0_R8Ki) then
            s = Tau(k)*(C(k, ic) + dot_product(C(k + 1:m, k), C(k + 1:m, ic)))
            C(k, ic) = C(k, ic) - s
            C(k + 1:m, ic) = C(k + 1:m, ic) - s*C(k + 1:m, k)
         end if
         ColNorm2(ic) = dot_product(C(k + 1:m, ic), C(k + 1:m, ic))
      end do
      Rank = k
   end do

   if (.not. Success) return

   ! U = H(1)*...*H(Rank)*[I; 0], accumulated backward
   allocate (U(m, Rank))
   U = 0.0_R8Ki
   do k = 1, Rank
      U(k, k) = 1.0_R8Ki
   end do
   do k = Rank, 1, -1
      if (Tau(k) == 0.0_R8Ki) cycle
      do ic = k, Rank
         s = Tau(k)*(U(k, ic) + dot_product(C(k + 1:m, k), U(k + 1:m, ic)))
         U(k, ic) = U(k, ic) - s
         U(k + 1:m, ic) = U(k + 1:m, ic) - s*C(k + 1:m, k)
      end do
   end do

   ! W = R(1:Rank, :) with columns returned to their original order
   allocate (W(Rank, nc))
   do ic = 1, nc
      j = min(ic, Rank)
      W(1:j, Jpvt(ic)) = C(1:j, ic)
      W(j + 1:Rank, Jpvt(ic)) = 0.0_R8Ki
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
