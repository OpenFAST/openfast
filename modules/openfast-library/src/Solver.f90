module Solver

use FAST_Solver
use FAST_ModTypes
use FAST_Eval
use NWTC_LAPACK
use ElastoDyn
use BeamDyn
use SubDyn
use AeroDyn
use AeroDyn14
use ServoDyn
use SC_DataEx

implicit none

private

! q matrix column indices (displacement, velocity, accel, algo accel)
integer(IntKi), parameter  :: COL_D = 1, COL_V = 2, COL_A = 3, COL_AA = 4

! Tight coupling modules
integer(IntKi), parameter  :: TC_Modules(*) = [Module_ED, Module_BD, Module_SD]

public Solver_Init, Solver_Step0, Solver_Step

! Transfer Module Flags
integer(IntKi), parameter  :: XM_x = 1, &
                              XM_dxdt = 2, &
                              XM_dYdx = 4, &
                              XM_dXdx = 8, &
                              XM_dYdu = 16, &
                              XM_dXdu = 32

integer(IntKi) :: munit

contains

subroutine Solver_Init(p, m, Mods, ErrStat, ErrMsg)

   type(TC_ParameterType), intent(inout)     :: p           !< Parameters
   type(TC_MiscVarType), intent(out)         :: m           !< Misc variables for optimization (not copied in glue code)
   type(ModDataType), intent(inout)          :: Mods(:)  !< Solution variables from modules
   integer(IntKi), intent(out)               :: ErrStat     !< Error status of the operation
   character(*), intent(out)                 :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   character(*), parameter                   :: RoutineName = 'Solver_Init'
   integer(IntKi)                            :: ErrStat2    ! local error status
   character(ErrMsgLen)                      :: ErrMsg2     ! local error message
   integer(IntKi)                            :: i, j, k, n
   integer(IntKi)                            :: NumX, NumU, NumY, NumQ, NumJac
   integer(IntKi), allocatable               :: iq(:), modIDs(:), vec1(:), vec2(:)
   logical                                   :: isLoad
   logical, allocatable                      :: isLoadTight(:), isLoadOption1(:)
   type(TC_MappingType)                      :: MeshMap

   !----------------------------------------------------------------------------
   ! Module ordering for solve
   !----------------------------------------------------------------------------

   ! Get array of module IDs
   modIDs = [(Mods(i)%ID, i=1, size(Mods))]

   ! Ordering array for all modules
   p%iModAll = [pack([(i, i=1, size(Mods))], ModIDs == Module_SrvD), &
                pack([(i, i=1, size(Mods))], ModIDs == Module_AD), &
                pack([(i, i=1, size(Mods))], ModIDs == Module_ED), &
                pack([(i, i=1, size(Mods))], ModIDs == Module_BD), &
                pack([(i, i=1, size(Mods))], ModIDs == Module_SD)]

   ! Ordering array for tight coupling modules
   p%iModTC = [pack([(i, i=1, size(Mods))], ModIDs == Module_ED), &
               pack([(i, i=1, size(Mods))], ModIDs == Module_BD), &
               pack([(i, i=1, size(Mods))], ModIDs == Module_SD)]

   ! Ordering array for Option 2 solve
   p%iModOpt2 = [pack([(i, i=1, size(Mods))], ModIDs == Module_ED), &
                 pack([(i, i=1, size(Mods))], ModIDs == Module_BD), &
                 pack([(i, i=1, size(Mods))], ModIDs == Module_SD)]

   ! Ordering array for Option 1 solve
   p%iModOpt1 = [pack([(i, i=1, size(Mods))], ModIDs == Module_ED), &
                 pack([(i, i=1, size(Mods))], ModIDs == Module_BD), &
                 pack([(i, i=1, size(Mods))], ModIDs == Module_SD), &
                 pack([(i, i=1, size(Mods))], ModIDs == Module_ExtPtfm), &
                 pack([(i, i=1, size(Mods))], ModIDs == Module_HD), &
                 pack([(i, i=1, size(Mods))], ModIDs == Module_Orca)]

   ! Ordering array for Option 1 modules that were not in Option 2
   ! These modules need to do update states and calc output before Option 1 solve
   p%iModOpt1US = [pack([(i, i=1, size(Mods))], ModIDs == Module_ExtPtfm), &
                   pack([(i, i=1, size(Mods))], ModIDs == Module_HD), &
                   pack([(i, i=1, size(Mods))], ModIDs == Module_Orca)]

   !----------------------------------------------------------------------------
   ! Initialize mesh mappings (must be done before calculating global indices)
   !----------------------------------------------------------------------------

   call Solver_DefineMappings(m%Mappings, Mods, ErrStat2, ErrMsg2); if (Failed()) return

   !----------------------------------------------------------------------------
   ! Calculate Variable cateogries and global indices
   ! TODO: reorder to improve data locality
   !----------------------------------------------------------------------------

   NumX = 0
   do i = 1, size(Mods)
      vec1 = [(0, k=1, 0)]
      vec2 = [(0, k=1, 0)]
      do j = 1, size(Mods(i)%Vars%x)
         Mods(i)%Vars%x(j)%Cat = VarCategory(Mods(i)%ID, Mods(i)%Vars%x(j)%Field)
         if (Mods(i)%Vars%x(j)%Cat == VC_Tight) then
            Mods(i)%Vars%x(j)%iGblSol = [(NumX + k, k=1, Mods(i)%Vars%x(j)%Size)]
            vec1 = [vec1, Mods(i)%Vars%x(j)%iLoc]
            vec2 = [vec2, Mods(i)%Vars%x(j)%iGblSol]
            NumX = NumX + Mods(i)%Vars%x(j)%Size
         end if
      end do
      Mods(i)%ixs = reshape([vec1, vec2], [size(vec1), 2])
   end do

   NumU = 0
   do i = 1, size(Mods)
      vec1 = [(0, k=1, 0)]
      vec2 = [(0, k=1, 0)]
      do j = 1, size(Mods(i)%Vars%u)
         Mods(i)%Vars%u(j)%Cat = VarCategory(Mods(i)%ID, Mods(i)%Vars%u(j)%Field)
         if (iand(Mods(i)%Vars%u(j)%Flags, VF_Solve) > 0) then
            Mods(i)%Vars%u(j)%iGblSol = [(NumU + k, k=1, Mods(i)%Vars%u(j)%Size)]
            vec1 = [vec1, Mods(i)%Vars%u(j)%iLoc]
            vec2 = [vec2, Mods(i)%Vars%u(j)%iGblSol]
            NumU = NumU + Mods(i)%Vars%u(j)%Size
         end if
      end do
      Mods(i)%ius = reshape([vec1, vec2], [size(vec1), 2])
   end do

   NumY = 0
   do i = 1, size(Mods)
      vec1 = [(0, k=1, 0)]
      vec2 = [(0, k=1, 0)]
      do j = 1, size(Mods(i)%Vars%y)
         Mods(i)%Vars%y(j)%Cat = VarCategory(Mods(i)%ID, Mods(i)%Vars%y(j)%Field)
         if (iand(Mods(i)%Vars%y(j)%Flags, VF_Solve) > 0) then
            Mods(i)%Vars%y(j)%iGblSol = [(NumY + k, k=1, Mods(i)%Vars%y(j)%Size)]
            vec1 = [vec1, Mods(i)%Vars%y(j)%iLoc]
            vec2 = [vec2, Mods(i)%Vars%y(j)%iGblSol]
            NumY = NumY + Mods(i)%Vars%y(j)%Size
         end if
      end do
      Mods(i)%iys = reshape([vec1, vec2], [size(vec1), 2])
   end do

   !----------------------------------------------------------------------------
   ! Allocate storage
   !----------------------------------------------------------------------------

   ! State and State Derivative
   call AllocAry(m%x, NumX, "m%x", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%xn, NumX, "m%xn", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%dx, NumX, "m%dx", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%dxdt, NumX, "m%dxdt", ErrStat2, ErrMsg2); if (Failed()) return
   m%x = 0.0_R8Ki
   m%dx = 0.0_R8Ki
   m%dxdt = 0.0_R8Ki

   ! Inputs and Input Temporary
   call AllocAry(m%u, NumU, "m%u", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%un, NumU, "m%un", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%du, NumU, "m%du", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%u_tmp, NumU, "m%u_tmp", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%UDiff, NumU, "m%UDiff", ErrStat2, ErrMsg2); if (Failed()) return

   ! Outputs
   call AllocAry(m%y, NumY, "m%y", ErrStat2, ErrMsg2); if (Failed()) return
   m%y = 0.0_R8Ki

   !----------------------------------------------------------------------------
   ! Allocate q storage for generalized alpha algorithm
   ! This matrix stores equation state in an (N,4) array where:
   !  - N is the number of equations (rows)
   !  - Column 1 is position
   !  - Column 2 is velocity
   !  - Column 3 is acceleration
   !  - Column 4 is generalized alpha algorithmic acceleration
   !----------------------------------------------------------------------------

   ! Allocate storage for matrix mapping mapping x to q
   call AllocAry(p%ixqd, 3, NumX, "p%ixqd", ErrStat2, ErrMsg2); if (Failed()) return
   p%ixqd = 0

   ! Initialize number of q states (ignore derivatives)
   NumQ = 0
   n = 0

   ! Loop through modules
   do i = 1, size(Mods)

      ! Allocate iq to store q index for each variable
      call AllocAry(iq, size(Mods(i)%Vars%x), "iq", ErrStat2, ErrMsg2); if (Failed()) return
      iq = 0

      ! Skip modules that aren't in tight coupling
      if (all(Mods(i)%ID /= TC_Modules)) cycle

      ! Loop through state variables
      do j = 1, size(Mods(i)%Vars%x)

         ! Skip variables which already have a q index
         if (iq(j) /= 0) cycle

         ! Skip load variables (force and moment)
         if (any(Mods(i)%Vars%x(j)%Field == LoadFields)) cycle

         ! Set q index for variable and update number
         iq(j) = NumQ + 1
         NumQ = NumQ + Mods(i)%Vars%x(j)%Size

         ! Loop through remaining vars  if the names match
         do k = i + 1, size(Mods(i)%Vars%x)

            ! If names are different then they don't match, skip
            if (Mods(i)%Vars%x(j)%Name /= Mods(i)%Vars%x(k)%Name) cycle

            ! If field is not the same or a derivative of current field, skip
            select case (Mods(i)%Vars%x(j)%Field)
            case (VF_TransDisp, VF_TransVel, VF_TransAcc)
               if (all(Mods(i)%Vars%x(k)%Field /= TransFields)) cycle
            case (VF_Orientation, VF_AngularDisp, VF_AngularVel, VF_AngularAcc)
               if (all(Mods(i)%Vars%x(k)%Field /= AngularFields)) cycle
            case (VF_Force, VF_Moment)
               cycle
            end select

            ! Copy q index
            iq(k) = iq(j)

         end do
      end do

      ! Loop through state variables and build mapping array between x and q
      ! ixqd is 3xN where each row is [global x array index, q matrix row, q matrix col]
      do j = 1, size(Mods(i)%Vars%x)
         do k = 1, Mods(i)%Vars%x(j)%Size
            n = n + 1
            p%ixqd(:, n) = [Mods(i)%Vars%x(j)%iGblSol(k), iq(j) + k - 1, Mods(i)%Vars%x(j)%DerivOrder + 1]
         end do
      end do

      ! Deallocate iq for use by next module
      deallocate (iq)
   end do

   ! Remove x->q indicies that aren't set
   p%ixqd = p%ixqd(:, 1:n)

   ! Allocate/initialize global q storage
   call AllocAry(m%q, NumQ, 4, "m%q", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%qn, NumQ, 4, "m%qn", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%dq, NumQ, 4, "m%qn", ErrStat2, ErrMsg2); if (Failed()) return
   m%q = 0.0_R8Ki
   m%qn = 0.0_R8Ki
   m%dq = 0.0_R8Ki

   !----------------------------------------------------------------------------
   ! Mapping arrays
   !----------------------------------------------------------------------------

   ! Calculate index mapping arrays for X1Tight and X2Tight
   call AllocAry(p%iX1Tight, 0, "p%iX1Tight", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(p%iX2Tight, 0, "p%iX2Tight", ErrStat2, ErrMsg2); if (Failed()) return
   do i = 1, size(Mods)
      do j = 1, size(Mods(i)%Vars%x)
         if (Mods(i)%Vars%x(j)%Cat == VC_Tight) then
            select case (Mods(i)%Vars%x(j)%DerivOrder)
            case (0)
               p%iX1Tight = [p%iX1Tight, Mods(i)%Vars%x(j)%iGblSol]
            case (1)
               p%iX2Tight = [p%iX2Tight, Mods(i)%Vars%x(j)%iGblSol]
            end select
         end if
      end do
   end do

   ! Calculate index mapping arrays for U^Tight and U^Option1
   call AllocAry(p%iuLoad, 0, "p%iuLoad", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(p%iUTight, 0, "p%iUTight", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(p%iUOpt1, 0, "p%iUOpt1", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(isLoadTight, 0, "isLoadTight", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(isLoadOption1, 0, "isLoadOption1", ErrStat2, ErrMsg2); if (Failed()) return
   do i = 1, size(Mods)
      do j = 1, size(Mods(i)%Vars%u)
         associate (Var => Mods(i)%Vars%u(j))
            if (iand(Var%Flags, VF_Solve) == 0) cycle
            isLoad = any(LoadFields == Var%Field)
            select case (Var%Cat)
            case (VC_Tight)
               p%iUTight = [p%iUTight, Var%iGblSol]
               isLoadTight = [isLoadTight, spread(isLoad, 1, Var%Size)]
            case (VC_Option1)
               p%iUOpt1 = [p%iUOpt1, Var%iGblSol]
               isLoadOption1 = [isLoadOption1, spread(isLoad, 1, Var%Size)]
            end select
            if (isLoad) then
               p%iuLoad = [p%iuLoad, Var%iGblSol]
            end if
         end associate
      end do
   end do

   ! Calculate index mapping arrays for y^Tight and y^Option1
   call AllocAry(p%iyTight, 0, "p%iyTight", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(p%iyOpt1, 0, "p%iyOpt1", ErrStat2, ErrMsg2); if (Failed()) return
   do i = 1, size(Mods)
      do j = 1, size(Mods(i)%Vars%y)
         if (iand(Mods(i)%Vars%y(j)%Flags, VF_Solve) == 0) cycle
         select case (Mods(i)%Vars%y(j)%Cat)
         case (VC_Tight)
            p%iyTight = [p%iyTight, Mods(i)%Vars%y(j)%iGblSol]
         case (VC_Option1)
            p%iyOpt1 = [p%iyOpt1, Mods(i)%Vars%y(j)%iGblSol]
         end select
      end do
   end do

   !----------------------------------------------------------------------------
   ! Allocate Jacobian matrices
   !----------------------------------------------------------------------------

   call AllocAry(m%dYdx, NumY, NumX, "m%dYdx", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%dYdu, NumY, NumU, "m%dYdu", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%dXdx, NumX, NumX, "m%dXdx", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%dXdu, NumX, NumU, "m%dXdu", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%dUdu, NumU, NumU, "m%dUdu", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%dUdy, NumU, NumY, "m%dUdy", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%G, NumU, NumU, "m%G", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%IPIV2, NumU, "m%IPIV2", ErrStat2, ErrMsg2); if (Failed()) return

   ! Initialize Jacobian matrices
   m%dYdx = 0.0_R8Ki
   m%dYdu = 0.0_R8Ki
   m%dXdx = 0.0_R8Ki
   m%dXdu = 0.0_R8Ki
   m%dUdu = 0.0_R8Ki
   m%dUdy = 0.0_R8Ki

   !----------------------------------------------------------------------------
   ! Allocate Jacobian matrix, RHS, and Delta
   !----------------------------------------------------------------------------

   ! Initialize size of system
   NumJac = 0

   ! Allocate and initialize indices of q vars in Jacobian matrix
   p%iJX2 = [(NumJac + i, i=1, size(p%iX2Tight))]
   NumJac = NumJac + size(p%iJX2)

   ! Allocate and initialize indices of tight coupling vars in Jacobian matrix
   p%iJT = [(NumJac + i, i=1, size(p%iUTight))]
   NumJac = NumJac + size(p%iJT)

   ! Allocate and initialize indices of option 1 vars in Jacobian matrix
   p%iJ1 = [(NumJac + i, i=1, size(p%iUOpt1))]
   NumJac = NumJac + size(p%iJ1)

   ! Get Jacobian indicies containing loads
   p%iJL = [pack(p%iJT, isLoadTight), pack(p%iJ1, isLoadOption1)]

   ! Allocate Macobian matrix, RHS/X matrix, Pivot array
   call AllocAry(m%Jac, NumJac, NumJac, "m%Jac", ErrStat, ErrMsg); if (Failed()) return
   call AllocAry(m%XB, NumJac, 1, "m%XB", ErrStat, ErrMsg); if (Failed()) return
   call AllocAry(m%IPIV, NumJac, "m%IPIV", ErrStat, ErrMsg); if (Failed()) return
   m%Jac = 0.0_R8Ki

   !----------------------------------------------------------------------------
   ! Calculate generalized alpha parameters
   !----------------------------------------------------------------------------

   p%AccBlend = 1.0_R8Ki

   p%AlphaM = (2.0_R8Ki*p%RhoInf - 1.0_R8Ki)/(p%RhoInf + 1.0_R8Ki)
   p%AlphaF = p%RhoInf/(p%RhoInf + 1.0_R8Ki)
   p%Gamma = 0.5_R8Ki - p%AlphaM + p%AlphaF
   p%Beta = (1.0_R8Ki - p%AlphaM + p%AlphaF)**2/4.0_R8Ki

   p%C(1) = (1.0_R8Ki - p%AlphaF)/(1.0_R8Ki - p%AlphaM)
   p%C(2) = p%DT*p%Gamma*p%C(1)
   p%C(3) = p%DT**2*p%Beta*p%C(1)

   !----------------------------------------------------------------------------
   ! Debug
   !----------------------------------------------------------------------------

   munit = -1

   call GetNewUnit(m%DebugUnit, ErrStat2, ErrMsg2); if (Failed()) return
   call OpenFOutFile(m%DebugUnit, "solver.dbg", ErrStat2, ErrMsg2); if (Failed()) return

   write (m%DebugUnit, *) "NumX      = ", NumX
   write (m%DebugUnit, *) "NumU      = ", NumU
   write (m%DebugUnit, *) "NumY      = ", NumY
   write (m%DebugUnit, *) "NumJac    = ", NumJac
   write (m%DebugUnit, '(A,*(I4))') " p%iJX2    = ", p%iJX2
   write (m%DebugUnit, '(A,*(I4))') " p%iJT     = ", p%iJT
   write (m%DebugUnit, '(A,*(I4))') " p%iJ1     = ", p%iJ1
   write (m%DebugUnit, '(A,*(I4))') " p%iJL     = ", p%iJL
   write (m%DebugUnit, '(A,*(I4))') " p%iX2Tight = ", p%iX2Tight
   write (m%DebugUnit, '(A,*(I4))') " p%iX1Tight = ", p%iX1Tight
   write (m%DebugUnit, '(A,*(I4))') " p%iUTight = ", p%iUTight
   write (m%DebugUnit, '(A,*(I4))') " p%iUOpt1  = ", p%iUOpt1
   write (m%DebugUnit, '(A,*(I4))') " p%iyTight = ", p%iyTight
   write (m%DebugUnit, '(A,*(I4))') " p%iyOpt1  = ", p%iyOpt1
   write (m%DebugUnit, *) "shape(m%dYdx) = ", shape(m%dYdx)
   write (m%DebugUnit, *) "shape(m%dYdu) = ", shape(m%dYdu)
   write (m%DebugUnit, *) "shape(m%dXdx) = ", shape(m%dXdx)
   write (m%DebugUnit, *) "shape(m%dXdu) = ", shape(m%dXdu)
   write (m%DebugUnit, *) "shape(m%dUdu) = ", shape(m%dUdu)
   write (m%DebugUnit, *) "shape(m%dUdy) = ", shape(m%dUdy)

   do i = 1, size(Mods)
      write (m%DebugUnit, *) "Module   = ", Mods(i)%Abbr
      write (m%DebugUnit, *) "ModuleID = ", Mods(i)%ID
      do j = 1, size(Mods(i)%Vars%x)
         if (.not. allocated(Mods(i)%Vars%x(j)%iGblSol)) cycle
         write (m%DebugUnit, *) "Var = "//trim(Mods(i)%Abbr)//trim(Num2LStr(Mods(i)%Ins))//" X "//trim(Mods(i)%Vars%x(j)%Name)// &
            " ("//trim(MV_FieldString(Mods(i)%Vars%x(j)%Field))//")"
         write (m%DebugUnit, '(A,*(I4))') "  X iLoc    = ", Mods(i)%Vars%x(j)%iLoc
         write (m%DebugUnit, '(A,*(I4))') "  X iGblSol = ", Mods(i)%Vars%x(j)%iGblSol
      end do
      do j = 1, size(Mods(i)%Vars%u)
         if (.not. allocated(Mods(i)%Vars%u(j)%iGblSol)) cycle
         write (m%DebugUnit, *) "Var = "//trim(Mods(i)%Abbr)//trim(Num2LStr(Mods(i)%Ins))//" U "//trim(Mods(i)%Vars%u(j)%Name)// &
            " ("//trim(MV_FieldString(Mods(i)%Vars%u(j)%Field))//")"
         write (m%DebugUnit, '(A,*(I4))') "  U iLoc    = ", Mods(i)%Vars%u(j)%iLoc
         write (m%DebugUnit, '(A,*(I4))') "  U iGblSol = ", Mods(i)%Vars%u(j)%iGblSol
      end do
      do j = 1, size(Mods(i)%Vars%y)
         if (.not. allocated(Mods(i)%Vars%y(j)%iGblSol)) cycle
         write (m%DebugUnit, *) "Var = "//trim(Mods(i)%Abbr)//trim(Num2LStr(Mods(i)%Ins))//" Y "//trim(Mods(i)%Vars%y(j)%Name)// &
            " ("//trim(MV_FieldString(Mods(i)%Vars%y(j)%Field))//")"
         write (m%DebugUnit, '(A,*(I4))') "  Y iLoc    = ", Mods(i)%Vars%y(j)%iLoc
         write (m%DebugUnit, '(A,*(I4))') "  Y iGblSol = ", Mods(i)%Vars%y(j)%iGblSol
      end do
   end do

contains
   function Failed()
      logical :: Failed
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function

   pure function VarCategory(ModID, VarField) result(VarCat)
      integer(IntKi), intent(in) :: ModID, VarField
      integer(IntKi)             :: VarCat
      ! If tight coupling module, then category is Tight
      if (any(ModID == TC_Modules)) then
         VarCat = VC_Tight
         return
      end if
      select case (VarField)
      case (VF_Orientation, VF_TransDisp, VF_AngularDisp)   ! Position
         VarCat = VC_Option2
      case (VF_TransVel, VF_AngularVel)                     ! Velocity
         VarCat = VC_Option2
      case (VF_TransAcc, VF_AngularAcc)                     ! Acceleration
         VarCat = VC_Option1
      case default
         VarCat = VC_None
      end select
   end function

end subroutine

subroutine Solver_DefineMappings(Mappings, Mods, ErrStat, ErrMsg)
   type(TC_MappingType), allocatable, intent(inout)   :: Mappings(:)
   type(ModDataType), intent(inout)                   :: Mods(:)  !< Solution variables from modules
   integer(IntKi), intent(out)                        :: ErrStat     !< Error status of the operation
   character(*), intent(out)                          :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   character(*), parameter   :: RoutineName = 'Solver_DefineMappings'
   integer(IntKi)            :: ErrStat2    ! local error status
   character(ErrMsgLen)      :: ErrMsg2     ! local error message
   integer(IntKi)            :: iMap, iModOut, iModIn, i, j
   logical, allocatable      :: isActive(:)

   ErrStat = ErrID_None
   ErrMsg = ''

   !----------------------------------------------------------------------------
   ! Define mesh mappings between modules
   !----------------------------------------------------------------------------

   ! Define a list of all possible module mesh mappings between modules
   ! Note: the mesh names must map those defined in MV_AddMeshVar in the modules
   allocate (Mappings(0), stat=ErrStat2)
   if (ErrStat2 /= 0) then
      call SetErrStat(ErrID_Fatal, "Error allocating mappings", ErrStat, ErrMsg, RoutineName)
      return
   end if

   do iMap = 1, size(Mods)
      if (Mods(iMap)%ID == Module_BD) then
         iModOut = Mods(iMap)%Ins
         call AddMotionMapping(Key='ED BladeRoot -> BD RootMotion', &
                               SrcModID=Module_ED, SrcIns=1, SrcMeshName='Blade root '//trim(Num2LStr(iModOut)), &
                               DstModID=Module_BD, DstIns=iModOut, DstMeshName='RootMotion')
         call AddLoadMapping(Key='BD ReactionForce -> ED HubLoad', &
                             SrcModID=Module_BD, SrcIns=iModOut, SrcMeshName='ReactionForce', SrcDispMeshName='RootMotion', &
                             DstModID=Module_ED, DstIns=1, DstMeshName='Hub', DstDispMeshName='Hub')
      end if
   end do

   !----------------------------------------------------------------------------
   ! Get module indices in ModData and determine which mappings are active
   !----------------------------------------------------------------------------

   ! Allocate array to indicate if mapping is active and initialize to false
   call AllocAry(isActive, size(Mappings), "isActive", ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   isActive = .false.

   ! Loop through mappings
   do iMap = 1, size(Mappings)

      ! Loop modules, if module ID matches source module ID
      ! and module instance matches source module instance, exit loop
      do iModOut = 1, size(Mods)
         if ((Mappings(iMap)%SrcModID == Mods(iModOut)%ID) .and. &
             (Mappings(iMap)%SrcIns == Mods(iModOut)%Ins)) exit
      end do

      ! Loop through modules, if module ID matches destination module ID
      ! and module instance matches destinatino module instance, exit loop
      do iModIn = 1, size(Mods)
         if ((Mappings(iMap)%DstModID == Mods(iModIn)%ID) .and. &
             (Mappings(iMap)%DstIns == Mods(iModIn)%Ins)) exit
      end do

      ! If input or output module not found, mapping is not active, cycle
      if (iModOut > size(Mods) .or. iModIn > size(Mods)) cycle

      ! Mark mapping as active
      isActive(iMap) = .true.

      ! Module input/ouput IDs and instances found, populate mapping
      Mappings(iMap)%SrcModIdx = iModOut
      Mappings(iMap)%DstModIdx = iModIn

      associate (map => Mappings(iMap), &
                 SrcMod => Mods(Mappings(iMap)%SrcModIdx), &
                 DstMod => Mods(Mappings(iMap)%DstModIdx))

         ! TODO: Add logic to check if mesh exists, skip mapping if it doesn't exist

         ! If load mapping
         if (map%IsLoad) then

            ! Source mesh variable indices
            map%SrcVarIdx = [(MV_VarIndex(SrcMod%Vars%y, map%SrcMeshName, LoadFields(iModOut)), iModOut=1, size(LoadFields))]
            map%SrcVarIdx = pack(map%SrcVarIdx, map%SrcVarIdx > 0)

            ! Destination mesh variable indices
            map%DstVarIdx = [(MV_VarIndex(DstMod%Vars%u, map%DstMeshName, LoadFields(iModOut)), iModOut=1, size(LoadFields))]
            map%DstVarIdx = pack(map%DstVarIdx, map%DstVarIdx > 0)

            ! Source displacement mesh is in input of source module (only translation displacement needed)
            map%SrcDispVarIdx = MV_VarIndex(SrcMod%Vars%u, map%SrcDispMeshName, VF_TransDisp)

            ! Destination displacement mesh is in output of destination module (only translation displacement needed)
            map%DstDispVarIdx = MV_VarIndex(DstMod%Vars%y, map%DstDispMeshName, VF_TransDisp)

            ! Mark displacement variables with Solve flag
            call SetFlags(SrcMod%Vars%u(map%SrcDispVarIdx), VF_Solve)
            call SetFlags(DstMod%Vars%y(map%DstDispVarIdx), VF_Solve)

         else

            ! Source mesh motion field variables
            map%SrcVarIdx = [(MV_VarIndex(SrcMod%Vars%y, map%SrcMeshName, MotionFields(iModOut)), iModOut=1, size(MotionFields))]
            map%SrcVarIdx = pack(map%SrcVarIdx, map%SrcVarIdx > 0)

            ! Destination mesh motion field variables
            map%DstVarIdx = [(MV_VarIndex(DstMod%Vars%u, map%DstMeshName, MotionFields(iModOut)), iModOut=1, size(MotionFields))]
            map%DstVarIdx = pack(map%DstVarIdx, map%DstVarIdx > 0)

         end if

         ! Mark variables with Solve flag
         do iModOut = 1, size(map%SrcVarIdx)
            call SetFlags(SrcMod%Vars%y(map%SrcVarIdx(iModOut)), VF_Solve)
         end do
         do iModOut = 1, size(map%DstVarIdx)
            call SetFlags(DstMod%Vars%u(map%DstVarIdx(iModOut)), VF_Solve)
         end do

      end associate

   end do

   ! Remove inactive mappings
   Mappings = pack(Mappings, mask=isActive)

contains
   subroutine AddLoadMapping(Key, SrcModID, SrcIns, SrcMeshName, SrcDispMeshName, &
                             DstModID, DstIns, DstMeshName, DstDispMeshName)
      character(*), intent(in)                           :: Key
      integer(IntKi), intent(in)                         :: SrcModID, DstModID
      integer(IntKi), intent(in)                         :: SrcIns, DstIns
      character(*), intent(in)                           :: SrcMeshName, DstMeshName
      character(*), intent(in)                           :: SrcDispMeshName, DstDispMeshName
      if (.not. allocated(Mappings)) allocate (Mappings(0))
      Mappings = [Mappings, TC_MappingType(Key=Key, isLoad=.true., &
                                           SrcModID=SrcModID, SrcIns=SrcIns, SrcMeshName=SrcMeshName, SrcDispMeshName=SrcDispMeshName, &
                                           DstModID=DstModID, DstIns=DstIns, DstMeshName=DstMeshName, DstDispMeshName=DstDispMeshName)]
   end subroutine
   subroutine AddMotionMapping(Key, SrcModID, SrcIns, SrcMeshName, &
                               DstModID, DstIns, DstMeshName)
      character(*), intent(in) :: Key
      integer(IntKi), intent(in) :: SrcModID, DstModID
      integer(IntKi), intent(in) :: SrcIns, DstIns
      character(*), intent(in) :: SrcMeshName, DstMeshName
      if (.not. allocated(Mappings)) allocate (Mappings(0))
      Mappings = [Mappings, TC_MappingType(Key=Key, isLoad=.false., &
                                           SrcModID=SrcModID, SrcIns=SrcIns, SrcMeshName=SrcMeshName, &
                                           DstModID=DstModID, DstIns=DstIns, DstMeshName=DstMeshName)]
   end subroutine
end subroutine

subroutine Solver_TransferXtoQ(ixqd, x, q)
   integer(IntKi), intent(in) :: ixqd(:, :)
   real(R8Ki), intent(in)     :: x(:)
   real(R8Ki), intent(inout)  :: q(:, :)
   integer(IntKi)             :: i
   do i = 1, size(ixqd, dim=2)
      q(ixqd(2, i), ixqd(3, i)) = x(ixqd(1, i))
   end do
end subroutine

subroutine Solver_TransferQtoX(ixqd, q, x)
   integer(IntKi), intent(in) :: ixqd(:, :)
   real(R8Ki), intent(in)     :: q(:, :)
   real(R8Ki), intent(inout)  :: x(:)
   integer(IntKi)             :: i
   do i = 1, size(ixqd, dim=2)
      x(ixqd(1, i)) = q(ixqd(2, i), ixqd(3, i))
   end do
end subroutine

pure function NeedWriteOutput(n_t_global, t_global, t_initial, n_DT_Out) result(WriteNeeded)
   integer(IntKi), intent(in)    :: n_t_global  !< Current global time step
   real(DbKi), intent(in)        :: t_initial   !< Initial time
   real(DbKi), intent(in)        :: t_global    !< Current global time
   integer(IntKi), intent(in)    :: n_DT_Out    !< Write output every n steps
   logical                       :: WriteNeeded !< Function result; if true, WriteOutput values are needed on this time step

   ! note that if TStart isn't an multiple of DT_out, we will not necessarially start output to the file at TStart
   if (t_global >= t_initial) then
      WriteNeeded = MOD(n_t_global, n_DT_Out) == 0
   else
      WriteNeeded = .false.
   end if
end function NeedWriteOutput

subroutine Solver_Step0(p, m, ModData, Turbine, ErrStat, ErrMsg)
   type(TC_ParameterType), intent(in)        :: p           !< Parameters
   type(TC_MiscVarType), intent(inout)       :: m           !< Misc variables
   type(ModDataType), intent(in)             :: ModData(:)  !< Solution variables from modules
   type(FAST_TurbineType), intent(inout)     :: Turbine     !< Turbine type
   integer(IntKi), intent(out)               :: ErrStat
   character(*), intent(out)                 :: ErrMsg

   character(*), parameter    :: RoutineName = 'Solver_Step0'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i, j, k
   real(R8Ki), allocatable    :: accel(:)
   real(R8Ki)                 :: diff
   logical                    :: converged
   integer(IntKi), parameter  :: n_t_global = -1     ! loop counter
   integer(IntKi), parameter  :: n_t_global_next = 0 ! loop counter
   real(DbKi)                 :: t_initial           ! next simulation time
   real(DbKi)                 :: t_global_next       ! next simulation time

   ErrStat = ErrID_None
   ErrMsg = ''

   !----------------------------------------------------------------------------
   ! Initialize module mappings
   ! TODO: Move this into init
   !----------------------------------------------------------------------------

   call FAST_InitMappings(m%Mappings, Turbine, ErrStat2, ErrMsg2); if (Failed()) return

   !----------------------------------------------------------------------------
   ! Miscellaneous initial step
   !----------------------------------------------------------------------------

   t_initial = Turbine%m_FAST%t_global
   t_global_next = t_initial + n_t_global_next*p%DT
   Turbine%y_FAST%WriteThisStep = NeedWriteOutput(n_t_global_next, t_global_next, t_initial, Turbine%p_FAST%n_DT_Out)

   if (Turbine%p_FAST%WrSttsTime) then
      call SimStatus_FirstTime(Turbine%m_FAST%TiLstPrn, Turbine%m_FAST%PrevClockTime, &
                               Turbine%m_FAST%SimStrtTime, Turbine%m_FAST%UsrTime2, t_initial, &
                               Turbine%p_FAST%TMax, Turbine%p_FAST%TDesc)
   end if

   !----------------------------------------------------------------------------
   ! Calculate initial accelerations
   !----------------------------------------------------------------------------

   ! Transfer initial state from modules to solver
   call PackModuleStates(ModData(p%iModTC), STATE_CURR, Turbine, x=m%x)

   ! Transfer initial state to state q matrix
   call Solver_TransferXtoQ(p%ixqd, m%x, m%qn)

   ! Allocate acceleration array which will be used to check for convergence
   ! of initial acceleration. Transfer initial accelerations from q matrix
   call AllocAry(accel, size(m%qn, dim=2), "accel", ErrStat2, ErrMsg2); if (Failed()) return
   accel = m%qn(:, COL_A)

   ! Reset mappings updated flags
   m%Mappings%Updated = .false.

   ! Loop until initial accelerations are converged, or max iterations are reached.
   ! TODO: may need a separate variable for max initial acceleration convergence iterations
   converged = .false.
   k = 1
   do while ((.not. converged) .and. (k <= p%MaxConvIter))

      ! Transfer inputs and calculate outputs for all modules (use current state)
      do i = 1, size(p%iModAll)
         call FAST_UpdateInputs(ModData(p%iModAll(i)), m%Mappings, 1, &
                                Turbine, ErrStat2, ErrMsg2); if (Failed()) return
         call FAST_CalcOutput(ModData(p%iModAll(i)), t_initial, STATE_CURR, &
                              Turbine, ErrStat2, ErrMsg2); if (Failed()) return
         call FAST_MapOutputs(ModData(p%iModAll(i)), m%Mappings, &
                              Turbine, ErrStat2, ErrMsg2); if (Failed()) return
      end do

      ! Calculate continuous state derivatives for tight coupling modules (use current state)
      do i = 1, size(p%iModTC)
         call FAST_CalcContStateDeriv(ModData(p%iModTC(i)), t_initial, STATE_CURR, &
                                      Turbine, ErrStat2, ErrMsg2, dxdt=m%dxdt); if (Failed()) return
      end do

      ! Copy acceleration (derivative of velocity) to acceleration array.
      ! Loop through x->q index array and transfer derivative of
      ! velocity variables (deriv=COL_V) to the acceleration column
      do i = 1, size(p%ixqd, dim=2)
         if (p%ixqd(3, i) == COL_V) accel(p%ixqd(2, i)) = m%dxdt(p%ixqd(1, i))
      end do

      ! Calculate diff as L2 norm of current and new accelerations
      diff = TwoNorm(accel - m%qn(:, COL_A))

      ! If difference is less than converence tolerance, set flag and exit loop
      if ((k > 1) .and. (diff < p%ConvTol)) converged = .true.

      ! Update acceleration in q matrix
      m%qn(:, COL_A) = accel

      ! Increment iteration counter
      k = k + 1
   end do

   ! Print warning if not converged
   if (.not. converged) then
      call WrScr("Solver: initial accel not converged, diff="//Num2LStr(diff)// &
                 ", tol="//Num2LStr(p%ConvTol))
   end if

   ! Initialize algorithmic acceleration from actual acceleration
   m%qn(:, COL_AA) = m%qn(:, COL_A)

   !----------------------------------------------------------------------------
   ! Initialize module input and state arrays for interpolation/extrapolation
   !----------------------------------------------------------------------------

   ! Loop through all module index array
   do i = 1, size(p%iModAll)

      ! Initialize IO and states for all modules (also copies STATE_CURR to STATE_PRED)
      call FAST_InitIO(ModData(p%iModAll(i)), t_initial, p%DT, Turbine, ErrStat2, ErrMsg2); if (Failed()) return

      ! Reset the Remap flags for all modules
      call FAST_ResetRemapFlags(ModData(p%iModAll(i)), Turbine, ErrStat2, ErrMsg2); if (Failed()) return
   end do

contains
   function Failed()
      logical :: Failed
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine UpdateBeamDynGlobalReference(p, m, Mod, T, ErrStat, ErrMsg)
   type(TC_ParameterType), intent(in)     :: p        !< Parameters
   type(TC_MiscVarType), intent(inout)    :: m        !< Misc variables
   type(ModDataType), intent(in)          :: Mod
   type(FAST_TurbineType), intent(inout)  :: T        !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'UpdateBeamDynGlobalReference'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   real(R8Ki)                 :: GlbWM_old(3), GlbWM_new(3), GlbWM_diff(3)
   real(R8Ki)                 :: GlbRot_old(3, 3), GlbRot_new(3, 3), GlbRot_diff(3, 3)
   real(R8Ki)                 :: GlbPos_old(3), GlbPos_new(3), GlbPos_diff(3)
   real(R8Ki)                 :: pos(3), rot(3), trans_vel(3), rot_vel(3), uuN0(3)
   integer(IntKi)             :: i, j, temp_id, temp_id2

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Save old global position, rotation, and WM
   GlbPos_old = T%BD%p(Mod%Ins)%GlbPos
   GlbRot_old = T%BD%p(Mod%Ins)%GlbRot
   GlbWM_old = T%BD%p(Mod%Ins)%Glb_crv

   ! Calculate new global position, rotation, and WM from root motion
   GlbPos_new = T%BD%Input(1, Mod%Ins)%RootMotion%Position(:, 1) + &
                T%BD%Input(1, Mod%Ins)%RootMotion%TranslationDisp(:, 1)
   GlbRot_new = transpose(T%BD%Input(1, Mod%Ins)%RootMotion%Orientation(:, :, 1))
   GlbWM_new = wm_from_dcm(GlbRot_new)

   ! Update the module global values
   T%BD%p(Mod%Ins)%GlbPos = GlbPos_new
   T%BD%p(Mod%Ins)%GlbRot = GlbRot_new
   T%BD%p(Mod%Ins)%Glb_crv = GlbWM_new

   ! Calculate differences between old and new reference
   GlbRot_diff = matmul(transpose(GlbRot_old), GlbRot_new)
   GlbWM_diff = wm_compose(wm_inv(GlbWM_old), GlbWM_new)
   GlbPos_diff = GlbPos_old - GlbPos_new

   associate (x_BD => T%BD%x(Mod%Ins, STATE_PRED), p_BD => T%BD%p(Mod%Ins))

      x_BD%q(:, 1) = 0.0_R8Ki
      x_BD%dqdt(1:3, 1) = matmul(transpose(GlbRot_diff), T%BD%Input(1, Mod%Ins)%RootMotion%TranslationVel(:, 1))
      x_BD%dqdt(4:6, 1) = matmul(transpose(GlbRot_diff), T%BD%Input(1, Mod%Ins)%RootMotion%RotationVel(:, 1))

      do i = 1, p_BD%elem_total
         do j = 1, p_BD%nodes_per_elem

            temp_id = (i - 1)*(p_BD%nodes_per_elem - 1) + j ! The last node of the first element is used as the first node in the second element.
            temp_id2 = (i - 1)*p_BD%nodes_per_elem + j      ! Index to a node within element i

            ! Calculate displacements from new reference
            x_BD%q(1:3, temp_id) = matmul(transpose(GlbRot_new), &
                                          GlbPos_old + matmul(GlbRot_old, p_BD%uuN0(1:3, j, i) + x_BD%q(1:3, temp_id)) - &
                                          GlbPos_new - matmul(GlbRot_new, p_BD%uuN0(1:3, j, i)))

            ! Update the node orientation rotation of the node
            x_BD%q(4:6, temp_id) = wm_compose(wm_inv(GlbWM_diff), x_BD%q(4:6, temp_id))

            ! Update the translational velocity
            x_BD%dqdt(1:3, temp_id) = matmul(transpose(GlbRot_diff), x_BD%dqdt(1:3, temp_id))

            ! Update the rotational velocity
            x_BD%dqdt(4:6, temp_id) = matmul(transpose(GlbRot_diff), x_BD%dqdt(4:6, temp_id))

         end do
      end do

   end associate

   ! T%BD%x(Mod%Ins, STATE_PRED)%q = 0
   ! T%BD%x(Mod%Ins, STATE_PRED)%dqdt = 0

   call BD_PackStateValues(T%BD%p(Mod%Ins), T%BD%x(Mod%Ins, STATE_PRED), T%BD%m(Mod%Ins)%Vals%x)
   call XferLocToGbl1D(Mod%ixs, T%BD%m(Mod%Ins)%Vals%x, m%xn)
   call Solver_TransferXtoQ(p%ixqd, m%xn, m%qn)

end subroutine

subroutine Solver_Step(n_t_global, t_initial, p, m, Mods, Turbine, ErrStat, ErrMsg)
   integer(IntKi), intent(in)                :: n_t_global  !< global time step
   real(DbKi), intent(in)                    :: t_initial   !< Initial simulation time
   type(TC_ParameterType), intent(in)        :: p           !< Parameters
   type(TC_MiscVarType), intent(inout)       :: m           !< Misc variables
   type(ModDataType), intent(in)             :: Mods(:)  !< Solution variables from modules
   type(FAST_TurbineType), intent(inout)     :: Turbine     !< Turbine type
   integer(IntKi), intent(out)               :: ErrStat
   character(*), intent(out)                 :: ErrMsg

   character(*), parameter    :: RoutineName = 'Solver_Step'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: NumCorrections      ! number of corrections for this time step
   integer(IntKi)             :: iterConv, iterCorr
   logical                    :: update_jacobian
   real(ReKi)                 :: delta_norm
   real(DbKi)                 :: t_global_next       ! next simulation time (m_FAST%t_global + p_FAST%dt)
   integer(IntKi)             :: n_t_global_next     ! n_t_global + 1
   integer(IntKi)             :: i

   ErrStat = ErrID_None
   ErrMsg = ''

   !----------------------------------------------------------------------------
   ! Miscellaneous step updates
   !----------------------------------------------------------------------------

   ! Calculate the next global time step number and time
   n_t_global_next = n_t_global + 1
   t_global_next = t_initial + n_t_global_next*p%DT

   ! Determine if output should be written in this step
   Turbine%y_FAST%WriteThisStep = NeedWriteOutput(n_t_global_next, t_global_next, t_initial, Turbine%p_FAST%n_DT_Out)

   ! Decrement number of time steps before updating the Jacobian
   m%StepsUntilUJac = m%StepsUntilUJac - 1

   write (m%DebugUnit, *) "step = ", n_t_global_next

   !----------------------------------------------------------------------------
   ! Extrapolate/interpolate inputs for all modules
   !----------------------------------------------------------------------------

   ! Loop through all modules and extrap/interp inputs
   do i = 1, size(p%iModAll)
      call FAST_ExtrapInterp(Mods(p%iModAll(i)), t_global_next, Turbine, ErrStat2, ErrMsg2); if (Failed()) return
   end do

   !----------------------------------------------------------------------------
   ! Prediction - guess solution state variables at end of time step
   ! m%qn contains the states at the end of the last step.
   ! m%q contains the prediction which is copied to m%qn at the start of the
   ! correction loop.
   !----------------------------------------------------------------------------

   ! Acceleration for next step
   m%q(:, COL_A) = p%AccBlend*m%qn(:, COL_A)

   ! Algorithm acceleration for next step
   m%q(:, COL_AA) = ((1.0_R8Ki - p%AlphaF)*m%q(:, COL_A) + &
                     p%AlphaF*m%qn(:, COL_A) - &
                     p%AlphaM*m%qn(:, COL_AA))/(1 - p%AlphaM)

   ! Calculate change in position and velocities
   ! (position states include orientations which must be composed with deltas)
   m%dq = 0.0_R8Ki
   m%dq(:, COL_V) = p%DT*(1 - p%Gamma)*m%qn(:, COL_AA) + p%DT*p%Gamma*m%q(:, COL_AA)
   m%dq(:, COL_D) = p%DT*m%qn(:, COL_V) + p%DT**2*(0.5 - p%Beta)*m%qn(:, COL_AA) + p%DT**2*p%Beta*m%q(:, COL_AA)

   ! Transfer delta state matrix to delta x array
   m%dx = 0.0_R8Ki
   call Solver_TransferQtoX(p%ixqd, m%dq, m%dx)

   ! Add delta to x array to get new states (respect variable fields)
   ! required for orientation fields in states
   call AddDeltaToStates(Mods, p%iModTC, m%dx, m%x)

   ! Update state matrix with updated state values
   call Solver_TransferXtoQ(p%ixqd, m%x, m%q)

   ! Initialize delta arrays for iterations
   m%dq = 0.0_R8Ki
   m%dx = 0.0_R8Ki
   m%du = 0.0_R8Ki

   !----------------------------------------------------------------------------
   ! Correction Iterations
   !----------------------------------------------------------------------------

   ! Loop through correction iterations
   iterCorr = 0
   do while (iterCorr <= p%NumCrctn)

      write (m%DebugUnit, *) "iterCorr = ", iterCorr

      ! Copy state for correction step
      m%qn = m%q
      m%xn = m%x

      ! Reset mappings updated flags
      m%Mappings%Updated = .false.

      ! Loop through Option 2 modules
      do i = 1, size(p%iModOpt2)
         call FAST_UpdateInputs(Mods(p%iModOpt2(i)), m%Mappings, 1, &
                                Turbine, ErrStat2, ErrMsg2); if (Failed()) return
         call FAST_UpdateStates(Mods(p%iModOpt2(i)), t_initial, n_t_global, &
                                Turbine, ErrStat2, ErrMsg2, m%xn); if (Failed()) return
         call FAST_CalcOutput(Mods(p%iModOpt2(i)), t_global_next, STATE_PRED, &
                              Turbine, ErrStat2, ErrMsg2); if (Failed()) return
         if (i < 2) then
            call FAST_MapOutputs(Mods(p%iModOpt2(i)), m%Mappings, &
                                 Turbine, ErrStat2, ErrMsg2); if (Failed()) return
         end if
      end do

      ! Get inputs and update states for Option 1 modules not in Option 2
      do i = 1, size(p%iModOpt1US)
         call FAST_UpdateInputs(Mods(p%iModOpt1US(i)), m%Mappings, 1, &
                                Turbine, ErrStat2, ErrMsg2); if (Failed()) return
         call FAST_UpdateStates(Mods(p%iModOpt1US(i)), t_initial, n_t_global, &
                                Turbine, ErrStat2, ErrMsg2); if (Failed()) return
      end do

      ! Pack Option 1 inputs into u array
      call PackModuleInputs(Mods, p%iModOpt1, Turbine, u=m%un)

      !-------------------------------------------------------------------------
      ! Convergence Iterations
      !-------------------------------------------------------------------------

      ! Loop through convergence iterations
      do iterConv = 0, p%MaxConvIter

         ! Decrement number of iterations before updating the Jacobian
         m%IterUntilUJac = m%IterUntilUJac - 1

         !----------------------------------------------------------------------
         ! Option 1 calculate outputs
         !----------------------------------------------------------------------

         do i = 1, size(p%iModOpt1)
            call FAST_CalcOutput(Mods(p%iModOpt1(i)), t_global_next, STATE_PRED, &
                                 Turbine, ErrStat2, ErrMsg2); if (Failed()) return
         end do
         ! call PackModuleOutputs(Mods, p%iModOpt1, Turbine, m%y)

         !----------------------------------------------------------------------
         ! If iteration limit reached, exit loop
         !----------------------------------------------------------------------

         if (iterConv >= p%MaxConvIter) exit

         write (m%DebugUnit, *) "iterConv = ", iterConv

         write (m%DebugUnit, '(A,*(ES16.7))') " BD1-eps  = ", pack(Turbine%BD%m(1)%qp%E1(1:3, :, 1) - Turbine%BD%m(1)%qp%RR0(1:3, 3, :, 1), .true.)
         ! write (m%DebugUnit, '(A,*(ES16.7))') " BD2-eps  = ", pack(Turbine%BD%m(2)%qp%E1(1:3,:,1) - Turbine%BD%m(1)%qp%RR0(1:3,3,:,1), .true.)
         write (m%DebugUnit, '(A,*(ES16.7))') " BD1-kappa  = ", pack(Turbine%BD%m(1)%qp%kappa(1:3, :, 1), .true.)
         ! write (m%DebugUnit, '(A,*(ES16.7))') " BD2-kappa  = ", pack(Turbine%BD%m(2)%qp%kappa(1:3,:,1), .true.)
         write (m%DebugUnit, '(A,*(ES16.7))') " BD1-Nrrr  = ", pack(Turbine%BD%m(1)%Nrrr(1:3, :, 1), .true.)
         ! write (m%DebugUnit, '(A,*(ES16.7))') " BD2-Nrrr  = ", pack(Turbine%BD%m(2)%Nrrr(1:3,:,1), .true.)
         write (m%DebugUnit, '(A,*(ES16.7))') " BD1-Glb_crv  = ", Turbine%BD%p(1)%Glb_crv
         ! write (m%DebugUnit, '(A,*(ES16.7))') " BD2-Glb_crv  = ", Turbine%BD%p(2)%Glb_crv
         write (m%DebugUnit, '(A,*(ES16.7))') " BD1-RRoot  = ", wm_from_dcm(Turbine%BD%Input(1, 1)%RootMotion%Orientation(:, :, 1))
         ! write (m%DebugUnit, '(A,*(ES16.7))') " BD2-RRoot  = ", wm_from_dcm(Turbine%BD%Input(1,2)%RootMotion%Orientation(:,:,1))
         write (m%DebugUnit, '(A,*(ES16.7))') " BD1-RR  = ", wm_compose(wm_inv(Turbine%BD%p(1)%Glb_crv), wm_from_dcm(Turbine%BD%Input(1, 1)%RootMotion%Orientation(:, :, 1)))
         ! write (m%DebugUnit, '(A,*(ES16.7))') " BD2-RR  = ", wm_compose(wm_inv(Turbine%BD%p(2)%Glb_crv), wm_from_dcm(Turbine%BD%Input(1,2)%RootMotion%Orientation(:,:,1)))
         ! write (m%DebugUnit, '(A,*(ES16.7))') " BD1-RRoot-dcm  = ", pack(Turbine%BD%Input(1,1)%RootMotion%Orientation(:,:,1), .true.)
         ! write (m%DebugUnit, '(A,*(ES16.7))') " BD2-RRoot-dcm  = ", pack(Turbine%BD%Input(1,2)%RootMotion%Orientation(:,:,1), .true.)

         !----------------------------------------------------------------------
         ! Update Jacobian
         !----------------------------------------------------------------------

         ! If number of iterations or steps until Jacobian is to be updated
         ! is zero or less, then rebuild the Jacobian. Note: Solver_BuildJacobian
         ! resets these counters.
         if ((m%IterUntilUJac <= 0) .or. (m%StepsUntilUJac <= 0) .or. (n_t_global_next == 1)) then

            call Solver_BuildJacobian(p, m, Mods, t_global_next, n_t_global_next*100 + iterConv, &
                                      Turbine, ErrStat2, ErrMsg2); if (Failed()) return
         end if

         !----------------------------------------------------------------------
         ! Formulate right hand side (X_2^tight, U^tight, U^Option1)
         !----------------------------------------------------------------------

         ! Calculate continuous state derivatives for tight coupling modules
         do i = 1, size(p%iModTC)
            call FAST_CalcContStateDeriv(Mods(p%iModTC(i)), t_global_next, STATE_PRED, &
                                         Turbine, ErrStat2, ErrMsg2, dxdt=m%dxdt); if (Failed()) return
         end do

         ! Calculate difference between predicted and actual accelerations
         m%XB(p%iJX2, 1) = m%qn(:, COL_A) - m%dxdt(p%iX2Tight)

         ! Transfer Option 1 outputs to temporary inputs and collect into u_tmp
         do i = 1, size(p%iModOpt1)
            call FAST_MapOutputs(Mods(p%iModOpt1(i)), m%Mappings, &
                                 Turbine, ErrStat2, ErrMsg2); if (Failed()) return
            call FAST_UpdateInputs(Mods(p%iModOpt1(i)), m%Mappings, 2, &
                                   Turbine, ErrStat2, ErrMsg2); if (Failed()) return
         end do
         call PackModuleInputs(Mods, p%iModOpt1, Turbine, u_tmp=m%u_tmp)

         ! Calculate difference in U for all Option 1 modules (un - u_tmp)
         ! and add to RHS for TC and Option 1 modules
         call ComputeDiffU(Mods, p%iModOpt1, m%un, m%u_tmp, m%UDiff)
         m%XB(p%iJT, 1) = m%UDiff(p%iUTight)
         m%XB(p%iJ1, 1) = m%UDiff(p%iUOpt1)

         ! Apply conditioning factor to loads in RHS
         m%XB(p%iJL, 1) = m%XB(p%iJL, 1)/p%Scale_UJac

         !----------------------------------------------------------------------
         ! Solve for state and input perturbations
         !----------------------------------------------------------------------

         write (m%DebugUnit, '(A,*(ES16.7))') " XB = ", m%XB

         ! Solve Jacobian and RHS
         call LAPACK_getrs('N', size(m%Jac, 1), m%Jac, m%IPIV, m%XB, ErrStat2, ErrMsg2); if (Failed()) return

         !----------------------------------------------------------------------
         ! Check perturbations for convergence and exit if below tolerance
         !----------------------------------------------------------------------

         delta_norm = TwoNorm(m%XB(:, 1))/size(m%XB)

         write (m%DebugUnit, '(A,*(ES16.7))') " y  = ", m%y
         write (m%DebugUnit, '(A,*(ES16.7))') " u  = ", m%un
         write (m%DebugUnit, '(A,*(ES16.7))') " u_tmp = ", m%u_tmp
         write (m%DebugUnit, '(A,*(ES16.7))') " U  = ", m%UDiff
         write (m%DebugUnit, '(A,*(ES16.7))') " x  = ", m%xn
         write (m%DebugUnit, *) "delta_norm = ", delta_norm

         if (delta_norm < p%ConvTol) exit

         ! Remove conditioning
         m%XB(p%iJL, 1) = m%XB(p%iJL, 1)*p%Scale_UJac

         !----------------------------------------------------------------------
         ! Update State for Tight Coupling modules
         !----------------------------------------------------------------------

         ! Calculate change in state matrix
         m%dq(:, COL_D) = -p%C(3)*m%XB(p%iJX2, 1)
         m%dq(:, COL_V) = -p%C(2)*m%XB(p%iJX2, 1)
         m%dq(:, COL_A) = -m%XB(p%iJX2, 1)
         m%dq(:, COL_AA) = -p%C(1)*m%XB(p%iJX2, 1)

         ! Transfer change in q state matrix to change in x array
         call Solver_TransferQtoX(p%ixqd, m%dq, m%dx)

         ! Add delta to x array to get new states (respect variable fields)
         call AddDeltaToStates(Mods, p%iModTC, m%dx, m%xn)

         ! Update new state matrix with new state array values
         ! The transfer overwrites values that were changed in AddDeltaToStates
         m%qn = m%qn + m%dq
         call Solver_TransferXtoQ(p%ixqd, m%xn, m%qn)

         ! Transfer updated state to TC modules
         call UnpackModuleStates(Mods, p%iModTC, STATE_PRED, Turbine, x=m%xn)

         !----------------------------------------------------------------------
         ! Update inputs for Tight Coupling and Option 1 modules
         !----------------------------------------------------------------------

         ! Update change in inputs
         m%du(p%iUTight) = -m%XB(p%iJT, 1)
         m%du(p%iUOpt1) = -m%XB(p%iJ1, 1)

         ! Apply deltas to inputs, update modules
         call AddDeltaToInputs(Mods, p%iModOpt1, m%du, m%un)

         ! Transfer updated inputs to Option 1 modules
         call UnpackModuleInputs(Mods, p%iModOpt1, Turbine, u=m%un)

         !----------------------------------------------------------------------
         ! Transfer updated states and inputs to relevant modules
         !----------------------------------------------------------------------

         write (m%DebugUnit, '(A,*(ES16.7))') " du = ", m%du
         write (m%DebugUnit, '(A,*(ES16.7))') " dx = ", m%dx

      end do

      iterCorr = iterCorr + 1

      ! Reset the remap flags on the meshes
      do i = 1, size(p%iModAll)
         call FAST_ResetRemapFlags(Mods(p%iModAll(i)), Turbine, ErrStat2, ErrMsg2); if (Failed()) return
      end do

   end do

   !----------------------------------------------------------------------------
   ! Update states for next step
   !----------------------------------------------------------------------------

   ! Reset BeamDyn global reference and rescale BeamDyn states
   ! Loop through tight coupling modules, if module is BeamDyn, update global ref
   do i = 1, size(p%iModTC)
      if (Mods(p%iModTC(i))%ID == Module_BD) then
         call UpdateBeamDynGlobalReference(p, m, Mods(p%iModTC(i)), Turbine, ErrStat2, ErrMsg2); if (Failed()) return
      end if
   end do

   ! Copy the final predicted states from step t_global_next to actual states for that step
   do i = 1, size(p%iModAll)
      call FAST_SaveStates(Mods(p%iModAll(i)), Turbine, ErrStat2, ErrMsg2); if (Failed()) return
   end do

   ! Save new state
   m%x = m%xn

   ! Update the global time
   Turbine%m_FAST%t_global = t_global_next

   !----------------------------------------------------------------------------
   ! Display simulation status every SttsTime-seconds (i.e., n_SttsTime steps):
   !----------------------------------------------------------------------------

   if (Turbine%p_FAST%WrSttsTime) then
      if (MOD(n_t_global_next, Turbine%p_FAST%n_SttsTime) == 0) then
         call SimStatus(Turbine%m_FAST%TiLstPrn, Turbine%m_FAST%PrevClockTime, &
                        Turbine%m_FAST%t_global, Turbine%p_FAST%TMax, Turbine%p_FAST%TDesc)
      end if
   end if

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine ComputeDiffU(Mods, ModOrder, PosAry, NegAry, DiffAry)
   type(ModDataType), intent(in) :: Mods(:)      ! Array of variables
   integer(IntKi), intent(in)    :: ModOrder(:)
   real(R8Ki), intent(in)        :: PosAry(:)      ! Positive result array
   real(R8Ki), intent(in)        :: NegAry(:)      ! Negative result array
   real(R8Ki), intent(inout)     :: DiffAry(:)     ! Array containing difference
   integer(IntKi)                :: i, j, k, ind(3)
   real(R8Ki)                    :: DeltaWM(3), n(3), phi

   ! Loop through module index order
   do i = 1, size(ModOrder)

      ! Loop through input variables in module
      do j = 1, size(Mods(ModOrder(i))%Vars%u)

         associate (Var => Mods(ModOrder(i))%Vars%u(j))

            if (.not. allocated(Var%iGblSol)) cycle

            ! If variable field is orientation
            if (Var%Field == VF_Orientation) then

               ! Loop through nodes
               do k = 1, Var%Nodes

                  ! Get vector of indicies of WM rotation parameters in array
                  ind = Var%iGblSol(3*(k - 1) + 1:3*k)

                  ! Compose WM parameters to go from negative to positive array
                  DeltaWM = wm_compose(wm_inv(NegAry(ind)), PosAry(ind))

                  ! Calculate change in rotation in XYZ in radians
                  ! phi = TwoNorm(DeltaWM)
                  ! n = DeltaWM/phi
                  ! DiffAry(ind) = 4.0_R8Ki*atan(phi/4.0_R8Ki)*n
                  ! DiffAry(ind) = 4.0_R8Ki*atan(DeltaWM/4.0_R8Ki)
                  DiffAry(ind) = DeltaWM
               end do

            else

               ! Subtract negative array from positive array
               DiffAry(Var%iGblSol) = PosAry(Var%iGblSol) - NegAry(Var%iGblSol)
            end if
         end associate
      end do
   end do
end subroutine

subroutine Solver_BuildJacobian(p, m, Mods, this_time, iter, Turbine, ErrStat, ErrMsg)
   type(TC_ParameterType), intent(in)     :: p           !< Parameters
   type(TC_MiscVarType), intent(inOUT)    :: m           !< Misc variables
   type(ModDataType), intent(in)          :: Mods(:)     !< Array of module data
   real(DbKi), intent(in)                 :: this_time   !< Time
   integer(IntKi), intent(in)             :: iter
   type(FAST_TurbineType), intent(inout)  :: Turbine     !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter                :: RoutineName = 'Solver_BuildJacobian'
   integer(IntKi)                         :: ErrStat2
   character(ErrMsgLen)                   :: ErrMsg2
   integer(IntKi)                         :: i, j
   real(R8Ki), allocatable                :: tmp(:, :)
   real(R8Ki), dimension(3)               :: wm_b, wm_p, wm_n, wm_d, wm_pert, delta

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Reset Jacobian update countdown values
   m%IterUntilUJac = p%NIter_UJac
   m%StepsUntilUJac = p%NStep_UJac

   !----------------------------------------------------------------------------
   ! Get module Jacobians and assemble
   !----------------------------------------------------------------------------

   ! Initialize matrices
   m%dYdx = 0.0_R8Ki
   m%dXdx = 0.0_R8Ki
   m%dXdu = 0.0_R8Ki
   m%dYdu = 0.0_R8Ki
   m%dUdy = 0.0_R8Ki
   call Eye2D(m%dUdu, ErrStat2, ErrMsg2); if (Failed()) return

   ! Calculate dYdx, dXdx, dXdu for tight coupling modules
   do i = 1, size(p%iModTC)
      call FAST_CalcJacobian(Mods(p%iModTC(i)), this_time, STATE_PRED, Turbine, ErrStat2, ErrMsg2, &
                             dYdx=m%dYdx, dXdx=m%dXdx, dXdu=m%dXdu); if (Failed()) return
   end do

   ! Calculate dYdu Loop for Option 1 modules
   do i = 1, size(p%iModOpt1)
      call FAST_CalcJacobian(Mods(p%iModOpt1(i)), this_time, STATE_PRED, Turbine, ErrStat2, ErrMsg2, &
                             dYdu=m%dYdu); if (Failed()) return
   end do

   ! Calculate dUdu and dUdy for Option 1 meshes
   call FAST_LinearizeMappings(Mods, p%iModOpt1, m%Mappings, Turbine, ErrStat2, ErrMsg2, &
                               dUdu=m%dUdu, dUdy=m%dUdy); if (Failed()) return

   !----------------------------------------------------------------------------
   ! Form full system matrices
   !----------------------------------------------------------------------------

   !A: rows = x; columns = x (dXdx)
   !B: rows = x; columns = u (dXdu)
   !C: rows = y; columns = x (dYdx)
   !D: rows = y; columns = u (dYdu)

   ! G
   ! m%G = m%dUdu + matmul(m%dUdy, m%dYdu)
   m%G = m%dUdu
   call LAPACK_GEMM('N', 'N', 1.0_R8Ki, m%dUdy, m%dYdu, 1.0_R8Ki, m%G, ErrStat2, ErrMsg2); if (Failed()) return
   !----------------------------------------------------------------------------
   ! Assemble Jacobian
   !----------------------------------------------------------------------------

   ! Group (1,1)
   m%Jac(p%iJX2, p%iJX2) = -p%C(2)*m%dXdx(p%iX2Tight, p%iX2Tight) - p%C(3)*m%dXdx(p%iX2Tight, p%iX1Tight)
   do i = 1, size(p%iJX2)
      m%Jac(p%iJX2(i), p%iJX2(i)) = m%Jac(p%iJX2(i), p%iJX2(i)) + 1.0_R8Ki
   end do

   ! Group (2,1)
   m%Jac(p%iJT, p%iJX2) = p%C(2)*matmul(m%dUdy(p%iUTight, p%iyTight), m%dYdx(p%iyTight, p%iX2Tight)) + &
                          p%C(3)*matmul(m%dUdy(p%iUTight, p%iyTight), m%dYdx(p%iyTight, p%iX1Tight))

   ! Group (1,2)
   m%Jac(p%iJX2, p%iJT) = -m%dXdu(p%iX2Tight, p%iUTight)

   ! Group (2,2)
   m%Jac(p%iJT, p%iJT) = m%G(p%iUTight, p%iUTight)
   ! m%Jac(p%iJT, p%iJT) = m%dUdu(p%iUTight, p%iUTight) + &
   !                       matmul(m%dUdy(p%iUTight, p%iyTight), m%dYdu(p%iyTight, p%iUTight))

   ! If modules in option 1
   if (size(p%iJ1) > 0) then

      ! Group (3,2)
      m%Jac(p%iJ1, p%iJT) = m%dUdu(p%iUOpt1, p%iUTight) + &
                            matmul(m%dUdy(p%iUOpt1, p%iyTight), m%dYdu(p%iyTight, p%iUTight))
      ! Group (2,3)
      m%Jac(p%iJT, p%iJ1) = m%dUdu(p%iUTight, p%iUOpt1) + &
                            matmul(m%dUdy(p%iUTight, p%iyOpt1), m%dYdu(p%iyOpt1, p%iUOpt1))
      ! Group (3,3)
      m%Jac(p%iJ1, p%iJ1) = m%dUdu(p%iUOpt1, p%iUOpt1) + &
                            matmul(m%dUdy(p%iUOpt1, p%iyOpt1), m%dYdu(p%iyOpt1, p%iUOpt1))
   end if

   ! Condition jacobian matrix before factoring
   m%Jac(p%iJL, :) = m%Jac(p%iJL, :)/p%Scale_UJac
   m%Jac(:, p%iJL) = m%Jac(:, p%iJL)*p%Scale_UJac

   ! if (munit == -1) then
   !    call GetNewUnit(munit, ErrStat2, ErrMsg2); if (Failed()) return
   !    call DumpMatrix(munit, "dUdu.bin", m%dUdu, ErrStat2, ErrMsg2); if (Failed()) return
   !    call DumpMatrix(munit, "dUdy.bin", m%dUdy, ErrStat2, ErrMsg2); if (Failed()) return
   !    call DumpMatrix(munit, "dXdu.bin", m%dXdu, ErrStat2, ErrMsg2); if (Failed()) return
   !    call DumpMatrix(munit, "dXdx.bin", m%dXdx, ErrStat2, ErrMsg2); if (Failed()) return
   !    call DumpMatrix(munit, "dYdu.bin", m%dYdu, ErrStat2, ErrMsg2); if (Failed()) return
   !    call DumpMatrix(munit, "dYdx.bin", m%dYdx, ErrStat2, ErrMsg2); if (Failed()) return
   !    ! call DumpMatrix(munit, "ED-dXdu.bin", Turbine%ED%m%Vals%dXdu, ErrStat2, ErrMsg2); if (Failed()) return
   !    ! call DumpMatrix(munit, "ED-dXdx.bin", Turbine%ED%m%Vals%dXdx, ErrStat2, ErrMsg2); if (Failed()) return
   !    ! call DumpMatrix(munit, "ED-dYdu.bin", Turbine%ED%m%Vals%dYdu, ErrStat2, ErrMsg2); if (Failed()) return
   !    ! call DumpMatrix(munit, "ED-dYdx.bin", Turbine%ED%m%Vals%dYdx, ErrStat2, ErrMsg2); if (Failed()) return
   !    call DumpMatrix(munit, "G.bin", m%G, ErrStat2, ErrMsg2); if (Failed()) return
   ! end if
   ! call DumpMatrix(munit, "jacs/BD-dXdu-"//trim(num2lstr(iter))//".bin", Turbine%BD%m(1)%Vals%dXdu, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(munit, "jacs/BD-dXdx-"//trim(num2lstr(iter))//".bin", Turbine%BD%m(1)%Vals%dXdx, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(munit, "jacs/BD-dYdu-"//trim(num2lstr(iter))//".bin", Turbine%BD%m(1)%Vals%dYdu, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(munit, "jacs/BD-dYdx-"//trim(num2lstr(iter))//".bin", Turbine%BD%m(1)%Vals%dYdx, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(munit, "jacs/J-"//trim(num2lstr(iter))//".bin", m%Jac, ErrStat2, ErrMsg2); if (Failed()) return

   ! Factor jacobian matrix
   call LAPACK_getrf(size(m%Jac, 1), size(m%Jac, 1), m%Jac, m%IPIV, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine DumpMatrix(unit, filename, A, ErrStat, ErrMsg)
   integer(IntKi), intent(in)             :: unit
   character(*), intent(in)               :: filename
   real(R8Ki), intent(in)                 :: A(:, :)
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter                :: RoutineName = 'DumpMatrix'
   integer(IntKi)                         :: ErrStat2
   character(ErrMsgLen)                   :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ''

   call OpenBOutFile(unit, filename, ErrStat2, ErrMsg2)
   write (unit) int(shape(A), B4Ki)
   write (unit) pack(A, .true.)
   close (unit)
end subroutine

subroutine AddDeltaToStates(Mods, ModOrder, dx, x)
   type(ModDataType), intent(in)    :: Mods(:)  !< Solution variables from modules
   integer(IntKi), intent(in)       :: ModOrder(:)
   real(R8Ki), intent(in)           :: dx(:)
   real(R8Ki), intent(inout)        :: x(:)

   character(*), parameter          :: RoutineName = 'AddDeltaToStates'
   integer(IntKi)                   :: iMod, iIns
   integer(IntKi)                   :: i, j, k, ind(3)
   real(R8Ki)                       :: n(3), phi

   ! Loop through modules in order
   do i = 1, size(ModOrder)

      ! Loop through variables
      do j = 1, size(Mods(ModOrder(i))%Vars%x)

         associate (Var => Mods(ModOrder(i))%Vars%x(j))

            ! Select based on field type
            select case (Var%Field)
            case (VF_Force, VF_Moment, VF_TransDisp, VF_TransVel, VF_TransAcc, VF_AngularVel, VF_AngularAcc)
               ! Add delta x to x
               x(Var%iGblSol) = x(Var%iGblSol) + dx(Var%iGblSol)
            case (VF_AngularDisp)
               ! Add delta x to x and limit to between -2pi and 2pi
               ! x(ModData(i)%Vars%x(j)%iGblSol) = mod(x(ModData(i)%Vars%x(j)%iGblSol) + dx(ModData(i)%Vars%x(j)%iGblSol), TwoPi_R8)
               x(Var%iGblSol) = x(Var%iGblSol) + dx(Var%iGblSol)
            case (VF_Orientation)
               ! Compose WM components (dx is in radians)
               do k = 1, size(Var%iGblSol), 3
                  ind = Var%iGblSol(k:k + 2)
                  x(ind) = wm_compose(wm_from_xyz(dx(ind)), x(ind)) ! dx is in radians
               end do
            end select
         end associate
      end do
   end do

end subroutine

subroutine AddDeltaToInputs(Mods, ModOrder, du, u)
   type(ModDataType), intent(in)    :: Mods(:)  !< Solution variables from modules
   integer(IntKi), intent(in)       :: ModOrder(:)
   real(R8Ki), intent(in)           :: du(:)
   real(R8Ki), intent(inout)        :: u(:)

   character(*), parameter          :: RoutineName = 'AddDeltaToInputs'
   integer(IntKi)                   :: i, j, k, ind(3)

   ! Loop through modules in order
   do i = 1, size(ModOrder)

      ! Loop through variables
      do j = 1, size(Mods(ModOrder(i))%Vars%u)

         associate (Var => Mods(ModOrder(i))%Vars%u(j))

            ! If this is not a solve variable, cycle
            if (iand(Var%Flags, VF_Solve) == 0) cycle

            ! Select based on field type
            select case (Var%Field)
            case (VF_Force, VF_Moment, VF_TransDisp, VF_TransVel, VF_TransAcc, VF_AngularVel, VF_AngularAcc)
               ! Add delta u to u
               u(Var%iGblSol) = u(Var%iGblSol) + du(Var%iGblSol)
            case (VF_AngularDisp)
               ! Add delta u to u and limit to between -2pi and 2pi
               ! un(Var%iGblSol) = mod(u(Var%iGblSol) + du(Var%iGblSol), TwoPi_R8)
               u(Var%iGblSol) = u(Var%iGblSol) + du(Var%iGblSol)
            case (VF_Orientation)
               ! Compose WM components (du is in radians)
               do k = 1, size(Var%iGblSol), 3
                  ind = Var%iGblSol(k:k + 2)
                  u(ind) = wm_compose(wm_from_xyz(du(ind)), u(ind))  ! du is in radians
               end do
            end select

         end associate
      end do
   end do
end subroutine

subroutine PackModuleStates(ModData, this_state, T, x)
   type(ModDataType), intent(in)                      :: ModData(:)  !< Solution variables from modules
   integer(IntKi), intent(in)                         :: this_state  !< State index
   type(FAST_TurbineType), intent(inout)              :: T           !< Turbine type
   real(R8Ki), allocatable, optional, intent(inout)   :: x(:)
   integer(IntKi)                                     :: ii, j

   ! Must support all Tight Coupling modules
   if (present(x)) then
      do j = 1, size(ModData)
         ii = ModData(j)%Ins
         select case (ModData(j)%ID)
         case (Module_ED)
            call ED_PackStateValues(T%ED%p, T%ED%x(this_state), T%ED%m%Vals%x)
            call XferLocToGbl1D(ModData(j)%ixs, T%ED%m%Vals%x, x)
         case (Module_BD)
            call BD_PackStateValues(T%BD%p(ii), T%BD%x(ii, this_state), T%BD%m(ii)%Vals%x)
            call XferLocToGbl1D(ModData(j)%ixs, T%BD%m(ii)%Vals%x, x)
         case (Module_SD)
            ! call SD_PackStateValues(SD%p, SD%x(this_state), SD%m%Vals%x)
            ! x(SD%p%Vars%ix) = SD%m%Vals%x
         end select
      end do
   end if

end subroutine

subroutine UnpackModuleStates(ModData, ModOrder, this_state, T, x)
   type(ModDataType), intent(in)          :: ModData(:)  !< Solution variables from modules
   integer(IntKi), intent(in)             :: ModOrder(:)
   integer(IntKi), intent(in)             :: this_state  !< State index
   type(FAST_TurbineType), intent(inout)  :: T           !< Turbine type
   real(R8Ki), intent(inout)              :: x(:)
   integer(IntKi)                         :: j

   do j = 1, size(ModOrder)
      associate (Mod => ModData(ModOrder(j)))
         select case (Mod%ID)
         case (Module_ED)
            call ED_PackStateValues(T%ED%p, T%ED%x(this_state), T%ED%m%Vals%x)
            call XferGblToLoc1D(Mod%ixs, x, T%ED%m%Vals%x)
            call ED_UnpackStateValues(T%ED%p, T%ED%m%Vals%x, T%ED%x(this_state))
         case (Module_BD)
            call BD_PackStateValues(T%BD%p(Mod%Ins), T%BD%x(Mod%Ins, this_state), T%BD%m(Mod%Ins)%Vals%x)
            call XferGblToLoc1D(Mod%ixs, x, T%BD%m(Mod%Ins)%Vals%x)
            call BD_UnpackStateValues(T%BD%p(Mod%Ins), T%BD%m(Mod%Ins)%Vals%x, T%BD%x(Mod%Ins, this_state))
         case (Module_SD)
            ! call SD_UnpackStateValues(SD%p, x(SD%p%Vars%ix), SD%x(this_state))
         end select
      end associate
   end do
end subroutine

! PackModuleInputs packs input values from Option 1 modules
subroutine PackModuleInputs(ModData, ModOrder, T, u, u_tmp)
   type(ModDataType), intent(in)          :: ModData(:)  !< Solution variables from modules
   integer(IntKi), intent(in)             :: ModOrder(:)
   type(FAST_TurbineType), intent(inout)  :: T           !< Turbine type
   real(R8Ki), optional, intent(inout)    :: u(:), u_tmp(:)
   integer(IntKi)                         :: j

   if (present(u)) then
      do j = 1, size(ModOrder)
         associate (Mod => ModData(ModOrder(j)))
            select case (Mod%ID)
            case (Module_ED)
               call ED_PackInputValues(T%ED%p, T%ED%Input(1), T%ED%m%Vals%u)
               call XferLocToGbl1D(Mod%ius, T%ED%m%Vals%u, u)
            case (Module_BD)
               call BD_PackInputValues(T%BD%p(Mod%Ins), T%BD%Input(1, Mod%Ins), T%BD%m(Mod%Ins)%Vals%u)
               call XferLocToGbl1D(Mod%ius, T%BD%m(Mod%Ins)%Vals%u, u)
            case (Module_SD)
               ! call SD_PackInputValues(SD%p, SD%Input, SD%m%Vals%u)
               ! u(SD%p%Vars%iu) = SD%m%Vals%u
            end select
         end associate
      end do
   end if

   if (present(u_tmp)) then
      do j = 1, size(ModOrder)
         associate (Mod => ModData(ModOrder(j)))
            select case (Mod%ID)
            case (Module_ED)
               call ED_PackInputValues(T%ED%p, T%ED%u, T%ED%m%Vals%u)
               call XferLocToGbl1D(Mod%ius, T%ED%m%Vals%u, u_tmp)
            case (Module_BD)
               call BD_PackInputValues(T%BD%p(Mod%Ins), T%BD%u(Mod%Ins), T%BD%m(Mod%Ins)%Vals%u)
               call XferLocToGbl1D(Mod%ius, T%BD%m(Mod%Ins)%Vals%u, u_tmp)
            case (Module_SD)
               ! call SD_PackInputValues(SD%p, SD%Input, SD%m%Vals%u)
               ! u(SD%p%Vars%iu) = SD%m%Vals%u
            end select
         end associate
      end do
   end if
end subroutine

! PackModuleInputs packs input values from Option 1 modules
subroutine UnpackModuleInputs(ModData, ModOrder, T, u)
   type(ModDataType), intent(in)          :: ModData(:)  !< Solution variables from modules
   integer(IntKi), intent(in)             :: ModOrder(:)
   type(FAST_TurbineType), intent(inout)  :: T           !< Turbine type
   real(R8Ki), intent(in)                 :: u(:)
   integer(IntKi)                         :: j

   do j = 1, size(ModOrder)
      associate (Mod => ModData(ModOrder(j)))
         select case (Mod%ID)
         case (Module_ED)
            call ED_PackInputValues(T%ED%p, T%ED%Input(1), T%ED%m%Vals%u)
            call XferGblToLoc1D(Mod%ius, u, T%ED%m%Vals%u)
            call ED_UnpackInputValues(T%ED%p, T%ED%m%Vals%u, T%ED%Input(1))
         case (Module_BD)
            call BD_PackInputValues(T%BD%p(Mod%Ins), T%BD%Input(1, Mod%Ins), T%BD%m(Mod%Ins)%Vals%u)
            call XferGblToLoc1D(Mod%ius, u, T%BD%m(Mod%Ins)%Vals%u)
            call BD_UnpackInputValues(T%BD%p(Mod%Ins), T%BD%m(Mod%Ins)%Vals%u, T%BD%Input(1, Mod%Ins))
         case (Module_SD)
            ! call SD_UnpackInputValues(SD%p, u(SD%p%Vars%iu), SD%Input(1))
         end select
      end associate
   end do

end subroutine

! PackModuleOutputs packs output values from Option 1 modules
subroutine PackModuleOutputs(ModData, ModOrder, T, y)
   type(ModDataType), intent(in)          :: ModData(:)  !< Solution variables from modules
   integer(IntKi), intent(in)             :: ModOrder(:)
   type(FAST_TurbineType), intent(inout)  :: T           !< Turbine type
   real(R8Ki), intent(inout)              :: y(:)
   integer(IntKi)                         :: j

   do j = 1, size(ModOrder)
      associate (Mod => ModData(ModOrder(j)))
         select case (Mod%ID)
         case (Module_ED)
            call ED_PackOutputValues(T%ED%p, T%ED%y, T%ED%m%Vals%y)
            call XferLocToGbl1D(Mod%iys, T%ED%m%Vals%y, y)
         case (Module_BD)
            call BD_PackOutputValues(T%BD%p(Mod%Ins), T%BD%y(Mod%Ins), T%BD%m(Mod%Ins)%Vals%y)
            call XferLocToGbl1D(Mod%iys, T%BD%m(Mod%Ins)%Vals%y, y)
         end select
      end associate
   end do

end subroutine

end module
