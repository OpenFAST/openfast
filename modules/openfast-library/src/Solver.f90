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

contains

subroutine Solver_Init(p, m, ModData, ErrStat, ErrMsg)

   type(TC_ParameterType), intent(inout)     :: p           !< Parameters
   type(TC_MiscVarType), intent(out)         :: m           !< Misc variables for optimization (not copied in glue code)
   type(ModDataType), intent(inout)          :: ModData(:)  !< Solution variables from modules
   integer(IntKi), intent(out)               :: ErrStat     !< Error status of the operation
   character(*), intent(out)                 :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   character(*), parameter                   :: RoutineName = 'Solver_Init'
   integer(IntKi)                            :: ErrStat2    ! local error status
   character(ErrMsgLen)                      :: ErrMsg2     ! local error message
   integer(IntKi)                            :: i, j, k, n
   integer(IntKi)                            :: NumX, NumU, NumY, NumQ, NumJac
   integer(IntKi), allocatable               :: iq(:), modIDs(:)
   logical                                   :: isLoad
   logical, allocatable                      :: isLoadTight(:), isLoadOption1(:)

   !----------------------------------------------------------------------------
   ! Calculate generalized alpha parameters
   !----------------------------------------------------------------------------

   p%AccBlend = 1.0_R8Ki

   p%AlphaM = (2.0_R8Ki*p%RhoInf - 1.0_R8Ki)/(p%RhoInf + 1.0_R8Ki)
   p%AlphaF = p%RhoInf/(p%RhoInf + 1.0_R8Ki)
   p%Gamma = 0.5_ReKi - p%AlphaM + p%AlphaF
   p%Beta = (1.0_R8Ki - p%AlphaM + p%AlphaF)**2/4.0_R8Ki

   p%C(1) = (1.0_R8Ki - p%AlphaF)/(1.0_R8Ki - p%AlphaM)
   p%C(2) = p%DT*p%Gamma*p%C(1)
   p%C(3) = p%DT**2*p%Beta*p%C(1)

   !----------------------------------------------------------------------------
   ! Module substepping
   !----------------------------------------------------------------------------

   ! Loop through module data
   do i = 1, size(ModData)

      ! Calculate the number of substeps
      ModData(i)%SubSteps = nint(p%DT/ModData(i)%DT)

      ! If module time step is same as global time step, set substeps to 1
      if (EqualRealNos(ModData(i)%DT, p%DT)) cycle

      ! If the module time step is greater than the global time step, set error
      if (ModData(i)%DT > p%DT) then
         call SetErrStat(ErrID_Fatal, "The "//trim(ModData(i)%Abbr)// &
                         " module time step ("//trim(Num2LStr(ModData(i)%DT))//" s) "// &
                         "cannot be larger than FAST time step ("//trim(Num2LStr(p%DT))//" s).", &
                         ErrStat, ErrMsg, RoutineName)
         return
      end if

      ! If the module DT is not an exact integer divisor of the global time step, set error
      if (.not. EqualRealNos(p%DT, ModData(i)%DT*ModData(i)%SubSteps)) then
         call SetErrStat(ErrID_Fatal, "The "//trim(ModData(i)%Abbr)// &
                         " module time step ("//trim(Num2LStr(ModData(i)%DT))//" s) "// &
                         "must be an integer divisor of the FAST time step ("//trim(Num2LStr(p%DT))//" s).", &
                         ErrStat, ErrMsg, RoutineName)
         return
      end if
   end do

   !----------------------------------------------------------------------------
   ! Calculate Variable cateogries and global indices
   ! TODO: reorder to improve data locality
   !----------------------------------------------------------------------------

   NumX = 0
   do i = 1, size(ModData)
      do j = 1, size(ModData(i)%Vars%x)
         ModData(i)%Vars%x(j)%Cat = VarCategory(ModData(i)%ID, ModData(i)%Vars%x(j)%Field)
         ModData(i)%Vars%x(j)%iGbl = [(NumX + k, k=1, ModData(i)%Vars%x(j)%Size)]
         NumX = NumX + ModData(i)%Vars%x(j)%Size
      end do
      call MV_CollectGlobalIndices(ModData(i)%Vars%x, ModData(i)%Vars%ixg)
   end do

   NumU = 0
   do i = 1, size(ModData)
      do j = 1, size(ModData(i)%Vars%u)
         ModData(i)%Vars%u(j)%Cat = VarCategory(ModData(i)%ID, ModData(i)%Vars%u(j)%Field)
         ModData(i)%Vars%u(j)%iGbl = [(NumU + k, k=1, ModData(i)%Vars%u(j)%Size)]
         NumU = NumU + ModData(i)%Vars%u(j)%Size
      end do
      call MV_CollectGlobalIndices(ModData(i)%Vars%u, ModData(i)%Vars%iug)
   end do

   NumY = 0
   do i = 1, size(ModData)
      do j = 1, size(ModData(i)%Vars%y)
         ModData(i)%Vars%y(j)%Cat = VarCategory(ModData(i)%ID, ModData(i)%Vars%y(j)%Field)
         ModData(i)%Vars%y(j)%iGbl = [(NumY + k, k=1, ModData(i)%Vars%y(j)%Size)]
         NumY = NumY + ModData(i)%Vars%y(j)%Size
      end do
      call MV_CollectGlobalIndices(ModData(i)%Vars%y, ModData(i)%Vars%iyg)
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
   m%u = 0.0_R8Ki
   m%un = 0.0_R8Ki
   m%du = 0.0_R8Ki
   m%u_tmp = 0.0_R8Ki

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
   do i = 1, size(ModData)

      ! Allocate iq to store q index for each variable
      call AllocAry(iq, size(ModData(i)%Vars%x), "iq", ErrStat2, ErrMsg2); if (Failed()) return
      iq = 0

      ! Skip modules that aren't in tight coupling
      if (all(ModData(i)%ID /= TC_Modules)) cycle

      ! Loop through state variables
      do j = 1, size(ModData(i)%Vars%x)

         ! Skip variables which already have a q index
         if (iq(j) /= 0) cycle

         ! Skip load variables (force and moment)
         if (any(ModData(i)%Vars%x(j)%Field == LoadFields)) cycle

         ! Set q index for variable and update number
         iq(j) = NumQ + 1
         NumQ = NumQ + ModData(i)%Vars%x(j)%Size

         ! Loop through remaining vars  if the names match
         do k = i + 1, size(ModData(i)%Vars%x)

            ! If names are different then they don't match, skip
            if (ModData(i)%Vars%x(j)%Name /= ModData(i)%Vars%x(k)%Name) cycle

            ! If field is not the same or a derivative of current field, skip
            select case (ModData(i)%Vars%x(j)%Field)
            case (VF_Orientation)
               if (ModData(i)%Vars%x(k)%Field /= VF_Orientation) cycle
            case (VF_TransDisp, VF_TransVel, VF_TransAcc)
               if (all(ModData(i)%Vars%x(k)%Field /= TransFields)) cycle
            case (VF_AngularDisp, VF_AngularVel, VF_AngularAcc)
               if (all(ModData(i)%Vars%x(k)%Field /= AngularFields)) cycle
            case (VF_Force, VF_Moment)
               cycle
            end select

            ! Copy q index
            iq(k) = iq(j)

         end do
      end do

      ! Loop through state variables and build mapping array between x and q
      ! ixqd is 3xN where each row is [global x array index, q matrix row, q matrix col]
      do j = 1, size(ModData(i)%Vars%x)
         do k = 1, ModData(i)%Vars%x(j)%Size
            n = n + 1
            p%ixqd(:, n) = [ModData(i)%Vars%x(j)%iGbl(k), iq(j) + k - 1, ModData(i)%Vars%x(j)%DerivOrder + 1]
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
   do i = 1, size(ModData)
      do j = 1, size(ModData(i)%Vars%x)
         if (ModData(i)%Vars%x(j)%Cat == VC_Tight) then
            select case (ModData(i)%Vars%x(j)%DerivOrder)
            case (0)
               p%iX1Tight = [p%iX1Tight, ModData(i)%Vars%x(j)%iGbl]
            case (1)
               p%iX2Tight = [p%iX2Tight, ModData(i)%Vars%x(j)%iGbl]
            end select
         end if
      end do
   end do

   ! Calculate index mapping arrays for U^Tight and U^Option1
   call AllocAry(p%iUTight, 0, "p%iUTight", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(p%iUOpt1, 0, "p%iUOpt1", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(isLoadTight, 0, "isLoadTight", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(isLoadOption1, 0, "isLoadOption1", ErrStat2, ErrMsg2); if (Failed()) return
   ! do i = 1, size(ModData)
   !    do j = 1, size(ModData(i)%Vars%u)
   !       isLoad = any(LoadFields == ModData(i)%Vars%u(j)%Field)
   !       select case (ModData(i)%Vars%u(j)%Cat)
   !       case (VC_Tight)
   !          p%iUTight = [p%iUTight, ModData(i)%Vars%u(j)%iGbl]
   !          isLoadTight = [isLoadTight, spread(isLoad, 1, ModData(i)%Vars%u(j)%Size)]
   !       case (VC_Option1)
   !          p%iUOpt1 = [p%iUOpt1, ModData(i)%Vars%u(j)%iGbl]
   !          isLoadOption1 = [isLoadOption1, spread(isLoad, 1, ModData(i)%Vars%u(j)%Size)]
   !       end select
   !    end do
   ! end do

   ! Calculate index mapping arrays for y^Tight and y^Option1
   call AllocAry(p%iyTight, 0, "p%iyTight", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(p%iyOpt1, 0, "p%iyOpt1", ErrStat2, ErrMsg2); if (Failed()) return
   ! do i = 1, size(ModData)
   !    do j = 1, size(ModData(i)%Vars%u)
   !       select case (ModData(i)%Vars%u(j)%Cat)
   !       case (VC_Tight)
   !          p%iyTight = [p%iyTight, ModData(i)%Vars%u(j)%iGbl]
   !       case (VC_Option1)
   !          p%iyOpt1 = [p%iyOpt1, ModData(i)%Vars%u(j)%iGbl]
   !       end select
   !    end do
   ! end do

   !----------------------------------------------------------------------------
   ! Allocate Jacobian matrices
   !----------------------------------------------------------------------------

   call AllocAry(m%dYdx, NumY, NumX, "m%dYdx", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%dYdu, NumY, NumU, "m%dYdu", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%dXdx, NumX, NumX, "m%dXdx", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%dXdu, NumX, NumU, "m%dXdu", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%dUdu, NumU, NumU, "m%dUdu", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%dUdy, NumU, NumY, "m%dUdy", ErrStat2, ErrMsg2); if (Failed()) return

   m%dYdx = 0.0_R8Ki
   m%dYdu = 0.0_R8Ki
   m%dXdx = 0.0_R8Ki
   m%dXdu = 0.0_R8Ki
   m%dUdu = 0.0_R8Ki
   m%dUdy = 0.0_R8Ki

   ! Populate dUdu
   do i = 1, NumU
      m%dUdu(i, i) = 1.0_R8Ki
   end do

   ! TODO: Figure out how to calculate from mesh mapping
   ! Populate dUdy based on output/input transfer relationships
   ! do i = 1, size(p%Vars%u)
   !    do j = 1, size(p%Vars%u(i)%iy)
   !       do k = 1, size(p%Vars%u(i)%iGbl)
   !          m%dUdy(p%Vars%u(i)%iGbl(k), p%Vars%u(i)%iy(k, j)) = -1.0_R8Ki
   !       end do
   !    end do
   ! end do

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

   !----------------------------------------------------------------------------
   ! Module ordering for solve
   !----------------------------------------------------------------------------

   ! Get array of module IDs
   modIDs = [(ModData(i)%ID, i=1, size(ModData))]

   ! Allocate ordering array for all modules
   p%iModAll = [pack([(i, i=1, size(ModData))], ModIDs == Module_SrvD), &
                pack([(i, i=1, size(ModData))], ModIDs == Module_AD), &
                pack([(i, i=1, size(ModData))], ModIDs == Module_ED), &
                pack([(i, i=1, size(ModData))], ModIDs == Module_BD), &
                pack([(i, i=1, size(ModData))], ModIDs == Module_SD)]

   p%iModTC = [pack([(i, i=1, size(ModData))], ModIDs == Module_ED), &
               pack([(i, i=1, size(ModData))], ModIDs == Module_BD), &
               pack([(i, i=1, size(ModData))], ModIDs == Module_SD)]

   ! TODO: Add other modules
   p%iModOpt2 = [pack([(i, i=1, size(ModData))], ModIDs == Module_ED)]

   ! TODO: Add other modules
   p%iModOpt1 = [pack([(i, i=1, size(ModData))], ModIDs == Module_ED)]

   ! Option 1 modules that require input solve and update states
   ! TODO: Add aerodyn with hydrodynamics
   p%iModOpt1US = [pack([(i, i=1, size(ModData))], ModIDs == Module_ExtPtfm), &
                   pack([(i, i=1, size(ModData))], ModIDs == Module_HD), &
                   pack([(i, i=1, size(ModData))], ModIDs == Module_Orca)]

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

   ! TODO: Add SrvD_SetExternalInputs
   ! TODO: Add SeaSt_CalcOutput

   !----------------------------------------------------------------------------
   ! Calculate initial accelerations
   !----------------------------------------------------------------------------

   ! Transfer initial state from modules to solver
   call TransferFromModules(p%iModTC, ModData, STATE_CURR, Turbine, x=m%x)

   ! Transfer initial state to next state q matrix (qn)
   ! The qn matrix is being used because the solution step routine predicts
   ! the the step starting state, q, from qn.
   call Solver_TransferXtoQ(p%ixqd, m%x, m%qn)

   ! Allocate acceleration array which will be used to check for convergence
   ! of initial acceleration. Transfer initial accelerations from q matrix
   call AllocAry(accel, size(m%qn, dim=2), "accel", ErrStat2, ErrMsg2); if (Failed()) return
   accel = m%qn(:, COL_A)

   ! Loop until initial accelerations are converged, or max iterations are reached.
   ! TODO: may need a separate variable for max initial acceleration convergence iterations
   converged = .false.
   k = 1
   do while ((.not. converged) .and. (k <= p%MaxConvIter))

      ! Transfer inputs and calculate outputs for all modules (use current state)
      call FAST_EvalModules(t_initial, n_t_global, p%iModAll, ModData, STATE_CURR, &
                            EM_InputSolve + EM_CalcOutput, Turbine, ErrStat2, ErrMsg2); if (Failed()) return

      ! Calculate continuous state derivatives for tight coupling modules (use current state)
      call FAST_EvalModules(t_initial, n_t_global, p%iModTC, ModData, STATE_CURR, &
                            EM_CalcContStateDeriv, Turbine, ErrStat2, ErrMsg2, dxdt=m%dxdt); if (Failed()) return

      ! Copy acceleration (derivative of velocity) to acceleration array.
      ! Loop through x->q index array and transfer derivative of
      ! velocity variables (deriv=COL_V) to the acceleration column
      do i = 1, size(p%ixqd, dim=2)
         if (p%ixqd(3, i) == COL_V) accel(p%ixqd(2, i)) = m%dxdt(p%ixqd(1, i))
      end do

      ! Calculate diff as L2 norm of current and new accelerations
      diff = TwoNorm(accel - m%qn(:, COL_A))

      ! If difference is less than converence tolerance, set flag and exit loop
      if (diff < p%ConvTol) converged = .true.

      ! Update acceleration in q matrix
      m%qn(:, COL_A) = accel

      ! Increment iteration counter
      k = k + 1
   end do

   ! Print warning if not converged
   if (.not. converged) then
      call WrScr("Solver: initial accel not converged, diff="//Num2LStr(diff)//", tol="//Num2LStr(p%ConvTol))
   end if

   ! Initialize algorithmic acceleration from actual acceleration
   m%qn(:, COL_AA) = m%qn(:, COL_A)
   ! m%qn(:, COL_A) = 0.0_R8Ki

   !----------------------------------------------------------------------------
   ! Initialize module input and state arrays for interpolation/extrapolation
   !----------------------------------------------------------------------------

   call FAST_EvalModules(t_initial, n_t_global, p%iModAll, ModData, STATE_CURR, &
                         EM_InitIO, Turbine, ErrStat2, ErrMsg2); if (Failed()) return

   !----------------------------------------------------------------------------
   ! Calculate inititial Jacobian
   !----------------------------------------------------------------------------

   ! Calculate the Jacobians for TC modules
   call FAST_EvalModules(t_initial, n_t_global, p%iModTC, ModData, STATE_CURR, &
                         EM_JacobianPContState + EM_JacobianPInput, Turbine, ErrStat2, ErrMsg2, &
                         dYdx=m%dYdx, dXdx=m%dXdx, dYdu=m%dYdu, dXdu=m%dXdu); if (Failed()) return

   ! Calculate the Jacobian
   call CalcJacobian(p, m, ErrStat2, ErrMsg2); if (Failed()) return

   !----------------------------------------------------------------------------
   ! Supercontroller
   !----------------------------------------------------------------------------

   ! TODO: add SC_DX_SetInputs

   !----------------------------------------------------------------------------
   ! Write Output
   !----------------------------------------------------------------------------

   call WriteOutputToFile(n_t_global_next, t_initial, Turbine, ErrStat2, ErrMsg2); if (Failed()) return

contains
   function Failed()
      logical :: Failed
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine Solver_Step(n_t_global, t_initial, p, m, ModData, Turbine, ErrStat, ErrMsg)
   integer(IntKi), intent(in)                :: n_t_global  !< global time step
   real(DbKi), intent(in)                    :: t_initial   !< Initial simulation time
   type(TC_ParameterType), intent(in)        :: p           !< Parameters
   type(TC_MiscVarType), intent(inout)       :: m           !< Misc variables
   type(ModDataType), intent(in)             :: ModData(:)  !< Solution variables from modules
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

   ! Set number of corrections to be used for this time step
   if (Turbine%p_FAST%CompElast /= Module_BD) then
      ! Use input file value if BeamDyn is not used
      NumCorrections = Turbine%p_FAST%NumCrctn
   else if (n_t_global == 0) then
      ! BD accelerations have fewer spikes on the first several time steps with these corrections
      NumCorrections = max(Turbine%p_FAST%NumCrctn, 16)
   else if (n_t_global <= 2) then
      ! Use at least 1 correction on the first 3 steps (the 2 should probably be related to Turbine%p_FAST%InterpOrder)
      NumCorrections = max(Turbine%p_FAST%NumCrctn, 1)
   else
      ! Use input file value on subsequent steps
      NumCorrections = Turbine%p_FAST%NumCrctn
   end if

   ! Decrement number of time steps before updating the Jacobian
   m%StepsUntilUJac = m%StepsUntilUJac - 1

   !----------------------------------------------------------------------------
   ! Extrapolate/interpolate inputs for all modules
   !----------------------------------------------------------------------------

   call FAST_EvalModules(t_initial, n_t_global, p%iModAll, ModData, STATE_PRED, &
                         EM_ExtrapInterp, Turbine, ErrStat2, ErrMsg2)

   !----------------------------------------------------------------------------
   ! Prediction - guess solution state variables at end of time step
   !----------------------------------------------------------------------------

   ! Calculate displacements, velocities, and accelerations for next step
   associate (q => m%qn(:, 1), qd => m%qn(:, 2), qdd => m%qn(:, 3), aa => m%qn(:, 4), &
              nq => m%q(:, 1), nqd => m%q(:, 2), nqdd => m%q(:, 3), naa => m%q(:, 4))
      nqdd = qdd*p%AccBlend                                                         ! Acceleration
      naa = ((1.0_R8Ki - p%AlphaF)*nqdd + p%AlphaF*qdd - p%AlphaM*aa)/(1 - p%AlphaM)  ! Algorithmic acceleration
      nq = q + p%DT*qd + p%DT**2*(0.5 - p%Beta)*aa + p%DT**2*p%Beta*naa             ! Displacment
      nqd = qd + p%DT*(1 - p%Gamma)*aa + p%DT*p%Gamma*naa                           ! Velocity
   end associate

   ! Transfer algorithm matrix to solver states
   call Solver_TransferQtoX(p%ixqd, m%q, m%x)

   !----------------------------------------------------------------------------
   ! Correction Iterations
   !----------------------------------------------------------------------------

   ! Loop through correction iterations
   iterCorr = 0
   do while (iterCorr <= p%NumCrctn)

      ! Copy state for correction step
      m%qn = m%q
      m%xn = m%x

      ! Reset state matrix delta
      m%dq = 0.0_R8Ki

      ! Option 2 calculation:
      ! - Input Solve
      ! - Update states from current to predicted (set new states in TC modules)
      ! - Calculate output
      call FAST_EvalModules(t_initial, n_t_global, p%iModOpt2, ModData, STATE_PRED, &
                            EM_UpdateStates + EM_InputSolve + EM_CalcOutput, &
                            Turbine, ErrStat2, ErrMsg2, x=m%xn); if (Failed()) return

      ! Option 1 modules not in Option 2:
      ! - Input solve
      ! - Update states
      call FAST_EvalModules(t_initial, n_t_global, p%iModOpt1US, ModData, STATE_PRED, &
                            EM_InputSolve + EM_UpdateStates, &
                            Turbine, ErrStat2, ErrMsg2); if (Failed()) return

      ! Get all current inputs
      call TransferFromModules(p%iModAll, ModData, STATE_PRED, Turbine, u=m%un)

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

         call FAST_EvalModules(t_initial, n_t_global, p%iModOpt1, ModData, STATE_PRED, &
                               EM_CalcOutput, &
                               Turbine, ErrStat2, ErrMsg2); if (Failed()) return

         !----------------------------------------------------------------------
         ! If iteration limit reached, exit loop
         !----------------------------------------------------------------------

         if (iterConv >= p%MaxConvIter) exit

         !----------------------------------------------------------------------
         ! Update Jacobian
         !----------------------------------------------------------------------

         ! If number of iterations or steps until Jacobian is to be updated
         ! is zero or less, then update the Jacobian. Note: CalcJacobian
         ! resets these counters.
         if ((m%IterUntilUJac <= 0) .or. (m%StepsUntilUJac <= 0)) then

            ! Calculate Jacobians for TC modules
            call FAST_EvalModules(t_initial, n_t_global, p%iModTC, ModData, STATE_PRED, &
                                  EM_JacobianPContState + EM_JacobianPInput, &
                                  Turbine, ErrStat2, ErrMsg2, &
                                  dYdx=m%dYdx, dXdx=m%dXdx, dYdu=m%dYdu, dXdu=m%dXdu); if (Failed()) return

            ! Calculate the Jacobian
            call CalcJacobian(p, m, ErrStat2, ErrMsg2); if (Failed()) return
         end if

         !----------------------------------------------------------------------
         ! Formulate right hand side (X_2^tight, U^tight, U^Option1)
         !----------------------------------------------------------------------

         ! Calculate continuous state derivatives for tight coupling modules
         call FAST_EvalModules(t_initial, n_t_global, p%iModTC, ModData, STATE_PRED, &
                               EM_CalcContStateDeriv, Turbine, ErrStat2, ErrMsg2, &
                               dxdt=m%dxdt); if (Failed()) return

         ! Option 1 transfer outputs to inputs
         call FAST_EvalModules(t_initial, n_t_global, p%iModOpt1, ModData, STATE_PRED, &
                               EM_InputSolve, Turbine, ErrStat2, ErrMsg2); if (Failed()) return

         ! Get all updated inputs and store in temporary input variable
         call TransferFromModules(p%iModOpt1, ModData, STATE_PRED, Turbine, u=m%u_tmp)

         ! Calculate difference between predicted and actual acceleration and add to RHS
         m%XB(p%iJX2, 1) = m%qn(:, COL_A) - m%dxdt(p%iX2Tight)

         ! Calculate difference in input/output and add to RHS
         m%XB(p%iJT, 1) = m%un(p%iUTight) - m%u_tmp(p%iUTight)

         ! If Option 1 inputs, calc difference and add to RHS
         if (size(p%iJ1) > 0) m%XB(p%iJ1, 1) = m%un(p%iUOpt1) - m%u_tmp(p%iUOpt1)

         ! Apply conditioning factor to loads in RHS
         m%XB(p%iJL, 1) = m%XB(p%iJL, 1)/p%Scale_UJac

         !----------------------------------------------------------------------
         ! Solve for state and input perturbations
         !----------------------------------------------------------------------

         ! Solve Jacobian and RHS
         call LAPACK_getrs('N', size(m%Jac, 1), m%Jac, m%IPIV, m%XB, ErrStat2, ErrMsg2); if (Failed()) return

         ! Remove conditioning
         m%XB(p%iJL, 1) = m%XB(p%iJL, 1)*p%Scale_UJac

         !----------------------------------------------------------------------
         ! Check perturbations for convergence and exit if below tolerance
         !----------------------------------------------------------------------

         delta_norm = TwoNorm(m%XB(:, 1))
         if (delta_norm < p%ConvTol) exit

         !----------------------------------------------------------------------
         ! Update State for Tight Coupling modules
         !----------------------------------------------------------------------

         ! Update delta state matrix
         m%dq(:, COL_D) = m%dq(:, COL_D) - p%C(3)*m%XB(p%iJX2, 1)
         m%dq(:, COL_V) = m%dq(:, COL_V) - p%C(2)*m%XB(p%iJX2, 1)
         m%dq(:, COL_A) = m%dq(:, COL_A) - m%XB(p%iJX2, 1)
         m%dq(:, COL_AA) = m%dq(:, COL_AA) - p%C(1)*m%XB(p%iJX2, 1)

         ! Transfer change in q state matrix to change in x array
         call Solver_TransferQtoX(p%ixqd, m%dq, m%dx)

         ! Add delta to x array to get new states (respect variable fields)
         call AddDeltaToStates(ModData, p%iModTC, m%x, m%dx, m%xn)

         ! Calculate new state matrix as sum of previous state and deltas
         m%qn = m%q + m%dq

         ! Update new state matrix with new state array values
         call Solver_TransferXtoQ(p%ixqd, m%xn, m%qn)

         ! Transfer new state to modules
         call TransferToModules(p%iModTC, ModData, STATE_PRED, Turbine, x=m%xn)

         !----------------------------------------------------------------------
         ! Update inputs for Tight Coupling and Option 1 modules
         !----------------------------------------------------------------------

         ! Update change in inputs
         m%du(p%iUTight) = m%du(p%iUTight) - m%XB(p%iJT, 1)
         m%du(p%iUOpt1) = m%du(p%iUOpt1) - m%XB(p%iJ1, 1)

         ! Apply deltas to inputs, update modules
         call AddDeltaToInputs(ModData, p%iModOpt1, m%u, m%du, m%un)
         call TransferToModules(p%iModOpt1, ModData, STATE_PRED, Turbine, u=m%un)

      end do

      m%u = m%un
      iterCorr = iterCorr + 1

   end do

   !----------------------------------------------------------------------------
   ! Update for next step
   !----------------------------------------------------------------------------

   ! Copy the final predicted states from step t_global_next to actual states for that step
   call FAST_EvalModules(t_initial, n_t_global, p%iModAll, ModData, STATE_PRED, &
                         EM_SavePredStates, Turbine, ErrStat2, ErrMsg2); if (Failed()) return

   ! Save new state
   m%x = m%xn

   ! Update the global time
   Turbine%m_FAST%t_global = t_global_next

   !----------------------------------------------------------------------------
   ! Write output to file
   !----------------------------------------------------------------------------

   call WriteOutputToFile(n_t_global_next, t_global_next, Turbine, ErrStat2, ErrMsg2); if (Failed()) return

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
   function Failed()
      logical :: Failed
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine CalcJacobian(p, m, ErrStat, ErrMsg)

   type(TC_ParameterType), intent(in)  :: p           !< Parameters
   type(TC_MiscVarType), intent(inOUT) :: m           !< Misc variables
   integer(IntKi), intent(out)         :: ErrStat
   character(*), intent(out)           :: ErrMsg

   character(*), parameter             :: RoutineName = 'CalcJacobian'
   integer(IntKi)                      :: ErrStat2
   character(ErrMsgLen)                :: ErrMsg2
   integer(IntKi) :: i

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Group (1,1)
   m%Jac(p%iJX2, p%iJX2) = -p%C(2)*m%dXdx(p%iX2Tight, p%iX2Tight) - &
                           p%C(3)*m%dXdx(p%iX2Tight, p%iX1Tight)
   do i = 1, size(p%iJX2)
      m%Jac(i, i) = m%Jac(i, i) + 1.0_R8Ki
   end do

   ! Group (2,1)
   m%Jac(p%iJT, p%iJX2) = p%C(2)*matmul(m%dUdy(p%iUTight, p%iyTight), &
                                        m%dYdx(p%iyTight, p%iX2Tight)) + &
                          p%C(3)*matmul(m%dUdy(p%iUTight, p%iyTight), &
                                        m%dYdx(p%iyTight, p%iX1Tight))

   ! Group (1,2)
   m%Jac(p%iJX2, p%iJT) = -m%dXdu(p%iX2Tight, p%iUTight)

   ! Group (2,2)
   m%Jac(p%iJT, p%iJT) = m%dUdu(p%iUTight, p%iUTight) + &
                         matmul(m%dUdy(p%iUTight, p%iyTight), &
                                m%dYdu(p%iyTight, p%iUTight))

   ! If modules in option 1
   if (size(p%iJ1) > 0) then

      ! Group (3,2)
      m%Jac(p%iJ1, p%iJT) = m%dUdu(p%iUOpt1, p%iUTight) + &
                            matmul(m%dUdy(p%iUOpt1, p%iyTight), &
                                   m%dYdu(p%iyTight, p%iUTight))
      ! Group (2,3)
      m%Jac(p%iJT, p%iJ1) = m%dUdu(p%iUTight, p%iUOpt1) + &
                            matmul(m%dUdy(p%iUTight, p%iyOpt1), &
                                   m%dYdu(p%iyOpt1, p%iUOpt1))
      ! Group (3,3)
      m%Jac(p%iJ1, p%iJ1) = m%dUdu(p%iUOpt1, p%iUOpt1) + &
                            matmul(m%dUdy(p%iUOpt1, p%iyOpt1), &
                                   m%dYdu(p%iyOpt1, p%iUOpt1))
   end if

   ! Condition jacobian matrix before factoring
   m%Jac(p%iJL, :) = m%Jac(p%iJL, :)/p%Scale_UJac
   m%Jac(:, p%iJL) = m%Jac(:, p%iJL)*p%Scale_UJac

   ! Factor jacobian matrix
   call LAPACK_getrf(size(m%Jac, 1), size(m%Jac, 1), m%Jac, m%IPIV, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ! Reset Jacobian update countdown values
   m%IterUntilUJac = p%NIter_UJac
   m%StepsUntilUJac = p%NStep_UJac

end subroutine

subroutine AddDeltaToStates(ModData, ModOrder, x, dx, xn)
   use ModVar
   type(ModDataType), intent(in)   :: ModData(:)  !< Solution variables from modules
   integer(IntKi), intent(in)      :: ModOrder(:) !< Array of module indices to evaluate
   real(R8Ki), intent(in)          :: dx(:)
   real(R8Ki), intent(in)          :: x(:)
   real(R8Ki), intent(out)         :: xn(:)
   integer(IntKi)                  :: iMod, iIns
   integer(IntKi)                  :: i, j, k, iGbl(3)

   ! Loop through modules in order
   do i = 1, size(ModOrder)

      iMod = ModOrder(i)
      iIns = ModData(iMod)%Instance

      ! Loop through variables
      do j = 1, size(ModData(iMod)%Vars%x)

         ! Select based on field type
         select case (ModData(iMod)%Vars%x(j)%Field)
         case (VF_TransDisp, VF_TransVel, VF_TransAcc, VF_AngularVel, VF_AngularAcc)
            ! Add delta x to x
            xn(ModData(iMod)%Vars%x(j)%iGbl) = x(ModData(iMod)%Vars%x(j)%iGbl) + dx(ModData(iMod)%Vars%x(j)%iGbl)
         case (VF_AngularDisp)
            ! Add delta x to x and limit to between -2pi and 2pi
            ! xn(ModData(iMod)%Vars%x(j)%iGbl) = mod(x(ModData(iMod)%Vars%x(j)%iGbl) + dx(ModData(iMod)%Vars%x(j)%iGbl), TwoPi_R8)
            xn(ModData(iMod)%Vars%x(j)%iGbl) = x(ModData(iMod)%Vars%x(j)%iGbl) + dx(ModData(iMod)%Vars%x(j)%iGbl)
         case (VF_Orientation)
            ! Compose WM components
            do k = 1, size(ModData(iMod)%Vars%x(j)%iGbl), 3
               iGbl = ModData(iMod)%Vars%x(j)%iGbl(k:k + 2)
               xn(iGbl) = wm_compose(dx(iGbl), x(iGbl))
            end do
         end select

      end do
   end do
end subroutine

subroutine AddDeltaToInputs(ModData, ModOrder, u, du, un)
   use ModVar
   type(ModDataType), intent(in)    :: ModData(:)  !< Solution variables from modules
   integer(IntKi), intent(in)       :: ModOrder(:) !< Array of module indices to evaluate
   real(R8Ki), intent(in)           :: du(:)
   real(R8Ki), intent(in)           :: u(:)
   real(R8Ki), intent(out)          :: un(:)
   integer(IntKi)                   :: iMod, iIns
   integer(IntKi)                   :: i, j, k, iGbl(3)

   ! Loop through modules in order
   do i = 1, size(ModOrder)

      iMod = ModOrder(i)
      iIns = ModData(iMod)%Instance

      ! Loop through variables
      do j = 1, size(ModData(iMod)%Vars%u)

         ! Select based on field type
         select case (ModData(iMod)%Vars%u(j)%Field)
         case (VF_TransDisp, VF_TransVel, VF_TransAcc, VF_AngularVel, VF_AngularAcc)
            ! Add delta u to u
            un(ModData(iMod)%Vars%u(j)%iGbl) = u(ModData(iMod)%Vars%u(j)%iGbl) + du(ModData(iMod)%Vars%u(j)%iGbl)
         case (VF_AngularDisp)
            ! Add delta u to u and limit to between -2pi and 2pi
            ! un(ModData(iMod)%Vars%u(j)%iGbl) = mod(u(ModData(iMod)%Vars%u(j)%iGbl) + du(ModData(iMod)%Vars%u(j)%iGbl), TwoPi_R8)
            un(ModData(iMod)%Vars%u(j)%iGbl) = u(ModData(iMod)%Vars%u(j)%iGbl) + du(ModData(iMod)%Vars%u(j)%iGbl)
         case (VF_Orientation)
            ! Compose WM components
            do k = 1, size(ModData(iMod)%Vars%u(j)%iGbl), 3
               iGbl = ModData(iMod)%Vars%u(j)%iGbl(k:k + 2)
               un(iGbl) = wm_compose(du(iGbl), u(iGbl))
            end do
         end select

      end do
   end do
end subroutine

subroutine TransferFromModules(ModOrder, ModData, this_state, T, x, u, y)
   integer(IntKi), intent(in)                         :: ModOrder(:) !< Array of module indices to evaluate
   type(ModDataType), intent(in)                      :: ModData(:)  !< Solution variables from modules
   integer(IntKi), intent(in)                         :: this_state  !< State index
   type(FAST_TurbineType), intent(inout)              :: T           !< Turbine type
   real(R8Ki), allocatable, optional, intent(inout)   :: x(:), u(:), y(:)
   integer(IntKi)                                     :: i, iMod, iIns

   if (present(x)) then
      do i = 1, size(ModOrder)
         iMod = ModOrder(i)
         iIns = ModData(iMod)%Instance
         select case (ModData(iMod)%ID)
         case (Module_ED)
            call ED_PackStateValues(T%ED%p, T%ED%x(this_state), T%ED%m%Vals%x)
            x(T%ED%p%Vars%ixg) = T%ED%m%Vals%x
         case (Module_BD)
            ! BD_PackStateValues(BD%p(iIns), BD%x(this_state, iIns), BD%m(iIns)%Vals%x)
            ! x(BD%p%Vars%ix) = BD%m%Vals%x
         case (Module_SD)
            ! call SD_PackStateValues(SD%p, SD%x(this_state), SD%m%Vals%x)
            ! x(SD%p%Vars%ix) = SD%m%Vals%x
         end select
      end do
   end if

   if (present(u)) then
      do i = 1, size(ModOrder)
         iMod = ModOrder(i)
         iIns = ModData(iMod)%Instance
         select case (ModData(iMod)%ID)
         case (Module_ED)
            call ED_PackInputValues(T%ED%p, T%ED%Input(1), T%ED%m%Vals%u)
            u(T%ED%p%Vars%iug) = T%ED%m%Vals%u
         case (Module_BD)
            ! BD_PackInputValues(BD%p(iIns), BD%Input(iIns), BD%m(iIns)%Vals%u)
            ! u(BD%p%Vars%iu) = BD%m%Vals%u
         case (Module_SD)
            ! call SD_PackInputValues(SD%p, SD%Input, SD%m%Vals%u)
            ! u(SD%p%Vars%iu) = SD%m%Vals%u
         end select
      end do
   end if

   if (present(y)) then
      do i = 1, size(ModOrder)
         iMod = ModOrder(i)
         iIns = ModData(iMod)%Instance
         select case (ModData(iMod)%ID)
         case (Module_ED)
            call ED_PackOutputValues(T%ED%p, T%ED%y, T%ED%m%Vals%y)
            y(T%ED%p%Vars%iyg) = T%ED%m%Vals%y
         case (Module_BD)
            ! BD_PackOutputValues(BD%p(iIns), BD%y(iIns), BD%m(iIns)%Vals%y)
            ! y(BD%p%Vars%iy) = BD%m%Vals%y
         case (Module_SD)
            ! call SD_PackOutputValues(SD%p, SD%y, SD%m%Vals%y)
            ! y(SD%p%Vars%iy) = SD%m%Vals%y
         end select
      end do
   end if

   ! if (present(dxdt)) then
   !    do i = 1, size(ModOrder)
   !       iMod = ModOrder(i)
   !       iIns = ModData(iMod)%Instance
   !       select case (ModData(iMod)%ID)
   !       case (Module_ED)
   !          call ED_PackStateValues(T%ED%p, T%ED%dxdt, T%ED%m%Vals%dxdt)
   !          dxdt(T%ED%p%Vars%ixg) = T%ED%m%Vals%dxdt
   !       case (Module_BD)
   !          ! BD_PackStateValues(BD%p(iIns), BD%dxdt(iIns), BD%m(iIns)%Vals%dxdt)
   !          ! dxdt(BD%p%Vars%ix) = BD%m%Vals%dxdt
   !       case (Module_SD)
   !          ! call SD_PackStateValues(SD%p, SD%dxdt, SD%m%Vals%dxdt)
   !          ! dxdt(SD%p%Vars%ix) = SD%m%Vals%dxdt
   !       end select
   !    end do
   ! end if

   ! if (present(dYdx)) then
   !    do i = 1, size(ModOrder)
   !       iMod = ModOrder(i)
   !       iIns = ModData(iMod)%Instance
   !       select case (ModData(iMod)%ID)
   !       case (Module_ED)
   !          dYdx(T%ED%p%Vars%iyg, T%ED%p%Vars%ixg) = T%ED%m%Vals%dYdx
   !       case (Module_BD)
   !          ! dYdx(BD%p%Vars%iy, BD%p%Vars%ix) = BD%m%Vals%dYdx
   !       case (Module_SD)
   !          ! dYdx(SD%p%Vars%iy, SD%p%Vars%ix) = SD%m%Vals%dYdx
   !       end select
   !    end do
   ! end if

   ! if (present(dXdx)) then
   !    do i = 1, size(ModOrder)
   !       iMod = ModOrder(i)
   !       iIns = ModData(iMod)%Instance
   !       select case (ModData(iMod)%ID)
   !       case (Module_ED)
   !          dXdx(T%ED%p%Vars%ixg, T%ED%p%Vars%ixg) = T%ED%m%Vals%dXdx
   !       case (Module_BD)
   !          ! dXdx(BD%p%Vars%ix, BD%p%Vars%ix) = BD%m%Vals%dXdx
   !       case (Module_SD)
   !          ! dXdx(SD%p%Vars%ix, SD%p%Vars%ix) = SD%m%Vals%dXdx
   !       end select
   !    end do
   ! end if

   ! if (present(dYdu)) then
   !    do i = 1, size(ModOrder)
   !       iMod = ModOrder(i)
   !       iIns = ModData(iMod)%Instance
   !       select case (ModData(iMod)%ID)
   !       case (Module_ED)
   !          dYdu(T%ED%p%Vars%iyg, T%ED%p%Vars%iug) = T%ED%m%Vals%dYdu
   !       case (Module_BD)
   !          ! dYdu(BD%p%Vars%iy, BD%p%Vars%iu) = BD%m%Vals%dYdu
   !       case (Module_SD)
   !          ! dYdu(SD%p%Vars%iy, SD%p%Vars%iu) = SD%m%Vals%dYdu
   !       end select
   !    end do
   ! end if

   ! if (present(dXdu)) then
   !    do i = 1, size(ModOrder)
   !       iMod = ModOrder(i)
   !       iIns = ModData(iMod)%Instance
   !       select case (ModData(iMod)%ID)
   !       case (Module_ED)
   !          dXdu(T%ED%p%Vars%ixg, T%ED%p%Vars%iug) = T%ED%m%Vals%dXdu
   !       case (Module_BD)
   !          ! dXdu(BD%p%Vars%ix, BD%p%Vars%iu) = BD%m%Vals%dXdu
   !       case (Module_SD)
   !          ! dXdu(SD%p%Vars%ix, SD%p%Vars%iu) = SD%m%Vals%dXdu
   !       end select
   !    end do
   ! end if

end subroutine

subroutine TransferToModules(ModOrder, ModData, this_state, T, x, u)
   integer(IntKi), intent(in)                         :: ModOrder(:) !< Array of module indices to evaluate
   type(ModDataType), intent(in)                      :: ModData(:)  !< Solution variables from modules
   integer(IntKi), intent(in)                         :: this_state  !< State index
   type(FAST_TurbineType), intent(inout)              :: T           !< Turbine type
   real(R8Ki), allocatable, optional, intent(inout)   :: x(:), u(:)
   integer(IntKi)                                     :: i, iMod, iIns

   if (present(x)) then
      do i = 1, size(ModOrder)
         iMod = ModOrder(i)
         iIns = ModData(iMod)%Instance
         select case (ModData(iMod)%ID)
         case (Module_ED)
            call ED_UnpackStateValues(T%ED%p, x(T%ED%p%Vars%ixg), T%ED%x(this_state))
         case (Module_BD)
            ! BD_UnpackStateValues(BD%p(iIns), x(BD%p%Vars%ix), BD%x(this_state, iIns))
         case (Module_SD)
            ! call SD_UnpackStateValues(SD%p, x(SD%p%Vars%ix), SD%x(this_state))
         end select
      end do
   end if

   if (present(u)) then
      do i = 1, size(ModOrder)
         iMod = ModOrder(i)
         iIns = ModData(iMod)%Instance
         select case (ModData(iMod)%ID)
         case (Module_ED)
            call ED_UnpackInputValues(T%ED%p, u(T%ED%p%Vars%iug), T%ED%Input(1))
         case (Module_BD)
            ! BD_UnpackInputValues(BD%p(iIns), u(BD%p%Vars%iu), BD%Input(1, iIns))
         case (Module_SD)
            ! call SD_UnpackInputValues(SD%p, u(SD%p%Vars%iu), SD%Input(1))
         end select
      end do
   end if

end subroutine

subroutine WriteOutputToFile(n_t_global, t_global, T, ErrStat, ErrMsg)

   integer(IntKi), intent(in)             :: n_t_global  !< Current global time step
   real(DbKi), intent(in)                 :: t_global    !< Current global time
   type(FAST_TurbineType), intent(inout)  :: T           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat     !< Error status of the operation
   character(*), intent(out)              :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   character(*), parameter                :: RoutineName = 'WriteOutputToFile'
   integer(IntKi)                         :: ErrStat2
   character(ErrMsgLen)                   :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ""

   ! If output requested for this step, time-series channel data to glue-code output file
   if (T%y_FAST%WriteThisStep) then
      ! call WrOutputLine(t_global, T%p_FAST, T%y_FAST, T%IfW%y%WriteOutput, &
      !                   T%OpFM%y%WriteOutput, T%ED%y%WriteOutput, &
      !                   T%AD%y, T%SrvD%y%WriteOutput, T%SeaSt%y%WriteOutput, &
      !                   T%HD%y%WriteOutput, T%SD%y%WriteOutput, T%ExtPtfm%y%WriteOutput, &
      !                   T%MAP%y%WriteOutput, T%FEAM%y%WriteOutput, T%MD%y%WriteOutput, &
      !                   T%Orca%y%WriteOutput, T%IceF%y%WriteOutput, T%IceD%y, T%BD%y, ErrStat2, ErrMsg2)
      ! call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   end if

end subroutine

end module
