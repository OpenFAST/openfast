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

! Public functions
public Solver_Init, Solver_Step0, Solver_Step

! q matrix column indices (displacement, velocity, accel, algorithmic accel)
integer(IntKi), parameter  :: COL_D = 1, COL_V = 2, COL_A = 3, COL_AA = 4

! Tight coupling modules
integer(IntKi), parameter  :: TC_Modules(*) = [Module_ED, Module_BD, Module_SD]

! Debugging
logical, parameter         :: DebugSolver = .false.
integer(IntKi)             :: DebugUn = -1
character(*), parameter    :: DebugFile = 'solver.dbg'
logical, parameter         :: DebugJacobian = .false.
integer(IntKi)             :: MatrixUn = -1

contains

subroutine Solver_Init(p, m, Mods, ErrStat, ErrMsg)

   type(TC_ParameterType), intent(inout)     :: p           !< Parameters
   type(TC_MiscVarType), intent(out)         :: m           !< Misc variables for optimization (not copied in glue code)
   type(ModDataType), intent(inout)          :: Mods(:)     !< Module data
   integer(IntKi), intent(out)               :: ErrStat     !< Error status of the operation
   character(*), intent(out)                 :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   character(*), parameter                   :: RoutineName = 'Solver_Init'
   integer(IntKi)                            :: ErrStat2    ! local error status
   character(ErrMsgLen)                      :: ErrMsg2     ! local error message
   integer(IntKi)                            :: i, j, k, n
   integer(IntKi)                            :: NumX, NumQ, NumU, NumY, NumJ
   integer(IntKi), allocatable               :: modIDs(:), modInds(:)
   type(TC_MappingType)                      :: MeshMap

   !----------------------------------------------------------------------------
   ! Module ordering for solve
   !----------------------------------------------------------------------------

   ! Create array of indices for Mods array
   modInds = [(i, i=1, size(Mods))]

   ! Get array of module IDs
   modIDs = [(Mods(i)%ID, i=1, size(Mods))]

   ! Indicies of all modules
   p%iModAll = [pack(modInds, ModIDs == Module_SrvD), &
                pack(modInds, ModIDs == Module_IfW), &
                pack(modInds, ModIDs == Module_AD), &
                pack(modInds, ModIDs == Module_ED), &
                pack(modInds, ModIDs == Module_BD), &
                pack(modInds, ModIDs == Module_SD)]

   ! Indicies of tight coupling modules
   p%iModTC = [pack(modInds, ModIDs == Module_ED), &
               pack(modInds, ModIDs == Module_BD), &
               pack(modInds, ModIDs == Module_SD)]

   ! Indicies of Option 2 modules
   p%iModOpt2 = [pack(modInds, ModIDs == Module_SrvD), &
                 pack(modInds, ModIDs == Module_ED), &
                 pack(modInds, ModIDs == Module_BD), &
                 pack(modInds, ModIDs == Module_SD), &
                 pack(modInds, ModIDs == Module_IfW), &
                 pack(modInds, ModIDs == Module_AD), &
                 pack(modInds, ModIDs == Module_FEAM), &
                 pack(modInds, ModIDs == Module_IceD), &
                 pack(modInds, ModIDs == Module_IceF), &
                 pack(modInds, ModIDs == Module_MAP), &
                 pack(modInds, ModIDs == Module_MD)]

   ! Indicies of Option 1 modules
   p%iModOpt1 = [pack(modInds, ModIDs == Module_ED), &
                 pack(modInds, ModIDs == Module_BD), &
                 pack(modInds, ModIDs == Module_SD), &
                 pack(modInds, ModIDs == Module_ExtPtfm), &
                 pack(modInds, ModIDs == Module_HD), &
                 pack(modInds, ModIDs == Module_Orca)]

   ! Indicies of Option 1 modules that were not in Option 2
   ! These modules need to do update states and calc output before Option 1 solve
   p%iModOpt1US = [pack(modInds, ModIDs == Module_ExtPtfm), &
                   pack(modInds, ModIDs == Module_HD), &
                   pack(modInds, ModIDs == Module_Orca)]

   ! Indices of BeamDyn modules
   p%iModBD = [pack(modInds, ModIDs == Module_BD)]

   !----------------------------------------------------------------------------
   ! Initialize mesh mappings (must be done before calculating global indices)
   !----------------------------------------------------------------------------

   call DefineMappings(m%Mappings, Mods, ErrStat2, ErrMsg2)
   if (Failed()) return

   !----------------------------------------------------------------------------
   ! Calculate variable global indices and related index arrays
   !----------------------------------------------------------------------------

   call CalcVarGlobalIndices(p, Mods, NumX, NumU, NumY, NumQ, NumJ, ErrStat2, ErrMsg2)
   if (Failed()) return

   !----------------------------------------------------------------------------
   ! Allocate state, input, and output storage
   !----------------------------------------------------------------------------

   ! State arrays
   call AllocAry(m%x, NumX, "m%x", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%xn, NumX, "m%xn", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%dx, NumX, "m%dx", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%dxdt, NumX, "m%dxdt", ErrStat2, ErrMsg2); if (Failed()) return

   ! State matrices
   call AllocAry(m%q, NumQ, 4, "m%q", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%qn, NumQ, 4, "m%qn", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%dq, NumQ, 4, "m%dq", ErrStat2, ErrMsg2); if (Failed()) return

   ! Inputs array
   call AllocAry(m%u, NumU, "m%u", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%un, NumU, "m%un", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%du, NumU, "m%du", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%u_tmp, NumU, "m%u_tmp", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%UDiff, NumU, "m%UDiff", ErrStat2, ErrMsg2); if (Failed()) return

   ! Output arrays
   call AllocAry(m%y, NumY, "m%y", ErrStat2, ErrMsg2); if (Failed()) return

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

   !----------------------------------------------------------------------------
   ! Allocate solver Jacobian matrix and right hand side
   !----------------------------------------------------------------------------

   ! Allocate Jacobian matrix, RHS/X matrix, Pivot array
   call AllocAry(m%Jac, NumJ, NumJ, "m%Jac", ErrStat, ErrMsg); if (Failed()) return
   call AllocAry(m%XB, NumJ, 1, "m%XB", ErrStat, ErrMsg); if (Failed()) return
   call AllocAry(m%IPIV, NumJ, "m%IPIV", ErrStat, ErrMsg); if (Failed()) return
   m%Jac = 0.0_R8Ki

   ! Initialize iterations and steps until Jacobian is calculated to zero
   ! so it is calculated in first step
   m%IterUntilUJac = 0
   m%StepsUntilUJac = 0

   !----------------------------------------------------------------------------
   ! Calculate generalized alpha parameters
   !----------------------------------------------------------------------------

   ! Acceleration blending term between current and next value, used during
   ! state prediction (1.0 = constant acceleration, 0.0 = no acceleration).
   p%AccBlend = 1.0_R8Ki

   p%AlphaM = (2.0_R8Ki*p%RhoInf - 1.0_R8Ki)/(p%RhoInf + 1.0_R8Ki)
   p%AlphaF = p%RhoInf/(p%RhoInf + 1.0_R8Ki)
   p%Gamma = 0.5_R8Ki - p%AlphaM + p%AlphaF
   p%Beta = (1.0_R8Ki - p%AlphaM + p%AlphaF)**2.0_R8Ki/4.0_R8Ki

   ! Constants used in Jacobian construction and state prediction
   p%C(1) = (1.0_R8Ki - p%AlphaF)/(1.0_R8Ki - p%AlphaM)
   p%C(2) = p%DT*p%Gamma*p%C(1)
   p%C(3) = p%DT**2*p%Beta*p%C(1)
   p%C(4) = p%DT**2*(0.5_R8Ki - p%Beta)
   p%C(5) = p%DT**2*p%Beta
   p%C(6) = p%DT*(1.0_R8Ki - p%Gamma)
   p%C(7) = p%DT*p%Gamma

   !----------------------------------------------------------------------------
   ! Write debug info to file
   !----------------------------------------------------------------------------

   if (DebugSolver) then
      call GetNewUnit(DebugUn, ErrStat2, ErrMsg2); if (Failed()) return
      call OpenFOutFile(DebugUn, DebugFile, ErrStat2, ErrMsg2); if (Failed()) return
      call Solver_Init_Debug(p, m, Mods)
   end if

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine CalcVarGlobalIndices(p, Mods, NumX, NumU, NumY, NumQ, NumJ, ErrStat, ErrMsg)
   type(TC_ParameterType), intent(inout)  :: p           !< Parameters
   type(ModDataType), intent(inout)       :: Mods(:)     !< Module data
   integer(IntKi), intent(out)            :: NumX, NumU, NumY, NumQ, NumJ
   integer(IntKi), intent(out)            :: ErrStat     !< Error status of the operation
   character(*), intent(out)              :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   character(*), parameter                :: RoutineName = 'Solver_Init'
   integer(IntKi)                         :: ErrStat2    ! local error status
   character(ErrMsgLen)                   :: ErrMsg2     ! local error message
   integer(IntKi), allocatable            :: modIDs(:), vec1(:), vec2(:), iuLoad(:)
   integer(IntKi)                         :: i, j, k, n

   ErrStat = ErrID_None
   ErrMsg = ''

   !----------------------------------------------------------------------------
   ! Calculate state variable global indices and transfer arrays
   ! for tight coupling and option one modules
   !----------------------------------------------------------------------------

   ! Initialize number of state variables to zero
   NumX = 0

   ! Loop through tight coupling modules, set Solve flag
   do i = 1, size(p%iModTC)
      associate (Mod => Mods(p%iModTC(i)))
         do j = 1, size(Mod%Vars%x)
            call SetFlags(Mod%Vars%x(j), VF_Solve)
         end do
      end associate
   end do

   ! Loop through x displacement variables (DerivOrder == 0), set global index
   do i = 1, size(p%iModTC)
      associate (Mod => Mods(p%iModTC(i)))
         do j = 1, size(Mod%Vars%x)
            if (Mod%Vars%x(j)%DerivOrder == 0) then
               Mod%Vars%x(j)%iGblSol = [(NumX + k, k=1, Mod%Vars%x(j)%Num)]
               NumX = NumX + Mod%Vars%x(j)%Num
            end if
         end do
      end associate
   end do

   ! Start and end indices for first derivative of X
   p%iX1 = [1, NumX]

   ! Loop through x velocity variables (DerivOrder == 1), set global index
   do i = 1, size(p%iModTC)
      associate (Mod => Mods(p%iModTC(i)))
         do j = 1, size(Mod%Vars%x)
            if (Mod%Vars%x(j)%DerivOrder == 1) then
               Mod%Vars%x(j)%iGblSol = [(NumX + k, k=1, Mod%Vars%x(j)%Num)]
               NumX = NumX + Mod%Vars%x(j)%Num
            end if
         end do
      end associate
   end do

   ! Start and end indices for second derivative of X
   p%iX2 = [p%iX1(2) + 1, NumX]

   ! Loop through x variables, and collect global and local indices
   do i = 1, size(p%iModTC)
      allocate (vec1(0), vec2(0))
      associate (Mod => Mods(p%iModTC(i)))
         do j = 1, size(Mod%Vars%x)
            vec1 = [vec1, Mod%Vars%x(j)%iLoc]
            vec2 = [vec2, Mod%Vars%x(j)%iGblSol]
         end do
         call AllocAry(Mod%ixs, 2, size(vec1), "Mod%ixs", ErrStat2, ErrMsg2); if (Failed()) return
         Mod%ixs(1, :) = vec1
         Mod%ixs(2, :) = vec2
      end associate
      deallocate (vec1, vec2)
   end do

   !----------------------------------------------------------------------------
   ! Calculate input variable global indices and transfer arrays
   ! for tight coupling and option one modules
   !----------------------------------------------------------------------------

   allocate (iuLoad(0))

   ! Initialze number of input variables to zero
   NumU = 0

   ! Loop through tight coupling modules and calculate global indices
   do i = 1, size(p%iModTC)
      do j = 1, size(Mods(p%iModTC(i))%Vars%u)
         associate (Var => Mods(p%iModTC(i))%Vars%u(j))
            if ((.not. allocated(Var%iGblSol)) .and. (iand(Var%Flags, VF_Solve) > 0)) then
               Var%iGblSol = [(NumU + k, k=1, Var%Num)]
               NumU = NumU + Var%Num
            end if
         end associate
      end do
   end do

   ! Save number of tight coupling inputs
   p%iUT = [1, NumU]

   ! Loop through option 1 modules and calculate global indices
   do i = 1, size(p%iModOpt1)
      do j = 1, size(Mods(p%iModOpt1(i))%Vars%u)
         associate (Var => Mods(p%iModOpt1(i))%Vars%u(j))
            if ((.not. allocated(Var%iGblSol)) .and. (iand(Var%Flags, VF_Solve) > 0)) then
               Var%iGblSol = [(NumU + k, k=1, Var%Num)]
               NumU = NumU + Var%Num
            end if
         end associate
      end do
   end do

   ! Save number of option 1 inputs
   p%iU1 = [p%iUT(2) + 1, NumU]

   ! Loop through all modules
   do i = 1, size(Mods)
      allocate (vec1(0), vec2(0))
      do j = 1, size(Mods(i)%Vars%u)
         if (.not. allocated(Mods(i)%Vars%u(j)%iGblSol)) cycle
         vec1 = [vec1, Mods(i)%Vars%u(j)%iLoc]
         vec2 = [vec2, Mods(i)%Vars%u(j)%iGblSol]
         select case (Mods(i)%Vars%u(j)%Field)
         case (VF_Force, VF_Moment)
            iuLoad = [iuLoad, Mods(i)%Vars%u(j)%iGblSol]
         end select
      end do
      call AllocAry(Mods(i)%ius, 2, size(vec1), "Mods(i)%ius", ErrStat2, ErrMsg2); if (Failed()) return
      Mods(i)%ius(1, :) = vec1
      Mods(i)%ius(2, :) = vec2
      deallocate (vec1, vec2)
   end do

   !----------------------------------------------------------------------------
   ! Calculate output variable categories and indices
   !----------------------------------------------------------------------------

   ! Initialize the number of output variables
   NumY = 0

   ! Loop through tight coupling modules and calculate global indices
   do i = 1, size(p%iModTC)
      do j = 1, size(Mods(p%iModTC(i))%Vars%y)
         associate (Var => Mods(p%iModTC(i))%Vars%y(j))
            if ((.not. allocated(Var%iGblSol)) .and. (iand(Var%Flags, VF_Solve) > 0)) then
               Var%iGblSol = [(NumY + k, k=1, Var%Num)]
               NumY = NumY + Var%Num
            end if
         end associate
      end do
   end do

   ! Save number of tight coupling inputs
   p%iyT = [1, NumY]

   ! Loop through option 1 modules and calculate global indices
   do i = 1, size(p%iModOpt1)
      do j = 1, size(Mods(p%iModOpt1(i))%Vars%y)
         associate (Var => Mods(p%iModOpt1(i))%Vars%y(j))
            if ((.not. allocated(Var%iGblSol)) .and. (iand(Var%Flags, VF_Solve) > 0)) then
               Var%iGblSol = [(NumY + k, k=1, Var%Num)]
               NumY = NumY + Var%Num
            end if
         end associate
      end do
   end do

   ! Calculate number of option 1 inputs
   p%iy1 = [p%iyT(2) + 1, NumY]

   ! Loop through all modules
   do i = 1, size(Mods)
      allocate (vec1(0), vec2(0))
      do j = 1, size(Mods(i)%Vars%y)
         if (.not. allocated(Mods(i)%Vars%y(j)%iGblSol)) cycle
         vec1 = [vec1, Mods(i)%Vars%y(j)%iLoc]
         vec2 = [vec2, Mods(i)%Vars%y(j)%iGblSol]
      end do
      call AllocAry(Mods(i)%iys, 2, size(vec1), "Mods(i)%iys", ErrStat2, ErrMsg2); if (Failed()) return
      Mods(i)%iys(1, :) = vec1
      Mods(i)%iys(2, :) = vec2
      deallocate (vec1, vec2)
   end do

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

      ! Skip modules that aren't in tight coupling
      if (all(Mods(i)%ID /= TC_Modules)) cycle

      ! Loop through state variables
      do j = 1, size(Mods(i)%Vars%x)

         ! Skip variables which already have a q index
         if (allocated(Mods(i)%Vars%x(j)%iq)) cycle

         ! Skip load variables (force and moment)
         if (any(Mods(i)%Vars%x(j)%Field == LoadFields)) cycle

         ! Set q index for variable and update number
         Mods(i)%Vars%x(j)%iq = [(NumQ + k, k=1, Mods(i)%Vars%x(j)%Num)]
         NumQ = NumQ + Mods(i)%Vars%x(j)%Num

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

            ! Copy q row indices
            Mods(i)%Vars%x(k)%iq = Mods(i)%Vars%x(j)%iq

         end do
      end do

      ! Loop through state variables and build mapping array between x and q
      ! ixqd is 3xN where each row is [global x array index, q matrix row, q matrix col]
      do j = 1, size(Mods(i)%Vars%x)
         do k = 1, Mods(i)%Vars%x(j)%Num
            n = n + 1
            p%ixqd(:, n) = [Mods(i)%Vars%x(j)%iGblSol(k), Mods(i)%Vars%x(j)%iq(k), Mods(i)%Vars%x(j)%DerivOrder + 1]
         end do
      end do

   end do

   ! Remove unused x->q indicies
   p%ixqd = p%ixqd(:, 1:n)

   !----------------------------------------------------------------------------
   ! Jacobian indices and ranges
   !----------------------------------------------------------------------------

   ! Calculate size of Jacobian matrix
   NumJ = NumQ + NumU

   ! Get start and end indicies for state part of Jacobian
   p%iJX = [1, NumQ]

   ! Get start and end indicies for input part of Jacobian
   p%iJUT = p%iUT + NumQ
   p%iJU = [1, NumU] + NumQ

   ! Get Jacobian indicies containing loads
   p%iJL = iuLoad + NumQ

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine DefineMappings(Mappings, Mods, ErrStat, ErrMsg)
   type(TC_MappingType), allocatable, intent(inout)   :: Mappings(:)
   type(ModDataType), intent(inout)                   :: Mods(:)     !< Module data
   integer(IntKi), intent(out)                        :: ErrStat     !< Error status of the operation
   character(*), intent(out)                          :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   character(*), parameter   :: RoutineName = 'DefineMappings'
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

subroutine TransferXtoQ(ixqd, x, q)
   integer(IntKi), intent(in) :: ixqd(:, :)
   real(R8Ki), intent(in)     :: x(:)
   real(R8Ki), intent(inout)  :: q(:, :)
   integer(IntKi)             :: i
   do i = 1, size(ixqd, dim=2)
      q(ixqd(2, i), ixqd(3, i)) = x(ixqd(1, i))
   end do
end subroutine

subroutine TransferQtoX(ixqd, q, x)
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
end function

subroutine Solver_Step0(p, m, ModData, Turbine, ErrStat, ErrMsg)
   type(TC_ParameterType), intent(in)        :: p           !< Parameters
   type(TC_MiscVarType), intent(inout)       :: m           !< Misc variables
   type(ModDataType), intent(in)             :: ModData(:)  !< Module data
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
   call TransferXtoQ(p%ixqd, m%x, m%qn)

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
         call FAST_InputSolve(ModData(p%iModAll(i)), m%Mappings, 1, &
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
   m%q = 0.0_R8Ki

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
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine Solver_Step(n_t_global, t_initial, p, m, Mods, Turbine, ErrStat, ErrMsg)
   integer(IntKi), intent(in)                :: n_t_global  !< global time step
   real(DbKi), intent(in)                    :: t_initial   !< Initial simulation time
   type(TC_ParameterType), intent(in)        :: p           !< Parameters
   type(TC_MiscVarType), intent(inout)       :: m           !< Misc variables
   type(ModDataType), intent(in)             :: Mods(:)  !< Module data
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
   integer(IntKi)             :: i, j

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
                     p%AlphaM*m%qn(:, COL_AA))/(1.0_R8Ki - p%AlphaM)

   ! Calculate displacements for the next step
   m%q(:, COL_D) = m%qn(:, COL_D) + &
                   p%DT*m%qn(:, COL_V) + &
                   p%C(4)*m%qn(:, COL_AA) + &
                   p%C(5)*m%q(:, COL_AA)

   ! Calculate velocities for the next step
   m%q(:, COL_V) = m%qn(:, COL_V) + &
                   p%C(6)*m%qn(:, COL_AA) + &
                   p%C(7)*m%q(:, COL_AA)

   ! Calculate difference between new and old states
   m%dq = m%q - m%qn

   ! Transfer delta state matrix to delta x array
   m%dx = 0.0_R8Ki
   call TransferQtoX(p%ixqd, m%dq, m%dx)

   ! Add delta to x array to get new states (respect variable fields)
   ! required for orientation fields in states
   call AddDeltaToStates(Mods, p%iModTC, m%dx, m%x)

   ! Update state matrix with updated state values
   call TransferXtoQ(p%ixqd, m%x, m%q)

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

      ! Copy state for correction step
      m%qn = m%q
      m%xn = m%x

      ! Reset mappings updated flags
      m%Mappings%Updated = .false.

      !-------------------------------------------------------------------------
      ! Option 2 Solve
      !-------------------------------------------------------------------------

      ! Loop through Option 2 modules
      do i = 1, size(p%iModOpt2)
         call FAST_InputSolve(Mods(p%iModOpt2(i)), m%Mappings, 1, &
                              Turbine, ErrStat2, ErrMsg2); if (Failed()) return
         call FAST_UpdateStates(Mods(p%iModOpt2(i)), t_initial, n_t_global, m%xn, m%qn, &
                                Turbine, ErrStat2, ErrMsg2); if (Failed()) return
         call FAST_CalcOutput(Mods(p%iModOpt2(i)), t_global_next, STATE_PRED, &
                              Turbine, ErrStat2, ErrMsg2); if (Failed()) return
         if (i < 2) then
            call FAST_MapOutputs(Mods(p%iModOpt2(i)), m%Mappings, &
                                 Turbine, ErrStat2, ErrMsg2); if (Failed()) return
         end if
      end do

      !-------------------------------------------------------------------------
      ! Option 1 Solve
      !-------------------------------------------------------------------------

      ! Get inputs and update states for Option 1 modules not in Option 2
      do i = 1, size(p%iModOpt1US)
         call FAST_InputSolve(Mods(p%iModOpt1US(i)), m%Mappings, 1, &
                              Turbine, ErrStat2, ErrMsg2); if (Failed()) return
         call FAST_UpdateStates(Mods(p%iModOpt1US(i)), t_initial, n_t_global, m%xn, m%qn, &
                                Turbine, ErrStat2, ErrMsg2); if (Failed()) return
      end do

      ! Populate state matrix with latest values from state array in case it
      ! changed during FAST_UpdateStates calls
      call TransferXtoQ(p%ixqd, m%xn, m%qn)

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

         !----------------------------------------------------------------------
         ! If iteration limit reached, exit loop
         !----------------------------------------------------------------------

         if (iterConv >= p%MaxConvIter) exit

         !----------------------------------------------------------------------
         ! Update Jacobian
         !----------------------------------------------------------------------

         ! If number of iterations or steps until Jacobian is to be updated
         ! is zero or less, or first solution step, then rebuild the Jacobian.
         ! Note: BuildJacobian resets these counters.
         if ((m%IterUntilUJac <= 0) .or. (m%StepsUntilUJac <= 0) .or. (n_t_global_next == 1)) then
            call BuildJacobian(p, m, Mods, t_global_next, n_t_global_next*100 + iterConv, &
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
         m%XB(p%iJX(1):p%iJX(2), 1) = m%qn(:, COL_A) - m%dxdt(p%iX2(1):p%iX2(2))

         ! Transfer Option 1 outputs to temporary inputs and collect into u_tmp
         do i = 1, size(p%iModOpt1)
            call FAST_MapOutputs(Mods(p%iModOpt1(i)), m%Mappings, &
                                 Turbine, ErrStat2, ErrMsg2); if (Failed()) return
            call FAST_InputSolve(Mods(p%iModOpt1(i)), m%Mappings, 2, &
                                 Turbine, ErrStat2, ErrMsg2); if (Failed()) return
         end do
         call PackModuleUs(Mods, p%iModOpt1, Turbine, m%u_tmp)

         ! Calculate difference in U for all Option 1 modules (un - u_tmp)
         ! and add to RHS for TC and Option 1 modules
         call ComputeDeltaU(Mods, p%iModOpt1, m%un, m%u_tmp, m%UDiff)
         m%XB(p%iJU(1):p%iJU(2), 1) = m%UDiff

         ! Apply conditioning factor to loads in RHS
         m%XB(p%iJL, 1) = m%XB(p%iJL, 1)/p%Scale_UJac

         !----------------------------------------------------------------------
         ! Solve for state and input perturbations
         !----------------------------------------------------------------------

         ! Solve Jacobian and RHS
         call LAPACK_getrs('N', size(m%Jac, 1), m%Jac, m%IPIV, m%XB, ErrStat2, ErrMsg2); if (Failed()) return

         !----------------------------------------------------------------------
         ! Check perturbations for convergence and exit if below tolerance
         !----------------------------------------------------------------------

         delta_norm = TwoNorm(m%XB(:, 1))/size(m%XB)

         if (DebugSolver) then
            call Solver_Step_Debug(p, m, n_t_global_next, iterCorr, iterConv, delta_norm)
         end if

         if (delta_norm < p%ConvTol) exit

         ! Remove conditioning
         m%XB(p%iJL, 1) = m%XB(p%iJL, 1)*p%Scale_UJac

         !----------------------------------------------------------------------
         ! Update State for Tight Coupling modules
         !----------------------------------------------------------------------

         ! Calculate change in state matrix
         m%dq(:, COL_D) = -p%C(3)*m%XB(p%iJX(1):p%iJX(2), 1)
         m%dq(:, COL_V) = -p%C(2)*m%XB(p%iJX(1):p%iJX(2), 1)
         m%dq(:, COL_A) = -m%XB(p%iJX(1):p%iJX(2), 1)
         m%dq(:, COL_AA) = -p%C(1)*m%XB(p%iJX(1):p%iJX(2), 1)

         ! Update state matrix with deltas
         m%qn = m%qn + m%dq

         ! Transfer change in q state matrix to change in x array
         m%dx = 0.0_R8Ki
         call TransferQtoX(p%ixqd, m%dq, m%dx)

         ! Add delta to x array to get new states (respect variable fields)
         call AddDeltaToStates(Mods, p%iModTC, m%dx, m%xn)

         ! Overwrites state matrix values that were changed in AddDeltaToStates
         call TransferXtoQ(p%ixqd, m%xn, m%qn)

         ! Transfer updated state to TC modules
         call UnpackModuleStates(Mods, p%iModTC, STATE_PRED, Turbine, x=m%xn)

         !----------------------------------------------------------------------
         ! Update inputs for Tight Coupling and Option 1 modules
         !----------------------------------------------------------------------

         ! Update change in inputs
         m%du = -m%XB(p%iJU(1):p%iJU(2), 1)

         ! Apply deltas to inputs, update modules
         call AddDeltaToInputs(Mods, p%iModOpt1, m%du, m%un)

         ! Transfer updated inputs to Option 1 modules
         call UnpackModuleInputs(Mods, p%iModOpt1, Turbine, u=m%un)

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

   ! Loop through BeamDyn instances
   do i = 1, size(p%iModBD)
      associate (Mod => Mods(p%iModBD(i)), &
                 p_BD => Turbine%BD%p(Mods(p%iModBD(i))%Ins), &
                 m_BD => Turbine%BD%m(Mods(p%iModBD(i))%Ins), &
                 u_BD => Turbine%BD%Input(1, Mods(p%iModBD(i))%Ins), &
                 x_BD => Turbine%BD%x(Mods(p%iModBD(i))%Ins, STATE_PRED), &
                 os_BD => Turbine%BD%OtherSt(Mods(p%iModBD(i))%Ins, STATE_PRED))

         ! Update accelerations and algorithmic accelerations
         do j = 1, size(p_BD%Vars%x)
            select case (p_BD%Vars%x(j)%Field)
            case (VF_TransDisp)
               os_BD%acc(1:3, p_BD%Vars%x(j)%iUsr(1)) = m%qn(p_BD%Vars%x(j)%iq, COL_A)
               os_BD%xcc(1:3, p_BD%Vars%x(j)%iUsr(1)) = m%qn(p_BD%Vars%x(j)%iq, COL_AA)
            case (VF_Orientation)
               os_BD%acc(4:6, p_BD%Vars%x(j)%iUsr(1)) = m%qn(p_BD%Vars%x(j)%iq, COL_A)
               os_BD%xcc(4:6, p_BD%Vars%x(j)%iUsr(1)) = m%qn(p_BD%Vars%x(j)%iq, COL_AA)
            end select
         end do

         ! Update the global reference
         call BD_UpdateGlobalRef(u_BD, p_BD, x_BD, os_BD, ErrStat, ErrMsg); if (Failed()) return

         ! Update accelerations and algorithmic accelerations
         do j = 1, size(p_BD%Vars%x)
            select case (p_BD%Vars%x(j)%Field)
            case (VF_TransDisp)
               m%qn(p_BD%Vars%x(j)%iq, COL_A) = os_BD%acc(1:3, p_BD%Vars%x(j)%iUsr(1))
               m%qn(p_BD%Vars%x(j)%iq, COL_AA) = os_BD%xcc(1:3, p_BD%Vars%x(j)%iUsr(1))
            case (VF_Orientation)
               m%qn(p_BD%Vars%x(j)%iq, COL_A) = os_BD%acc(4:6, p_BD%Vars%x(j)%iUsr(1))
               m%qn(p_BD%Vars%x(j)%iq, COL_AA) = os_BD%xcc(4:6, p_BD%Vars%x(j)%iUsr(1))
            end select
         end do

         ! Transfer updated states to solver
         call BD_PackStateValues(p_BD, x_BD, m_BD%Vals%x)
         call XferLocToGbl1D(Mod%ixs, m_BD%Vals%x, m%xn)
      end associate
   end do

   ! Update state matrix from state array
   call TransferXtoQ(p%ixqd, m%xn, m%qn)

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

subroutine BuildJacobian(p, m, Mods, this_time, iter, Turbine, ErrStat, ErrMsg)
   type(TC_ParameterType), intent(in)     :: p           !< Parameters
   type(TC_MiscVarType), intent(inout)    :: m           !< Misc variables
   type(ModDataType), intent(in)          :: Mods(:)     !< Array of module data
   real(DbKi), intent(in)                 :: this_time   !< Time
   integer(IntKi), intent(in)             :: iter
   type(FAST_TurbineType), intent(inout)  :: Turbine     !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter                :: RoutineName = 'BuildJacobian'
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

   ! Initialize Jacobian matrices
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
   call LAPACK_GEMM('N', 'N', 1.0_R8Ki, m%dUdy, m%dYdu, 1.0_R8Ki, m%G, ErrStat2, ErrMsg2)
   if (Failed()) return

   !----------------------------------------------------------------------------
   ! Assemble Jacobian
   !----------------------------------------------------------------------------

   ! Group (1,1)
   m%Jac(p%iJX(1):p%iJX(2), p%iJX(1):p%iJX(2)) = &
      -p%C(2)*m%dXdx(p%iX2(1):p%iX2(2), p%iX2(1):p%iX2(2)) - &
      p%C(3)*m%dXdx(p%iX2(1):p%iX2(2), p%iX1(1):p%iX1(2))
   do i = p%iJX(1), p%iJX(2)
      m%Jac(i, i) = m%Jac(i, i) + 1.0_R8Ki
   end do

   ! Group (2,1)
   m%Jac(p%iJUT(1):p%iJUT(2), p%iJX(1):p%iJX(2)) = &
      p%C(2)*matmul(m%dUdy(p%iUT(1):p%iUT(2), p%iyT(1):p%iyT(2)), &
                    m%dYdx(p%iyT(1):p%iyT(2), p%iX2(1):p%iX2(2))) + &
      p%C(3)*matmul(m%dUdy(p%iUT(1):p%iUT(2), p%iyT(1):p%iyT(2)), &
                    m%dYdx(p%iyT(1):p%iyT(2), p%iX1(1):p%iX1(2)))

   ! Group (1,2)
   m%Jac(p%iJX(1):p%iJX(2), p%iJUT(1):p%iJUT(2)) = -m%dXdu(p%iX2(1):p%iX2(2), p%iUT(1):p%iUT(2))

   ! Group (2,2) - Inputs
   m%Jac(p%iJU(1):p%iJU(2), p%iJU(1):p%iJU(2)) = m%G

   ! Condition jacobian matrix before factoring
   m%Jac(p%iJL, :) = m%Jac(p%iJL, :)/p%Scale_UJac
   m%Jac(:, p%iJL) = m%Jac(:, p%iJL)*p%Scale_UJac

   ! Write debug matrices if requested
   if (DebugJacobian) then
      call BuildJacobian_Debug(m, Turbine, ErrStat2, ErrMsg2)
      if (Failed()) return
   end if

   ! Factor jacobian matrix
   call LAPACK_getrf(size(m%Jac, 1), size(m%Jac, 1), m%Jac, m%IPIV, ErrStat2, ErrMsg2)
   if (Failed()) return

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine AddDeltaToStates(Mods, ModOrder, dx, x)
   type(ModDataType), intent(in)    :: Mods(:)  !< Module data
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
   type(ModDataType), intent(in)    :: Mods(:)  !< Module data
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
               ! Compose WM components (change in orientation with orientation)
               do k = 1, size(Var%iGblSol), 3
                  ind = Var%iGblSol(k:k + 2)
                  u(ind) = wm_compose(wm_from_xyz(du(ind)), u(ind))
               end do
            end select

         end associate
      end do
   end do
end subroutine

subroutine ComputeDeltaU(Mods, ModOrder, PosAry, NegAry, DiffAry)
   type(ModDataType), intent(in) :: Mods(:)        ! Array of modules
   integer(IntKi), intent(in)    :: ModOrder(:)    ! Array of module indices
   real(R8Ki), intent(in)        :: PosAry(:)      ! Positive result array
   real(R8Ki), intent(in)        :: NegAry(:)      ! Negative result array
   real(R8Ki), intent(inout)     :: DiffAry(:)     ! Array containing difference
   integer(IntKi)                :: i, j, k, ind(3)
   real(R8Ki)                    :: DeltaWM(3)

   ! Loop through module index order
   do i = 1, size(ModOrder)

      ! Loop through input variables in module
      do j = 1, size(Mods(ModOrder(i))%Vars%u)

         associate (Var => Mods(ModOrder(i))%Vars%u(j))

            if (.not. allocated(Var%iGblSol)) cycle

            ! If variable field is orientation
            if (Var%Field == VF_Orientation) then

               ! Loop through nodes
               do k = 1, Var%Nodes, 3

                  ! Get vector of indicies of WM rotation parameters in array
                  ind = Var%iGblSol(k:k + 2)

                  ! Compose WM parameters to go from negative to positive array
                  ! then store change in radians
                  DiffAry(ind) = wm_to_xyz(wm_compose(wm_inv(NegAry(ind)), PosAry(ind)))
               end do

            else

               ! Subtract negative array from positive array
               DiffAry(Var%iGblSol) = PosAry(Var%iGblSol) - NegAry(Var%iGblSol)
            end if
         end associate
      end do
   end do
end subroutine

subroutine PackModuleStates(ModData, this_state, T, x)
   type(ModDataType), intent(in)          :: ModData(:)  !< Module data
   integer(IntKi), intent(in)             :: this_state  !< State index
   type(FAST_TurbineType), intent(inout)  :: T           !< Turbine type
   real(R8Ki), intent(inout)              :: x(:)
   integer(IntKi)                         :: ii, j

   ! Must support all Tight Coupling modules
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

end subroutine

subroutine UnpackModuleStates(ModData, ModOrder, this_state, T, x)
   type(ModDataType), intent(in)          :: ModData(:)  !< Module data
   integer(IntKi), intent(in)             :: ModOrder(:)
   integer(IntKi), intent(in)             :: this_state  !< State index
   type(FAST_TurbineType), intent(inout)  :: T           !< Turbine type
   real(R8Ki), intent(inout)              :: x(:)
   integer(IntKi)                         :: j

   do j = 1, size(ModOrder)
      associate (Mod => ModData(ModOrder(j)))
         select case (Mod%ID)
         case (Module_ED)
            call XferGblToLoc1D(Mod%ixs, x, T%ED%m%Vals%x)
            call ED_UnpackStateValues(T%ED%p, T%ED%m%Vals%x, T%ED%x(this_state))
         case (Module_BD)
            call XferGblToLoc1D(Mod%ixs, x, T%BD%m(Mod%Ins)%Vals%x)
            call BD_UnpackStateValues(T%BD%p(Mod%Ins), T%BD%m(Mod%Ins)%Vals%x, T%BD%x(Mod%Ins, this_state))
         case (Module_SD)
            ! call SD_UnpackStateValues(SD%p, x(SD%p%Vars%ix), SD%x(this_state))
         end select
      end associate
   end do
end subroutine

! PackModuleInputs packs input values from Option 1 modules
subroutine PackModuleInputs(ModData, ModOrder, T, u)
   type(ModDataType), intent(in)          :: ModData(:)  !< Module data
   integer(IntKi), intent(in)             :: ModOrder(:)
   type(FAST_TurbineType), intent(inout)  :: T           !< Turbine type
   real(R8Ki), optional, intent(inout)    :: u(:)
   integer(IntKi)                         :: j

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

end subroutine

! PackModuleUs packs U values from Option 1 modules
subroutine PackModuleUs(ModData, ModOrder, T, u)
   type(ModDataType), intent(in)          :: ModData(:)  !< Module data
   integer(IntKi), intent(in)             :: ModOrder(:)
   type(FAST_TurbineType), intent(inout)  :: T           !< Turbine type
   real(R8Ki), intent(inout)              :: u(:)
   integer(IntKi)                         :: j

   do j = 1, size(ModOrder)
      associate (Mod => ModData(ModOrder(j)))
         select case (Mod%ID)
         case (Module_ED)
            call ED_PackInputValues(T%ED%p, T%ED%u, T%ED%m%Vals%u)
            call XferLocToGbl1D(Mod%ius, T%ED%m%Vals%u, u)
         case (Module_BD)
            call BD_PackInputValues(T%BD%p(Mod%Ins), T%BD%u(Mod%Ins), T%BD%m(Mod%Ins)%Vals%u)
            call XferLocToGbl1D(Mod%ius, T%BD%m(Mod%Ins)%Vals%u, u)
         case (Module_SD)
            ! call SD_PackInputValues(SD%p, SD%Input, SD%m%Vals%u)
            ! u(SD%p%Vars%iu) = SD%m%Vals%u
         end select
      end associate
   end do

end subroutine

! UnpackModuleInputs unpacks input values from Option 1 modules
subroutine UnpackModuleInputs(ModData, ModOrder, T, u)
   type(ModDataType), intent(in)          :: ModData(:)  !< Module data
   integer(IntKi), intent(in)             :: ModOrder(:) !< Module index array
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

!-------------------------------------------------------------------------------
! Debugging routines
!-------------------------------------------------------------------------------

subroutine Solver_Init_Debug(p, m, Mods)
   type(TC_ParameterType), intent(in)  :: p           !< Parameters
   type(TC_MiscVarType), intent(in)    :: m           !< Misc variables
   type(ModDataType), intent(in)       :: Mods(:)     !< Module data
   integer(IntKi)                      :: i, j

   write (DebugUn, '(A,*(I4))') " p%iJX2  = ", p%iJX
   write (DebugUn, '(A,*(I4))') " p%iJUT  = ", p%iJUT
   write (DebugUn, '(A,*(I4))') " p%iJU   = ", p%iJU
   write (DebugUn, '(A,*(I4))') " p%iJL   = ", p%iJL
   write (DebugUn, '(A,*(I4))') " p%iX2   = ", p%iX2
   write (DebugUn, '(A,*(I4))') " p%iX1   = ", p%iX1
   write (DebugUn, '(A,*(I4))') " p%iUT   = ", p%iUT
   write (DebugUn, '(A,*(I4))') " p%iU1   = ", p%iU1
   write (DebugUn, '(A,*(I4))') " p%iyT   = ", p%iyT
   write (DebugUn, '(A,*(I4))') " p%iy1   = ", p%iy1
   write (DebugUn, *) "shape(m%dYdx) = ", shape(m%dYdx)
   write (DebugUn, *) "shape(m%dYdu) = ", shape(m%dYdu)
   write (DebugUn, *) "shape(m%dXdx) = ", shape(m%dXdx)
   write (DebugUn, *) "shape(m%dXdu) = ", shape(m%dXdu)
   write (DebugUn, *) "shape(m%dUdu) = ", shape(m%dUdu)
   write (DebugUn, *) "shape(m%dUdy) = ", shape(m%dUdy)

   do i = 1, size(Mods)
      write (DebugUn, *) "Module   = ", Mods(i)%Abbr
      write (DebugUn, *) "ModuleID = ", Mods(i)%ID
      do j = 1, size(Mods(i)%Vars%x)
         if (.not. allocated(Mods(i)%Vars%x(j)%iGblSol)) cycle
         write (DebugUn, *) "Var = "//trim(Mods(i)%Abbr)//trim(Num2LStr(Mods(i)%Ins))//" X "//trim(Mods(i)%Vars%x(j)%Name)// &
            " ("//trim(MV_FieldString(Mods(i)%Vars%x(j)%Field))//")"
         write (DebugUn, '(A,*(I4))') "  X iLoc    = ", Mods(i)%Vars%x(j)%iLoc
         write (DebugUn, '(A,*(I4))') "  X iGblSol = ", Mods(i)%Vars%x(j)%iGblSol
      end do
      do j = 1, size(Mods(i)%Vars%u)
         if (.not. allocated(Mods(i)%Vars%u(j)%iGblSol)) cycle
         write (DebugUn, *) "Var = "//trim(Mods(i)%Abbr)//trim(Num2LStr(Mods(i)%Ins))//" U "//trim(Mods(i)%Vars%u(j)%Name)// &
            " ("//trim(MV_FieldString(Mods(i)%Vars%u(j)%Field))//")"
         write (DebugUn, '(A,*(I4))') "  U iLoc    = ", Mods(i)%Vars%u(j)%iLoc
         write (DebugUn, '(A,*(I4))') "  U iGblSol = ", Mods(i)%Vars%u(j)%iGblSol
      end do
      do j = 1, size(Mods(i)%Vars%y)
         if (.not. allocated(Mods(i)%Vars%y(j)%iGblSol)) cycle
         write (DebugUn, *) "Var = "//trim(Mods(i)%Abbr)//trim(Num2LStr(Mods(i)%Ins))//" Y "//trim(Mods(i)%Vars%y(j)%Name)// &
            " ("//trim(MV_FieldString(Mods(i)%Vars%y(j)%Field))//")"
         write (DebugUn, '(A,*(I4))') "  Y iLoc    = ", Mods(i)%Vars%y(j)%iLoc
         write (DebugUn, '(A,*(I4))') "  Y iGblSol = ", Mods(i)%Vars%y(j)%iGblSol
      end do
   end do
end subroutine

subroutine Solver_Step_Debug(p, m, step, iterCorr, iterConv, delta_norm)
   type(TC_ParameterType), intent(in)  :: p           !< Parameters
   type(TC_MiscVarType), intent(in)    :: m           !< Misc variables
   integer(IntKi), intent(in)          :: step, iterCorr, iterConv
   real(R8Ki), intent(in)              :: delta_norm

   write (DebugUn, *) "step = ", step
   write (DebugUn, *) "iterCorr = ", iterCorr
   write (DebugUn, *) "iterConv = ", iterConv
   write (DebugUn, '(A,*(ES16.7))') " du = ", m%du
   write (DebugUn, '(A,*(ES16.7))') " dx = ", m%dx
   write (DebugUn, '(A,*(ES16.7))') " u  = ", m%un
   write (DebugUn, '(A,*(ES16.7))') " u_tmp = ", m%u_tmp
   write (DebugUn, '(A,*(ES16.7))') " U  = ", m%UDiff
   write (DebugUn, '(A,*(ES16.7))') " x  = ", m%xn
   write (DebugUn, *) "delta_norm = ", delta_norm
end subroutine

subroutine BuildJacobian_Debug(m, T, ErrStat, ErrMsg)
   type(TC_MiscVarType), intent(in)    :: m        !< Misc variables
   type(FAST_TurbineType), intent(in)  :: T        !< Turbine
   integer(IntKi), intent(out)         :: ErrStat
   character(*), intent(out)           :: ErrMsg

   character(*), parameter    :: RoutineName = 'BuildJacobian_Debug'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2

   if (MatrixUn == -1) then
      call GetNewUnit(MatrixUn, ErrStat2, ErrMsg2); if (Failed()) return
   end if

   call DumpMatrix(MatrixUn, "dUdu.bin", m%dUdu, ErrStat2, ErrMsg2); if (Failed()) return
   call DumpMatrix(MatrixUn, "dUdy.bin", m%dUdy, ErrStat2, ErrMsg2); if (Failed()) return
   call DumpMatrix(MatrixUn, "dXdu.bin", m%dXdu, ErrStat2, ErrMsg2); if (Failed()) return
   call DumpMatrix(MatrixUn, "dXdx.bin", m%dXdx, ErrStat2, ErrMsg2); if (Failed()) return
   call DumpMatrix(MatrixUn, "dYdu.bin", m%dYdu, ErrStat2, ErrMsg2); if (Failed()) return
   call DumpMatrix(MatrixUn, "dYdx.bin", m%dYdx, ErrStat2, ErrMsg2); if (Failed()) return
   call DumpMatrix(MatrixUn, "ED-dXdu.bin", T%ED%m%Vals%dXdu, ErrStat2, ErrMsg2); if (Failed()) return
   call DumpMatrix(MatrixUn, "ED-dXdx.bin", T%ED%m%Vals%dXdx, ErrStat2, ErrMsg2); if (Failed()) return
   call DumpMatrix(MatrixUn, "ED-dYdu.bin", T%ED%m%Vals%dYdu, ErrStat2, ErrMsg2); if (Failed()) return
   call DumpMatrix(MatrixUn, "ED-dYdx.bin", T%ED%m%Vals%dYdx, ErrStat2, ErrMsg2); if (Failed()) return
   call DumpMatrix(MatrixUn, "G.bin", m%G, ErrStat2, ErrMsg2); if (Failed()) return
   call DumpMatrix(MatrixUn, "BD-dXdu.bin", T%BD%m(1)%Vals%dXdu, ErrStat2, ErrMsg2); if (Failed()) return
   call DumpMatrix(MatrixUn, "BD-dXdx.bin", T%BD%m(1)%Vals%dXdx, ErrStat2, ErrMsg2); if (Failed()) return
   call DumpMatrix(MatrixUn, "BD-dYdu.bin", T%BD%m(1)%Vals%dYdu, ErrStat2, ErrMsg2); if (Failed()) return
   call DumpMatrix(MatrixUn, "BD-dYdx.bin", T%BD%m(1)%Vals%dYdx, ErrStat2, ErrMsg2); if (Failed()) return
   call DumpMatrix(MatrixUn, "J.bin", m%Jac, ErrStat2, ErrMsg2); if (Failed()) return

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

end module
