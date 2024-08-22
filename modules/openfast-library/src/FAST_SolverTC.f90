module FAST_SolverTC

use NWTC_LAPACK
use FAST_Solver
use FAST_ModTypes
use FAST_Mapping
use FAST_ModGlue
use FAST_Funcs
use ElastoDyn
use BeamDyn
use SubDyn
use AeroDyn
use ServoDyn
use SC_DataEx

implicit none

private

! Public functions
public Solver_Init, Solver_Step0, Solver_Step

! Debugging
logical, parameter         :: DebugSolver = .false.
integer(IntKi)             :: DebugUn = -1
character(*), parameter    :: DebugFile = 'solver.dbg'
logical, parameter         :: DebugJacobian = .false.
integer(IntKi)             :: MatrixUn = -1

contains

subroutine Solver_Init(p_FAST, p, m, GlueModData, GlueModMaps, Turbine, ErrStat, ErrMsg)
   type(FAST_ParameterType), intent(in)      :: p_FAST         !< FAST parameters
   type(Glue_TCParam), intent(inout)         :: p              !< Glue Parameters
   type(Glue_TCMisc), intent(out)            :: m              !< Glue miscellaneous variables
   type(ModDataType), intent(inout)          :: GlueModData(:) !< Glue module data
   type(MappingType), intent(inout)          :: GlueModMaps(:) !< Module mappings
   type(FAST_TurbineType), intent(inout)     :: Turbine        !< all data for one instance of a turbine
   integer(IntKi), intent(out)               :: ErrStat        !< Error status of the operation
   character(*), intent(out)                 :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   character(*), parameter                   :: RoutineName = 'Solver_Init'
   integer(IntKi)                            :: ErrStat2
   character(ErrMsgLen)                      :: ErrMsg2
   integer(IntKi)                            :: i, j, k
   integer(IntKi), allocatable               :: modIDs(:), modInds(:), iMod(:)

   !----------------------------------------------------------------------------
   ! Initialize data in TC structure
   !----------------------------------------------------------------------------

   ! Solver time step
   p%h = p_FAST%DT

   ! Number of steps to advance between Jacobian updates
   if (p_FAST%DT_UJac/p_FAST%DT + 1 < huge(1_IntKi)) then
      p%NStep_UJac = ceiling(p_FAST%DT_UJac/p_FAST%DT, IntKi)
   else
      p%NStep_UJac = huge(1_IntKi)
   end if

   ! Number of convergence iterations between Jacobian updates
   ! TODO: read from input file
   p%NIter_UJac = 10000

   ! Generalized alpha damping coefficient
   ! TODO: read from input file
   p%RhoInf = 1.0_R8Ki

   ! Max number of convergence iterations
   ! TODO: read from input file
   p%MaxConvIter = 6

   ! Convergence tolerance
   ! TODO: read from input file
   p%ConvTol = 1.0e-4_R8Ki

   ! Jacobian conditioning
   p%Scale_UJac = p_FAST%UJacSclFact

   ! Generalized alpha integration constants
   p%AlphaM = (2.0_R8Ki*p%RhoInf - 1.0_R8Ki)/(p%RhoInf + 1.0_R8Ki)
   p%AlphaF = p%RhoInf/(p%RhoInf + 1.0_R8Ki)
   p%Gamma = 0.5_R8Ki - p%AlphaM + p%AlphaF
   p%Beta = (1.0_R8Ki - p%AlphaM + p%AlphaF)**2.0_R8Ki/4.0_R8Ki

   ! Precalculate some coefficients
   p%BetaPrime = p%h*p%h*p%Beta*(1.0_R8Ki - p%AlphaF)/(1.0_R8Ki - p%AlphaM)
   p%GammaPrime = p%h*p%Gamma*(1.0_R8Ki - p%AlphaF)/(1.0_R8Ki - p%AlphaM)

   !----------------------------------------------------------------------------
   ! Module ordering for solve
   !----------------------------------------------------------------------------

   ! Create array of indices for Mods array
   modInds = [(i, i=1, size(GlueModData))]

   ! Get array of module IDs
   modIDs = [(GlueModData(i)%ID, i=1, size(GlueModData))]

   ! Indices of all modules in Step 0 initialization order
   p%iModInit = [pack(modInds, ModIDs == Module_ED), &
                 pack(modInds, ModIDs == Module_BD), &
                 pack(modInds, ModIDs == Module_SD), &
                 pack(modInds, ModIDs == Module_IfW), &
                 pack(modInds, ModIDs == Module_AD), &
                 pack(modInds, ModIDs == Module_SrvD)] ! ServoDyn is last

   ! Indices of tight coupling modules
   p%iModTC = [pack(modInds, ModIDs == Module_ED), &
               pack(modInds, ModIDs == Module_BD), &
               pack(modInds, ModIDs == Module_SD)]

   ! Indices of Option 1 modules
   p%iModOpt1 = [pack(modInds, ModIDs == Module_ExtPtfm), &
                 pack(modInds, ModIDs == Module_HD), &
                 pack(modInds, ModIDs == Module_Orca)]

   ! Indices of Option 2 modules
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

   ! Indices of modules to perform InputSolves after the Option 1 solve
   p%iModPost = [pack(modInds, ModIDs == Module_SrvD), &
                 pack(modInds, ModIDs == Module_MD), &
                 pack(modInds, ModIDs == Module_ExtInfw)]

   !----------------------------------------------------------------------------
   ! Set solve flags and combine relevant modules into TC module
   !----------------------------------------------------------------------------

   ! Set VF_Solve flag on Jacobian variables use by the tight coupling solver
   call SetVarSolveFlags()

   ! Combination of TC and Option 1 module indices
   iMod = [p%iModTC, p%iModOpt1]

   ! Build tight coupling module using solve variables from TC and Option 1 modules
   call Glue_CombineModules(m%Mod, GlueModData, GlueModMaps, iMod, &
                            VF_Solve, .true., ErrStat2, ErrMsg2)
   if (Failed()) return

   !----------------------------------------------------------------------------
   ! Recalculate glue variable locations to simplify Jacobian construction
   !----------------------------------------------------------------------------

   call CalcVarGlobalIndices(p, m%Mod, p%NumQ, p%NumJ, ErrStat2, ErrMsg2)
   if (Failed()) return

   !----------------------------------------------------------------------------
   ! Allocate state, input, and output storage
   !----------------------------------------------------------------------------

   ! Calculated inputs array
   call AllocAry(m%uCalc, m%Mod%Vars%Nu, "m%uCalc", ErrStat2, ErrMsg2); if (Failed()) return

   ! Generalized alpha state arrays
   call AllocAry(m%State%q_prev, p%NumQ, "m%State%q_prev", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%State%x, p%NumQ, "m%State%q_delta", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%State%q, p%NumQ, "m%State%q", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%State%v, p%NumQ, "m%State%v", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%State%vd, p%NumQ, "m%State%vd", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%State%a, p%NumQ, "m%State%a", ErrStat2, ErrMsg2); if (Failed()) return
   m%State%q_prev = 0.0_R8Ki
   m%State%x = 0.0_R8Ki
   m%State%q = 0.0_R8Ki
   m%State%v = 0.0_R8Ki
   m%State%vd = 0.0_R8Ki
   m%State%a = 0.0_R8Ki

   !----------------------------------------------------------------------------
   ! Allocate solver Jacobian matrix and right hand side
   !----------------------------------------------------------------------------

   ! Allocate Jacobian matrix, RHS/X matrix, Pivot array
   call AllocAry(m%Mod%Lin%J, p%NumJ, p%NumJ, "m%J", ErrStat, ErrMsg); if (Failed()) return
   call AllocAry(m%T, p%NumJ, 3, "m%T", ErrStat, ErrMsg); if (Failed()) return
   call AllocAry(m%XB, p%NumJ, 1, "m%XB", ErrStat, ErrMsg); if (Failed()) return
   call AllocAry(m%IPIV, p%NumJ, "m%IPIV", ErrStat, ErrMsg); if (Failed()) return
   m%Mod%Lin%J = 0.0_R8Ki

   ! Initialize iterations and steps until Jacobian is calculated to zero
   ! so it is calculated in first step
   m%IterUntilUJac = 0
   m%StepsUntilUJac = 0

   !----------------------------------------------------------------------------
   ! Write debug info to file
   !----------------------------------------------------------------------------

   if (DebugSolver) then
      call GetNewUnit(DebugUn, ErrStat2, ErrMsg2); if (Failed()) return
      call OpenFOutFile(DebugUn, DebugFile, ErrStat2, ErrMsg2); if (Failed()) return
      call Solver_Init_Debug(p, m, GlueModData, GlueModMaps)
   end if

contains

   subroutine SetVarSolveFlags()
      ! Loop through tight coupling modules and add VF_Solve flag to
      do i = 1, size(p%iModTC)
         associate (ModData => GlueModData(p%iModTC(i)))
            do j = 1, size(ModData%Vars%x)
               call MV_SetFlags(ModData%Vars%x(j), VF_Solve)   ! Continuous state variables
            end do
         end associate
      end do

      ! Loop through module mappings
      do j = 1, size(GlueModMaps)
         associate (Mapping => GlueModMaps(j), &
                    SrcMod => GlueModData(GlueModMaps(j)%iModSrc), &
                    DstMod => GlueModData(GlueModMaps(j)%iModDst))

            ! Skip custom mapping types
            if (Mapping%MapType == Map_Variable .or. Mapping%MapType == Map_Custom) cycle

            ! Skip mappings where source and destination are not in tight coupling
            if (all(SrcMod%ID /= TC_Modules) .and. all(DstMod%ID /= TC_Modules)) cycle

            ! If source module is in tight coupling
            if (any(SrcMod%ID == TC_Modules)) then

               ! Set mapping flag on source variables
               do i = 1, size(SrcMod%Vars%y)
                  associate (Var => SrcMod%Vars%y(i))
                     if (MV_EqualDL(Mapping%SrcDL, Var%DL)) call MV_SetFlags(Var, VF_Solve)
                  end associate
               end do

               ! Set mapping flag on source displacement mesh variables
               if (Mapping%MapType == Map_LoadMesh) then
                  do i = 1, size(SrcMod%Vars%u)
                     associate (Var => SrcMod%Vars%u(i))
                        if (MV_EqualDL(Mapping%SrcDispDL, Var%DL)) call MV_SetFlags(Var, VF_Solve)
                     end associate
                  end do
               end if
            end if

            ! If source module is in option 1
            if (any(SrcMod%ID == O1_Modules)) then

               ! Set mapping flag on source variables
               do i = 1, size(SrcMod%Vars%y)
                  associate (Var => SrcMod%Vars%y(i))
                     if (.not. MV_EqualDL(Mapping%SrcDL, Var%DL)) cycle
                     select case (Var%Field)
                     case (FieldForce, FieldMoment)
                        call MV_SetFlags(Var, VF_Solve)
                     end select
                  end associate
               end do

               ! Set mapping flag on source displacement mesh variables
               if (Mapping%MapType == Map_LoadMesh) then
                  do i = 1, size(SrcMod%Vars%u)
                     associate (Var => SrcMod%Vars%u(i))
                        if (.not. MV_EqualDL(Mapping%SrcDispDL, Var%DL)) cycle
                        select case (Var%Field)
                        case (FieldForce, FieldMoment)
                           call MV_SetFlags(Var, VF_Solve)
                        end select
                     end associate
                  end do
               end if
            end if

            ! If destination module is in tight coupling
            if (any(DstMod%ID == TC_Modules) ) then

               ! Set mapping flag on destination variables
               do i = 1, size(DstMod%Vars%u)
                  associate (Var => DstMod%Vars%u(i))
                     if (MV_EqualDL(Mapping%DstDL, Var%DL)) call MV_SetFlags(Var, VF_Solve)
                  end associate
               end do

               ! Set mapping flag on destination displacement mesh variables
               if (Mapping%MapType == Map_LoadMesh) then
                  do i = 1, size(DstMod%Vars%y)
                     associate (Var => DstMod%Vars%y(i))
                        if (MV_EqualDL(Mapping%DstDispDL, Var%DL)) call MV_SetFlags(Var, VF_Solve)
                     end associate
                  end do
               end if
            end if

            ! If destination module is in option 1
            if (any(DstMod%ID == O1_Modules)) then

               ! Set mapping flag on destination variables
               do i = 1, size(DstMod%Vars%u)
                  associate (Var => DstMod%Vars%u(i))
                     if (.not. MV_EqualDL(Mapping%DstDL, Var%DL)) cycle
                     select case (Var%Field)
                     case (FieldTransAcc, FieldAngularAcc)
                        call MV_SetFlags(Var, VF_Solve)
                     end select
                  end associate
               end do

               ! Set mapping flag on destination displacement mesh variables
               if (Mapping%MapType == Map_LoadMesh) then
                  do i = 1, size(DstMod%Vars%y)
                     associate (Var => DstMod%Vars%y(i))
                        if (.not. MV_EqualDL(Mapping%DstDispDL, Var%DL)) cycle
                        select case (Var%Field)
                        case (FieldTransAcc, FieldAngularAcc)
                           call MV_SetFlags(Var, VF_Solve)
                        end select
                     end associate
                  end do
               end if
            end if
         end associate
      end do
   end subroutine

   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine CalcVarGlobalIndices(p, ModTC, NumQ, NumJ, ErrStat, ErrMsg)
   type(Glue_TCParam), intent(inout)   :: p        !< Parameters
   type(ModGlueType), intent(inout)    :: ModTC    !< Module data
   integer(IntKi), intent(out)         :: NumJ     !< Number of rows in Jacobian
   integer(IntKi), intent(out)         :: NumQ     !< Number of rows in state matrix
   integer(IntKi), intent(out)         :: ErrStat  !< Error status of the operation
   character(*), intent(out)           :: ErrMsg   !< Error message if ErrStat /= ErrID_None

   character(*), parameter             :: RoutineName = 'CalcVarGlobalIndices'
   integer(IntKi)                      :: ErrStat2    ! local error status
   character(ErrMsgLen)                :: ErrMsg2     ! local error message
   integer(IntKi)                      :: i, j, k, num, iGlu
   integer(IntKi)                      :: ix, iu, iy

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Initialize indices to zero
   p%iX1 = 0
   p%iX2 = 0
   p%iUT = 0
   p%iU1 = 0
   p%iUL = 0
   p%iyT = 0
   p%iy1 = 0
   p%iJX = 0
   p%iJU = 0
   p%iJUT = 0
   p%iJL = 0

   ! Loop through modules in data array and zero glue locations
   do i = 1, size(ModTC%ModData)
      associate (Vars => ModTC%ModData(i)%Vars)
         if (allocated(Vars%x)) then
            do j = 1, size(Vars%x)
               Vars%x(j)%iGlu = 0
            end do
         end if
         if (allocated(Vars%u)) then
            do j = 1, size(Vars%u)
               Vars%u(j)%iGlu = 0
            end do
         end if
         if (allocated(Vars%y)) then
            do j = 1, size(Vars%y)
               Vars%y(j)%iGlu = 0
            end do
         end if
      end associate
   end do

   !----------------------------------------------------------------------------
   ! Calculate TC state glue locations (displacements then velocities)
   !----------------------------------------------------------------------------

   ! Initialize glue index
   iGlu = 0

   ! Set indices for displacement variables
   do i = 1, size(ModTC%ModData)
      associate (Vars => ModTC%ModData(i)%Vars)
         if (.not. allocated(Vars%x)) cycle
         do j = 1, size(Vars%x)
            if (Vars%x(j)%DerivOrder == 0) then
               Vars%x(j)%iGlu = [iGlu + 1, iGlu + Vars%x(j)%Num]
               iGlu = Vars%x(j)%iGlu(2)
            end if
         end do
      end associate
   end do

   ! Start and end indices of displacement variables
   if (iGlu > 0) p%iX1 = [1, iGlu]

   ! Set indices for velocity variables
   do i = 1, size(ModTC%ModData)
      associate (Vars => ModTC%ModData(i)%Vars)
         if (.not. allocated(Vars%x)) cycle
         do j = 1, size(Vars%x)
            if (Vars%x(j)%DerivOrder == 1) then
               Vars%x(j)%iGlu = [iGlu + 1, iGlu + Vars%x(j)%Num]
               iGlu = Vars%x(j)%iGlu(2)
            end if
         end do
      end associate
   end do

   ! Start and end indices of velocity variables
   if (iGlu > p%iX1(2)) p%iX2 = [p%iX1(2) + 1, iGlu]

   !----------------------------------------------------------------------------
   ! Calculate input variable glue locations (group load and non-load)
   !----------------------------------------------------------------------------

   ! Initialize glue index
   iGlu = 0

   ! Set indices of Tight Coupling input variables (non-load)
   do i = 1, size(p%iModTC)
      associate (Vars => ModTC%ModData(i)%Vars)
         if (.not. allocated(Vars%u)) cycle
         do j = 1, size(Vars%u)
            if (.not. MV_IsLoad(Vars%u(j))) then
               Vars%u(j)%iGlu = [iGlu + 1, iGlu + Vars%u(j)%Num]
               iGlu = Vars%u(j)%iGlu(2)
            end if
         end do
      end associate
   end do

   ! Set start index of load values
   p%iUL(1) = iGlu + 1

   ! Set indices of Tight Coupling input variables (load)
   do i = 1, size(p%iModTC)
      associate (Vars => ModTC%ModData(i)%Vars)
         if (.not. allocated(Vars%u)) cycle
         do j = 1, size(Vars%u)
            if (MV_IsLoad(Vars%u(j))) then
               Vars%u(j)%iGlu = [iGlu + 1, iGlu + Vars%u(j)%Num]
               iGlu = Vars%u(j)%iGlu(2)
            end if
         end do
      end associate
   end do

   ! Set start/end indices for tight coupling inputs
   if (iGlu > 0) p%iUT = [1, iGlu]

   ! Set indices of Option 1 input variables (load)
   do i = size(p%iModTC) + 1, size(ModTC%ModData)
      associate (Vars => ModTC%ModData(i)%Vars)
         if (.not. allocated(Vars%u)) cycle
         do j = 1, size(Vars%u)
            if (MV_IsLoad(Vars%u(j))) then
               Vars%u(j)%iGlu = [iGlu + 1, iGlu + Vars%u(j)%Num]
               iGlu = Vars%u(j)%iGlu(2)
            end if
         end do
      end associate
   end do

   ! Set end index of load values
   if (iGlu >= p%iUL(1)) p%iUL(2) = iGlu

   ! Set indices of Option 1 input variables (non-load)
   do i = size(p%iModTC) + 1, size(ModTC%ModData)
      associate (Vars => ModTC%ModData(i)%Vars)
         if (.not. allocated(Vars%u)) cycle
         do j = 1, size(Vars%u)
            if (.not. MV_IsLoad(Vars%u(j))) then
               Vars%u(j)%iGlu = [iGlu + 1, iGlu + Vars%u(j)%Num]
               iGlu = Vars%u(j)%iGlu(2)
            end if
         end do
      end associate
   end do

   ! Set start/end indices for Option 1 inputs
   if (iGlu > p%iUT(2)) p%iU1 = [p%iUT(2) + 1, iGlu]

   !----------------------------------------------------------------------------
   ! Calculate output variable categories and indices
   !----------------------------------------------------------------------------

   ! Initialize glue index
   iGlu = 0

   ! Set indices of Tight Coupling output variables
   do i = 1, size(p%iModTC)
      associate (Vars => ModTC%ModData(i)%Vars)
         if (.not. allocated(Vars%y)) cycle
         do j = 1, size(Vars%y)
            Vars%y(j)%iGlu = [iGlu + 1, iGlu + Vars%y(j)%Num]
            iGlu = Vars%y(j)%iGlu(2)
         end do
      end associate
   end do

   ! Save number of tight coupling inputs
   if (iGlu > 0) p%iyT = [1, iGlu]

   ! Set indices of Option 1 output variables
   do i = size(p%iModTC) + 1, size(ModTC%ModData)
      associate (Vars => ModTC%ModData(i)%Vars)
         if (.not. allocated(Vars%y)) cycle
         do j = 1, size(Vars%y)
            Vars%y(j)%iGlu = [iGlu + 1, iGlu + Vars%y(j)%Num]
            iGlu = Vars%y(j)%iGlu(2)
         end do
      end associate
   end do

   ! Calculate number of option 1 outputs
   if (iGlu > p%iyT(2)) p%iy1 = [p%iyT(2) + 1, iGlu]

   !----------------------------------------------------------------------------
   ! Allocate q storage for generalized alpha algorithm
   ! This matrix stores equation state in an (N,4) array where:
   !  - N is the number of equations (rows)
   !  - Column 1 is position
   !  - Column 2 is velocity
   !  - Column 3 is acceleration
   !  - Column 4 is generalized alpha algorithmic acceleration
   !----------------------------------------------------------------------------

   ! Initialize number of q states (ignore derivatives)
   NumQ = 0

   ! Loop through tight coupling modules in glue module
   do i = 1, size(p%iModTC)

      associate (xVars => ModTC%ModData(i)%Vars%x)

         ! Loop through state variables
         do j = 1, size(xVars)

            ! Skip variables which already have a q index
            if (xVars(j)%iq(1) > 0) cycle

            ! Set q index for variable and update number
            xVars(j)%iq = [NumQ + 1, NumQ + xVars(j)%Num]
            NumQ = NumQ + xVars(j)%Num

            ! Loop through remaining vars  if the names match
            do k = j + 1, size(xVars)

               ! If names are different then they don't match, skip
               if (xVars(j)%Name /= xVars(k)%Name) cycle

               ! If field is not the same or a derivative of current field, skip
               select case (xVars(j)%Field)
               case (FieldTransDisp, FieldTransVel, FieldTransAcc)
                  if (all(xVars(k)%Field /= TransFields)) cycle
               case (FieldOrientation, FieldAngularDisp, FieldAngularVel, FieldAngularAcc)
                  if (all(xVars(k)%Field /= AngularFields)) cycle
               case (FieldForce, FieldMoment)
                  cycle
               end select

               ! Copy q row indices
               xVars(k)%iq = xVars(j)%iq

            end do
         end do
      end associate
   end do

   !----------------------------------------------------------------------------
   ! Populate combined variable arrays
   !----------------------------------------------------------------------------

   ix = 0; iu = 0; iy = 0
   do i = 1, size(ModTC%ModData)
      associate (ModData => ModTC%ModData(i))

         if (allocated(ModData%Vars%x)) then
            do j = 1, size(ModData%Vars%x)
               ix = ix + 1
               ModTC%Vars%x(ix)%iLoc = ModData%Vars%x(j)%iGlu
               ModTC%Vars%x(ix)%iq = ModData%Vars%x(j)%iq
            end do
         end if

         if (allocated(ModData%Vars%u)) then
            do j = 1, size(ModData%Vars%u)
               iu = iu + 1
               ModTC%Vars%u(iu)%iLoc = ModData%Vars%u(j)%iGlu
            end do
         end if

         if (allocated(ModData%Vars%y)) then
            do j = 1, size(ModData%Vars%y)
               iy = iy + 1
               ModTC%Vars%y(iy)%iLoc = ModData%Vars%y(j)%iGlu
            end do
         end if

      end associate
   end do

   !----------------------------------------------------------------------------
   ! Jacobian indices and ranges
   !----------------------------------------------------------------------------

   ! Calculate size of Jacobian matrix
   NumJ = NumQ + ModTC%Vars%Nu

   ! Get start and end indices for state part of Jacobian
   if (NumQ > 0) p%iJX = [1, NumQ]

   ! Get start and end indices for tight coupling input part of Jacobian
   if (p%iUT(1) > 0) p%iJUT = NumQ + p%iUT

   ! Get start and end indices for input part of Jacobian
   if (p%iUT(1) > 0 .or. p%iU1(2) > 0) p%iJU = NumQ + [1, max(p%iUT(2), p%iU1(2))]

   ! Get Jacobian indices containing loads
   if (p%iUL(1) > 0) p%iJL = NumQ + p%iUL

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine Solver_Step0(p, m, GlueModData, GlueModMaps, Turbine, ErrStat, ErrMsg)
   type(Glue_TCParam), intent(in)            :: p              !< Parameters
   type(Glue_TCMisc), intent(inout)          :: m              !< Misc variables
   type(ModDataType), intent(inout)          :: GlueModData(:) !< Glue module data
   type(MappingType), intent(inout)          :: GlueModMaps(:) !< Module mappings
   type(FAST_TurbineType), intent(inout)     :: Turbine        !< Turbine type
   integer(IntKi), intent(out)               :: ErrStat
   character(*), intent(out)                 :: ErrMsg

   character(*), parameter    :: RoutineName = 'Solver_Step0'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   logical, parameter         :: IsSolve = .true.
   integer(IntKi)             :: i, j, k
   real(R8Ki)                 :: ConvError
   logical                    :: converged
   integer(IntKi), parameter  :: n_t_global = -1     ! loop counter
   integer(IntKi), parameter  :: n_t_global_next = 0 ! loop counter
   real(DbKi)                 :: t_initial           ! next simulation time
   real(DbKi)                 :: t_global_next       ! next simulation time

   ErrStat = ErrID_None
   ErrMsg = ''

   !----------------------------------------------------------------------------
   ! Miscellaneous initial step setup
   !----------------------------------------------------------------------------

   t_initial = Turbine%m_FAST%t_global
   t_global_next = t_initial + n_t_global_next*p%h
   Turbine%y_FAST%WriteThisStep = NeedWriteOutput(n_t_global_next, t_initial, Turbine%p_FAST%TStart, Turbine%p_FAST%n_DT_Out)

   ! Set flag to warn about convergence errors
   m%ConvWarn = .true.

   !----------------------------------------------------------------------------
   ! Calculate initial accelerations
   !----------------------------------------------------------------------------

   ! Transfer initial state from modules to solver
   do i = 1, size(m%Mod%ModData)
      associate (ModData => m%Mod%ModData(i))

         ! Get continuous state operating points
         call FAST_GetOP(ModData, t_initial, INPUT_CURR, STATE_CURR, Turbine, ErrStat2, ErrMsg2, &
                         x_op=ModData%Lin%x, x_glue=m%Mod%Lin%x)
         if (Failed()) return

         ! Transfer initial module state to GA state
         call TransferXtoQ(ModData, m%Mod%Lin%x, m%State)
      end associate
   end do

   ! Reset mapping ready for transfer flag
   GlueModMaps%Ready = .false.

   ! Loop until initial accelerations are converged, or max iterations are reached.
   ! TODO: may need a separate variable for max initial acceleration convergence iterations
   converged = .false.
   k = 0
   do while ((.not. converged) .and. (k <= 2*p%MaxConvIter))

      ! Transfer inputs and calculate outputs for all modules (use current state)
      do i = 1, size(p%iModInit)
         call FAST_InputSolve(p%iModInit(i), GlueModData, GlueModMaps, INPUT_CURR, &
                              Turbine, ErrStat2, ErrMsg2)
         if (Failed()) return
         call FAST_CalcOutput(GlueModData(p%iModInit(i)), GlueModMaps, t_initial, INPUT_CURR, STATE_CURR, &
                              Turbine, ErrStat2, ErrMsg2)
         if (Failed()) return
      end do

      ! Calculate continuous state derivatives for tight coupling modules (use current state)
      do i = 1, size(m%Mod%ModData)
         call FAST_GetOP(m%Mod%ModData(i), t_initial, INPUT_CURR, STATE_CURR, &
                         Turbine, ErrStat2, ErrMsg2, &
                         dx_op=m%Mod%ModData(i)%Lin%dx, dx_glue=m%Mod%Lin%dx)
         if (Failed()) return
      end do

      ! Copy acceleration (derivative of velocity) to acceleration array.
      associate (Vars => m%Mod%Vars)
         do i = 1, size(Vars%x)
            select case (Vars%x(i)%DerivOrder)
            case (1) ! Velocity
               m%State%vd(Vars%x(i)%iq(1):Vars%x(i)%iq(2)) = m%Mod%Lin%dx(Vars%x(i)%iLoc(1):Vars%x(i)%iLoc(2))
            end select
         end do
      end associate

      ! Calculate convergence error as L2 norm of diff between current and new accelerations
      ConvError = TwoNorm(m%State%vd - m%State%a)

      ! Update algorithmic acceleration
      m%State%a = m%State%vd

      ! Transfer acceleration to BeamDyn if module is active
      do i = 1, size(p%iModTC)
         associate (ModData => m%Mod%ModData(i))
            if (ModData%ID == Module_BD) then
               call SetBDAccel(ModData, m%State, Turbine%BD%OtherSt(ModData%Ins, STATE_CURR))
            end if
         end associate
      end do

      ! If difference is less than convergence tolerance, set flag to exit loop
      if ((k > 1) .and. (ConvError < p%ConvTol)) converged = .true.

      ! Increment iteration counter
      k = k + 1
   end do

   ! Print warning if not converged
   if (.not. converged) then
      call WrScr("Solver: initial accel not converged, diff="// &
                 trim(Num2LStr(ConvError))//", tol="//trim(Num2LStr(p%ConvTol)))
   end if

   ! Copy current state to previous state
   call Glue_CopyTC_State(m%State, m%StatePrev, MESH_NEWCOPY, ErrStat2, ErrMsg2)
   if (Failed()) return

   ! Initialize IO and states for all modules (also copies STATE_CURR to STATE_PRED)
   call FAST_InitIO(GlueModData, t_initial, p%h, Turbine, ErrStat, ErrMsg)
   if (ErrStat >= AbortErrLev) return

   ! Reset the Remap flags for all modules
   call FAST_ResetRemapFlags(GlueModData, GlueModMaps, Turbine, ErrStat, ErrMsg)
   if (ErrStat >= AbortErrLev) return

   !----------------------------------------------------------------------------
   ! Set Outputs
   !----------------------------------------------------------------------------

   Turbine%y_FAST%DriverWriteOutput(1) = real(k, ReKi)         ! ConvIter
   Turbine%y_FAST%DriverWriteOutput(2) = real(ConvError, ReKi) ! ConvError
   Turbine%y_FAST%DriverWriteOutput(3) = 0.0_ReKi              ! NumUJac

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine Solver_Step(n_t_global, t_initial, p, m, GlueModData, GlueModMaps, Turbine, ErrStat, ErrMsg)
   integer(IntKi), intent(in)                :: n_t_global     !< global time step
   real(DbKi), intent(in)                    :: t_initial      !< Initial simulation time
   type(Glue_TCParam), intent(in)            :: p              !< Parameters
   type(Glue_TCMisc), intent(inout)          :: m              !< Misc variables
   type(ModDataType), intent(inout)          :: GlueModData(:) !< Glue module data
   type(MappingType), intent(inout)          :: GlueModMaps(:) !< Module mappings
   type(FAST_TurbineType), intent(inout)     :: Turbine        !< Turbine type
   integer(IntKi), intent(out)               :: ErrStat
   character(*), intent(out)                 :: ErrMsg

   character(*), parameter    :: RoutineName = 'Solver_Step'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   logical, parameter         :: IsSolve = .true.
   integer(IntKi)             :: iterConv, iterCorr, iterTotal
   integer(IntKi)             :: NumUJac, NumCorrections
   real(R8Ki)                 :: delta_norm
   real(DbKi)                 :: t_global_next     ! next simulation time (m_FAST%t_global + p_FAST%dt)
   integer(IntKi)             :: n_t_global_next   ! n_t_global + 1
   integer(IntKi)             :: i, j, k
   integer(IntKi)             :: iMod
   logical                    :: ConvUJac          ! Jacobian updated for convergence
   real(R8Ki)                 :: RotDiff(3, 3)

   ErrStat = ErrID_None
   ErrMsg = ''

   !----------------------------------------------------------------------------
   ! Miscellaneous step updates
   !----------------------------------------------------------------------------

   ! Calculate the next global time step number and time
   n_t_global_next = n_t_global + 1
   t_global_next = t_initial + n_t_global_next*p%h

   ! Determine if output should be written in this step
   Turbine%y_FAST%WriteThisStep = NeedWriteOutput(n_t_global_next, t_global_next, Turbine%p_FAST%TStart, Turbine%p_FAST%n_DT_Out)

   ! Decrement number of time steps before updating the Jacobian
   m%StepsUntilUJac = m%StepsUntilUJac - 1

   ! Set Jacobian updated for convergence flag to false
   ConvUJac = .false.

   ! Init counters for number of Jacobian updates and number of convergence iterations
   NumUJac = 0
   iterTotal = 0

   !----------------------------------------------------------------------------
   ! Extrapolate/interpolate inputs for all modules
   !----------------------------------------------------------------------------

   ! Loop through all modules and extrapolate inputs
   do i = 1, size(GlueModData)
      associate (ModData => GlueModData(i))
         call FAST_ExtrapInterp(ModData, t_global_next, Turbine, ErrStat2, ErrMsg2)
         if (Failed()) return
      end associate
   end do

   !----------------------------------------------------------------------------
   ! Reset BD States
   !----------------------------------------------------------------------------

   ! Perform additional state manipulation on a per-module basis
   do i = 1, size(p%iModTC)
      associate (ModData => m%Mod%ModData(i))
         select case (ModData%ID)
         case (Module_ED)

            ! Update the azimuth angle
            call ED_UpdateAzimuth(Turbine%ED%p, Turbine%ED%x(STATE_CURR), ModData%DT)

         case (Module_BD)

            ! Transfer acceleration from TC state to BeamDyn
            call SetBDAccel(ModData, m%StatePrev, Turbine%BD%OtherSt(ModData%Ins, STATE_CURR))

            ! Reset BeamDyn states so they are relative to the root node
            call BD_UpdateGlobalRef(Turbine%BD%Input(INPUT_CURR, ModData%Ins), &
                                    Turbine%BD%p(ModData%Ins), &
                                    Turbine%BD%x(ModData%Ins, STATE_CURR), &
                                    Turbine%BD%OtherSt(ModData%Ins, STATE_CURR), &
                                    ErrStat2, ErrMsg2)
            if (Failed()) return

            ! Transfer acceleration from BeamDyn to state
            call GetBDAccel(ModData, Turbine%BD%OtherSt(ModData%Ins, STATE_CURR), m%StatePrev)

         case default
            cycle
         end select

         ! Collect updated states
         call FAST_GetOP(ModData, t_global_next, INPUT_CURR, STATE_CURR, Turbine, ErrStat2, ErrMsg2, &
                         x_op=ModData%Lin%x, x_glue=m%Mod%Lin%x)
         if (Failed()) return

         ! Transfer current states to linearization array
         call TransferXtoQ(ModData, m%Mod%Lin%x, m%StatePrev)
      end associate
   end do

   !----------------------------------------------------------------------------
   ! Prediction - guess solution state variables at end of time step
   !----------------------------------------------------------------------------

   call PredictNextState(p, m%StatePrev, m%Mod%Vars)

   !----------------------------------------------------------------------------
   ! Correction Iterations
   !----------------------------------------------------------------------------

   ! Loop through correction iterations
   iterCorr = 0
   NumCorrections = p%NumCrctn
   do while (iterCorr <= NumCorrections)

      ! Reset mapping ready flags
      GlueModMaps%Ready = .false.

      ! Copy TC solver states from previous to current
      call Glue_CopyTC_State(m%StatePrev, m%State, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
      if (Failed()) return

      ! Loop through tight coupling modules
      do i = 1, size(p%iModTC)
         associate (ModData => m%Mod%ModData(i))

            ! Transfer current states to linearization array
            call TransferQtoX(ModData, m%State, m%Mod%Lin%x)

            ! Copy state from current to predicted
            call FAST_CopyStates(ModData, Turbine, STATE_CURR, STATE_PRED, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
            if (Failed()) return

            ! Transfer solver states to module
            call FAST_SetOP(ModData, INPUT_CURR, STATE_PRED, Turbine, ErrStat2, ErrMsg2, &
                            x_op=ModData%Lin%x, x_glue=m%Mod%Lin%x)
            if (Failed()) return
         end associate
      end do

      !-------------------------------------------------------------------------
      ! Option 2 Solve
      !-------------------------------------------------------------------------

      ! Loop through Option 2 modules
      do i = 1, size(p%iModOpt2)
         associate (ModData => GlueModData(p%iModOpt2(i)))

            ! Solve for inputs
            call FAST_InputSolve(p%iModOpt2(i), GlueModData, GlueModMaps, INPUT_CURR, Turbine, ErrStat2, ErrMsg2)
            if (Failed()) return

            ! Update states
            call FAST_UpdateStates(ModData, t_initial, n_t_global, Turbine, ErrStat2, ErrMsg2)
            if (Failed()) return

            ! Calculate outputs
            call FAST_CalcOutput(ModData, GlueModMaps, t_global_next, INPUT_CURR, STATE_PRED, Turbine, ErrStat2, ErrMsg2)
            if (Failed()) return
         end associate
      end do

      ! -4.9345743494052119       -9.4942092091626126E-003   91.498545759647968
      ! -4.9345900476124651       -1.8983439772237952E-002   91.498449503167748
      ! write (*, *) Turbine%BD%y(1)%BldMotion%Position(:, 1) + Turbine%BD%y(1)%BldMotion%TranslationDisp(:, 1)

      !-------------------------------------------------------------------------
      ! Option 1 Solve
      !-------------------------------------------------------------------------

      ! Get inputs and update states for Option 1 modules
      do i = 1, size(p%iModOpt1)
         call FAST_InputSolve(p%iModOpt1(i), GlueModData, GlueModMaps, INPUT_CURR, Turbine, ErrStat2, ErrMsg2)
         if (Failed()) return
         call FAST_UpdateStates(GlueModData(p%iModOpt1(i)), t_initial, n_t_global, Turbine, ErrStat2, ErrMsg2)
         if (Failed()) return
      end do

      !-------------------------------------------------------------------------
      ! Pack inputs and modify states
      !-------------------------------------------------------------------------

      ! Pack TC and Option 1 inputs into u array
      do i = 1, size(m%Mod%ModData)
         associate (ModData => m%Mod%ModData(i))
            call FAST_GetOP(ModData, t_global_next, INPUT_CURR, STATE_PRED, Turbine, ErrStat2, ErrMsg2, &
                            u_op=ModData%Lin%u, u_glue=m%Mod%Lin%u)
            if (Failed()) return
         end associate
      end do

      !-------------------------------------------------------------------------
      ! Convergence Iterations
      !-------------------------------------------------------------------------

      ! Loop through convergence iterations
      do iterConv = 0, p%MaxConvIter

         ! Increment total number of convergence iterations in step
         iterTotal = iterTotal + 1

         ! Decrement number of iterations before updating the Jacobian
         m%IterUntilUJac = m%IterUntilUJac - 1

         !----------------------------------------------------------------------
         ! Calculate outputs for TC & Opt1 modules
         !----------------------------------------------------------------------

         do i = 1, size(m%Mod%ModData)
            associate (ModData => m%Mod%ModData(i))
               call FAST_CalcOutput(ModData, GlueModMaps, t_global_next, INPUT_CURR, STATE_PRED, &
                                    Turbine, ErrStat2, ErrMsg2)
               if (Failed()) return
            end associate
         end do

         !----------------------------------------------------------------------
         ! Convergence iteration check
         !----------------------------------------------------------------------

         ! If convergence iteration has reached or exceeded limit
         if (iterConv >= p%MaxConvIter) then

            ! If Jacobian has not been updated for convergence
            if (.not. ConvUJac) then

               ! Set counter to trigger a Jacobian update on next convergence iteration
               m%IterUntilUJac = 0

               ! If at the maximum number of correction iterations,
               ! increase limit to retry the step after the Jacobian is updated
               if (iterCorr == NumCorrections) NumCorrections = NumCorrections + 1

               ! Set flag indicating that the jacobian has been updated for convergence
               ConvUJac = .true.

            else

               ! Otherwise, correction iteration with Jacobian update has been tried,
               ! display warning that convergence failed and move to next step
               call SetErrStat(ErrID_Warn, "Failed to converge in "//trim(Num2LStr(p%MaxConvIter))// &
                               " iterations on step "//trim(Num2LStr(n_t_global_next))// &
                               " (error="//trim(Num2LStr(delta_norm))// &
                               ", tolerance="//trim(Num2LStr(p%ConvTol))//").", &
                               ErrStat, ErrMsg, RoutineName)
            end if

            ! Exit convergence loop to next correction iteration or next step
            exit
         end if

         !----------------------------------------------------------------------
         ! Update Jacobian
         !----------------------------------------------------------------------

         ! If number of iterations or steps until Jacobian is to be updated
         ! is zero or less, or first solution step, then rebuild the Jacobian.
         ! Note: BuildJacobian resets these counters.
         if ((m%IterUntilUJac <= 0) .or. &
             (m%StepsUntilUJac <= 0) .or. &
             (n_t_global_next == 1)) then
            NumUJac = NumUJac + 1
            call BuildJacobian(p, m, GlueModMaps, t_global_next, Turbine, ErrStat2, ErrMsg2)
            if (Failed()) return
         end if

         !----------------------------------------------------------------------
         ! Formulate right hand side (X_2^tight, U^tight, U^Option1)
         !----------------------------------------------------------------------

         ! Calculate continuous state derivatives for tight coupling modules
         do i = 1, size(m%Mod%ModData)
            call FAST_GetOP(m%Mod%ModData(i), t_global_next, INPUT_CURR, STATE_PRED, Turbine, ErrStat2, ErrMsg2, &
                            dx_op=m%Mod%ModData(i)%Lin%dx, dx_glue=m%Mod%Lin%dx)
            if (Failed()) return
         end do

         ! Input solve for tight coupling modules
         do i = 1, size(p%iModTC)
            associate (ModData => GlueModData(p%iModTC(i)))
               call FAST_InputSolve(p%iModTC(i), GlueModData, GlueModMaps, INPUT_TEMP, Turbine, ErrStat2, ErrMsg2)
               if (Failed()) return
            end associate
         end do

         ! Input solve for Option 1 modules
         do i = 1, size(p%iModOpt1)
            associate (ModData => GlueModData(p%iModOpt1(i)))
               call FAST_InputSolve(p%iModOpt1(i), GlueModData, GlueModMaps, INPUT_TEMP, Turbine, ErrStat2, ErrMsg2)
               if (Failed()) return
            end associate
         end do

         ! Transfer collect inputs into uCalc
         do i = 1, size(m%Mod%ModData)
            call FAST_GetOP(m%Mod%ModData(i), t_global_next, INPUT_TEMP, STATE_PRED, Turbine, ErrStat2, ErrMsg2, &
                            u_op=m%Mod%ModData(i)%Lin%u, u_glue=m%uCalc)
            if (Failed()) return
         end do

         !----------------------------------------------------------------------
         ! Populate residual vector and apply conditioning to loads
         !----------------------------------------------------------------------

         ! Calculate difference between calculated and predicted accelerations
         if (p%iJX(1) > 0) m%XB(p%iJX(1):p%iJX(2), 1) = m%Mod%Lin%dx(p%iX2(1):p%iX2(2)) - m%State%vd

         ! Calculate difference in U for all Option 1 modules (un - u_tmp)
         ! and add to RHS for TC and Option 1 modules
         if (p%iJU(1) > 0) call MV_ComputeDiff(m%Mod%Vars%u, m%uCalc, m%Mod%Lin%u, m%XB(p%iJU(1):p%iJU(2), 1))

         ! Apply conditioning factor to loads in RHS
         do i = p%iJL(1), p%iJL(2)
            m%XB(i, 1) = m%XB(i, 1)/p%Scale_UJac
         end do

         !----------------------------------------------------------------------
         ! Solve for state and input perturbations
         !----------------------------------------------------------------------

         ! Solve Jacobian and RHS
         call LAPACK_getrs('N', size(m%Mod%Lin%J, 1), m%Mod%Lin%J, m%IPIV, m%XB, ErrStat2, ErrMsg2)
         if (Failed()) return

         !----------------------------------------------------------------------
         ! Check perturbations for convergence and exit if below tolerance
         !----------------------------------------------------------------------

         ! Calculate average L2 norm of change in states and inputs
         delta_norm = TwoNorm(m%XB(:, 1))/size(m%XB)

         ! Write step debug info if requested
         if (DebugSolver) call Solver_Step_Debug(p, m, n_t_global_next, iterCorr, iterConv, delta_norm)

         ! If at least one convergence iteration has been done and
         ! the RHS norm is less than convergence tolerance, exit loop
         if ((iterConv > 0) .and. (delta_norm < p%ConvTol)) exit

         ! Remove conditioning
         do i = p%iJL(1), p%iJL(2)
            m%XB(i, 1) = m%XB(i, 1)*p%Scale_UJac
         end do

         !----------------------------------------------------------------------
         ! Update State for Tight Coupling modules
         !----------------------------------------------------------------------

         if (p%iJX(1) > 0) call UpdateStatePrediction(p, m%Mod%Vars, m%XB(p%iJX(1):p%iJX(2), 1), m%State)

         !----------------------------------------------------------------------
         ! Update inputs for Tight Coupling and Option 1 modules
         !----------------------------------------------------------------------

         ! Add change in inputs
         if (p%iJU(1) > 0) call MV_AddDelta(m%Mod%Vars%u, m%XB(p%iJU(1):p%iJU(2), 1), m%Mod%Lin%u)

         !----------------------------------------------------------------------
         ! Transfer updated TC and Option 1 states and inputs to modules
         !----------------------------------------------------------------------

         do i = 1, size(m%Mod%ModData)
            associate (ModData => m%Mod%ModData(i))

               ! Transfer States to linearization array
               call TransferQtoX(ModData, m%State, m%Mod%Lin%x)

               ! Transfer states and inputs to modules
               call FAST_SetOP(ModData, INPUT_CURR, STATE_PRED, Turbine, ErrStat2, ErrMsg2, &
                               x_op=ModData%Lin%x, x_glue=m%Mod%Lin%x, &
                               u_op=ModData%Lin%u, u_glue=m%Mod%Lin%u)
               if (Failed()) return
            end associate
         end do

         ! Transfer acceleration to BD other states
         do i = 1, size(p%iModTC)
            if (m%Mod%ModData(i)%ID == Module_BD) then
               call SetBDAccel(m%Mod%ModData(i), m%State, Turbine%BD%OtherSt(m%Mod%ModData(i)%Ins, STATE_PRED))
            end if
         end do
      end do

      ! Increment correction iteration counter
      iterCorr = iterCorr + 1

      ! Perform input solve for modules post Option 1 convergence
      do i = 1, size(p%iModPost)
         call FAST_InputSolve(p%iModPost(i), GlueModData, GlueModMaps, INPUT_CURR, Turbine, ErrStat2, ErrMsg2)
         if (Failed()) return
      end do

      ! Reset mesh remap
      call FAST_ResetRemapFlags(GlueModData, GlueModMaps, Turbine, ErrStat2, ErrMsg2)
      if (Failed()) return
   end do

   !----------------------------------------------------------------------------
   ! Update states for next step
   !----------------------------------------------------------------------------

   ! Update algorithmic acceleration
   m%State%a = m%State%a + (1.0_R8Ki - p%AlphaF)/(1.0_R8Ki - p%AlphaM)*m%State%vd

   ! Copy current state to previous state
   call Glue_CopyTC_State(m%State, m%StatePrev, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
   if (Failed()) return

   ! Copy the final predicted states from step t_global_next to actual states for that step
   do i = 1, size(GlueModData)
      call FAST_SaveStates(GlueModData(i), Turbine, ErrStat2, ErrMsg2)
      if (Failed()) return
   end do

   !----------------------------------------------------------------------------
   ! Set Outputs
   !----------------------------------------------------------------------------

   Turbine%y_FAST%DriverWriteOutput(1) = real(iterTotal, ReKi)    ! ConvIter
   Turbine%y_FAST%DriverWriteOutput(2) = real(delta_norm, ReKi)   ! ConvError
   Turbine%y_FAST%DriverWriteOutput(3) = real(NumUJac, ReKi)      ! NumUJac

   !----------------------------------------------------------------------------
   ! Update the global time
   !----------------------------------------------------------------------------

   Turbine%m_FAST%t_global = t_global_next

contains
   logical function Failed()
      if (ErrStat2 /= ErrID_None) call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine BuildJacobian(p, m, GlueModMaps, ThisTime, Turbine, ErrStat, ErrMsg)
   type(Glue_TCParam), intent(in)         :: p              !< Parameters
   type(Glue_TCMisc), intent(inout)       :: m              !< Misc variables
   type(MappingType), intent(inout)       :: GlueModMaps(:) !< Module mappings at glue level
   real(DbKi), intent(in)                 :: ThisTime       !< Time
   type(FAST_TurbineType), intent(inout)  :: Turbine        !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter                :: RoutineName = 'BuildJacobian'
   integer(IntKi)                         :: ErrStat2
   character(ErrMsgLen)                   :: ErrMsg2
   real(R8Ki)                             :: phi, rv(3), T(3, 3), tmp1, tmp2, T2(3, 3)
   integer(IntKi)                         :: i, j, k, idx

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Reset Jacobian update countdown values
   m%IterUntilUJac = p%NIter_UJac
   m%StepsUntilUJac = p%NStep_UJac

   if (size(m%Mod%Lin%J) == 0) return

   !----------------------------------------------------------------------------
   ! Get module Jacobians and assemble
   ! A: rows = x; columns = x (dXdx)
   ! B: rows = x; columns = u (dXdu)
   ! C: rows = y; columns = x (dYdx)
   ! D: rows = y; columns = u (dYdu)
   !----------------------------------------------------------------------------

   ! Initialize Jacobian matrices
   if (allocated(m%Mod%Lin%dYdx)) m%Mod%Lin%dYdx = 0.0_R8Ki
   if (allocated(m%Mod%Lin%dXdx)) m%Mod%Lin%dXdx = 0.0_R8Ki
   if (allocated(m%Mod%Lin%dXdu)) m%Mod%Lin%dXdu = 0.0_R8Ki
   if (allocated(m%Mod%Lin%dYdu)) m%Mod%Lin%dYdu = 0.0_R8Ki
   if (allocated(m%Mod%Lin%dUdy)) m%Mod%Lin%dUdy = 0.0_R8Ki
   if (allocated(m%Mod%Lin%dUdu)) then
      call Eye2D(m%Mod%Lin%dUdu, ErrStat2, ErrMsg2)
      if (Failed()) return
   end if

   ! Loop through modules tight coupling modules
   do i = 1, size(p%iModTC)
      associate (ModData => m%Mod%ModData(i))

         ! Calculate dYdx, dXdx for tight coupling modules
         call FAST_JacobianPContState(ModData, ThisTime, INPUT_CURR, STATE_PRED, Turbine, ErrStat2, ErrMsg2, &
                                      dXdx=ModData%Lin%dXdx, dXdx_glue=m%Mod%Lin%dXdx, &
                                      dYdx=ModData%Lin%dYdx, dYdx_glue=m%Mod%Lin%dYdx)
         if (Failed()) return

         ! Calculate Jacobians wrt inputs
         call FAST_JacobianPInput(ModData, ThisTime, INPUT_CURR, STATE_PRED, Turbine, ErrStat2, ErrMsg2, &
                                  dXdu=ModData%Lin%dXdu, dXdu_glue=m%Mod%Lin%dXdu, &
                                  dYdu=ModData%Lin%dYdu, dYdu_glue=m%Mod%Lin%dYdu)
         if (Failed()) return
      end associate
   end do

   ! Loop through Option 1 modules and calculate dYdu
   do i = size(p%iModTC) + 1, size(m%Mod%ModData)
      associate (ModData => m%Mod%ModData(i))
         call FAST_JacobianPInput(ModData, ThisTime, INPUT_CURR, STATE_PRED, Turbine, ErrStat2, ErrMsg2, &
                                  dYdu=ModData%Lin%dYdu, dYdu_glue=m%Mod%Lin%dYdu)
         if (Failed()) return
      end associate
   end do

   ! Calculate dUdu and dUdy for TC and Option 1 modules
   if (allocated(m%Mod%Lin%dUdy) .and. allocated(m%Mod%Lin%dUdu)) then
      call FAST_LinearizeMappings(m%Mod, GlueModMaps, Turbine, ErrStat2, ErrMsg2)
      if (Failed()) return
   end if

   !----------------------------------------------------------------------------
   ! Assemble Jacobian
   !----------------------------------------------------------------------------

   ! If states in Jacobian
   if (p%iJX(1) > 0) then

      ! Group (1,1)
      associate (J11 => m%Mod%Lin%J(p%iJX(1):p%iJX(2), p%iJX(1):p%iJX(2)), &
               dX2dx2 => m%Mod%Lin%dXdx(p%iX2(1):p%iX2(2), p%iX2(1):p%iX2(2)), &
               dX2dx1 => m%Mod%Lin%dXdx(p%iX2(1):p%iX2(2), p%iX1(1):p%iX1(2)))
         J11 = -p%GammaPrime*dX2dx2 - p%BetaPrime*dX2dx1
         do i = p%iJX(1), p%iJX(2)
            J11(i, i) = J11(i, i) + 1.0_R8Ki
         end do
      end associate

      ! Group (2,1)
      if (p%iyT(1) > 0 .and. p%iUT(1) > 0) then
         associate (J21 => m%Mod%Lin%J(p%iJUT(1):p%iJUT(2), p%iJX(1):p%iJX(2)), &
                  dUTdyT => m%Mod%Lin%dUdy(p%iUT(1):p%iUT(2), p%iyT(1):p%iyT(2)), &
                  dYTdx2 => m%Mod%Lin%dYdx(p%iyT(1):p%iyT(2), p%iX2(1):p%iX2(2)), &
                  dYTdx1 => m%Mod%Lin%dYdx(p%iyT(1):p%iyT(2), p%iX1(1):p%iX1(2)))
            ! J21 = C1*matmul(dUTdyT, dYTdx2) + C2*matmul(dUTdyT, dYTdx1)
            call LAPACK_GEMM('N', 'N', p%GammaPrime, dUTdyT, dYTdx2, 0.0_R8Ki, J21, ErrStat2, ErrMsg2); if (Failed()) return
            call LAPACK_GEMM('N', 'N', p%BetaPrime, dUTdyT, dYTdx1, 1.0_R8Ki, J21, ErrStat2, ErrMsg2); if (Failed()) return
         end associate
      end if

      ! Group (1,2)
      if (p%iUT(1) > 0) then
         associate (J12 => m%Mod%Lin%J(p%iJX(1):p%iJX(2), p%iJUT(1):p%iJUT(2)), &
                  dX2duT => m%Mod%Lin%dXdu(p%iX2(1):p%iX2(2), p%iUT(1):p%iUT(2)))
            J12 = -dX2duT
         end associate
      end if

   end if

   ! Group (2,2) - Inputs = dUdu + matmul(dUdy, dYdu)
   if (p%iJU(1) > 0) then
      associate (J22 => m%Mod%Lin%J(p%iJU(1):p%iJU(2), p%iJU(1):p%iJU(2)))
         ! J22 = m%Mod%Lin%dUdu + matmul(m%Mod%Lin%dUdy, m%Mod%Lin%dYdu)
         J22 = m%Mod%Lin%dUdu
         call LAPACK_GEMM('N', 'N', 1.0_R8Ki, m%Mod%Lin%dUdy, m%Mod%Lin%dYdu, 1.0_R8Ki, J22, ErrStat2, ErrMsg2); if (Failed()) return
      end associate
   end if

   ! Write debug matrices if requested
   if (DebugJacobian) then

      ! Get module outputs
      do i = 1, size(m%Mod%ModData)
         associate (ModData => m%Mod%ModData(i))
            call FAST_GetOP(ModData, ThisTime, INPUT_CURR, STATE_PRED, Turbine, ErrStat2, ErrMsg2, &
                            y_op=ModData%Lin%y, y_glue=m%Mod%Lin%y)
            if (Failed()) return
         end associate
      end do

      ! Write debug info
      call BuildJacobian_Debug(m, Turbine, ErrStat2, ErrMsg2)
      if (Failed()) return
   end if

   ! Condition jacobian matrix before factoring
   if (p%iJL(1) > 0) then
      m%Mod%Lin%J(p%iJL(1):p%iJL(2), :) = m%Mod%Lin%J(p%iJL(1):p%iJL(2), :)/p%Scale_UJac
      m%Mod%Lin%J(:, p%iJL(1):p%iJL(2)) = m%Mod%Lin%J(:, p%iJL(1):p%iJL(2))*p%Scale_UJac
   end if

   ! Factor jacobian matrix
   call LAPACK_getrf(size(m%Mod%Lin%J, 1), size(m%Mod%Lin%J, 2), m%Mod%Lin%J, m%IPIV, ErrStat2, ErrMsg2)
   if (Failed()) return

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

!-------------------------------------------------------------------------------
! Utility functions
!-------------------------------------------------------------------------------

subroutine PredictNextState(p, State, Vars)
   type(Glue_TCParam), intent(in)   :: p
   type(TC_State), intent(inout)    :: State
   type(ModVarsType), intent(in)    :: Vars
   real(R8Ki)                       :: v_p, vd_p, a_p
   integer(IntKi)                   :: i

   ! Loop through values and calculate acceleration, algo acceleration, velocity, and delta displacement
   do i = 1, size(State%q)

      ! Store previous velocity, acceleration, and algorithmic acceleration
      v_p = State%v(i)
      vd_p = State%vd(i)
      a_p = State%a(i)

      ! Set acceleration to zero
      State%vd(i) = 0.0_R8Ki

      ! Calculate new algorithmic acceleration
      State%a(i) = (p%AlphaF*vd_p - p%AlphaM*a_p)/(1.0_R8Ki - p%AlphaM)

      ! Calculate new velocity
      State%v(i) = v_p + p%h*(1.0_R8Ki - p%Gamma)*a_p + p%Gamma*p%h*State%a(i)

      ! Copy current displacement to previous displacement
      State%q_prev(i) = State%q(i)

      ! Predict change in displacement
      State%x(i) = p%h*v_p + p%h*p%h*(0.5_R8Ki - p%Beta)*a_p + p%Beta*p%h*p%h*State%a(i)
   end do

   ! Calculate new displacements from delta
   call CalculateStateQ(State, Vars, p%h)
end subroutine

subroutine CalculateStateQ(State, Vars, h)
   type(TC_State), intent(inout) :: State
   type(ModVarsType), intent(in) :: Vars
   real(R8Ki), intent(in)        :: h
   integer(IntKi)                :: i, j, iq
   real(R8Ki)                    :: quat_prev(3), quat_delta(3), quat_new(3)

   ! Calculate new displacement (valid for all states except orientation)
   State%q = State%q_prev + State%x

   ! Loop through variables and compose rotations
   do i = 1, size(Vars%x)
      select case (Vars%x(i)%Field)
      case (FieldOrientation)
         iq = Vars%x(i)%iq(1)
         do j = 1, Vars%x(i)%Nodes
            quat_delta = rvec_to_quat(State%x(iq:iq + 2))
            quat_prev = State%q_prev(iq:iq + 2)
            quat_new = quat_compose(quat_prev, quat_delta)
            State%q(iq:iq + 2) = quat_new
            iq = iq + 3
         end do
      end select
   end do
end subroutine

subroutine UpdateStatePrediction(p, Vars, delta_vd, State)
   type(Glue_TCParam), intent(in)   :: p
   type(ModVarsType), intent(in)    :: Vars
   real(R8Ki), intent(in)           :: delta_vd(:)
   type(TC_State), intent(inout)    :: State

   ! Update x by delta x
   State%x = State%x + p%BetaPrime*delta_vd

   ! Update velocity
   State%v = State%v + p%GammaPrime*delta_vd

   ! Update acceleration
   State%vd = State%vd + delta_vd

   ! Update displacement calculation
   call CalculateStateQ(State, Vars, p%h)

end subroutine

subroutine TransferXtoQ(ModData, x, State)
   type(ModDataType), intent(in)       :: ModData
   real(R8Ki), intent(in)              :: x(:)
   type(TC_State), intent(inout)       :: State
   integer(IntKi)                      :: i
   do i = 1, size(ModData%Vars%x)
      associate (Var => ModData%Vars%x(i))
         select case (Var%DerivOrder)
         case (0) ! Displacement
            State%q(Var%iq(1):Var%iq(2)) = x(Var%iGlu(1):Var%iGlu(2))
         case (1) ! Velocity
            State%v(Var%iq(1):Var%iq(2)) = x(Var%iGlu(1):Var%iGlu(2))
         end select
      end associate
   end do
end subroutine

subroutine TransferQtoX(ModData, State, x)
   type(ModDataType), intent(in) :: ModData
   type(TC_State), intent(in)    :: State
   real(R8Ki), intent(inout)     :: x(:)
   integer(IntKi)                :: i
   do i = 1, size(ModData%Vars%x)
      associate (Var => ModData%Vars%x(i))
         select case (Var%DerivOrder)
         case (0) ! Displacement
            x(Var%iGlu(1):Var%iGlu(2)) = State%q(Var%iq(1):Var%iq(2))
         case (1) ! Velocity
            x(Var%iGlu(1):Var%iGlu(2)) = State%v(Var%iq(1):Var%iq(2))
         end select
      end associate
   end do
end subroutine

subroutine SetBDAccel(ModData, State, BD_OtherSt)
   type(ModDataType), intent(in)          :: ModData
   type(TC_State), intent(in)             :: State
   type(BD_OtherStateType), intent(inout) :: BD_OtherSt
   integer(IntKi)                         :: i
   do i = 1, size(ModData%Vars%x)
      associate (Var => ModData%Vars%x(i))
         select case (Var%Field)
         case (FieldTransVel, FieldAngularVel)
            BD_OtherSt%acc(Var%iLB:Var%iUB, Var%j) = State%vd(Var%iq(1):Var%iq(2))
            BD_OtherSt%xcc(Var%iLB:Var%iUB, Var%j) = State%a(Var%iq(1):Var%iq(2))
         end select
      end associate
   end do
end subroutine

subroutine GetBDAccel(ModData, BD_OtherSt, State)
   type(ModDataType), intent(in)       :: ModData
   type(BD_OtherStateType), intent(in) :: BD_OtherSt
   type(TC_State), intent(inout)       :: State
   integer(IntKi)                      :: i
   do i = 1, size(ModData%Vars%x)
      associate (Var => ModData%Vars%x(i))
         select case (Var%Field)
         case (FieldTransVel, FieldAngularVel)
            State%vd(Var%iq(1):Var%iq(2)) = BD_OtherSt%acc(Var%iLB:Var%iUB, Var%j)
            State%a(Var%iq(1):Var%iq(2)) = BD_OtherSt%xcc(Var%iLB:Var%iUB, Var%j)
         end select
      end associate
   end do
end subroutine

pure function NeedWriteOutput(n_t_global, t_global, t_initial, n_DT_Out) result(WriteNeeded)
   integer(IntKi), intent(in)    :: n_t_global  !< Current global time step
   real(DbKi), intent(in)        :: t_initial   !< Initial time
   real(DbKi), intent(in)        :: t_global    !< Current global time
   integer(IntKi), intent(in)    :: n_DT_Out    !< Write output every n steps
   logical                       :: WriteNeeded !< Function result; if true, WriteOutput values are needed on this time step

   ! note that if TStart isn't an multiple of DT_out, we will not necessarily start output to the file at TStart
   if (t_global >= t_initial) then
      WriteNeeded = MOD(n_t_global, n_DT_Out) == 0
   else
      WriteNeeded = .false.
   end if
end function

!-------------------------------------------------------------------------------
! Debugging routines
!-------------------------------------------------------------------------------

subroutine Solver_Init_Debug(p, m, GlueModData, GlueModMaps)
   type(Glue_TCParam), intent(in)   :: p              !< Parameters
   type(Glue_TCMisc), intent(in)    :: m              !< Misc variables
   type(ModDataType), intent(in)    :: GlueModData(:)  !< Module data
   type(MappingType), intent(in)    :: GlueModMaps(:) !< Module mappings at glue level
   integer(IntKi)                   :: i, j

   write (DebugUn, '(A,*(I6))') " p%iJX2  = ", p%iJX
   write (DebugUn, '(A,*(I6))') " p%iJUT  = ", p%iJUT
   write (DebugUn, '(A,*(I6))') " p%iJU   = ", p%iJU
   write (DebugUn, '(A,*(I6))') " p%iJL   = ", p%iJL
   write (DebugUn, '(A,*(I6))') " p%iX2   = ", p%iX2
   write (DebugUn, '(A,*(I6))') " p%iX1   = ", p%iX1
   write (DebugUn, '(A,*(I6))') " p%iUT   = ", p%iUT
   write (DebugUn, '(A,*(I6))') " p%iU1   = ", p%iU1
   write (DebugUn, '(A,*(I6))') " p%iyT   = ", p%iyT
   write (DebugUn, '(A,*(I6))') " p%iy1   = ", p%iy1
   write (DebugUn, *) "shape(m%dYdx) = ", shape(m%Mod%Lin%dYdx)
   write (DebugUn, *) "shape(m%dYdu) = ", shape(m%Mod%Lin%dYdu)
   write (DebugUn, *) "shape(m%dXdx) = ", shape(m%Mod%Lin%dXdx)
   write (DebugUn, *) "shape(m%dXdu) = ", shape(m%Mod%Lin%dXdu)
   write (DebugUn, *) "shape(m%dUdu) = ", shape(m%Mod%Lin%dUdu)
   write (DebugUn, *) "shape(m%dUdy) = ", shape(m%Mod%Lin%dUdy)

   do j = 1, size(m%Mod%Vars%x)
      write (DebugUn, *) "Var = X "//trim(m%Mod%Vars%x(j)%Name)// &
         " ("//trim(MV_FieldString(m%Mod%Vars%x(j)%Field))//")"
      write (DebugUn, '(A,*(I6))') "  X iLoc = ", m%Mod%Vars%x(j)%iLoc
      write (DebugUn, '(A,*(I6))') "  X iq   = ", m%Mod%Vars%x(j)%iGlu
   end do
   do j = 1, size(m%Mod%Vars%u)
      write (DebugUn, *) "Var = U "//trim(m%Mod%Vars%u(j)%Name)// &
         " ("//trim(MV_FieldString(m%Mod%Vars%u(j)%Field))//")"
      write (DebugUn, '(A,*(I6))') "  U iLoc = ", m%Mod%Vars%u(j)%iLoc
   end do
   do j = 1, size(m%Mod%Vars%y)
      write (DebugUn, *) "Var = Y "//trim(m%Mod%Vars%y(j)%Name)// &
         " ("//trim(MV_FieldString(m%Mod%Vars%y(j)%Field))//")"
      write (DebugUn, '(A,*(I6))') "  Y iLoc = ", m%Mod%Vars%y(j)%iLoc
   end do

   do i = 1, size(GlueModMaps)
      associate (SrcMod => GlueModData(GlueModMaps(i)%iModSrc), &
                 DstMod => GlueModData(GlueModMaps(i)%iModDst))
         write (DebugUn, *) "Mapping = "//GlueModMaps(i)%Desc
         write (DebugUn, *) "   Src = "//trim(SrcMod%Abbr)//' Ins:'//trim(num2lstr(SrcMod%Ins))//' iMod:'//trim(num2lstr(SrcMod%iMod))
         write (DebugUn, *) "   Dst = "//trim(DstMod%Abbr)//' Ins:'//trim(num2lstr(DstMod%Ins))//' iMod:'//trim(num2lstr(DstMod%iMod))
      end associate
   end do
end subroutine

subroutine Solver_Step_Debug(p, m, step, iterCorr, iterConv, delta_norm)
   type(Glue_TCParam), intent(in)   :: p           !< Parameters
   type(Glue_TCMisc), intent(in)    :: m           !< Misc variables
   integer(IntKi), intent(in)       :: step
   integer(IntKi), intent(in)       :: iterCorr
   integer(IntKi), intent(in)       :: iterConv
   real(R8Ki), intent(in)           :: delta_norm

   write (DebugUn, *) "step = ", step
   write (DebugUn, *) "iterCorr = ", iterCorr
   write (DebugUn, *) "iterConv = ", iterConv
   if (p%iJX(1) > 0) write (DebugUn, '(A,*(ES16.7))') " delta_x = ", m%XB(p%iJX(1):p%iJX(2), 1)
   if (p%iJU(1) > 0) write (DebugUn, '(A,*(ES16.7))') " delta_u = ", m%XB(p%iJU(1):p%iJU(2), 1)
   if (allocated(m%uCalc)) write (DebugUn, '(A,*(ES16.7))') " uCalc = ", m%uCalc
   if (allocated(m%Mod%Lin%x)) write (DebugUn, '(A,*(ES16.7))') " x = ", m%Mod%Lin%x
   if (allocated(m%Mod%Lin%u)) write (DebugUn, '(A,*(ES16.7))') " u = ", m%Mod%Lin%u
   write (DebugUn, *) "delta_norm = ", delta_norm
end subroutine

subroutine BuildJacobian_Debug(m, T, ErrStat, ErrMsg)
   type(Glue_TCMisc), intent(inout)    :: m        !< Misc variables
   type(FAST_TurbineType), intent(in)  :: T        !< Turbine
   integer(IntKi), intent(out)         :: ErrStat
   character(*), intent(out)           :: ErrMsg

   character(*), parameter    :: RoutineName = 'BuildJacobian_Debug'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i

   if (MatrixUn == -1) then
      call GetNewUnit(MatrixUn, ErrStat2, ErrMsg2); if (Failed()) return
   end if

   ! Write module matrices to file
   do i = 1, size(m%Mod%ModData)
      associate (ModData => m%Mod%ModData(i))
         call CalcWriteLinearMatrices(ModData%Vars, ModData%Lin, T%p_FAST, T%y_FAST, 0.0_R8Ki, MatrixUn, "SolverTC", VF_None, ErrStat2, ErrMsg2, &
                                      CalcGlue=.false., ModSuffix=ModData%Abbr, FullOutput=.true.)
         if (Failed()) return
      end associate
   end do

   ! Write glue code matrices to file
   call CalcWriteLinearMatrices(m%Mod%Vars, m%Mod%Lin, T%p_FAST, T%y_FAST, 0.0_R8Ki, MatrixUn, "SolverTC", VF_None, ErrStat2, ErrMsg2, CalcGlue=.false., FullOutput=.true.)
   if (Failed()) return

   ! call DumpMatrix(MatrixUn, "dUdu.bin", m%Mod%Lin%dUdu, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(MatrixUn, "dUdy.bin", m%Mod%Lin%dUdy, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(MatrixUn, "dXdu.bin", m%Mod%Lin%dXdu, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(MatrixUn, "dXdx.bin", m%Mod%Lin%dXdx, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(MatrixUn, "dYdu.bin", m%Mod%Lin%dYdu, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(MatrixUn, "dYdx.bin", m%Mod%Lin%dYdx, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(MatrixUn, "ED-dXdu.bin", T%ED%m%Vals%dXdu, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(MatrixUn, "ED-dXdx.bin", T%ED%m%Vals%dXdx, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(MatrixUn, "ED-dYdu.bin", T%ED%m%Vals%dYdu, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(MatrixUn, "ED-dYdx.bin", T%ED%m%Vals%dYdx, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(MatrixUn, "BD-dXdu.bin", T%BD%m(1)%Vals%dXdu, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(MatrixUn, "BD-dXdx.bin", T%BD%m(1)%Vals%dXdx, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(MatrixUn, "BD-dYdu.bin", T%BD%m(1)%Vals%dYdu, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(MatrixUn, "BD-dYdx.bin", T%BD%m(1)%Vals%dYdx, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(MatrixUn, "J.bin", m%Mod%Lin%J, ErrStat2, ErrMsg2); if (Failed()) return

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

end module