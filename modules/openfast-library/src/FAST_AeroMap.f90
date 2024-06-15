!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2024 Envision Energy USA, National Renewable Energy Laboratory
!
!    This file is part of FAST.
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
!> This module contains the routines used by FAST to solve input-output equations and to advance states.

module FAST_AeroMap

use FAST_ModData
use FAST_ModTypes
use FAST_Types
use FAST_Funcs
use FAST_Mapping

use FAST_Subs

implicit none

! Define array of module IDs used in AeroMap
integer(IntKi), parameter  :: AeroMapModIDs(*) = [Module_ED, Module_BD, Module_AD]

real(DbKi), parameter      :: SS_t_global = 0.0_DbKi
real(DbKi), parameter      :: UJacSclFact_x = 1.0d3

logical, parameter         :: output_debugging = .true.

contains

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! DRIVER ROUTINE (runs + ends simulation)
! Put here so that we can call from either stand-alone code or from the ENFAST executable.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine FAST_RunSteadyStateDriver(Turbine)
   type(FAST_TurbineType), intent(inout)  :: Turbine        !< all data for one instance of a turbine

   integer(IntKi)          :: ErrStat        !< Error status of the operation
   character(ErrMsgLen)    :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   ProgName = TRIM(FAST_Ver%Name)//' Aero Map'
   FAST_Ver%Name = ProgName

   call FAST_AeroMapDriver(Turbine%m_Glue, Turbine%p_FAST, Turbine%m_FAST, Turbine%y_FAST, Turbine, ErrStat, ErrMsg)
   call CheckError(ErrStat, ErrMsg, 'FAST_AeroMapDriver')

   call ExitThisProgram_T(Turbine, ErrID_None, .true.)

contains
   subroutine CheckError(ErrID, Msg, SimMsg)
      integer(IntKi), intent(in) :: ErrID       ! The error identifier (ErrStat)
      character(*), intent(in)   :: Msg         ! The error message (ErrMsg)
      character(*), intent(in)   :: SimMsg      ! a message describing the location of the error
      if (ErrID /= ErrID_None) then
         call WrScr(NewLine//TRIM(Msg)//NewLine)
         if (ErrID >= AbortErrLev) then
            call ExitThisProgram_T(Turbine, ErrID, .true., SimMsg)
         end if
      end if
   end subroutine CheckError
end subroutine

subroutine FAST_AeroMapDriver(m, p_FAST, m_FAST, y_FAST, T, ErrStat, ErrMsg)
   use InflowWind_IO, only: IfW_SteadyFlowField_Init
   type(Glue_MiscVarType), intent(inout)     :: m                   !< MiscVars for the glue code
   type(FAST_ParameterType), intent(in)      :: p_FAST              !< Parameters for the glue code
   type(FAST_OutputFileType), intent(inout)  :: y_FAST              !< Output variables for the glue code
   type(FAST_MiscVarType), intent(inout)     :: m_FAST
   type(FAST_TurbineType), intent(inout)     :: T              !< all data for one instance of a turbine
   character(*), intent(out)                 :: ErrMsg         !< Error message if ErrStat /= ErrID_None
   integer(IntKi), intent(out)               :: ErrStat        !< Error status of the operation

   character(*), parameter                   :: RoutineName = 'FAST_AeroMapDriver'
   character(ErrMsgLen)                      :: ErrMsg2
   integer(IntKi)                            :: ErrStat2
   logical, parameter                        :: CompAeroMaps = .true.
   real(DbKi), parameter                     :: t_initial = 0.0_DbKi
   integer(IntKi), allocatable               :: modIDs(:), modIdx(:), iModOrder(:)
   integer(IntKi)                            :: i
   integer(IntKi)                            :: JacSize
   integer(IntKi)                            :: n_case         !< loop counter
   real(DbKi)                                :: n_global
   real(ReKi), allocatable                   :: UnusedAry(:)
   type(AeroMapCase)                         :: CaseDataTmp  ! tsr, windSpeed, pitch, and rotor speed for this case (to try a different operating point first)
   integer(IntKi)                            :: NStatus
   character(MaxWrScrLen), parameter         :: BlankLine = " "

   ErrStat = ErrID_None
   ErrMsg = ''

   !----------------------------------------------------------------------------
   ! Initialization
   !----------------------------------------------------------------------------

   ! Set Turbine ID
   T%TurbID = 1

   ! Initialize linearization file number (will be incremented before use)
   m%AM%LinFileNum = 0

   ! Standard Turbine initialization
   call FAST_InitializeAll(t_initial, T%m_Glue, T%p_FAST, T%y_FAST, T%m_FAST, &
                           T%ED, T%BD, T%SrvD, T%AD14, T%AD, &
                           T%ExtLd, T%IfW, T%ExtInfw, T%SC_DX, &
                           T%SeaSt, T%HD, T%SD, T%ExtPtfm, T%MAP, &
                           T%FEAM, T%MD, T%Orca, T%IceF, T%IceD, &
                           T%MeshMapData, CompAeroMaps, ErrStat2, ErrMsg2)
   if (Failed()) return

   ! TODO: Move into FAST_InitializeAll
   ! Initialize module data transfer mappings
   call FAST_InitMappings(m%Modules, m%Mappings, T, ErrStat2, ErrMsg2)
   if (Failed()) return

   ! Initialize steady flow field in AeroDyn
   call IfW_SteadyFlowField_Init(T%AD%p%FlowField, &
                                 RefHt=100.0_ReKi, HWindSpeed=8.0_ReKi, PLExp=0.0_ReKi, &
                                 ErrStat=ErrStat2, ErrMsg=ErrMsg2)
   if (Failed()) return

   !----------------------------------------------------------------------------
   ! Modules
   !----------------------------------------------------------------------------

   ! Initialize module indices
   m%AM%iModED = 0
   m%AM%iModBD = 0
   m%AM%iModAD = 0

   ! Get indices of modules that are used by Aero Mapping (first instance only)
   call GetModuleOrder(m%Modules, AeroMapModIDs, m%AM%iModOrder)
   do i = 1, size(m%AM%iModOrder)
      associate (ModData => m%Modules(m%AM%iModOrder(i)))
         if (ModData%Ins == 1) then
            select case (ModData%ID)
            case (Module_ED)
               m%AM%iModED = i
            case (Module_BD)
               m%AM%iModBD = i
            case (Module_AD)
               m%AM%iModAD = i
            end select
         end if
      end associate
   end do

   ! If BeamDyn is active
   if (m%AM%iModBD > 0) then
      m%AM%iModED = 0
      m%AM%iModOrder = [m%AM%iModBD, m%AM%iModAD]
   else if (m%AM%iModED > 0) then
      m%AM%iModOrder = [m%AM%iModED, m%AM%iModAD]
   end if

   ! Loop through module indices, copy states and inputs
   do i = 1, size(m%AM%iModOrder)
      associate (ModData => m%Modules(m%AM%iModOrder(i)))

         ! Copy current state to predicted state
         call FAST_CopyStates(ModData, T, STATE_CURR, STATE_PRED, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         if (Failed()) return

         ! Copy current inputs to previous inputs
         call FAST_CopyInput(ModData, T, INPUT_CURR, INPUT_PREV, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         if (Failed()) return

         ! If linearization is enabled, set lin file module abbreviation for file name
         ! If module is BeamDyn or more than one instance, append instance number to abbreviation
         if ((ModData%ID == Module_BD) .or. (count(m%Modules%ID == ModData%ID) > 1)) then
            ModData%Lin%Abbr = trim(ModData%Abbr)//Num2LStr(ModData%Ins)
         else
            ModData%Lin%Abbr = ModData%Abbr
         end if

      end associate
   end do

   !----------------------------------------------------------------------------
   ! Build AeroMap module
   !----------------------------------------------------------------------------

   ! Generate index for variables with AeroMap flag
   call ModD_CombineModules(m%Modules, m%AM%iModOrder, VF_AeroMap, .true., m%AM%Mod, ErrStat2, ErrMsg2)
   if (Failed()) return

   !----------------------------------------------------------------------------
   ! Allocation
   !----------------------------------------------------------------------------

   ! Allocate components of the Jacobian matrix
   call AllocAry(m%AM%Jac11, m%AM%Mod%Vars%Nx, m%AM%Mod%Vars%Nx, 'Jac11', ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%AM%Jac12, m%AM%Mod%Vars%Nx, m%AM%Mod%Vars%Nu, 'Jac12', ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%AM%Jac21, m%AM%Mod%Vars%Nu, m%AM%Mod%Vars%Nx, 'Jac21', ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%AM%Jac22, m%AM%Mod%Vars%Nu, m%AM%Mod%Vars%Nu, 'Jac22', ErrStat2, ErrMsg2); if (Failed()) return

   ! Jacobian size is number of states plus number of inputs
   JacSize = m%AM%Mod%Vars%Nx + m%AM%Mod%Vars%Nu

   ! Allocate Jacobian pivot vector
   call AllocAry(m%AM%JacPivot, JacSize, 'Pivot array for Jacobian LU decomposition', ErrStat2, ErrMsg2); if (Failed()) return

   ! Storage for residual and solution delta
   call AllocAry(m%AM%Residual, JacSize, 'Residual', ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%AM%SolveDelta, JacSize, 'SolveDelta', ErrStat2, ErrMsg2); if (Failed()) return

   ! Allocate Jacobian matrix
   call AllocAry(m%AM%Mod%Lin%J, JacSize, JacSize, 'J', ErrStat2, ErrMsg2); if (Failed()) return

   ! Allocate Idx Jacobian storage
   call AllocAry(m%AM%Mod%Lin%dYdu, m%AM%Mod%Vars%Ny, m%AM%Mod%Vars%Nu, 'dYdu', ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%AM%Mod%Lin%dXdu, m%AM%Mod%Vars%Nx, m%AM%Mod%Vars%Nu, 'dXdu', ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%AM%Mod%Lin%dYdx, m%AM%Mod%Vars%Ny, m%AM%Mod%Vars%Nx, 'dYdx', ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%AM%Mod%Lin%dXdx, m%AM%Mod%Vars%Nx, m%AM%Mod%Vars%Nx, 'dXdx', ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%AM%Mod%Lin%dXdy, m%AM%Mod%Vars%Nx, m%AM%Mod%Vars%Ny, 'dXdy', ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%AM%Mod%Lin%dUdu, m%AM%Mod%Vars%Nu, m%AM%Mod%Vars%Nu, "dUdu", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%AM%Mod%Lin%dUdy, m%AM%Mod%Vars%Nu, m%AM%Mod%Vars%Ny, "dUdy", ErrStat2, ErrMsg2); if (Failed()) return

   ! Allocate operating point arrays
   if (output_debugging) then
      call AllocAry(m%AM%Mod%Lin%x, m%AM%Mod%Vars%Nx, 'x', ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(m%AM%Mod%Lin%dx, m%AM%Mod%Vars%Nx, 'dx', ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(m%AM%Mod%Lin%u, m%AM%Mod%Vars%Nu, 'u', ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(m%AM%Mod%Lin%y, m%AM%Mod%Vars%Ny, 'y', ErrStat2, ErrMsg2); if (Failed()) return
   end if

   ! Allocate arrays to store inputs
   call AllocAry(m%AM%u1, m%AM%Mod%Vars%Nu, 'u1', ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(m%AM%u2, m%AM%Mod%Vars%Nu, 'u2', ErrStat2, ErrMsg2); if (Failed()) return

   ! Move hub orientation matrices to AeroMap structure
   call move_alloc(T%MeshMapData%HubOrient, m%AM%HubOrientation)

   ! Allocate cases
   allocate (m%AM%Cases(p_FAST%NumSSCases), stat=ErrStat2)
   if (ErrStat2 /= 0) then
      call SetErrStat(ErrID_Fatal, "Error allocating AeroMap cases", ErrStat, ErrMsg, RoutineName)
      return
   end if

   ! Populate case data
   do n_case = 1, p_FAST%NumSSCases
      if (p_FAST%WindSpeedOrTSR == 1) then
         m%AM%Cases(n_case)%WindSpeed = p_FAST%WS_TSR(n_case)
         m%AM%Cases(n_case)%TSR = p_FAST%RotSpeed(n_case)*T%AD%p%rotors(1)%BEMT%rTipFixMax/m%AM%Cases(n_case)%WindSpeed
      else
         m%AM%Cases(n_case)%TSR = p_FAST%WS_TSR(n_case)
         m%AM%Cases(n_case)%WindSpeed = p_FAST%RotSpeed(n_case)*T%AD%p%rotors(1)%BEMT%rTipFixMax/m%AM%Cases(n_case)%TSR
      end if
      m%AM%Cases(n_case)%Pitch = p_FAST%Pitch(n_case)
      m%AM%Cases(n_case)%RotSpeed = p_FAST%RotSpeed(n_case)
   end do

   !----------------------------------------------------------------------------
   ! Calculate steady-state solution for each case
   !----------------------------------------------------------------------------

   ! how often do we inform the user which case we are on?
   NStatus = min(100, p_FAST%NumSSCases/100 + 1) ! at least 100 every 100 cases or 100 times per simulation
   call WrScr(NewLine)

   ! Loop through Aero Map cases
   do n_case = 1, p_FAST%NumSSCases

      ! If status should be written to screen
      if (n_case == 1 .or. n_case == p_FAST%NumSSCases .or. mod(n_case, NStatus) == 0) then
         call WrOver(' Case '//trim(Num2LStr(n_case))//' of '//trim(Num2LStr(p_FAST%NumSSCases)))
      end if

      ! Call steady-state solve for this pitch and rotor speed
      call SolveSteadyState(m, m%AM%Cases(n_case), p_FAST, y_FAST, m_FAST, T%MeshMapData, T, ErrStat2, ErrMsg2)

      ! we didn't converge; let's try a different operating point and see if that helps:
      if (ErrStat2 >= ErrID_Severe) then

         ! Create copy of case data for second attempt
         CaseDataTmp = m%AM%Cases(n_case)

         ! Modify pitch, TSR, and WindSpeed
         CaseDataTmp%Pitch = CaseDataTmp%Pitch*0.5_ReKi
         CaseDataTmp%TSR = CaseDataTmp%TSR*0.5_ReKi
         CaseDataTmp%WindSpeed = CaseDataTmp%WindSpeed*0.5_ReKi

         ! Write message about retrying case
         call WrScr('Retrying case '//trim(Num2LStr(n_case))//', first trying to get a better initial guess. Average error is '// &
                    trim(Num2LStr(y_FAST%DriverWriteOutput(SS_Indx_Err)))//'.')

         ! call SolveSteadyState(m, CaseDataTmp, p_FAST, y_FAST, m_FAST, T%ED, T%BD, T%AD, T%MeshMapData, T, ErrStat2, ErrMsg2)

         ! if that worked, try the real case again:
         if (ErrStat2 < AbortErrLev) then
            ! call SolveSteadyState(m, m%AM%Cases(n_case), p_FAST, y_FAST, m_FAST, T%ED, T%BD, T%AD, T%MeshMapData, T, ErrStat2, ErrMsg2)
            call WrOver(BlankLine)
         end if

      end if

      if (ErrStat2 > ErrID_None) then
         ErrMsg2 = trim(ErrMsg2)//" case "//trim(Num2LStr(n_case))// &
                   ' (tsr='//trim(Num2LStr(m%AM%Cases(n_case)%tsr))// &
                   ', wind speed='//trim(Num2LStr(m%AM%Cases(n_case)%windSpeed))//' m/s'// &
                   ', pitch='//trim(num2lstr(m%AM%Cases(n_case)%pitch*R2D))//' deg'// &
                   ', rotor speed='//trim(num2lstr(m%AM%Cases(n_case)%RotSpeed*RPS2RPM))//' rpm)'
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      end if

      !-------------------------------------------------------------------------
      ! Write results to file
      !-------------------------------------------------------------------------

      n_global = real(n_case, DbKi) ! n_global is double-precision so that we can reuse existing code.

      call WrOutputLine(n_global, p_FAST, y_FAST, UnusedAry, UnusedAry, T%ED%y%WriteOutput, &
                        T%AD%y, UnusedAry, UnusedAry, UnusedAry, UnusedAry, UnusedAry, UnusedAry, &
                        UnusedAry, UnusedAry, UnusedAry, UnusedAry, T%IceD%y, T%BD%y, ErrStat2, ErrMsg2)
      if (Failed()) return

      !-------------------------------------------------------------------------
      ! Write errors to screen
      !-------------------------------------------------------------------------

      if (ErrStat > ErrID_None) then
         call WrScr(trim(ErrMsg))
         call WrScr("")
         ErrStat = ErrID_None
         ErrMsg = ""
      end if

   end do

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine performs the Input-Output solve for the steady-state solver.
!! Note that this has been customized for the physics in the problems and is not a general solution.
subroutine SolveSteadyState(m, caseData, p_FAST, y_FAST, m_FAST, MeshMapData, T, ErrStat, ErrMsg)
   type(Glue_MiscVarType), intent(inout)     :: m              !< Miscellaneous variables
   type(AeroMapCase), intent(in)             :: caseData       !< tsr, windSpeed, pitch, and rotor speed for this case
   type(FAST_ParameterType), intent(in)      :: p_FAST         !< Glue-code simulation parameters
   type(FAST_OutputFileType), intent(inout)  :: y_FAST         !< Glue-code output file values
   type(FAST_MiscVarType), intent(inout)     :: m_FAST         !< Miscellaneous variables
   type(FAST_TurbineType), intent(inout)     :: T              !< Turbine type
   type(FAST_ModuleMapType), intent(inout)   :: MeshMapData    !< data for mapping meshes between modules
   integer(IntKi), intent(out)               :: ErrStat        !< Error status of the operation
   character(*), intent(out)                 :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   character(*), parameter :: RoutineName = 'SolveSteadyState'
   integer(IntKi)          :: ErrStat2                  ! temporary Error status of the operation
   character(ErrMsgLen)    :: ErrMsg2                   ! temporary Error message if ErrStat /= ErrID_None

   !bjj: store these so that we don't reallocate every time?
   real(R8Ki)                                        :: err
   real(R8Ki)                                        :: err_prev
   real(R8Ki), parameter                             :: reduction_factor = 0.1_R8Ki

   integer(IntKi)                                    :: nb                        ! loop counter (blade number)
   integer(IntKi)                                    :: MaxIter                   ! maximum number of iterations
   integer(IntKi)                                    :: K                         ! Input-output-solve iteration counter
   integer(IntKi)                                    :: i, j

   logical                                           :: GetWriteOutput            ! flag to determine if we need WriteOutputs from this call to CalcOutput

   ! Note: p_FAST%UJacSclFact is a scaling factor that gets us similar magnitudes between loads and accelerations...

   !bjj: note, that this routine may have a problem if there is remapping done

   ErrStat = ErrID_None
   ErrMsg = ""

   !----------------------------------------------------------------------------
   ! Some record keeping stuff:
   !----------------------------------------------------------------------------

   ! Set the rotor speed in ElastoDyn
   T%ED%x(STATE_CURR)%QDT(p_FAST%GearBox_Index) = caseData%RotSpeed

   ! Update module inputs
   call SetPrescribedInputs(caseData, p_FAST, y_FAST, m_FAST, T%ED, T%BD, T%AD)
   do i = 1, size(m%AM%iModOrder)
      associate (ModData => m%Modules(m%AM%iModOrder(i)))
         call FAST_CopyInput(ModData, T, INPUT_CURR, INPUT_PREV, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      end associate
   end do

   K = 0
   err = 1.0E3
   err_prev = err

   y_FAST%DriverWriteOutput(SS_Indx_Err) = -1
   y_FAST%DriverWriteOutput(SS_Indx_Iter) = 0
   y_FAST%DriverWriteOutput(SS_Indx_TSR) = caseData%tsr
   y_FAST%DriverWriteOutput(SS_Indx_WS) = caseData%windSpeed
   y_FAST%DriverWriteOutput(SS_Indx_Pitch) = caseData%Pitch*R2D
   y_FAST%DriverWriteOutput(SS_Indx_RotSpeed) = caseData%RotSpeed*RPS2RPM

   MaxIter = p_FAST%KMax + 1 ! adding 1 here so that we get the error calculated correctly when we hit the max iteration
   do

      !-------------------------------------------------------------------------
      ! Calculate outputs, based on inputs at this time
      !-------------------------------------------------------------------------

      ! Set GetWriteOutput flag true if not the first iteration
      GetWriteOutput = K > 0

      !-----------------------------------------
      ! Caclulate ElastoDyn / BeamDyn output
      !-----------------------------------------

      ! If BeamDyn is active
      if (m%AM%iModBD > 0) then

         ! Calculate BeamDyn output
         call FAST_CalcOutput(m%Modules(m%AM%iModBD), m%Mappings, SS_t_global, INPUT_CURR, STATE_CURR, T, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      else

         ! Calculate ElastoDyn output
         call FAST_CalcOutput(m%Modules(m%AM%iModED), m%Mappings, SS_t_global, INPUT_CURR, STATE_CURR, T, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      end if

      !-----------------------------------------
      ! AeroDyn InputSolve
      !-----------------------------------------

      ! If first iteration
      if (K == 0) then

         ! Perform AeroDyn input solve to get initial guess from structural module
         ! (this ensures that the pitch is accounted for in the fixed aero-map solve:):
         call SS_AD_InputSolve(m, INPUT_CURR, T, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         call SS_AD_InputSolve_OtherBlades(m, INPUT_CURR, T)

         ! set up x-u vector, using local initial guesses
         call SS_GetInputs(m, m%AM%u1, p_FAST%UJacSclFact, INPUT_CURR, STATE_CURR, T, ErrStat2, ErrMsg2)

      end if

      !-----------------------------------------
      ! Calculate AeroDyn Output
      !-----------------------------------------

      call FAST_CalcOutput(m%Modules(m%AM%iModAD), m%Mappings, SS_t_global, INPUT_CURR, STATE_CURR, T, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) then
         call ResetInputsAndStates()
         return
      end if

      ! If iteration is at or above maximum iteration, exit loop
      if (K >= MaxIter) exit

      !-------------------------------------------------------------------------------------------------
      ! Calculate residual and the Jacobian
      ! (note that we don't want to change module%Input(1), here)
      ! Also, the residual uses values from y_FAST, so do this before calculating the jacobian
      !-------------------------------------------------------------------------------------------------

      call SS_BuildResidual(caseData, m, T, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) then
         call ResetInputsAndStates()
         return
      end if

      ! If Jacobian needs to be recalculated
      if (mod(K, p_FAST%N_UJac) == 0) then
         call SS_BuildJacobian(m, caseData, p_FAST, y_FAST, m_FAST, T, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) then
            call ResetInputsAndStates()
            return
         end if
      end if

      !-------------------------------------------------------------------------
      ! Solve for delta u: Jac*u_delta = - Fn_U_Resid
      !  using the LAPACK routine
      !-------------------------------------------------------------------------

      ! Copy negative of residual into solve
      m%AM%SolveDelta = -m%AM%Residual

      ! Solve for changes in states and inputs
      call LAPACK_getrs(TRANS="N", N=size(m%AM%Mod%Lin%J, 1), A=m%AM%Mod%Lin%J, &
                        IPIV=m%AM%JacPivot, B=m%AM%SolveDelta, ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return

      !-------------------------------------------------------------------------
      ! Check for error, update inputs if necessary, and iterate again
      !-------------------------------------------------------------------------

      ! Save previous error
      err_prev = err

      ! Calculate new error
      err = dot_product(m%AM%SolveDelta, m%AM%SolveDelta)

      ! Store normalized error in output
      y_FAST%DriverWriteOutput(SS_Indx_Err) = sqrt(err)/size(m%AM%Mod%Lin%J, 1)

      ! If error is below tolerance
      if (err <= p_FAST%TolerSquared) then
         if (K == 0) then ! the error will be incorrect in this instance, but the outputs will be better
            MaxIter = K
         else
            exit
         end if
      end if

      if (K >= p_FAST%KMax) exit
      if (K > 5 .and. err > 1.0E35) exit ! this is obviously not converging. Let's try something else.

      !-------------------------------------------------------------------------
      ! Modify inputs and states for next iteration
      !-------------------------------------------------------------------------

      ! If current error is greater than previous error (solution diverging),
      ! reduce delta (take a smaller step)
      if (err > err_prev) then
         m%AM%SolveDelta = m%AM%SolveDelta*reduction_factor 
         err_prev = err_prev*reduction_factor
      end if

      ! ! call Add_SteadyState_delta(p_FAST, y_FAST, u_delta, AD, ED, BD, MeshMapData)

      ! !u = u + u_delta
      ! call SS_GetInputs(m, T, p_FAST%UJacSclFact, INPUT_CURR, STATE_CURR, ErrStat2, ErrMsg2)

      ! K = K + 1
      ! y_FAST%DriverWriteOutput(SS_Indx_Iter) = k

   end do ! K

   ! if (p_FAST%CompElast == Module_BD) then
   !    ! this doesn't actually get the correct hub point load from BD, but we'll get some outputs:
   !    call ED_CalcOutput(SS_t_global, ED%Input(1), ED%p, ED%x(STATE_CURR), ED%xd(STATE_CURR), ED%z(STATE_CURR), ED%OtherSt(STATE_CURR), ED%y, ED%m, ErrStat2, ErrMsg2)
   !    call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   ! end if

   call ResetInputsAndStates()

contains
   subroutine ResetInputsAndStates()

      if (err > p_FAST%TolerSquared) then

         call SetErrStat(ErrID_Severe, 'Steady-state solver did not converge.', ErrStat, ErrMsg, RoutineName)

         if (err > 100.0) then
            ! if we didn't get close on the solution, we should reset the states and inputs because they very well could
            ! lead to numerical issues on the next iteration. Here, set the initial values to 0:

            ! because loads occasionally get very large when it fails, manually set these to zero (otherwise
            ! roundoff can lead to non-zero values with the method below, which is most useful for states)
            if (p_FAST%CompElast == Module_BD) then
               do K = 1, p_FAST%nBeams
                  T%BD%Input(1, k)%DistrLoad%Force = 0.0_ReKi
                  T%BD%Input(1, k)%DistrLoad%Moment = 0.0_ReKi
               end do
            end if

            call SS_GetInputs(m, m%AM%u1, p_FAST%UJacSclFact, INPUT_CURR, STATE_CURR, T, ErrStat2, ErrMsg2)     ! find the values we have been modifying (in u... continuous states and inputs)
            ! call Add_SteadyState_delta(p_FAST, y_FAST, -u, AD, ED, BD, MeshMapData) ! and reset them to 0 (by adding -u)

         end if
      end if
   end subroutine ResetInputsAndStates

end subroutine SolveSteadyState

subroutine SS_BuildJacobian(m, caseData, p_FAST, y_FAST, m_FAST, T, ErrStat, ErrMsg)
   type(Glue_MiscVarType), intent(inout)     :: m              !< Miscellaneous variables
   type(AeroMapCase), intent(in)             :: caseData       !< tsr, windSpeed, pitch, and rotor speed for this case
   type(FAST_ParameterType), intent(IN)      :: p_FAST         !< Parameters for the glue code
   type(FAST_OutputFileType), intent(INOUT)  :: y_FAST         !< Output variables for the glue code
   type(FAST_MiscVarType), intent(INOUT)     :: m_FAST         !< Miscellaneous variables
   type(FAST_TurbineType), intent(inout)     :: T              !< Turbine type
   integer(IntKi), intent(OUT)               :: ErrStat        !< Error status of the operation
   character(*), intent(OUT)                 :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   character(*), parameter   :: RoutineName = 'SS_BuildJacobian'
   integer(IntKi)            :: ErrStat2
   character(ErrMsgLen)      :: ErrMSg2
   character(1024)           :: LinRootName
   integer(IntKi)            :: i, j, k, c, r, iRow(2), iCol(2)
   integer(IntKi)            :: nx        ! Number of states
   integer(IntKi)            :: Un

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Set number of states
   nx = m%AM%Mod%Vars%Nx

   ! If output debugging is requested
   if (output_debugging) then

      ! Get unit number for output files
      call GetNewUnit(Un, ErrStat2, ErrMsg2)
      if (Failed()) return

      ! Build linearization root name
      m%AM%LinFileNum = m%AM%LinFileNum + 1
      LinRootName = trim(p_FAST%OutFileRoot)//'.'//trim(Num2LStr(m%AM%LinFileNum))

      ! These values get printed in the linearization output files, so we'll set them here:
      y_FAST%Lin%WindSpeed = caseData%WindSpeed
      y_FAST%Lin%RotSpeed = caseData%RotSpeed
      y_FAST%Lin%Azimuth = 0.0_ReKi
   end if

   ! Initialize Jacobian
   m%AM%Mod%Lin%J = 0.0_R8Ki

   !----------------------------------------------------------------------------
   ! dXdy
   !----------------------------------------------------------------------------

   m%AM%Mod%Lin%dXdy = 0.0_R8Ki

   !----------------------------------------------------------------------------
   ! Module Jacobians
   !----------------------------------------------------------------------------

   ! Loop through modules
   do i = 1, size(m%AM%iModOrder)
      associate (ModData => m%Modules(m%AM%iModOrder(i)), iMod => m%AM%iModOrder(i))

         ! Calculate dYdu and dXdu
         call FAST_JacobianPInput(ModData, SS_t_global, STATE_CURR, T, ErrStat, ErrMsg, &
                                  FlagFilter=VF_AeroMap, dYdu=ModData%Lin%dYdu, dXdu=ModData%Lin%dXdu)
         if (Failed()) return

         ! Calculate dYdx and dXdx
         call FAST_JacobianPContState(ModData, SS_t_global, STATE_CURR, T, ErrStat, ErrMsg, &
                                      FlagFilter=VF_AeroMap, dYdx=ModData%Lin%dYdx, dXdx=ModData%Lin%dXdx)
         if (Failed()) return

         ! If output debugging requested
         if (output_debugging) then

            ! Calculate operating point values
            call FAST_GetOP(ModData, SS_t_global, INPUT_CURR, STATE_CURR, T, ErrStat2, ErrMsg2, &
                            u_op=ModData%Lin%u, y_op=ModData%Lin%y, x_op=ModData%Lin%x, dx_op=ModData%Lin%dx)
            if (Failed()) return

            ! Write linearization matrices
            call CalcWriteLinearMatrices(ModData, p_FAST, y_FAST, SS_t_global, Un, LinRootName, VF_AeroMap, .false., ErrStat2, ErrMsg2)
            if (Failed()) return

            ! Pack values into module
            if (allocated(ModData%Lin%x)) call ModD_PackAry(m%AM%Mod%Xfr(iMod)%x, ModData%Lin%x, m%AM%Mod%Lin%x)
            if (allocated(ModData%Lin%dx)) call ModD_PackAry(m%AM%Mod%Xfr(iMod)%x, ModData%Lin%dx, m%AM%Mod%Lin%dx)
            if (allocated(ModData%Lin%u)) call ModD_PackAry(m%AM%Mod%Xfr(iMod)%u, ModData%Lin%u, m%AM%Mod%Lin%u)
            if (allocated(ModData%Lin%y)) call ModD_PackAry(m%AM%Mod%Xfr(iMod)%y, ModData%Lin%y, m%AM%Mod%Lin%y)

         end if

         ! If this module is BeamDyn, calculate dxdotdy
         if (ModData%ID == Module_BD) then

            ! TODO: implement beamdyn
            ! NOTE that this implies that the FEA nodes (states) are the same as the output nodes!!!! (note that we have overlapping nodes at the element end points)
            ! r = 1
            ! do i = 2, BD%p(k)%node_total ! the first node isn't technically a state
            !    c = (BD%p(k)%NdIndx(i) - 1)*3 + 1 ! BldMeshNode = BD%p(k)%NdIndx(i)

            !    !dxdotdy(r:r+2,c:c+2) = SkewSymMat( [p_FAST%RotSpeed, 0.0_ReKi, 0.0_ReKi] )
            !    dxdotdy(r + 2, c + 1) = caseData%RotSpeed
            !    dxdotdy(r + 1, c + 2) = -caseData%RotSpeed

            !    ! derivative
            !    dxdotdy(r + nx + 1, c + 1) = -OmegaSquared
            !    dxdotdy(r + nx + 2, c + 2) = -OmegaSquared

            !    r = r + BD%p(k)%dof_node
            ! end do
         end if

         ! Add module Jacobians to global Jacobians
         if (allocated(ModData%Lin%dYdu)) call ModD_PackMatrix(m%AM%Mod%Xfr(iMod)%y, m%AM%Mod%Xfr(iMod)%u, ModData%Lin%dYdu, m%AM%Mod%Lin%dYdu)
         if (allocated(ModData%Lin%dXdu)) call ModD_PackMatrix(m%AM%Mod%Xfr(iMod)%x, m%AM%Mod%Xfr(iMod)%u, ModData%Lin%dXdu, m%AM%Mod%Lin%dXdu)
         if (allocated(ModData%Lin%dYdx)) call ModD_PackMatrix(m%AM%Mod%Xfr(iMod)%y, m%AM%Mod%Xfr(iMod)%x, ModData%Lin%dYdx, m%AM%Mod%Lin%dYdx)
         if (allocated(ModData%Lin%dXdx)) call ModD_PackMatrix(m%AM%Mod%Xfr(iMod)%x, m%AM%Mod%Xfr(iMod)%x, ModData%Lin%dXdx, m%AM%Mod%Lin%dXdx)

      end associate
   end do

   !----------------------------------------------------------------------------
   ! Glue Jacobians
   !----------------------------------------------------------------------------

   m%AM%Mod%Lin%dUdy = 0.0_R8Ki
   call Eye2D(m%AM%Mod%Lin%dUdu, ErrStat2, ErrMsg2); if (Failed()) return
   call FAST_LinearizeMappings(T, m%Modules, m%Mappings, m%AM%iModOrder, m%AM%Mod%Xfr, ErrStat2, ErrMsg2, &
                               m%AM%Mod%Lin%dUdu, m%AM%Mod%Lin%dUdy)
   if (Failed()) return

   !----------------------------------------------------------------------------
   ! Form Jacobian matrix
   !----------------------------------------------------------------------------

   ! Calculate Jacobian block 11 = dX/dx - dX/dy * dY/dx
   m%AM%Jac11 = m%AM%Mod%Lin%dXdx
   call LAPACK_GEMM('N', 'N', -1.0_R8Ki, m%AM%Mod%Lin%dXdy, m%AM%Mod%Lin%dYdx, 1.0_R8Ki, m%AM%Jac11, ErrStat2, ErrMsg2)
   if (Failed()) return

   ! Calculate Jacobian block 12 = dX/du - dX/dy * dY/du
   m%AM%Jac12 = m%AM%Mod%Lin%dXdu
   call LAPACK_GEMM('N', 'N', -1.0_R8Ki, m%AM%Mod%Lin%dXdy, m%AM%Mod%Lin%dYdu, 1.0_R8Ki, m%AM%Jac12, ErrStat2, ErrMsg2)
   if (Failed()) return

   ! Calculate Jacobian block 21 = dU/dy * dY/dx
   call LAPACK_GEMM('N', 'N', 1.0_R8Ki, m%AM%Mod%Lin%dUdy, m%AM%Mod%Lin%dYdx, 0.0_R8Ki, m%AM%Jac21, ErrStat2, ErrMsg2)
   if (Failed()) return

   ! Calculate Jacobian block 22 = dU/du + dU/dy * dY/du
   m%AM%Jac22 = m%AM%Mod%Lin%dUdu
   call LAPACK_GEMM('N', 'N', 1.0_R8Ki, m%AM%Mod%Lin%dUdy, m%AM%Mod%Lin%dYdu, 1.0_R8Ki, m%AM%Jac22, ErrStat2, ErrMsg2)
   if (Failed()) return

   ! Assemble blocks to form full Jacobian
   m%AM%Mod%Lin%J(:nx, :nx) = m%AM%Jac11
   m%AM%Mod%Lin%J(:nx, nx + 1:) = m%AM%Jac12
   m%AM%Mod%Lin%J(nx + 1:, :nx) = m%AM%Jac21
   m%AM%Mod%Lin%J(nx + 1:, nx + 1:) = m%AM%Jac22

   ! If output debugging is enabled, write combined matrices and Jacobian
   if (output_debugging) then
      call CalcWriteLinearMatrices(m%AM%Mod, p_FAST, y_FAST, SS_t_global, Un, LinRootName, VF_AeroMap, .false., ErrStat2, ErrMsg2)
      if (Failed()) return
   end if

   !----------------------------------------------------------------------------
   ! Condition Jacobian matrix
   !----------------------------------------------------------------------------

   ! Loop through inputs
   do c = 1, size(m%AM%Mod%Vars%u)

      iCol = m%AM%Mod%Vars%u(c)%iLoc + nx

      ! If column is a load
      if (MV_IsLoad(m%AM%Mod%Vars%u(c))) then

         ! Column is a load, state rows are not loads
         m%AM%Mod%Lin%J(1:nx, iCol(1):iCol(2)) = &
            m%AM%Mod%Lin%J(1:nx, iCol(1):iCol(2))*p_FAST%UJacSclFact

         ! Loop through rows
         do r = 1, size(m%AM%Mod%Vars%u)
            ! If column is load, but row is a motion
            if (.not. MV_IsLoad(m%AM%Mod%Vars%u(r))) then
               iRow = m%AM%Mod%Vars%u(r)%iLoc + nx
               m%AM%Mod%Lin%J(iRow(1):iRow(2), iCol(1):iCol(2)) = &
                  m%AM%Mod%Lin%J(iRow(1):iRow(2), iCol(1):iCol(2))*p_FAST%UJacSclFact
            end if
         end do

      else

         ! Loop through rows
         do r = 1, size(m%AM%Mod%Vars%u)
            ! Column is a motion, but row is a load
            if (MV_IsLoad(m%AM%Mod%Vars%u(r))) then
               iRow = m%AM%Mod%Vars%u(r)%iLoc + nx
               m%AM%Mod%Lin%J(iRow(1):iRow(2), iCol(1):iCol(2)) = &
                  m%AM%Mod%Lin%J(iRow(1):iRow(2), iCol(1):iCol(2))/p_FAST%UJacSclFact
            end if
         end do

      end if
   end do

   !----------------------------------------------------------------------------
   ! Factor Jacobian matrix
   ! Get the LU decomposition of this matrix using a LAPACK routine:
   ! The result is of the form Jmat = P * L * U
   !----------------------------------------------------------------------------

   call LAPACK_getrf(M=size(m%AM%Mod%Lin%J, 1), N=size(m%AM%Mod%Lin%J, 2), &
                     A=m%AM%Mod%Lin%J, IPIV=m%AM%JacPivot, ErrStat=ErrStat2, ErrMsg=ErrMsg2)
   if (Failed()) return

contains

   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
      if (Failed) call Cleanup()
   end function

   subroutine Cleanup()
      if (Un > 0) close (Un)
   end subroutine Cleanup

end subroutine SS_BuildJacobian

!----------------------------------------------------------------------------------------------------------------------------------
subroutine SS_BuildResidual(caseData, m, T, ErrStat, ErrMsg)
   type(AeroMapCase), intent(IN)             :: caseData    !< tsr, windSpeed, pitch, and rotor speed for this case
   type(Glue_MiscVarType), intent(INOUT)     :: m           !< Miscellaneous variables
   type(FAST_TurbineType), intent(INOUT)     :: T           !< Turbine type
   integer(IntKi), intent(OUT)               :: ErrStat     !< Error status of the operation
   character(*), intent(OUT)                 :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   character(*), parameter                   :: RoutineName = 'SS_BuildResidual'
   integer(IntKi)                            :: ErrStat2
   character(ErrMsgLen)                      :: ErrMsg2
   integer(IntKi)                            :: i, j, iVarMod(2), iVarGbl(2)

   integer, parameter :: InputIndex = INPUT_PREV
   integer, parameter :: StateIndex = STATE_PRED

   ErrStat = ErrID_None
   ErrMsg = ""

   !note: prescribed inputs are already set in both InputIndex=1 and InputIndex=2 so we can ignore them here
   ! Use current inputs to calculate CCSD in STATE_PRED
   call SteadyStateCCSD(m, caseData, InputIndex, T, ErrStat2, ErrMsg2)
   call SetErrStat(Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ! Store state accelerations in residual
   if (m%AM%iModBD > 0) then
      call ModD_PackAry(m%AM%Mod%Xfr(m%AM%iModBD)%x, m%Modules(m%AM%iModBD)%Lin%dx, m%AM%Residual)
   else if (m%AM%iModED > 0) then
      call ModD_PackAry(m%AM%Mod%Xfr(m%AM%iModED)%x, m%Modules(m%AM%iModED)%Lin%dx, m%AM%Residual)
   end if

   ! note that we don't need to calculate the inputs on more than p_FAST%NumBl_Lin blades because we are only using them to compute the SS_GetInputs
   call SteadyStateCalculatedInputs(m, InputIndex, T, ErrStat2, ErrMsg2) ! calculate new inputs and store in InputIndex=2
   call SetErrStat(Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ! Pack the output "residual vector" with these state derivatives and new inputs:
   call SS_GetInputs(m, m%AM%u2, T%p_FAST%UJacSclFact, InputIndex, StateIndex, T, ErrStat2, ErrMsg2)

   ! Store difference in inputs
   m%AM%Residual(m%AM%Mod%Vars%Nx + 1:) = m%AM%u1 - m%AM%u2

end subroutine SS_BuildResidual

!-------------------------------------------------------------------------------

!> SS_BD_InputSolve sets the blade load inputs required for BD.
subroutine SS_BD_InputSolve(m, InputIndex, T, ErrStat, ErrMsg)
   type(Glue_MiscVarType), intent(INOUT)     :: m           !< Miscellaneous variables
   integer(IntKi), intent(in)                :: InputIndex  !< Input index to transfer
   type(FAST_TurbineType), intent(INOUT)     :: T           !< Turbine type
   integer(IntKi), intent(OUT)               :: ErrStat     !< Error status of the operation
   character(*), intent(OUT)                 :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   character(*), parameter                   :: RoutineName = 'SS_BD_InputSolve'
   integer(IntKi)                            :: ErrStat2
   character(ErrMsgLen)                      :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ""

   call FAST_InputSolve(m%Modules(m%AM%iModBD), m%Modules, m%Mappings, InputIndex, T, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
end subroutine SS_BD_InputSolve

!> SS_BD_InputSolve_OtherBlades sets the blade-load ElastoDyn inputs from blade 1 to the other blades.
subroutine SS_BD_InputSolve_OtherBlades(m, InputIndex, T)
   type(Glue_MiscVarType), intent(INOUT)     :: m           !< Miscellaneous variables
   integer(IntKi), intent(in)                :: InputIndex  !< Input index to transfer
   type(FAST_TurbineType), intent(INOUT)     :: T           !< Turbine type
   integer(IntKi)                            :: j, k
   do k = 2, T%p_FAST%nBeams
      do j = 1, T%BD%Input(InputIndex, k)%DistrLoad%NNodes
         T%BD%Input(InputIndex, k)%DistrLoad%Force(:, j) = MATMUL(T%BD%Input(InputIndex, 1)%DistrLoad%Force(:, j), m%AM%HubOrientation(:, :, k))
         T%BD%Input(InputIndex, k)%DistrLoad%Moment(:, j) = MATMUL(T%BD%Input(InputIndex, 1)%DistrLoad%Moment(:, j), m%AM%HubOrientation(:, :, k))
      end do
   end do
end subroutine SS_BD_InputSolve_OtherBlades

!> This routine sets the blade load inputs required for ED.
subroutine SS_ED_InputSolve(m, InputIndex, T, ErrStat, ErrMsg)
   type(Glue_MiscVarType), intent(INOUT)     :: m           !< Miscellaneous variables
   integer(IntKi), intent(in)                :: InputIndex  !< Input index to transfer
   type(FAST_TurbineType), intent(INOUT)     :: T           !< Turbine type
   integer(IntKi), intent(OUT)               :: ErrStat     !< Error status of the operation
   character(*), intent(OUT)                 :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   character(*), parameter                   :: RoutineName = 'SS_ED_InputSolve'
   integer(IntKi)                            :: ErrStat2
   character(ErrMsgLen)                      :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ""

   call FAST_InputSolve(m%Modules(m%AM%iModED), m%Modules, m%Mappings, InputIndex, T, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
end subroutine SS_ED_InputSolve

!> SS_ED_InputSolve_OtherBlades sets the blade-load ElastoDyn inputs from blade 1 to the other blades.
subroutine SS_ED_InputSolve_OtherBlades(m, InputIndex, T)
   type(Glue_MiscVarType), intent(INOUT)     :: m           !< Miscellaneous variables
   integer(IntKi), intent(in)                :: InputIndex  !< Input index to transfer
   type(FAST_TurbineType), intent(INOUT)     :: T           !< Turbine type

   integer(IntKi)                            :: j, k

   associate (BladePtLoads => T%ED%Input(InputIndex)%BladePtLoads)
      do k = 2, size(BladePtLoads, 1)
         do j = 1, BladePtLoads(k)%NNodes
            BladePtLoads(k)%Force(:, j) = MATMUL(BladePtLoads(1)%Force(:, j), m%AM%HubOrientation(:, :, k))
            BladePtLoads(k)%Moment(:, j) = MATMUL(BladePtLoads(1)%Moment(:, j), m%AM%HubOrientation(:, :, k))
         end do
      end do
   end associate
end subroutine SS_ED_InputSolve_OtherBlades

!> SS_AD_InputSolve sets the blade-motion AeroDyn inputs for Blade 1.
subroutine SS_AD_InputSolve(m, InputIndex, T, ErrStat, ErrMsg)
   type(Glue_MiscVarType), intent(INOUT)     :: m           !< Miscellaneous variables
   integer(IntKi), intent(in)                :: InputIndex  !< Input index to transfer
   type(FAST_TurbineType), intent(INOUT)     :: T           !< Turbine type
   integer(IntKi), intent(OUT)               :: ErrStat     !< Error status of the operation
   character(*), intent(OUT)                 :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   character(*), parameter                   :: RoutineName = 'SS_AD_InputSolve'
   integer(IntKi)                            :: ErrStat2
   character(ErrMsgLen)                      :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Get blade motion inputs
   call FAST_InputSolve(m%Modules(m%AM%iModAD), m%Modules, m%Mappings, InputIndex, T, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ! Set prescribed values for first blade
   T%AD%Input(InputIndex)%rotors(1)%BladeMotion(1)%RotationVel = 0.0_ReKi
   T%AD%Input(InputIndex)%rotors(1)%BladeMotion(1)%TranslationAcc = 0.0_ReKi

end subroutine SS_AD_InputSolve

!> SS_AD_InputSolve_OtherBlades sets the blade-motion AeroDyn inputs.
subroutine SS_AD_InputSolve_OtherBlades(m, InputIndex, T)
   type(Glue_MiscVarType), intent(INOUT)     :: m           !< Miscellaneous variables
   integer(IntKi), intent(in)                :: InputIndex  !< Input index to transfer
   type(FAST_TurbineType), intent(INOUT)     :: T           !< Turbine type

   integer(IntKi)                            :: j, k

   associate (BladeMotion => T%AD%Input(InputIndex)%rotors(1)%BladeMotion)
      do k = 2, size(BladeMotion, 1)
         do j = 1, BladeMotion(k)%NNodes
            BladeMotion(k)%TranslationDisp(:, j) = MATMUL(BladeMotion(1)%TranslationDisp(:, j), m%AM%HubOrientation(:, :, k))
            BladeMotion(k)%Orientation(:, :, j) = MATMUL(BladeMotion(1)%Orientation(:, :, j), m%AM%HubOrientation(:, :, k))
            BladeMotion(k)%TranslationVel(:, j) = MATMUL(BladeMotion(1)%TranslationVel(:, j), m%AM%HubOrientation(:, :, k))
         end do
      end do
   end associate
end subroutine SS_AD_InputSolve_OtherBlades

subroutine SteadyStateCCSD(m, caseData, InputIndex, T, ErrStat, ErrMsg)
   type(Glue_MiscVarType), intent(INOUT)     :: m                   !< Miscellaneous variables
   type(AeroMapCase), intent(IN)             :: caseData            !< tsr, windSpeed, pitch, and rotor speed for this case
   integer(IntKi), intent(IN)                :: InputIndex          !< Index into input array
   type(FAST_TurbineType), intent(INOUT)     :: T                   !< Turbine type
   integer(IntKi), intent(OUT)               :: ErrStat             !< Error status
   character(*), intent(OUT)                 :: ErrMsg              !< Error message

   character(*), parameter                   :: RoutineName = 'SteadyStateCCSD'
   integer(IntKi)                            :: ErrStat2            ! temporary Error status of the operation
   character(ErrMsgLen)                      :: ErrMsg2             ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)                            :: i, k
   integer(IntKi)                            :: BldMeshNode
   real(R8Ki)                                :: Omega_Hub(3)
   real(R8Ki)                                :: position(3)
   real(R8Ki)                                :: omega_cross_position(3)

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Select based on which module is simulating the blades
   select case (T%p_FAST%CompElast)

   case (Module_ED) ! ElastoDyn

      call FAST_GetOP(m%Modules(m%AM%iModED), SS_t_global, InputIndex, STATE_CURR, T, ErrStat2, ErrMsg2, dx_op=m%Modules(m%AM%iModED)%Lin%dx)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   case (Module_BD) ! BeamDyn

      ! Set hub rotation speed
      Omega_Hub = [caseData%RotSpeed, 0.0_R8Ki, 0.0_R8Ki]

      ! TODO: Make this work for BeamDyn
      ! do K = 1, T%p_FAST%nBeams

      !    call BD_CalcContStateDeriv(SS_t_global, BD%Input(InputIndex, k), BD%p(k), BD%x(k, STATE_CURR), BD%xd(k, STATE_CURR), BD%z(k, STATE_CURR), &
      !                               BD%OtherSt(k, STATE_CURR), BD%m(k), BD%x(k, STATE_PRED), ErrStat2, ErrMsg2)
      !    call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      !    ! subtract xdot(y) here:
      !    ! note that this only works when the BldMotion mesh is on the FE nodes
      !    do i = 2, BD%p(k)%node_total ! the first node isn't technically a state
      !       BldMeshNode = BD%p(k)%NdIndx(i)
      !       position = BD%y(k)%BldMotion%Position(:, BldMeshNode) + BD%y(k)%BldMotion%TranslationDisp(:, BldMeshNode)
      !       omega_cross_position = cross_product(Omega_Hub, position)

      !       BD%x(k, STATE_PRED)%q(1:3, i) = BD%x(k, STATE_PRED)%q(1:3, i) - omega_cross_position
      !       BD%x(k, STATE_PRED)%q(4:6, i) = BD%x(k, STATE_PRED)%q(4:6, i) - Omega_Hub
      !       BD%x(k, STATE_PRED)%dqdt(1:3, i) = BD%x(k, STATE_PRED)%dqdt(1:3, i) - cross_product(Omega_Hub, omega_cross_position)
      !    end do

      ! end do

   end select

end subroutine SteadyStateCCSD

!----------------------------------------------------------------------------------------------------------------------------------
subroutine SteadyStateCalculatedInputs(m, InputIndex, T, ErrStat, ErrMsg)
   type(Glue_MiscVarType), intent(INOUT)     :: m                   !< Miscellaneous variables
   integer(IntKi), intent(IN)                :: InputIndex          !< Index into input array
   type(FAST_TurbineType), intent(INOUT)     :: T                   !< Turbine type
   integer(IntKi), intent(OUT)               :: ErrStat             !< Error status
   character(*), intent(OUT)                 :: ErrMsg              !< Error message

   character(*), parameter                   :: RoutineName = 'SteadyStateCalculatedInputs'
   integer(IntKi)                            :: ErrStat2            ! temporary Error status of the operation
   character(ErrMsgLen)                      :: ErrMsg2             ! temporary Error message if ErrStat /= ErrID_None

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Transfer motions to AeroDyn first
   call SS_AD_InputSolve(m, InputIndex, T, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ! Transfer loads to structural solver next
   if (m%AM%iModBD > 0) then
      call SS_BD_InputSolve(m, InputIndex, T, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   else if (m%AM%iModED > 0) then
      call SS_ED_InputSolve(m, InputIndex, T, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   end if

end subroutine SteadyStateCalculatedInputs

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine adds u_delta to the corresponding mesh field and scales it as appropriate
! subroutine Add_SteadyState_delta(p_FAST, y_FAST, u_delta, AD, ED, BD, MeshMapData)
!    !..................................................................................................................................
!    type(FAST_ParameterType), intent(IN) :: p_FAST           !< Glue-code simulation parameters
!    type(FAST_OutputFileType), intent(IN) :: y_FAST           !< Output variables for the glue code
!    real(R8Ki), intent(IN) :: u_delta(:)       !< The delta amount to add to the appropriate mesh fields
!    type(ElastoDyn_Data), intent(INOUT) :: ED               !< ElastoDyn data
!    type(BeamDyn_Data), intent(INOUT) :: BD               !< BeamDyn data
!    type(AeroDyn_Data), intent(INOUT) :: AD               !< AeroDyn data
!    type(FAST_ModuleMapType), intent(IN) :: MeshMapData      !< data for mapping meshes between modules

!    ! local variables
!    integer                                             :: n
!    integer                                             :: fieldIndx
!    integer                                             :: node
!    integer                                             :: indx, indx_last
!    integer                                             :: i, j, k
!    integer                                             :: nx, nStates

!    real(R8Ki)                                          :: orientation(3, 3)
!    real(R8Ki)                                          :: rotation(3, 3)

!    integer(IntKi)                                      :: ErrStat2
!    character(ErrMsgLen)                                :: ErrMsg2

!    nx = y_FAST%Lin%Glue%SizeLin(LIN_ContSTATE_COL)

!    ! structural code states:
!    if (p_FAST%CompElast == Module_ED) then
!       nStates = nx

!       do j = 1, nStates

!          do k = 1, ED%p%NActvDOF_Stride ! transfer these states to the other blades (this means that the original states MUST be set the same for all blades!!!)
!             indx = ED%p%DOFs%PS((j - 1)*ED%p%NActvDOF_Stride + k)

!             ED%x(STATE_CURR)%QT(indx) = ED%x(STATE_CURR)%QT(indx) + u_delta(j)
!             ED%x(STATE_CURR)%QDT(indx) = 0.0_R8Ki !ED%x( STATE_CURR)%QDT(indx)  + u_delta(j+nStates)
!          end do

!       end do

!    elseif (p_FAST%CompElast == Module_BD) then
!       nStates = nx/2

!       ! see BD's Perturb_x function:

!       do k = 1, p_FAST%nBeams
!          indx = 1
!          do i = 2, BD%p(k)%node_total
!             indx_last = indx + BD%p(k)%dof_node - 1
!             BD%x(k, STATE_CURR)%dqdt(:, i) = BD%x(k, STATE_CURR)%dqdt(:, i) + u_delta(nStates + indx:indx_last + nStates)
!             BD%x(k, STATE_CURR)%q(1:3, i) = BD%x(k, STATE_CURR)%q(1:3, i) + u_delta(indx:indx + 2)

!             ! w-m parameters
!             call BD_CrvMatrixR(BD%x(k, STATE_CURR)%q(4:6, i), rotation) ! returns the rotation matrix (transpose of DCM) that was stored in the state as a w-m parameter
!             orientation = transpose(rotation)

!             call PerturbOrientationMatrix(Orientation, Perturbations=u_delta(indx + 3:indx_last))

!             rotation = transpose(orientation)
!             call BD_CrvExtractCrv(rotation, BD%x(k, STATE_CURR)%q(4:6, i), ErrStat2, ErrMsg2) ! return the w-m parameters of the new orientation

!             indx = indx_last + 1
!          end do
!       end do
!    end if !CompElast

!    ! inputs:
!    ! we are at u_delta(nx+1 : end)
!    n = nx + 1
!    if (p_FAST%CompElast == Module_ED) then

!       do K = 1, p_FAST%NumBl_Lin !we don't need all blades here: SIZE(ED%Input(1)%BladePtLoads,1) ! Loop through all blades

!          do node = 1, ED%Input(1)%BladePtLoads(k)%NNodes
!             do fieldIndx = 1, 3
!                ED%Input(1)%BladePtLoads(k)%Force(fieldIndx, node) = ED%Input(1)%BladePtLoads(k)%Force(fieldIndx, node) + u_delta(n)*p_FAST%UJacSclFact
!                n = n + 1
!             end do
!          end do

!          do node = 1, ED%Input(1)%BladePtLoads(k)%NNodes
!             do fieldIndx = 1, 3
!                ED%Input(1)%BladePtLoads(k)%Moment(fieldIndx, node) = ED%Input(1)%BladePtLoads(k)%Moment(fieldIndx, node) + u_delta(n)*p_FAST%UJacSclFact
!                n = n + 1
!             end do
!          end do

!       end do

!       call SS_ED_InputSolve_OtherBlades(p_FAST, ED%Input(1), MeshMapData)

!    elseif (p_FAST%CompElast == Module_BD) then

!       do K = 1, p_FAST%NumBl_Lin !we don't need all blades here: p_FAST%nBeams ! Loop through all blades

!          do node = 1, BD%Input(1, k)%DistrLoad%NNodes
!             do fieldIndx = 1, 3
!                BD%Input(1, k)%DistrLoad%Force(fieldIndx, node) = BD%Input(1, k)%DistrLoad%Force(fieldIndx, node) + u_delta(n)*p_FAST%UJacSclFact
!                n = n + 1
!             end do
!          end do

!          do node = 1, BD%Input(1, k)%DistrLoad%NNodes
!             do fieldIndx = 1, 3
!                BD%Input(1, k)%DistrLoad%Moment(fieldIndx, node) = BD%Input(1, k)%DistrLoad%Moment(fieldIndx, node) + u_delta(n)*p_FAST%UJacSclFact
!                n = n + 1
!             end do
!          end do

!       end do

!       call SS_BD_InputSolve_OtherBlades(p_FAST, BD, MeshMapData, 1) ! 1 is for the input index (i.e., Input(1,Blades2-end)

!    end if !CompElast

!    ! AeroDyn
!    do k = 1, p_FAST%NumBl_Lin !we don't need all blades here: SIZE(AD%Input(1)%BladeMotion)
!       do node = 1, AD%Input(1)%rotors(1)%BladeMotion(k)%NNodes
!          do fieldIndx = 1, 3
!             AD%Input(1)%rotors(1)%BladeMotion(k)%TranslationDisp(fieldIndx, node) = AD%Input(1)%rotors(1)%BladeMotion(k)%TranslationDisp(fieldIndx, node) + u_delta(n)
!             n = n + 1
!          end do
!       end do

!       do node = 1, AD%Input(1)%rotors(1)%BladeMotion(k)%NNodes
!          call PerturbOrientationMatrix(AD%Input(1)%rotors(1)%BladeMotion(k)%Orientation(:, :, node), Perturbations=u_delta(n:n + 2))
!          n = n + 3
!       end do

!       do node = 1, AD%Input(1)%rotors(1)%BladeMotion(k)%NNodes
!          AD%Input(1)%rotors(1)%BladeMotion(k)%TranslationVel(:, node) = AD%Input(1)%rotors(1)%BladeMotion(k)%TranslationVel(:, node) + u_delta(n:n + 2)

!          n = n + 3
!       end do

!    end do

!    ! now update the inputs on other blades:
!    call SS_AD_InputSolve_OtherBlades(p_FAST, AD%Input(1), MeshMapData) ! transfer results from blade 1 to other blades

! end subroutine Add_SteadyState_delta

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine basically packs the relevant parts of the modules' inputs and states for use in the steady-state solver.
subroutine SS_GetInputs(m, u_vec, ScaleFactor, InputIndex, StateIndex, T, ErrStat, ErrMsg)
   type(Glue_MiscVarType), intent(inout)  :: m           !< Glue-code simulation parameters
   real(R8Ki), intent(inout)              :: u_vec(:)    !< Array of input packed values
   real(R8Ki), intent(in)                 :: ScaleFactor !< Jacobian scaling factor
   integer(IntKi), intent(in)             :: InputIndex  !< Input array index
   integer(IntKi), intent(in)             :: StateIndex  !< State array index
   type(FAST_TurbineType), intent(inout)  :: T           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat     !< Error status of the operation
   character(*), intent(out)              :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   character(*), parameter :: RoutineName = 'SolveSteadyState'
   integer(IntKi)          :: ErrStat2                  ! temporary Error status of the operation
   character(ErrMsgLen)    :: ErrMsg2                   ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)          :: i, j, k, ieMod, ieGbl
   integer(IntKi)          :: iMod(3), iVarMod(2), iVarGbl(2)

   iMod = [m%AM%iModED, m%AM%iModBD, m%AM%iModAD]

   ! Loop through modules
   do i = 1, size(iMod)

      ! Skip inactive modules
      if (iMod(i) == 0) cycle

      associate (ModData => m%Modules(iMod(i)))

         ! Get states and outputs
         call FAST_GetOP(ModData, SS_t_global, InputIndex, StateIndex, T, ErrStat2, ErrMsg2, &
                         u_op=ModData%Lin%u)
         if (Failed()) return

         ! Transfer selected input data from module to RHS based on Idx
         if (allocated(ModData%Lin%u)) then
            do j = 1, size(ModData%Vars%u)

               ! Get module and global variable indices from Idx, skip if not used
               if (.not. ModD_GetValLoc(m%AM%Mod%Xfr(iMod(i))%u, j, iVarMod, iVarGbl)) cycle

               ! Convert or store based on field type
               select case (ModData%Vars%u(j)%Field)

               case (FieldForce, FieldMoment)
                  ! If field is a force or moment, scale by scale factor
                  u_vec(iVarGbl(1):iVarGbl(2)) = ModData%Lin%u(iVarMod(1):iVarMod(2))/ScaleFactor

               case (FieldOrientation)
                  ! Convert orientations to rotation vectors
                  ieMod = iVarMod(1)
                  ieGbl = iVarGbl(1)
                  do k = 1, ModData%Vars%u(j)%Nodes
                     u_vec(ieGbl:ieGbl + 2) = -quat_to_rvec(ModData%Lin%u(ieMod:ieMod + 2))
                     ieMod = ieMod + 3
                     ieGbl = ieGbl + 3
                  end do

               case default
                  u_vec(iVarGbl(1):iVarGbl(2)) = ModData%Lin%u(iVarMod(1):iVarMod(2))
               end select

            end do
         end if

      end associate

   end do

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine SS_GetInputs

!----------------------------------------------------------------------------------------------------------------------------------
subroutine SetPrescribedInputs(caseData, p_FAST, y_FAST, m_FAST, ED, BD, AD)
   type(AeroMapCase), intent(in)             :: caseData            !< tsr, windSpeed, pitch, and rotor speed for this case
   type(FAST_ParameterType), intent(in)      :: p_FAST              !< Parameters for the glue code
   type(FAST_OutputFileType), intent(inout)  :: y_FAST              !< Output variables for the glue code
   type(FAST_MiscVarType), intent(inout)     :: m_FAST              !< Miscellaneous variables

   type(ElastoDyn_Data), intent(inout)       :: ED                  !< ElastoDyn data
   type(BeamDyn_Data), intent(inout)         :: BD                  !< BeamDyn data
   type(AeroDyn_Data), intent(inout)         :: AD                  !< AeroDyn data

   integer(IntKi)                              :: k
   real(R8Ki)                                  :: theta(3)

   ! Set prescribed inputs for all of the modules in the steady-state solve

   ED%Input(1)%TwrAddedMass = 0.0_ReKi
   ED%Input(1)%PtfmAddedMass = 0.0_ReKi

   ED%Input(1)%TowerPtLoads%Force = 0.0
   ED%Input(1)%TowerPtLoads%Moment = 0.0
   ED%Input(1)%NacelleLoads%Force = 0.0
   ED%Input(1)%NacelleLoads%Moment = 0.0
   ED%Input(1)%HubPtLoad%Force = 0.0      ! these are from BD, but they don't affect the ED calculations for aeromaps, so set them to 0
   ED%Input(1)%HubPtLoad%Moment = 0.0     ! these are from BD, but they don't affect the ED calculations for aeromaps, so set them to 0

   ED%Input(1)%BlPitchCom = caseData%Pitch
   ED%Input(1)%YawMom = 0.0
   ED%Input(1)%HSSBrTrqC = 0.0
   ED%Input(1)%GenTrq = 0.0

   ! BeamDyn
   if (p_FAST%CompElast == Module_BD) then

      !CALL ED_CalcOutput( 0.0_DbKi, ED%Input(1), ED%p, ED%x(STATE_CURR), ED%xd(STATE_CURR), ED%z(STATE_CURR), ED%OtherSt(STATE_CURR), ED%y, ED%m, ErrStat2, ErrMsg2 )
      !   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )

      do k = 1, p_FAST%nBeams
         BD%Input(1, k)%RootMotion%TranslationDisp = 0.0_ReKi

         theta = EulerExtract(BD%Input(1, k)%RootMotion%RefOrientation(:, :, 1))
         theta(3) = -caseData%Pitch
         BD%Input(1, k)%RootMotion%Orientation(:, :, 1) = EulerConstruct(theta)

         BD%Input(1, k)%RootMotion%RotationVel(1, 1) = caseData%RotSpeed !BD%Input(1,k)%RootMotion%RotationVel = ED%y_interp%BladeRootMotion(k)%RotationVel
         BD%Input(1, k)%RootMotion%RotationVel(2:3, 1) = 0.0_ReKi

         BD%Input(1, k)%RootMotion%TranslationVel(:, 1) = cross_product(BD%Input(1, k)%RootMotion%RotationVel(:, 1), BD%Input(1, k)%RootMotion%Position(:, 1) - AD%Input(1)%rotors(1)%HubMotion%Position(:, 1)) ! ED%y_interp%BladeRootMotion(k)%TranslationVel
         BD%Input(1, k)%RootMotion%TranslationAcc(:, 1) = cross_product(BD%Input(1, k)%RootMotion%RotationVel(:, 1), BD%Input(1, k)%RootMotion%TranslationVel(:, 1)) ! ED%y_interp%BladeRootMotion(k)%TranslationAcc

         BD%Input(1, k)%RootMotion%RotationAcc = 0.0_ReKi
      end do ! k=p_FAST%nBeams

   end if  ! BeamDyn
   !BeamDyn's first "state" is not actually the state. So, do we need to do something with that?????

   !AeroDyn
   !note: i'm skipping the (unused) TowerMotion mesh
   AD%Input(1)%rotors(1)%HubMotion%TranslationDisp = 0.0
   AD%Input(1)%rotors(1)%HubMotion%Orientation = AD%Input(1)%rotors(1)%HubMotion%RefOrientation
   AD%Input(1)%rotors(1)%HubMotion%RotationVel(1, :) = caseData%RotSpeed
   AD%Input(1)%rotors(1)%HubMotion%RotationVel(2:3, :) = 0.0_ReKi

   do k = 1, size(AD%Input(1)%rotors(1)%BladeRootMotion, 1)
      theta = EulerExtract(AD%Input(1)%rotors(1)%BladeRootMotion(k)%RefOrientation(:, :, 1))
      theta(3) = -caseData%Pitch
      AD%Input(1)%rotors(1)%BladeRootMotion(k)%Orientation(:, :, 1) = EulerConstruct(theta) !AD%Input(1)%BladeRootMotion(k)%RefOrientation

      AD%Input(1)%rotors(1)%BladeMotion(k)%RotationVel = 0.0_ReKi
      !AD%Input(1)%rotors(1)%BladeMotion(k)%RotationAcc = 0.0_ReKi
      AD%Input(1)%rotors(1)%BladeMotion(k)%TranslationAcc = 0.0_ReKi
   end do

   ! Set FlowField information -- AD calculates everything from the data stored in the FlowField pointer
   AD%p%FlowField%Uniform%VelH(:) = caseData%WindSpeed
   AD%p%FlowField%Uniform%LinShrV(:) = 0.0_ReKi
   AD%p%FlowField%Uniform%AngleH(:) = 0.0_ReKi
   AD%p%FlowField%PropagationDir = 0.0_ReKi

   AD%Input(1)%rotors(1)%UserProp = 0.0_ReKi

end subroutine SetPrescribedInputs

end module
