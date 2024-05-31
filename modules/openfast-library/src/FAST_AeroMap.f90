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

use FAST_Types
use FAST_ModTypes
use FAST_Funcs
use FAST_Idx

use FAST_Subs

implicit none

! Define array of module IDs used in AeroMap
integer(IntKi), parameter :: AeroMapModIDs(*) = [Module_ED, Module_BD, Module_AD]

real(DbKi), parameter    :: SS_t_global = 0.0_DbKi
real(DbKi), parameter    :: UJacSclFact_x = 1.0d3

logical, parameter    :: output_debugging = .false.

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

   call FAST_AeroMapDriver(Turbine, ErrStat, ErrMsg)
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

subroutine FAST_AeroMapDriver(Turbine, ErrStat, ErrMsg)
   use InflowWind_IO, only: IfW_SteadyFlowField_Init
   type(FAST_TurbineType), intent(inout)  :: Turbine        !< all data for one instance of a turbine
   character(*), intent(out)              :: ErrMsg         !< Error message if ErrStat /= ErrID_None
   integer(IntKi), intent(out)            :: ErrStat        !< Error status of the operation

   character(*), parameter                :: RoutineName = 'FAST_AeroMapDriver'
   character(ErrMsgLen)                   :: ErrMsg2
   integer(IntKi)                         :: ErrStat2
   logical, parameter                     :: CompAeroMaps = .true.
   real(DbKi), parameter                  :: t_initial = 0.0_DbKi
   integer(IntKi), allocatable            :: modIDs(:), modIdx(:), iModOrder(:)
   integer(IntKi)                         :: i
   integer(IntKi)                         :: JacSize
   real(R8Ki), allocatable                :: Jmat(:, :)
   integer(IntKi), allocatable            :: JacPivot(:)
   type(VarsIdxType)                      :: AeroMapIdx

   ErrStat = ErrID_None
   ErrMsg = ''

   !----------------------------------------------------------------------------
   ! Initialization
   !----------------------------------------------------------------------------

   Turbine%TurbID = 1

   ! Standard Turbine initialization
   call FAST_InitializeAll(t_initial, Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
                           Turbine%ED, Turbine%BD, Turbine%SrvD, Turbine%AD14, Turbine%AD, &
                           Turbine%ExtLd, Turbine%IfW, Turbine%ExtInfw, Turbine%SC_DX, &
                           Turbine%SeaSt, Turbine%HD, Turbine%SD, Turbine%ExtPtfm, Turbine%MAP, &
                           Turbine%FEAM, Turbine%MD, Turbine%Orca, Turbine%IceF, Turbine%IceD, &
                           Turbine%MeshMapData, CompAeroMaps, ErrStat, ErrMsg)

   ! Initialize steady flow field in AeroDyn
   call IfW_SteadyFlowField_Init(Turbine%AD%p%FlowField, &
                                 RefHt=100.0_ReKi, HWindSpeed=8.0_ReKi, PLExp=0.0_ReKi, &
                                 ErrStat=ErrStat2, ErrMsg=ErrMsg2); if (Failed()) return

   ! Get indices of Mods that are used by Aero Mapping
   call GetModuleOrder(Turbine%y_FAST%Modules, AeroMapModIDs, iModOrder)

   ! Loop through module data indices
   do i = 1, size(iModOrder)

      ! Copy current state to predicted state
      call FAST_CopyStates(Turbine%y_FAST%Modules(iModOrder(i)), Turbine, STATE_CURR, STATE_PRED, MESH_NEWCOPY, &
                           ErrStat2, ErrMsg2); if (Failed()) return

      ! Copy inputs to second index
      call FAST_CopyInput(Turbine%y_FAST%Modules(iModOrder(i)), Turbine, 0.0_DbKi, iSrc=1, iDst=2, CtrlCode=MESH_NEWCOPY, &
                           ErrStat=ErrStat2, ErrMsg=ErrMsg2); if (Failed()) return
   end do

   ! Generate index for variables with AeroMap flag
   call Idx_Init(Turbine%y_FAST%Modules, iModOrder, AeroMapIdx, VF_AeroMap, ErrStat2, ErrMsg2); if (Failed()) return

   ! Jacobian size is number of states plus number of inputs
   JacSize = AeroMapIdx%Nx + AeroMapIdx%Nu

   ! Allocate Jacobian matrix
   call AllocAry(Jmat, JacSize, JacSize, 'Jmat', ErrStat2, ErrMsg2); if (Failed()) return

   ! Allocate Jacobian pivot vector
   call AllocAry(JacPivot, JacSize, 'Pivot array for Jacobian LU decomposition', ErrStat2, ErrMsg2); if (Failed()) return

   !----------------------------------------------------------------------------
   ! Calculate steady-state solutions:
   !----------------------------------------------------------------------------

   call FAST_SteadyState(Turbine%y_FAST%Modules, iModOrder, Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
                         Turbine, ErrStat2, ErrMsg2); if (Failed()) return

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that calls FAST_Solution for one instance of a Turbine data structure. This is a separate subroutine so that the FAST
!! driver programs do not need to change or operate on the individual module level.
subroutine FAST_SteadyState(Mods, ModOrder, p_FAST, y_FAST, m_FAST, T, ErrStat, ErrMsg)
   type(ModDataType), intent(inout)          :: Mods(:)
   integer(IntKi), intent(in)                :: ModOrder(:)
   type(FAST_ParameterType), intent(IN)      :: p_FAST              !< Parameters for the glue code
   type(FAST_OutputFileType), intent(INOUT)  :: y_FAST              !< Output variables for the glue code
   type(FAST_MiscVarType), intent(INOUT)     :: m_FAST
   type(FAST_TurbineType), intent(inout)     :: T              !< all data for one instance of a turbine
   integer(IntKi), intent(out)               :: ErrStat        !< Error status of the operation
   character(*), intent(out)                 :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   integer(IntKi)                         :: n_case         !< loop counter
   real(DbKi)                             :: n_global
   real(ReKi), allocatable                :: UnusedAry(:)
   real(R8Ki), allocatable                :: Jmat(:, :)
   type(FAST_SS_CaseType)                 :: CaseData       ! tsr, windSpeed, pitch, and rotor speed for this case
   type(FAST_SS_CaseType)                 :: caseData_try2  ! tsr, windSpeed, pitch, and rotor speed for this case (to try a different operating point first)

   character(*), parameter                :: RoutineName = 'FAST_SteadyState'
   integer(IntKi)                         :: ErrStat2
   character(ErrMsgLen)                   :: ErrMSg2
   integer(IntKi)                         :: NStatus
   character(MaxWrScrLen), parameter      :: BlankLine = " "

   ErrStat = ErrID_None
   ErrMsg = ""

   ! how often do we inform the user which case we are on?
   NStatus = min(100, p_FAST%NumSSCases/100 + 1) ! at least 100 every 100 cases or 100 times per simulation
   call WrScr(NewLine)

   ! Loop through Aero Map cases
   do n_case = 1, p_FAST%NumSSCases

      ! If status should be written to screen
      if (n_case == 1 .or. n_case == p_FAST%NumSSCases .or. mod(n_case, NStatus) == 0) then
         call WrOver(' Case '//trim(Num2LStr(n_case))//' of '//trim(Num2LStr(p_FAST%NumSSCases)))
      end if

      ! Populate case data
      if (p_FAST%WindSpeedOrTSR == 1) then
         CaseData%WindSpeed = p_FAST%WS_TSR(n_case)
         CaseData%TSR = p_FAST%RotSpeed(n_case)*T%AD%p%rotors(1)%BEMT%rTipFixMax/CaseData%WindSpeed
      else
         CaseData%TSR = p_FAST%WS_TSR(n_case)
         CaseData%WindSpeed = p_FAST%RotSpeed(n_case)*T%AD%p%rotors(1)%BEMT%rTipFixMax/CaseData%TSR
      end if
      CaseData%Pitch = p_FAST%Pitch(n_case)
      CaseData%RotSpeed = p_FAST%RotSpeed(n_case)

      ! Call steady-state solve for this pitch and rotor speed
      call SolveSteadyState(Mods, ModOrder, caseData, Jmat, p_FAST, y_FAST, m_FAST, T%ED, T%BD, T%AD, T%MeshMapData, T, ErrStat2, ErrMsg2)

      ! we didn't converge; let's try a different operating point and see if that helps:
      if (ErrStat2 >= ErrID_Severe) then

         ! Create copy of case data for second attempt
         caseData_try2 = CaseData

         ! Modify pitch, TSR, and WindSpeed
         caseData_try2%Pitch = caseData_try2%Pitch*0.5_ReKi
         caseData_try2%TSR = caseData_try2%TSR*0.5_ReKi
         caseData_try2%WindSpeed = caseData_try2%WindSpeed*0.5_ReKi

         ! Write message about retrying case
         call WrScr('Retrying case '//trim(Num2LStr(n_case))//', first trying to get a better initial guess. Average error is '// &
                    trim(Num2LStr(y_FAST%DriverWriteOutput(SS_Indx_Err)))//'.')

         call SolveSteadyState(Mods, ModOrder, caseData_try2, Jmat, p_FAST, y_FAST, m_FAST, T%ED, T%BD, T%AD, T%MeshMapData, T, ErrStat2, ErrMsg2)

         ! if that worked, try the real case again:
         if (ErrStat2 < AbortErrLev) then
            call SolveSteadyState(Mods, ModOrder, caseData, Jmat, p_FAST, y_FAST, m_FAST, T%ED, T%BD, T%AD, T%MeshMapData, T, ErrStat2, ErrMsg2)
            call WrOver(BlankLine)
         end if

      end if

      if (ErrStat2 > ErrID_None) then
         ErrMsg2 = trim(ErrMsg2)//" case "//trim(Num2LStr(n_case))// &
                   ' (tsr='//trim(Num2LStr(CaseData%tsr))// &
                   ', wind speed='//trim(Num2LStr(CaseData%windSpeed))//' m/s'// &
                   ', pitch='//trim(num2lstr(CaseData%pitch*R2D))//' deg'// &
                   ', rotor speed='//trim(num2lstr(CaseData%RotSpeed*RPS2RPM))//' rpm)'
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
end subroutine FAST_SteadyState

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine performs the Input-Output solve for the steady-state solver.
!! Note that this has been customized for the physics in the problems and is not a general solution.
subroutine SolveSteadyState(Mods, ModOrder, caseData, Jmat, p_FAST, y_FAST, m_FAST, ED, BD, AD, MeshMapData, Turbine, ErrStat, ErrMsg)
   type(ModDataType), intent(inout)          :: Mods(:)
   integer(IntKi), intent(in)                :: ModOrder(:)
   type(FAST_SS_CaseType), intent(IN)        :: caseData       !< tsr, windSpeed, pitch, and rotor speed for this case
   real(R8Ki), intent(INOUT)                 :: Jmat(:, :)     !< temporary storage space for jacobian matrix
   type(FAST_ParameterType), intent(IN)      :: p_FAST         !< Glue-code simulation parameters
   type(FAST_OutputFileType), intent(INOUT)  :: y_FAST         !< Glue-code output file values
   type(FAST_MiscVarType), intent(INOUT)     :: m_FAST         !< Miscellaneous variables
   type(ElastoDyn_Data), intent(INOUT)       :: ED             !< ElastoDyn data
   type(BeamDyn_Data), intent(INOUT)         :: BD             !< BeamDyn data
   type(AeroDyn_Data), intent(INOUT)         :: AD             !< AeroDyn data
   type(FAST_TurbineType), intent(inout)     :: Turbine
   type(FAST_ModuleMapType), intent(INOUT)   :: MeshMapData    !< data for mapping meshes between modules
   integer(IntKi), intent(OUT)               :: ErrStat        !< Error status of the operation
   character(*), intent(OUT)                 :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   character(*), parameter :: RoutineName = 'SolveSteadyState'
   integer(IntKi)                                    :: ErrStat2                  ! temporary Error status of the operation
   character(ErrMsgLen)                              :: ErrMsg2                   ! temporary Error message if ErrStat /= ErrID_None

   !bjj: store these so that we don't reallocate every time?
   real(R8Ki)                                        :: u(p_FAST%SizeJac_Opt1(1))   ! size of loads/accelerations passed between the 6 modules
   real(R8Ki)                                        :: u_delta(p_FAST%SizeJac_Opt1(1))   ! size of loads/accelerations passed between the 6 modules
   real(R8Ki)                                        :: Fn_U_Resid(p_FAST%SizeJac_Opt1(1))   ! Residual of U
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
   ED%x(STATE_CURR)%QDT(p_FAST%GearBox_Index) = caseData%RotSpeed

   call SteadyStatePrescribedInputs(caseData, p_FAST, y_FAST, m_FAST, ED, BD, AD)
   call CopyStatesInputs(p_FAST, ED, BD, AD, ErrStat2, ErrMsg2, MESH_UPDATECOPY) ! COPY the inputs to the temp copy (so we get updated input values)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

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

      ! Loop through modules in order
      do i = 1, size(ModOrder)
         associate (ModData => Mods(ModOrder(i)))

            !-------------------------------------------------------------------
            ! ElastoDyn / BeamDyn CalcOutput
            !-------------------------------------------------------------------

            ! If ElastoDyn blades and module is ED or BeamDyn Blades and module is BD and 1st blade, calculate output
            if (((p_FAST%CompElast == Module_ED) .and. (ModData%ID == Module_ED)) .or. &
                ((p_FAST%CompElast == Module_BD) .and. (ModData%ID == Module_BD) .and. (ModData%Ins == 1))) then

               call FAST_CalcOutput(ModData, m_FAST%ModLin%Mappings, SS_t_global, STATE_CURR, Turbine, ErrStat2, ErrMsg2)
               call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            end if

            !-------------------------------------------------------------------
            ! AeroDyn InputSolve
            !-------------------------------------------------------------------

            ! If module is AD (assumes AD comes after ED/BD in ModOrder)
            if (ModData%ID == Module_AD) then

               ! If first iteration 
               if (K == 0) then

                  ! Perform AeroDyn input solve to get initial guess from structural module
                  ! (this ensures that the pitch is accounted for in the fixed aero-map solve:):
                  call FAST_InputSolve(ModData, Mods, m_FAST%ModLin%Mappings, Turbine, ErrStat2, ErrMsg2)
                  call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

                  ! Prescribe AeroDyn blade motion inputs on first blade:
                  Turbine%AD%u%rotors(1)%BladeMotion(1)%RotationVel = 0.0_ReKi
                  Turbine%AD%u%rotors(1)%BladeMotion(1)%TranslationAcc = 0.0_ReKi

                  ! Initialize AeroDyn blade motion from blade 1 to remaining blades
                  ! adjusting for hub orientation
                  do k = 2, size(Turbine%AD%u%rotors(1)%BladeMotion)
                     do j = 1, Turbine%AD%u%rotors(1)%BladeMotion(k)%NNodes
                        Turbine%AD%u%rotors(1)%BladeMotion(k)%TranslationDisp(:, j) = matmul(Turbine%AD%u%rotors(1)%BladeMotion(1)%TranslationDisp(:, j), MeshMapData%HubOrient(:, :, k))
                        Turbine%AD%u%rotors(1)%BladeMotion(k)%Orientation(:, :, j) = matmul(Turbine%AD%u%rotors(1)%BladeMotion(1)%Orientation(:, :, j), MeshMapData%HubOrient(:, :, k))
                        Turbine%AD%u%rotors(1)%BladeMotion(k)%TranslationVel(:, j) = matmul(Turbine%AD%u%rotors(1)%BladeMotion(1)%TranslationVel(:, j), MeshMapData%HubOrient(:, :, k))
                     end do
                  end do

                  !----------------------------------------------------------------------------------------------------
                  ! set up x-u vector, using local initial guesses:
                  !----------------------------------------------------------------------------------------------------
                  call Create_SS_Vector(p_FAST, y_FAST, u, AD, ED, BD, 1, STATE_CURR)

               end if   ! K == 0

               ! Calculate AeroDyn Output
               call FAST_CalcOutput(ModData, m_FAST%ModLin%Mappings, SS_t_global, STATE_CURR, Turbine, ErrStat2, ErrMsg2)
               call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               if (ErrStat >= AbortErrLev) then
                  call ResetInputsAndStates()
                  return
               end if

            end if   ! ModData%ID == Module_AD

         end associate
      end do

      ! If iteration is at or above maximum iteration, exit loop
      if (K >= MaxIter) exit

      !-------------------------------------------------------------------------------------------------
      ! Calculate residual and the Jacobian:
      ! (note that we don't want to change module%Input(1), here)
      ! Also, the residual uses values from y_FAST, so do this before calculating the jacobian
      !-------------------------------------------------------------------------------------------------

      call SteadyStateSolve_Residual(caseData, p_FAST, y_FAST, m_FAST, ED, BD, AD, MeshMapData, u, Fn_U_Resid, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) then
         call ResetInputsAndStates()
         return
      end if

      ! If Jacobian needs to be recalculated
      if (mod(K, p_FAST%N_UJac) == 0) then


         call FormSteadyStateJacobian(caseData, Jmat, p_FAST, y_FAST, m_FAST, ED, BD, AD, MeshMapData, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

         call Precondition_Jmat(p_FAST, y_FAST, Jmat)

         ! Get the LU decomposition of this matrix using a LAPACK routine:
         ! The result is of the form Jmat = P * L * U

         call LAPACK_getrf(M=size(Jmat, 1), N=size(Jmat, 2), &
                           A=Jmat, IPIV=MeshMapData%Jacobian_pivot, &
                           ErrStat=ErrStat2, ErrMsg=ErrMsg2)
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

      u_delta = -Fn_U_Resid
      call LAPACK_getrs(TRANS="N", N=SIZE(Jmat, 1), A=Jmat, &
                        IPIV=MeshMapData%Jacobian_pivot, B=u_delta, ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return

      !-------------------------------------------------------------------------
      ! check for error, update inputs if necessary, and iterate again
      !-------------------------------------------------------------------------
      err_prev = err
      err = DOT_PRODUCT(u_delta, u_delta)
      y_FAST%DriverWriteOutput(SS_Indx_Err) = sqrt(err)/p_FAST%SizeJac_Opt1(1)

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
      ! modify inputs and states for next iteration
      !-------------------------------------------------------------------------
      if (err > err_prev) then
         u_delta = u_delta*reduction_factor ! don't take a full step if we're getting farther from the solution!
         err_prev = err_prev*reduction_factor
      end if

      call Add_SteadyState_delta(p_FAST, y_FAST, u_delta, AD, ED, BD, MeshMapData)

      !u = u + u_delta
      call Create_SS_Vector(p_FAST, y_FAST, u, AD, ED, BD, 1, STATE_CURR)

      K = K + 1
      y_FAST%DriverWriteOutput(SS_Indx_Iter) = k

   end do ! K

   if (p_FAST%CompElast == Module_BD) then
      ! this doesn't actually get the correct hub point load from BD, but we'll get some outputs:
      call ED_CalcOutput(SS_t_global, ED%Input(1), ED%p, ED%x(STATE_CURR), ED%xd(STATE_CURR), ED%z(STATE_CURR), ED%OtherSt(STATE_CURR), ED%y, ED%m, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   end if

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
                  BD%Input(1, k)%DistrLoad%Force = 0.0_ReKi
                  BD%Input(1, k)%DistrLoad%Moment = 0.0_ReKi
               end do

            end if

            call Create_SS_Vector(p_FAST, y_FAST, u, AD, ED, BD, 1, STATE_CURR)     ! find the values we have been modifying (in u... continuous states and inputs)
            call Add_SteadyState_delta(p_FAST, y_FAST, -u, AD, ED, BD, MeshMapData) ! and reset them to 0 (by adding -u)

         end if
      end if
   end subroutine ResetInputsAndStates

end subroutine SolveSteadyState

!----------------------------------------------------------------------------------------------------------------------------------
subroutine SteadyStatePrescribedInputs(caseData, p_FAST, y_FAST, m_FAST, ED, BD, AD)
   type(FAST_SS_CaseType), intent(IN) :: caseData            !< tsr, windSpeed, pitch, and rotor speed for this case
   type(FAST_ParameterType), intent(IN) :: p_FAST              !< Parameters for the glue code
   type(FAST_OutputFileType), intent(INOUT) :: y_FAST              !< Output variables for the glue code
   type(FAST_MiscVarType), intent(INOUT) :: m_FAST              !< Miscellaneous variables

   type(ElastoDyn_Data), intent(INOUT) :: ED                  !< ElastoDyn data
   type(BeamDyn_Data), intent(INOUT) :: BD                  !< BeamDyn data
   type(AeroDyn_Data), intent(INOUT) :: AD                  !< AeroDyn data

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

end subroutine SteadyStatePrescribedInputs

end module
