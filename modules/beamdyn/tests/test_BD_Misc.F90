module test_BD_Misc

use test_tools
use BeamDyn_Subs
use BeamDyn

implicit none

private
public :: test_BD_Misc_suite

contains

!> Collect all exported unit tests
subroutine test_BD_Misc_suite(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
               new_unittest("test_BD_DistrLoadCopy", test_BD_DistrLoadCopy), &
               new_unittest("test_BD_InputGlobalLocal", test_BD_InputGlobalLocal), &
               new_unittest("test_BD_GravityForce", test_BD_GravityForce), &
               new_unittest("test_BD_QPData_mEta_rho", test_BD_QPData_mEta_rho) &
               ]
end subroutine

subroutine test_BD_DistrLoadCopy(error)
   type(error_type), allocatable, intent(out) :: error

   ! branches to test
   ! - the 2D array is correctly stored in the 3D array

   integer                    :: i, j
   type(BD_ParameterType)     :: parametertype
   type(BD_InputType)         :: inputtype
   type(BD_MiscVarType)       :: miscvartype
   integer(IntKi)             :: ErrStat
   character                  :: ErrMsg
   character(1024)            :: testname

   ! --------------------------------------------------------------------------
   testname = "static simple beam under gravity:"

   ! build the parametertype, inputtype, and miscvartype
   parametertype = simpleParameterType(1, 16, 16, 0, 1)
   miscvartype = simpleMiscVarType(parametertype%nqp, parametertype%dof_node, parametertype%elem_total, parametertype%nodes_per_elem)
   inputtype = simpleInputType(parametertype%nqp, parametertype%elem_total)

   call BD_DistrLoadCopy(parametertype, inputtype, miscvartype)

   do j = 1, parametertype%elem_total
      do i = 1, parametertype%nqp
         call check_array(error, real([9 * (j - 1) + 3 * (i - 1) + 1, 9 * (j - 1) + 3 * (i - 1) + 2, 9 * (j - 1) + 3 * (i - 1) + 3], R8Ki), miscvartype%DistrLoad_QP(1:3, i, j)); if (allocated(error)) return
         call check_array(error, real([-9 * (j - 1) - 3 * (i - 1) - 1, -9 * (j - 1) - 3 * (i - 1) - 2, -9 * (j - 1) - 3 * (i - 1) - 3], R8Ki), miscvartype%DistrLoad_QP(4:6, i, j)); if (allocated(error)) return
      end do
   end do

   call BD_DestroyParam(parametertype, ErrStat, ErrMsg)
end subroutine

subroutine test_BD_InputGlobalLocal(error)
   type(error_type), allocatable, intent(out) :: error

   ! branches to test
   ! - a simple rotation does the rotation

   ! Check the following quanities are actually rotated
    !!  1 Displacements                -> u%RootMotion%TranslationDisp(:,1)
    !!  2 Linear/Angular velocities    -> u%RootMotion%TranslationVel(:,1), u%RootMotion%RotationVel(:,1)
    !!  3 Linear/Angular accelerations -> u%RootMotion%TranslationAcc(:,1), u%RootMotion%RotationAcc(:,1)
    !!  4 Point forces/moments         -> u%PointLoad%Force(1:3,i), u%PointLoad%Moment(1:3,i)
    !!  5 Distributed forces/moments   -> u%DistrLoad%Force(1:3,i), u%DistrLoad%Moment(1:3,i)

   ! Verify the DCM is transposed
    !! u%RootMotion%Orientation(:,:,1)

   integer                    :: i, totalnodes
   type(BD_ParameterType)     :: parametertype
   type(BD_OtherStateType)    :: otherstate
   type(BD_InputType)         :: inputtype
   real(BDKi), dimension(3)   :: vectorInit, vectorAfterRotation, rotationaxis
   integer(IntKi)             :: ErrStat
   character                  :: ErrMsg
   character(1024)            :: testname

   ! --------------------------------------------------------------------------
   testname = "test_BD_InputGlobalLocal:"

   totalnodes = 2
   vectorInit = [0.0, 0.0, 1.0]
   vectorAfterRotation = [0.0, 0.0, -1.0]
   rotationaxis = [1.0, 0.0, 0.0]

   ! build the parameter type
   parametertype%node_total = totalnodes
   otherstate = simpleOtherState()
   otherstate%GlbRot = calcRotationMatrix(real(Pi, BDKi), rotationaxis)

   ! build the inputs
   call AllocAry(inputtype%RootMotion%TranslationDisp, 3, 1, 'TranslationDisp', ErrStat, ErrMsg)
   call AllocAry(inputtype%RootMotion%TranslationVel, 3, 1, 'TranslationVel', ErrStat, ErrMsg)
   call AllocAry(inputtype%RootMotion%RotationVel, 3, 1, 'RotationVel', ErrStat, ErrMsg)
   call AllocAry(inputtype%RootMotion%TranslationAcc, 3, 1, 'TranslationAcc', ErrStat, ErrMsg)
   call AllocAry(inputtype%RootMotion%RotationAcc, 3, 1, 'RotationAcc', ErrStat, ErrMsg)
   inputtype%RootMotion%TranslationDisp(:, 1) = vectorInit
   inputtype%RootMotion%TranslationVel(:, 1) = vectorInit
   inputtype%RootMotion%RotationVel(:, 1) = vectorInit
   inputtype%RootMotion%TranslationAcc(:, 1) = vectorInit
   inputtype%RootMotion%RotationAcc(:, 1) = vectorInit

   call AllocAry(inputtype%PointLoad%Force, 3, totalnodes, 'PointLoad%Force', ErrStat, ErrMsg)
   call AllocAry(inputtype%PointLoad%Moment, 3, totalnodes, 'PointLoad%Moment', ErrStat, ErrMsg)
   do i = 1, parametertype%node_total
      inputtype%PointLoad%Force(1:3, i) = vectorInit
      inputtype%PointLoad%Moment(1:3, i) = vectorInit
   end do

   inputtype%DistrLoad%Nnodes = totalnodes
   call AllocAry(inputtype%DistrLoad%Force, 3, totalnodes, 'DistrLoad%Force', ErrStat, ErrMsg)
   call AllocAry(inputtype%DistrLoad%Moment, 3, totalnodes, 'DistrLoad%Moment', ErrStat, ErrMsg)
   do i = 1, inputtype%DistrLoad%Nnodes
      inputtype%DistrLoad%Force(1:3, i) = vectorInit
      inputtype%DistrLoad%Moment(1:3, i) = vectorInit
   end do

   call AllocAry(inputtype%RootMotion%Orientation, 3, 3, totalnodes, 'RootMotion%Orientation', ErrStat, ErrMsg)
   inputtype%RootMotion%Orientation(:, :, 1) = otherstate%GlbRot

   ! call the subroutine to test
   call BD_InputGlobalLocal(parametertype, otherstate, inputtype)

   ! test the values
   call check_array(error, vectorAfterRotation, real(inputtype%RootMotion%TranslationDisp(:, 1), BDKi), testname, tolerance); if (allocated(error)) return
   call check_array(error, vectorAfterRotation, real(inputtype%RootMotion%TranslationVel(:, 1), BDKi), testname, tolerance); if (allocated(error)) return
   call check_array(error, vectorAfterRotation, real(inputtype%RootMotion%RotationVel(:, 1), BDKi), testname, tolerance); if (allocated(error)) return
   call check_array(error, vectorAfterRotation, real(inputtype%RootMotion%TranslationAcc(:, 1), BDKi), testname, tolerance); if (allocated(error)) return
   call check_array(error, vectorAfterRotation, real(inputtype%RootMotion%RotationAcc(:, 1), BDKi), testname, tolerance); if (allocated(error)) return

   do i = 1, parametertype%node_total
      call check_array(error, vectorAfterRotation, real(inputtype%PointLoad%Force(1:3, i), BDKi), testname, tolerance); if (allocated(error)) return
      call check_array(error, vectorAfterRotation, real(inputtype%PointLoad%Moment(1:3, i), BDKi), testname, tolerance); if (allocated(error)) return
   end do

   inputtype%DistrLoad%Nnodes = totalnodes
   do i = 1, inputtype%DistrLoad%Nnodes
      call check_array(error, vectorAfterRotation, real(inputtype%DistrLoad%Force(1:3, i), BDKi), testname, tolerance); if (allocated(error)) return
      call check_array(error, vectorAfterRotation, real(inputtype%DistrLoad%Moment(1:3, i), BDKi), testname, tolerance); if (allocated(error)) return
   end do

   call check_array(error, transpose(otherstate%GlbRot), inputtype%RootMotion%Orientation(:, :, 1), testname, tolerance); if (allocated(error)) return

end subroutine

subroutine test_BD_GravityForce(error)
   type(error_type), allocatable, intent(out) :: error
   integer                    :: i, j
   real(BDKi)                 :: gravity(3)
   type(BD_ParameterType)     :: parametertype
   type(BD_MiscVarType)       :: miscvartype
   real(BDKi)                 :: baseline(6)
   integer(IntKi)             :: ErrStat
   character                  :: ErrMsg
   character(1024)            :: testname

   ! --------------------------------------------------------------------------
   testname = "static simple beam under gravity:"
   baseline(1:3) = getGravityInZ()
   baseline(4:6) = [0.0, 0.0, 0.0]

   ! allocate and build the custom types
   parametertype = simpleParameterType(1, 16, 16, 0, 1)
   miscvartype = simpleMiscVarType(parametertype%nqp, parametertype%dof_node, parametertype%elem_total, parametertype%nodes_per_elem)

   gravity = getGravityInZ()

   ! call the subroutine to test
   call BD_GravityForce(1, parametertype, miscvartype, gravity)

   ! test the values
   call check_array(error, baseline, miscvartype%qp%Fg(:, 1, 1), testname, tolerance); if (allocated(error)) return

   call BD_DestroyParam(parametertype, ErrStat, ErrMsg)

end subroutine

subroutine test_BD_QPData_mEta_rho(error)
   type(error_type), allocatable, intent(out) :: error
   integer                    :: i, j
   type(BD_MiscVarType)       :: miscvartype
   type(BD_ParameterType)     :: parametertype
   real(BDKi)                 :: baselineRho(3, 3), baselineRR0mEta(3)
   character(1024)            :: testname
   integer(IntKi)             :: ErrStat
   character                  :: ErrMsg

   ! --------------------------------------------------------------------------
   testname = "static simple beam under gravity:"

   baselineRho(1, :) = [1.0, 0.0, 0.0]
   baselineRho(2, :) = [0.0, 1.0, 0.0]
   baselineRho(3, :) = [0.0, 0.0, 2.0]

   baselineRR0mEta = [0.0, 0.0, 0.0]

   ! allocate and build the custom input types
   parametertype = simpleParameterType(1, 16, 16, 0, 1)
   miscvartype = simpleMiscVarType(parametertype%nqp, parametertype%dof_node, parametertype%elem_total, parametertype%nodes_per_elem)

   ! allocate the results
   call BD_QPData_mEta_rho(parametertype, miscvartype)

   do j = 1, parametertype%elem_total
      do i = 1, parametertype%nqp
         call check_array(error, baselineRho, miscvartype%qp%rho(:, :, i, j), testname, tolerance); if (allocated(error)) return
         call check_array(error, baselineRR0mEta, miscvartype%qp%RR0mEta(:, i, j), testname, tolerance); if (allocated(error)) return
      end do
   end do
   call BD_DestroyParam(parametertype, ErrStat, ErrMsg)
end subroutine

end module
