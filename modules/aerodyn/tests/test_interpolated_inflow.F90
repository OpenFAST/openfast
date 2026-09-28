program test_interpolated_inflow
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: check, error_type, new_testsuite, run_testsuite, testsuite_type, unittest_type, new_unittest
   use AeroDyn_Driver_Subs, only: dat, Dvr_Init, Dvr_InitCase, Dvr_CleanUp
   use AeroDyn, only: AD_UpdateStates
   use AeroDyn_Types, only: AD_InputType, AD_CopyInput, AD_DestroyInput
   use IfW_FlowField_Types, only: FlowFieldType
   use NWTC_Library, only: ReKi, DbKi, IntKi, AbortErrLev, ErrMsgLen, MESH_NEWCOPY
   implicit none

   type(testsuite_type), allocatable :: suites(:)
   character(16) :: scenario
   integer :: status
   real(ReKi), parameter :: translation_speed = 20.0_ReKi

   call get_environment_variable('AD_INFLOW_CASE', scenario)
   if (scenario /= 'time' .and. scenario /= 'time_inputs' .and. &
       scenario /= 'steady' .and. scenario /= 'null') error stop 'Unknown inflow test case'
   suites = [new_testsuite(trim(scenario), collect)]
   status = 0
   call run_testsuite(suites(1)%collect, error_unit, status, parallel=.false.)
   if (status /= 0) error stop 'AeroDyn inflow test failed'

contains

   subroutine collect(tests)
      type(unittest_type), allocatable, intent(out) :: tests(:)
      tests = [new_unittest(trim(scenario), run_case)]
   end subroutine collect

   subroutine run_case(error)
      type(error_type), allocatable, intent(out) :: error
      type(AD_InputType) :: inputs(2)
      type(FlowFieldType), pointer :: flow_field
      real(DbKi) :: input_times(2), sample_times(2)
      real(ReKi) :: expected(3), wind(3), velocity(3), z, max_input_error, max_bemt_error
      real(ReKi) :: fixture_wind_delta, fixture_height_delta
      integer(IntKi) :: err_stat, blade, node, slot, source, component
      character(ErrMsgLen) :: err_msg
      logical :: initialized, remove_flow_field, time_varying

      initialized = .false.
      remove_flow_field = scenario == 'null'
      time_varying = scenario == 'time' .or. scenario == 'time_inputs'
      call Dvr_Init(dat%dvr, dat%ADI, dat%FED, dat%SeaSt, err_stat, err_msg)
      call check(error, err_stat < AbortErrLev, trim(err_msg)); if (allocated(error)) return
      call Dvr_InitCase(1, dat%dvr, dat%ADI, dat%FED, dat%SeaSt, err_stat, err_msg)
      call check(error, err_stat < AbortErrLev, trim(err_msg)); if (allocated(error)) return
      initialized = .true.

      input_times = [real(dat%dvr%dt, DbKi), 0.0_DbKi]
      sample_times = [0.0_DbKi, real(dat%dvr%dt, DbKi)]
      flow_field => dat%ADI%p%AD%FlowField
      call check(error, associated(flow_field), 'Fixture must provide a flow field')
      if (allocated(error)) return
      call check(error, flow_field%Uniform%DataSize >= 2, 'Fixture must provide wind samples')
      if (allocated(error)) return
      do slot = 1, 2
         call AD_CopyInput(dat%ADI%u(1)%AD, inputs(slot), MESH_NEWCOPY, err_stat, err_msg)
         call check(error, err_stat < AbortErrLev, trim(err_msg)); if (allocated(error)) return
         call set_translation(inputs(slot), real(input_times(slot), ReKi))
      end do
      fixture_height_delta = abs(inputs(1)%rotors(1)%BladeMotion(1)%TranslationDisp(3,1) &
                                 - inputs(2)%rotors(1)%BladeMotion(1)%TranslationDisp(3,1))
      fixture_wind_delta = abs(flow_field%Uniform%VelH(2) - flow_field%Uniform%VelH(1))
      call check(error, fixture_height_delta > 0.1_ReKi, 'Fixture must sample distinct heights')
      if (allocated(error)) return
      if (time_varying) then
         call check(error, fixture_wind_delta > 0.5_ReKi, 'Fixture must sample distinct wind times')
         if (allocated(error)) return
         call check(error, all(abs(flow_field%Uniform%ShrV - 1.0_ReKi) < 1.0e-6_ReKi), &
                    'Fixture must use unit power-law shear')
         if (allocated(error)) return
      else if (scenario == 'steady') then
         call check(error, fixture_wind_delta < 1.0e-6_ReKi, 'Steady fixture must have constant wind')
         if (allocated(error)) return
      end if

      if (remove_flow_field) then
         do slot = 1, 2
            do blade = 1, size(dat%ADI%m%AD%Inflow(slot)%RotInflow(1)%Blade)
               dat%ADI%m%AD%Inflow(slot)%RotInflow(1)%Blade(blade)%InflowVel(1,:) = 2.0_ReKi + slot
               dat%ADI%m%AD%Inflow(slot)%RotInflow(1)%Blade(blade)%InflowVel(2,:) = 0.5_ReKi * slot
               dat%ADI%m%AD%Inflow(slot)%RotInflow(1)%Blade(blade)%InflowVel(3,:) = 0.25_ReKi * slot
            end do
            dat%ADI%m%AD%Inflow(slot)%RotInflow(1)%InflowOnHub(:,1) = [2.0_ReKi + slot, 0.5_ReKi * slot, 0.25_ReKi * slot]
         end do
         nullify(dat%ADI%p%AD%FlowField)
      end if

      call AD_UpdateStates(0.0_DbKi, 0_IntKi, inputs, input_times, dat%ADI%p%AD, dat%ADI%x(1)%AD, &
                           dat%ADI%xd(1)%AD, dat%ADI%z(1)%AD, dat%ADI%OtherState(1)%AD, dat%ADI%m%AD, err_stat, err_msg)
      if (remove_flow_field) dat%ADI%p%AD%FlowField => flow_field
      call check(error, err_stat < AbortErrLev, trim(err_msg)); if (allocated(error)) return

      max_input_error = 0.0_ReKi
      max_bemt_error = 0.0_ReKi
      do slot = 1, 2
         source = 3 - slot
         do blade = 1, size(inputs(source)%rotors(1)%BladeMotion)
            do node = 1, inputs(source)%rotors(1)%BladeMotion(blade)%NNodes
               z = inputs(source)%rotors(1)%BladeMotion(blade)%Position(3,node) &
                   + inputs(source)%rotors(1)%BladeMotion(blade)%TranslationDisp(3,node)
               wind = 0.0_ReKi
               if (scenario == 'steady') then
                  wind(1) = 8.0_ReKi
               else if (time_varying) then
                  wind(1) = (8.0_ReKi + 100.0_ReKi * real(sample_times(slot), ReKi)) * z / 140.0_ReKi
               else
                  wind = [3.0_ReKi, 0.5_ReKi, 0.25_ReKi]
               end if
               velocity = inputs(source)%rotors(1)%BladeMotion(blade)%TranslationVel(:,node)
               do component = 1, 3
                  expected(component) = dot_product(wind - velocity, &
                       dat%ADI%m%AD%rotors(1)%orientationAnnulus(component,:,node,blade))
               end do
               max_bemt_error = max(max_bemt_error, &
                    abs(dat%ADI%m%AD%rotors(1)%BEMT_u(slot)%Vx(node,blade) - expected(1)), &
                    abs(dat%ADI%m%AD%rotors(1)%BEMT_u(slot)%Vy(node,blade) - expected(2)), &
                    abs(dat%ADI%m%AD%rotors(1)%BEMT_u(slot)%Vz(node,blade) - expected(3)))
            end do
         end do
      end do

      do slot = 1, 2
         do blade = 1, size(inputs(slot)%rotors(1)%BladeMotion)
            do node = 1, inputs(slot)%rotors(1)%BladeMotion(blade)%NNodes
               z = inputs(slot)%rotors(1)%BladeMotion(blade)%Position(3,node) &
                   + inputs(slot)%rotors(1)%BladeMotion(blade)%TranslationDisp(3,node)
               if (scenario == 'steady') then
                  expected = [8.0_ReKi, 0.0_ReKi, 0.0_ReKi]
               else if (time_varying) then
                  expected = [((8.0_ReKi + 100.0_ReKi * real(input_times(slot), ReKi)) * z / 140.0_ReKi), &
                              0.0_ReKi, 0.0_ReKi]
               else
                  expected = [2.0_ReKi + slot, 0.5_ReKi * slot, 0.25_ReKi * slot]
               end if
               max_input_error = max(max_input_error, &
                    maxval(abs(dat%ADI%m%AD%Inflow(slot)%RotInflow(1)%Blade(blade)%InflowVel(:,node) - expected)))
            end do
         end do
      end do

      do slot = 1, 2
         call AD_DestroyInput(inputs(slot), err_stat, err_msg)
      end do
      if (initialized) call Dvr_CleanUp(dat%dvr, dat%ADI, dat%FED, initialized, err_stat, err_msg)

      write(error_unit, '(a,1x,a,2(a,es12.4))') 'CASE', trim(scenario), &
           ' max_input_error=', max_input_error, ' max_bemt_error=', max_bemt_error
      if (scenario == 'time_inputs') then
         call check(error, max_input_error < 2.0e-3_ReKi, 'Original input inflow differs from analytic wind')
         return
      end if
      call check(error, max_bemt_error < 2.0e-3_ReKi, 'BEMT slot differs from analytic relative wind')
      if (allocated(error)) return
      call check(error, max_input_error < 2.0e-3_ReKi, 'Original input inflow differs from analytic wind')
   end subroutine run_case

   subroutine set_translation(input, time)
      type(AD_InputType), intent(inout) :: input
      real(ReKi), intent(in) :: time
      integer :: blade
      input%rotors(1)%HubMotion%TranslationDisp(3,:) = translation_speed * time
      input%rotors(1)%HubMotion%TranslationVel(3,:) = translation_speed
      input%rotors(1)%NacelleMotion%TranslationDisp(3,:) = translation_speed * time
      input%rotors(1)%NacelleMotion%TranslationVel(3,:) = translation_speed
      do blade = 1, size(input%rotors(1)%BladeMotion)
         input%rotors(1)%BladeRootMotion(blade)%TranslationDisp(3,:) = translation_speed * time
         input%rotors(1)%BladeRootMotion(blade)%TranslationVel(3,:) = translation_speed
         input%rotors(1)%BladeMotion(blade)%TranslationDisp(3,:) = translation_speed * time
         input%rotors(1)%BladeMotion(blade)%TranslationVel(3,:) = translation_speed
      end do
   end subroutine set_translation

end program test_interpolated_inflow
