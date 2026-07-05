module test_NWTC_CheckInput

use testdrive, only: new_unittest, unittest_type, error_type, check
use NWTC_Library, only: IntKi, ErrMsgLen, NewLine, ErrID_None, ErrID_Info, ErrID_Warn, ErrID_Severe, ErrID_Fatal, AbortErrLev
use NWTC_CheckInput   ! NWTC_CheckInput has a PRIVATE default and does NOT re-export NWTC_Library names

implicit none

private
public :: test_NWTC_CheckInput_suite

contains

subroutine test_NWTC_CheckInput_suite(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
               new_unittest("test_collect_splits_multiline_message", test_collect_splits_multiline_message), &
               new_unittest("test_collect_no_error_is_passed",       test_collect_no_error_is_passed), &
               new_unittest("test_collect_warning_does_not_fail",    test_collect_warning_does_not_fail), &
               new_unittest("test_collect_status_override",          test_collect_status_override), &
               new_unittest("test_skipped_status",                   test_skipped_status), &
               new_unittest("test_failed_status_sticky",             test_failed_status_sticky), &
               new_unittest("test_unavailable_status_sticky",        test_unavailable_status_sticky), &
               new_unittest("test_component_status_unknown",         test_component_status_unknown), &
               new_unittest("test_exit_code",                        test_exit_code), &
               new_unittest("test_empty_message_fatal",              test_empty_message_fatal), &
               new_unittest("test_driverrecord_fatal",                test_driverrecord_fatal), &
               new_unittest("test_closereport_lastresort",             test_closereport_lastresort) &
               ]
end subroutine

subroutine test_collect_splits_multiline_message(error)
   type(error_type), allocatable, intent(out) :: error
   type(CheckInputCollectorType) :: collector
   character(ErrMsgLen) :: msg

   ! Exactly the shape NWTC_Base::SetErrStat produces for two severe errors in sequence:
   !   ErrMess = TRIM(ErrMess)//new_line('a')//TRIM(RoutineName)//':'//TRIM(ErrMessLcl)
   msg = 'SD_Init:Number of joints must be at least 2.'//NewLine// &
         'SD_Init:Member 3 references an undefined joint.'

   call CkIn_Collect(collector, 'SubDyn', ErrID_Severe, msg)

   call check(error, collector%NumMsgs, 2); if (allocated(error)) return
   call check(error, collector%NumErrors, 2); if (allocated(error)) return
   call check(error, collector%NumWarnings, 0); if (allocated(error)) return
   call check(error, trim(collector%Msgs(1)%Source), 'SD_Init'); if (allocated(error)) return
   call check(error, trim(collector%Msgs(1)%Text), 'Number of joints must be at least 2.'); if (allocated(error)) return
   call check(error, trim(collector%Msgs(2)%Text), 'Member 3 references an undefined joint.'); if (allocated(error)) return
   call check(error, CkIn_ComponentStatus(collector, 'SubDyn'), CkIn_St_Failed)
end subroutine

subroutine test_collect_no_error_is_passed(error)
   type(error_type), allocatable, intent(out) :: error
   type(CheckInputCollectorType) :: collector
   call CkIn_Collect(collector, 'AeroDyn', ErrID_None, '')
   call check(error, collector%NumMsgs, 0); if (allocated(error)) return
   call check(error, CkIn_ComponentStatus(collector, 'AeroDyn'), CkIn_St_Passed)
end subroutine

subroutine test_collect_warning_does_not_fail(error)
   type(error_type), allocatable, intent(out) :: error
   type(CheckInputCollectorType) :: collector
   call CkIn_Collect(collector, 'ServoDyn', ErrID_Warn, 'ServoDyn_Init:Using default gains.')
   call check(error, collector%NumWarnings, 1); if (allocated(error)) return
   call check(error, collector%NumErrors, 0); if (allocated(error)) return
   call check(error, CkIn_ComponentStatus(collector, 'ServoDyn'), CkIn_St_Passed)
end subroutine

subroutine test_collect_status_override(error)
   type(error_type), allocatable, intent(out) :: error
   type(CheckInputCollectorType) :: collector
   call CkIn_Collect(collector, 'HydroDyn', ErrID_None, '', Status='not_used')
   call check(error, CkIn_ComponentStatus(collector, 'HydroDyn'), CkIn_St_NotUsed)
end subroutine

subroutine test_skipped_status(error)
   type(error_type), allocatable, intent(out) :: error
   type(CheckInputCollectorType) :: collector
   call CkIn_Collect(collector, 'FAST_InitMappings', ErrID_Info, &
        'blocked by upstream module failure(s)', Status='skipped')
   call check(error, CkIn_ComponentStatus(collector, 'FAST_InitMappings'), CkIn_St_Skipped); if (allocated(error)) return
   call check(error, collector%NumErrors, 0)   ! info-level note is not an error
end subroutine

subroutine test_failed_status_sticky(error)
   type(error_type), allocatable, intent(out) :: error
   type(CheckInputCollectorType) :: collector
   ! A fatal collect followed by a clean collect for the SAME component must stay failed:
   call CkIn_Collect(collector, 'ElastoDyn', ErrID_Fatal, 'ED_Init:Blade file not found.')
   call CkIn_Collect(collector, 'ElastoDyn', ErrID_None, '')
   call check(error, CkIn_ComponentStatus(collector, 'ElastoDyn'), CkIn_St_Failed)
end subroutine

subroutine test_unavailable_status_sticky(error)
   type(error_type), allocatable, intent(out) :: error
   type(CheckInputCollectorType) :: collector
   ! An 'unavailable' mark (attempted only against fabricated/tainted upstream data) must not be
   ! silently laundered into Passed by a later benign collect for the same component:
   call CkIn_Collect(collector, 'AeroDyn', ErrID_Info, &
        'ElastoDyn initialization failed; attempted with stubbed ElastoDyn interface data', &
        Status='unavailable')
   call CkIn_Collect(collector, 'AeroDyn', ErrID_None, '')
   call check(error, CkIn_ComponentStatus(collector, 'AeroDyn'), CkIn_St_Unavailable); if (allocated(error)) return
   ! But a real failure still beats the taint marker:
   call CkIn_Collect(collector, 'AeroDyn', ErrID_Fatal, 'AD_Init:msg')
   call check(error, CkIn_ComponentStatus(collector, 'AeroDyn'), CkIn_St_Failed)
end subroutine

subroutine test_component_status_unknown(error)
   type(error_type), allocatable, intent(out) :: error
   type(CheckInputCollectorType) :: collector
   logical :: found
   call check(error, CkIn_ComponentStatus(collector, 'DoesNotExist', found), CkIn_St_Unavailable); if (allocated(error)) return
   call check(error, found, .false.)
end subroutine

subroutine test_exit_code(error)
   type(error_type), allocatable, intent(out) :: error
   type(CheckInputCollectorType) :: collector
   call check(error, CkIn_ExitCode(collector), 0); if (allocated(error)) return
   call CkIn_Collect(collector, 'BeamDyn', ErrID_Fatal, 'BD_Init:Blade input file not found.')
   call check(error, CkIn_ExitCode(collector), 1)
end subroutine

subroutine test_empty_message_fatal(error)
   type(error_type), allocatable, intent(out) :: error
   type(CheckInputCollectorType) :: collector
   ! A fatal collect with an empty message: CkIn_Collect still marks the component failed (that
   ! branch does not depend on ErrMsg), but it does NOT increment NumErrors (only the message-split
   ! loop does that, and it is gated on LEN_TRIM(ErrMsg) > 0). This is exactly the divergence
   ! CkIn_CloseReport's overall_status must not reproduce -- it now derives Overall from
   ! CkIn_ExitCode (component status), not from NumErrors, so this empty-message fatal is not
   ! silently reported as passed.
   call CkIn_Collect(collector, 'X', ErrID_Fatal, '')
   call check(error, collector%NumErrors, 0); if (allocated(error)) return
   call check(error, CkIn_ComponentStatus(collector, 'X'), CkIn_St_Failed); if (allocated(error)) return
   call check(error, CkIn_ExitCode(collector), 1)
end subroutine

subroutine test_driverrecord_fatal(error)
   type(error_type), allocatable, intent(out) :: error
   type(CheckInputCollectorType) :: collector
   ! No CkIn_OpenReport call: CkIn_ReportComponent's "not open" severe must be swallowed by
   ! CkIn_DriverRecord (collector state is still updated) rather than crashing or propagating.
   call CkIn_DriverRecord(collector, 'HydroDyn', ErrID_Fatal, 'HD_Init:WAMIT file not found.')
   call check(error, CkIn_ComponentStatus(collector, 'HydroDyn'), CkIn_St_Failed); if (allocated(error)) return
   call check(error, CkIn_ExitCode(collector), 1)
end subroutine

subroutine test_closereport_lastresort(error)
   type(error_type), allocatable, intent(out) :: error
   type(CheckInputCollectorType) :: collector
   integer(IntKi) :: ErrStat, Un, IOS
   character(ErrMsgLen) :: ErrMsg
   character(*), parameter :: FileName = 'checkinput.verify.yaml'
   character(256) :: Line
   logical :: FileExists, FoundOverallFailed, FoundDriverComp

   ! Never opened: driver failed before it knew its RootName. Clean up any leftover file from a
   ! previous aborted run of this test first.
   inquire(file=FileName, exist=FileExists)
   if (FileExists) then
      open(newunit=Un, file=FileName, status='old')
      close(Un, status='delete')
   end if

   call CkIn_Collect(collector, 'Driver', ErrID_Fatal, 'X:bad arg')
   call CkIn_CloseReport(collector, ErrStat, ErrMsg)

   call check(error, ErrStat < AbortErrLev); if (allocated(error)) return

   inquire(file=FileName, exist=FileExists)
   call check(error, FileExists); if (allocated(error)) return

   FoundOverallFailed = .false.
   FoundDriverComp    = .false.
   open(newunit=Un, file=FileName, status='old', action='read')
   do
      read(Un, '(A)', iostat=IOS) Line
      if (IOS /= 0) exit
      if (index(Line, 'overall_status: failed') > 0) FoundOverallFailed = .true.
      if (index(Line, '- name: Driver') > 0)         FoundDriverComp    = .true.
   end do
   close(Un, status='delete')

   call check(error, FoundOverallFailed); if (allocated(error)) return
   call check(error, FoundDriverComp)
end subroutine

end module
