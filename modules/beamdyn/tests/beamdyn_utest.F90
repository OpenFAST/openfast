program beamdyn_utest
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type

   use test_BD_Crv, only: test_BD_Crv_suite
   use test_BD_diffmtc, only: test_BD_diffmtc_suite
   use test_BD_InitializeNodalLocations, only: test_BD_InitializeNodalLocations_suite
   use test_BD_MemberEta, only: test_BD_MemberEta_suite
   use test_BD_Misc, only: test_BD_Misc_suite
   use test_BD_QuadraturePointData, only: test_BD_QuadraturePointData_suite
   use test_BD_ShapeFuncs, only: test_BD_ShapeFuncs_suite
   use test_BD_TrapezoidalPointWeight, only: test_BD_TrapezoidalPointWeight_suite
   use NWTC_Num
   
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'
   
   stat = 0
   
   call SetConstants()
   
   testsuites = [ &
                new_testsuite("Crv", test_BD_Crv_suite), &
                new_testsuite("diffmtc", test_BD_diffmtc_suite), &
                new_testsuite("InitializeNodalLocations", test_BD_InitializeNodalLocations_suite), &
                new_testsuite("MemberEta", test_BD_MemberEta_suite), &
                new_testsuite("Misc", test_BD_Misc_suite), &
                new_testsuite("QuadraturePointData", test_BD_QuadraturePointData_suite), &
                new_testsuite("ShapeFuncs", test_BD_ShapeFuncs_suite), &
                new_testsuite("TrapezoidalPointWeight", test_BD_TrapezoidalPointWeight_suite) &
                ]
   
   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   
   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if

   write (error_unit, fmt) "All tests PASSED"
   
   end program
   