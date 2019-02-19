!  DBEMT_Test.f90 
!
!  FUNCTIONS:
!  DBEMT_Test - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: DBEMT_Test
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

program DBEMT_Test
   use DBEMT
   use NWTC_Library
   
   implicit none

   ! Variables
   integer(IntKi)                  :: errStat       !< Error status of the operation
   character(ErrMsgLen)            :: errMsg        !< Error message if ErrStat /= ErrID_None
   logical                         :: result
      ! Initialize the NWTC Subroutine Library

   call NWTC_Init( EchoLibVer=.FALSE. )
   
   ! Body of DBEMT_Test
   print *, 'Beginning DBEMT Testing'
   errStat = ErrID_None
   errMsg  = ""
   
   result =  DBEMT_TEST01(errStat, errMsg)
   print *, 'Test01: ', result
#ifdef PRINT_ERRMSG
   print *, trim(errMsg)
   print *
#endif

   result =  DBEMT_TEST02(errStat, errMsg)
   print *, 'Test02: ', result
#ifdef PRINT_ERRMSG
   print *, trim(errMsg)
   print *
#endif

result =  DBEMT_TEST03(errStat, errMsg)
   print *, 'Test03: ', result
#ifdef PRINT_ERRMSG
   print *, trim(errMsg)
   print *
#endif

   
   result =  DBEMT_TEST04(errStat, errMsg)
   print *, 'Test04: ', result
#ifdef PRINT_ERRMSG
   print *, trim(errMsg)
   print *
#endif

   
   result =  DBEMT_TEST05(errStat, errMsg)
   print *, 'Test05: ', result
#ifdef PRINT_ERRMSG
   print *, trim(errMsg)
   print *
#endif

   
   result =  DBEMT_TEST06(errStat, errMsg)
   print *, 'Test06: ', result
#ifdef PRINT_ERRMSG
   print *, trim(errMsg)
   print *
#endif

   
   result =  DBEMT_TEST07(errStat, errMsg)
   print *, 'Test07: ', result
#ifdef PRINT_ERRMSG
   print *, trim(errMsg)
   print *
#endif
   
   result =  DBEMT_TEST08(errStat, errMsg)
   print *, 'Test08: ', result
#ifdef PRINT_ERRMSG
   print *, trim(errMsg)
   print *
#endif

result =  DBEMT_TEST09(errStat, errMsg)
   print *, 'Test09: ', result
#ifdef PRINT_ERRMSG
   print *, trim(errMsg)
   print *
#endif

result =  DBEMT_TEST10(errStat, errMsg)
   print *, 'Test10: ', result
#ifdef PRINT_ERRMSG
   print *, trim(errMsg)
   print *
#endif

result =  DBEMT_TEST11(errStat, errMsg)
   print *, 'Test11: ', result
#ifdef PRINT_ERRMSG
   print *, trim(errMsg)
   print *
#endif

contains

logical function DBEMT_TEST01(errStat, errMsg)

   
   integer(IntKi), INTENT(out)                  :: errStat       !< Error status of the operation
   character(ErrMsgLen), INTENT(out)                    :: errMsg        !< Error message if ErrStat /= ErrID_None
   type(DBEMT_InitInputType)       :: InitInp       !< Input data for initialization routine
   type(DBEMT_InputType)           :: u             !< An initial guess for the input; input mesh must be defined
   type(DBEMT_ParameterType)       :: p             !< Parameters
   type(DBEMT_ContinuousStateType) :: x             !< Initial continuous states
   type(DBEMT_MiscVarType)         :: m             !< Initial misc/optimization variables
   real(DbKi)                      :: interval      !< Coupling interval in seconds: the rate that
                                                     !!   (1) DBEMT_UpdateStates() is called in loose coupling &
                                                     !!   (2) DBEMT_UpdateDiscState() is called in tight coupling.
                                                     !!   Input is the suggested time from the glue code;
                                                     !!   Output is the actual coupling interval that will be used
                                                     !!   by the glue code.
   type(DBEMT_InitOutputType)      :: InitOut       !< Output for initialization routine
   
   
   
   DBEMT_TEST01 = .false.
   
   ! This test will evaluate the initialization of the module 
   
   errStat = ErrID_None
   errMsg  = ""

   InitInp%DBEMT_Mod = -1
   
   call DBEMT_Init(InitInp, u, p, x, m, Interval, InitOut, errStat, errMsg)

   if (errStat == AbortErrLev) DBEMT_TEST01 = .true.

end function DBEMT_TEST01

logical function DBEMT_TEST02(errStat, errMsg)

   
   integer(IntKi), INTENT(out)                  :: errStat       !< Error status of the operation
   character(ErrMsgLen), INTENT(out)                    :: errMsg        !< Error message if ErrStat /= ErrID_None
   type(DBEMT_InitInputType)       :: InitInp       !< Input data for initialization routine
   type(DBEMT_InputType)           :: u(2)             !< An initial guess for the input; input mesh must be defined
   type(DBEMT_ParameterType)       :: p             !< Parameters
   type(DBEMT_ContinuousStateType) :: x             !< Initial continuous states
   type(DBEMT_MiscVarType)         :: m             !< Initial misc/optimization variables
   real(DbKi)                      :: interval      !< Coupling interval in seconds: the rate that
                                                     !!   (1) DBEMT_UpdateStates() is called in loose coupling &
                                                     !!   (2) DBEMT_UpdateDiscState() is called in tight coupling.
                                                     !!   Input is the suggested time from the glue code;
                                                     !!   Output is the actual coupling interval that will be used
                                                     !!   by the glue code.
   type(DBEMT_InitOutputType)      :: InitOut       !< Output for initialization routine
   
   integer(IntKi)                              :: errStat2      ! temporary error status of the operation
   character(ErrMsgLen)                        :: errMsg2       ! temporary error message 
   character(*), parameter                     :: RoutineName = 'DBEMT_TEST02'
   
     ! Initialize variables for this routine

   errStat2 = ErrID_None
   errMsg2  = ""

   
   DBEMT_TEST02 = .false.
   
   ! This test will evaluate the initialization of the module 
   
   errStat = ErrID_None
   errMsg  = ""

   InitInp%DBEMT_Mod = 0
   call DBEMT_Init(InitInp, u(1), p, x, m, Interval, InitOut, errStat2, errMsg2)
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName)
      
   call DBEMT_UpdateStates(1, 1, interval, u,  p, x, m, errStat2, errMsg2 )  
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName)
      
      ! Since DBEMT is turned off, their should be no errors, even though other initialization data is missing
   if (errStat == ErrID_None) DBEMT_TEST02 = .true.

end function DBEMT_TEST02

logical function DBEMT_TEST03(errStat, errMsg)

   
   integer(IntKi), INTENT(out)                  :: errStat       !< Error status of the operation
   character(ErrMsgLen), INTENT(out)                    :: errMsg        !< Error message if ErrStat /= ErrID_None
   type(DBEMT_InitInputType)       :: InitInp       !< Input data for initialization routine
   type(DBEMT_InputType)           :: u             !< An initial guess for the input; input mesh must be defined
   type(DBEMT_ParameterType)       :: p             !< Parameters
   type(DBEMT_ContinuousStateType) :: x             !< Initial continuous states
   type(DBEMT_MiscVarType)         :: m             !< Initial misc/optimization variables
   real(DbKi)                      :: interval      !< Coupling interval in seconds: the rate that
                                                     !!   (1) DBEMT_UpdateStates() is called in loose coupling &
                                                     !!   (2) DBEMT_UpdateDiscState() is called in tight coupling.
                                                     !!   Input is the suggested time from the glue code;
                                                     !!   Output is the actual coupling interval that will be used
                                                     !!   by the glue code.
   type(DBEMT_InitOutputType)      :: InitOut       !< Output for initialization routine
   
   
   
   DBEMT_TEST03 = .false.
   
   ! This test will evaluate the initialization of the module 
   
   errStat = ErrID_None
   errMsg  = ""
   
   InitInp%DBEMT_Mod = 1
   
      ! The InitInp%numBlades has not been set.  This should produce a fatal error.
      ! The InitInp%numNodes has not been set.  This should produce a fatal error.
      ! The InitInp%rlocal array has not been allocated.  This should produce a fatal error.

   call DBEMT_Init(InitInp, u, p, x, m, Interval, InitOut, errStat, errMsg)

   if (errStat == AbortErrLev) DBEMT_TEST03 = .true.

end function DBEMT_TEST03

logical function DBEMT_TEST04(errStat, errMsg)

   
   integer(IntKi), INTENT(out)                  :: errStat       !< Error status of the operation
   character(ErrMsgLen), INTENT(out)                    :: errMsg        !< Error message if ErrStat /= ErrID_None
   type(DBEMT_InitInputType)       :: InitInp       !< Input data for initialization routine
   type(DBEMT_InputType)           :: u             !< An initial guess for the input; input mesh must be defined
   type(DBEMT_ParameterType)       :: p             !< Parameters
   type(DBEMT_ContinuousStateType) :: x             !< Initial continuous states
   type(DBEMT_MiscVarType)         :: m             !< Initial misc/optimization variables
   real(DbKi)                      :: interval      !< Coupling interval in seconds: the rate that
                                                     !!   (1) DBEMT_UpdateStates() is called in loose coupling &
                                                     !!   (2) DBEMT_UpdateDiscState() is called in tight coupling.
                                                     !!   Input is the suggested time from the glue code;
                                                     !!   Output is the actual coupling interval that will be used
                                                     !!   by the glue code.
   type(DBEMT_InitOutputType)      :: InitOut       !< Output for initialization routine
   
   
   
   DBEMT_TEST04 = .false.
   
   ! This test will evaluate the initialization of the module 
   
   errStat = ErrID_None
   errMsg  = ""

   InitInp%DBEMT_Mod = 1
   InitInp%numBlades = 1
   call DBEMT_Init(InitInp, u, p, x, m, Interval, InitOut, errStat, errMsg)

      ! The InitInp%numNodes has not been set.  This should produce a fatal error.
      ! The InitInp%rlocal array has not been allocated.  This should produce a fatal error.
   if (errStat == AbortErrLev) DBEMT_TEST04 = .true.

end function DBEMT_TEST04

logical function DBEMT_TEST05(errStat, errMsg)

   
   integer(IntKi), INTENT(out)                  :: errStat       !< Error status of the operation
   character(ErrMsgLen), INTENT(out)                    :: errMsg        !< Error message if ErrStat /= ErrID_None
   type(DBEMT_InitInputType)       :: InitInp       !< Input data for initialization routine
   type(DBEMT_InputType)           :: u             !< An initial guess for the input; input mesh must be defined
   type(DBEMT_ParameterType)       :: p             !< Parameters
   type(DBEMT_ContinuousStateType) :: x             !< Initial continuous states
   type(DBEMT_MiscVarType)         :: m             !< Initial misc/optimization variables
   real(DbKi)                      :: interval      !< Coupling interval in seconds: the rate that
                                                     !!   (1) DBEMT_UpdateStates() is called in loose coupling &
                                                     !!   (2) DBEMT_UpdateDiscState() is called in tight coupling.
                                                     !!   Input is the suggested time from the glue code;
                                                     !!   Output is the actual coupling interval that will be used
                                                     !!   by the glue code.
   type(DBEMT_InitOutputType)      :: InitOut       !< Output for initialization routine
   
   
   
   DBEMT_TEST05 = .false.
   
   ! This test will evaluate the initialization of the module 
   
   errStat = ErrID_None
   errMsg  = ""

   InitInp%DBEMT_Mod = 1
   InitInp%numBlades = 1
   InitInp%numNodes  = 3
   call DBEMT_Init(InitInp, u, p, x, m, Interval, InitOut, errStat, errMsg)

      ! The InitInp%rlocal array has not been allocated.  This should produce a fatal error.
   if (errStat == AbortErrLev) DBEMT_TEST05 = .true.

end function DBEMT_TEST05

logical function DBEMT_TEST06(errStat, errMsg)

   
   integer(IntKi), INTENT(out)                  :: errStat       !< Error status of the operation
   character(ErrMsgLen), INTENT(out)                    :: errMsg        !< Error message if ErrStat /= ErrID_None
   type(DBEMT_InitInputType)       :: InitInp       !< Input data for initialization routine
   type(DBEMT_InputType)           :: u             !< An initial guess for the input; input mesh must be defined
   type(DBEMT_ParameterType)       :: p             !< Parameters
   type(DBEMT_ContinuousStateType) :: x             !< Initial continuous states
   type(DBEMT_MiscVarType)         :: m             !< Initial misc/optimization variables
   real(DbKi)                      :: interval      !< Coupling interval in seconds: the rate that
                                                     !!   (1) DBEMT_UpdateStates() is called in loose coupling &
                                                     !!   (2) DBEMT_UpdateDiscState() is called in tight coupling.
                                                     !!   Input is the suggested time from the glue code;
                                                     !!   Output is the actual coupling interval that will be used
                                                     !!   by the glue code.
   type(DBEMT_InitOutputType)      :: InitOut       !< Output for initialization routine
   
   
   
   DBEMT_TEST06 = .false.
   
   ! This test will evaluate the initialization of the module 
   
   errStat = ErrID_None
   errMsg  = ""

   InitInp%DBEMT_Mod = 1
   InitInp%numBlades = 1
   InitInp%numNodes  = 2
   allocate(InitInp%rlocal(InitInp%numNodes,InitInp%numBlades))
   
   ! The rlocal values have been allocated, but not initialized.  This should produce a fatal error, but the values 
   ! that are set are indeterminate, so we will never know if an error is triggered or not.
   
   call DBEMT_Init(InitInp, u, p, x, m, Interval, InitOut, errStat, errMsg)

   if (errStat == AbortErrLev) DBEMT_TEST06 = .true.

end function DBEMT_TEST06

logical function DBEMT_TEST07(errStat, errMsg)

   
   integer(IntKi), INTENT(out)                  :: errStat       !< Error status of the operation
   character(ErrMsgLen), INTENT(out)                    :: errMsg        !< Error message if ErrStat /= ErrID_None
   type(DBEMT_InitInputType)       :: InitInp       !< Input data for initialization routine
   type(DBEMT_InputType)           :: u             !< An initial guess for the input; input mesh must be defined
   type(DBEMT_ParameterType)       :: p             !< Parameters
   type(DBEMT_ContinuousStateType) :: x             !< Initial continuous states
   type(DBEMT_MiscVarType)         :: m             !< Initial misc/optimization variables
   real(DbKi)                      :: interval      !< Coupling interval in seconds: the rate that
                                                     !!   (1) DBEMT_UpdateStates() is called in loose coupling &
                                                     !!   (2) DBEMT_UpdateDiscState() is called in tight coupling.
                                                     !!   Input is the suggested time from the glue code;
                                                     !!   Output is the actual coupling interval that will be used
                                                     !!   by the glue code.
   type(DBEMT_InitOutputType)      :: InitOut       !< Output for initialization routine
   
   integer(IntKi)                  :: i,j
   
   DBEMT_TEST07 = .false.
   
   ! This test will evaluate the initialization of the module 
   
   errStat = ErrID_None
   errMsg  = ""

   InitInp%DBEMT_Mod = 1
   InitInp%numBlades = 1
   InitInp%numNodes  = 2
   allocate(InitInp%rlocal(InitInp%numNodes,InitInp%numBlades))
   do j=1,InitInp%numBlades
      do i=1,InitInp%numNodes
         InitInp%rlocal(i,j) = real(i,ReKi)
      end do
   end do
   
   call DBEMT_Init(InitInp, u, p, x, m, Interval, InitOut, errStat, errMsg)

   if (errStat == AbortErrLev) DBEMT_TEST07 = .true.

end function DBEMT_TEST07

logical function DBEMT_TEST08(errStat, errMsg)

   
   integer(IntKi), INTENT(out)                  :: errStat       !< Error status of the operation
   character(ErrMsgLen), INTENT(out)                    :: errMsg        !< Error message if ErrStat /= ErrID_None
   type(DBEMT_InitInputType)       :: InitInp       !< Input data for initialization routine
   type(DBEMT_InputType)           :: u(2)             !< An initial guess for the input; input mesh must be defined
   type(DBEMT_ParameterType)       :: p             !< Parameters
   type(DBEMT_ContinuousStateType) :: x             !< Initial continuous states
   type(DBEMT_MiscVarType)         :: m             !< Initial misc/optimization variables
   real(DbKi)                      :: interval      !< Coupling interval in seconds: the rate that
                                                     !!   (1) DBEMT_UpdateStates() is called in loose coupling &
                                                     !!   (2) DBEMT_UpdateDiscState() is called in tight coupling.
                                                     !!   Input is the suggested time from the glue code;
                                                     !!   Output is the actual coupling interval that will be used
                                                     !!   by the glue code.
   type(DBEMT_InitOutputType)      :: InitOut       !< Output for initialization routine
   
   integer(IntKi)                  :: i,j
   integer(IntKi)                              :: errStat2      ! temporary error status of the operation
   character(ErrMsgLen)                        :: errMsg2       ! temporary error message 
   character(*), parameter                     :: RoutineName = 'DBEMT_TEST08'
   
     ! Initialize variables for this routine

   errStat2 = ErrID_None
   errMsg2  = ""
   DBEMT_TEST08 = .false.
   
   ! This test will evaluate the initialization of the module 
   
   errStat = ErrID_None
   errMsg  = ""

   InitInp%DBEMT_Mod = 1
   InitInp%numBlades = 1
   InitInp%numNodes  = 2
   allocate(InitInp%rlocal(InitInp%numNodes,InitInp%numBlades))
   do j=1,InitInp%numBlades
      do i=1,InitInp%numNodes
         InitInp%rlocal(i,j) = real(i,ReKi)
      end do
   end do
   Interval = 0.1_ReKi
   InitInp%tau1_const = 0.0_ReKi
   
   call DBEMT_Init(InitInp, u(1), p, x, m, Interval, InitOut, errStat2, errMsg2)
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName)
     
   u(1)%vind_s(1) = 10.0_ReKi
   u(1)%vind_s(2) =  0.0_ReKi
   
   u(2)%vind_s(1) = 12.0_ReKi
   u(2)%vind_s(2) =  0.0_ReKi
   
   call DBEMT_UpdateStates(1, 1, 0.0_DbKi, u,  p, x, m, errStat2, errMsg2 )  
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName)
   
   if ( (errStat == ErrID_None) .and. (EqualRealNos(u(1)%vind_s(1),x%vind(1,1,1)) ) &
         .and. (EqualRealNos(u(1)%vind_s(1),x%vind_1(1,1,1)) )                       &
         .and. (EqualRealNos(u(1)%vind_s(2),x%vind(2,1,1)) )                         &
         .and. (EqualRealNos(u(1)%vind_s(2),x%vind_1(2,1,1)) )) then
      DBEMT_TEST08 = .true.
   else
      DBEMT_TEST08 = .false.
   end if
   
   call DBEMT_End(u, p, x, m, ErrStat, ErrMsg )
   
   InitInp%tau1_const = 1.0_ReKi
   call DBEMT_Init(InitInp, u(1), p, x, m, Interval, InitOut, errStat2, errMsg2)
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName)
      
   u(1)%vind_s(1) = 10.0_ReKi
   u(1)%vind_s(2) =  0.0_ReKi
   
   u(2)%vind_s(1) = 12.0_ReKi
   u(2)%vind_s(2) =  0.0_ReKi
   
   call DBEMT_UpdateStates(1, 1, 0.0_DbKi, u,  p, x, m, errStat2, errMsg2 )  
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName)
   
   if ( (errStat == ErrID_None) .and. (EqualRealNos(u(1)%vind_s(1),x%vind(1,1,1)) ) &
         .and. (EqualRealNos(u(1)%vind_s(1),x%vind_1(1,1,1)) )                       &
         .and. (EqualRealNos(u(1)%vind_s(2),x%vind(2,1,1)) )                         &
         .and. (EqualRealNos(u(1)%vind_s(2),x%vind_1(2,1,1)) )) then
      DBEMT_TEST08 = (.true. .and. DBEMT_TEST08)
   else
      DBEMT_TEST08 = .false.
      
   end if
   
   call DBEMT_End(u, p, x, m, ErrStat, ErrMsg )
   
end function DBEMT_TEST08

logical function DBEMT_TEST09(errStat, errMsg)

   
   integer(IntKi), INTENT(out)                  :: errStat       !< Error status of the operation
   character(ErrMsgLen), INTENT(out)                    :: errMsg        !< Error message if ErrStat /= ErrID_None
   type(DBEMT_InitInputType)       :: InitInp       !< Input data for initialization routine
   type(DBEMT_InputType)           :: u(2)             !< An initial guess for the input; input mesh must be defined
   type(DBEMT_ParameterType)       :: p             !< Parameters
   type(DBEMT_ContinuousStateType) :: x             !< Initial continuous states
   type(DBEMT_MiscVarType)         :: m             !< Initial misc/optimization variables
   real(DbKi)                      :: interval      !< Coupling interval in seconds: the rate that
                                                     !!   (1) DBEMT_UpdateStates() is called in loose coupling &
                                                     !!   (2) DBEMT_UpdateDiscState() is called in tight coupling.
                                                     !!   Input is the suggested time from the glue code;
                                                     !!   Output is the actual coupling interval that will be used
                                                     !!   by the glue code.
   type(DBEMT_InitOutputType)      :: InitOut       !< Output for initialization routine
   
   integer(IntKi)                  :: i,j,n
   real(DbKi)                      :: t
   integer(IntKi)                              :: errStat2      ! temporary error status of the operation
   character(ErrMsgLen)                        :: errMsg2       ! temporary error message 
   character(*), parameter                     :: RoutineName = 'DBEMT_TEST09'
   
     ! Initialize variables for this routine

   errStat2 = ErrID_None
   errMsg2  = ""
   DBEMT_TEST09 = .true.
   
   ! This test will evaluate the initialization of the module 
   
   errStat = ErrID_None
   errMsg  = ""

   InitInp%DBEMT_Mod = 1
   InitInp%numBlades = 1
   InitInp%numNodes  = 2
   allocate(InitInp%rlocal(InitInp%numNodes,InitInp%numBlades))
   do j=1,InitInp%numBlades
      do i=1,InitInp%numNodes
         InitInp%rlocal(i,j) = real(i,ReKi)
      end do
   end do
   Interval = 0.1_ReKi
   InitInp%tau1_const = 0.0_ReKi
   
   call DBEMT_Init(InitInp, u(1), p, x, m, Interval, InitOut, errStat2, errMsg2)
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName)
     
   u(1)%vind_s(1) = 10.0_ReKi
   u(1)%vind_s(2) =  0.0_ReKi
   
   u(2)%vind_s(1) = 10.0_ReKi
   u(2)%vind_s(2) =  0.0_ReKi
   do n=1,10
      t = n*Interval
      call DBEMT_UpdateStates(1, 1, t, u,  p, x, m, errStat2, errMsg2 )  
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName)
   
      if ( (errStat == ErrID_None) .and. (EqualRealNos(u(1)%vind_s(1),x%vind(1,1,1)) ) &
            .and. (EqualRealNos(u(1)%vind_s(1),x%vind_1(1,1,1)) )                       &
            .and. (EqualRealNos(u(1)%vind_s(2),x%vind(2,1,1)) )                         &
            .and. (EqualRealNos(u(1)%vind_s(2),x%vind_1(2,1,1)) )) then
         DBEMT_TEST09 = .true. .and. DBEMT_TEST09
      else
         DBEMT_TEST09 = .false.
      end if
   end do
   call DBEMT_End(u, p, x, m, ErrStat, ErrMsg )
   
   
   
end function DBEMT_TEST09


logical function DBEMT_TEST10(errStat, errMsg)

   
   integer(IntKi), INTENT(out)                  :: errStat       !< Error status of the operation
   character(ErrMsgLen), INTENT(out)                    :: errMsg        !< Error message if ErrStat /= ErrID_None
   type(DBEMT_InitInputType)       :: InitInp       !< Input data for initialization routine
   type(DBEMT_InputType)           :: u(2)             !< An initial guess for the input; input mesh must be defined
   type(DBEMT_ParameterType)       :: p             !< Parameters
   type(DBEMT_ContinuousStateType) :: x             !< Initial continuous states
   type(DBEMT_MiscVarType)         :: m             !< Initial misc/optimization variables
   real(DbKi)                      :: interval      !< Coupling interval in seconds: the rate that
                                                     !!   (1) DBEMT_UpdateStates() is called in loose coupling &
                                                     !!   (2) DBEMT_UpdateDiscState() is called in tight coupling.
                                                     !!   Input is the suggested time from the glue code;
                                                     !!   Output is the actual coupling interval that will be used
                                                     !!   by the glue code.
   type(DBEMT_InitOutputType)      :: InitOut       !< Output for initialization routine
   
   integer(IntKi)                  :: i,j,n
   real(DbKi)                      :: t
   integer(IntKi)                              :: errStat2      ! temporary error status of the operation
   character(ErrMsgLen)                        :: errMsg2       ! temporary error message 
   character(*), parameter                     :: RoutineName = 'DBEMT_TEST10'
   
     ! Initialize variables for this routine

   errStat2 = ErrID_None
   errMsg2  = ""
   DBEMT_TEST10 = .true.
   
   ! This test will evaluate the initialization of the module 
   
   errStat = ErrID_None
   errMsg  = ""

   InitInp%DBEMT_Mod = 1
   InitInp%numBlades = 1
   InitInp%numNodes  = 2
   allocate(InitInp%rlocal(InitInp%numNodes,InitInp%numBlades))
   do j=1,InitInp%numBlades
      do i=1,InitInp%numNodes
         InitInp%rlocal(i,j) = real(i,ReKi)
      end do
   end do
   Interval = 0.1_ReKi
   InitInp%tau1_const = 0.0_ReKi
   
   call DBEMT_Init(InitInp, u(1), p, x, m, Interval, InitOut, errStat2, errMsg2)
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName)
     
   u(1)%vind_s(1) = 10.0_ReKi
   u(1)%vind_s(2) = 0.0_ReKi
   
   u(2)%vind_s(1) = 12.0_ReKi
   u(2)%vind_s(2) =  1.0_ReKi
   do n=0,10
      t = n*Interval
      call DBEMT_UpdateStates(1, 1, t, u,  p, x, m, errStat2, errMsg2 )  
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName)
   
      if ( (errStat == ErrID_None) .and. (EqualRealNos(u(1)%vind_s(1),x%vind(1,1,1)) ) &
            .and. (EqualRealNos(u(1)%vind_s(1),x%vind_1(1,1,1)) )                       &
            .and. (EqualRealNos(u(1)%vind_s(2),x%vind(2,1,1)) )                         &
            .and. (EqualRealNos(u(1)%vind_s(2),x%vind_1(2,1,1)) )) then
         DBEMT_TEST10 = .true. .and. DBEMT_TEST10
      else
         DBEMT_TEST10 = .false.
      end if
      
      u(1)%vind_s(1) = 12.0_ReKi
      u(1)%vind_s(2) = 1.0_ReKi
   
   end do
   call DBEMT_End(u, p, x, m, ErrStat, ErrMsg )
   
   
   
end function DBEMT_TEST10

logical function DBEMT_TEST11(errStat, errMsg)

   
   integer(IntKi), INTENT(out)                  :: errStat       !< Error status of the operation
   character(ErrMsgLen), INTENT(out)                    :: errMsg        !< Error message if ErrStat /= ErrID_None
   type(DBEMT_InitInputType)       :: InitInp       !< Input data for initialization routine
   type(DBEMT_InputType)           :: u(2)             !< An initial guess for the input; input mesh must be defined
   type(DBEMT_ParameterType)       :: p             !< Parameters
   type(DBEMT_ContinuousStateType) :: x             !< Initial continuous states
   type(DBEMT_MiscVarType)         :: m             !< Initial misc/optimization variables
   real(DbKi)                      :: interval      !< Coupling interval in seconds: the rate that
                                                     !!   (1) DBEMT_UpdateStates() is called in loose coupling &
                                                     !!   (2) DBEMT_UpdateDiscState() is called in tight coupling.
                                                     !!   Input is the suggested time from the glue code;
                                                     !!   Output is the actual coupling interval that will be used
                                                     !!   by the glue code.
   type(DBEMT_InitOutputType)      :: InitOut       !< Output for initialization routine
   
   integer(IntKi)                  :: i,j,n
   real(DbKi)                      :: t
   integer(IntKi)                              :: errStat2      ! temporary error status of the operation
   character(ErrMsgLen)                        :: errMsg2       ! temporary error message 
   character(*), parameter                     :: RoutineName = 'DBEMT_TEST11'
   real(ReKi)          :: maxVind1
     ! Initialize variables for this routine

   errStat2 = ErrID_None
   errMsg2  = ""
   DBEMT_TEST11 = .true.
   maxVind1 = 0.0
   ! This test will evaluate the initialization of the module 
   
   errStat = ErrID_None
   errMsg  = ""

   InitInp%DBEMT_Mod = 1
   InitInp%numBlades = 1
   InitInp%numNodes  = 2
   allocate(InitInp%rlocal(InitInp%numNodes,InitInp%numBlades))
   do j=1,InitInp%numBlades
      do i=1,InitInp%numNodes
         InitInp%rlocal(i,j) = real(i,ReKi)
      end do
   end do
   Interval = 0.1_ReKi
   InitInp%tau1_const = .01000_ReKi
   
   call DBEMT_Init(InitInp, u(1), p, x, m, Interval, InitOut, errStat2, errMsg2)
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName)
     
   u(1)%vind_s(1) = 10.0_ReKi
   u(1)%vind_s(2) = 0.0_ReKi
   
   u(2)%vind_s(1) = 12.0_ReKi
   u(2)%vind_s(2) =  1.0_ReKi
   do n=0,10000
      t = n*Interval
      call DBEMT_UpdateStates(1, 1, t, u,  p, x, m, errStat2, errMsg2 )  
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName)
   
      if ( x%vind(1,1,1) > maxVind1 ) maxVind1 = x%vind(1,1,1)
      
      u(1)%vind_s(1) = 12.0_ReKi
      u(1)%vind_s(2) = 1.0_ReKi
   
   end do
   if ( (errStat == ErrID_None) .and. (EqualRealNos(u(1)%vind_s(1),x%vind(1,1,1)) ) &
            .and. (EqualRealNos(u(1)%vind_s(1),x%vind_1(1,1,1)) )                       &
            .and. (EqualRealNos(u(1)%vind_s(2),x%vind(2,1,1)) )                         &
            .and. (EqualRealNos(u(1)%vind_s(2),x%vind_1(2,1,1)) )) then
         DBEMT_TEST11 = (.true. .and. DBEMT_TEST11)
      else
         DBEMT_TEST11 = .false.
      end if
   call DBEMT_End(u, p, x, m, ErrStat, ErrMsg )
   
   
   
end function DBEMT_TEST11


   end program DBEMT_Test

