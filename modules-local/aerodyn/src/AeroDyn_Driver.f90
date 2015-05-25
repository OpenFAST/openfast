program AeroDyn_Driver

   use AeroDyn_Driver_Subs   
    
   implicit none   
   
      ! Program variables

   real(DbKi)                                     :: time                 !< Variable for storing time, in seconds 
   real(DbKi)                                     :: dT_Dvr               !< copy of DT, to make sure AD didn't change it
                                                    
   type(Dvr_SimData)                              :: DvrData              ! The data required for running the AD driver
   type(AeroDyn_Data)                             :: AD                   ! AeroDyn data 
                                                  
   integer(IntKi)                                 :: n                    ! loop counter (for time step)
   integer(IntKi)                                 :: numSteps             ! number of time steps in the simulation
   integer(IntKi)                                 :: errStat              ! Status of error message
   character(ErrMsgLen)                           :: errMsg               ! Error message if ErrStat /= ErrID_None

   !integer                                        :: StrtTime (8)                            ! Start time of simulation (including intialization)
   !integer                                        :: SimStrtTime (8)                         ! Start time of simulation (after initialization)
   !real(ReKi)                                     :: PrevClockTime                           ! Clock time at start of simulation in seconds
   !real                                           :: UsrTime1                                ! User CPU time for simulation initialization
   !real                                           :: UsrTime2                                ! User CPU time for simulation (without intialization)
   !real                                           :: UsrTimeDiff                             ! Difference in CPU time from start to finish of program execution
   !real(DbKi)                                     :: TiLstPrn                                ! The simulation time of the last print
   !real(DbKi)                                     :: SttsTime                                ! Amount of time between screen status messages (sec)
   !integer                                        :: n_SttsTime                              ! Number of time steps between screen status messages (-)
   integer                                        :: nt, iCase
   logical                                        :: AD_Initialized
   
                            

   errStat     = ErrID_None
   errMsg      = ''
   AD_Initialized = .false.
   
   time        = 0.0 ! seconds
      
            
      ! Get the current time
   !call date_and_time ( Values=StrtTime )                               ! Let's time the whole simulation
   !call cpu_time ( UsrTime1 )                                           ! Initial time (this zeros the start time when used as a MATLAB function)
   
   
      ! initialize this driver:
   call Dvr_Init( DvrData, ErrStat, ErrMsg)
      call CheckError()
   
     ! figure out how many time steps we should go before writing screen output:      
    !n_SttsTime = MAX( 1, NINT( SttsTime / dT ) )
                      
   
   !call SimStatus_FirstTime( TiLstPrn, PrevClockTime, SimStrtTime, UsrTime2, time, TMax )
   !   if (ErrStat >= AbortErrLev) call WTP_DvrCleanup
   
   do iCase = 1, DvrData%NumCases
      call WrScr( 'Running case '//trim(num2lstr(iCase))//' of '//trim(num2lstr(DvrData%NumCases))//'.' )
      
      !dT = TwoPi/DvrData%Cases(iCase)%RotSpeed / DvrData%NumSect ! sec
      
      numSteps = ceiling( DvrData%Cases(iCase)%TMax / DvrData%Cases(iCase)%dT)      
      dT_Dvr   = DvrData%Cases(iCase)%dT
      
         ! Set the Initialization input data for AeroDyn based on the Driver input file data, and initialize AD
         !bjj: I'm doing this inside the loop because DT is changing (and don't we have to reinitialize states?); 
         ! TODO: check with GJH/JMJ/RRD about this change
      call Init_AeroDyn(DvrData, AD, dT_Dvr, errStat, errMsg)
         call CheckError()
         AD_Initialized = .true.
         
         if (.not. EqualRealNos( dT_Dvr, DvrData%Cases(iCase)%dT ) ) then
            ErrStat = ErrID_Fatal
            ErrMsg = 'AeroDyn changed the time step for case '//trim(num2lstr(iCase))//'. Change DTAero to "default".'
            call CheckError()
         end if

         ! set TSR and/or WS here (RotorRad comes from Init_AeroDyn):
      IF ( DvrData%InputTSR )  THEN
         DvrData%Cases(iCase)%WndSpeed = DvrData%RotorRad*DvrData%Cases(iCase)%RotSpeed/DvrData%Cases(iCase)%TSR
      ELSE
         DvrData%Cases(iCase)%TSR      = DvrData%RotorRad*DvrData%Cases(iCase)%RotSpeed/DvrData%Cases(iCase)%WndSpeed
      ENDIF
                        
      call Dvr_InitializeOutputFile( iCase, DvrData%Cases(iCase), DvrData%OutFileData, errStat, errMsg)
         call CheckError()
      
      
      do nt = 1, numSteps
         
         !...............................
         ! set AD inputs
         !...............................
         
         call Set_AD_Inputs(iCase,nt,time,DvrData,AD,errStat,errMsg)
            call CheckError()
   
            ! Get state variables at next step: INPUT at step n, OUTPUT at step n + 1

         call AD_UpdateStates( time, n, AD%u, AD%InputTime, AD%p, AD%x, AD%xd, AD%z, AD%OtherState, errStat, errMsg )
            call CheckError()
      
      
            ! Calculate outputs at n

         call AD_CalcOutput( time, AD%u(1), AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%y, errStat, errMsg )
            call CheckError()
   
         call Dvr_WriteOutputLine(DvrData%OutFileData, time, AD%y%WriteOutput, errStat, errMsg)
            call CheckError()
            
      end do !nt=1,numSteps
      
      call AD_End( AD%u(1), AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%y, errStat, errMsg )
         AD_Initialized = .false.         
         call CheckError()
         close( DvrData%OutFileData%unOutFile )
               
   end do !iCase = 1, DvrData%NumCases
   
   
   call Dvr_End()
   
contains
!................................   
   subroutine CheckError()
   
      if (ErrStat /= ErrID_None) then
         call WrScr(TRIM(ErrMsg))
         
         if (ErrStat >= AbortErrLev) then
            call Dvr_End()
         end if
      end if
         
   end subroutine CheckError
!................................   
   subroutine Dvr_End()
   
         ! Local variables
      character(ErrMsgLen)                          :: errMsg2                 ! temporary Error message if ErrStat /= ErrID_None
      integer(IntKi)                                :: errStat2                ! temporary Error status of the operation
      
      character(*), parameter                       :: RoutineName = 'Dvr_End'
         ! Close the output file
      if (DvrData%OutFileData%unOutFile > 0) close(DvrData%OutFileData%unOutFile)
            
      if ( AD_Initialized ) then
         call AD_End( AD%u(1), AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%y, errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      end if
           
      call AD_Dvr_DestroyDvr_SimData( DvrData, ErrStat2, ErrMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

      call AD_Dvr_DestroyAeroDyn_Data( AD, ErrStat2, ErrMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
               
      if (ErrStat >= AbortErrLev) then      
         CALL ProgAbort( 'AeroDyn Driver encountered simulation error level: '&
             //TRIM(GetErrStr(ErrStat)), TrapErrors=.FALSE., TimeWait=3._ReKi )  ! wait 3 seconds (in case they double-clicked and got an error)
      else
         call NormStop()
      end if
      
      
   end subroutine Dvr_End
!................................   
end program AeroDyn_Driver
   