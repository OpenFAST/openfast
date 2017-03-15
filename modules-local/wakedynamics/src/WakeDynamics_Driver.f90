!  WakeDynamics_Driver.f90 
!
!  FUNCTIONS:
!  WakeDynamics_Driver - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: WakeDynamics_Driver
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

program WakeDynamics_Driver
   use WakeDynamics
   use WD_Driver_Subs
   use NWTC_Library
      
   implicit none

   

   character(1024)             :: WD_Dvr_InputFile                   !< Name of the file containing the driver input data
   type(WD_InitInputType)   :: WD_InitInp                         !< Input data for initialization routine
   type(WD_InitOutputType)     :: WD_InitOut                         !< Input data for initialization routine
   type(WD_InputType)       :: WD_u                               !< Input data for initialization routine
   type(WD_ParameterType)   :: WD_p                               !< Input data for initialization routine
   type(WD_DiscreteStateType)   :: WD_xd                              !< Input data for initialization routine
   type(WD_OutputType)      :: WD_y                               !< Input data for initialization routine
                                                                      
   integer(IntKi)           :: errStat                            !< Error status
   character(ErrMsgLen)     :: errMsg                             !< Error message
   
   errStat     = ErrID_None
   errMsg      = ''
   
      ! Initialize the NWTC Subroutine Library
   call NWTC_Init( EchoLibVer=.FALSE. )

      ! Initialize the Driver and the WD module
  ! call WD_TEST_Init_BadData(errStat, ErrMsg)
  ! CALL CheckError( ErrStat, ErrMsg  )
   
  ! call WD_TEST_Init_GoodData(errStat, ErrMsg)
  ! CALL CheckError( ErrStat, ErrMsg  )
   !call WD_Dvr_Init( WD_InitInp,  WD_InitOut, WD_u,WD_p, WD_xd,  WD_y, errStat, errMsg)
   call WD_TEST_UpdateStates(errStat, errMsg) 
   CALL CheckError( ErrStat, ErrMsg  )
      ! Run the time marching loop
  ! call WD_Dvr_Time_Marching()
    
      ! Cleanup 
 !  call WD_Dvr_Cleanup()
 CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg,ErrLocMsg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN)           :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN)           :: Msg         ! The error message (ErrMsg)
      CHARACTER(*),   INTENT(IN), OPTIONAL :: ErrLocMsg   ! an optional message describing the location of the error

     ! CHARACTER(1024)                      :: SimMsg      
     ! integer(IntKi)                       :: i_turb2
      
      
      IF ( ErrID /= ErrID_None ) THEN
         CALL WrScr( NewLine//TRIM(Msg)//NewLine )
         !IF ( ErrID >= AbortErrLev ) THEN
         !   
         !   IF (PRESENT(ErrLocMsg)) THEN
         !      SimMsg = ErrLocMsg
         !   ELSE
         !      SimMsg = 'at simulation time '//TRIM(Num2LStr(Turbine(1)%m_FAST%t_global))//' of '//TRIM(Num2LStr(Turbine(1)%p_FAST%TMax))//' seconds'
         !   END IF
         !   
         !   DO i_turb2 = 1,NumTurbines
         !      CALL ExitThisProgram_T( Turbine(i_turb2), ErrID, i_turb2==NumTurbines, SimMsg )
         !   END DO
         !               
         !END IF
         
      END IF


   END SUBROUTINE CheckError      
end program WakeDynamics_Driver

