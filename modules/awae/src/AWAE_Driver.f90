!  AWAE_Driver.f90 
!
!  FUNCTIONS:
!  AWAE_Driver - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: AWAE_Driver
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

program AWAE_Driver
   use AWAE
   use AWAE_Driver_Subs
   use NWTC_Library
      
   implicit none

   

   character(1024)             :: AWAE_Dvr_InputFile                   !< Name of the file containing the driver input data
   type(AWAE_InitInputType)   :: AWAE_InitInp                         !< Input data for initialization routine
   type(AWAE_InitOutputType)     :: AWAE_InitOut                         !< Input data for initialization routine
   type(AWAE_InputType)       :: AWAE_u                               !< Input data for initialization routine
   type(AWAE_ParameterType)   :: AWAE_p                               !< Input data for initialization routine
   type(AWAE_DiscreteStateType)   :: AWAE_xd                              !< Input data for initialization routine
   type(AWAE_OutputType)      :: AWAE_y                               !< Input data for initialization routine
                                                                      
   integer(IntKi)           :: errStat                            !< Error status
   character(ErrMsgLen)     :: errMsg                             !< Error message
   
   errStat     = ErrID_None
   errMsg      = ''
   
      ! Initialize the NWTC Subroutine Library
   call NWTC_Init( EchoLibVer=.FALSE. )
   !call AWAE_TEST_ExtractSlice(errStat, errMsg)
   call AWAE_TEST_LowResGridCalcs(errStat, errMsg)
!   call AWAE_Dvr_Tests(1, errStat, errMsg)
   call CheckError( errStat, errMsg  )
      ! Initialize the Driver and the WD module
   !call AWAE_TEST_Init_BadData(errStat, ErrMsg)
   !call CheckError( ErrStat, ErrMsg  )
   
  ! call AWAE_TEST_Init_GoodData(errStat, ErrMsg)
  ! call CheckError( ErrStat, ErrMsg  )
   !call AWAE_Dvr_Init( AWAE_InitInp,  AWAE_InitOut, AWAE_u,AWAE_p, AWAE_xd,  AWAE_y, errStat, errMsg)
  ! call AWAE_TEST_CalcOutput(errStat, errMsg) 
  ! call CheckError( ErrStat, ErrMsg  )
      ! Run the time marching loop
  ! call AWAE_Dvr_Time_Marching()
    
      ! Cleanup 
 !  call AWAE_Dvr_Cleanup()
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
end program AWAE_Driver

