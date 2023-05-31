!**********************************************************************************************************************************
! SS_Excitation_Driver: This code tests the SS_Excitation module
!..................................................................................................................................
! LICENSING
! Copyright (C) 2018  National Renewable Energy Laboratory
!
!    This file is part of SS_Excitation.
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
!    
!**********************************************************************************************************************************
PROGRAM SS_Excitation_Driver

   USE NWTC_Library
   USE SS_Excitation
   USE SS_Excitation_Types

   IMPLICIT NONE

      ! Program variables

   REAL(DbKi)                                         :: Time                 ! Variable for storing time, in seconds
   REAL(DbKi)                                         :: waveDT
   !REAL(DbKi)                                         :: Time2(145201,1)                 ! Variable for storing time, in seconds
   !REAL(DbKi)                                         :: tdq(145201,7)           ! Variable for storing time and body velocities, in m/s or rad/s
   !REAL(DbKi)                                         :: dq(145201,6)           ! Variable for storing body velocities, in m/s or rad/s
   REAL(DbKi)                                         :: TimeInterval         ! Interval between time steps, in seconds
   !INTEGER(B1Ki), ALLOCATABLE                         :: SaveAry(:)           ! Array to store packed data structure

   TYPE(SS_Exc_InitInputType)                        :: InitInData           ! Input data for initialization
   TYPE(SS_Exc_InitOutputType)                       :: InitOutData          ! Output data from initialization

   TYPE(SS_Exc_ContinuousStateType)                  :: x                    ! Continuous states
   TYPE(SS_Exc_ContinuousStateType)                  :: x_new                ! Continuous states at updated time
   TYPE(SS_Exc_DiscreteStateType)                    :: xd                   ! Discrete states
   TYPE(SS_Exc_DiscreteStateType)                    :: xd_new               ! Discrete states at updated time
   TYPE(SS_Exc_ConstraintStateType)                  :: z                    ! Constraint states
   TYPE(SS_Exc_ConstraintStateType)                  :: z_residual           ! Residual of the constraint state equations (Z)
   TYPE(SS_Exc_OtherStateType)                       :: OtherState           ! Other states

   TYPE(SS_Exc_ParameterType)                        :: p                    ! Parameters
   TYPE(SS_Exc_InputType)                            :: u(1)                 ! System inputs
   REAL(DbKi)                                        :: InputTimes(1)        ! System input times
   TYPE(SS_Exc_OutputType)                           :: y                    ! System outputs
   TYPE(SS_Exc_MiscVarType)                          :: m                    ! misc/optimization variables

   TYPE(SS_Exc_ContinuousStateType)                  :: dxdt                 ! First time derivatives of the continuous states



    !Local Variables
   INTEGER(IntKi)                                     :: n                    ! Loop counter (for time step)
   INTEGER(IntKi)                                     :: I                    ! Loop counter (for time step)
   INTEGER(IntKi)                                     :: J                    ! Loop counter (for time step)
   REAL(SiKi)                                         :: ElevData
   INTEGER(IntKi)                                     :: UnWvEl               ! Input file identifier
   INTEGER(IntKi)                                     :: Outputy              ! Output file identifier
   INTEGER(IntKi)                                     :: ErrStat, ErrStat2    ! Status of error message
   CHARACTER(1024)                                    :: ErrMsg, ErrMsg2      ! Error message if ErrStat /= ErrID_None
   INTEGER                                            :: Sttus                ! Error in reading input file
   REAL(ReKi)                                         :: Start                ! CPU Time at start of the program
   REAL(ReKi)                                         :: Finnish              ! CPU Time at the end of the program
   REAL(ReKi)                                         :: UsrTime 
   REAL(ReKi)                                         :: Tratio 
   REAL(ReKi)                                         :: Factor
   CHARACTER(8)                                       :: TimePer
   INTEGER(4)                                         :: EndTimes (8)         ! An array holding the ending clock time of the simulation.
   INTEGER(4)                                         :: StrtTime (8)         ! An array holding the starting clock time of the simulation.
   REAL(ReKi)                                         :: ClckTime               
   INTEGER                                            :: len                  ! Number of input arguments
   CHARACTER(1024)                                    :: waveFile

   !...............................................................................................................................
   ! Routines called in initialization
   !...............................................................................................................................

   ErrStat = ErrID_None
   ErrMsg  = ''
   
   call NWTC_Init()
      
   ! Call Time 
   !call cpu_time(start)
   !call DATE_AND_TIME ( Values=StrtTime ) 
   

      
      ! Populate the InitInData data structure
     
   
      ! This file name should be the WAMIT file name without extension!
   InitInData%InputFile = 'C:\Dev\Envision\all-changes\Test_Models\5MW_Baseline\HydroData\barge'
   InitInData%WaveDir   = 0.0_ReKi 
   InitInData%NStepWave = 14520
   waveDT = 0.25
   allocate ( InitInData%WaveElev0(0:InitInData%NStepWave) , STAT=ErrStat2 )
   allocate ( InitInData%WaveTime (0:InitInData%NStepWave) , STAT=ErrStat2 )
   
      ! Construct the wave times array 
   do i = 0,InitInData%NStepWave
      InitInData%WaveTime(i) = waveDT*i
   end do
   
      ! Need to read in the wave elevation data to pass in as initialization data
   waveFile = 'C:\Dev\Envision\all-changes\Test_Models\5MW_ITIBarge_DLL_WTurb_WavesIrr\barge.Elev'
   call GetNewUnit ( UnWvEl, ErrStat, ErrMsg )
   call OpenFInpFile ( UnWvEl,  trim(waveFile), ErrStat, ErrMsg   )  ! Open wave elevation file.
      if ( ErrStat /= 0 ) then
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Could not open wave elevation file.'
         print*, ( ErrMsg )
      end if
      
   call ReadCom ( UnWvEl, trim(waveFile), 'Header',ErrStat2, ErrMsg2  )! Reads the first entire line (Title header)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Excitation_Driver')  
         
   do i = 0,InitInData%NStepWave - 1
      call ReadVar( UnWvEl,trim(waveFile), InitInData%WaveElev0(i), 'InitInData%WaveElev0(i)', 'Wave elevation',ErrStat2, ErrMsg2) ! Reads in the third line, containing the number of states
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Excitation_Driver')        
   end do  
    
   close ( UnWvEl ) !Close dq input file
    
      ! Now set the last element of the Wave elevation array to match the initial elevation for wrapping
    InitInData%WaveElev0(InitInData%NStepWave) = InitInData%WaveElev0(0)
    
   
   
      ! Set the driver's request for time interval here: This should be the Rdtn DT defined in the hydrodyn input file
   TimeInterval = 0.005

   CALL SS_Exc_Init( InitInData, u(1), p,  x, xd, z, OtherState, y, m, TimeInterval, InitOutData, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
      END IF
   
  

      ! Initialize output file
   call GetNewUnit ( Outputy, ErrStat, ErrMsg )
   CALL OpenFOutFile ( Outputy, (TRIM(InitInData%InputFile)//'.out'), ErrStat, ErrMsg)
      IF ( ErrStat /= 0 ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error opening output file.'
         CALL WrScr( ErrMsg )
      END IF
                      
   WRITE(Outputy,*,IOSTAT=Sttus)  InitOutData%WriteOutputHdr
      IF ( Sttus /= 0 )  THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error writing output file.'
         CALL WrScr( ErrMsg )
      ENDIF
      
   WRITE(Outputy,*,IOSTAT=Sttus)  InitOutData%WriteOutputUnt  
         IF ( Sttus /= 0 )  THEN
            ErrStat = ErrID_Fatal
            ErrMsg  = ' Error writing output file.'
            CALL WrScr( ErrMsg )
         ENDIF
             
   !...............................................................................................................................
   ! Routines called in loose coupling -- the glue code may implement this in various ways
   !...............................................................................................................................
  
   CALL WrScr( 'Runnig SS_Excitation in Loose Coupling using a Adams-Bashforth-Moulton Method'  ) 

   CALL SS_Exc_CopyDiscState( xd, xd_new, MESH_NEWCOPY, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
   END IF
   

   CALL SS_Exc_CopyContState( x, x_new, MESH_NEWCOPY, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
   END IF
   !

   DO n = 0,InitInData%NStepWave-1
   
      Time = n*TimeInterval
      InputTimes(1) = Time

         ! Get state variables at next step: constraint states (z) at step n, continuous and discrete states at step n + 1
      CALL SS_Exc_UpdateStates( Time, n, u, InputTimes, p, x_new, xd_new, z, OtherState, m, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
      END IF
   !print*, x%x
         ! Calculate outputs at n
   
      CALL SS_Exc_CalcOutput( Time, u(1), p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
      END IF
   
         ! Update x and xd with continuous and discrete states at n + 1
         ! Note that the constraint state guess at n+1 is the value of the constraint state at n (so it doesn't need updating here)
     
      CALL SS_Exc_CopyContState( x_new, x, MESH_UPDATECOPY, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
      END IF

      CALL SS_Exc_CopyDiscState( xd_new, xd, MESH_UPDATECOPY, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
      END IF
   
      !Write Output to file
      WRITE(Outputy,'(7(e16.6))',IOSTAT=Sttus)  y%WriteOutput
      IF ( Sttus /= 0 )  THEN
        ErrStat = ErrID_Fatal
        ErrMsg  = ' Error writing output file.'
        CALL WrScr( ErrMsg )
        print*, ErrMsg
      ENDIF
   END DO
   
   
   CALL SS_Exc_DestroyDiscState( xd_new, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
   END IF
   
   CALL SS_Exc_DestroyContState( x_new, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
   END IF
   
   CALL SS_Exc_DestroyInitInput( InitInData, ErrStat, ErrMsg, DEALLOCATEpointers = .true. ) ! pointers were allocated in this data type, so we need to deallocate them here, too
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
   END IF

   !...............................................................................................................................
   ! Routine to terminate program execution
   !...............................................................................................................................
   CALL SS_Exc_End( u(1), p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN
      CALL WrScr( ErrMsg )
   END IF
   

   !!!! GREG: This is also to ouput values (dont need it)
   !CALL DATE_AND_TIME ( VALUES=EndTimes )
   !CALL cpu_time(finnish)
   !
   !ClckTime =  0.001*( EndTimes(8) - StrtTime(8) ) + ( EndTimes(7) - StrtTime(7) ) + 60.0*( EndTimes(6) - StrtTime(6) ) &
   !      + 3600.0*( EndTimes(5) - StrtTime(5) ) + 86400.0*( EndTimes(3) - StrtTime(3) )  
   !
   !UsrTime = finnish-start
   !
   !IF ( UsrTime /= 0.0 )  THEN
   !
   !TRatio = Time / UsrTime
   !
   !IF     ( UsrTime > 86400.0 )  THEN
   !   Factor = 1.0/86400.0
   !   TimePer = ' days'
   !ELSEIF ( UsrTime >  3600.0 )  THEN
   !   Factor = 1.0/3600.0
   !   TimePer = ' hours'
   !ELSEIF ( UsrTime >    60.0 )  THEN
   !   Factor = 1.0/60.0
   !   TimePer = ' minutes'
   !ELSE
   !   Factor = 1.0
   !   TimePer = ' seconds'
   !ENDIF
   !
   !CALL WrScr ( ' Total Real Time:       '//TRIM( Flt2LStr( Factor*ClckTime      ) )//TRIM( TimePer ) )
   !CALL WrScr ( ' Total CPU Time:        '//TRIM( Flt2LStr( Factor*UsrTime       ) )//TRIM( TimePer ) )
   !CALL WrScr ( ' Simulated Time:        '//TRIM( Flt2LStr( Factor*REAL( Time ) ) )//TRIM( TimePer ) )
   !CALL WrScr ( ' Time Ratio (Sim/CPU):  '//TRIM( Flt2LStr( TRatio ) ) )
   !
   !ENDIF
   
   
    !!Write Output to file
    !  WRITE(Outputy,'(1(e16.6))',IOSTAT=Sttus)  TRatio
    !     ! Ending routines
   
   CLOSE( Outputy )



END PROGRAM SS_Excitation_Driver

