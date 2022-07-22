!**********************************************************************************************************************************
! SS_Radiation_DriverCode: This code tests the template modules
!..................................................................................................................................
! LICENSING
! Copyright (C) 2012  National Renewable Energy Laboratory
!
!    This file is part of SS_Radiation.
!
!    SS_Radiation is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with SS_Radiation.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
PROGRAM SS_Radiation_Driver

   USE NWTC_Library
   USE SS_Radiation
   USE SS_Radiation_Types

   IMPLICIT NONE

      ! Program variables

   REAL(DbKi)                                         :: Time                 ! Variable for storing time, in seconds
   REAL(DbKi)                                         :: Time2(145201,1)                 ! Variable for storing time, in seconds
   REAL(DbKi)                                         :: tdq(145201,7)           ! Variable for storing time and body velocities, in m/s or rad/s
   REAL(DbKi)                                         :: dq(145201,6)           ! Variable for storing body velocities, in m/s or rad/s
   REAL(DbKi)                                         :: TimeInterval         ! Interval between time steps, in seconds
   INTEGER(B1Ki), ALLOCATABLE                         :: SaveAry(:)           ! Array to store packed data structure

   TYPE(SS_Rad_InitInputType)                        :: InitInData           ! Input data for initialization
   TYPE(SS_Rad_InitOutputType)                       :: InitOutData          ! Output data from initialization

   TYPE(SS_Rad_ContinuousStateType)                  :: x                    ! Continuous states
   TYPE(SS_Rad_ContinuousStateType)                  :: x_new                ! Continuous states at updated time
   TYPE(SS_Rad_DiscreteStateType)                    :: xd                   ! Discrete states
   TYPE(SS_Rad_DiscreteStateType)                    :: xd_new               ! Discrete states at updated time
   TYPE(SS_Rad_ConstraintStateType)                  :: z                    ! Constraint states
   TYPE(SS_Rad_ConstraintStateType)                  :: z_residual           ! Residual of the constraint state equations (Z)
   TYPE(SS_Rad_OtherStateType)                       :: OtherState           ! Other states

   TYPE(SS_Rad_ParameterType)                        :: p                    ! Parameters
   TYPE(SS_Rad_InputType)                            :: u                    ! System inputs
   TYPE(SS_Rad_OutputType)                           :: y                    ! System outputs
   TYPE(SS_Rad_MiscVarType)                          :: m                    ! misc/optimization variables

   TYPE(SS_Rad_ContinuousStateType)                  :: dxdt                 ! First time derivatives of the continuous states



    !Local Variables
   INTEGER(IntKi)                                     :: n                    ! Loop counter (for time step)
   INTEGER(IntKi)                                     :: I                    ! Loop counter (for time step)
   INTEGER(IntKi)                                     :: J                    ! Loop counter (for time step)
   INTEGER(IntKi)                                     :: Inputdq              ! Input file identifier
   INTEGER(IntKi)                                     :: Outputy              ! Output file identifier
   INTEGER(IntKi)                                     :: ErrStat              ! Status of error message
   CHARACTER(1024)                                    :: ErrMsg               ! Error message if ErrStat /= ErrID_None
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

   !...............................................................................................................................
   ! Routines called in initialization
   !...............................................................................................................................

   ! Call Time 
   CALL cpu_time(start)
   CALL DATE_AND_TIME ( Values=StrtTime )  
   
   ! Populate the InitInData data structure here:

    InitInData%InputFile = 'C:\Users\tduarte\Documents\SS_Module\Comparisons\FAST_output_freq\spar_IMP_097'
    !!! GREG !!!: This file name should be the WAMIT file name without extension!


   InitInData%Dofs = 1
    !!! GREG: This is a vector of [1x6] containing 0 and 1 if each of the 6 dofs is enabled or not (as we discussed today in the meeting)
   
   
   ! Set the driver's request for time interval here:
   TimeInterval = 0.025                     ! Glue code's request for delta time (likely based on information from other modules)
    !!! GREG: This should be the Rdtn DT defined in the platform input file@

   CALL SS_Rad_Init( InitInData, u, p,  x, xd, z, OtherState, y, m, TimeInterval, InitOutData, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
   END IF
   
   !!! GREG: This version reads in the desired file containing the platform velocities. You don't need this in your case.
   CALL CheckArgs( InitInData%InputFile )
   
   CALL Get_Arg_Num (len )
     
   ! Read the time dependent input vector dq
   CALL OpenFInpFile ( Inputdq,  (TRIM(InitInData%InputFile)//'.txt'), ErrStat   )  ! Open motion file.
    IF ( ErrStat /= 0 ) THEN
        ErrStat = ErrID_Fatal
        ErrMsg  = ' Error allocating memory for the dq array.'
        print*, ( ErrMsg )
    END IF

    
    
    DO I = 1,145201 !Read dq Matrix
        READ (Inputdq,*,IOSTAT=Sttus) (tdq (I,J), J=1,7)   
    ENDDO  
    
    CLOSE ( Inputdq ) !Close dq input file
    
    Time2(:,1) = tdq(:,1)
    dq = tdq(:,2:7)
      
    !!!GREG: here the output file is opened, you should not need this
      !Initialize output file
    CALL OpenFOutFile ( Outputy, (TRIM(InitInData%InputFile)//'.out'), ErrStat)
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
  
   CALL WrScr( 'Runnig SS_Radiation in Loose Coupling using a Adams-Bashforth-Moulton Method'  ) 

   CALL SS_Rad_CopyDiscState( xd, xd_new, MESH_NEWCOPY, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
   END IF
   

   CALL SS_Rad_CopyContState( x, x_new, MESH_NEWCOPY, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
   END IF
   !
!CALL cpu_time(T1)
   DO n = 0,145200
   
      Time = n*TimeInterval
   
         ! Modify u (likely from the outputs of another module or a set of test conditions) here:
   
        u%dq(1,1) = dq (n+1,1)
        u%dq(2,1) = dq (n+1,2)
        u%dq(3,1) = dq (n+1,3)
        u%dq(4,1) = dq (n+1,4)
        u%dq(5,1) = dq (n+1,5)
        u%dq(6,1) = dq (n+1,6)

         ! Get state variables at next step: constraint states (z) at step n, continuous and discrete states at step n + 1
   
      CALL SS_Rad_UpdateStates( Time, u, p, x_new, xd_new, z, OtherState, m, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
      END IF
   !print*, x%x
         ! Calculate outputs at n
   
      CALL SS_Rad_CalcOutput( Time, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
      END IF
   
         ! Update x and xd with continuous and discrete states at n + 1
         ! Note that the constraint state guess at n+1 is the value of the constraint state at n (so it doesn't need updating here)
     
      CALL SS_Rad_CopyContState( x_new, x, MESH_UPDATECOPY, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
      END IF

      CALL SS_Rad_CopyDiscState( xd_new, xd, MESH_UPDATECOPY, ErrStat, ErrMsg )
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
   
   
   CALL SS_Rad_DestroyDiscState( xd_new, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
   END IF
   
   CALL SS_Rad_DestroyContState( x_new, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
   END IF
   

   !...............................................................................................................................
   ! Routine to terminate program execution
   !...............................................................................................................................
   CALL SS_Rad_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN
      CALL WrScr( ErrMsg )
   END IF
   

   !!! GREG: This is also to ouput values (dont need it)
   CALL DATE_AND_TIME ( VALUES=EndTimes )
   CALL cpu_time(finnish)
   
   ClckTime =  0.001*( EndTimes(8) - StrtTime(8) ) + ( EndTimes(7) - StrtTime(7) ) + 60.0*( EndTimes(6) - StrtTime(6) ) &
         + 3600.0*( EndTimes(5) - StrtTime(5) ) + 86400.0*( EndTimes(3) - StrtTime(3) )  

   UsrTime = finnish-start

   IF ( UsrTime /= 0.0 )  THEN

   TRatio = Time / UsrTime

   IF     ( UsrTime > 86400.0 )  THEN
      Factor = 1.0/86400.0
      TimePer = ' days'
   ELSEIF ( UsrTime >  3600.0 )  THEN
      Factor = 1.0/3600.0
      TimePer = ' hours'
   ELSEIF ( UsrTime >    60.0 )  THEN
      Factor = 1.0/60.0
      TimePer = ' minutes'
   ELSE
      Factor = 1.0
      TimePer = ' seconds'
   ENDIF

   CALL WrScr ( ' Total Real Time:       '//TRIM( Flt2LStr( Factor*ClckTime      ) )//TRIM( TimePer ) )
   CALL WrScr ( ' Total CPU Time:        '//TRIM( Flt2LStr( Factor*UsrTime       ) )//TRIM( TimePer ) )
   CALL WrScr ( ' Simulated Time:        '//TRIM( Flt2LStr( Factor*REAL( Time ) ) )//TRIM( TimePer ) )
   CALL WrScr ( ' Time Ratio (Sim/CPU):  '//TRIM( Flt2LStr( TRatio ) ) )

   ENDIF
   
   
    !Write Output to file
      WRITE(Outputy,'(1(e16.6))',IOSTAT=Sttus)  TRatio
         ! Ending routines
   CLOSE( Outputy )



END PROGRAM SS_Radiation_Driver

