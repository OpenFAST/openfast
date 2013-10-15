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
   TYPE(SS_Rad_OtherStateType)                       :: OtherState           ! Other/optimization states

   TYPE(SS_Rad_ParameterType)                        :: p                    ! Parameters
   TYPE(SS_Rad_InputType)                            :: u                    ! System inputs
   TYPE(SS_Rad_OutputType)                           :: y                    ! System outputs

   TYPE(SS_Rad_ContinuousStateType)                  :: dxdt                 ! First time derivatives of the continuous states

   !TYPE(SS_Rad_PartialOutputPInputType)              :: dYdu                 ! Partial derivatives of the output equations
                                                                              !  (Y) with respect to the inputs (u)
  ! TYPE(SS_Rad_PartialContStatePInputType)           :: dXdu                 ! Partial derivatives of the continuous state
                                                                              !  equations (X) with respect to the inputs (u)
  ! TYPE(SS_Rad_PartialDiscStatePInputType)           :: dXddu                ! Partial derivatives of the discrete state
                                                                              !  equations (Xd) with respect to the inputs (u)
   !TYPE(SS_Rad_PartialConstrStatePInputType)         :: dZdu                 ! Partial derivatives of the constraint state
                                                                              !  equations (Z) with respect to the inputs (u)
   !TYPE(SS_Rad_PartialOutputPContStateType)          :: dYdx                 ! Partial derivatives of the output equations (Y)
                                                                              !  with respect to the continuous states (x)
!   TYPE(SS_Rad_PartialContStatePContStateType)       :: dXdx                 ! Partial derivatives of the continuous state equat-
                                                                              !  ions (X) with respect to the continuous states (x)
   !TYPE(SS_Rad_PartialDiscStatePContStateType)       :: dXddx                ! Partial derivatives of the discrete state equat-
                                                                              !  ions (Xd) with respect to continuous states (x)
   !TYPE(SS_Rad_PartialConstrStatePContStateType)     :: dZdx                 ! Partial derivatives of the constraint state equat-
                                                                              !  ions (Z) with respect to the continuous states (x)
   !TYPE(SS_Rad_PartialOutputPDiscStateType)          :: dYdxd                ! Partial derivatives of the output equations (Y)
                                                                              !  with respect to the discrete states (xd)
  ! TYPE(SS_Rad_PartialContStatePDiscStateType)       :: dXdxd                ! Partial derivatives of the continuous state equat-
                                                                              !  ions (X) with respect to the discrete states (xd)
 !  TYPE(SS_Rad_PartialDiscStatePDiscStateType)       :: dXddxd               ! Partial derivatives of the discrete state equat-
                                                                              !  ions (Xd) with respect to the discrete states (xd)
  ! TYPE(SS_Rad_PartialConstrStatePDiscStateType)     :: dZdxd                ! Partial derivatives of the constraint state equat-
                                                                              !  ions (Z) with respect to the discrete states (xd)
   !TYPE(SS_Rad_PartialOutputPConstrStateType)        :: dYdz                 ! Partial derivatives of the output equations (Y)
                                                                              !  with respect to the constraint states (z)
  ! TYPE(SS_Rad_PartialContStatePConstrStateType)     :: dXdz                 ! Partial derivatives of the continuous state equat-
                                                                              !  ions (X) with respect to the constraint states (z)
 !  TYPE(SS_Rad_PartialDiscStatePConstrStateType)     :: dXddz                ! Partial derivatives of the discrete state equat-
                                                                              !  ions (Xd) with respect to constraint states (z)
!   TYPE(SS_Rad_PartialConstrStatePConstrStateType)   :: dZdz                 ! Partial derivatives of the constraint state equat-
                                                                              !  ions (Z) with respect to the constraint states (z)

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

   CALL SS_Rad_Init( InitInData, u, p,  x, xd, z, OtherState, y, TimeInterval, InitOutData, ErrStat, ErrMsg )
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
   
      CALL SS_Rad_UpdateStates( Time, u, p, x_new, xd_new, z, OtherState, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
      END IF
   !print*, x%x
         ! Calculate outputs at n
   
      CALL SS_Rad_CalcOutput( Time, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
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
   

!   !...............................................................................................................................
!   ! Routines called in tight coupling -- time marching only
!   !...............................................................................................................................
!
!   CALL WrScr( 'Runnig SS_Radiation in Loose Coupling using a Adams-Bashforth-Moulton Method'  ) 
!   DO n=0,145200
!
!      Time = n * TimeInterval   ! Note that the discrete states must be updated only at the TimeInterval defined in initialization
!
!         ! set inputs (u) here:
!        u%dq(1,1) = dq (n+1,1)
!        u%dq(2,1) = dq (n+1,2)
!        u%dq(3,1) = dq (n+1,3)
!        u%dq(4,1) = dq (n+1,4)
!        u%dq(5,1) = dq (n+1,5)
!        u%dq(6,1) = dq (n+1,6)
!
!         ! Update constraint states at Time
!      !
!      !CALL SS_Rad_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_residual, ErrStat, ErrMsg )
!      !IF ( ErrStat /= ErrID_None ) THEN      ! Check if there was an error and do something about it if necessary
!      !   CALL WrScr( ErrMsg )
!      !END IF
!
!      ! DO WHILE ( z_residual% > tolerance )
!
!         ! z =
!
!!         CALL SS_Rad_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_residual, ErrStat, ErrMsg )
!!         IF ( ErrStat /= ErrID_None ) THEN      ! Check if there was an error and do something about it if necessary
!!            CALL WrScr( ErrMsg )
!!         END IF
!
!      ! END DO
!
!
!
!         ! Calculate the outputs at Time
!
!      CALL SS_Rad_CalcOutput( Time, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
!      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
!         CALL WrScr( ErrMsg )
!      END IF
!
!
!         ! Calculate the continuous state derivatives at Time
!
!      CALL SS_Rad_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )
!      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
!         CALL WrScr( ErrMsg )
!      END IF
!
!
!         ! Update the discrete state from step n to step n+1
!         ! Note that the discrete states must be updated only at the TimeInterval defined in initialization
!
!      !CALL SS_Rad_UpdateDiscState( Time, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!      !IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
!      !   CALL WrScr( ErrMsg )
!      !END IF
!
!
!         ! Driver should integrate (update) continuous states here:
!
!     ! 4th order solver
!      
!     !IF ( EqualRealNos(OtherState%LastTime + p%DT,Time)) THEN !Time must have been reduced TD: Function EqualRealNos only works for Single!
!     IF ( OtherState%LastTime + p%DT/=Time) THEN !Time must have been reduced
!         
!         OtherState%dxdt = 0 ! Remove previous history and start from zero with a
!                              ! runge-Kutta method
!         OtherState%Step = 0
!     ENDIF
!     
!    !Update time step
!     OtherState%Step = OtherState%Step + 1 
!    !Update the OtherStates matrices, with the previous dXdt Values
!        OtherState%dxdt (:,4) = OtherState%dxdt (:,3)
!        OtherState%dxdt (:,3) = OtherState%dxdt (:,2)
!        OtherState%dxdt (:,2) = OtherState%dxdt (:,1)
!        OtherState%dxdt (:,1) = dxdt%x (:,1)
!    
!     Call Solver (Time, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg)
!       
!     !Update LastTime
!     OtherState%LastTime = Time
!        
!
!
!         ! Jacobians required:
!
!      CALL SS_Rad_JacobianPInput( Time, u, p, x, xd, z, OtherState, dYdu=dYdu, dZdu=dZdu, ErrStat=ErrStat, ErrMsg=ErrMsg )
!      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
!         CALL WrScr( ErrMsg )
!      END IF
!
!      CALL SS_Rad_JacobianPConstrState( Time, u, p, x, xd, z, OtherState, dYdz=dYdz, dZdz=dZdz, ErrStat=ErrStat, ErrMsg=ErrMsg )
!      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
!         CALL WrScr( ErrMsg )
!      END IF
!
!
!   END DO
!
!
!      ! Destroy z_residual and dxdt because they are not necessary anymore
!
!   CALL SS_Rad_DestroyConstrState( z_residual, ErrStat, ErrMsg )
!   IF ( ErrStat /= ErrID_None ) THEN   ! Check if there was an error and do something about it if necessary
!      CALL WrScr( ErrMsg )
!   END IF
!
!   CALL SS_Rad_DestroyContState( dxdt, ErrStat, ErrMsg )
!   IF ( ErrStat /= ErrID_None ) THEN   ! Check if there was an error and do something about it if necessary
!      CALL WrScr( ErrMsg )
!   END IF
!
!   !...............................................................................................................................
!   ! Jacobian routines called in tight coupling
!   !...............................................................................................................................
!
!   CALL SS_Rad_JacobianPInput( Time, u, p, x, xd, z, OtherState, dYdu, dXdu, dXddu, dZdu, ErrStat, ErrMsg )
!   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
!      CALL WrScr( ErrMsg )
!   END IF
!
!   CALL SS_Rad_JacobianPContState( Time, u, p, x, xd, z, OtherState, dYdx, dXdx, dXddx, dZdx, ErrStat, ErrMsg )
!   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
!      CALL WrScr( ErrMsg )
!   END IF
!
!   CALL SS_Rad_JacobianPDiscState( Time, u, p, x, xd, z, OtherState, dYdxd, dXdxd, dXddxd, dZdxd, ErrStat, ErrMsg )
!   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
!      CALL WrScr( ErrMsg )
!   END IF
!
!   CALL SS_Rad_JacobianPConstrState( Time, u, p, x, xd, z, OtherState, dYdz, dXdz, dXddz, dZdz, ErrStat, ErrMsg )
!   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
!      CALL WrScr( ErrMsg )
!   END IF

!
!   !!...............................................................................................................................
!   !! Routines to pack data (to restart later)
!   !!...............................................................................................................................
!   !CALL SS_Rad_Pack(SaveAry, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg)
!   !IF ( ErrStat /= ErrID_None ) THEN
!   !   CALL WrScr( ErrMsg )
!   !END IF
!
!
   !...............................................................................................................................
   ! Routine to terminate program execution
   !...............................................................................................................................
   CALL SS_Rad_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
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

!   !...............................................................................................................................
!   ! Routines to retreive packed data (unpack for restart)
!   !...............................................................................................................................
!   !CALL SS_Rad_Unpack( SaveAry, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
!   !
!   !IF ( ErrStat /= ErrID_None ) THEN
!   !   CALL WrScr( ErrMsg )
!   !END IF
!   !
!
!   !...............................................................................................................................
!   ! Routines to copy data (not already tested)
!   !...............................................................................................................................
!
!
!
!   !...............................................................................................................................
!   ! Routines to destroy data (not already tested)
!   !...............................................................................................................................
!
!   IF ( ALLOCATED( SaveAry ) ) DEALLOCATE( SaveAry )
!
!!   CALL SS_Rad_DestroyPartialOutputPInput ( )  % Jacobian Routine not yet implemented
!
!
!
!   !...............................................................................................................................
!   ! Routine to terminate program execution (again)
!   !...............................................................................................................................
!
!   CALL SS_Rad_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
!   IF ( ErrStat /= ErrID_None ) THEN
!      CALL WrScr( ErrMsg )
!   END IF


END PROGRAM SS_Radiation_Driver

