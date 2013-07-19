!**********************************************************************************************************************************
! SubDyn_DriverCode: This code tests the SubDyn modules
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
!
!    This file is part of SubDyn.
!
!    SubDyn is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with SubDyn.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
PROGRAM TestSubDyn

   USE NWTC_Library
   USE SubDyn
   USE TempMod
   USE SubDyn_Types
   USE SubDyn_Output

   IMPLICIT NONE

   INTEGER(IntKi), PARAMETER                          :: NumInp = 1           ! Number of inputs sent to SD_UpdateStates
   
      ! Program variables

   REAL(DbKi)                                         :: Time                 ! Variable for storing time, in seconds
   REAL(DbKi)                                         :: TimeInterval         ! Interval between time steps, in seconds
   REAL(DbKi)                                         :: InputTime(NumInp)    ! Variable for storing time associated with inputs, in seconds
   
   TYPE(SD_InitInputType)                        :: InitInData           ! Input data for initialization
   TYPE(SD_InitOutputType)                       :: InitOutData          ! Output data from initialization

   TYPE(SD_ContinuousStateType)                  :: x                    ! Continuous states
   TYPE(SD_DiscreteStateType)                    :: xd                   ! Discrete states
   TYPE(SD_ConstraintStateType)                  :: z, z_next                    ! Constraint states
   TYPE(SD_ConstraintStateType)                  :: Z_residual           ! Residual of the constraint state functions (Z)
   TYPE(SD_OtherStateType)                       :: OtherState           ! Other/optimization states

   TYPE(SD_ParameterType)                        :: p                    ! Parameters
   TYPE(SD_InputType)                            :: u(NumInp)            ! System inputs
   TYPE(SD_OutputType)                           :: y                    ! System outputs

   TYPE(SD_ContinuousStateType)                  :: dxdt                 ! First time derivatives of the continuous states

   TYPE(SD_PartialOutputPInputType)              :: dYdu                 ! Partial derivatives of the output functions
                                                                              !  (Y) with respect to the inputs (u)
   TYPE(SD_PartialContStatePInputType)           :: dXdu                 ! Partial derivatives of the continuous state
                                                                              !  functions (X) with respect to the inputs (u)
   TYPE(SD_PartialDiscStatePInputType)           :: dXddu                ! Partial derivatives of the discrete state
                                                                              !  functions (Xd) with respect to the inputs (u)
   TYPE(SD_PartialConstrStatePInputType)         :: dZdu                 ! Partial derivatives of the constraint state
                                                                              !  functions (Z) with respect to the inputs (u)
   TYPE(SD_PartialOutputPContStateType)          :: dYdx                 ! Partial derivatives of the output functions (Y)
                                                                              !  with respect to the continuous states (x)
   TYPE(SD_PartialContStatePContStateType)       :: dXdx                 ! Partial derivatives of the continuous state funct-
                                                                              !  ions (X) with respect to the continuous states (x)
   TYPE(SD_PartialDiscStatePContStateType)       :: dXddx                ! Partial derivatives of the discrete state funct-
                                                                              !  ions (Xd) with respect to continuous states (x)
   TYPE(SD_PartialConstrStatePContStateType)     :: dZdx                 ! Partial derivatives of the constraint state funct-
                                                                              !  ions (Z) with respect to the continuous states (x)
   TYPE(SD_PartialOutputPDiscStateType)          :: dYdxd                ! Partial derivatives of the output functions (Y)
                                                                              !  with respect to the discrete states (xd)
   TYPE(SD_PartialContStatePDiscStateType)       :: dXdxd                ! Partial derivatives of the continuous state funct-
                                                                              !  ions (X) with respect to the discrete states (xd)
   TYPE(SD_PartialDiscStatePDiscStateType)       :: dXddxd               ! Partial derivatives of the discrete state funct-
                                                                              !  ions (Xd) with respect to the discrete states (xd)
   TYPE(SD_PartialConstrStatePDiscStateType)     :: dZdxd                ! Partial derivatives of the constraint state funct-
                                                                              !  ions (Z) with respect to the discrete states (xd)
   TYPE(SD_PartialOutputPConstrStateType)        :: dYdz                 ! Partial derivatives of the output functions (Y)
                                                                              !  with respect to the constraint states (z)
   TYPE(SD_PartialContStatePConstrStateType)     :: dXdz                 ! Partial derivatives of the continuous state funct-
                                                                              !  ions (X) with respect to the constraint states (z)
   TYPE(SD_PartialDiscStatePConstrStateType)     :: dXddz                ! Partial derivatives of the discrete state funct-
                                                                              !  ions (Xd) with respect to constraint states (z)
   TYPE(SD_PartialConstrStatePConstrStateType)   :: dZdz                 ! Partial derivatives of the constraint state funct-
                                                                              !  ions (Z) with respect to the constraint states (z)


   INTEGER(IntKi)                                     :: n                    ! Loop counter (for time step)
   INTEGER(IntKi)                                     :: ErrStat, ErrStat1, ErrStat2, ErrStat3              ! Status of error message
   CHARACTER(1024)                                    :: ErrMsg, ErrMsg1, ErrMsg2, ErrMsg3              ! Error message if ErrStat /= ErrID_None


   REAL(ReKi), ALLOCATABLE                            :: Re_SaveAry  (:)      ! Array to store reals in packed data structure
   REAL(DbKi), ALLOCATABLE                            :: Db_SaveAry  (:)      ! Array to store doubles in packed data structure
   INTEGER(IntKi), ALLOCATABLE                        :: Int_SaveAry (:)      ! Array to store integers in packed data structure

   !...............................................................................................................................
   ! Routines called in initialization
   !...............................................................................................................................

         ! Populate the InitInData data structure here:

 ! InitInData%SDInputFile = '..\BeamFEM\IOFiles\TestBeam2.txt'
  InitInData%SDInputFile = '..\MergedSubDyn\IOFiles\TestBeam3.txt'
 ! InitInData%SDInputFile = '..\BeamFEM\IOFiles\TestFrame.txt'
   InitInData%g =  9.80665
   !InitInData%TP_RefPoint = (/0.0, 0.0, 100.0/)  !testbeam2
   InitInData%TP_RefPoint = (/50.0, 0.0, 50.0/)  !testbeam3
   !InitInData%TP_RefPoint = (/0.0, 0.0, 40.0/)  !testframe
         ! Set the driver's request for time interval here:
   TimeInterval = 0.001 ! Glue code's request for delta time (likely based on information from other modules)
   
   
         ! Initialize the module
   
   CALL SubDyn_Init( InitInData, u(1), p,  x, xd, z, OtherState, y, TimeInterval, InitOutData, ErrStat1, ErrMsg1 )
   


         ! Destroy initialization data

   CALL SD_DestroyInitInput(  InitInData,  ErrStat2, ErrMsg2 )
   CALL SD_DestroyInitOutput( InitOutData, ErrStat3, ErrMsg3 )

   
      ! Handle the initialization error after destroying the data structures
   
   IF ( ErrStat1 /= ErrID_None .OR. ErrStat2 /=0 .OR. ErrStat3 /= 0) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg1 )
      STOP
   END IF
   
   IF ( ErrStat2 /=0 .OR. ErrStat3 /= 0) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( 'Error destroying SubDyn intialization data' )
      STOP
   END IF

   !...............................................................................................................................
   ! Routines called in loose coupling -- the glue code may implement this in various ways
   !...............................................................................................................................

   ! Force the displacement of the interface node in the global Z direction to be the sag of the column under it's own weight

   !u(1)%UFL(3) = -0.001821207  !-0.001821235   !This is for testbeam.txt
   ! u(1)%UFL(3)=-12.958  !this is for testbeam3
    
   DO n = 0,600

      Time = n*TimeInterval
      InputTime(1) = Time

         ! Modify u (likely from the outputs of another module or a set of test conditions) here:

         
         ! Calculate outputs at n
      
      CALL SubDyn_CalcOutput( Time, u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
      END IF

         
         ! Get state variables at next step: INPUT at step n, OUTPUT at step n + 1
                                     
      CALL SubDyn_UpdateStates( Time, n, u, InputTime, p, x, xd, z, z_next, OtherState, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
      END IF     
      
      
   END DO


   !...............................................................................................................................
   ! Routines called in tight coupling -- time marching only
   !...............................................................................................................................

!   DO n = 0,2
!
!      Time = n * TimeInterval   ! Note that the discrete states must be updated only at the TimeInterval defined in initialization
!
!         ! set inputs (u) here:
!!      u =
!
!         ! Update constraint states at Time
!
!      ! DO 
!
!   !      CALL SD_CalcConstrStateResidual( Time, u(1), p, x, xd, z, OtherState, Z_residual, ErrStat, ErrMsg )
!   !      IF ( ErrStat /= ErrID_None ) THEN      ! Check if there was an error and do something about it if necessary
!   !         CALL WrScr( ErrMsg )
!   !      END IF
!
!         ! z =
!
!      ! END DO
!
!
!
!         ! Calculate the outputs at Time
!
!      CALL SubDyn_CalcOutput( Time, u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
!      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
!         CALL WrScr( ErrMsg )
!      END IF
!
!
!         ! Calculate the continuous state derivatives at Time
!
!   !   CALL SD_CalcContStateDeriv( Time, u(1), p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )
!   !   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
!   !      CALL WrScr( ErrMsg )
!   !   END IF
!
!
!         ! Update the discrete state from step n to step n+1
!         ! Note that the discrete states must be updated only at the TimeInterval defined in initialization
!
!!      CALL SD_UpdateDiscState( Time, n, u(1), p, x, xd, z, OtherState, ErrStat, ErrMsg )
!!      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
!!         CALL WrScr( ErrMsg )
!!      END IF
!
!
!         ! Driver should integrate (update) continuous states here:
!
!      !x = function of dxdt, x
!
!
!         ! Jacobians required:
!                              
!      CALL SD_JacobianPInput( Time, u(1), p, x, xd, z, OtherState, dYdu, dXdu, dXddu, dZdu, ErrStat, ErrMsg )
!      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
!         CALL WrScr( ErrMsg )
!      END IF
!                                   
!      CALL SD_JacobianPConstrState( Time, u(1), p, x, xd, z, OtherState, dYdz,dXdz, dXddz, dZdz, ErrStat, ErrMsg  )
!      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
!         CALL WrScr( ErrMsg )
!      END IF
!
!
!   END DO


      ! Destroy Z_residual and dxdt because they are not necessary anymore

   CALL SD_DestroyConstrState( Z_residual, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN   ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
   END IF

   CALL SD_DestroyContState( dxdt, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN   ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
   END IF

   !...............................................................................................................................
   ! Jacobian routines called in tight coupling
   !...............................................................................................................................

   CALL SD_JacobianPInput( Time, u(1), p, x, xd, z, OtherState, dYdu, dXdu, dXddu, dZdu, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
   END IF

   CALL SD_JacobianPContState( Time, u(1), p, x, xd, z, OtherState, dYdx, dXdx, dXddx, dZdx, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
   END IF

   CALL SD_JacobianPDiscState( Time, u(1), p, x, xd, z, OtherState, dYdxd, dXdxd, dXddxd, dZdxd, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
   END IF

   CALL SD_JacobianPConstrState( Time, u(1), p, x, xd, z, OtherState, dYdz, dXdz, dXddz, dZdz, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
   END IF


   !...............................................................................................................................
   ! Routines to pack data (to restart later)
   !...............................................................................................................................  
   CALL SD_Pack(Re_SaveAry, Db_SaveAry, Int_SaveAry, u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg) 
     
   IF ( ErrStat /= ErrID_None ) THEN
      CALL WrScr( ErrMsg )
   END IF


   !...............................................................................................................................
   ! Routine to terminate program execution
   !...............................................................................................................................
   CALL SubDyn_End( u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
   
   IF ( ErrStat /= ErrID_None ) THEN
      CALL WrScr( ErrMsg )
   END IF


   !...............................................................................................................................
   ! Routines to retreive packed data (unpack for restart)
   !...............................................................................................................................
   ! TODO:  BUG with Unpack and the added meshes?  GJH 6/12/13
 !  CALL SD_Unpack( Re_SaveAry, Db_SaveAry, Int_SaveAry, u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      CALL WrScr( ErrMsg )
   END IF


   !...............................................................................................................................
   ! Routines to copy data (not already tested)
   !...............................................................................................................................



   !...............................................................................................................................
   ! Routines to destroy data (not already tested)
   !...............................................................................................................................

   IF ( ALLOCATED( Re_SaveAry  ) ) DEALLOCATE( Re_SaveAry )
   IF ( ALLOCATED( Db_SaveAry  ) ) DEALLOCATE( Db_SaveAry )
   IF ( ALLOCATED( Int_SaveAry ) ) DEALLOCATE( Int_SaveAry )

!   CALL SD_DestroyPartialOutputPInput ( )  ! Jacobian Routine not yet implemented


   !...............................................................................................................................
   ! Routine to terminate program execution (again)
   !...............................................................................................................................

   CALL SubDyn_End( u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN
      CALL WrScr( ErrMsg )
   END IF


END PROGRAM TestSubDyn
