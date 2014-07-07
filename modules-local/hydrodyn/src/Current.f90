!**********************************************************************************************************************************
! The Current and Current_Types modules make up a template for creating user-defined calculations in the FAST Modularization 
! Framework. Currents_Types will be auto-generated based on a description of the variables for the module.
!
! "Current" should be replaced with the name of your module. Example: HydroDyn
! "Current" (in Current_*) should be replaced with the module name or an abbreviation of it. Example: HD
!..................................................................................................................................
! LICENSING
! Copyright (C) 2012  National Renewable Energy Laboratory
!
!    This file is part of Current.
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
! limitations under the License.
!    
!**********************************************************************************************************************************
! File last committed: $Date: 2014-06-18 12:55:01 -0600 (Wed, 18 Jun 2014) $
! (File) Revision #: $Rev: 426 $
! URL: $HeadURL: https://windsvn.nrel.gov/HydroDyn/branches/HydroDyn_Modularization/Source/Current.f90 $
!**********************************************************************************************************************************
MODULE Current

   USE Current_Types   
   USE NWTC_Library
      
   IMPLICIT NONE
   
   PRIVATE

!   INTEGER(IntKi), PARAMETER            :: DataFormatID = 1   ! Update this value if the data types change (used in Current_Pack)
   TYPE(ProgDesc), PARAMETER            :: Current_ProgDesc = ProgDesc( 'Current', '(v1.00.01, 19-October-2012)', '05-Mar-2013' )

   
      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: Current_Init                           ! Initialization routine
   PUBLIC :: Current_End                            ! Ending routine (includes clean up)
   
   PUBLIC :: Current_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating 
                                                    !   continuous states, and updating discrete states
   PUBLIC :: Current_CalcOutput                     ! Routine for computing outputs
   
   PUBLIC :: Current_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: Current_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: Current_UpdateDiscState                ! Tight coupling routine for updating discrete states
      
   !PUBLIC :: Current_JacobianPInput                 ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
   !                                                 !   (Xd), and constraint-state (Z) equations all with respect to the inputs (u)
   !PUBLIC :: Current_JacobianPContState             ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
   !                                                 !   (Xd), and constraint-state (Z) equations all with respect to the continuous 
   !                                                 !   states (x)
   !PUBLIC :: Current_JacobianPDiscState             ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
   !                                                 !   (Xd), and constraint-state (Z) equations all with respect to the discrete 
   !                                                 !   states (xd)
   !PUBLIC :: Current_JacobianPConstrState           ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
                                                    !   (Xd), and constraint-state (Z) equations all with respect to the constraint 
                                                    !   states (z)
   
   
CONTAINS

!=======================================================================
!JASON: MOVE THIS USER-DEFINED ROUTINE (UserCurrent) TO THE UserSubs.f90 OF HydroDyn WHEN THE PLATFORM LOADING FUNCTIONALITY HAS BEEN DOCUMENTED!!!!!
SUBROUTINE UserCurrent ( zi, WtrDpth, DirRoot, CurrVxi, CurrVyi )

         ! This is a dummy routine for holding the place of a user-specified
         ! current profile.  Modify this code to create your own profile.


      IMPLICIT                        NONE


         ! Passed Variables:

      REAL(ReKi), INTENT(OUT)      :: CurrVxi                                         ! xi-component of the current velocity at elevation zi, m/s.
      REAL(ReKi), INTENT(OUT)      :: CurrVyi                                         ! yi-component of the current velocity at elevation zi, m/s.
      REAL(ReKi), INTENT(IN )      :: WtrDpth                                         ! Water depth ( WtrDpth       >  0 ), meters.
      REAL(ReKi), INTENT(IN )      :: zi                                              ! Elevation   (-WtrDpth <= zi <= 0 ), meters.
      
      CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.



      CurrVxi = 0.0
      CurrVyi = 0.0



      RETURN
      
END SUBROUTINE UserCurrent
      
      
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Calc_Current( InitInp, z, h , DirRoot, CurrVxi, CurrVyi )
! This routine computes the x and y current components for a given water elevation
!----------------------------------------------------------------------------------------------------------------------------------

         ! This routine is used to initialize the variables associated with
         ! current.



      IMPLICIT                        NONE


         ! Passed Variables:

      REAL(ReKi),              INTENT(OUT) :: CurrVxi         ! xi-component of the current velocity at elevation z (m/s)
      REAL(ReKi),              INTENT(OUT) :: CurrVyi         ! yi-component of the current velocity at elevation z (m/s)

      REAL(ReKi),              INTENT(IN ) :: h               ! Water depth (meters)  This quantity must be positive-valued
      REAL(ReKi),              INTENT(IN ) :: z               ! Elevation relative to the mean sea level (meters)
      CHARACTER(1024),         INTENT(IN ) :: DirRoot         ! The name of the root file including the full path to the current working directory.  
                                                              ! This may be useful if you want this routine to write a permanent record of what it does 
                                                              ! to be stored with the simulation results: the results should be stored in a file whose name 
                                                              ! (including path) is generated by appending any suitable extension to DirRoot.
      TYPE(Current_InitInputType), INTENT(IN ) :: InitInp         ! Initialization data for the current module 


         ! Local Variables:

      REAL(ReKi)                           :: CurrSSV         ! Magnitude of sub -surface current velocity at elevation z (m/s)
      REAL(ReKi)                           :: CurrNSV         ! Magnitude of near-surface current velocity at elevation z (m/s)



         ! If elevation z lies between the seabed and the mean sea level, compute the
         !   xi- and yi-components of the current (which depends on which current
         !   profile model is selected), else set CurrVxi and CurrVyi to zero:

      IF ( ( z < -h ) .OR. ( z > 0.0 ) )  THEN  ! .TRUE. if elevation z lies below the seabed or above mean sea level (exclusive)


            CurrVxi = 0.0  ! Set both the xi- and yi-direction
            CurrVyi = 0.0  ! current velocities to zero


      ELSE                                      ! Elevation z must lie between the seabed and the mean sea level (inclusive)


         SELECT CASE ( InitInp%CurrMod ) ! Which current profile model are we using?

         CASE ( 0 )              ! None!

            CurrVxi = 0.0  ! Set both the xi- and yi-direction
            CurrVyi = 0.0  ! current velocities to zero


         CASE ( 1 )              ! Standard (using inputs from PtfmFile).

            CurrSSV =      InitInp%CurrSSV0*( ( z + h                      )/h                      )**(1.0/7.0)
            CurrNSV = MAX( InitInp%CurrNSV0*( ( z + InitInp%CurrNSRef )/InitInp%CurrNSRef )           , 0.0_ReKi )

            CurrVxi = InitInp%CurrDIV*COS( D2R*InitInp%CurrDIDir ) + CurrSSV*COS( D2R*InitInp%CurrSSDir ) + &
                                   CurrNSV*COS( D2R*InitInp%CurrNSDir )

            CurrVyi = InitInp%CurrDIV*SIN( D2R*InitInp%CurrDIDir ) + CurrSSV*SIN( D2R*InitInp%CurrSSDir ) + &
                                   CurrNSV*SIN( D2R*InitInp%CurrNSDir )


         CASE ( 2 )              ! User-defined current profile model.

            CALL UserCurrent ( z, h, DirRoot, CurrVxi, CurrVyi )

         ENDSELECT

      END IF

   RETURN
   
END SUBROUTINE Calc_Current


!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Current_Init( InitInp, u, p, x, xd, z, OtherState, y, Interval, InitOut, ErrStat, ErrMsg )
! This routine is called at the start of the simulation to perform initialization steps. 
! The parameters are set here and not changed during the simulation.
! The initial states and initial guess for the input are defined.
!..................................................................................................................................

   TYPE(Current_InitInputType),       INTENT(IN   )  :: InitInp     ! Input data for initialization routine
   TYPE(Current_InputType),           INTENT(  OUT)  :: u           ! An initial guess for the input; input mesh must be defined
   TYPE(Current_ParameterType),       INTENT(  OUT)  :: p           ! Parameters      
   TYPE(Current_ContinuousStateType), INTENT(  OUT)  :: x           ! Initial continuous states
   TYPE(Current_DiscreteStateType),   INTENT(  OUT)  :: xd          ! Initial discrete states
   TYPE(Current_ConstraintStateType), INTENT(  OUT)  :: z           ! Initial guess of the constraint states
   TYPE(Current_OtherStateType),      INTENT(  OUT)  :: OtherState  ! Initial other/optimization states            
   TYPE(Current_OutputType),          INTENT(  OUT)  :: y           ! Initial system outputs (outputs are not calculated; 
                                                                    !   only the output mesh is initialized)
   REAL(DbKi),                        INTENT(INOUT)  :: Interval    ! Coupling interval in seconds: the rate that 
                                                                    !   (1) Current_UpdateStates() is called in loose coupling &
                                                                    !   (2) Current_UpdateDiscState() is called in tight coupling.
                                                                    !   Input is the suggested time from the glue code; 
                                                                   !   Output is the actual coupling interval that will be used 
                                                                     !   by the glue code.
   TYPE(Current_InitOutputType),      INTENT(  OUT)  :: InitOut     ! Output for initialization routine
   INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


     
      

      ! Local Variables:

   REAL(ReKi)                   :: CurrVxi                          ! xi-component of the current velocity at elevation z (m/s)     
   REAL(ReKi)                   :: CurrVyi                          ! yi-component of the current velocity at elevation z (m/s)
   REAL(ReKi)                   :: CurrVxi0                         ! xi-component of the current velocity at zi =  0.0 meters            (m/s)
   REAL(ReKi)                   :: CurrVyi0                         ! yi-component of the current velocity at zi =  0.0 meters            (m/s)
   REAL(ReKi)                   :: CurrVxiS                         ! xi-component of the current velocity at zi = -SmllNmbr meters       (m/s)
   REAL(ReKi)                   :: CurrVyiS                         ! yi-component of the current velocity at zi = -SmllNmbr meters       (m/s)
   REAL(ReKi), PARAMETER        :: SmllNmbr  = 9.999E-4             ! A small number representing epsilon for taking numerical derivatives.
   
   INTEGER                      :: I                                ! Generic index
      
      
      ! Initialize ErrStat
         
   ErrStat = ErrID_None         
   ErrMsg  = ""               
      
      
      ! Initialize the NWTC Subroutine Library
         
   CALL NWTC_Init(  )

  
   
      ! IF there are Morison elements, then compute the current components at each morison node elevation
      
   IF ( InitInp%NMorisonNodes > 0 ) THEN    
         
      ALLOCATE ( InitOut%CurrVxi( InitInp%NMorisonNodes ) , STAT=ErrStat )
      IF ( ErrStat /= ErrID_None )  THEN
         ErrMsg = ' Error allocating memory for the CurrVxi array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      
      ALLOCATE ( InitOut%CurrVyi( InitInp%NMorisonNodes ) , STAT=ErrStat )
      IF ( ErrStat /= ErrID_None )  THEN
         ErrMsg = ' Error allocating memory for the CurrVyi array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      
      
         ! Loop over all of the points where current information is required
      
      DO I = 1, InitInp%NMorisonNodes
         
         CALL Calc_Current( InitInp, InitInp%MorisonNodezi(I), InitInp%WtrDpth, InitInp%DirRoot, CurrVxi, CurrVyi )        
         
         InitOut%CurrVxi(I) = CurrVxi
         InitOut%CurrVyi(I) = CurrVyi
       
      END DO     
     
   END IF   
      

      ! Compute the partial derivative for wave stretching
   CALL    Calc_Current( InitInp,  0.0_ReKi, InitInp%WtrDpth, InitInp%DirRoot, CurrVxi0, CurrVyi0 )
   CALL    Calc_Current( InitInp, -SmllNmbr, InitInp%WtrDpth, InitInp%DirRoot, CurrVxiS, CurrVyiS )

   InitOut%PCurrVxiPz0 = ( CurrVxi0 - CurrVxiS )/SmllNmbr                    ! xi-direction
   InitOut%PCurrVyiPz0 = ( CurrVyi0 - CurrVyiS )/SmllNmbr                    ! yi-direction
   
   
   u%DummyInput = 0.0
   p%DT = Interval
   x%DummyContState = 0.0
   xd%DummyDiscState = 0.0
   z%DummyConstrState = 0.0
   OtherState%DummyOtherState = 0.0
   y%DummyOutput = 0.0
   
END SUBROUTINE Current_Init


!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Current_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
! This routine is called at the end of the simulation.
!..................................................................................................................................

      TYPE(Current_InputType),           INTENT(INOUT)  :: u           ! System inputs
      TYPE(Current_ParameterType),       INTENT(INOUT)  :: p           ! Parameters     
      TYPE(Current_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states
      TYPE(Current_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Discrete states
      TYPE(Current_ConstraintStateType), INTENT(INOUT)  :: z           ! Constraint states
      TYPE(Current_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states            
      TYPE(Current_OutputType),          INTENT(INOUT)  :: y           ! System outputs
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Place any last minute operations or calculations here:


         ! Close files here:     
                  
                  

         ! Destroy the input data:
         
      CALL Current_DestroyInput( u, ErrStat, ErrMsg )


         ! Destroy the parameter data:
         
      CALL Current_DestroyParam( p, ErrStat, ErrMsg )


         ! Destroy the state data:
         
      CALL Current_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL Current_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL Current_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL Current_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )
         

         ! Destroy the output data:
         
      CALL Current_DestroyOutput( y, ErrStat, ErrMsg )


      

END SUBROUTINE Current_End
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Current_UpdateStates( Time, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
! Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
! Constraint states are solved for input Time; Continuous and discrete states are updated for Time + Interval
!..................................................................................................................................
   
      REAL(DbKi),                         INTENT(IN   ) :: Time        ! Current simulation time in seconds
      TYPE(Current_InputType),            INTENT(IN   ) :: u           ! Inputs at Time                    
      TYPE(Current_ParameterType),        INTENT(IN   ) :: p           ! Parameters                              
      TYPE(Current_ContinuousStateType),  INTENT(INOUT) :: x           ! Input: Continuous states at Time; 
                                                                       !   Output: Continuous states at Time + Interval
      TYPE(Current_DiscreteStateType),    INTENT(INOUT) :: xd          ! Input: Discrete states at Time; 
                                                                       !   Output: Discrete states at Time  + Interval
      TYPE(Current_ConstraintStateType),  INTENT(INOUT) :: z           ! Input: Initial guess of constraint states at Time;
                                                                       !   Output: Constraint states at Time
      TYPE(Current_OtherStateType),       INTENT(INOUT) :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                     INTENT(  OUT) :: ErrStat     ! Error status of the operation     
      CHARACTER(*),                       INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None

         ! Local variables
         
      TYPE(Current_ContinuousStateType)                 :: dxdt        ! Continuous state derivatives at Time
      TYPE(Current_ConstraintStateType)                 :: z_Residual  ! Residual of the constraint state equations (Z)
         
      INTEGER(IntKi)                                    :: ErrStat2    ! Error status of the operation (occurs after initial error)
      CHARACTER(LEN(ErrMsg))                            :: ErrMsg2     ! Error message if ErrStat2 /= ErrID_None
                        
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
           
      
         ! Solve for the constraint states (z) here:
                           
         ! Check if the z guess is correct and update z with a new guess.
         ! Iterate until the value is within a given tolerance. 
                                    
      CALL Current_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_Residual, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN      
         CALL Current_DestroyConstrState( z_Residual, ErrStat2, ErrMsg2)
         ErrMsg = TRIM(ErrMsg)//' '//TRIM(ErrMsg2)
         RETURN      
      END IF
         
      ! DO WHILE ( z_Residual% > tolerance )
      !
      !  z = 
      !
      !  CALL Current_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_Residual, ErrStat, ErrMsg )
      !  IF ( ErrStat >= AbortErrLev ) THEN      
      !     CALL Current_DestroyConstrState( z_Residual, ErrStat2, ErrMsg2)
      !     ErrMsg = TRIM(ErrMsg)//' '//TRIM(ErrMsg2)
      !     RETURN      
      !  END IF
      !           
      ! END DO         
      
      
         ! Destroy z_Residual because it is not necessary for the rest of the subroutine:
            
      CALL Current_DestroyConstrState( z_Residual, ErrStat, ErrMsg)
      IF ( ErrStat >= AbortErrLev ) RETURN      
         
         
         
         ! Get first time derivatives of continuous states (dxdt):
      
      CALL Current_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN      
         CALL Current_DestroyContState( dxdt, ErrStat2, ErrMsg2)
         ErrMsg = TRIM(ErrMsg)//' '//TRIM(ErrMsg2)
         RETURN
      END IF
               
               
         ! Update discrete states:
         !   Note that xd [discrete state] is changed in Current_UpdateDiscState(), so Current_CalcOutput(),  
         !   Current_CalcContStateDeriv(), and Current_CalcConstrStates() must be called first (see above).
      
      CALL Current_UpdateDiscState(Time, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )   
      IF ( ErrStat >= AbortErrLev ) THEN      
         CALL Current_DestroyContState( dxdt, ErrStat2, ErrMsg2)
         ErrMsg = TRIM(ErrMsg)//' '//TRIM(ErrMsg2)
         RETURN      
      END IF
         
         
         ! Integrate (update) continuous states (x) here:
         
      !x = function of dxdt and x


         ! Destroy dxdt because it is not necessary for the rest of the subroutine
            
      CALL Current_DestroyContState( dxdt, ErrStat, ErrMsg)
      IF ( ErrStat >= AbortErrLev ) RETURN      
     
   
      
END SUBROUTINE Current_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Current_CalcOutput( Time, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )   
! Routine for computing outputs, used in both loose and tight coupling.
!..................................................................................................................................
   
      REAL(DbKi),                        INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(Current_InputType),           INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(Current_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(Current_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(Current_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(Current_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(Current_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(Current_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at Time (Input only so that mesh con-
                                                                       !   nectivity information does not have to be recalculated)
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Compute outputs here:
      y%DummyOutput    = 2.0_ReKi

   
               

END SUBROUTINE Current_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Current_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )  
! Tight coupling routine for computing derivatives of continuous states
!..................................................................................................................................
   
      REAL(DbKi),                        INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(Current_InputType),           INTENT(IN   )  :: u           ! Inputs at Time                    
      TYPE(Current_ParameterType),       INTENT(IN   )  :: p           ! Parameters                             
      TYPE(Current_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(Current_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(Current_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(Current_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states                    
      TYPE(Current_ContinuousStateType), INTENT(  OUT)  :: dxdt        ! Continuous state derivatives at Time
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation     
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

               
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Compute the first time derivatives of the continuous states here:
      
      dxdt%DummyContState = 0.0
         

END SUBROUTINE Current_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Current_UpdateDiscState( Time, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )   
! Tight coupling routine for updating discrete states
!..................................................................................................................................
   
      REAL(DbKi),                        INTENT(IN   )  :: Time        ! Current simulation time in seconds   
      TYPE(Current_InputType),           INTENT(IN   )  :: u           ! Inputs at Time                       
      TYPE(Current_ParameterType),       INTENT(IN   )  :: p           ! Parameters                                 
      TYPE(Current_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(Current_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Input: Discrete states at Time; 
                                                                       !   Output: Discrete states at Time + Interval
      TYPE(Current_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(Current_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states           
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

               
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Update discrete states here:
      
      ! StateData%DiscState = 

END SUBROUTINE Current_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Current_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_residual, ErrStat, ErrMsg )   
! Tight coupling routine for solving for the residual of the constraint state equations
!..................................................................................................................................
   
      REAL(DbKi),                        INTENT(IN   )  :: Time        ! Current simulation time in seconds   
      TYPE(Current_InputType),           INTENT(IN   )  :: u           ! Inputs at Time                       
      TYPE(Current_ParameterType),       INTENT(IN   )  :: p           ! Parameters                           
      TYPE(Current_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(Current_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(Current_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time (possibly a guess)
      TYPE(Current_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states                    
      TYPE(Current_ConstraintStateType), INTENT(  OUT)  :: z_residual  ! Residual of the constraint state equations using  
                                                                       !     the input values described above      
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

               
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Solve for the constraint states here:
      
      z_residual%DummyConstrState = 0

END SUBROUTINE Current_CalcConstrStateResidual
!!----------------------------------------------------------------------------------------------------------------------------------
!SUBROUTINE Current_JacobianPInput( Time, u, p, x, xd, z, OtherState, dYdu, dXdu, dXddu, dZdu, ErrStat, ErrMsg )   
!! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) equations 
!! with respect to the inputs (u). The partial derivatives dY/du, dX/du, dXd/du, and DZ/du are returned.
!!..................................................................................................................................
!   
!      REAL(DbKi),                                INTENT(IN   )           :: Time       ! Current simulation time in seconds   
!      TYPE(Current_InputType),                   INTENT(IN   )           :: u          ! Inputs at Time                       
!      TYPE(Current_ParameterType),               INTENT(IN   )           :: p          ! Parameters                           
!      TYPE(Current_ContinuousStateType),         INTENT(IN   )           :: x          ! Continuous states at Time
!      TYPE(Current_DiscreteStateType),           INTENT(IN   )           :: xd         ! Discrete states at Time
!      TYPE(Current_ConstraintStateType),         INTENT(IN   )           :: z          ! Constraint states at Time
!      TYPE(Current_OtherStateType),              INTENT(INOUT)           :: OtherState ! Other/optimization states                    
!      TYPE(Current_PartialOutputPInputType),     INTENT(  OUT), OPTIONAL :: dYdu       ! Partial derivatives of output equations
!                                                                                       !   (Y) with respect to the inputs (u)
!      TYPE(Current_PartialContStatePInputType),  INTENT(  OUT), OPTIONAL :: dXdu       ! Partial derivatives of continuous state
!                                                                                       !   equations (X) with respect to inputs (u)
!      TYPE(Current_PartialDiscStatePInputType),  INTENT(  OUT), OPTIONAL :: dXddu      ! Partial derivatives of discrete state 
!                                                                                       !   equations (Xd) with respect to inputs (u)
!      TYPE(Current_PartialConstrStatePInputType),INTENT(  OUT), OPTIONAL :: dZdu       ! Partial derivatives of constraint state 
!                                                                                       !   equations (Z) with respect to inputs (u)
!      INTEGER(IntKi),                            INTENT(  OUT)           :: ErrStat    ! Error status of the operation
!      CHARACTER(*),                              INTENT(  OUT)           :: ErrMsg     ! Error message if ErrStat /= ErrID_None
!
!               
!         ! Initialize ErrStat
!         
!      ErrStat = ErrID_None         
!      ErrMsg  = ""               
!      
!      
!      IF ( PRESENT( dYdu ) ) THEN
!      
!         ! Calculate the partial derivative of the output equations (Y) with respect to the inputs (u) here:
!
!         dYdu%DummyOutput%DummyInput = 0
!
!      END IF
!      
!      IF ( PRESENT( dXdu ) ) THEN
!      
!         ! Calculate the partial derivative of the continuous state equations (X) with respect to the inputs (u) here:
!      
!         dXdu%DummyContState%DummyInput = 0
!
!      END IF
!      
!      IF ( PRESENT( dXddu ) ) THEN
!
!         ! Calculate the partial derivative of the discrete state equations (Xd) with respect to the inputs (u) here:
!
!         dXddu%DummyDiscState%DummyInput = 0
!
!      END IF
!      
!      IF ( PRESENT( dZdu ) ) THEN
!
!         ! Calculate the partial derivative of the constraint state equations (Z) with respect to the inputs (u) here:
!      
!         dZdu%DummyConstrState%DummyInput = 0
!
!      END IF
!
!
!END SUBROUTINE Current_JacobianPInput
!!----------------------------------------------------------------------------------------------------------------------------------
!SUBROUTINE Current_JacobianPContState( Time, u, p, x, xd, z, OtherState, dYdx, dXdx, dXddx, dZdx, ErrStat, ErrMsg )   
!! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) equations
!! with respect to the continuous states (x). The partial derivatives dY/dx, dX/dx, dXd/dx, and DZ/dx are returned.
!!..................................................................................................................................
!   
!      REAL(DbKi),                                    INTENT(IN   )           :: Time       ! Current simulation time in seconds   
!      TYPE(Current_InputType),                       INTENT(IN   )           :: u          ! Inputs at Time                       
!      TYPE(Current_ParameterType),                   INTENT(IN   )           :: p          ! Parameters                           
!      TYPE(Current_ContinuousStateType),             INTENT(IN   )           :: x          ! Continuous states at Time
!      TYPE(Current_DiscreteStateType),               INTENT(IN   )           :: xd         ! Discrete states at Time
!      TYPE(Current_ConstraintStateType),             INTENT(IN   )           :: z          ! Constraint states at Time
!      TYPE(Current_OtherStateType),                  INTENT(INOUT)           :: OtherState ! Other/optimization states                    
!      TYPE(Current_PartialOutputPContStateType),     INTENT(  OUT), OPTIONAL :: dYdx       ! Partial derivatives of output equations
!                                                                                           !   (Y) with respect to the continuous 
!                                                                                           !   states (x)
!      TYPE(Current_PartialContStatePContStateType),  INTENT(  OUT), OPTIONAL :: dXdx       ! Partial derivatives of continuous state
!                                                                                           !   equations (X) with respect to 
!                                                                                           !   the continuous states (x)
!      TYPE(Current_PartialDiscStatePContStateType),  INTENT(  OUT), OPTIONAL :: dXddx      ! Partial derivatives of discrete state 
!                                                                                           !   equations (Xd) with respect to 
!                                                                                           !   the continuous states (x)
!      TYPE(Current_PartialConstrStatePContStateType),INTENT(  OUT), OPTIONAL :: dZdx       ! Partial derivatives of constraint state
!                                                                                           !   equations (Z) with respect to 
!                                                                                           !   the continuous states (x)
!      INTEGER(IntKi),                                INTENT(  OUT)           :: ErrStat    ! Error status of the operation
!      CHARACTER(*),                                  INTENT(  OUT)           :: ErrMsg     ! Error message if ErrStat /= ErrID_None
!
!               
!         ! Initialize ErrStat
!         
!      ErrStat = ErrID_None         
!      ErrMsg  = ""               
!      
!      
!     
!      IF ( PRESENT( dYdx ) ) THEN
!
!         ! Calculate the partial derivative of the output equations (Y) with respect to the continuous states (x) here:
!
!         dYdx%DummyOutput%DummyContState = 0
!
!      END IF
!      
!      IF ( PRESENT( dXdx ) ) THEN
!      
!         ! Calculate the partial derivative of the continuous state equations (X) with respect to the continuous states (x) here:
!      
!         dXdx%DummyContState%DummyContState = 0
!
!      END IF
!      
!      IF ( PRESENT( dXddx ) ) THEN
!
!         ! Calculate the partial derivative of the discrete state equations (Xd) with respect to the continuous states (x) here:
!
!         dXddx%DummyDiscState%DummyContState = 0
!         
!      END IF
!      
!      IF ( PRESENT( dZdx ) ) THEN
!
!
!         ! Calculate the partial derivative of the constraint state equations (Z) with respect to the continuous states (x) here:
!      
!         dZdx%DummyConstrState%DummyContState = 0
!
!      END IF
!      
!
!   END SUBROUTINE Current_JacobianPContState
!!----------------------------------------------------------------------------------------------------------------------------------
!SUBROUTINE Current_JacobianPDiscState( Time, u, p, x, xd, z, OtherState, dYdxd, dXdxd, dXddxd, dZdxd, ErrStat, ErrMsg )   
!! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) equations
!! with respect to the discrete states (xd). The partial derivatives dY/dxd, dX/dxd, dXd/dxd, and DZ/dxd are returned.
!!..................................................................................................................................
!
!      REAL(DbKi),                                    INTENT(IN   )           :: Time       ! Current simulation time in seconds   
!      TYPE(Current_InputType),                       INTENT(IN   )           :: u          ! Inputs at Time                       
!      TYPE(Current_ParameterType),                   INTENT(IN   )           :: p          ! Parameters                           
!      TYPE(Current_ContinuousStateType),             INTENT(IN   )           :: x          ! Continuous states at Time
!      TYPE(Current_DiscreteStateType),               INTENT(IN   )           :: xd         ! Discrete states at Time
!      TYPE(Current_ConstraintStateType),             INTENT(IN   )           :: z          ! Constraint states at Time
!      TYPE(Current_OtherStateType),                  INTENT(INOUT)           :: OtherState ! Other/optimization states                    
!      TYPE(Current_PartialOutputPDiscStateType),     INTENT(  OUT), OPTIONAL :: dYdxd      ! Partial derivatives of output equations
!                                                                                           !  (Y) with respect to the discrete 
!                                                                                           !  states (xd)
!      TYPE(Current_PartialContStatePDiscStateType),  INTENT(  OUT), OPTIONAL :: dXdxd      ! Partial derivatives of continuous state
!                                                                                           !   equations (X) with respect to the 
!                                                                                           !   discrete states (xd)
!      TYPE(Current_PartialDiscStatePDiscStateType),  INTENT(  OUT), OPTIONAL :: dXddxd     ! Partial derivatives of discrete state 
!                                                                                           !   equations (Xd) with respect to the
!                                                                                           !   discrete states (xd)
!      TYPE(Current_PartialConstrStatePDiscStateType),INTENT(  OUT), OPTIONAL :: dZdxd      ! Partial derivatives of constraint state
!                                                                                           !   equations (Z) with respect to the 
!                                                                                           !   discrete states (xd)
!      INTEGER(IntKi),                                INTENT(  OUT)           :: ErrStat    ! Error status of the operation
!      CHARACTER(*),                                  INTENT(  OUT)           :: ErrMsg     ! Error message if ErrStat /= ErrID_None
!
!               
!         ! Initialize ErrStat
!         
!      ErrStat = ErrID_None         
!      ErrMsg  = ""               
!      
!      
!      IF ( PRESENT( dYdxd ) ) THEN
!      
!         ! Calculate the partial derivative of the output equations (Y) with respect to the discrete states (xd) here:
!
!         dYdxd%DummyOutput%DummyDiscState = 0
!
!      END IF
!      
!      IF ( PRESENT( dXdxd ) ) THEN
!
!         ! Calculate the partial derivative of the continuous state equations (X) with respect to the discrete states (xd) here:
!      
!         dXdxd%DummyContState%DummyDiscState = 0
!
!      END IF
!      
!      IF ( PRESENT( dXddxd ) ) THEN
!
!         ! Calculate the partial derivative of the discrete state equations (Xd) with respect to the discrete states (xd) here:
!
!         dXddxd%DummyDiscState%DummyDiscState = 0
!
!      END IF
!      
!      IF ( PRESENT( dZdxd ) ) THEN
!
!         ! Calculate the partial derivative of the constraint state equations (Z) with respect to the discrete states (xd) here:
!      
!         dZdxd%DummyConstrState%DummyDiscState = 0
!
!      END IF
!      
!
!
!END SUBROUTINE Current_JacobianPDiscState
!!----------------------------------------------------------------------------------------------------------------------------------    
!SUBROUTINE Current_JacobianPConstrState( Time, u, p, x, xd, z, OtherState, dYdz, dXdz, dXddz, dZdz, ErrStat, ErrMsg )   
!! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) equations
!! with respect to the constraint states (z). The partial derivatives dY/dz, dX/dz, dXd/dz, and DZ/dz are returned.
!!..................................................................................................................................
!   
!      REAL(DbKi),                                      INTENT(IN   )           :: Time       ! Current simulation time in seconds   
!      TYPE(Current_InputType),                         INTENT(IN   )           :: u          ! Inputs at Time                       
!      TYPE(Current_ParameterType),                     INTENT(IN   )           :: p          ! Parameters                           
!      TYPE(Current_ContinuousStateType),               INTENT(IN   )           :: x          ! Continuous states at Time
!      TYPE(Current_DiscreteStateType),                 INTENT(IN   )           :: xd         ! Discrete states at Time
!      TYPE(Current_ConstraintStateType),               INTENT(IN   )           :: z          ! Constraint states at Time
!      TYPE(Current_OtherStateType),                    INTENT(INOUT)           :: OtherState ! Other/optimization states                    
!      TYPE(Current_PartialOutputPConstrStateType),     INTENT(  OUT), OPTIONAL :: dYdz       ! Partial derivatives of output 
!                                                                                             !  equations (Y) with respect to the 
!                                                                                             !  constraint states (z)
!      TYPE(Current_PartialContStatePConstrStateType),  INTENT(  OUT), OPTIONAL :: dXdz       ! Partial derivatives of continuous
!                                                                                             !  state equations (X) with respect to 
!                                                                                             !  the constraint states (z)
!      TYPE(Current_PartialDiscStatePConstrStateType),  INTENT(  OUT), OPTIONAL :: dXddz      ! Partial derivatives of discrete state
!                                                                                             !  equations (Xd) with respect to the 
!                                                                                             !  constraint states (z)
!      TYPE(Current_PartialConstrStatePConstrStateType),INTENT(  OUT), OPTIONAL :: dZdz       ! Partial derivatives of constraint 
!                                                                                             ! state equations (Z) with respect to 
!                                                                                             !  the constraint states (z)
!      INTEGER(IntKi),                                  INTENT(  OUT)           :: ErrStat    ! Error status of the operation
!      CHARACTER(*),                                    INTENT(  OUT)           :: ErrMsg     ! Error message if ErrStat /= ErrID_None
!
!               
!         ! Initialize ErrStat
!         
!      ErrStat = ErrID_None         
!      ErrMsg  = ""               
!      
!      IF ( PRESENT( dYdz ) ) THEN
!      
!            ! Calculate the partial derivative of the output equations (Y) with respect to the constraint states (z) here:
!        
!         dYdz%DummyOutput%DummyConstrState = 0
!         
!      END IF
!      
!      IF ( PRESENT( dXdz ) ) THEN
!      
!            ! Calculate the partial derivative of the continuous state equations (X) with respect to the constraint states (z) here:
!         
!         dXdz%DummyContState%DummyConstrState = 0
!
!      END IF
!      
!      IF ( PRESENT( dXddz ) ) THEN
!
!            ! Calculate the partial derivative of the discrete state equations (Xd) with respect to the constraint states (z) here:
!
!         dXddz%DummyDiscState%DummyConstrState = 0
!
!      END IF
!      
!      IF ( PRESENT( dZdz ) ) THEN
!
!            ! Calculate the partial derivative of the constraint state equations (Z) with respect to the constraint states (z) here:
!         
!         dZdz%DummyConstrState%DummyConstrState = 0
!
!      END IF
!      
!
!END SUBROUTINE Current_JacobianPConstrState

!----------------------------------------------------------------------------------------------------------------------------------
   
END MODULE Current
!**********************************************************************************************************************************
