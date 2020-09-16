!**********************************************************************************************************************************
! The SS_Excitation and SS_Excitation_Types modules make up a template for creating user-defined calculations in the FAST Modularization 
! Framework. SS_Excitations_Types will be auto-generated based on a description of the variables for the module.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2012, 2018  National Renewable Energy Laboratory
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
MODULE SS_Excitation

   USE SS_Excitation_Types   
   USE NWTC_Library
      
   IMPLICIT NONE
   
   PRIVATE

   TYPE(ProgDesc), PARAMETER  :: SS_Exc_ProgDesc = ProgDesc( 'SS_Excitation', '', '' )

   
      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: SS_Exc_Init                           ! Initialization routine
   PUBLIC :: SS_Exc_End                            ! Ending routine (includes clean up)
   
   PUBLIC :: SS_Exc_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating 
                                                   !   continuous states, and updating discrete states
   PUBLIC :: SS_Exc_CalcOutput                     ! Routine for computing outputs
   
   PUBLIC :: SS_Exc_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: SS_Exc_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: SS_Exc_UpdateDiscState                ! Tight coupling routine for updating discrete states
         
   
CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the start of the simulation to perform initialization steps. 
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
SUBROUTINE SS_Exc_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )
!..................................................................................................................................

    TYPE(SS_Exc_InitInputType),       INTENT(IN   )  :: InitInp     !< Input data for initialization routine
    TYPE(SS_Exc_InputType),           INTENT(  OUT)  :: u           !< An initial guess for the input; input mesh must be defined
    TYPE(SS_Exc_ParameterType),       INTENT(  OUT)  :: p           !< Parameters      
    TYPE(SS_Exc_ContinuousStateType), INTENT(  OUT)  :: x           !< Initial continuous states
    TYPE(SS_Exc_DiscreteStateType),   INTENT(  OUT)  :: xd          !< Initial discrete states
    TYPE(SS_Exc_ConstraintStateType), INTENT(  OUT)  :: z           !< Initial guess of the constraint states
    TYPE(SS_Exc_OtherStateType),      INTENT(  OUT)  :: OtherState  !< Initial other states            
    TYPE(SS_Exc_OutputType),          INTENT(  OUT)  :: y           !< Initial system outputs (outputs are not calculated; 
                                                                    !!   only the output mesh is initialized)
    TYPE(SS_Exc_MiscVarType),         INTENT(  OUT)  :: m           !< Initial misc/optimization variables            
    REAL(DbKi),                       INTENT(INOUT)  :: Interval    !< Coupling interval in seconds: the rate that 
                                                                    !!   (1) SS_Exc_UpdateStates() is called in loose coupling &
                                                                    !!   (2) SS_Exc_UpdateDiscState() is called in tight coupling.
                                                                    !!   Input is the suggested time from the glue code; 
                                                                    !!   Output is the actual coupling interval that will be used 
                                                                    !!   by the glue code.
    TYPE(SS_Exc_InitOutputType),      INTENT(  OUT)  :: InitOut     !< Output for initialization routine
    INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     !< Error status of the operation
    CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

    ! Local Variables:
         
    INTEGER                                :: I                                    ! Generic index
!    INTEGER                                :: J                                    ! Generic index  
    INTEGER                                :: xx (1,6)                             ! Active DOF's on the input file .ss
    INTEGER                                :: Nlines                               ! Number of lines in the input file, used to determine N
    INTEGER                                :: UnSS                                 ! I/O unit number for the WAMIT output file with the .ss extension; this file contains the state-space matrices.
    INTEGER                                :: Sttus                                ! Error in reading .ssexctn file
    !CHARACTER                              :: Line                                 ! Temp line of file
    real(ReKi)                             :: WaveDir                              ! Temp wave direction angle (deg)
    integer                                :: ErrStat2
    character(ErrMsgLen)                   :: ErrMsg2
    
    ! Initialize ErrStat   
    ErrStat = ErrID_None
    ErrMsg  = ""
      
    u%DummyInput = 0.0_ReKi
      
    UnSS  = -1
    p%N     =  0
      
    ! Open the .ss input file!
    CALL GetNewUnit( UnSS )
    CALL OpenFInpFile ( UnSS, TRIM(InitInp%InputFile)//'.ssexctn', ErrStat2, ErrMsg2 )  ! Open file.
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Exc_Init')
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

    ! Determine the number of states and size of the matrices
    Nlines = 1
    
    CALL ReadCom ( UnSS, TRIM(InitInp%InputFile)//'.ssexctn', 'Header',ErrStat2, ErrMsg2  )! Reads the first entire line (Title header)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Exc_Init')    
   
    CALL ReadVar( UnSS,TRIM(InitInp%InputFile)//'.ssexctn', WaveDir, 'WaveDir', 'Wave direction (deg)',ErrStat2, ErrMsg2) ! Reads in the second line, containing the wave direction
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Exc_Init')
         
   ! Check that excitation state-space file Beta angle (in degrees) matches the HydroDyn input file angle
   if ( .not. EqualRealNos(InitInp%WaveDir, WaveDir) ) call SetErrStat(ErrID_FATAL,'HydroDyn Wave direction does not match the wave excitation wave direction',ErrStat,ErrMsg,'SS_Exc_Init')

   CALL ReadVar( UnSS,TRIM(InitInp%InputFile)//'.ssexctn', p%Tc, 'p%Tc', 'Time offset (s)',ErrStat2, ErrMsg2) ! Reads in the third line, containing the number of states
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Exc_Init')
            
    CALL ReadVar( UnSS,TRIM(InitInp%InputFile)//'.ssexctn', p%N, 'p%N', 'Number of states',ErrStat2, ErrMsg2) ! Reads in the third line, containing the number of states
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Exc_Init')
            
   CALL ReadAry( UnSS,TRIM(InitInp%InputFile)//'.ssexctn', p%spDOF, 6, 'p%spDOF', 'States per DOF',ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Exc_Init')
          
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF
      
    DO !Loop through all the lines of the file
        CALL ReadCom ( UnSS, TRIM(InitInp%InputFile)//'.ssexctn', 'Header',Sttus,ErrMsg2  )! Reads the first entire line (Title header)
        IF ( Sttus == ErrID_None )  THEN ! .TRUE. when data is read in successfully                    
            Nlines=Nlines+1                    
        ELSE !We must have reached the end of the file
            EXIT
        END IF
    END DO

    ! The input file contains the matrices A [NxN], B [Nx1] and C [6xN], so
    !p%N = ( Nlines - 1 ) / 2 ! this is the number of states
    
    !Verifications on the input file
    IF ( ( Nlines - 6 ) / 2 /= p%N) THEN
      CALL SetErrStat(ErrID_Severe,'Error in the input file .ssexctn: The size of the matrices does not correspond to the number of states!',ErrStat,ErrMsg,'SS_Exc_Init')
    END IF
        
    
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF
    
    ! Now we can allocate the temporary matrices A, B and C
    
    CALL AllocAry( p%A, p%N,    p%N,    'p%A', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Exc_Init')
    CALL AllocAry( p%B, p%N,            'p%B', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Exc_Init')
    CALL AllocAry( p%C,   6,    p%N,    'p%C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Exc_Init')
    
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF
    
        
    REWIND (UNIT=UnSS)   ! REWIND the file so we can read it in a second time.

    ! Skip the first 4 lines:  (NOTE: no error handling here because we would have caught it the first time through)
    CALL ReadCom ( UnSS, TRIM(InitInp%InputFile)//'.ssexctn', 'Header', ErrStat2, ErrMsg2  )! Reads the first entire line (Title header)
    CALL ReadCom ( UnSS, TRIM(InitInp%InputFile)//'.ssexctn', 'Wave direction (deg)', ErrStat2, ErrMsg2  )! Reads the first entire line (Title header)
    CALL ReadCom ( UnSS, TRIM(InitInp%InputFile)//'.ssexctn', 'Time offset (s)', ErrStat2, ErrMsg2  )! Reads the first entire line (Title header)
    CALL ReadCom ( UnSS, TRIM(InitInp%InputFile)//'.ssexctn', 'Number of Excitation States', ErrStat2, ErrMsg2  )! Reads the first entire line (Title header)
    CALL ReadCom ( UnSS, TRIM(InitInp%InputFile)//'.ssexctn', 'Number of states per dofs', ErrStat2, ErrMsg2  )! Reads the first entire line (Title header)   
    
    DO I = 1,p%N !Read A MatriX
        CALL ReadAry( UnSS,TRIM(InitInp%InputFile)//'.ssexctn', p%A(I,:), p%N, 'p%A', 'A_Matrix',ErrStat2, ErrMsg2)
          CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Exc_Init')
    END DO
    
    DO I = 1,p%N !Read B Matrix
        CALL ReadVar( UnSS, TRIM(InitInp%InputFile)//'.ssexctn', p%B(I), 'p%B', 'B_Matrix',ErrStat2, ErrMsg2) 
          CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Exc_Init')
    END DO
    
   DO I = 1,6 !Read C Matrix
      CALL ReadAry( UnSS, TRIM(InitInp%InputFile)//'.ssexctn', p%C(I,:), p%N, 'p%C', 'C_Matrix',ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Exc_Init')
   END DO 
    CLOSE ( UnSS ) !Close .ss input file
    UnSS = -1        ! Indicate the file is closed
    
    
    CALL WrScr1 ( 'Using SS_Excitation Module, with '//TRIM( Num2LStr(p%N ))//' excitation states' )
  
    ! Define parameters here:
         
      p%DT  = Interval
      
         ! Allocate Wave-elevation related arrays
      p%NStepWave = InitInp%NStepWave
      allocate ( p%WaveElev0(0:p%NStepWave) , STAT=ErrStat2 )
      allocate ( p%WaveTime (0:p%NStepWave) , STAT=ErrStat2 )
!TODO: Error Handling
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF
      p%WaveTime = InitInp%WaveTime
      p%WaveElev0 = InitInp%WaveElev0
       
    ! Define initial system states here:
    CALL AllocAry( x%x, p%N,  'x%x', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Exc_Init')      
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

      x%x = 0
     
      xd%DummyDiscState          = 0 !TD: SS doesn't have disc states
      z%DummyConstrState         = 0 !TD: SS doesn't have constr states
      
    ! Define other States: 
      DO I=1,SIZE(OtherState%xdot)
         CALL SS_Exc_CopyContState( x, OtherState%xdot(i), MESH_NEWCOPY, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Exc_Init')
      END DO
      OtherState%n = -1

   ! misc vars:

      
     ! Inputs     
     ! no inputs

         ! Define system output initializations (set up mesh) here:
         
      y%y = 0         
      y%WriteOutput = 0
      
         
         ! Define initialization-routine output here:
         
      InitOut%WriteOutputHdr = (/ 'Time', 'FX  ' , 'FY  ' , 'FZ  ' , 'MX  ' , 'MY  ' , 'MZ  ' /)
      InitOut%WriteOutputUnt = (/ '(s) ', '(N) ' , '(N) ' , '(N) ' , '(Nm)' , '(Nm)' , '(Nm)' /)     
      
       
   CALL CleanUp() ! deallocate local arrays

CONTAINS
   SUBROUTINE CleanUp()
   
      IF (UnSS > 0 ) CLOSE ( UnSS )
      
   END SUBROUTINE CleanUp
       
END SUBROUTINE SS_Exc_Init
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
SUBROUTINE SS_Exc_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(SS_Exc_InputType),           INTENT(INOUT)  :: u           !< System inputs
      TYPE(SS_Exc_ParameterType),       INTENT(INOUT)  :: p           !< Parameters     
      TYPE(SS_Exc_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states
      TYPE(SS_Exc_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Discrete states
      TYPE(SS_Exc_ConstraintStateType), INTENT(INOUT)  :: z           !< Constraint states
      TYPE(SS_Exc_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states            
      TYPE(SS_Exc_OutputType),          INTENT(INOUT)  :: y           !< System outputs
      TYPE(SS_Exc_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables            
      INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Place any last minute operations or calculations here:
         ! Destroy the input data:
         
      CALL SS_Exc_DestroyInput( u, ErrStat, ErrMsg )


         ! Destroy the parameter data:
         
      CALL SS_Exc_DestroyParam( p, ErrStat, ErrMsg )


         ! Destroy the state data:
         
      CALL SS_Exc_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL SS_Exc_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL SS_Exc_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL SS_Exc_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )
         
         ! Destroy misc vars:
      CALL SS_Exc_DestroyMisc(  m,  ErrStat, ErrMsg )
      
      
         ! Destroy the output data:
         
      CALL SS_Exc_DestroyOutput( y, ErrStat, ErrMsg )


      

END SUBROUTINE SS_Exc_End
!----------------------------------------------------------------------------------------------------------------------------------
!> Loose coupling routine for solving constraint states, integrating continuous states, and updating discrete states.
!! Continuous, constraint, and discrete states are updated to values at t + Interval.
SUBROUTINE SS_Exc_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   ) :: t               !< Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   ) :: n               !< Current step of the simulation: t = n*Interval
      TYPE(SS_Exc_InputType),             INTENT(INOUT) :: Inputs(:)       !< Inputs at InputTimes
      REAL(DbKi),                         INTENT(IN   ) :: InputTimes(:)   !< Times in seconds associated with Inputs
      TYPE(SS_Exc_ParameterType),         INTENT(IN   ) :: p               !< Parameters
      TYPE(SS_Exc_ContinuousStateType),   INTENT(INOUT) :: x               !< Input: Continuous states at t;
                                                                           !!   Output: Continuous states at t + Interval
      TYPE(SS_Exc_DiscreteStateType),     INTENT(INOUT) :: xd              !< Input: Discrete states at t;
                                                                           !!   Output: Discrete states at t + Interval
      TYPE(SS_Exc_ConstraintStateType),   INTENT(INOUT) :: z               !< Input: Constraint states at t;
                                                                           !!   Output: Constraint states at t + Interval
      TYPE(SS_Exc_OtherStateType),        INTENT(INOUT) :: OtherState      !< Input: Other states at t;
                                                                           !!   Output: Other states at t + Interval
      TYPE(SS_Exc_MiscVarType),           INTENT(INOUT) :: m               !< Initial misc/optimization variables            
      INTEGER(IntKi),                     INTENT(  OUT) :: ErrStat         !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT) :: ErrMsg          !< Error message if ErrStat /= ErrID_None

      INTEGER, PARAMETER :: IntegrationMethod = 3   
      
                                   
      SELECT CASE ( IntegrationMethod )
         
      CASE (1) ! RK4
      
         CALL SS_Exc_RK4( t, n, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
         
      CASE (2) ! AB4
      
         CALL SS_Exc_AB4( t, n, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
      
      CASE (3) ! ABM4
      
         CALL SS_Exc_ABM4( t, n, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
         
      CASE DEFAULT  !bjj: we already checked this at initialization, but for completeness:
         
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in SS_Exc_UpdateStates: method must be 1 (RK4), 2 (AB4), or 3 (ABM4)'
         RETURN
         
      END SELECT
      
     
END SUBROUTINE SS_Exc_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
SUBROUTINE SS_Exc_CalcOutput( Time, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )   
!..................................................................................................................................
   
      REAL(DbKi),                       INTENT(IN   )   :: Time        !< Current simulation time in seconds
      TYPE(SS_Exc_InputType),           INTENT(IN   )   :: u           !< Inputs at Time
      TYPE(SS_Exc_ParameterType),       INTENT(IN   )   :: p           !< Parameters
      TYPE(SS_Exc_ContinuousStateType), INTENT(IN   )   :: x           !< Continuous states at Time
      TYPE(SS_Exc_DiscreteStateType),   INTENT(IN   )   :: xd          !< Discrete states at Time
      TYPE(SS_Exc_ConstraintStateType), INTENT(IN   )   :: z           !< Constraint states at Time
      TYPE(SS_Exc_OtherStateType),      INTENT(IN   )   :: OtherState  !< Other states at Time
      TYPE(SS_Exc_OutputType),          INTENT(INOUT)   :: y           !< Outputs computed at Time (Input only so that mesh con-
                                                                       !!   nectivity information does not have to be recalculated)
      TYPE(SS_Exc_MiscVarType),         INTENT(INOUT)   :: m           !< Initial misc/optimization variables            
      INTEGER(IntKi),                   INTENT(  OUT)   :: ErrStat     !< Error status of the operation
      CHARACTER(*),                     INTENT(  OUT)   :: ErrMsg      !< Error message if ErrStat /= ErrID_None
!      REAL(DbKi)  :: test(6,1)
      
      ! Initialize ErrStat    
      ErrStat = ErrID_None         
      ErrMsg  = ""                   

      ! Calc outputs of system, based on system states 
      ! [y] = [C]*[xr]

      y%y = matmul(p%C,x%x)    
      
      ! Compute outputs here:
      
      y%WriteOutput(1)   = REAL(Time,ReKi)
      y%WriteOutput(2:7) = y%y
                   
END SUBROUTINE SS_Exc_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for computing derivatives of continuous states
SUBROUTINE SS_Exc_CalcContStateDeriv( Time, waveElev0, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )  
!..................................................................................................................................
   
      REAL(DbKi),                        INTENT(IN   )  :: Time        !< Current simulation time in seconds
      REAL(SiKi),                        INTENT(IN   )  :: waveElev0   !< Wave elevation at origin at time: Time (m)                  
      TYPE(SS_Exc_ParameterType),        INTENT(IN   )  :: p           !< Parameters                             
      TYPE(SS_Exc_ContinuousStateType),  INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(SS_Exc_DiscreteStateType),    INTENT(IN   )  :: xd          !< Discrete states at Time
      TYPE(SS_Exc_ConstraintStateType),  INTENT(IN   )  :: z           !< Constraint states at Time
      TYPE(SS_Exc_OtherStateType),       INTENT(IN   )  :: OtherState  !< Other states                    
      TYPE(SS_Exc_MiscVarType),          INTENT(INOUT)  :: m           !< Initial misc/optimization variables            
      TYPE(SS_Exc_ContinuousStateType),  INTENT(  OUT)  :: dxdt        !< Continuous state derivatives at Time
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     !< Error status of the operation     
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
      CALL AllocAry( dxdt%x, p%N, 'SS_Exc_CalcContStateDeriv:dxdt%x', ErrStat, ErrMsg)
      IF ( ErrStat >= AbortErrLev) RETURN
            
      ! Compute the first time derivatives of the continuous states here:
      
      !Calc dxdt of a state space system
      ! [dxdt] = [A]*[xr]+B*[q]
      
      dxdt%x =matmul(p%A,x%x) +  p%B * waveElev0
        
END SUBROUTINE SS_Exc_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for updating discrete states
SUBROUTINE SS_Exc_UpdateDiscState( Time, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )   
!..................................................................................................................................
   
      REAL(DbKi),                       INTENT(IN   )  :: Time        !< Current simulation time in seconds   
      TYPE(SS_Exc_InputType),           INTENT(IN   )  :: u           !< Inputs at Time                       
      TYPE(SS_Exc_ParameterType),       INTENT(IN   )  :: p           !< Parameters                                 
      TYPE(SS_Exc_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(SS_Exc_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Input: Discrete states at Time; 
                                                                      !!   Output: Discrete states at Time + Interval
      TYPE(SS_Exc_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time
      TYPE(SS_Exc_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other/optimization states           
      TYPE(SS_Exc_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables            
      INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

               
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
         ! Update discrete states here:
      
      ! StateData%DiscState = 

END SUBROUTINE SS_Exc_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for solving for the residual of the constraint state equations
SUBROUTINE SS_Exc_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, m, z_residual, ErrStat, ErrMsg )   
!..................................................................................................................................
   
      REAL(DbKi),                       INTENT(IN   )  :: Time        !< Current simulation time in seconds   
      TYPE(SS_Exc_InputType),           INTENT(IN   )  :: u           !< Inputs at Time                       
      TYPE(SS_Exc_ParameterType),       INTENT(IN   )  :: p           !< Parameters                           
      TYPE(SS_Exc_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(SS_Exc_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at Time
      TYPE(SS_Exc_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time (possibly a guess)
      TYPE(SS_Exc_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other/optimization states                    
      TYPE(SS_Exc_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables            
      TYPE(SS_Exc_ConstraintStateType), INTENT(  OUT)  :: z_residual  !< Residual of the constraint state equations using  
                                                                      !!     the input values described above      
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat    !< Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg     !< Error message if ErrStat /= ErrID_None

               
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Solve for the constraint states here:
      
      z_residual%DummyConstrState = 0

END SUBROUTINE SS_Exc_CalcConstrStateResidual
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine implements the fourth-order Runge-Kutta Method (RK4) for numerically integrating ordinary differential equations:
!!
!!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!!   Define constants k1, k2, k3, and k4 as 
!!        k1 = dt * f(t        , x_t        )
!!        k2 = dt * f(t + dt/2 , x_t + k1/2 )
!!        k3 = dt * f(t + dt/2 , x_t + k2/2 ), and
!!        k4 = dt * f(t + dt   , x_t + k3   ).
!!   Then the continuous states at t = t + dt are
!!        x_(t+dt) = x_t + k1/6 + k2/3 + k3/3 + k4/6 + O(dt^5)
!!
!! For details, see:
!! Press, W. H.; Flannery, B. P.; Teukolsky, S. A.; and Vetterling, W. T. "Runge-Kutta Method" and "Adaptive Step Size Control for 
!!   Runge-Kutta." �16.1 and 16.2 in Numerical Recipes in FORTRAN: The Art of Scientific Computing, 2nd ed. Cambridge, England: 
!!   Cambridge University Press, pp. 704-716, 1992.
!!
SUBROUTINE SS_Exc_RK4( t, n, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                       INTENT(IN   )  :: t           !< Current simulation time in seconds
      INTEGER(IntKi),                   INTENT(IN   )  :: n           !< time step number
      REAL(DbKi),                       INTENT(IN   )  :: utimes(:)   !< times of input
      TYPE(SS_Exc_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(SS_Exc_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
      TYPE(SS_Exc_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(SS_Exc_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t (possibly a guess)
      TYPE(SS_Exc_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other/optimization states
      TYPE(SS_Exc_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables            
      INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! local variables
         
      TYPE(SS_Exc_ContinuousStateType)                 :: xdot        ! time derivatives of continuous states      
      TYPE(SS_Exc_ContinuousStateType)                 :: k1          ! RK4 constant; see above
      TYPE(SS_Exc_ContinuousStateType)                 :: k2          ! RK4 constant; see above 
      TYPE(SS_Exc_ContinuousStateType)                 :: k3          ! RK4 constant; see above 
      TYPE(SS_Exc_ContinuousStateType)                 :: k4          ! RK4 constant; see above 
      TYPE(SS_Exc_ContinuousStateType)                 :: x_tmp       ! Holds temporary modification to x
      real(SiKi)                                       :: waveElev0   ! interpolated value of the wave elevation at the origin
      INTEGER(IntKi)                                   :: ErrStat2    ! local error status
      CHARACTER(ErrMsgLen)                             :: ErrMsg2     ! local error message (ErrMsg)
      
      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      CALL SS_Exc_CopyContState( x, k1, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
      CALL SS_Exc_CopyContState( x, k2, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
      CALL SS_Exc_CopyContState( x, k3, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
      CALL SS_Exc_CopyContState( x, k4,    MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
      CALL SS_Exc_CopyContState( x, x_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN
                     
      ! find waveElev0 for time, t+p%Tc
      
      waveElev0 = InterpWrappedStpReal ( REAL(t+p%Tc, SiKi), p%WaveTime(:), p%WaveElev0(:), m%LastIndWave, p%NStepWave + 1 )        
      ! find xdot at t
      CALL SS_Exc_CalcContStateDeriv( t, waveElev0, p, x, xd, z, OtherState, m, xdot, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      k1%x     = p%dt * xdot%x  
      x_tmp%x  = x%x  + 0.5 * k1%x

      ! find waveElev0 for time, t + p%Tc + dt/2
      waveElev0 = InterpWrappedStpReal ( REAL(t+p%Tc+p%DT/2.0, SiKi), p%WaveTime(:), p%WaveElev0(:), m%LastIndWave, p%NStepWave + 1 ) 

      ! find xdot at t  + dt/2
      CALL SS_Exc_CalcContStateDeriv( t + 0.5*p%dt, waveElev0, p, x_tmp, xd, z, OtherState, m, xdot, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      k2%x     = p%dt * xdot%x
      x_tmp%x  = x%x  + 0.5 * k2%x

      ! find xdot at t + dt/2
      CALL SS_Exc_CalcContStateDeriv( t + 0.5*p%dt, waveElev0, p, x_tmp, xd, z, OtherState, m, xdot, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      k3%x     = p%dt * xdot%x 
      x_tmp%x  = x%x  + k3%x

      ! find waveElev0 for time, (t + p%Tc + dt)
      waveElev0 = InterpWrappedStpReal ( REAL(t+p%Tc+p%DT, SiKi), p%WaveTime(:), p%WaveElev0(:), m%LastIndWave, p%NStepWave + 1 )   
      

      ! find xdot at t + dt
      CALL SS_Exc_CalcContStateDeriv( t + p%dt, waveElev0, p, x_tmp, xd, z, OtherState, m, xdot, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      k4%x  = p%dt * xdot%x
      x%x   = x%x  +  ( k1%x  + 2. * k2%x  + 2. * k3%x  + k4%x  ) / 6.      

         ! clean up local variables:
      CALL ExitThisRoutine(  )
         
CONTAINS      
   !...............................................................................................................................
   SUBROUTINE ExitThisRoutine()
   ! This subroutine destroys all the local variables
   !...............................................................................................................................

         ! local variables
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)
   
   
      CALL SS_Exc_DestroyContState( xdot,     ErrStat3, ErrMsg3 )
      CALL SS_Exc_DestroyContState( k1,       ErrStat3, ErrMsg3 )
      CALL SS_Exc_DestroyContState( k2,       ErrStat3, ErrMsg3 )
      CALL SS_Exc_DestroyContState( k3,       ErrStat3, ErrMsg3 )
      CALL SS_Exc_DestroyContState( k4,       ErrStat3, ErrMsg3 )
      CALL SS_Exc_DestroyContState( x_tmp,    ErrStat3, ErrMsg3 )

         
   END SUBROUTINE ExitThisRoutine      
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)

         ! local variables
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)

      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'SS_Exc_RK4:'//TRIM(Msg)         
         ErrStat = MAX(ErrStat,ErrID)
         
         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         
         IF ( ErrStat >= AbortErrLev ) CALL ExitThisRoutine( )                  
                  
         
      END IF

   END SUBROUTINE CheckError                    
      
END SUBROUTINE SS_Exc_RK4
!-----------------------------------------------------------------------------
!! This subroutine implements the fourth-order Adams-Bashforth Method (RK4) for numerically integrating ordinary differential 
!! equations:
!!
!!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!!
!!   x(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
!!
!!  See, e.g.,
!!  http://en.wikipedia.org/wiki/Linear_multistep_method
!!
!!  or
!!
!!  K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
!!
SUBROUTINE SS_Exc_AB4( t, n, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   )  :: t           !< Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   )  :: n           !< time step number
      REAL(DbKi),                         INTENT(IN   )  :: utimes(:)   !< times of input
      TYPE(SS_Exc_ParameterType),         INTENT(IN   )  :: p           !< Parameters
      TYPE(SS_Exc_ContinuousStateType),   INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
      TYPE(SS_Exc_DiscreteStateType),     INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(SS_Exc_ConstraintStateType),   INTENT(IN   )  :: z           !< Constraint states at t (possibly a guess)
      TYPE(SS_Exc_OtherStateType),        INTENT(INOUT)  :: OtherState  !< Other/optimization states
      TYPE(SS_Exc_MiscVarType),           INTENT(INOUT)  :: m           !< Initial misc/optimization variables            
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


      ! local variables

      real(SiKi)                                         :: waveElev0   
      INTEGER(IntKi)                                     :: ErrStat2    ! local error status
      CHARACTER(ErrMsgLen)                               :: ErrMsg2     ! local error message (ErrMsg)

      
      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 
      
      
      if (OtherState%n .lt. n) then

         OtherState%n = n
                    
         CALL SS_Exc_CopyContState( OtherState%xdot ( 3 ), OtherState%xdot ( 4 ), MESH_UPDATECOPY, ErrStat2, ErrMsg )
            CALL CheckError(ErrStat2,ErrMsg2)
            IF ( ErrStat >= AbortErrLev ) RETURN
         CALL SS_Exc_CopyContState( OtherState%xdot ( 2 ), OtherState%xdot ( 3 ), MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            CALL CheckError(ErrStat2,ErrMsg2)
            IF ( ErrStat >= AbortErrLev ) RETURN
         CALL SS_Exc_CopyContState( OtherState%xdot ( 1 ), OtherState%xdot ( 2 ), MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            CALL CheckError(ErrStat2,ErrMsg2)
            IF ( ErrStat >= AbortErrLev ) RETURN
         
         
            
      elseif (OtherState%n .gt. n) then
 
         CALL CheckError( ErrID_Fatal, ' Backing up in time is not supported with a multistep method.')
         RETURN

      endif        

      ! find waveElev at  t + Tc
      waveElev0 = InterpWrappedStpReal ( REAL(t+p%Tc, SiKi), p%WaveTime(:), p%WaveElev0(:), m%LastIndWave, p%NStepWave + 1 ) 
         
      CALL SS_Exc_CalcContStateDeriv( t, waveElev0, p, x, xd, z, OtherState, m, OtherState%xdot ( 1 ), ErrStat2, ErrMsg2 ) ! initializes OtherState%xdot ( 1 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

                                                    
      if (n .le. 2) then
                                               
         CALL SS_Exc_RK4(t, n, utimes, p, x, xd, z, OtherState, m, ErrStat2, ErrMsg2 )
            CALL CheckError(ErrStat2,ErrMsg2)
            IF ( ErrStat >= AbortErrLev ) RETURN

      else

         x%x  = x%x  + p%DT/24.0  * ( 55.*OtherState%xdot(1)%x  - 59.*OtherState%xdot(2)%x   &
                                    + 37.*OtherState%xdot(3)%x   - 9.*OtherState%xdot(4)%x )


      endif

      
CONTAINS      
  
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)

         ! local variables
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)

      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'SS_Exc_AB4:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
              
         
      END IF

   END SUBROUTINE CheckError            
         
END SUBROUTINE SS_Exc_AB4
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine implements the fourth-order Adams-Bashforth-Moulton Method (RK4) for numerically integrating ordinary 
!! differential equations:
!!
!!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!!
!!   Adams-Bashforth Predictor:
!!   x^p(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
!!
!!   Adams-Moulton Corrector:
!!   x(t+dt) = x(t)  + (dt / 24.) * ( 9.*f(t+dt,x^p) + 19.*f(t,x) - 5.*f(t-dt,x) + 1.*f(t-2.*dt,x) )
!!
!!  See, e.g.,
!!  http://en.wikipedia.org/wiki/Linear_multistep_method
!!
!!  or
!!
!!  K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
SUBROUTINE SS_Exc_ABM4( t, n, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   )  :: t           !< Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   )  :: n           !< time step number
      REAL(DbKi),                         INTENT(IN   )  :: utimes(:)   !< times of input
      TYPE(SS_Exc_ParameterType),         INTENT(IN   )  :: p           !< Parameters
      TYPE(SS_Exc_ContinuousStateType),   INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
      TYPE(SS_Exc_DiscreteStateType),     INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(SS_Exc_ConstraintStateType),   INTENT(IN   )  :: z           !< Constraint states at t (possibly a guess)
      TYPE(SS_Exc_OtherStateType),        INTENT(INOUT)  :: OtherState  !< Other/optimization states
      TYPE(SS_Exc_MiscVarType),           INTENT(INOUT)  :: m           !< Initial misc/optimization variables            
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! local variables

      TYPE(SS_Exc_ContinuousStateType)                   :: x_pred      ! Continuous states at t
      TYPE(SS_Exc_ContinuousStateType)                   :: xdot_pred   ! Derivative of continuous states at t
      real(SiKi)                                         :: waveElev0
      INTEGER(IntKi)                                     :: ErrStat2    ! local error status
      CHARACTER(ErrMsgLen)                               :: ErrMsg2     ! local error message (ErrMsg)
      
      
      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      CALL SS_Exc_CopyContState(x, x_pred, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      CALL SS_Exc_AB4( t, n, utimes, p, x_pred, xd, z, OtherState, m, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      if (n .gt. 2_IntKi) then

         waveElev0 = InterpWrappedStpReal ( REAL(t+p%Tc+p%DT, SiKi), p%WaveTime(:), p%WaveElev0(:), m%LastIndWave, p%NStepWave + 1 ) 
         CALL SS_Exc_CalcContStateDeriv(t + p%dt, waveElev0, p, x_pred, xd, z, OtherState, m, xdot_pred, ErrStat2, ErrMsg2 )
            CALL CheckError(ErrStat2,ErrMsg2)
            IF ( ErrStat >= AbortErrLev ) RETURN

         x%x  = x%x  + p%DT/24. * ( 9. * xdot_pred%x +  19. * OtherState%xdot(1)%x &
                                                       - 5. * OtherState%xdot(2)%x &
                                                       + 1. * OtherState%xdot(3)%x )
                                
      else

         x%x  = x_pred%x

      endif
      
      
         ! clean up local variables:
      CALL ExitThisRoutine()
      
CONTAINS      
   !...............................................................................................................................
   SUBROUTINE ExitThisRoutine()
   ! This subroutine destroys all the local variables
   !...............................................................................................................................

         ! local variables
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)
   
   
      CALL SS_Exc_DestroyContState( xdot_pred,  ErrStat3, ErrMsg3 )
      CALL SS_Exc_DestroyContState( x_pred,     ErrStat3, ErrMsg3 )             
      
   END SUBROUTINE ExitThisRoutine    
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)

         ! local variables
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)

      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'SS_Exc_ABM4:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) CALL ExitThisRoutine( )                  
         
      END IF

   END SUBROUTINE CheckError                 

END SUBROUTINE SS_Exc_ABM4
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE SS_Excitation
!**********************************************************************************************************************************
