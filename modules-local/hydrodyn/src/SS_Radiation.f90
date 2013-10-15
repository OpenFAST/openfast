!**********************************************************************************************************************************
! The SS_Radiation and SS_Radiation_Types modules make up a template for creating user-defined calculations in the FAST Modularization 
! Framework. SS_Radiations_Types will be auto-generated based on a description of the variables for the module.
!
! "SS_Radiation" should be replaced with the name of your module. Example: HydroDyn
! "SS_Rad" (in SS_Rad_*) should be replaced with the module name or an abbreviation of it. Example: HD
!..................................................................................................................................
! LICENSING
! Copyright (C) 2012  National Renewable Energy Laboratory
!
!    This file is part of SS_Radiation.
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
! File last committed: $Date: 2013-10-03 21:39:29 -0600 (Thu, 03 Oct 2013) $
! (File) Revision #: $Rev: 266 $
! URL: $HeadURL: https://windsvn.nrel.gov/HydroDyn/branches/HydroDyn_Modularization/Source/SS_Radiation.f90 $
!**********************************************************************************************************************************
MODULE SS_Radiation

   USE SS_Radiation_Types   
   USE NWTC_Library
      
   IMPLICIT NONE
   
   PRIVATE

!   INTEGER(IntKi), PARAMETER            :: DataFormatID = 1   ! Update this value if the data types change (used in SS_Rad_Pack)
   TYPE(ProgDesc), PARAMETER            :: SS_Rad_ProgDesc = ProgDesc( 'SS_Radiation', '(v1.00.01, 19-October-2012)', '05-Mar-2013' )

   
      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: SS_Rad_Init                           ! Initialization routine
   PUBLIC :: SS_Rad_End                            ! Ending routine (includes clean up)
   
   PUBLIC :: SS_Rad_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating 
                                                    !   continuous states, and updating discrete states
   PUBLIC :: SS_Rad_CalcOutput                     ! Routine for computing outputs
   
   PUBLIC :: SS_Rad_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: SS_Rad_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: SS_Rad_UpdateDiscState                ! Tight coupling routine for updating discrete states
      
  ! PUBLIC :: SS_Rad_JacobianPInput                 ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
                                                    !   (Xd), and constraint-state (Z) equations all with respect to the inputs (u)
   !PUBLIC :: SS_Rad_JacobianPContState             ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
                                                    !   (Xd), and constraint-state (Z) equations all with respect to the continuous 
                                                    !   states (x)
 !  PUBLIC :: SS_Rad_JacobianPDiscState             ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
                                                    !   (Xd), and constraint-state (Z) equations all with respect to the discrete 
                                                    !   states (xd)
!   PUBLIC :: SS_Rad_JacobianPConstrState           ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
                                                    !   (Xd), and constraint-state (Z) equations all with respect to the constraint 
                                                    !   states (z)
   !PUBLIC :: Solver
   
! Note that the following routines will be updated with new definitions of arrays returned (no longer one-byte arrays)
   !PUBLIC :: SS_Rad_Pack                           ! Routine to pack (save) data into one array of bytes 
   !PUBLIC :: SS_Rad_Unpack                         ! Routine to unpack an array of bytes into data structures usable by the module
   
CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SS_Rad_Init( InitInp, u, p, x, xd, z, OtherState, y, Interval, InitOut, ErrStat, ErrMsg )
! This routine is called at the start of the simulation to perform initialization steps. 
! The parameters are set here and not changed during the simulation.
! The initial states and initial guess for the input are defined.
!..................................................................................................................................

    TYPE(SS_Rad_InitInputType),       INTENT(IN   )  :: InitInp     ! Input data for initialization routine
    TYPE(SS_Rad_InputType),           INTENT(  OUT)  :: u           ! An initial guess for the input; input mesh must be defined
    TYPE(SS_Rad_ParameterType),       INTENT(  OUT)  :: p           ! Parameters      
    TYPE(SS_Rad_ContinuousStateType), INTENT(  OUT)  :: x           ! Initial continuous states
    TYPE(SS_Rad_DiscreteStateType),   INTENT(  OUT)  :: xd          ! Initial discrete states
    TYPE(SS_Rad_ConstraintStateType), INTENT(  OUT)  :: z           ! Initial guess of the constraint states
    TYPE(SS_Rad_OtherStateType),      INTENT(  OUT)  :: OtherState  ! Initial other/optimization states            
    TYPE(SS_Rad_OutputType),          INTENT(  OUT)  :: y           ! Initial system outputs (outputs are not calculated; 
                                                                    !   only the output mesh is initialized)
    REAL(DbKi),                       INTENT(INOUT)  :: Interval    ! Coupling interval in seconds: the rate that 
                                                                    !   (1) SS_Rad_UpdateStates() is called in loose coupling &
                                                                    !   (2) SS_Rad_UpdateDiscState() is called in tight coupling.
                                                                    !   Input is the suggested time from the glue code; 
                                                                    !   Output is the actual coupling interval that will be used 
                                                                    !   by the glue code.
    TYPE(SS_Rad_InitOutputType),      INTENT(  OUT)  :: InitOut     ! Output for initialization routine
    INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     ! Error status of the operation
    CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

    ! Local Variables:
         
    REAL(ReKi), ALLOCATABLE                :: Rad_A (:,:)                          ! A matrix of the radiation state-space system on the input file ss
    REAL(ReKi), ALLOCATABLE                :: Rad_B (:,:)                          ! B matrix of the radiation state-space system on the input file ss
    REAL(ReKi), ALLOCATABLE                :: Rad_C (:,:)                          ! C matrix of the radiation state-space system on the input file ss

    INTEGER                                :: I                                    ! Generic index
    INTEGER                                :: J                                    ! Generic index  
    INTEGER                                :: xx (1,6)                             ! Active DOF's on the input file .ss
    INTEGER(IntKi)                         :: spdof (1,6)                          ! States per dof  
    INTEGER                                :: DOFs                                 ! Number of DOFS  
    INTEGER                                :: N                                    ! Number of states
    INTEGER                                :: Nlines                               ! Number of lines in the input file, used to determine N
    INTEGER                                :: UnSS       = 34                      ! I/O unit number for the WAMIT output file with the .ss extension; this file contains the state-space matrices.
    INTEGER                                :: Sttus                                ! Error in reading .ss file
    CHARACTER                              :: Line                                 ! Temp line of file
    
    ! Initialize ErrStat   
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
    ! Open the .ss input file!
    CALL OpenFInpFile ( UnSS, TRIM(InitInp%InputFile)//'.ss', Sttus )  ! Open file.
    IF (Sttus==1) THEN
        ErrStat = ErrID_Fatal
        ErrMsg  = 'Unable to open the input file .ss' 
    END IF        

    ! Determine the number of states and size of the matrices
    Nlines = 1
    
    CALL ReadCom ( UnSS, TRIM(InitInp%InputFile)//'.ss', 'Header',ErrStat  )! Reads the first entire line (Title header)
    CALL ReadAry( UnSS,TRIM(InitInp%InputFile)//'.ss', xx(1,:), 6, 'xx', 'xx vector containing the enabled dofs',ErrStat, ErrMsg) ! Reads in the second line, containing the active dofs vector
    CALL ReadVar( UnSS,TRIM(InitInp%InputFile)//'.ss', N, 'N', 'Number of Dofs',ErrStat, ErrMsg) ! Reads in the third line, containing the number of states
    CALL ReadAry( UnSS,TRIM(InitInp%InputFile)//'.ss', spdof(1,:), 6, 'spdof', 'spdof vector containing the number of states per dofs',ErrStat, ErrMsg) ! Reads in the forth line, containing the state per dofs vector

    DO !Loop through all the lines of the file
        CALL ReadCom ( UnSS, TRIM(InitInp%InputFile)//'.ss', 'Header',Sttus,ErrMsg  )! Reads the first entire line (Title header)
        IF ( Sttus == 0 )  THEN ! .TRUE. when data is read in successfully                    
            Nlines=Nlines+1                    
        ELSE !We must have reach the end of the file
            EXIT
        END IF
    END DO

    ! The input file contains the matrices A [NxN], B [Nx6] and C [6xN], so
    !p%N = ( Nlines - 6 ) / 2 ! this is the number of states
    
    !Verifications on the input file
    IF ( ( Nlines - 6 ) / 2 /= N) THEN
        ErrStat = ErrID_Severe 
        ErrMsg  = 'Error in the input file .ss: The size of the matrices does not correspond to the number of states!' 
    END IF
    
    IF ( N /= SUM(spdof)) THEN
        ErrStat = ErrID_Severe 
        ErrMsg  = 'Error in the input file .ss: The size of the matrices does not correspond to the number of states!' 
    END IF        
    
    !Verify if the DOFs active in the input file correspond to the ones active by FAST in this run
    DO I=1,6 !Loop through all 6 DOFs           
        IF ( InitInp%DOFs (1,I) == 1)  THEN !  True when the current DOF is active in FAST                   
            IF ( xx (1,I) /= 1) THEN ! True if a DOF enabled by FAST is not available in the INPUT File
                ErrStat = ErrID_Severe 
                ErrMsg  = 'Error in the input file .ss: The enabled DOFs in the current FAST Simulation don`t match the ones on the input file .ss!' 
            END IF           
        END IF
    END DO
    
    DOFs = SUM (xx) !Number of DOFS in the input file
    
    ! Now we can allocate the temporary matrices A, B and C
    
    ALLOCATE (Rad_A  (N,N) , STAT=Sttus )
    IF ( Sttus /= 0 )  THEN
        ErrStat = ErrID_Fatal
        ErrMsg  = ' Error allocating memory for the Rad_A array.'
    END IF
    
    ALLOCATE ( Rad_B (N,DOFs) , STAT=Sttus )
    IF ( Sttus /= 0 )  THEN
        ErrStat = ErrID_Fatal
        ErrMsg  = ' Error allocating memory for the Rad_B array.'
    END IF
    
    ALLOCATE ( Rad_C  (DOFs,N) , STAT=Sttus )
    IF ( Sttus /= 0 )  THEN
        ErrStat = ErrID_Fatal
        ErrMsg  = ' Error allocating memory for the Rad_C array.'
    END IF   
    
    REWIND (UNIT=UnSS)   ! REWIND the file so we can read it in a second time.

    ! Skip the first 4 lines:
    CALL ReadCom ( UnSS, TRIM(InitInp%InputFile)//'.ss', 'Header'  )! Reads the first entire line (Title header)
    CALL ReadCom ( UnSS, TRIM(InitInp%InputFile)//'.ss', 'Enabled dofs'  )! Reads the first entire line (Title header)
    CALL ReadCom ( UnSS, TRIM(InitInp%InputFile)//'.ss', 'N'  )! Reads the first entire line (Title header)
    CALL ReadCom ( UnSS, TRIM(InitInp%InputFile)//'.ss', 'N per dofs'  )! Reads the first entire line (Title header)   
    
    DO I = 1,N !Read A MatriX
        CALL ReadAry( UnSS,TRIM(InitInp%InputFile)//'.ss', Rad_A(I,:), N, 'Rad_A', 'A_Matrix',ErrStat, ErrMsg)   
    END DO
    
    DO I = 1,N !Read B Matrix
        CALL ReadAry( UnSS, TRIM(InitInp%InputFile)//'.ss', Rad_B(I,:), 6, 'Rad_B', 'B_Matrix',ErrStat, ErrMsg)     
    END DO
    
    DO I = 1,6 !Read C Matrix
        CALL ReadAry( UnSS, TRIM(InitInp%InputFile)//'.ss', Rad_C(I,:), N, 'Rad_C', 'C_Matrix',ErrStat, ErrMsg)
    END DO
    
    CLOSE ( UnSS ) !Close .ss input file
    
    !Now we are ready to reduce the matrices to the correspondent active dofs in FAST
    p%N=0
    DO I=1,6 !For each state
        IF ( InitInp%DOFs (1,I) == 1)  THEN !  True when the current DOF is active in FAST          
            p%N = p%N + spdof(1,I) !Add the correspondent number of states to the vector
        END IF
    END DO
    
    CALL WrScr1 ( 'Using SS_Radiation Module, with '//TRIM( Num2LStr(p%N ))//' of '//TRIM( Num2LStr(N ))// ' radiation states' )
    
    !Now we can allocate the final size of the SS matrices
    ALLOCATE (p%A  (p%N,p%N) , STAT=Sttus )
    IF ( Sttus /= 0 )  THEN
        ErrStat = ErrID_Fatal
        ErrMsg  = ' Error allocating memory for the A array.'
    END IF
    
    ALLOCATE ( p%B (p%N,6) , STAT=Sttus )
    IF ( Sttus /= 0 )  THEN
        ErrStat = ErrID_Fatal
        ErrMsg  = ' Error allocating memory for the B array.'
    END IF
    
    ALLOCATE ( p%C  (6,p%N) , STAT=Sttus )
    IF ( Sttus /= 0 )  THEN
        ErrStat = ErrID_Fatal
        ErrMsg  = ' Error allocating memory for the C array.'
    END IF   
            
    !Finaly we write the ss matrices, based on the ones on the input file and on the active dofs
        p%A = 0
        p%B = 0
        p%C = 0
        
    IF ( p%N == N ) THEN !The matrices are the same
        
        p%A = Rad_A
        p%B = Rad_B
        p%C = Rad_C
            
    ELSE !We need to cut some of the lines and columns
        N=1 !Use as number of active states introduced
        
        DO I=1,6 !For each dof...
            IF ( InitInp%DOFs (1,I) == 1 .AND. sum(spdof(1,1:I))<size(Rad_A(:,1)))  THEN !  That is enabled in FAST
    
                p%A (N:N+spdof(1,I),N:N+spdof(1,I)) = Rad_A (sum(spdof(1,1:I-1))+1:sum(spdof(1,1:I)),sum(spdof(1,1:I-1))+1:sum(spdof(1,1:I)))
                p%B (N:N+spdof(1,I),:)= Rad_B (sum(spdof(1,1:I-1))+1:sum(spdof(1,1:I)),:)
                p%C (:,N:N+spdof(1,I))= Rad_C (:,sum(spdof(1,1:I-1))+1:sum(spdof(1,1:I)))
                
                N = N + spdof(1,I) !Number of lines added to the A and B Matrix and columns to the C Matrix
            END IF
        END DO
    END IF
            
    ! Define parameters here:
         
      p%DT  = Interval
       
    ! Define initial system states here:

    ALLOCATE ( x%x  (p%N,1) , STAT=Sttus )
    IF ( Sttus /= 0 )  THEN
        ErrStat = ErrID_Fatal
        ErrMsg  = ' Error allocating memory for the State vector.'
    END IF         
      
      x%x = 0
     
      xd%DummyDiscState          = 0 !TD: SS doesn't have disc states
      z%DummyConstrState         = 0 !TD: SS doesn't have constr states
      
    ! Define other States:    
    ALLOCATE ( OtherState%dxdt  (p%N,4) , STAT=Sttus )
    IF ( Sttus /= 0 )  THEN
        ErrStat = ErrID_Fatal
        ErrMsg  = ' Error allocating memory for the OtherState%dxdt array.'
    END IF
     OtherState%dxdt = 0   
     OtherState%Step = 0
     OtherState%LastTime = 0    ! Define initial guess for the system inputs here:

     !Inputs     
      u%dq = 0 !6 DoF's velocities

         ! Define system output initializations (set up mesh) here:
         
      y%y = 0         
      y%WriteOutput = 0
      
         
         ! Define initialization-routine output here:
         
      InitOut%WriteOutputHdr = (/ 'Time', 'F1  ' , 'F2  ' , 'F3  ' , 'F4  ' , 'F5  ' , 'F6  ' /)
      InitOut%WriteOutputUnt = (/ '(s) ',  '(N) '   ,  '(N) '   ,  '(N) '   ,  '(Nm)'   ,  '(Nm)',  '(Nm)'        /)     
      
         ! If you want to choose your own rate instead of using what the glue code suggests, tell the glue code the rate at which
         !   this module must be called here:
         
       !p%DT=Interval

END SUBROUTINE SS_Rad_Init
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SS_Rad_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
! This routine is called at the end of the simulation.
!..................................................................................................................................

      TYPE(SS_Rad_InputType),           INTENT(INOUT)  :: u           ! System inputs
      TYPE(SS_Rad_ParameterType),       INTENT(INOUT)  :: p           ! Parameters     
      TYPE(SS_Rad_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states
      TYPE(SS_Rad_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Discrete states
      TYPE(SS_Rad_ConstraintStateType), INTENT(INOUT)  :: z           ! Constraint states
      TYPE(SS_Rad_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states            
      TYPE(SS_Rad_OutputType),          INTENT(INOUT)  :: y           ! System outputs
      INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Place any last minute operations or calculations here:
         ! Destroy the input data:
         
      CALL SS_Rad_DestroyInput( u, ErrStat, ErrMsg )


         ! Destroy the parameter data:
         
      CALL SS_Rad_DestroyParam( p, ErrStat, ErrMsg )


         ! Destroy the state data:
         
      CALL SS_Rad_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL SS_Rad_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL SS_Rad_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL SS_Rad_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )
         

         ! Destroy the output data:
         
      CALL SS_Rad_DestroyOutput( y, ErrStat, ErrMsg )


      

END SUBROUTINE SS_Rad_End
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SS_Rad_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
! Loose coupling routine for solving constraint states, integrating continuous states, and updating discrete states.
! Continuous, constraint, and discrete states are updated to values at t + Interval.
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   ) :: t               ! Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   ) :: n               ! Current step of the simulation: t = n*Interval
      TYPE(SS_Rad_InputType),            INTENT(INOUT   ) :: Inputs(:)       ! Inputs at InputTimes
      REAL(DbKi),                         INTENT(IN   ) :: InputTimes(:)   ! Times in seconds associated with Inputs
      TYPE(SS_Rad_ParameterType),        INTENT(IN   ) :: p               ! Parameters
      TYPE(SS_Rad_ContinuousStateType),  INTENT(INOUT) :: x               ! Input: Continuous states at t;
                                                                           !   Output: Continuous states at t + Interval
      TYPE(SS_Rad_DiscreteStateType),    INTENT(INOUT) :: xd              ! Input: Discrete states at t;
                                                                           !   Output: Discrete states at t + Interval
      TYPE(SS_Rad_ConstraintStateType),  INTENT(INOUT) :: z               ! Input: Constraint states at t;
                                                                           !   Output: Constraint states at t + Interval
      TYPE(SS_Rad_OtherStateType),       INTENT(INOUT) :: OtherState      ! Other/optimization states
      INTEGER(IntKi),                     INTENT(  OUT) :: ErrStat         ! Error status of the operation
      CHARACTER(*),                       INTENT(  OUT) :: ErrMsg          ! Error message if ErrStat /= ErrID_None

         
         ! Local variables     
      TYPE(SS_Rad_ContinuousStateType)                 :: dxdt        ! Continuous state derivatives at Time
      TYPE(SS_Rad_ConstraintStateType)                 :: z_Residual  ! Residual of the constraint state equations 
      TYPE(SS_Rad_InputType)                          :: u               ! Instantaneous inputs
      INTEGER(IntKi)                                   :: I           ! loop indice
      INTEGER(IntKi)                                   :: ErrStat2    ! Error status of the operation (occurs after initial error)
      CHARACTER(LEN(ErrMsg))                           :: ErrMsg2     ! Error message if ErrStat2 /= ErrID_None
                        
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""             
      
      
       ! Get the inputs at time t, based on the array of values sent by the glue code:
         
      CALL SS_Rad_Input_ExtrapInterp( Inputs, InputTimes, u, t, ErrStat, ErrMsg )  
      IF ( ErrStat >= AbortErrLev ) RETURN
      
    ! Integrate (update) continuous states (x) here
    
    !TD First order solver:
    ! Get first time derivatives of continuous states (dxdt): 
    
    CALL SS_Rad_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )
    IF ( ErrStat >= AbortErrLev ) THEN        
         CALL SS_Rad_DestroyContState( dxdt, ErrStat2, ErrMsg2)
         ErrMsg = TRIM(ErrMsg)//' '//TRIM(ErrMsg2)
         RETURN
    END IF         
     !x%x  = x%x + dxdt%x * p%DT
         
     ! 4th order solver
      
     !IF ( EqualRealNos(OtherState%LastTime + p%DT,Time)) THEN !Time must have been reduced TD: Function EqualRealNos only works for Single!
     IF ( OtherState%LastTime + p%DT/=t) THEN !Time must have been reduced
         
         OtherState%dxdt = 0 ! Remove previous history and start from zero with a
                              ! runge-Kutta method
         OtherState%Step = 0
     END IF
     
    !Update time step
     OtherState%Step = OtherState%Step + 1 
    !Update the OtherStates matrices, with the previous dXdt Values
        OtherState%dxdt (:,4) = OtherState%dxdt (:,3)
        OtherState%dxdt (:,3) = OtherState%dxdt (:,2)
        OtherState%dxdt (:,2) = OtherState%dxdt (:,1)
        OtherState%dxdt (:,1) = dxdt%x (:,1)
    
    Call Solver (t, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg)
       
     !Update LastTime
     OtherState%LastTime = t
        
    !Destroy dxdt because it is not necessary for the rest of the subroutine
      CALL SS_Rad_DestroyContState( dxdt, ErrStat, ErrMsg)
      
   
      IF ( ErrStat >= AbortErrLev ) RETURN      
     
END SUBROUTINE SS_Rad_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SS_Rad_CalcOutput( Time, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )   
! Routine for computing outputs, used in both loose and tight coupling.
!..................................................................................................................................
   
      REAL(DbKi),                        INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(SS_Rad_InputType),           INTENT(IN   )   :: u           ! Inputs at Time
      TYPE(SS_Rad_ParameterType),       INTENT(IN   )   :: p           ! Parameters
      TYPE(SS_Rad_ContinuousStateType), INTENT(IN   )   :: x           ! Continuous states at Time
      TYPE(SS_Rad_DiscreteStateType),   INTENT(IN   )   :: xd          ! Discrete states at Time
      TYPE(SS_Rad_ConstraintStateType), INTENT(IN   )   :: z           ! Constraint states at Time
      TYPE(SS_Rad_OtherStateType),      INTENT(INOUT)   :: OtherState  ! Other/optimization states
      TYPE(SS_Rad_OutputType),          INTENT(INOUT)   :: y           ! Outputs computed at Time (Input only so that mesh con-
                                                                       !   nectivity information does not have to be recalculated)
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
      REAL(DbKi)  :: test(6,1)
      
      ! Initialize ErrStat    
      ErrStat = ErrID_None         
      ErrMsg  = ""                   

      ! Calc outputs of system, based on system states 
      ! [y] = [C]*[xr]

      y%y = matmul(p%C,x%x)    
      
      ! Compute outputs here:
      
      y%WriteOutput(1,1)   = REAL(Time,ReKi)
      y%WriteOutput(1,2:7) = y%y(:,1)
                   
END SUBROUTINE SS_Rad_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SS_Rad_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )  
! Tight coupling routine for computing derivatives of continuous states
!..................................................................................................................................
   
      REAL(DbKi),                        INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(SS_Rad_InputType),           INTENT(IN   )  :: u           ! Inputs at Time                    
      TYPE(SS_Rad_ParameterType),       INTENT(IN   )  :: p           ! Parameters                             
      TYPE(SS_Rad_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(SS_Rad_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(SS_Rad_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(SS_Rad_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states                    
      TYPE(SS_Rad_ContinuousStateType), INTENT(  OUT)  :: dxdt        ! Continuous state derivatives at Time
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation     
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      ALLOCATE ( dxdt%x  (p%N,1) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
          ErrStat = ErrID_Fatal
          ErrMsg  = ' Error allocating memory for the State vector.'
      END IF         
      
      dxdt%x = 0
      
      ! Compute the first time derivatives of the continuous states here:
      
      !Calc dxdt of a state space system
      ! [dxdt] = [A]*[xr]+B*[q]
      
      dxdt%x =matmul(p%A,x%x) + matmul( p%B, u%dq)
        
END SUBROUTINE SS_Rad_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SS_Rad_UpdateDiscState( Time, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )   
! Tight coupling routine for updating discrete states
!..................................................................................................................................
   
      REAL(DbKi),                        INTENT(IN   )  :: Time        ! Current simulation time in seconds   
      TYPE(SS_Rad_InputType),           INTENT(IN   )  :: u           ! Inputs at Time                       
      TYPE(SS_Rad_ParameterType),       INTENT(IN   )  :: p           ! Parameters                                 
      TYPE(SS_Rad_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(SS_Rad_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Input: Discrete states at Time; 
                                                                       !   Output: Discrete states at Time + Interval
      TYPE(SS_Rad_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(SS_Rad_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states           
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

               
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
         ! Update discrete states here:
      
      ! StateData%DiscState = 

END SUBROUTINE SS_Rad_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SS_Rad_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_residual, ErrStat, ErrMsg )   
! Tight coupling routine for solving for the residual of the constraint state equations
!..................................................................................................................................
   
      REAL(DbKi),                        INTENT(IN   )  :: Time        ! Current simulation time in seconds   
      TYPE(SS_Rad_InputType),           INTENT(IN   )  :: u           ! Inputs at Time                       
      TYPE(SS_Rad_ParameterType),       INTENT(IN   )  :: p           ! Parameters                           
      TYPE(SS_Rad_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(SS_Rad_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(SS_Rad_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time (possibly a guess)
      TYPE(SS_Rad_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states                    
      TYPE(SS_Rad_ConstraintStateType), INTENT(  OUT)  :: z_residual  ! Residual of the constraint state equations using  
                                                                       !     the input values described above      
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

               
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Solve for the constraint states here:
      
      z_residual%DummyConstrState = 0

END SUBROUTINE SS_Rad_CalcConstrStateResidual
!----------------------------------------------------------------------------------------------------------------------------------
!SUBROUTINE SS_Rad_JacobianPInput( Time, u, p, x, xd, z, OtherState, dYdu, dXdu, dXddu, dZdu, ErrStat, ErrMsg )   
!! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) equations 
!! with respect to the inputs (u). The partial derivatives dY/du, dX/du, dXd/du, and DZ/du are returned.
!!..................................................................................................................................
!   
!      REAL(DbKi),                                INTENT(IN   )           :: Time       ! Current simulation time in seconds   
!      TYPE(SS_Rad_InputType),                   INTENT(IN   )           :: u          ! Inputs at Time                       
!      TYPE(SS_Rad_ParameterType),               INTENT(IN   )           :: p          ! Parameters                           
!      TYPE(SS_Rad_ContinuousStateType),         INTENT(IN   )           :: x          ! Continuous states at Time
!      TYPE(SS_Rad_DiscreteStateType),           INTENT(IN   )           :: xd         ! Discrete states at Time
!      TYPE(SS_Rad_ConstraintStateType),         INTENT(IN   )           :: z          ! Constraint states at Time
!      TYPE(SS_Rad_OtherStateType),              INTENT(INOUT)           :: OtherState ! Other/optimization states                    
!      TYPE(SS_Rad_PartialOutputPInputType),     INTENT(  OUT), OPTIONAL :: dYdu       ! Partial derivatives of output equations
!                                                                                       !   (Y) with respect to the inputs (u)
!      TYPE(SS_Rad_PartialContStatePInputType),  INTENT(  OUT), OPTIONAL :: dXdu       ! Partial derivatives of continuous state
!                                                                                       !   equations (X) with respect to inputs (u)
!      TYPE(SS_Rad_PartialConstrStatePInputType ),  INTENT(  OUT), OPTIONAL :: dXddu      ! Partial derivatives of discrete state 
!                                                                                       !   equations (Xd) with respect to inputs (u)
!      TYPE(SS_Rad_PartialConstrStatePInputType),INTENT(  OUT), OPTIONAL :: dZdu       ! Partial derivatives of constraint state 
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
!         dYdu%dYdu = 0
!
!      END IF
!      
!      IF ( PRESENT( dXdu ) ) THEN
!      
!         ! Calculate the partial derivative of the continuous state equations (X) with respect to the inputs (u) here:
!      
!         dXdu%dXdu = p%B
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
!END SUBROUTINE SS_Rad_JacobianPInput
!!----------------------------------------------------------------------------------------------------------------------------------
!SUBROUTINE SS_Rad_JacobianPContState( Time, u, p, x, xd, z, OtherState, dYdx, dXdx, dXddx, dZdx, ErrStat, ErrMsg )   
!! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) equations
!! with respect to the continuous states (x). The partial derivatives dY/dx, dX/dx, dXd/dx, and DZ/dx are returned.
!!..................................................................................................................................
!   
!      REAL(DbKi),                                    INTENT(IN   )           :: Time       ! Current simulation time in seconds   
!      TYPE(SS_Rad_InputType),                       INTENT(IN   )           :: u          ! Inputs at Time                       
!      TYPE(SS_Rad_ParameterType),                   INTENT(IN   )           :: p          ! Parameters                           
!      TYPE(SS_Rad_ContinuousStateType),             INTENT(IN   )           :: x          ! Continuous states at Time
!      TYPE(SS_Rad_DiscreteStateType),               INTENT(IN   )           :: xd         ! Discrete states at Time
!      TYPE(SS_Rad_ConstraintStateType),             INTENT(IN   )           :: z          ! Constraint states at Time
!      TYPE(SS_Rad_OtherStateType),                  INTENT(INOUT)           :: OtherState ! Other/optimization states                    
!      TYPE(SS_Rad_PartialOutputPContStateType),     INTENT(  OUT), OPTIONAL :: dYdx       ! Partial derivatives of output equations
!                                                                                           !   (Y) with respect to the continuous 
!                                                                                           !   states (x)
!      TYPE(SS_Rad_PartialContStatePContStateType),  INTENT(  OUT), OPTIONAL :: dXdx       ! Partial derivatives of continuous state
!                                                                                           !   equations (X) with respect to 
!                                                                                           !   the continuous states (x)
!      TYPE(SS_Rad_PartialDiscStatePContStateType),  INTENT(  OUT), OPTIONAL :: dXddx      ! Partial derivatives of discrete state 
!                                                                                           !   equations (Xd) with respect to 
!                                                                                           !   the continuous states (x)
!      TYPE(SS_Rad_PartialConstrStatePContStateType),INTENT(  OUT), OPTIONAL :: dZdx       ! Partial derivatives of constraint state
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
!         dYdx%dYdx = p%C
!
!      END IF
!      
!      IF ( PRESENT( dXdx ) ) THEN
!      
!         ! Calculate the partial derivative of the continuous state equations (X) with respect to the continuous states (x) here:
!      
!         dXdx%dXdx = p%A
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
!   END SUBROUTINE SS_Rad_JacobianPContState
!!----------------------------------------------------------------------------------------------------------------------------------
!SUBROUTINE SS_Rad_JacobianPDiscState( Time, u, p, x, xd, z, OtherState, dYdxd, dXdxd, dXddxd, dZdxd, ErrStat, ErrMsg )   
!! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) equations
!! with respect to the discrete states (xd). The partial derivatives dY/dxd, dX/dxd, dXd/dxd, and DZ/dxd are returned.
!!..................................................................................................................................
!
!      REAL(DbKi),                                    INTENT(IN   )           :: Time       ! Current simulation time in seconds   
!      TYPE(SS_Rad_InputType),                       INTENT(IN   )           :: u          ! Inputs at Time                       
!      TYPE(SS_Rad_ParameterType),                   INTENT(IN   )           :: p          ! Parameters                           
!      TYPE(SS_Rad_ContinuousStateType),             INTENT(IN   )           :: x          ! Continuous states at Time
!      TYPE(SS_Rad_DiscreteStateType),               INTENT(IN   )           :: xd         ! Discrete states at Time
!      TYPE(SS_Rad_ConstraintStateType),             INTENT(IN   )           :: z          ! Constraint states at Time
!      TYPE(SS_Rad_OtherStateType),                  INTENT(INOUT)           :: OtherState ! Other/optimization states                    
!      TYPE(SS_Rad_PartialOutputPDiscStateType),     INTENT(  OUT), OPTIONAL :: dYdxd      ! Partial derivatives of output equations
!                                                                                           !  (Y) with respect to the discrete 
!                                                                                           !  states (xd)
!      TYPE(SS_Rad_PartialContStatePDiscStateType),  INTENT(  OUT), OPTIONAL :: dXdxd      ! Partial derivatives of continuous state
!                                                                                           !   equations (X) with respect to the 
!                                                                                           !   discrete states (xd)
!      TYPE(SS_Rad_PartialDiscStatePDiscStateType),  INTENT(  OUT), OPTIONAL :: dXddxd     ! Partial derivatives of discrete state 
!                                                                                           !   equations (Xd) with respect to the
!                                                                                           !   discrete states (xd)
!      TYPE(SS_Rad_PartialConstrStatePDiscStateType),INTENT(  OUT), OPTIONAL :: dZdxd      ! Partial derivatives of constraint state
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
!END SUBROUTINE SS_Rad_JacobianPDiscState
!!----------------------------------------------------------------------------------------------------------------------------------    
!SUBROUTINE SS_Rad_JacobianPConstrState( Time, u, p, x, xd, z, OtherState, dYdz, dXdz, dXddz, dZdz, ErrStat, ErrMsg )   
!! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) equations
!! with respect to the constraint states (z). The partial derivatives dY/dz, dX/dz, dXd/dz, and DZ/dz are returned.
!!..................................................................................................................................
!   
!      REAL(DbKi),                                      INTENT(IN   )           :: Time       ! Current simulation time in seconds   
!      TYPE(SS_Rad_InputType),                         INTENT(IN   )           :: u          ! Inputs at Time                       
!      TYPE(SS_Rad_ParameterType),                     INTENT(IN   )           :: p          ! Parameters                           
!      TYPE(SS_Rad_ContinuousStateType),               INTENT(IN   )           :: x          ! Continuous states at Time
!      TYPE(SS_Rad_DiscreteStateType),                 INTENT(IN   )           :: xd         ! Discrete states at Time
!      TYPE(SS_Rad_ConstraintStateType),               INTENT(IN   )           :: z          ! Constraint states at Time
!      TYPE(SS_Rad_OtherStateType),                    INTENT(INOUT)           :: OtherState ! Other/optimization states                    
!      TYPE(SS_Rad_PartialOutputPConstrStateType),     INTENT(  OUT), OPTIONAL :: dYdz       ! Partial derivatives of output 
!                                                                                             !  equations (Y) with respect to the 
!                                                                                             !  constraint states (z)
!      TYPE(SS_Rad_PartialContStatePConstrStateType),  INTENT(  OUT), OPTIONAL :: dXdz       ! Partial derivatives of continuous
!                                                                                             !  state equations (X) with respect to 
!                                                                                             !  the constraint states (z)
!      TYPE(SS_Rad_PartialDiscStatePConstrStateType),  INTENT(  OUT), OPTIONAL :: dXddz      ! Partial derivatives of discrete state
!                                                                                             !  equations (Xd) with respect to the 
!                                                                                             !  constraint states (z)
!      TYPE(SS_Rad_PartialConstrStatePConstrStateType),INTENT(  OUT), OPTIONAL :: dZdz       ! Partial derivatives of constraint 
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
!END SUBROUTINE SS_Rad_JacobianPConstrState
!----------------------------------------------------------------------------------------------------------------------------------
!     
SUBROUTINE Solver(Time, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )
    ! Solver solves the equations of motion by marching in time using a
    !   predictor-corrector scheme.  Fourth order Runge-Kutta is used to
    !   get the first 4 points from the initial degrees of freedom and
    !   velocities.

IMPLICIT                        NONE

      REAL(DbKi),                        INTENT(IN   ) :: Time        ! Current simulation time in seconds
      TYPE(SS_Rad_InputType),            INTENT(IN   ) :: u           ! Inputs at Time                    
      TYPE(SS_Rad_ParameterType),        INTENT(IN   ) :: p           ! Parameters                              
      TYPE(SS_Rad_ContinuousStateType),  INTENT(INOUT) :: x           ! Input: Continuous states at Time; 
                                                                      !   Output: Continuous states at Time + Interval
      TYPE(SS_Rad_DiscreteStateType),    INTENT(INOUT) :: xd          ! Input: Discrete states at Time; 
                                                                      !   Output: Discrete states at Time  + Interval
      TYPE(SS_Rad_ConstraintStateType),  INTENT(INOUT) :: z           ! Input: Initial guess of constraint states at Time;
                                                                      !   Output: Constraint states at Time
      TYPE(SS_Rad_OtherStateType),       INTENT(INOUT) :: OtherState  ! Other/optimization states
      TYPE(SS_Rad_ContinuousStateType),  INTENT(INOUT) :: dxdt        ! Continuous state derivatives at Time
      INTEGER(IntKi),                    INTENT(  OUT) :: ErrStat     ! Error status of the operation     
      CHARACTER(*),                      INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None
        
        
        ! Local variables       
        INTEGER(IntKi)                          :: I                  ! Loops through all DOFs
        INTEGER(IntKi)                          :: Sttus              ! Status returned from an attempt to allocate an array.
        REAL(DbKi), ALLOCATABLE                 ::    ZK1 (:)         ! Temporary values of dq                
        REAL(DbKi), ALLOCATABLE                 ::    ZK2 (:)         ! Temporary values of dq                
        REAL(DbKi), ALLOCATABLE                 ::    ZK3 (:)         ! Temporary values of dq                
        REAL(DbKi), ALLOCATABLE                 ::    ZK4 (:)         ! Temporary values of dq                
        INTEGER(IntKi)                          :: ErrStat2           ! Error status of the operation (occurs after initial error)
        CHARACTER(LEN(ErrMsg))                  :: ErrMsg2            ! Error message if ErrStat2 /= ErrID_None
        TYPE(SS_Rad_ContinuousStateType)        :: xt        ! Dummy Continuous states at Time+1, created by predictor
       
      ErrStat = ErrID_None
      ErrMsg  = ''
        
IF ( OtherState%Step < 3 )  THEN   ! Use Runge-Kutta integration at the the start of the simulation (first 3 steps).


   ! Allocate arrays that vary with the number of DOFs..
   Sttus = 0

   IF (.NOT. ALLOCATED(ZK1)) ALLOCATE ( ZK1(p%N) , STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = ' Error allocating memory for the ZK1 array.'
   END IF

   IF (.NOT. ALLOCATED(ZK2)) ALLOCATE ( ZK2(p%N) , STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = ' Error allocating memory for the ZK2 array.'
   END IF

   IF (.NOT. ALLOCATED(ZK3)) ALLOCATE ( ZK3(p%N) , STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = ' Error allocating memory for the ZK3 array.'
   END IF

   IF (.NOT. ALLOCATED(ZK4)) ALLOCATE ( ZK4(p%N) , STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = ' Error allocating memory for the ZK4 array.'
   END IF
   
   ! First call to dynamics routine:
   CALL SS_Rad_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN      
         CALL SS_Rad_DestroyContState( dxdt, ErrStat2, ErrMsg2)
         ErrMsg = TRIM(ErrMsg)//' '//TRIM(ErrMsg2)
         RETURN
      END IF  
      
   ! Compute intermediate functions to estimate next Q and QD.
   DO I = 1,p%N  ! Loop through all DOFs
      ZK1 (I) = p%DT * dxdt%x(I,1)   
   END DO          ! I - All DOFs

   ! Get first time derivatives of continuous states (dxdt):  
   CALL SS_Rad_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN      
         CALL SS_Rad_DestroyContState( dxdt, ErrStat2, ErrMsg2)
         ErrMsg = TRIM(ErrMsg)//' '//TRIM(ErrMsg2)
         RETURN
      END IF  

   ! Repeat above steps for each ZK, ZKD:
   DO I = 1,p%N  ! Loop through all DOFs
      ZK2 (I) = p%DT * dxdt%x(I,1)   
   END DO          ! I - All DOFs

   ! Get first time derivatives of continuous states (dxdt):  
      CALL SS_Rad_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN      
         CALL SS_Rad_DestroyContState( dxdt, ErrStat2, ErrMsg2)
         ErrMsg = TRIM(ErrMsg)//' '//TRIM(ErrMsg2)
         RETURN
      END IF  

   DO I = 1,p%N  ! Loop through all DOFs
      ZK3 (I) = p%DT * dxdt%x(I,1)   
   END DO          ! I - All DOFs

   ! Get first time derivatives of continuous states (dxdt):     
      CALL SS_Rad_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN      
         CALL SS_Rad_DestroyContState( dxdt, ErrStat2, ErrMsg2)
         ErrMsg = TRIM(ErrMsg)//' '//TRIM(ErrMsg2)
         RETURN
      END IF  

   ! Compute best estimate for Q, QD at next time step using
   !   the intermediate functions (Runge-Kutta).
   ! IC(NMX) locates the i + 1 value of Q, QD.
   DO I = 1,p%N  ! Loop through all DOFs
      ZK4 (I) = p%DT*dxdt%x(I,1)
      x%x  (I,1) = x%x  (I,1) + ( ZK1 (I) + 2.0*ZK2 (I) + 2.0*ZK3 (I) + ZK4 (I) ) / 6.0
   END DO          ! I - All DOFs

   IF (ALLOCATED(ZK1) ) DEALLOCATE ( ZK1  )
   IF (ALLOCATED(ZK2) ) DEALLOCATE ( ZK2  )
   IF (ALLOCATED(ZK3) ) DEALLOCATE ( ZK3  )
   IF (ALLOCATED(ZK4) ) DEALLOCATE ( ZK4  )
   
ELSE ! Use Adams-Bashforth predictor and Adams-Moulton corrector integration scheme for all other time steps.

    ! Allocate size to XT
      IF (SIZE(xt%x) /= SIZE( x%x)) THEN
        ALLOCATE ( xt%x  (size(x%x),1) , STAT=ErrStat )
        IF ( ErrStat /= 0 )  THEN
              ErrStat = ErrID_Fatal
              ErrMsg  = ' Error allocating memory for the ZK1 array.'
        END IF
    END IF

 
   ! PREDICTOR (Adams-Bashforth)
   ! Compute predictor from current (IC(1)) and 3 previous values of
   !   Q, QD.  OtherState%dxdt (I,1) = i, ...,2) = i-1, ...,3) = i-2 etc... 

   DO I = 1,p%N  ! Loop through all DOFs
      xt%x(I,1) = x%x(I,1) + p%DT/24*( 55.0*OtherState%dxdt (I,1) - 59.0*OtherState%dxdt (I,2) + 37.0*OtherState%dxdt (I,3) - 9.0*OtherState%dxdt (I,4) )
   END DO          ! I - All DOFs

   ! Get first time derivatives of continuous states (dxdt):
      
      CALL SS_Rad_CalcContStateDeriv( Time, u, p, xt, xd, z, OtherState, dxdt, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN      
         CALL SS_Rad_DestroyContState( dxdt, ErrStat2, ErrMsg2)
         ErrMsg = TRIM(ErrMsg)//' '//TRIM(ErrMsg2)
         RETURN
      END IF  

   ! CORRECTOR (Adams-Moulton)
   ! Compute corrector from predictor value of Q, QD and 3
   !   previous values of Q, QD.  OtherState%dxdt (I,1) = i, IC(2) = i-1,
   !   IC(3) = i-2 etc...

   DO I = 1,p%N  ! Loop through all DOFs
      x%x(I,1) = x%x(I,1) + p%DT/24*( 9.0*dxdt%x (I,1) + 19.0*OtherState%dxdt (I,1) - 5.0*OtherState%dxdt (I,2) + OtherState%dxdt (I,3) )
   END DO          ! I - All DOFs
END IF

END SUBROUTINE Solver
   
END MODULE SS_Radiation
!**********************************************************************************************************************************
