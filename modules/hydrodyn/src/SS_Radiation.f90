!**********************************************************************************************************************************
! The SS_Radiation and SS_Radiation_Types modules make up a template for creating user-defined calculations in the FAST Modularization 
! Framework. SS_Radiations_Types will be auto-generated based on a description of the variables for the module.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2012, 2015  National Renewable Energy Laboratory
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
MODULE SS_Radiation

   USE SS_Radiation_Types   
   USE NWTC_Library
      
   IMPLICIT NONE
   
   PRIVATE

   TYPE(ProgDesc), PARAMETER  :: SS_Rad_ProgDesc = ProgDesc( 'SS_Radiation', '', '' )

   
      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: SS_Rad_Init                           ! Initialization routine
   PUBLIC :: SS_Rad_End                            ! Ending routine (includes clean up)
   
   PUBLIC :: SS_Rad_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating 
                                                   !   continuous states, and updating discrete states
   PUBLIC :: SS_Rad_CalcOutput                     ! Routine for computing outputs
   
   PUBLIC :: SS_Rad_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: SS_Rad_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: SS_Rad_UpdateDiscState                ! Tight coupling routine for updating discrete states
         
   
   CONTAINS
   
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine transforms  the State Space input file data from a local (heading-angle, based) coordinate system to the global system. 
!> NOTE: This routine ONLY works if all the DOFs are enabled!!!!!!!!!!
subroutine TransformStateSpaceMatrices( NBody, RotZ, B, C )
!..................................................................................................................................
   integer(IntKi), intent( in    ) :: NBody   ! Number of WAMIT bodies in this WAMIT object ( = 1 if NBodyMod > 1)
   real(R8Ki),     intent( in    ) :: RotZ(:) ! NBody heading angles (radians)
   real(ReKi),     intent( inout ) :: B(:,:)  ! Matrix data to be transformed, if NBodyMOD = 1 and NBody > 1 then we will be transforming the individual sub 6x6 matrices
   real(ReKi),     intent( inout ) :: C(:,:)  ! Matrix data to be transformed, if NBodyMOD = 1 and NBody > 1 then we will be transforming the individual sub 6x6 matrices
      
   integer(IntKi)   :: i,j,indx
   real(R8Ki)       :: R(3,3)
   real(R8Ki)       :: Rt(3,3)
      
   !do j = 1, NBody
   !   Rj(1,:) = (/ cos(RotZ(j)), sin(RotZ(j)), 0.0_R8Ki/)
   !   Rj(2,:) = (/-sin(RotZ(j)), cos(RotZ(j)), 0.0_R8Ki/)
   !   Rj(3,:) = (/ 0.0_R8Ki    , 0.0_R8Ki    , 1.0_R8Ki/)
      do i = 1, NBody
         if ( .not. EqualRealNos(RotZ(i), 0.0_R8Ki)  ) then
            R(1,:) = (/ cos(RotZ(i)), sin(RotZ(i)), 0.0_R8Ki/)
            R(2,:) = (/-sin(RotZ(i)), cos(RotZ(i)), 0.0_R8Ki/)
            R(3,:) = (/ 0.0_R8Ki    , 0.0_R8Ki    , 1.0_R8Ki/)
            Rt     = transpose(R)
            do j = 1,2  ! Need to do this twice, since a single R (3x3) matrix is used to transform all 6 DOFs associated with the ith Body data
               indx = (i-1)*6 + (j-1)*3 + 1 
               ! Create sub matrix which is all rows of B but only necessary columns for transformation work, NOTE: B is numStates X (6*NBody)
               B(:,indx:indx+2) = matmul(     B(:,indx:indx+2), R ) 
               ! Create sub matrix which is all columns of C but only necessary rows for transformation work, NOTE: c is (6*NBody) X numStates 
               C(indx:indx+2,:) = matmul( Rt, C(indx:indx+2,:)    ) 
            end do 
         end if
      end do

  ! end do
end subroutine TransformStateSpaceMatrices
   
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the start of the simulation to perform initialization steps. 
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
SUBROUTINE SS_Rad_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )
!..................................................................................................................................

    TYPE(SS_Rad_InitInputType),       INTENT(IN   )  :: InitInp     !< Input data for initialization routine
    TYPE(SS_Rad_InputType),           INTENT(  OUT)  :: u           !< An initial guess for the input; input mesh must be defined
    TYPE(SS_Rad_ParameterType),       INTENT(  OUT)  :: p           !< Parameters      
    TYPE(SS_Rad_ContinuousStateType), INTENT(  OUT)  :: x           !< Initial continuous states
    TYPE(SS_Rad_DiscreteStateType),   INTENT(  OUT)  :: xd          !< Initial discrete states
    TYPE(SS_Rad_ConstraintStateType), INTENT(  OUT)  :: z           !< Initial guess of the constraint states
    TYPE(SS_Rad_OtherStateType),      INTENT(  OUT)  :: OtherState  !< Initial other states            
    TYPE(SS_Rad_OutputType),          INTENT(  OUT)  :: y           !< Initial system outputs (outputs are not calculated; 
                                                                    !!   only the output mesh is initialized)
    TYPE(SS_Rad_MiscVarType),         INTENT(  OUT)  :: m           !< Initial misc/optimization variables            
    REAL(DbKi),                       INTENT(INOUT)  :: Interval    !< Coupling interval in seconds: the rate that 
                                                                    !!   (1) SS_Rad_UpdateStates() is called in loose coupling &
                                                                    !!   (2) SS_Rad_UpdateDiscState() is called in tight coupling.
                                                                    !!   Input is the suggested time from the glue code; 
                                                                    !!   Output is the actual coupling interval that will be used 
                                                                    !!   by the glue code.
    TYPE(SS_Rad_InitOutputType),      INTENT(  OUT)  :: InitOut     !< Output for initialization routine
    INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     !< Error status of the operation
    CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

    ! Local Variables:
         
    REAL(ReKi), ALLOCATABLE                :: Rad_A (:,:)                          ! A matrix of the radiation state-space system on the input file ss
    REAL(ReKi), ALLOCATABLE                :: Rad_B (:,:)                          ! B matrix of the radiation state-space system on the input file ss
    REAL(ReKi), ALLOCATABLE                :: Rad_C (:,:)                          ! C matrix of the radiation state-space system on the input file ss

    INTEGER                                :: I                                    ! Generic index
    INTEGER, allocatable                   :: xx (:)                               ! Active DOF's on the input file .ss
    INTEGER                                :: numDOFs                              ! Number of DOFS  
    INTEGER                                :: numStates                            ! Number of states
    integer(IntKi)                         :: N                                    ! Counter
    INTEGER                                :: Nlines                               ! Number of lines in the input file, used to determine N
    INTEGER                                :: UnSS                                 ! I/O unit number for the WAMIT output file with the .ss extension; this file contains the state-space matrices.
    INTEGER                                :: Sttus                                ! Error in reading .ss file
    character(3)                           :: bodystr
    integer                                :: ErrStat2
    character(ErrMsgLen)                   :: ErrMsg2
    
    ! Initialize ErrStat   
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
    UnSS  = -1
    N     =  0
    numStates = 0
    p%NBody = InitInp%NBody  ! Number of WAMIT bodies: =1 if WAMIT is using NBodyMod > 1,  >=1 if NBodyMod=1
    
    ! Open the .ss input file!
    CALL GetNewUnit( UnSS )
    CALL OpenFInpFile ( UnSS, TRIM(InitInp%InputFile)//'.ss', ErrStat2, ErrMsg2 )  ! Open file.
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Rad_Init')
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

    ! Determine the number of states and size of the matrices
    Nlines = 1
    
    CALL ReadCom ( UnSS, TRIM(InitInp%InputFile)//'.ss', 'Header',ErrStat2, ErrMsg2  )! Reads the first entire line (Title header)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Rad_Init')  
    call AllocAry( xx, 6*p%NBody, 'xx', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Rad_Init')
    CALL ReadAry( UnSS,TRIM(InitInp%InputFile)//'.ss', xx(:), 6*p%NBody, 'xx', 'xx vector containing the enabled dofs',ErrStat2, ErrMsg2) ! Reads in the second line, containing the active dofs vector
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Rad_Init')
    CALL ReadVar( UnSS,TRIM(InitInp%InputFile)//'.ss', numStates, 'numStates', 'Number of States',ErrStat2, ErrMsg2) ! Reads in the third line, containing the number of states
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Rad_Init')
    call AllocAry( p%spdof, 6*p%NBody, 'p%spdof', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Rad_Init')
    CALL ReadAry( UnSS,TRIM(InitInp%InputFile)//'.ss', p%spdof, 6*p%NBody, 'spdof', 'spdof vector containing the number of states per dofs',ErrStat2, ErrMsg2) ! Reads in the forth line, containing the state per dofs vector
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Rad_Init')
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF
      
    DO !Loop through all the lines of the file
        CALL ReadCom ( UnSS, TRIM(InitInp%InputFile)//'.ss', 'Header',Sttus,ErrMsg2  )! Reads the first entire line (Title header)
        IF ( Sttus == ErrID_None )  THEN ! .TRUE. when data is read in successfully                    
            Nlines=Nlines+1                    
        ELSE !We must have reach the end of the file
            EXIT
        END IF
    END DO

    ! The input file contains the matrices A [NxN], B [N x 6*NBody] and C [6*NBody x N], so
    !p%numStates = ( Nlines - 6*p%NBody ) / 2 ! this is the number of states
    
    !Verifications on the input file
    IF ( ( Nlines - 6*p%NBody ) / 2 /= numStates) THEN
      CALL SetErrStat(ErrID_Severe,'Error in the input file .ss: The size of the matrices does not correspond to the number of states!',ErrStat,ErrMsg,'SS_Rad_Init')
    END IF
    
    IF ( numStates /= SUM(p%spdof)) THEN
      CALL SetErrStat(ErrID_Severe,'Error in the input file .ss: The size of the matrices does not correspond to the number of states!',ErrStat,ErrMsg,'SS_Rad_Init')
    END IF        
    
    !Verify if the DOFs active in the input file correspond to the ones active by FAST in this run
    DO I=1,6*p%NBody !Loop through all 6 DOFs           
        IF ( InitInp%enabledDOFs (I) == 1)  THEN !  True when the current DOF is active in FAST                   
            IF ( xx (I) /= 1) THEN ! True if a DOF enabled by FAST is not available in the INPUT File
               CALL SetErrStat(ErrID_Severe,'Error in the input file .ss: The enabled DOFs in the current FAST Simulation don`t match the ones on the input file .ss!',ErrStat,ErrMsg,'SS_Rad_Init')
            END IF           
        END IF
    END DO
    
    numDOFs = SUM (xx) !Number of DOFS in the input file
    
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF
    
    ! Now we can allocate the temporary matrices A, B and C
    
    CALL AllocAry( Rad_A, numStates,    numStates,    'Rad_A', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Rad_Init')
    CALL AllocAry( Rad_B, numStates,    numDOFs, 'Rad_B', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Rad_Init')
    CALL AllocAry( Rad_C, numDOFs, numStates,    'Rad_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Rad_Init')
    
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF
    
        
    REWIND (UNIT=UnSS)   ! REWIND the file so we can read it in a second time.

    ! Skip the first 4 lines:  (NOTE: no error handling here because we would have caught it the first time through)
    CALL ReadCom ( UnSS, TRIM(InitInp%InputFile)//'.ss', 'Header', ErrStat2, ErrMsg2  )! Reads the first entire line (Title header)
    CALL ReadCom ( UnSS, TRIM(InitInp%InputFile)//'.ss', 'Enabled dofs', ErrStat2, ErrMsg2  )! Reads the first entire line (Title header)
    CALL ReadCom ( UnSS, TRIM(InitInp%InputFile)//'.ss', 'numStates', ErrStat2, ErrMsg2  )! Reads the first entire line (Title header)
    CALL ReadCom ( UnSS, TRIM(InitInp%InputFile)//'.ss', 'numStates per dofs', ErrStat2, ErrMsg2  )! Reads the first entire line (Title header)   
    
    DO I = 1,numStates !Read A MatriX
        CALL ReadAry( UnSS,TRIM(InitInp%InputFile)//'.ss', Rad_A(I,:), numStates, 'Rad_A', 'A_Matrix',ErrStat2, ErrMsg2)
          CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Rad_Init')
    END DO
    
    DO I = 1,numStates !Read B Matrix
        CALL ReadAry( UnSS, TRIM(InitInp%InputFile)//'.ss', Rad_B(I,:), 6*p%NBody, 'Rad_B', 'B_Matrix',ErrStat2, ErrMsg2) 
          CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Rad_Init')
    END DO
    
    DO I = 1,6*p%NBody !Read C Matrix
        CALL ReadAry( UnSS, TRIM(InitInp%InputFile)//'.ss', Rad_C(I,:), numStates, 'Rad_C', 'C_Matrix',ErrStat2, ErrMsg2)
          CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Rad_Init')
    END DO
    
    CLOSE ( UnSS ) !Close .ss input file
    UnSS = -1        ! Indicate the file is closed
    
    ! Transform the SS matrices using the heading angles
    ! NOTE: This transformation routine ONLY works if all the DOFs are enabled so we will do it on Rad_B and Rad_C which have all DOFs !!!!!!!!!!
         call TransformStateSpaceMatrices( p%NBody, InitInp%PtfmRefztRot, Rad_B, Rad_C )    
         
    !Now we are ready to reduce the matrices to the correspondent active dofs in FAST
    p%numStates=0
    DO I=1,6*p%NBody !For each state
        IF ( InitInp%enabledDOFs (I) == 1)  THEN !  True when the current DOF is active in FAST          
            p%numStates = p%numStates + p%spdof(I) !Add the correspondent number of states to the vector
        END IF
    END DO
    
    CALL WrScr1 ( 'Using SS_Radiation Module, with '//TRIM( Num2LStr(p%numStates ))//' of '//TRIM( Num2LStr(numStates ))// ' radiation states' )
    
    !Now we can allocate the final size of the SS matrices
    CALL AllocAry( p%A, p%numStates, p%numStates,    'p%A', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Rad_Init')
    CALL AllocAry( p%B, p%numStates, 6*p%NBody,      'p%B', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Rad_Init')
    CALL AllocAry( p%C, 6*p%NBody,   p%numStates,    'p%C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Rad_Init')
      
    ! if these arrays weren't allocated, return before a seg fault occurs:      
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF
    
                    
    !Finaly we write the ss matrices, based on the ones on the input file and on the active dofs
        
    IF ( p%numStates == numStates ) THEN !The matrices are the same
        
        p%A = Rad_A
        p%B = Rad_B
        p%C = Rad_C
        
    ELSE !We need to cut some of the lines and columns
       
        p%A = 0
        p%B = 0
        p%C = 0
       
       
        N=1 !Use as number of active states introduced
        
        DO I=1,6*p%NBody !For each dof...
            IF ( InitInp%enabledDOFs (I) == 1 .AND. sum(p%spdof(1:I))<size(Rad_A(:,1)))  THEN !  That is enabled in FAST
    
                p%A (N:N+p%spdof(I),N:N+p%spdof(I)) = Rad_A (sum(p%spdof(1:I-1))+1:sum(p%spdof(1:I)),sum(p%spdof(1:I-1))+1:sum(p%spdof(1:I)))
                p%B (N:N+p%spdof(I),:)= Rad_B (sum(p%spdof(1:I-1))+1:sum(p%spdof(1:I)),:)
                p%C (:,N:N+p%spdof(I))= Rad_C (:,sum(p%spdof(1:I-1))+1:sum(p%spdof(1:I)))
                
                N = N + p%spdof(I) !Number of lines added to the A and B Matrix and columns to the C Matrix
            END IF
        END DO
    END IF
            
    ! Define parameters here:
         
      p%DT  = Interval
       
    ! Define initial system states here:
    CALL AllocAry( x%x, p%numStates,  'x%x', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Rad_Init')      
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

   x%x = 0
     
   xd%DummyDiscState          = 0 !TD: SS doesn't have disc states
   z%DummyConstrState         = 0 !TD: SS doesn't have constr states
      
   ! Define other States: 
   DO I=1,SIZE(OtherState%xdot)
      CALL SS_Rad_CopyContState( x, OtherState%xdot(i), MESH_NEWCOPY, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Rad_Init')
   END DO
   OtherState%n = -1

! misc vars:
   m%DummyMiscVar = 0
      
   !Inputs    
   call AllocAry( u%dq, 6*p%NBody, 'u%dq', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Rad_Init')
   u%dq = 0 !All DoF's velocities

      ! Define system output initializations (set up mesh) here:
   call AllocAry( y%y,           6*p%NBody+1, 'y%y',           ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Rad_Init')
   call AllocAry( y%WriteOutput, 6*p%NBody+1, 'y%WriteOutput', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Rad_Init')
   y%y = 0         
   y%WriteOutput = 0
      
         
   ! Define initialization-routine output here:
   
   ! This output channels are only used by the stand-alone driver program and not by the OpenFAST coupled version.  
   !  For OpenFAST, these outputs are attached (via HydroDyn) to the Radiation Force/Moment channels within HydroDyn
   
   call AllocAry( InitOut%WriteOutputHdr, 6*p%NBody+1, 'InitOut%WriteOutputHdr', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Rad_Init')
   call AllocAry( InitOut%WriteOutputUnt, 6*p%NBody+1, 'InitOut%WriteOutputUnt', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SS_Rad_Init')

   InitOut%WriteOutputHdr(1) = 'Time'
   InitOut%WriteOutputUnt(1) = '(s) '
   do i = 1, p%NBody
      bodystr = 'B'//trim(num2lstr(i))
      InitOut%WriteOutputHdr( (i-1)*6+2: (i-1)*6+7 ) = (/ trim(bodystr)//'FX  ' , trim(bodystr)//'FY  ' , trim(bodystr)//'FZ  ' , trim(bodystr)//'MX  ' , trim(bodystr)//'MY  ' , trim(bodystr)//'MZ  ' /)
      InitOut%WriteOutputUnt( (i-1)*6+2: (i-1)*6+7 ) = (/ '(N) ' , '(N) ' , '(N) ' , '(Nm)' , '(Nm)' , '(Nm)' /)     
   end do   
      ! If you want to choose your own rate instead of using what the glue code suggests, tell the glue code the rate at which
      !   this module must be called here:
         
      !p%DT=Interval
   
   CALL CleanUp() ! deallocate local arrays

CONTAINS
   SUBROUTINE CleanUp()
   
      IF (UnSS > 0 ) CLOSE ( UnSS )
      
      IF ( ALLOCATED( Rad_A ) ) DEALLOCATE( Rad_A )
      IF ( ALLOCATED( Rad_B ) ) DEALLOCATE( Rad_B )
      IF ( ALLOCATED( Rad_C ) ) DEALLOCATE( Rad_C )
   
   END SUBROUTINE CleanUp
       
END SUBROUTINE SS_Rad_Init
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
SUBROUTINE SS_Rad_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(SS_Rad_InputType),           INTENT(INOUT)  :: u           !< System inputs
      TYPE(SS_Rad_ParameterType),       INTENT(INOUT)  :: p           !< Parameters     
      TYPE(SS_Rad_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states
      TYPE(SS_Rad_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Discrete states
      TYPE(SS_Rad_ConstraintStateType), INTENT(INOUT)  :: z           !< Constraint states
      TYPE(SS_Rad_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states            
      TYPE(SS_Rad_OutputType),          INTENT(INOUT)  :: y           !< System outputs
      TYPE(SS_Rad_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables            
      INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None



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
         
         ! Destroy misc vars:
      CALL SS_Rad_DestroyMisc(  m,  ErrStat, ErrMsg )
      
      
         ! Destroy the output data:
         
      CALL SS_Rad_DestroyOutput( y, ErrStat, ErrMsg )


      

END SUBROUTINE SS_Rad_End
!----------------------------------------------------------------------------------------------------------------------------------
!> Loose coupling routine for solving constraint states, integrating continuous states, and updating discrete states.
!! Continuous, constraint, and discrete states are updated to values at t + Interval.
SUBROUTINE SS_Rad_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   ) :: t               !< Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   ) :: n               !< Current step of the simulation: t = n*Interval
      TYPE(SS_Rad_InputType),             INTENT(INOUT) :: Inputs(:)       !< Inputs at InputTimes
      REAL(DbKi),                         INTENT(IN   ) :: InputTimes(:)   !< Times in seconds associated with Inputs
      TYPE(SS_Rad_ParameterType),         INTENT(IN   ) :: p               !< Parameters
      TYPE(SS_Rad_ContinuousStateType),   INTENT(INOUT) :: x               !< Input: Continuous states at t;
                                                                           !!   Output: Continuous states at t + Interval
      TYPE(SS_Rad_DiscreteStateType),     INTENT(INOUT) :: xd              !< Input: Discrete states at t;
                                                                           !!   Output: Discrete states at t + Interval
      TYPE(SS_Rad_ConstraintStateType),   INTENT(INOUT) :: z               !< Input: Constraint states at t;
                                                                           !!   Output: Constraint states at t + Interval
      TYPE(SS_Rad_OtherStateType),        INTENT(INOUT) :: OtherState      !< Input: Other states at t;
                                                                           !!   Output: Other states at t + Interval
      TYPE(SS_Rad_MiscVarType),           INTENT(INOUT) :: m               !< Initial misc/optimization variables            
      INTEGER(IntKi),                     INTENT(  OUT) :: ErrStat         !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT) :: ErrMsg          !< Error message if ErrStat /= ErrID_None

      INTEGER, PARAMETER :: IntegrationMethod = 3   
      
                                   
      SELECT CASE ( IntegrationMethod )
         
      CASE (1) ! RK4
      
         CALL SS_Rad_RK4( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
         
      CASE (2) ! AB4
      
         CALL SS_Rad_AB4( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
      
      CASE (3) ! ABM4
      
         CALL SS_Rad_ABM4( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
         
      CASE DEFAULT  !bjj: we already checked this at initialization, but for completeness:
         
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in SS_Rad_UpdateStates: method must be 1 (RK4), 2 (AB4), or 3 (ABM4)'
         RETURN
         
      END SELECT
      
     
END SUBROUTINE SS_Rad_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
SUBROUTINE SS_Rad_CalcOutput( Time, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )   
!..................................................................................................................................
   
      REAL(DbKi),                       INTENT(IN   )   :: Time        !< Current simulation time in seconds
      TYPE(SS_Rad_InputType),           INTENT(IN   )   :: u           !< Inputs at Time
      TYPE(SS_Rad_ParameterType),       INTENT(IN   )   :: p           !< Parameters
      TYPE(SS_Rad_ContinuousStateType), INTENT(IN   )   :: x           !< Continuous states at Time
      TYPE(SS_Rad_DiscreteStateType),   INTENT(IN   )   :: xd          !< Discrete states at Time
      TYPE(SS_Rad_ConstraintStateType), INTENT(IN   )   :: z           !< Constraint states at Time
      TYPE(SS_Rad_OtherStateType),      INTENT(IN   )   :: OtherState  !< Other states at Time
      TYPE(SS_Rad_OutputType),          INTENT(INOUT)   :: y           !< Outputs computed at Time (Input only so that mesh con-
                                                                       !!   nectivity information does not have to be recalculated)
      TYPE(SS_Rad_MiscVarType),         INTENT(INOUT)   :: m           !< Initial misc/optimization variables            
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
      y%WriteOutput(2:6*p%NBody+1) = y%y
                   
END SUBROUTINE SS_Rad_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for computing derivatives of continuous states
SUBROUTINE SS_Rad_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )  
!..................................................................................................................................
   
      REAL(DbKi),                        INTENT(IN   )  :: Time        !< Current simulation time in seconds
      TYPE(SS_Rad_InputType),            INTENT(IN   )  :: u           !< Inputs at Time                    
      TYPE(SS_Rad_ParameterType),        INTENT(IN   )  :: p           !< Parameters                             
      TYPE(SS_Rad_ContinuousStateType),  INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(SS_Rad_DiscreteStateType),    INTENT(IN   )  :: xd          !< Discrete states at Time
      TYPE(SS_Rad_ConstraintStateType),  INTENT(IN   )  :: z           !< Constraint states at Time
      TYPE(SS_Rad_OtherStateType),       INTENT(IN   )  :: OtherState  !< Other states                    
      TYPE(SS_Rad_MiscVarType),          INTENT(INOUT)  :: m           !< Initial misc/optimization variables            
      TYPE(SS_Rad_ContinuousStateType),  INTENT(  OUT)  :: dxdt        !< Continuous state derivatives at Time
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     !< Error status of the operation     
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
     
      CALL AllocAry( dxdt%x, p%numStates, 'SS_Rad_CalcContStateDeriv:dxdt%x', ErrStat, ErrMsg)
      IF ( ErrStat >= AbortErrLev) RETURN
            
      ! Compute the first time derivatives of the continuous states here:
      
      !Calc dxdt of a state space system
      ! [dxdt] = [A]*[xr]+B*[q]
      
      dxdt%x =matmul(p%A,x%x) + matmul( p%B, u%dq)
        
END SUBROUTINE SS_Rad_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for updating discrete states
SUBROUTINE SS_Rad_UpdateDiscState( Time, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )   
!..................................................................................................................................
   
      REAL(DbKi),                       INTENT(IN   )  :: Time        !< Current simulation time in seconds   
      TYPE(SS_Rad_InputType),           INTENT(IN   )  :: u           !< Inputs at Time                       
      TYPE(SS_Rad_ParameterType),       INTENT(IN   )  :: p           !< Parameters                                 
      TYPE(SS_Rad_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(SS_Rad_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Input: Discrete states at Time; 
                                                                      !!   Output: Discrete states at Time + Interval
      TYPE(SS_Rad_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time
      TYPE(SS_Rad_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other/optimization states           
      TYPE(SS_Rad_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables            
      INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

               
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
         ! Update discrete states here:
      
      ! StateData%DiscState = 

END SUBROUTINE SS_Rad_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for solving for the residual of the constraint state equations
SUBROUTINE SS_Rad_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, m, z_residual, ErrStat, ErrMsg )   
!..................................................................................................................................
   
      REAL(DbKi),                       INTENT(IN   )  :: Time        !< Current simulation time in seconds   
      TYPE(SS_Rad_InputType),           INTENT(IN   )  :: u           !< Inputs at Time                       
      TYPE(SS_Rad_ParameterType),       INTENT(IN   )  :: p           !< Parameters                           
      TYPE(SS_Rad_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(SS_Rad_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at Time
      TYPE(SS_Rad_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time (possibly a guess)
      TYPE(SS_Rad_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other/optimization states                    
      TYPE(SS_Rad_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables            
      TYPE(SS_Rad_ConstraintStateType), INTENT(  OUT)  :: z_residual  !< Residual of the constraint state equations using  
                                                                      !!     the input values described above      
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat    !< Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg     !< Error message if ErrStat /= ErrID_None

               
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Solve for the constraint states here:
      
      z_residual%DummyConstrState = 0

END SUBROUTINE SS_Rad_CalcConstrStateResidual
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
!!   Runge-Kutta."Sections 16.1 and 16.2 in Numerical Recipes in FORTRAN: The Art of Scientific Computing, 2nd ed. Cambridge, England: 
!!   Cambridge University Press, pp. 704-716, 1992.
!!
SUBROUTINE SS_Rad_RK4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                       INTENT(IN   )  :: t           !< Current simulation time in seconds
      INTEGER(IntKi),                   INTENT(IN   )  :: n           !< time step number
      TYPE(SS_Rad_InputType),           INTENT(INOUT)  :: u(:)        !< Inputs at t (out only for mesh record-keeping in ExtrapInterp routine)
      REAL(DbKi),                       INTENT(IN   )  :: utimes(:)   !< times of input
      TYPE(SS_Rad_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(SS_Rad_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
      TYPE(SS_Rad_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(SS_Rad_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t (possibly a guess)
      TYPE(SS_Rad_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other/optimization states
      TYPE(SS_Rad_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables            
      INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! local variables
         
      TYPE(SS_Rad_ContinuousStateType)                 :: xdot        ! time derivatives of continuous states      
      TYPE(SS_Rad_ContinuousStateType)                 :: k1          ! RK4 constant; see above
      TYPE(SS_Rad_ContinuousStateType)                 :: k2          ! RK4 constant; see above 
      TYPE(SS_Rad_ContinuousStateType)                 :: k3          ! RK4 constant; see above 
      TYPE(SS_Rad_ContinuousStateType)                 :: k4          ! RK4 constant; see above 
      TYPE(SS_Rad_ContinuousStateType)                 :: x_tmp       ! Holds temporary modification to x
      TYPE(SS_Rad_InputType)                           :: u_interp    ! interpolated value of inputs 

      INTEGER(IntKi)                                   :: ErrStat2    ! local error status
      CHARACTER(ErrMsgLen)                             :: ErrMsg2     ! local error message (ErrMsg)
      
      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      CALL SS_Rad_CopyContState( x, k1, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
      CALL SS_Rad_CopyContState( x, k2, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
      CALL SS_Rad_CopyContState( x, k3, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
      CALL SS_Rad_CopyContState( x, k4,    MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
      CALL SS_Rad_CopyContState( x, x_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN


      CALL SS_Rad_CopyInput( u(1), u_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN
                     
      ! interpolate u to find u_interp = u(t)
      CALL SS_Rad_Input_ExtrapInterp( u, utimes, u_interp, t, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      ! find xdot at t
      CALL SS_Rad_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, m, xdot, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      k1%x     = p%dt * xdot%x  
      x_tmp%x  = x%x  + 0.5 * k1%x

      ! interpolate u to find u_interp = u(t + dt/2)
      CALL SS_Rad_Input_ExtrapInterp(u, utimes, u_interp, t+0.5*p%dt, ErrStat2, ErrMsg2)
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      ! find xdot at t + dt/2
      CALL SS_Rad_CalcContStateDeriv( t + 0.5*p%dt, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      k2%x     = p%dt * xdot%x
      x_tmp%x  = x%x  + 0.5 * k2%x

      ! find xdot at t + dt/2
      CALL SS_Rad_CalcContStateDeriv( t + 0.5*p%dt, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      k3%x     = p%dt * xdot%x 
      x_tmp%x  = x%x  + k3%x

      ! interpolate u to find u_interp = u(t + dt)
      CALL SS_Rad_Input_ExtrapInterp(u, utimes, u_interp, t + p%dt, ErrStat2, ErrMsg2)
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      ! find xdot at t + dt
      CALL SS_Rad_CalcContStateDeriv( t + p%dt, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, ErrStat2, ErrMsg2 )
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
   
   
      CALL SS_Rad_DestroyContState( xdot,     ErrStat3, ErrMsg3 )
      CALL SS_Rad_DestroyContState( k1,       ErrStat3, ErrMsg3 )
      CALL SS_Rad_DestroyContState( k2,       ErrStat3, ErrMsg3 )
      CALL SS_Rad_DestroyContState( k3,       ErrStat3, ErrMsg3 )
      CALL SS_Rad_DestroyContState( k4,       ErrStat3, ErrMsg3 )
      CALL SS_Rad_DestroyContState( x_tmp,    ErrStat3, ErrMsg3 )

      CALL SS_Rad_DestroyInput(     u_interp, ErrStat3, ErrMsg3 )
         
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
         ErrMsg = TRIM(ErrMsg)//'SS_Rad_RK4:'//TRIM(Msg)         
         ErrStat = MAX(ErrStat,ErrID)
         
         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         
         IF ( ErrStat >= AbortErrLev ) CALL ExitThisRoutine( )                  
                  
         
      END IF

   END SUBROUTINE CheckError                    
      
END SUBROUTINE SS_Rad_RK4
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
SUBROUTINE SS_Rad_AB4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   )  :: t           !< Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   )  :: n           !< time step number
      TYPE(SS_Rad_InputType),             INTENT(INOUT)  :: u(:)        !< Inputs at t (out only for mesh record-keeping in ExtrapInterp routine)
      REAL(DbKi),                         INTENT(IN   )  :: utimes(:)   !< times of input
      TYPE(SS_Rad_ParameterType),         INTENT(IN   )  :: p           !< Parameters
      TYPE(SS_Rad_ContinuousStateType),   INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
      TYPE(SS_Rad_DiscreteStateType),     INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(SS_Rad_ConstraintStateType),   INTENT(IN   )  :: z           !< Constraint states at t (possibly a guess)
      TYPE(SS_Rad_OtherStateType),        INTENT(INOUT)  :: OtherState  !< Other/optimization states
      TYPE(SS_Rad_MiscVarType),           INTENT(INOUT)  :: m           !< Initial misc/optimization variables            
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


      ! local variables
      TYPE(SS_Rad_InputType)                             :: u_interp
         
      INTEGER(IntKi)                                     :: ErrStat2    ! local error status
      CHARACTER(ErrMsgLen)                               :: ErrMsg2     ! local error message (ErrMsg)

      
      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 
      
      
      if (OtherState%n .lt. n) then

         OtherState%n = n
                    
         CALL SS_Rad_CopyContState( OtherState%xdot ( 3 ), OtherState%xdot ( 4 ), MESH_UPDATECOPY, ErrStat2, ErrMsg )
            CALL CheckError(ErrStat2,ErrMsg2)
            IF ( ErrStat >= AbortErrLev ) RETURN
         CALL SS_Rad_CopyContState( OtherState%xdot ( 2 ), OtherState%xdot ( 3 ), MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            CALL CheckError(ErrStat2,ErrMsg2)
            IF ( ErrStat >= AbortErrLev ) RETURN
         CALL SS_Rad_CopyContState( OtherState%xdot ( 1 ), OtherState%xdot ( 2 ), MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            CALL CheckError(ErrStat2,ErrMsg2)
            IF ( ErrStat >= AbortErrLev ) RETURN
         
         
            
      elseif (OtherState%n .gt. n) then
 
         CALL CheckError( ErrID_Fatal, ' Backing up in time is not supported with a multistep method.')
         RETURN

      endif        
      
      
      ! Allocate the input arrays
      CALL SS_Rad_CopyInput( u(1), u_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      
      ! need xdot at t
      CALL SS_Rad_Input_ExtrapInterp(u, utimes, u_interp, t, ErrStat2, ErrMsg2)
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN
         
      CALL SS_Rad_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, m, OtherState%xdot ( 1 ), ErrStat2, ErrMsg2 ) ! initializes OtherState%xdot ( 1 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

                                                    
      if (n .le. 2) then
                                               
         CALL SS_Rad_RK4(t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat2, ErrMsg2 )
            CALL CheckError(ErrStat2,ErrMsg2)
            IF ( ErrStat >= AbortErrLev ) RETURN

      else

         x%x  = x%x  + p%DT/24.0  * ( 55.*OtherState%xdot(1)%x  - 59.*OtherState%xdot(2)%x   &
                                    + 37.*OtherState%xdot(3)%x   - 9.*OtherState%xdot(4)%x )


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
   
   
      CALL SS_Rad_DestroyInput(     u_interp, ErrStat3, ErrMsg3 )
         
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
         ErrMsg = TRIM(ErrMsg)//'SS_Rad_AB4:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         
         IF ( ErrStat >= AbortErrLev ) CALL ExitThisRoutine( )                  
         
      END IF

   END SUBROUTINE CheckError            
         
END SUBROUTINE SS_Rad_AB4
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
SUBROUTINE SS_Rad_ABM4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   )  :: t           !< Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   )  :: n           !< time step number
      TYPE(SS_Rad_InputType),             INTENT(INOUT)  :: u(:)        !< Inputs at t (out only for mesh record-keeping in ExtrapInterp routine)
      REAL(DbKi),                         INTENT(IN   )  :: utimes(:)   !< times of input
      TYPE(SS_Rad_ParameterType),         INTENT(IN   )  :: p           !< Parameters
      TYPE(SS_Rad_ContinuousStateType),   INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
      TYPE(SS_Rad_DiscreteStateType),     INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(SS_Rad_ConstraintStateType),   INTENT(IN   )  :: z           !< Constraint states at t (possibly a guess)
      TYPE(SS_Rad_OtherStateType),        INTENT(INOUT)  :: OtherState  !< Other/optimization states
      TYPE(SS_Rad_MiscVarType),           INTENT(INOUT)  :: m           !< Initial misc/optimization variables            
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! local variables

      TYPE(SS_Rad_InputType)                             :: u_interp    ! Inputs at t
      TYPE(SS_Rad_ContinuousStateType)                   :: x_pred      ! Continuous states at t
      TYPE(SS_Rad_ContinuousStateType)                   :: xdot_pred   ! Derivative of continuous states at t

      INTEGER(IntKi)                                     :: ErrStat2    ! local error status
      CHARACTER(ErrMsgLen)                               :: ErrMsg2     ! local error message (ErrMsg)
      
      
      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      CALL SS_Rad_CopyContState(x, x_pred, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      CALL SS_Rad_AB4( t, n, u, utimes, p, x_pred, xd, z, OtherState, m, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      if (n .gt. 2_IntKi) then
            ! allocate the arrays in u_interp
         CALL SS_Rad_CopyInput( u(1), u_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
            CALL CheckError(ErrStat2,ErrMsg2)
            IF ( ErrStat >= AbortErrLev ) RETURN
         
         CALL SS_Rad_Input_ExtrapInterp(u, utimes, u_interp, t + p%dt, ErrStat2, ErrMsg2)
            CALL CheckError(ErrStat2,ErrMsg2)
            IF ( ErrStat >= AbortErrLev ) RETURN

         CALL SS_Rad_CalcContStateDeriv(t + p%dt, u_interp, p, x_pred, xd, z, OtherState, m, xdot_pred, ErrStat2, ErrMsg2 )
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
   
   
      CALL SS_Rad_DestroyContState( xdot_pred,  ErrStat3, ErrMsg3 )
      CALL SS_Rad_DestroyContState( x_pred,     ErrStat3, ErrMsg3 )
      CALL SS_Rad_DestroyInput(     u_interp,   ErrStat3, ErrMsg3 )               
      
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
         ErrMsg = TRIM(ErrMsg)//'SS_Rad_ABM4:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) CALL ExitThisRoutine( )                  
         
      END IF

   END SUBROUTINE CheckError                 

END SUBROUTINE SS_Rad_ABM4
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE SS_Radiation
!**********************************************************************************************************************************
