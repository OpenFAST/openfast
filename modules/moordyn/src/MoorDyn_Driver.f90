!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2020 National Renewable Energy Laboratory
! Copyright (C) 2020 Matthew Hall
!
!    This file is part of MoorDyn.
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
PROGRAM MoorDyn_Driver

   USE MoorDyn_Types
   USE MoorDyn
   USE NWTC_Library 

   IMPLICIT NONE 

   INTEGER(IntKi)                        :: ErrStat          ! Status of error message   
   CHARACTER(1024)                       :: ErrMsg           ! Error message if ErrStat /= ErrID_None

   TYPE (MD_InitInputType)               :: MD_InitInput    
   TYPE (MD_ParameterType)               :: MD_Parameter
   TYPE (MD_ContinuousStateType)         :: MD_ContinuousState
   TYPE (MD_InitOutputType)              :: MD_InitOutput    
   TYPE (MD_DiscreteStateType)           :: MD_DiscreteState
   TYPE (MD_ConstraintStateType)         :: MD_ConstraintState
   TYPE (MD_OtherStateType)              :: MD_OtherState
   TYPE (MD_MiscVarType)                 :: MD_MiscVar

   TYPE (MD_InputType),      ALLOCATABLE :: MD_Input(:)
   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: MD_InputTimes

   TYPE (MD_OutputType)                  :: MD_Output        ! Output file identifier

   INTEGER(IntKi)                        :: UnPtfmMotIn      ! platform motion input file identifier
   CHARACTER(100)                        :: Line             ! String to temporarially hold value of read line
   REAL(ReKi), ALLOCATABLE               :: PtfmMotIn(:,:)   ! Variable for storing time, and 6DOF platform position from motion input file
   REAL(ReKi), ALLOCATABLE               :: PtfmMot(:,:)     ! Variable for storing interpolated 6DOF platform position
   INTEGER(IntKi)                        :: ntIn             ! number of time steps read from motion file
   INTEGER(IntKi)                        :: nt               ! number of coupling time steps to use in simulation

   REAL(DbKi)                            :: t                ! current time (s)
   REAL(DbKi)                            :: tMax             ! sim end time (s)
   REAL(DbKi)                            :: dtC              ! fixed/constant global time step
   REAL(DbKi)                            :: frac             ! fraction used in interpolation
         
   INTEGER(IntKi)                        :: MD_interp_order     ! order of interpolation/extrapolation

   ! Local variables
   Integer(IntKi)                        :: i                    ! counter for various loops
   Integer(IntKi)                        :: j                    ! counter for various loops
   Integer(IntKi)                        :: k                    ! counter for various loops
   Integer(IntKi)                        :: iIn
   integer(intKi)                        :: Un
   
  
   ! -------------------------------------------------------------------------
   ! Read in prescribed motions from text file if available 
   ! (single 6DOF platform for now, plus one active tensioning command)
   ! (to be updated for versatile coupling in future)
   ! -------------------------------------------------------------------------
   
   CALL GetNewUnit( UnPtfmMotIn ) 
   
   CALL OpenFInpFile ( UnPtfmMotIn, "PtfmMotion.txt", ErrStat, ErrMsg ) 
   
   IF (ErrStat == 0 ) THEN
   
      print *, "Reading platform motion input data from PtfmMotion.txt"
         
      ! Read through length of file to find its length
      i = 1  ! start counter
      DO
         READ(UnPtfmMotIn,'(A)',IOSTAT=ErrStat) Line      !read into a line
         
         
         IF (ErrStat /= 0) EXIT            
         
         print *, TRIM(Line)
         
         i = i+1
      END DO

      ! rewind to start of input file to re-read things now that we know how long it is
      REWIND(UnPtfmMotIn)      

      ntIn = i-1     ! save number of lines of file
      

      ! allocate space for input motion array
      ALLOCATE ( PtfmMotIn(ntIn, 8), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = '  Error allocating space for PtfmMotIn array.'
         CALL WrScr( ErrMsg )
      END IF 

      ! read the data in from the file
      DO i = 1, ntIn
         READ (UnPtfmMotIn,*,IOSTAT=ErrStat) (PtfmMotIn (i,J), J=1,8)
            
         IF ( ErrStat /= 0 ) THEN
            ErrMsg = '  Error reading the input time-series file. '
            CALL WrScr( ErrMsg )
         END IF 
      END DO  

         ! Close the inputs file 
      CLOSE ( UnPtfmMotIn ) 
      
      print *, "read ", ntIn, " time steps from input file."
      print *, PtfmMotIn
   
   ELSE
      ntIn = 0   ! flag to indicate no motion input file
   END IF

  
   ! ----------------------- specify stepping details -----------------------

   dtC = 0.01                   ! desired coupling time step size for communicating with MoorDyn

   IF (ntIn > 0) THEN
      tMax = PtfmMotIn(ntIn, 1)    ! save last time step as total sim time
   ELSE
      tMax = 60
   END IF
   

   nt = tMax/dtC - 1            ! number of coupling time steps

   CALL WrScr(" ")
   print *, "Tmax - ", tMax, " and nt=", nt
   CALL WrScr(" ")
   
   ! allocate space for processed motion array
   ALLOCATE ( PtfmMot(nt, 7), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = '  Error allocating space for PtfmMot array.'
      CALL WrScr( ErrMsg )
   END IF 


   ! go through and interpolate inputs to new regular time steps (if nt=0 this array should be left as zeros)
   IF (ntIn > 0) THEN
      DO i = 1,nt
         
         t = dtC*(i-1)
         
         ! interpolation routine
         DO iIn = 1,ntIn-1      
            IF (PtfmMotIn(iIn+1, 1) > t) THEN
               frac = (t - PtfmMotIn(iIn, 1) )/( PtfmMotIn(iIn+1, 1) - PtfmMotIn(iIn, 1) )
               
 !              print *, "t=", t, ", iIn=", iIn, ", frac=", frac
               
               DO J=1,7
                  PtfmMot(i, J) = PtfmMotIn(iIn, J+1) + frac*(PtfmMotIn(iIn+1, J+1) - PtfmMotIn(iIn, J+1))
               END DO
               
               EXIT
            END IF
         END DO
      
 !        print *, t, "s", PtfmMot(i,:)
      
      END DO
      
      
   ELSE 
      PtfmMot = 0.0_Reki
   END IF

   ! -------------------------------------------------------------------------
   ! Initialize simulation
   ! -------------------------------------------------------------------------

  
  
  
   MD_interp_order = 0
  
   ! MAP: allocate Input and Output arrays; used for interpolation and extrapolation
   Allocate(MD_InputTimes(MD_interp_order + 1)) 
  
   ! @bonnie : This is in the FAST developers glue code example, but it's probably not needed here. 
   Allocate(MD_Input(MD_interp_order + 1))
     
   ! set the input file name and other environment terms.
   MD_InitInput%FileName    = "MoorDyn.dat" 
   !MD_InitInput%NStepWave   = 1        ! an arbitrary number > 0 (to set the size of the wave data, which currently contains all zero values)     
   MD_InitInput%g           = 9.81     ! This need to be according to g used in ElastoDyn 
   MD_InitInput%rhoW        = 1025     ! This needs to be set according to seawater density in HydroDyn      
   MD_InitInput%PtfmInit    = 0.0
   MD_InitInput%RootName    = "MoorDyn.MD"
   
   CALL GetNewUnit( Un )
   OPEN(Unit=Un,FILE='MD.out',STATUS='UNKNOWN')
  
   ! call the initialization routine
   CALL MD_Init( MD_InitInput      , &
                   MD_Input(1)        , & 
                   MD_Parameter       , &
                   MD_ContinuousState , &
                   MD_DiscreteState   , &
                   MD_ConstraintState , & 
                   MD_OtherState      , &
                   MD_Output          , &
                   MD_MiscVar         , &
                   dtC                , &
                   MD_InitOutput      , &
                   ErrStat            , &
                   ErrMsg )  
      IF ( ErrStat .NE. ErrID_None ) THEN
         IF (ErrStat >=AbortErrLev) CALL ProgAbort(ErrMsg)
         CALL WrScr( ErrMsg )
      END IF
   
   CALL MD_DestroyInitInput  ( MD_InitInput  , ErrStat, ErrMsg )
   CALL MD_DestroyInitOutput ( MD_InitOutput , ErrStat, ErrMsg )
      
   CALL DispNVD( MD_InitOutput%Ver ) 
  
   ! set the initial input values
   MD_Input(1)%PtFairleadDisplacement%TranslationDisp = 0.0_ReKi   ! zero the displacements to start with
   MD_Input(1)%DeltaL = 0.0_ReKi
   MD_Input(1)%DeltaLdot = 0.0_ReKi
   
   IF (size(MD_Input(1)%DeltaL) .NE. 2) then
      PRINT *, "Error: this MoorDyn Driver expects there to be two active tensioning channels"
      STOP
   END IF
   
   DO i = 2, MD_interp_order + 1  
      CALL MD_CopyInput( MD_Input(1), MD_Input(i), MESH_NEWCOPY, ErrStat, ErrMsg )
   END DO
  
   DO i = 1, MD_interp_order + 1  
       MD_InputTimes(i) = -(i - 1) * dtC
   ENDDO


   t = 0

   CALL MD_CalcOutput(  t                  , &
                        MD_Input(1)        , &
                        MD_Parameter       , &
                        MD_ContinuousState , &
                        MD_DiscreteState   , &
                        MD_ConstraintState , &
                        MD_OtherState      , &
                        MD_Output          , &
                        MD_MiscVar         , &        
                        ErrStat              , &
                        ErrMsg )
   IF ( ErrStat .NE. ErrID_None ) THEN
      IF (ErrStat >=AbortErrLev) CALL ProgAbort(ErrMsg)
      CALL WrScr( ErrMsg )
   END IF
  
  
  ! -------------------------------------------------------------------------
  ! BEGIN time marching
  ! -------------------------------------------------------------------------

   print *,"Doing time marching now..."

   DO i = 1,nt

      t = dtC*(i-1)

      !==========   NOTE   ======     <-----------------------------------------+
      MD_InputTimes(1) = t + dtC
      !MD_InputTimes(2) = MD_InputTimes(1) - dtC 
      !MD_InputTimes(3) = MD_InputTimes(2) - dtC

      ! apply platform translations (neglecting rotations for now)
      MD_Input(1)%PtFairleadDisplacement%TranslationDisp(1,1) = PtfmMot(i, 1)  
      MD_Input(1)%PtFairleadDisplacement%TranslationDisp(1,2) = PtfmMot(i, 2)  
      MD_Input(1)%PtFairleadDisplacement%TranslationDisp(1,3) = PtfmMot(i, 3)  

      !MD_Input(2)%PtFairleadDisplacement%TranslationDisp(1,1) = .001*n_t_global  
      !MD_Input(3)%PtFairleadDisplacement%TranslationDisp(1,1) = .001*n_t_global  
      
      ! what about velocities??
      
      ! also provide the active tensioning command on two channels (second negative) for testing
      MD_Input(1)%DeltaL(1) =  PtfmMot(i, 7)  
      MD_Input(1)%DeltaL(2) = -PtfmMot(i, 7)  
      ! and the rate of change
      IF (i>1) then
         MD_Input(1)%DeltaLdot(1) = (PtfmMot(i, 7) - PtfmMot(i-1, 7))/dtC
         MD_Input(1)%DeltaLdot(2) = -MD_Input(1)%DeltaLdot(1)
      ELSE
         MD_Input(1)%DeltaLdot(1) = 0.0_ReKi
         MD_Input(1)%DeltaLdot(2) = 0.0_ReKi
      END IF
      
      !===========================================================================


      ! @bonnie & @jason: the FAST glue code will update the new fairlead position 
      !                   based on the new platform position in the global frame.
      CALL  MD_UpdateStates(  t                  , &
                             nt                 , &
                             MD_Input           , &
                             MD_InputTimes      , &
                             MD_Parameter       , &
                             MD_ContinuousState , &
                             MD_DiscreteState   , &
                             MD_ConstraintState , &
                             MD_OtherState      , &
                             MD_MiscVar         , &        
                             ErrStat              , &
                             ErrMsg )    
      IF ( ErrStat .NE. ErrID_None ) THEN
         IF (ErrStat >=AbortErrLev) CALL ProgAbort(ErrMsg)
         CALL WrScr( ErrMsg )
         EXIT
      END IF
  
      ! update the global time step by one delta t               <<<< ??? why?
      t = t + dtC
     
      CALL MD_CalcOutput(  t                  , &
                          MD_Input(1)        , &
                          MD_Parameter       , &
                          MD_ContinuousState , &
                          MD_DiscreteState   , &
                          MD_ConstraintState , &
                          MD_OtherState      , &
                          MD_Output          , &
                          MD_MiscVar         , &        
                          ErrStat              , &
                          ErrMsg )
      IF ( ErrStat .NE. ErrID_None ) THEN
         IF (ErrStat >=AbortErrLev) CALL ProgAbort(ErrMsg)
         CALL WrScr( ErrMsg )
      END IF
  

      WRITE(Un,100) t, MD_Input(1)%PtFairleadDisplacement%TranslationDisp(1,1), &
      ((MD_Output%PtFairleadLoad%Force(k,j), k=1,3),j=1,3)
      !WRITE(*,*) t_global
     
   END DO
   ! -------------------------------------------------------------------------
   ! END time marching
   ! -------------------------------------------------------------------------
   
   ! Destroy all objects
   CALL MD_End( MD_Input(1)       , &
                 MD_Parameter       , &
                 MD_ContinuousState , &
                 MD_DiscreteState   , &
                 MD_ConstraintState , & 
                 MD_OtherState      , &
                 MD_Output          , &
                 MD_MiscVar         , &
                 ErrStat              , &
                 ErrMsg )  
   IF ( ErrStat .NE. ErrID_None ) THEN
      IF (ErrStat >=AbortErrLev) CALL ProgAbort(ErrMsg)
      CALL WrScr( ErrMsg )
   END IF  
   
   do j = 2,MD_interp_order+1
      call MD_DestroyInput( MD_Input(j), ErrStat, ErrMsg)
   end do  
   
   DEALLOCATE(MD_Input)
   DEALLOCATE(MD_InputTimes)
   
   IF (ALLOCATED(PtfmMot)  ) DEALLOCATE(PtfmMot)
   IF (ALLOCATED(PtfmMotIn)) DEALLOCATE(PtfmMotIn)
   
   CALL WrScr( "Program has ended" )
   close (un) 
  
100 FORMAT(2(1X,F8.3),9(1X,E12.5))
     
 END PROGRAM