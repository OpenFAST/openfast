PROGRAM Main

   ! A Driver Program for the MoorDyn module (first draft)

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
   REAL(ReKi), ALLOCATABLE               :: PtfmMotIn(:,:)   ! Variable for storing time, and 6DOF platform position from motion input file
   REAL(ReKi), ALLOCATABLE               :: PtfmMot(:,:)     ! Variable for storing interpolated 6DOF platform position
   INTEGER(IntKi)                        :: ntIn             ! number of time steps read from motion file
   INTEGER(IntKi)                        :: nt               ! number of coupling time steps to use in simulation

   REAL(DbKi)                            :: dtC              ! fixed/constant global time step
                                         
   INTEGER(IntKi)                        :: MD_interp_order     ! order of interpolation/extrapolation

   ! Local variables
   Integer(IntKi)                        :: i                    ! counter for various loops
   Integer(IntKi)                        :: j                    ! counter for various loops
   integer(intKi)                        :: Un

  
  
   ! -------------------------------------------------------------------------
   ! Read in prescribed motions from text file if available 
   ! (single 6DOF platform for now, to be updated for versatile coupling in future)
   ! -------------------------------------------------------------------------
   
   CALL GetNewUnit( UnPtfmMotIn ) 
   
   CALL OpenFInpFile ( UnPtfmMotIn, "PtfmMotion.txt", ErrStat, ErrMsg ) 
   
   IF (ErrStat == 0 ) THEN
         
      ! Read through length of file to find its length
      i = 1  ! start counter
      DO
         READ(UnPtfmMotIn,'(A)',IOSTAT=ErrStat) Line      !read into a line
         IF (ErrStat2 > 0) EXIT            
         i = i+1
      END DO

      ! rewind to start of input file to re-read things now that we know how long it is
      REWIND(UnPtfmMotIn)      

      ntIn = i     ! save number of lines of file

      ! allocate space for input motion array
      ALLOCATE ( PtfmMotIn(ntIn, 7), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = '  Error allocating space for PtfmMotIn array.'
         CALL WrScr( ErrMsg )
      END IF 

      ! read the data in from the file
      DO i = 1, ntIn
         READ (UnPtfmMotIn,*,IOSTAT=ErrStat) (PtfmMotIn (i,J), J=1,7)
            
         IF ( ErrStat /= 0 ) THEN
            ErrMsg = '  Error reading the input time-series file. '
            CALL WrScr( ErrMsg )
         END IF 
      END DO  

         ! Close the inputs file 
      CLOSE ( UnPtfmMotIn ) 
   
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

   ! allocate space for processed motion array
   ALLOCATE ( PtfmMot(nt, 6), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = '  Error allocating space for PtfmMot array.'
      CALL WrScr( ErrMsg )
   END IF 

   ! go through and interpolate inputs to new regular time steps (if nt=0 this array should be left as zeros)
   DO i = 0,nt
      
      t = i*dtC
      
      ! interpolation routine
      DO iIn = 1,ntIn-1      
         IF (PtfmMotIn(i+1, 1) > t) THEN
            frac = (t - PtfmMotIn(iIn, 1) )/( PtfmMotIn(iIn+1, 1) - PtfmMotIn(iIn, 1) )
            
            DO J=1,6
               PtfmMot(i, J) = PtfmMotIn(iIn, J+1) + frac*(PtfmMotIn(iIn+1, J+1) - PtfmMotIn(iIn, J+1))
            END DO
         END IF
      END DO
   END DO
   

   ! -------------------------------------------------------------------------
   ! Initialize simulation
   ! -------------------------------------------------------------------------

  
  
  
  MD_interp_order = 0

  ! MAP: allocate Input and Output arrays; used for interpolation and extrapolation
  Allocate(MD_InputTimes(MD_interp_order + 1)) 

  ! @bonnie : This is in the FAST developers glue code example, but it's probably not needed here. 
  Allocate(MD_Input(MD_interp_order + 1))
    
  ! set the input file name and other environment terms.
  MD_InitInput%InputFile    = "MoorDyn.dat" 
  MD_InitInput%NStepWave   = 1        ! an arbitrary number > 0 (to set the size of the wave data, which currently contains all zero values)     
  MD_InitInput%gravity     = 9.81     ! This need to be according to g used in ElastoDyn 
  MD_InitInput%WtrDens     = 1025     ! This needs to be set according to seawater density in HydroDyn      
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

     
  MD_Input(1)%PtFairleadDisplacement%TranslationDisp = 0.0   !@mhall what's this????  <<< zeros the whole array?
  
  DO i = 2, MD_interp_order + 1  
     CALL MD_CopyInput( MD_Input(1), MD_Input(i), MESH_NEWCOPY, ErrStat, ErrMsg )
  END DO

  DO i = 1, MD_interp_order + 1  
      MD_InputTimes(i) = -(i - 1) * dtC
  ENDDO

   CALL MD_CalcOutput( t_global            , &
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

  DO i = 1,nt

     t = dtC*(i-1)

     !==========   NOTE   ======     <-----------------------------------------+
     MD_InputTimes(1) = t + dtC
     !MD_InputTimes(2) = MD_InputTimes(1) - dtC 
     !MD_InputTimes(3) = MD_InputTimes(2) - dtC
     
     ! apply platform translations (neglecting rotations for now)
     MD_Input(1)%PtFairleadDisplacement%TranslationDisp(1,1) = PtfmMot(i, 1))  
     MD_Input(1)%PtFairleadDisplacement%TranslationDisp(1,2) = PtfmMot(i, 2))  
     MD_Input(1)%PtFairleadDisplacement%TranslationDisp(1,3) = PtfmMot(i, 3))  

     !MD_Input(2)%PtFairleadDisplacement%TranslationDisp(1,1) = .001*n_t_global  
     !MD_Input(3)%PtFairleadDisplacement%TranslationDisp(1,1) = .001*n_t_global  
     !===========================================================================


     ! @bonnie & @jason: the FAST glue code will update the new fairlead position 
     !                   based on the new platform position in the global frame.
     CALL  MD_UpdateStates( t            , &
                             n_t_global           , &
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
     END IF
  
     CALL MD_CalcOutput( t            , &
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
  
     ! update the global time step by one delta t               <<<< ??? why?
     t = t + dtC

     WRITE(Un,100) t, MD_Input(1)%PtFairleadDisplacement%TranslationDisp(1,1), &
     ((MD_Output%PtFairleadLoad%Force(i,j), i=1,3),j=1,3)
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
  
  DEALLOCATE(PtfmMot)
  DEALLOCATE(PtfmMotIn)
  
  CALL WrScr( "Program has ended" )
  close (un) 
  
100 FORMAT(2(1X,F8.3),9(1X,E12.5))
     
END PROGRAM Main
