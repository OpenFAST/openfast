PROGRAM Main
  
  USE FEAMooring_Types
  USE FEAMooring

  USE NWTC_Library 

  IMPLICIT NONE 

  INTEGER(IntKi)                         :: ErrStat          ! Status of error message   
  CHARACTER(1024)                        :: ErrMsg           ! Error message if ErrStat /= ErrID_None
                                         
  REAL(DbKi)                             :: dt_global        ! fixed/constant global time step
  REAL(DbKi)                             :: t_initial        ! time at initialization
  REAL(DbKi)                             :: t_final          ! time at simulation end 
  REAL(DbKi)                             :: t_global         ! global-loop time marker
                                         
  INTEGER(IntKi)                         :: n_t_final        ! total number of time steps
  INTEGER(IntKi)                         :: n_t_global       ! global-loop time counter

  TYPE (FEAM_InitInputType)               :: FEAM_InitInput    
  TYPE (FEAM_ParameterType)               :: FEAM_Parameter
  TYPE (FEAM_ContinuousStateType)         :: FEAM_ContinuousState
  TYPE (FEAM_InitOutputType)              :: FEAM_InitOutput    
  TYPE (FEAM_DiscreteStateType)           :: FEAM_DiscreteState
  TYPE (FEAM_ConstraintStateType)         :: FEAM_ConstraintState
  TYPE (FEAM_OtherStateType)              :: FEAM_OtherState
  TYPE (FEAM_MiscVarType)                 :: FEAM_MiscVar

  TYPE (FEAM_InputType),      ALLOCATABLE :: FEAM_Input(:)
  REAL(DbKi) , DIMENSION(:), ALLOCATABLE :: FEAM_InputTimes

  TYPE (FEAM_OutputType)                  :: FEAM_Output

  INTEGER(IntKi)                         :: FEAM_interp_order     ! order of interpolation/extrapolation
  
  ! Local variables
  Integer(IntKi)                         :: i                    ! counter for various loops
  Integer(IntKi)                         :: j                    ! counter for various loops
  integer(intKi)                         :: Un
  
  ! -------------------------------------------------------------------------
  ! Initialization of glue-code time-step variables
  ! -------------------------------------------------------------------------
  
  t_initial = 0.
  t_final   = 50 !500
  
  ! specify time increment; currently, all modules will be time integrated with this increment size
  dt_global = 0.01
  n_t_final = ((t_final - t_initial) / dt_global ) - 1  
  t_global = t_initial  
  
  ! @bonnie : is this right? What's a good interp order?
  ! @marco: the interp order is specified in the FAST input file. It can be 0, 1, or 2
  FEAM_interp_order = 0

  ! MAP: allocate Input and Output arrays; used for interpolation and extrapolation
  Allocate(FEAM_InputTimes(FEAM_interp_order + 1)) 

  ! @bonnie : This is in the FAST developers glue code example, but it's probably not needed here. 
  Allocate(FEAM_Input(FEAM_interp_order + 1))  
    
  ! set the FEAM input file name and other environment terms.
  FEAM_InitInput%InputFile    = "FE_Mooring.dat"  ! @bonnie : This needs to be set according to what is in the FAST input file. 
  FEAM_InitInput%NStepWave   = 1                          ! an arbitrary number > 0 (to set the size of the wave data, which currently contains all zero values)     
  FEAM_InitInput%gravity     = 9.81     ! This need to be according to g used in ElastoDyn 
  FEAM_InitInput%WtrDens     = 1025     ! This needs to be set according to seawater density in HydroDyn      
  FEAM_InitInput%PtfmInit    = 0.0
  FEAM_InitInput%RootName    = "FE_Mooring.FEAM"
  
  CALL GetNewUnit( Un )
  OPEN(Unit=Un,FILE='FEAM.out',STATUS='UNKNOWN')
  !OPEN(Unit=2,FILE='FEAM_TTN1.out',STATUS='UNKNOWN')
  !OPEN(Unit=3,FILE='FEAM_TTN2.out',STATUS='UNKNOWN')
  !OPEN(Unit=4,FILE='FEAM_TTN3.out',STATUS='UNKNOWN')

  ! call the initialization routine
  CALL FEAM_Init( FEAM_InitInput      , &
                  FEAM_Input(1)        , & 
                  FEAM_Parameter       , &
                  FEAM_ContinuousState , &
                  FEAM_DiscreteState   , &
                  FEAM_ConstraintState , & 
                  FEAM_OtherState      , &
                  FEAM_Output          , &
                  FEAM_MiscVar         , &
                  dt_global            , &
                  FEAM_InitOutput      , &
                  ErrStat              , &
                  ErrMsg )  
     IF ( ErrStat .NE. ErrID_None ) THEN
        IF (ErrStat >=AbortErrLev) CALL ProgAbort(ErrMsg)
        CALL WrScr( ErrMsg )
     END IF
  
  CALL FEAM_DestroyInitInput  ( FEAM_InitInput  , ErrStat, ErrMsg )
  CALL FEAM_DestroyInitOutput ( FEAM_InitOutput , ErrStat, ErrMsg )
     
  CALL DispNVD( FEAM_InitOutput%Ver ) 

     
  FEAM_Input(1)%PtFairleadDisplacement%TranslationDisp = 0.0 
  
  DO i = 2, FEAM_interp_order + 1  
     CALL FEAM_CopyInput( FEAM_Input(1), FEAM_Input(i), MESH_NEWCOPY, ErrStat, ErrMsg )
  END DO

  DO i = 1, FEAM_interp_order + 1  
      FEAM_InputTimes(i) = t_initial - (i - 1) * dt_global
  ENDDO


  
  ! -------------------------------------------------------------------------
  ! BEGIN time marching
  ! -------------------------------------------------------------------------
   CALL FEAM_CalcOutput( t_global            , &
                        FEAM_Input(1)        , &
                        FEAM_Parameter       , &
                        FEAM_ContinuousState , &
                        FEAM_DiscreteState   , &
                        FEAM_ConstraintState , &
                        FEAM_OtherState      , &
                        FEAM_Output          , &
                        FEAM_MiscVar         , &        
                        ErrStat              , &
                        ErrMsg )
   IF ( ErrStat .NE. ErrID_None ) THEN
      IF (ErrStat >=AbortErrLev) CALL ProgAbort(ErrMsg)
      CALL WrScr( ErrMsg )
   END IF
     
  DO n_t_global = 0, n_t_final

     t_global =  t_initial + dt_global*n_t_global

     !==========   NOTE   ======     <-----------------------------------------+
     FEAM_InputTimes(1) = t_global + dt_global
     !FEAM_InputTimes(2) = FEAM_InputTimes(1) - dt_global 
     !FEAM_InputTimes(3) = FEAM_InputTimes(2) - dt_global
     
     ! Sinusoidal platform motion 
     FEAM_Input(1)%PtFairleadDisplacement%TranslationDisp(1,1) = 20*SIN(0.0005*n_t_global)  
     FEAM_Input(1)%PtFairleadDisplacement%TranslationDisp(1,2) = 20*SIN(0.0005*n_t_global)  
     FEAM_Input(1)%PtFairleadDisplacement%TranslationDisp(1,3) = 20*SIN(0.0005*n_t_global)  

     !FEAM_Input(2)%PtFairleadDisplacement%TranslationDisp(1,1) = .001*n_t_global  
     !FEAM_Input(3)%PtFairleadDisplacement%TranslationDisp(1,1) = .001*n_t_global  
     !===========================================================================


     ! @bonnie & @jason: the FAST glue code will update the new fairlead position 
     !                   based on the new platform position in the global frame.
     CALL  FEAM_UpdateStates( t_global            , &
                             n_t_global           , &
                             FEAM_Input           , &
                             FEAM_InputTimes      , &
                             FEAM_Parameter       , &
                             FEAM_ContinuousState , &
                             FEAM_DiscreteState   , &
                             FEAM_ConstraintState , &
                             FEAM_OtherState      , &
                             FEAM_MiscVar         , &        
                             ErrStat              , &
                             ErrMsg )    
     IF ( ErrStat .NE. ErrID_None ) THEN
        IF (ErrStat >=AbortErrLev) CALL ProgAbort(ErrMsg)
        CALL WrScr( ErrMsg )
     END IF
  
     CALL FEAM_CalcOutput( t_global            , &
                          FEAM_Input(1)        , &
                          FEAM_Parameter       , &
                          FEAM_ContinuousState , &
                          FEAM_DiscreteState   , &
                          FEAM_ConstraintState , &
                          FEAM_OtherState      , &
                          FEAM_Output          , &
                          FEAM_MiscVar         , &        
                          ErrStat              , &
                          ErrMsg )
     IF ( ErrStat .NE. ErrID_None ) THEN
        IF (ErrStat >=AbortErrLev) CALL ProgAbort(ErrMsg)
        CALL WrScr( ErrMsg )
     END IF
  
     ! update the global time step by one delta t
     t_global = ( n_t_global + 1 )* dt_global + t_initial

     WRITE(Un,100) t_global, FEAM_Input(1)%PtFairleadDisplacement%TranslationDisp(1,1), &
     ((FEAM_Output%PtFairleadLoad%Force(i,j), i=1,3),j=1,3)
     !WRITE(*,*) t_global
     
  END DO
  ! -------------------------------------------------------------------------
  ! END time marching
  ! -------------------------------------------------------------------------
  
  ! Destroy all objects
  CALL FEAM_End( FEAM_Input(1)       , &
                FEAM_Parameter       , &
                FEAM_ContinuousState , &
                FEAM_DiscreteState   , &
                FEAM_ConstraintState , & 
                FEAM_OtherState      , &
                FEAM_Output          , &
                FEAM_MiscVar         , &
                ErrStat              , &
                ErrMsg )  
  IF ( ErrStat .NE. ErrID_None ) THEN
     IF (ErrStat >=AbortErrLev) CALL ProgAbort(ErrMsg)
     CALL WrScr( ErrMsg )
  END IF  

  do j = 2,FEAM_interp_order+1
     call FEAM_DestroyInput( FEAM_Input(j), ErrStat, ErrMsg)
  end do  
  
  DEALLOCATE(FEAM_Input)
  DEALLOCATE(FEAM_InputTimes)
  
  CALL WrScr( "Program has ended" )
  close (un) 
  
100 FORMAT(2(1X,F8.3),9(1X,E12.5))
     
END PROGRAM Main
