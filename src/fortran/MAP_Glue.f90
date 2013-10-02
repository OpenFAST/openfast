PROGRAM Main
  
  USE MAP_Types
  USE MAP

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

  TYPE (MAP_InitInputType)               :: MAP_InitInput    
  TYPE (MAP_ParameterType)               :: MAP_Parameter
  TYPE (MAP_ContinuousStateType)         :: MAP_ContinuousState
  TYPE (MAP_ContinuousStateType)         :: MAP_ContinuousStateDeriv
  TYPE (MAP_InitOutputType)              :: MAP_InitOutput    
  TYPE (MAP_DiscreteStateType)           :: MAP_DiscreteState
  TYPE (MAP_ConstraintStateType)         :: MAP_ConstraintState
  TYPE (MAP_OtherStateType)              :: MAP_OtherState

  TYPE (MAP_InputType),      ALLOCATABLE :: MAP_Input(:)
  REAL(DbKi) , DIMENSION(:), ALLOCATABLE :: MAP_InputTimes

  TYPE (MAP_OutputType)                  :: MAP_Output
  REAL(DbKi) , DIMENSION(:), ALLOCATABLE :: MAP_OutputTimes

  INTEGER(IntKi)                         :: MAP_interp_order     ! order of interpolation/extrapolation
  
  ! Local variables
  Integer(IntKi)                         :: i                    ! counter for various loops
  
  
  ! -------------------------------------------------------------------------
  ! Initialization of glue-code time-step variables
  ! -------------------------------------------------------------------------
  
  t_initial = 0.
  t_final   = 5.0
  
  ! specify time increment; currently, all modules will be time integrated with this increment size
  dt_global = 0.5
  n_t_final = ((t_final - t_initial) / dt_global ) - 1  
  t_global = t_initial  
  
  ! @bonnie : is this right? What's a good interp order?
  ! @marco: the interp order is specified in the FAST input file. It can be 0, 1, or 2
  MAP_interp_order = 0

  ! MAP: allocate Input and Output arrays; used for interpolation and extrapolation
  Allocate(MAP_OutputTimes(MAP_interp_order + 1)) 
  Allocate(MAP_InputTimes(MAP_interp_order + 1)) 

  ! @bonnie : This is in the FAST developers glue code example, but it's probably not needed here. 
  Allocate(MAP_Input(MAP_interp_order + 1))  
    
  ! set the MAP input file name and other environment terms.
  MAP_InitInput%filename    = "NRELOffshrBsline5MW_Platform_ITIBarge4.map"  ! @bonnie : This needs to be set according to what is in the FAST input file. 
  MAP_InitInput%gravity     = 9.81          ! @bonnie : This need to be according to g used in FAST
  MAP_InitInput%sea_density = 1025          ! @bonnie : This needs to be set according to seawater density in FAST
  MAP_InitInput%depth       = 150           ! @bonnie : This need to be set according to the water depth in FAST
  MAP_Parameter%dt = dt_global              ! @bonnie : This is for the glue code to set

  ! call the initialization routine
  CALL MAP_Init( MAP_InitInput       , &
                 MAP_Input(1)        , & 
                 MAP_Parameter       , &
                 MAP_ContinuousState , &
                 MAP_DiscreteState   , &
                 MAP_ConstraintState , & 
                 MAP_OtherState      , &
                 MAP_Output          , &
                 dt_global           , &
                 MAP_InitOutput      , &
                 ErrStat             , &
                 ErrMsg )  
  IF ( ErrStat .NE. 0 ) THEN
     CALL WrScr(ErrMsg) 
  END IF  
  
  CALL DispNVD( MAP_InitOutput%Ver ) 

  DO i = 1, MAP_interp_order + 1  
      MAP_InputTimes(i) = t_initial - (i - 1) * dt_global
      MAP_OutputTimes(i) = t_initial - (i - 1) * dt_global
   ENDDO
   
  DO i = 2, MAP_interp_order + 1  
     CALL MAP_CopyInput( MAP_Input(1), MAP_Input(i), MESH_NEWCOPY, ErrStat, ErrMsg )
  END DO


  ! should probably delete this at the end bc save/retrieve will need it
  !................
  !@marco, some questions:
  ! 1) why will save/retrieve need this? Unless you're using the data type later, I don't know why we'd pack the InitInput/Output data.
  ! 2) I don't want to have to call MAP_InitInput_Destroy in my glue code; I don't want the glue code to act any differently for C or Fortran code.
  !    Can we call MAP_InitInput_Destroy from MAP_DestroyInitInput?
  ! 3) Same as (2) for MAP_InitOutput_Destroy.
  !................
!  CALL MAP_InitInput_Destroy ( MAP_InitInput%C_obj%object )  
  CALL MAP_DestroyInitInput  ( MAP_InitInput , ErrStat, ErrMsg )

!  CALL MAP_InitOutput_Destroy( MAP_InitOutput%C_obj%object )  
  CALL MAP_DestroyInitOutput ( MAP_InitOutput , ErrStat, ErrMsg )

  ! @bonnie : don't we need to initialize the messhes once?
  IF (MAP_interp_order .EQ. 2) THEN
     DO i = 1,MAP_Input(3)%PtFairleadDisplacement%NNodes
        MAP_Input(3)%PtFairleadDisplacement%TranslationDisp(1,i) = 0.0  
        MAP_Input(3)%PtFairleadDisplacement%TranslationDisp(2,i) = 0.0
        MAP_Input(3)%PtFairleadDisplacement%TranslationDisp(3,i) = 0.0
     END DO
     DO i = 1,MAP_Input(2)%PtFairleadDisplacement%NNodes
        MAP_Input(2)%PtFairleadDisplacement%TranslationDisp(1,i) = 0.0  
        MAP_Input(2)%PtFairleadDisplacement%TranslationDisp(2,i) = 0.0
        MAP_Input(2)%PtFairleadDisplacement%TranslationDisp(3,i) = 0.0
     END DO
  ELSE IF (MAP_interp_order .EQ. 1) THEN
     DO i = 1,MAP_Input(2)%PtFairleadDisplacement%NNodes
        MAP_Input(2)%PtFairleadDisplacement%TranslationDisp(1,i) = 0.0  
        MAP_Input(2)%PtFairleadDisplacement%TranslationDisp(2,i) = 0.0
        MAP_Input(2)%PtFairleadDisplacement%TranslationDisp(3,i) = 0.0
     END DO
  END IF
  DO i = 1,MAP_Input(1)%PtFairleadDisplacement%NNodes
     MAP_Input(1)%PtFairleadDisplacement%TranslationDisp(1,i) = 0.0  
     MAP_Input(1)%PtFairleadDisplacement%TranslationDisp(2,i) = 0.0
     MAP_Input(1)%PtFairleadDisplacement%TranslationDisp(3,i) = 0.0
  END DO

  
  ! -------------------------------------------------------------------------
  ! BEGIN time marching
  ! -------------------------------------------------------------------------
  DO n_t_global = 0, n_t_final

     t_global =  t_initial + dt_global*n_t_global

     !==========   NOTE   ======     <-----------------------------------------+
     ! @bonnie : I am assuming this MAP_InputTimes{:} and MAP_Input{:} 
     !           will be assigned by the glue code   

     !MAP_InputTimes(1) = t_global + dt_global
     !MAP_InputTimes(2) = MAP_InputTimes(1) - dt_global 
     !MAP_InputTimes(3) = MAP_InputTimes(2) - dt_global

     !MAP_Input(1)%PtFairleadDisplacement%TranslationDisp(1,1) = .001*n_t_global  
     !MAP_Input(2)%PtFairleadDisplacement%TranslationDisp(1,1) = .001*n_t_global  
     !MAP_Input(3)%PtFairleadDisplacement%TranslationDisp(1,1) = .001*n_t_global  
     !===========================================================================


     ! @bonnie & @jason: the FAST glue code will update the new fairlead position 
     !                   based on the new platform position in the global frame.
     CALL  MAP_UpdateStates( t_global            , &
                             n_t_global          , &
                             MAP_Input           , &
                             MAP_InputTimes      , &
                             MAP_Parameter       , &
                             MAP_ContinuousState , &
                             MAP_DiscreteState   , &
                             MAP_ConstraintState , &
                             MAP_OtherState      , &
                             ErrStat             , &
                             ErrMsg )    
     IF ( ErrStat .NE. 0 ) THEN
        CALL WrScr(ErrMsg) 
    END IF
  
     CALL MAP_CalcOutput( t_global            , &
                          MAP_Input(1)        , &
                          MAP_Parameter       , &
                          MAP_ContinuousState , &
                          MAP_DiscreteState   , &
                          MAP_ConstraintState , &
                          MAP_OtherState      , &
                          MAP_Output          , &
                          ErrStat             , &
                          ErrMsg )
     IF ( ErrStat .NE. 0 ) THEN
        CALL WrScr(ErrMsg) 
     END IF
  
     ! update the global time step by one delta t
     t_global = ( n_t_global + 1 )* dt_global + t_initial
  END DO
  ! -------------------------------------------------------------------------
  ! END time marching
  ! -------------------------------------------------------------------------


  !==========   NOTE   ======     <-----------------------------------------+
  ! @bonnie : I am assuming the glue code will do this
  IF (MAP_interp_order .EQ. 1) THEN  
     CALL MeshDestroy(MAP_Input(2)%PtFairleadDisplacement, ErrStat,ErrMsg)
  ELSE IF (MAP_interp_order .EQ. 2) THEN
     CALL MeshDestroy(MAP_Input(2)%PtFairleadDisplacement, ErrStat,ErrMsg)
     CALL MeshDestroy(MAP_Input(3)%PtFairleadDisplacement, ErrStat,ErrMsg)
  END IF
  !===========================================================================

  
  ! Destroy all objects
  CALL MAP_End( MAP_Input(1)        , &
                MAP_Parameter       , &
                MAP_ContinuousState , &
                MAP_DiscreteState   , &
                MAP_ConstraintState , & 
                MAP_OtherState      , &
                MAP_Output          , &
                ErrStat             , &
                ErrMsg )  
  IF ( ErrStat .NE. 0 ) THEN
     WRITE(*,*) ErrMsg 
  END IF  

  DEALLOCATE(MAP_InputTimes)
  DEALLOCATE(MAP_OutputTimes)
  
  WRITE(*,*) "Program has ended"
    
END PROGRAM Main
