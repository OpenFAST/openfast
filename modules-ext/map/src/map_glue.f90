! /****************************************************************
!  *   Copyright (C) 2014 mdm                                     *
!  *   map[dot]plus[dot]plus[dot]help[at]gmail                     *
!  *                                                              *
!  * Licensed to the Apache Software Foundation (ASF) under one   *
!  * or more contributor license agreements.  See the NOTICE file *
!  * distributed with this work for additional information        *
!  * regarding copyright ownership.  The ASF licenses this file   *
!  * to you under the Apache License, Version 2.0 (the            *
!  * "License"); you may not use this file except in compliance   *
!  * with the License.  You may obtain a copy of the License at   *
!  *                                                              *
!  *   http://www.apache.org/licenses/LICENSE-2.0                 *
!  *                                                              *
!  * Unless required by applicable law or agreed to in writing,   *
!  * software distributed under the License is distributed on an  *
!  * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY       *
!  * KIND, either express or implied.  See the License for the    *
!  * specific language governing permissions and limitations      *      
!  * under the License.                                           *  
!  ****************************************************************/

PROGRAM Main
  
  USE MAP_Types
  USE MAP
  USE NWTC_Library 

  IMPLICIT NONE 

  INTEGER(IntKi)                         :: i,j, j_ss        ! generic loop counter
  INTEGER(IntKi)                         :: ErrStat          ! Status of error message   
  CHARACTER(1024)                        :: ErrMsg           ! Error message if ErrStat /= ErrID_None
                                         
  REAL(DbKi)                             :: dt_global        ! fixed/constant global time step
  REAL(DbKi)                             :: t_initial        ! time at initialization
  REAL(DbKi)                             :: t_final          ! time at simulation end 
  REAL(DbKi)                             :: t_global         ! global-loop time marker
                                         
  INTEGER(IntKi)                         :: n_t_final        ! total number of time steps
  INTEGER(IntKi)                         :: n_t_global       ! global-loop time counter

  TYPE (MAP_InitInputType)               :: InitInData_MAP    
  TYPE (MAP_InitOutputType)              :: InitOutData_MAP    
  TYPE (MAP_InputType)                   :: u_MAP
  TYPE (MAP_ParameterType)               :: p_MAP
  TYPE (MAP_ContinuousStateType)         :: x_MAP
  TYPE (MAP_DiscreteStateType)           :: xd_MAP
  TYPE (MAP_ConstraintStateType)         :: z_MAP
  TYPE (MAP_OtherStateType)              :: other_MAP


  TYPE(MAP_ContinuousStateType)         :: x_MAP_pred                              ! Predicted continuous states
  TYPE(MAP_DiscreteStateType)           :: xd_MAP_pred                             ! Predicted discrete states
  TYPE(MAP_ConstraintStateType)         :: z_MAP_pred                              ! Predicted constraint states
  TYPE(MAP_OtherStateType)              :: other_MAP_old                           ! Other/optimization states (copied for the case of subcycling)

  TYPE (MAP_InputType),      ALLOCATABLE :: MAP_Input(:)
  REAL(DbKi) , DIMENSION(:), ALLOCATABLE :: MAP_InputTimes(:)

  TYPE (MAP_OutputType)                  :: y_MAP
!  REAL(DbKi) , DIMENSION(:), ALLOCATABLE :: MAP_OutputTimes

  INTEGER(IntKi)                         :: MAP_interp_order     ! order of interpolation/extrapolation

  
  ! -------------------------------------------------------------------------
  ! Initialization of glue-code time-step variables
  ! -------------------------------------------------------------------------
  
  t_initial = 0.
  t_final   = 10.0
  
  ! specify time increment; currently, all modules will be time integrated with this increment size
  dt_global = 0.5
  n_t_final = ((t_final - t_initial) / dt_global ) - 1  
  t_global = t_initial  
  
  ! @bonnie : is this right? What's a good interp order?
  ! @marco: the interp order is specified in the FAST input file. It can be 0, 1, or 2
  MAP_interp_order = 0

  ! MAP: allocate Input and Output arrays; used for interpolation and extrapolation
!  Allocate(MAP_InputTimes(MAP_interp_order + 1)) 

  ! @bonnie : This is in the FAST developers glue code example, but it's probably not needed here. 
  ALLOCATE(MAP_Input(MAP_interp_order + 1), MAP_InputTimes(MAP_interp_order + 1), STAT=ErrStat )
  IF (ErrStat /= 0) CALL CheckError(ErrID_Fatal,"Error allocating MAP_Input and MAP_InputTimes.") 

    
  ! set the MAP input file name and other environment terms.
  InitInData_MAP%file_name = "baseline.map"      
  InitInData_MAP%summary_file_name = "baseline.sum.map"  
  InitInData_MAP%gravity = 9.81       ! @bonnie : This need to be according to g used in FAST. Positive value
  InitInData_MAP%sea_density = 1025    ! @bonnie : This needs to be set according to seawater density in FAST. Positive value
  InitInData_MAP%depth = -320         ! @bonnie : This need to be set according to the water depth in FAST. Negative value
  ! p_MAP%dt = dt_global     ! @bonnie : This is for the glue code to set

  ! call the initialization routine
  ! CALL MAP_Init( InitInData_MAP, MAP_Input(1), p_MAP,  x_MAP, xd_MAP, z_MAP, OtherSt_MAP, y_MAP, p_FAST%dt_module( MODULE_MAP ), InitOutData_MAP, ErrStat, ErrMsg )
  CALL MAP_Init(InitInData_MAP , &
                MAP_Input(1)   , & 
                p_MAP          , &
                x_MAP          , &
                xd_MAP         , &
                z_MAP          , &
                other_MAP      , &
                y_MAP          , &
                dt_global      , &
                InitOutData_MAP, &
                ErrStat        , &
                ErrMsg )
  IF (ErrStat.NE.0) THEN
     CALL WrScr(ErrMsg) 
  END IF

  CALL DispNVD(InitOutData_MAP%Ver) 


  CALL MAP_DestroyInitInput(InitInData_MAP, ErrStat, ErrMsg)    
  CALL MAP_DestroyInitOutput(InitOutData_MAP, ErrStat, ErrMsg)


  DO j = 1, MAP_interp_order + 1
     MAP_InputTimes(j) = t_initial - (j - 1) * dt_global
  END DO

  DO j = 2, MAP_interp_order + 1
     CALL MAP_CopyInput (MAP_Input(1),  MAP_Input(j),  MESH_NEWCOPY, Errstat, ErrMsg)
     CALL CheckError( ErrStat, 'Message from MAP_CopyInput (MAP_Input): '//NewLine//ErrMsg )
  END DO
  CALL MAP_CopyInput (MAP_Input(1),  u_MAP,  MESH_NEWCOPY, Errstat, ErrMsg) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
  CALL CheckError( ErrStat, 'Message from MAP_CopyInput (u_MAP): '//NewLine//ErrMsg )
   
  ! Initialize predicted states for j_pc loop:
  CALL MAP_CopyContState   ( x_MAP,  x_MAP_pred, MESH_NEWCOPY, Errstat, ErrMsg)
  CALL CheckError( ErrStat, 'Message from MAP_CopyContState (init): '//NewLine//ErrMsg )
  CALL MAP_CopyDiscState   (xd_MAP, xd_MAP_pred, MESH_NEWCOPY, Errstat, ErrMsg)  
  CALL CheckError( ErrStat, 'Message from MAP_CopyDiscState (init): '//NewLine//ErrMsg )
  CALL MAP_CopyConstrState ( z_MAP,  z_MAP_pred, MESH_NEWCOPY, Errstat, ErrMsg)
  CALL CheckError( ErrStat, 'Message from MAP_CopyConstrState (init): '//NewLine//ErrMsg )
   
  ! IF ( p_FAST%n_substeps( MODULE_MAP ) > 1 ) THEN
  !    CALL MAP_CopyOtherState( other_MAP, other_MAP_old, MESH_NEWCOPY, Errstat, ErrMsg)
  !    CALL CheckError( ErrStat, 'Message from MAP_CopyOtherState (init): '//NewLine//ErrMsg )   
  ! END IF
     
   
  ! -------------------------------------------------------------------------
  ! BEGIN time marching
  ! -------------------------------------------------------------------------
  DO n_t_global = 0, n_t_final
     t_global =  t_initial + dt_global*n_t_global
     
     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     ! Step 1.a: Extrapolate Inputs (gives predicted values at t+dt
     ! 
     ! a) Extrapolate inputs (and outputs -- bjj: output extrapolation not necessary, yet) 
     !    to t + dt (i.e., t_global_next); will only be used by modules with an implicit dependence on input data.
     ! b) Shift "window" of the ModName_Input and ModName_Output
     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
     ! MAP
     CALL MAP_Input_ExtrapInterp(MAP_Input, MAP_InputTimes, u_MAP, t_global, ErrStat, ErrMsg)
      CALL CheckError(ErrStat,'Message from MAP_Input_ExtrapInterp (FAST): '//NewLine//ErrMsg )              
     DO j = MAP_interp_order, 1, -1
        CALL MAP_CopyInput (MAP_Input(j),  MAP_Input(j+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
        MAP_InputTimes(j+1) = MAP_InputTimes(j)
     END DO
   
     CALL MAP_CopyInput (u_MAP,  MAP_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
     MAP_InputTimes(1) = t_global          
   
   
     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     ! Step 1.b: Advance states (yield state and constraint values at t_global_next)
     !
     ! x, xd, and z contain val0ues at t_global;
     ! values at t_global_next are stored in the *_pred variables.
     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     !----------------------------------------------------------------------------------------
     ! copy the states at step t_global and get prediction for step t_global_next
     ! (note that we need to copy the states because UpdateStates updates the values
     ! and we need to have the old values [at t_global] for the next j_pc step)
     !----------------------------------------------------------------------------------------
     ! ElastoDyn: get predicted states
     ! AeroDyn: get predicted states
     ! ServoDyn: get predicted states
     ! HydroDyn: get predicted states
     ! SubDyn: get predicted states
     ! MAP/FEAM: get predicted states
   
     CALL MAP_CopyContState   ( x_MAP,  x_MAP_pred, MESH_UPDATECOPY, Errstat, ErrMsg)
     CALL MAP_CopyDiscState   (xd_MAP, xd_MAP_pred, MESH_UPDATECOPY, Errstat, ErrMsg)  
     CALL MAP_CopyConstrState ( z_MAP,  z_MAP_pred, MESH_UPDATECOPY, Errstat, ErrMsg)
     
     ! IF ( p_FAST%n_substeps( Module_MAP ) > 1 ) THEN
     !    CALL MAP_CopyOtherState( OtherSt_MAP, OtherSt_MAP_old, MESH_UPDATECOPY, Errstat, ErrMsg)
     ! END IF
         
     DO j_ss = 1, 1 !p_FAST%n_substeps( Module_MAP )
        ! n_t_module = n_t_global*p_FAST%n_substeps( Module_MAP ) + j_ss - 1
        ! t_module   = n_t_module*p_FAST%dt_module( Module_MAP )           
        CALL  MAP_UpdateStates(t_global       , &
                               n_t_global     , &
                               MAP_Input      , &
                               MAP_InputTimes , &
                               p_MAP          , &
                               x_MAP_pred     , &
                               xd_MAP_pred    , &
                               z_MAP_pred     , &
                               other_MAP      , &
                               ErrStat        , &
                               ErrMsg )    
        IF (ErrStat.NE.0) THEN
           CALL WrScr(trim(ErrMsg)) 
        END IF
        ! CALL CheckError( ErrStat, 'Message from MAP_UpdateStates: '//NewLine//ErrMsg )
     END DO !j_ss
        
     
     ! !==========   NOTE   ======     <-----------------------------------------+
     ! ! @bonnie : I am assuming this MAP_InputTimes{:} and MAP_Input{:} 
     ! !           will be assigned by the glue code   
     ! 
     ! MAP_InputTimes(1) = t_global + dt_global
     ! ! MAP_InputTimes(2) = MAP_InputTimes(1) - dt_global 
     ! ! MAP_InputTimes(3) = MAP_InputTimes(2) - dt_global
     ! 
     MAP_Input(1)%PtFairDisplacement%TranslationDisp(1,1) = .05*n_t_global  
     MAP_Input(1)%PtFairDisplacement%TranslationDisp(1,2) = .05*n_t_global  
     ! ! MAP_Input(3)%PtFairDisplacement%TranslationDisp(1,1) = 1*n_t_global  
     ! !===========================================================================
          
     ! @bonnie & @jason: the FAST glue code will update the new fairlead position 
     !                   based on the new platform position in the global frame.
     ! CALL MAP_CalcOutput( this_time, MAP_Input(1), p_MAP, x_MAP_this, xd_MAP_this, z_MAP_this, OtherSt_MAP, y_MAP, ErrStat, ErrMsg )
     CALL MAP_CalcOutput(t_global     , &
                         MAP_Input(1) , &
                         p_MAP        , &
                         x_MAP_pred   , &
                         xd_MAP_pred  , &
                         z_MAP_pred   , &
                         other_MAP    , &
                         y_MAP        , &
                         ErrStat      , &
                         ErrMsg )
     IF (ErrStat.NE.0) THEN
        CALL WrScr(ErrMsg) 
     END IF
   
   
   
     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     ! Step 2: Correct (continue in loop) 
     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     ! IF ( j_pc /= p_FAST%NumCrctn)  THEN          ! Don't copy these on the last loop iteration...                  
     !    IF ( p_FAST%n_substeps( Module_MAP ) > 1 ) THEN
     !       CALL MAP_CopyOtherState( OtherSt_MAP_old, OtherSt_MAP, MESH_UPDATECOPY, Errstat, ErrMsg)
     !    ELSEIF ( p_FAST%n_substeps( Module_FEAM ) > 1 ) THEN
     !       CALL FEAM_CopyOtherState( OtherSt_FEAM_old, OtherSt_FEAM, MESH_UPDATECOPY, Errstat, ErrMsg)
     !    END IF
     ! END IF
     
     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     ! Predictor-corrector ends here!
     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   
     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     ! Step 3: Save all final variables (advance to next time)
     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     ! MAP: copy final predictions to actual states
     CALL MAP_CopyContState   ( x_MAP_pred,  x_MAP, MESH_UPDATECOPY, Errstat, ErrMsg)
     CALL MAP_CopyDiscState   (xd_MAP_pred, xd_MAP, MESH_UPDATECOPY, Errstat, ErrMsg)  
     CALL MAP_CopyConstrState ( z_MAP_pred,  z_MAP, MESH_UPDATECOPY, Errstat, ErrMsg)
  END DO
  ! -------------------------------------------------------------------------
  ! END time marching
  ! -------------------------------------------------------------------------


   !==========   NOTE   ======     <-----------------------------------------+
   ! @bonnie : I am assuming the glue code will do this
   IF (MAP_interp_order .EQ. 1) THEN  
      CALL MeshDestroy(MAP_Input(2)%PtFairDisplacement, ErrStat,ErrMsg)
   ELSE IF (MAP_interp_order .EQ. 2) THEN
      CALL MeshDestroy(MAP_Input(2)%PtFairDisplacement, ErrStat,ErrMsg)
      CALL MeshDestroy(MAP_Input(3)%PtFairDisplacement, ErrStat,ErrMsg)
   END IF
   !===========================================================================  

   ! Destroy all objects
   ! CALL MAP_End(    MAP_Input(1),   p_MAP,   x_MAP,   xd_MAP,   z_MAP,   OtherSt_MAP,   y_MAP,   ErrStat2, ErrMsg2)
   CALL MAP_End(MAP_Input(1), &
                p_MAP       , &
                x_MAP       , &
                xd_MAP      , &
                z_MAP       , &
                other_MAP   , &
                y_MAP       , &
                ErrStat     , &
                ErrMsg )  
   IF (ErrStat.NE.0) THEN
      WRITE(*,*) ErrMsg 
   END IF  
   
   CALL MAP_DestroyInput(u_MAP, ErrStat, ErrMsg)
   IF (ErrStat/=ErrID_None) CALL WrScr(TRIM(ErrMsg))  
   IF (ALLOCATED(MAP_Input)) THEN
      DO j = 2,MAP_interp_order+1 !note that SD_Input(1) was destroyed in MAP_End       
         CALL MAP_DestroyInput(MAP_Input(j), ErrStat, ErrMsg)
         IF (ErrStat/=ErrID_None) CALL WrScr(TRIM(ErrMsg))
      END DO
      DEALLOCATE(MAP_Input)
   END IF
   IF (ALLOCATED(MAP_Input)) DEALLOCATE(MAP_Input)      
   IF (ALLOCATED(MAP_InputTimes)) DEALLOCATE(MAP_InputTimes)
   
   CALL MAP_DestroyContState(x_MAP_pred, ErrStat, ErrMsg); IF (ErrStat/= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   CALL MAP_DestroyDiscState(xd_MAP_pred, ErrStat, ErrMsg); IF (ErrStat/= ErrID_None) CALL WrScr(TRIM(ErrMsg))  
   CALL MAP_DestroyConstrState(z_MAP_pred, ErrStat, ErrMsg); IF (ErrStat/= ErrID_None) CALL WrScr(TRIM(ErrMsg))  
   CALL MAP_DestroyOtherState(other_MAP_old, ErrStat, ErrMsg); IF (ErrStat/= ErrID_None) CALL WrScr(TRIM(ErrMsg))  

 CONTAINS 

    SUBROUTINE CheckError(ErrID,Msg)
      ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)
      IF ( ErrID /= ErrID_None ) THEN
         CALL WrScr( NewLine//TRIM(Msg)//NewLine )
         ! IF ( ErrID >= AbortErrLev ) CALL ExitThisProgram( Error=.TRUE., ErrLev=ErrID )
      END IF
   END SUBROUTINE CheckError   

END PROGRAM Main


