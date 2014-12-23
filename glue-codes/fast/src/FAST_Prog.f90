!**********************************************************************************************************************************
! The FAST_Prog.f90, FAST_IO.f90, and FAST_Mods.f90 make up the FAST glue code in the FAST Modularization Framework.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2014  National Renewable Energy Laboratory
!
!    This file is part of FAST.
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
! File last committed: $Date$
! (File) Revision #: $Rev$
! URL: $HeadURL$
!**********************************************************************************************************************************
PROGRAM FAST
! This program models 2- or 3-bladed turbines of a standard configuration.
!
! noted compilation switches:
!   SOLVE_OPTION_1_BEFORE_2 (uses a different order for solving input-output relationships)
!   OUTPUT_ADDEDMASS        (outputs a file called "<RootName>.AddedMass" that contains HydroDyn's added-mass matrix.
!   OUTPUT_JACOBIAN
!   FPE_TRAP_ENABLED        (use with gfortran when checking for floating point exceptions)
!.................................................................................................


   USE FAST_IO_Subs   ! all of the ModuleName and ModuleName_types modules are inherited from FAST_IO_Subs
                       
IMPLICIT  NONE
   
   ! Local variables:
   REAL(DbKi), PARAMETER     :: t_initial = 0.0_DbKi                    ! Initial time

   
   ! Data for the glue code:
TYPE(FAST_ParameterType)              :: p_FAST                                  ! Parameters for the glue code (bjj: made global for now)
TYPE(FAST_OutputType)                 :: y_FAST                                  ! Output variables for the glue code
TYPE(FAST_MiscVarType)                :: m_FAST                                  ! Miscellaneous variables

TYPE(FAST_ModuleMapType)              :: MeshMapData                             ! Data for mapping between modules
   
TYPE(ElastoDyn_Data)                  :: ED                                      ! Data for the ElastoDyn module
TYPE(ServoDyn_Data)                   :: SrvD                                    ! Data for the ServoDyn module
TYPE(AeroDyn_Data)                    :: AD                                      ! Data for the AeroDyn module
TYPE(InflowWind_Data)                 :: IfW                                     ! Data for InflowWind module
TYPE(HydroDyn_Data)                   :: HD                                      ! Data for the HydroDyn module
TYPE(SubDyn_Data)                     :: SD                                      ! Data for the SubDyn module
TYPE(MAP_Data)                        :: MAPp                                    ! Data for the MAP (Mooring Analysis Program) module
TYPE(FEAMooring_Data)                 :: FEAM                                    ! Data for the FEAMooring module
TYPE(IceFloe_Data)                    :: IceF                                    ! Data for the IceFloe module
TYPE(IceDyn_Data)                     :: IceD                                    ! Data for the IceDyn module

   ! Other/Misc variables

REAL(DbKi)                            :: t_global_next                           ! next simulation time (m_FAST%t_global + p_FAST%dt)


INTEGER(IntKi)                        :: I,J                                     ! generic loop counter
INTEGER(IntKi)                        :: n_t_global                              ! simulation time step, loop counter for global (FAST) simulation
INTEGER(IntKi)                        :: j_pc                                    ! predictor-corrector loop counter 
INTEGER(IntKi)                        :: ErrStat                                 ! Error status
CHARACTER(1024)                       :: ErrMsg                                  ! Error message


   !...............................................................................................................................
   ! initialization
   !...............................................................................................................................

   CALL FAST_InitializeAll( t_initial, p_FAST, y_FAST, m_FAST, ED, SrvD, AD, IfW, HD, SD, MAPp, FEAM, IceF, IceD, MeshMapData, ErrStat, ErrMsg )
      CALL CheckError( ErrStat, ErrMsg, 'during module initialization' )
                                                  
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! loose coupling
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   !...............................................................................................................................
   ! Initialization: (calculate outputs based on states at t=t_initial as well as guesses of inputs and constraint states)
   !...............................................................................................................................     
   CALL FAST_Solution0(p_FAST, y_FAST, m_FAST, ED, SrvD, AD, IfW, HD, SD, MAPp, FEAM, IceF, IceD, MeshMapData, ErrStat, ErrMsg )
      CALL CheckError( ErrStat, ErrMsg, 'during simulation initialization'  )
              
   !...............................................................................................................................
   ! Time Stepping:
   !...............................................................................................................................         
   
   DO n_t_global = 0, m_FAST%n_TMax_m1
      ! this takes data from n_t_global and gets values at n_t_global + 1
  
      t_global_next = t_initial + (n_t_global+1)*p_FAST%DT  ! = m_FAST%t_global + p_FAST%dt
                       
         ! determine if the Jacobian should be calculated this time
      IF ( m_FAST%calcJacobian ) THEN ! this was true (possibly at initialization), so we'll advance the time for the next calculation of the Jacobian
         m_FAST%NextJacCalcTime = m_FAST%t_global + p_FAST%DT_UJac         
      END IF
      
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Step 1.a: Extrapolate Inputs (gives predicted values at t+dt)
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL FAST_ExtrapInterpMods( t_global_next, p_FAST, y_FAST, m_FAST, ED, SrvD, AD, IfW, HD, SD, MAPp, FEAM, IceF, IceD, ErrStat, ErrMsg )
         CALL CheckError( ErrStat, ErrMsg  )

      
      ! predictor-corrector loop:
      DO j_pc = 0, p_FAST%NumCrctn
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Step 1.b: Advance states (yield state and constraint values at t_global_next)
      !
      ! x, xd, and z contain values at m_FAST%t_global;
      ! values at t_global_next are stored in the *_pred variables.
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         !----------------------------------------------------------------------------------------
         ! copy the states at step m_FAST%t_global and get prediction for step t_global_next
         ! (note that we need to copy the states because UpdateStates updates the values
         ! and we need to have the old values [at m_FAST%t_global] for the next j_pc step)
         !----------------------------------------------------------------------------------------
      
         CALL FAST_AdvanceStates( t_initial, n_t_global, p_FAST, y_FAST, m_FAST, ED, SrvD, AD, IfW, HD, SD, MAPp, FEAM, IceF, IceD, ErrStat, ErrMsg )               
            CALL CheckError( ErrStat, ErrMsg )
         
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Step 1.c: Input-Output Solve      
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            CALL CalcOutputs_And_SolveForInputs( n_t_global, m_FAST%t_global,  STATE_PRED, m_FAST%calcJacobian, m_FAST%NextJacCalcTime, &
                        p_FAST, ED, SrvD, AD, IfW, HD, SD, MAPp, FEAM, IceF, IceD, MeshMapData, ErrStat, ErrMsg )
            
               CALL CheckError( ErrStat, ErrMsg )
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Step 2: Correct (continue in loop) 
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         IF ( j_pc /= p_FAST%NumCrctn)  THEN          ! Don't copy these on the last loop iteration...
                  
            IF ( p_FAST%n_substeps( Module_ED ) > 1 ) THEN
               CALL ED_CopyOtherState( ED%OtherSt_old, ED%OtherSt, MESH_UPDATECOPY, Errstat, ErrMsg)
            END IF
            
            IF ( p_FAST%n_substeps( Module_AD ) > 1 ) THEN
               CALL AD_CopyOtherState( AD%OtherSt_old, AD%OtherSt, MESH_UPDATECOPY, Errstat, ErrMsg)
            END IF
            
            IF ( p_FAST%n_substeps( Module_SrvD ) > 1 ) THEN
               CALL SrvD_CopyOtherState( SrvD%OtherSt_old, SrvD%OtherSt, MESH_UPDATECOPY, Errstat, ErrMsg)
            END IF
            
            IF ( p_FAST%n_substeps( Module_HD ) > 1 ) THEN
               CALL HydroDyn_CopyOtherState( HD%OtherSt_old, HD%OtherSt, MESH_UPDATECOPY, Errstat, ErrMsg)
            END IF
            
            IF ( p_FAST%n_substeps( Module_SD ) > 1 ) THEN
               CALL SD_CopyOtherState( SD%OtherSt_old, SD%OtherSt, MESH_UPDATECOPY, Errstat, ErrMsg)
            END IF

            IF ( p_FAST%n_substeps( Module_MAP ) > 1 ) THEN
               CALL MAP_CopyOtherState( MAPp%OtherSt_old, MAPp%OtherSt, MESH_UPDATECOPY, Errstat, ErrMsg)
            ELSEIF ( p_FAST%n_substeps( Module_FEAM ) > 1 ) THEN
               CALL FEAM_CopyOtherState( FEAM%OtherSt_old, FEAM%OtherSt, MESH_UPDATECOPY, Errstat, ErrMsg)
            END IF
         
            IF ( p_FAST%n_substeps( Module_IceF ) > 1 ) THEN
               CALL IceFloe_CopyOtherState( IceF%OtherSt_old, IceF%OtherSt, MESH_UPDATECOPY, Errstat, ErrMsg)
            ELSEIF ( p_FAST%n_substeps( Module_IceD ) > 1 ) THEN
               DO i=1,p_FAST%numIceLegs
                  CALL IceD_CopyOtherState( IceD%OtherSt_old(i), IceD%OtherSt(i), MESH_UPDATECOPY, Errstat, ErrMsg)
               END DO
            END IF
            
         END IF
                              
      enddo ! j_pc
      
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Step 3: Save all final variables (advance to next time)
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      !----------------------------------------------------------------------------------------
      ! copy the final predicted states from step t_global_next to actual states for that step
      !----------------------------------------------------------------------------------------
      
      ! ElastoDyn: copy final predictions to actual states
      CALL ED_CopyContState   (ED%x( STATE_PRED), ED%x( STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg)
      CALL ED_CopyDiscState   (ED%xd(STATE_PRED), ED%xd(STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg)  
      CALL ED_CopyConstrState (ED%z( STATE_PRED), ED%z( STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg)      
      
      
      ! AeroDyn: copy final predictions to actual states; copy current outputs to next 
      IF ( p_FAST%CompAero == Module_AD ) THEN
         CALL AD_CopyContState   (AD%x( STATE_PRED), AD%x( STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg)
         CALL AD_CopyDiscState   (AD%xd(STATE_PRED), AD%xd(STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg)  
         CALL AD_CopyConstrState (AD%z( STATE_PRED), AD%z( STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg)      
      END IF
            
      
      ! ServoDyn: copy final predictions to actual states; copy current outputs to next 
      IF ( p_FAST%CompServo == Module_SrvD ) THEN
         CALL SrvD_CopyContState   (SrvD%x( STATE_PRED), SrvD%x( STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg)
         CALL SrvD_CopyDiscState   (SrvD%xd(STATE_PRED), SrvD%xd(STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg)  
         CALL SrvD_CopyConstrState (SrvD%z( STATE_PRED), SrvD%z( STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg)      
      END IF
      
      
      ! HydroDyn: copy final predictions to actual states
      IF ( p_FAST%CompHydro == Module_HD ) THEN         
         CALL HydroDyn_CopyContState   (HD%x( STATE_PRED), HD%x( STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg)
         CALL HydroDyn_CopyDiscState   (HD%xd(STATE_PRED), HD%xd(STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg)  
         CALL HydroDyn_CopyConstrState (HD%z( STATE_PRED), HD%z( STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg)
      END IF
            
            
      ! SubDyn: copy final predictions to actual states
      IF ( p_FAST%CompSub == Module_SD ) THEN
         CALL SD_CopyContState   (SD%x( STATE_PRED), SD%x( STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg)
         CALL SD_CopyDiscState   (SD%xd(STATE_PRED), SD%xd(STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg)  
         CALL SD_CopyConstrState (SD%z( STATE_PRED), SD%z( STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg)
      END IF
         
      
      ! MAP: copy final predictions to actual states
      IF (p_FAST%CompMooring == Module_MAP) THEN
         CALL MAP_CopyContState   (MAPp%x( STATE_PRED), MAPp%x( STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg)
         CALL MAP_CopyDiscState   (MAPp%xd(STATE_PRED), MAPp%xd(STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg)  
         CALL MAP_CopyConstrState (MAPp%z( STATE_PRED), MAPp%z( STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg)
      ELSEIF (p_FAST%CompMooring == Module_FEAM) THEN
         CALL FEAM_CopyContState   (FEAM%x( STATE_PRED), FEAM%x( STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg)
         CALL FEAM_CopyDiscState   (FEAM%xd(STATE_PRED), FEAM%xd(STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg)  
         CALL FEAM_CopyConstrState (FEAM%z( STATE_PRED), FEAM%z( STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg)
      END IF
             
            ! IceFloe: copy final predictions to actual states
      IF ( p_FAST%CompIce == Module_IceF ) THEN
         CALL IceFloe_CopyContState   (IceF%x( STATE_PRED), IceF%x( STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg)
         CALL IceFloe_CopyDiscState   (IceF%xd(STATE_PRED), IceF%xd(STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg)  
         CALL IceFloe_CopyConstrState (IceF%z( STATE_PRED), IceF%z( STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg)
      ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
         DO i=1,p_FAST%numIceLegs
            CALL IceD_CopyContState   (IceD%x( i,STATE_PRED), IceD%x( i,STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg)
            CALL IceD_CopyDiscState   (IceD%xd(i,STATE_PRED), IceD%xd(i,STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg)  
            CALL IceD_CopyConstrState (IceD%z( i,STATE_PRED), IceD%z( i,STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg)
         END DO
      END IF

            
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! We've advanced everything to the next time step: 
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                       
      
      ! update the global time 
  
      m_FAST%t_global = t_global_next 
      
      
      !----------------------------------------------------------------------------------------
      ! Check to see if we should output data this time step:
      !----------------------------------------------------------------------------------------

      CALL WriteOutputToFile(m_FAST%t_global, p_FAST, y_FAST, ED, AD, IfW, HD, SD, SrvD, MAPp, FEAM, IceF, IceD, ErrStat, ErrMsg)
      
      !----------------------------------------------------------------------------------------
      ! Display simulation status every SttsTime-seconds (i.e., n_SttsTime steps):
      !----------------------------------------------------------------------------------------   
      
      IF ( MOD( n_t_global + 1, p_FAST%n_SttsTime ) == 0 ) THEN

         CALL SimStatus( m_FAST%TiLstPrn, m_FAST%PrevClockTime, m_FAST%t_global, p_FAST%TMax )

      ENDIF
      
            
   END DO ! n_t_global
  
  
   !...............................................................................................................................
   !  Write simulation times and stop
   !...............................................................................................................................
   
   CALL ExitThisProgram( p_FAST, y_FAST, m_FAST, ED, SrvD, AD, IfW, HD, SD, MAPp, FEAM, IceF, IceD, MeshMapData, ErrID_None )


CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg,ErrLocMsg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN)           :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN)           :: Msg         ! The error message (ErrMsg)
      CHARACTER(*),   INTENT(IN), OPTIONAL :: ErrLocMsg   ! an optional message describing the location of the error

      CHARACTER(1024)                      :: SimMsg      
      
      IF ( ErrID /= ErrID_None ) THEN
         CALL WrScr( NewLine//TRIM(Msg)//NewLine )
         IF ( ErrID >= AbortErrLev ) THEN
            
            IF (PRESENT(ErrLocMsg)) THEN
               SimMsg = ErrLocMsg
            ELSE
               SimMsg = 'at simulation time '//TRIM(Num2LStr(m_FAST%t_global))//' of '//TRIM(Num2LStr(p_FAST%TMax))//' seconds'
            END IF
            
            CALL ExitThisProgram( p_FAST, y_FAST, m_FAST, ED, SrvD, AD, IfW, HD, SD, MAPp, FEAM, IceF, IceD, MeshMapData, ErrID, SimMsg )
            
         END IF
         
      END IF


   END SUBROUTINE CheckError   
   !...............................................................................................................................
END PROGRAM FAST
!=======================================================================
