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
!   OUTPUT_INPUTMESHES
!.................................................................................................


   USE FAST_IO_Subs   ! all of the ModuleName and ModuleName_types modules are inherited from FAST_IO_Subs
                       
IMPLICIT  NONE
   
   ! Local variables:

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
REAL                                  :: UsrTimeDiff                             ! Difference in CPU time from start to finish of program execution


INTEGER(IntKi)                        :: I,J                                     ! generic loop counter
INTEGER(IntKi)                        :: n_t_global                              ! simulation time step, loop counter for global (FAST) simulation
INTEGER(IntKi)                        :: j_pc                                    ! predictor-corrector loop counter 
INTEGER(IntKi)                        :: ErrStat                                 ! Error status
CHARACTER(1024)                       :: ErrMsg                                  ! Error message


   !...............................................................................................................................
   ! initialization
   !...............................................................................................................................

   CALL FAST_InitializeAll( p_FAST, y_FAST, m_FAST, ED, SrvD, AD, IfW, HD, SD, MAPp, FEAM, IceF, IceD, MeshMapData, ErrStat, ErrMsg )
      CALL CheckError( ErrStat, ErrMsg )
                                                  
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! loose coupling
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !CALL WrScr1 ( '' )

   
   !...............................................................................................................................
   ! Initialization: (calculate outputs based on states at t=t_initial as well as guesses of inputs and constraint states)
   !...............................................................................................................................
   
   m_FAST%t_global   = t_initial
   n_t_global = -1  ! initialize here because CalcOutputs_And_SolveForInputs uses it
   j_PC       = -1
   m_FAST%n_TMax_m1  = CEILING( ( (p_FAST%TMax - t_initial) / p_FAST%DT ) ) - 1 ! We're going to go from step 0 to n_TMax (thus the -1 here)
     
  
   CALL SimStatus_FirstTime( m_FAST%TiLstPrn, m_FAST%PrevClockTime, m_FAST%SimStrtTime, m_FAST%UsrTime2, m_FAST%t_global, p_FAST%TMax )

   ! Solve input-output relations; this section of code corresponds to Eq. (35) in Gasmi et al. (2013)
   ! This code will be specific to the underlying modules
   
#ifdef SOLVE_OPTION_1_BEFORE_2
! used for Option 1 before Option 2:

   IF ( p_FAST%CompSub == Module_SD .OR. p_FAST%CompHydro == Module_HD ) THEN
   ! Because SubDyn needs a better initial guess from ElastoDyn, we'll add an additional call to ED_CalcOutput to get them:
   ! (we'll do the same for HydroDyn, though I'm not sure it's as critical)
   
      CALL ED_CalcOutput( m_FAST%t_global, ED%Input(1), ED%p, ED%x(STATE_CURR), ED%xd(STATE_CURR), ED%z(STATE_CURR), ED%OtherSt, &
                          ED%Output(1), ErrStat, ErrMsg )
         CALL CheckError( ErrStat, ErrMsg  )    
      
      CALL Transfer_ED_to_HD_SD_Mooring( p_FAST, ED%Output(1), HD%Input(1), SD%Input(1), MAPp%Input(1), FEAM%Input(1), MeshMapData, ErrStat, ErrMsg )         
         CALL CheckError( ErrStat, ErrMsg  )    
               
   END IF   
#endif   

   CALL CalcOutputs_And_SolveForInputs(  n_t_global, m_FAST%t_global,  STATE_CURR, m_FAST%calcJacobian, m_FAST%NextJacCalcTime, &
                        p_FAST, ED, SrvD, AD, IfW, HD, SD, MAPp, FEAM, IceF, IceD, MeshMapData, ErrStat, ErrMsg )
         
      CALL CheckError( ErrStat, ErrMsg  )
   
#ifdef OUTPUT_INPUTMESHES
   CALL WriteInputMeshesToFile( ED%Input(1), SD%Input(1), HD%Input(1), MAPp%Input(1), AD%Input(1), TRIM(p_FAST%OutFileRoot)//'.InputMeshes.bin', ErrStat, ErrMsg) 
#else
      IF (p_FAST%WrGraphics) THEN
         CALL WriteInputMeshesToFile( ED%Input(1), SD%Input(1), HD%Input(1), MAPp%Input(1), AD%Input(1), TRIM(p_FAST%OutFileRoot)//'.InputMeshes.bin', ErrStat, ErrMsg) 
      END IF 
#endif

      !----------------------------------------------------------------------------------------
      ! Check to see if we should output data this time step:
      !----------------------------------------------------------------------------------------

      CALL WriteOutputToFile(m_FAST%t_global, p_FAST, y_FAST, ED, AD, IfW, HD, SD, SrvD, MAPp, FEAM, IceF, IceD, ErrStat, ErrMsg)   
         CALL CheckError( ErrStat, ErrMsg  )
   !...............
   ! Copy values of these initial guesses for interpolation/extrapolation and 
   ! initialize predicted states for j_pc loop (use MESH_NEWCOPY here so we can use MESH_UPDATE copy later)
   !...............
         
   ! Initialize Input-Output arrays for interpolation/extrapolation:

   CALL FAST_InitIOarrays( t_initial, p_FAST, y_FAST, m_FAST, ED, SrvD, AD, IfW, HD, SD, MAPp, FEAM, IceF, IceD, ErrStat, ErrMsg )
      CALL CheckError( ErrStat, ErrMsg  )

           
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
         CALL CheckError( ErrStat, ErrMsg )

      
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
   n_t_global =  m_FAST%n_TMax_m1 + 1               ! set this for the message in ProgAbort, if necessary
   CALL ExitThisProgram( Error=.FALSE. )


CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


      IF ( ErrID /= ErrID_None ) THEN
         CALL WrScr( NewLine//TRIM(Msg)//NewLine )
         IF ( ErrID >= AbortErrLev ) CALL ExitThisProgram( Error=.TRUE., ErrLev=ErrID )
      END IF


   END SUBROUTINE CheckError   
   !...............................................................................................................................  
   SUBROUTINE ExitThisProgram( Error, ErrLev )
   ! This subroutine is called when FAST exits. It calls all the modules' end routines and cleans up variables declared in the
   ! main program. If there was an error, it also aborts. Otherwise, it prints the run times and performs a normal exit.
   !...............................................................................................................................

         ! Passed arguments
      LOGICAL,        INTENT(IN)           :: Error        ! flag to determine if this is an abort or normal stop
      INTEGER(IntKi), INTENT(IN), OPTIONAL :: ErrLev       ! Error level when Error == .TRUE. (required when Error is .TRUE.)

         ! Local arguments:
      INTEGER(IntKi)                       :: ErrStat2                                    ! Error status
      CHARACTER(LEN(ErrMsg))               :: ErrMsg2                                     ! Error message

      
      
      !...............................................................................................................................
      ! Clean up modules (and write binary FAST output file), destroy any other variables
      !...............................................................................................................................
!bjj: if any of these operations produces an error >= AbortErrLev, we should also set Error = TRUE and update ErrLev appropriately.

      CALL FAST_End( p_FAST, y_FAST, ErrStat2, ErrMsg2 )
      IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

      IF ( p_FAST%ModuleInitialized(Module_ED) ) THEN
         CALL ED_End(   ED%Input(1),   ED%p,   ED%x(STATE_CURR),   ED%xd(STATE_CURR),   ED%z(STATE_CURR),   ED%OtherSt,   &
                        ED%Output(1),   ErrStat2, ErrMsg2 )
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      END IF

      IF ( p_FAST%ModuleInitialized(Module_AD) ) THEN
         CALL AD_End(   AD%Input(1),   AD%p,   AD%x(STATE_CURR),   AD%xd(STATE_CURR),   AD%z(STATE_CURR),   AD%OtherSt,   &
                        AD%y,   ErrStat2, ErrMsg2 )
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      END IF
      
      IF ( p_FAST%ModuleInitialized(Module_SrvD) ) THEN
         CALL SrvD_End( SrvD%Input(1), SrvD%p, SrvD%x(STATE_CURR), SrvD%xd(STATE_CURR), SrvD%z(STATE_CURR), SrvD%OtherSt, &
                        SrvD%y, ErrStat2, ErrMsg2 )
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      END IF

      IF ( p_FAST%ModuleInitialized(Module_HD) ) THEN
         CALL HydroDyn_End(    HD%Input(1),   HD%p,   HD%x(STATE_CURR),   HD%xd(STATE_CURR),   HD%z(STATE_CURR),   HD%OtherSt,   &
                               HD%y,   ErrStat2, ErrMsg2)
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      END IF

      IF ( p_FAST%ModuleInitialized(Module_SD) ) THEN
         CALL SD_End(    SD%Input(1),   SD%p,   SD%x(STATE_CURR),   SD%xd(STATE_CURR),   SD%z(STATE_CURR),   SD%OtherSt,   &
                         SD%y,   ErrStat2, ErrMsg2)
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      END IF
      
      IF ( p_FAST%ModuleInitialized(Module_MAP) ) THEN
         CALL MAP_End(    MAPp%Input(1),   MAPp%p,   MAPp%x(STATE_CURR),   MAPp%xd(STATE_CURR),   MAPp%z(STATE_CURR),   MAPp%OtherSt,   &
                          MAPp%y,   ErrStat2, ErrMsg2)
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      ELSEIF ( p_FAST%ModuleInitialized(Module_FEAM) ) THEN
         CALL FEAM_End(   FEAM%Input(1),  FEAM%p,  FEAM%x(STATE_CURR),  FEAM%xd(STATE_CURR),  FEAM%z(STATE_CURR),  FEAM%OtherSt,  &
                          FEAM%y,  ErrStat2, ErrMsg2)
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      END IF
      
      IF ( p_FAST%ModuleInitialized(Module_IceF) ) THEN
         CALL IceFloe_End(IceF%Input(1),  IceF%p,  IceF%x(STATE_CURR),  IceF%xd(STATE_CURR),  IceF%z(STATE_CURR),  IceF%OtherSt, &
                          IceF%y,  ErrStat2, ErrMsg2)
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      ELSEIF ( p_FAST%ModuleInitialized(Module_IceD) ) THEN
         
         DO i=1,p_FAST%numIceLegs                     
            CALL IceD_End(IceD%Input(1,i),  IceD%p(i),  IceD%x(i,STATE_CURR),  IceD%xd(i,STATE_CURR),  IceD%z(i,STATE_CURR), &
                          IceD%OtherSt(i),  IceD%y(i),  ErrStat2, ErrMsg2)
            IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )            
         END DO
         
      END IF
      
                  
      ! -------------------------------------------------------------------------
      ! Deallocate/Destroy structures associated with mesh mapping
      ! -------------------------------------------------------------------------

      CALL FAST_DestroyModuleMapType( MeshMapData, ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))
                              
      ! -------------------------------------------------------------------------
      ! variables for ExtrapInterp:
      ! -------------------------------------------------------------------------

      ! ElastoDyn
      CALL FAST_DestroyElastoDyn_Data( ED, ErrStat2, ErrMsg2 )
      IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )      
      
      ! ServoDyn     
      IF ( ALLOCATED(SrvD%Input)      ) THEN
         
         IF ( p_FAST%CompServo == Module_SrvD ) THEN
         
            CALL SrvD_DestroyOutput( SrvD%y_prev, ErrStat2, ErrMsg2)
               IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
                  
            CALL SrvD_DestroyInput( SrvD%u, ErrStat2, ErrMsg2 )
            IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

            DO j = 2,p_FAST%InterpOrder+1 !note that SrvD%Input(1) was destroyed in SrvD_End
               CALL SrvD_DestroyInput( SrvD%Input(j), ErrStat2, ErrMsg2 )
               IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
            END DO
         END IF
         
         DEALLOCATE( SrvD%Input )
      END IF

      IF ( ALLOCATED(SrvD%InputTimes) ) DEALLOCATE( SrvD%InputTimes )
                           
         
      ! AeroDyn
      IF ( p_FAST%CompAero == Module_AD ) THEN
         CALL AD_DestroyInput( AD%u, ErrStat2, ErrMsg2 )
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

         IF ( ALLOCATED(AD%Input)      )  THEN
            DO j = 2,p_FAST%InterpOrder+1 !note that AD%Input(1) was destroyed in AD_End
               CALL AD_DestroyInput( AD%Input(j), ErrStat2, ErrMsg2 )
               IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
            END DO
            DEALLOCATE( AD%Input )
         END IF
      ELSE
         IF ( ALLOCATED(AD%Input)      ) DEALLOCATE( AD%Input )         
      END IF

      IF ( ALLOCATED(AD%InputTimes) ) DEALLOCATE( AD%InputTimes )
                  
      
      ! HydroDyn
      IF ( p_FAST%CompHydro == Module_HD ) THEN                  
         CALL HydroDyn_DestroyInput( HD%u, ErrStat2, ErrMsg2 )
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

         IF ( ALLOCATED(HD%Input)      )  THEN
            DO j = 2,p_FAST%InterpOrder+1 !note that HD%Input(1) was destroyed in HydroDyn_End
               CALL HydroDyn_DestroyInput( HD%Input(j), ErrStat2, ErrMsg2 )
               IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
            END DO
            DEALLOCATE( HD%Input )
         END IF
      ELSE
         IF ( ALLOCATED(HD%Input)      ) DEALLOCATE( HD%Input )         
      END IF

      IF ( ALLOCATED(HD%InputTimes) ) DEALLOCATE( HD%InputTimes )

      ! SubDyn
      IF ( p_FAST%CompSub == Module_SD ) THEN
         CALL SD_DestroyInput( SD%u, ErrStat2, ErrMsg2 )
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

         IF ( ALLOCATED(SD%Input)      ) THEN
            DO j = 2,p_FAST%InterpOrder+1 !note that SD%Input(1) was destroyed in SD_End
               CALL SD_DestroyInput( SD%Input(j), ErrStat2, ErrMsg2 )
               IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
            END DO
            DEALLOCATE( SD%Input )
         END IF
      ELSE
         IF ( ALLOCATED(SD%Input)      ) DEALLOCATE( SD%Input )
      END IF

      IF ( ALLOCATED(SD%InputTimes) ) DEALLOCATE( SD%InputTimes )
      
      ! MAP      
      IF ( p_FAST%ModuleInitialized(Module_MAP)  ) THEN        
         CALL MAP_DestroyInput( MAPp%u, ErrStat2, ErrMsg2 )
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

         IF ( ALLOCATED(MAPp%Input)      ) THEN
            DO j = 2,p_FAST%InterpOrder+1 !note that SD%Input(1) was destroyed in MAP_End
               CALL MAP_DestroyInput( MAPp%Input(j), ErrStat2, ErrMsg2 )
               IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
            END DO
            DEALLOCATE( MAPp%Input )
         END IF
      ELSE
         IF ( ALLOCATED(MAPp%Input)      ) DEALLOCATE( MAPp%Input )
      END IF

      IF ( ALLOCATED(MAPp%InputTimes) ) DEALLOCATE( MAPp%InputTimes )
      
      
      ! FEAM      
      IF ( p_FAST%ModuleInitialized(Module_FEAM)  ) THEN        
         CALL FEAM_DestroyInput( FEAM%u, ErrStat2, ErrMsg2 )
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

         IF ( ALLOCATED(FEAM%Input)      ) THEN
            DO j = 2,p_FAST%InterpOrder+1 !note that SD%Input(1) was destroyed in MAP_End
               CALL FEAM_DestroyInput( FEAM%Input(j), ErrStat2, ErrMsg2 )
               IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
            END DO
            DEALLOCATE( FEAM%Input )
         END IF
      ELSE
         IF ( ALLOCATED(FEAM%Input)      ) DEALLOCATE( FEAM%Input )
      END IF

      IF ( ALLOCATED(FEAM%InputTimes) ) DEALLOCATE( FEAM%InputTimes )
      
      ! IceFloe
      IF ( p_FAST%CompIce == Module_IceF ) THEN
         CALL IceFloe_DestroyInput( IceF%u, ErrStat2, ErrMsg2 )
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

         IF ( ALLOCATED(IceF%Input)      ) THEN
            DO j = 2,p_FAST%InterpOrder+1 !note that IceF%Input(1) was destroyed in IceFloe_End
               CALL IceFloe_DestroyInput( IceF%Input(j), ErrStat2, ErrMsg2 )
               IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
            END DO
            DEALLOCATE( IceF%Input )
         END IF
      ELSE
         IF ( ALLOCATED(IceF%Input)      ) DEALLOCATE( IceF%Input )
      END IF

      IF ( ALLOCATED(IceF%InputTimes) ) DEALLOCATE( IceF%InputTimes )      
      
      
      ! IceDyn
      IF ( p_FAST%CompIce == Module_IceD ) THEN
                                          
         IF ( ALLOCATED(IceD%Input)  ) THEN
            DO i=1,p_FAST%numIceLegs
               DO j = 2,p_FAST%InterpOrder+1 !note that IceD%Input(1,:) was destroyed in ID_End
                  CALL IceD_DestroyInput( IceD%Input(j,i), ErrStat2, ErrMsg2 );   IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
               END DO
            END DO
         END IF
                              
         IF ( ALLOCATED(IceD%OtherSt_old) ) THEN  ! and all the others that need to be allocated...
            DO i=1,p_FAST%numIceLegs
               CALL IceD_DestroyContState(  IceD%x(          i,STATE_CURR), ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
               CALL IceD_DestroyDiscState(  IceD%xd(         i,STATE_CURR), ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
               CALL IceD_DestroyConstrState(IceD%z(          i,STATE_CURR), ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
               CALL IceD_DestroyContState(  IceD%x(          i,STATE_PRED), ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
               CALL IceD_DestroyDiscState(  IceD%xd(         i,STATE_PRED), ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
               CALL IceD_DestroyConstrState(IceD%z(          i,STATE_PRED), ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
               CALL IceD_DestroyOtherState( IceD%OtherSt(    i), ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
               CALL IceD_DestroyOtherState( IceD%OtherSt_old(i), ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )                              
               CALL IceD_DestroyParam(      IceD%p(          i), ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
               CALL IceD_DestroyInput(      IceD%u(          i), ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
               CALL IceD_DestroyOutput(     IceD%y(          i), ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )                              
            END DO                        
         END IF         
         
      END IF 
      
      IF ( ALLOCATED(IceD%Input      ) ) DEALLOCATE( IceD%Input      )      
      IF ( ALLOCATED(IceD%InputTimes ) ) DEALLOCATE( IceD%InputTimes )      
      IF ( ALLOCATED(IceD%x          ) ) DEALLOCATE( IceD%x          )      
      IF ( ALLOCATED(IceD%xd         ) ) DEALLOCATE( IceD%xd         )      
      IF ( ALLOCATED(IceD%z          ) ) DEALLOCATE( IceD%z          )      
      IF ( ALLOCATED(IceD%OtherSt    ) ) DEALLOCATE( IceD%OtherSt    )      
      IF ( ALLOCATED(IceD%p          ) ) DEALLOCATE( IceD%p          )      
      IF ( ALLOCATED(IceD%u          ) ) DEALLOCATE( IceD%u          )      
      IF ( ALLOCATED(IceD%y          ) ) DEALLOCATE( IceD%y          )      
      IF ( ALLOCATED(IceD%OtherSt_old) ) DEALLOCATE( IceD%OtherSt_old)      
                        
      ! -------------------------------------------------------------------------
      ! predicted state variables:
      ! -------------------------------------------------------------------------
                                                          
      CALL AD_DestroyContState   ( AD%x( STATE_PRED),      ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      CALL AD_DestroyDiscState   ( AD%xd(STATE_PRED),      ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL AD_DestroyConstrState ( AD%z( STATE_PRED),      ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL AD_DestroyOtherState  ( AD%OtherSt_old,         ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
                                                          
      CALL SrvD_DestroyContState   ( SrvD%x( STATE_PRED),  ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      CALL SrvD_DestroyDiscState   ( SrvD%xd(STATE_PRED),  ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL SrvD_DestroyConstrState ( SrvD%z( STATE_PRED),  ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL SrvD_DestroyOtherState  ( SrvD%OtherSt_old,     ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
                                                          
      CALL HydroDyn_DestroyContState   ( HD%x( STATE_PRED),ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      CALL HydroDyn_DestroyDiscState   ( HD%xd(STATE_PRED),ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL HydroDyn_DestroyConstrState ( HD%z( STATE_PRED),ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL HydroDyn_DestroyOtherState  ( HD%OtherSt_old,   ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
                                                        
      CALL SD_DestroyContState   ( SD%x( STATE_PRED),      ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      CALL SD_DestroyDiscState   ( SD%xd(STATE_PRED),      ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL SD_DestroyConstrState ( SD%z( STATE_PRED),      ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL SD_DestroyOtherState  ( SD%OtherSt_old,         ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
                                                          
      CALL MAP_DestroyContState   ( MAPp%x( STATE_PRED),   ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      CALL MAP_DestroyDiscState   ( MAPp%xd(STATE_PRED),   ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL MAP_DestroyConstrState ( MAPp%z( STATE_PRED),   ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL MAP_DestroyOtherState  ( MAPp%OtherSt_old,      ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  

!TODO:
!BJJ: do I have to call these other routines for MAP's c_obj stuff, like Marco indicates in his glue code?
!CALL MAP_InitInput_Destroy ( MAP_InitInput%C_obj%object )              
      
      CALL FEAM_DestroyContState   ( FEAM%x( STATE_PRED),  ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      CALL FEAM_DestroyDiscState   ( FEAM%xd(STATE_PRED),  ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL FEAM_DestroyConstrState ( FEAM%z( STATE_PRED),  ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL FEAM_DestroyOtherState  ( FEAM%OtherSt_old,     ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      
      CALL IceFloe_DestroyContState   ( IceF%x( STATE_PRED),ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      CALL IceFloe_DestroyDiscState   ( IceF%xd(STATE_PRED),ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL IceFloe_DestroyConstrState ( IceF%z( STATE_PRED),ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL IceFloe_DestroyOtherState  ( IceF%OtherSt_old,   ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      
      
      
      !............................................................................................................................
      ! Set exit error code if there was an error;
      !............................................................................................................................
      IF (Error) THEN !This assumes PRESENT(ErrID) is also .TRUE. :
         IF ( m_FAST%t_global < t_initial ) THEN
            ErrMsg = 'at initialization'
         ELSEIF ( n_t_global > m_FAST%n_TMax_m1 ) THEN
            ErrMsg = 'after computing the solution'
         ELSE            
            ErrMsg = 'at simulation time '//TRIM(Num2LStr(m_FAST%t_global))//' of '//TRIM(Num2LStr(p_FAST%TMax))//' seconds'
         END IF
                    
         
         CALL ProgAbort( 'FAST encountered an error '//TRIM(ErrMsg)//'.'//NewLine//' Simulation error level: '&
                         //TRIM(GetErrStr(ErrLev)), TrapErrors=.FALSE., TimeWait=3._ReKi )  ! wait 3 seconds (in case they double-clicked and got an error)
      END IF
      
      !............................................................................................................................
      !  Write simulation times and stop
      !............................................................................................................................

      CALL RunTimes( m_FAST%StrtTime, m_FAST%UsrTime1, m_FAST%SimStrtTime, m_FAST%UsrTime2, m_FAST%t_global, UsrTimeDiff )

      CALL NormStop( )


   END SUBROUTINE ExitThisProgram
   !...............................................................................................................................

END PROGRAM FAST
!=======================================================================
