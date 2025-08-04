!**********************************************************************************************************************************
! FAST_Solver.f90, FAST_Subs.f90, FAST_Lin.f90, and FAST_Mods.f90 make up the FAST glue code in the FAST Modularization Framework.
! FAST_Prog.f90, FAST_Library.f90, FAST_Prog.c are different drivers for this code.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
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
!**********************************************************************************************************************************
MODULE FAST_SS_Subs

   USE FAST_SS_Solver
   
   IMPLICIT NONE

   
CONTAINS
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! DRIVER ROUTINE (runs + ends simulation)
! Put here so that we can call from either stand-alone code or from the ENFAST executable.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE FAST_RunSteadyStateDriver( Turbine )
   TYPE(FAST_TurbineType),            INTENT(INOUT) :: Turbine        !< all data for one instance of a turbine
   
   INTEGER(IntKi)                                   :: ErrStat        !< Error status of the operation
   CHARACTER(ErrMsgLen)                             :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   ProgName = TRIM(FAST_Ver%Name)//' Steady State'
   FAST_Ver%Name = ProgName

   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! initialization
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         
   CALL FAST_InitializeSteadyState_T( Turbine, ErrStat, ErrMsg )
   CALL CheckError( ErrStat, ErrMsg, 'during module initialization' )
                        
      
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Calculate steady-state solutions:
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   CALL FAST_SteadyState_T( Turbine, ErrStat, ErrMsg )
   CALL CheckError( ErrStat, ErrMsg, 'during steady-state solve' )
   
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !  Clean up and stop
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   CALL ExitThisProgram_T( Turbine, ErrID_None, .true. )
   
   CONTAINS
      !...............................................................................................................................
      SUBROUTINE CheckError(ErrID,Msg,SimMsg)
      ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
      !...............................................................................................................................

            ! Passed arguments
         INTEGER(IntKi), INTENT(IN)           :: ErrID       ! The error identifier (ErrStat)
         CHARACTER(*),   INTENT(IN)           :: Msg         ! The error message (ErrMsg)
         CHARACTER(*),   INTENT(IN)           :: SimMsg      ! a message describing the location of the error

         IF ( ErrID /= ErrID_None ) THEN
            CALL WrScr( NewLine//TRIM(Msg)//NewLine )

            IF ( ErrID >= AbortErrLev ) THEN
               CALL ExitThisProgram_T( Turbine, ErrID, .true., SimMsg )
            END IF
         
         END IF

      END SUBROUTINE CheckError
END SUBROUTINE FAST_RunSteadyStateDriver
!----------------------------------------------------------------------------------------------------------------------------------

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! INITIALIZATION ROUTINES
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE FAST_InitializeSteadyState_T( Turbine, ErrStat, ErrMsg )
   TYPE(FAST_TurbineType),            INTENT(INOUT) :: Turbine        !< all data for one instance of a turbine
   INTEGER(IntKi),                    INTENT(  OUT) :: ErrStat        !< Error status of the operation
   CHARACTER(*),                      INTENT(  OUT) :: ErrMsg         !< Error message if ErrStat /= ErrID_None
   
   LOGICAL,                  PARAMETER              :: CompAeroMaps = .true.
   REAL(DbKi),               PARAMETER              :: t_initial = 0.0_DbKi

   Turbine%TurbID = 1

      CALL FAST_InitializeAll( t_initial, Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
                     Turbine%ED, Turbine%SED, Turbine%BD, Turbine%SrvD, Turbine%AD, Turbine%ADsk, Turbine%ExtLd, Turbine%IfW, Turbine%ExtInfw, &
                     Turbine%SeaSt, Turbine%HD, Turbine%SD, Turbine%ExtPtfm, Turbine%MAP, Turbine%FEAM, Turbine%MD, Turbine%Orca, &
                     Turbine%IceF, Turbine%IceD, Turbine%MeshMapData, CompAeroMaps, ErrStat, ErrMsg )

      call InitFlowField()

contains
   !> AD15 now directly accesses FlowField data from IfW.  Since we don't use IfW, we need to manually set the FlowField data
   !! NOTE: we deallocate(AD%p%FlowField) at the end of the simulation if CompAeroMaps is true
   subroutine InitFlowField()
      use InflowWind_IO, only: IfW_SteadyWind_Init
      use InflowWind_IO_Types, only: InflowWind_IO_DestroySteady_InitInputType, InflowWind_IO_DestroyWindFileDat
      type(Steady_InitInputType) :: InitInp
      integer(IntKi)             :: SumFileUnit = -1
      type(WindFileDat)          :: WFileDat          ! throw away data returned form init
      integer(IntKi)             :: ErrStat2
      character(ErrMsgLen)       :: ErrMsg2

      allocate(Turbine%AD%p%FlowField)
      Turbine%AD%p%FlowField%FieldType = 1  ! Steady wind, init below.
      InitInp%RefHt      = 100.0_ReKi  ! Any value will do here.  No exponent, so this doesn't matter
      InitInp%HWindSpeed = 8.0_ReKi    ! This gets overwritten later before used
      InitInp%PLExp      = 0.0_ReKi    ! no shear used
      call IfW_SteadyWind_Init(InitInp, SumFileUnit, Turbine%AD%p%FlowField%Uniform, WFileDat, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'FAST_InitializeSteadyState_T:InitFlowField')
      if (ErrStat >= AbortErrLev) deallocate(Turbine%AD%p%FlowField)

      call InflowWind_IO_DestroySteady_InitInputType(InitInp, ErrStat2, ErrMsg2)  ! ignore errors here because I'm lazy
      call InflowWind_IO_DestroyWindFileDat(WFileDat,  ErrStat2, ErrMsg2)  ! ignore errors here because I'm lazy
   end subroutine
END SUBROUTINE FAST_InitializeSteadyState_T
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that calls FAST_Solution for one instance of a Turbine data structure. This is a separate subroutine so that the FAST
!! driver programs do not need to change or operate on the individual module level. 
SUBROUTINE FAST_SteadyState_T( Turbine, ErrStat, ErrMsg )

   TYPE(FAST_TurbineType),   INTENT(INOUT) :: Turbine             !< all data for one instance of a turbine
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None
   
      CALL FAST_SteadyState( Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
                  Turbine%ED, Turbine%BD, Turbine%AD, Turbine%MeshMapData, ErrStat, ErrMsg )
                  
END SUBROUTINE FAST_SteadyState_T
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine takes data from n_t_global and gets values at n_t_global + 1
SUBROUTINE FAST_SteadyState(p_FAST, y_FAST, m_FAST, ED, BD, AD, MeshMapData, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         !< Data for mapping between modules
      
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None
   
   INTEGER(IntKi)                          :: n_case              !< loop counter
   REAL(DbKi)                              :: n_global
   REAL(ReKi), ALLOCATABLE                 :: UnusedAry(:)
   REAL(R8Ki), ALLOCATABLE                 :: Jmat(:,:)
   TYPE(FAST_SS_CaseType)                  :: caseData            ! tsr, windSpeed, pitch, and rotor speed for this case
   TYPE(FAST_SS_CaseType)                  :: caseData_try2       ! tsr, windSpeed, pitch, and rotor speed for this case (to try a different operating point first)
   
   INTEGER(IntKi)                          :: NStatus
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   TYPE(IceD_OutputType),    ALLOCATABLE   :: y_IceD (:)         !< IceDyn outputs (WriteOutput values are subset)
   CHARACTER(MaxWrScrLen), PARAMETER       :: BlankLine = " "
   
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_SteadyState'
   
   ErrStat = ErrID_None
   ErrMsg = ""
   
  CALL InitSSVariables(p_FAST, y_FAST, m_FAST, ED, BD, AD, MeshMapData, JMat, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if
   
      ! how often do we inform the user which case we are on?
   NStatus = min( 100, p_FAST%NumSSCases/100 + 1) ! at least 100 every 100 cases or 100 times per simulation
   call WrScr(NewLine)

   DO n_case = 1, p_FAST%NumSSCases
   
      if (mod(n_case,NStatus) == 0 .or. n_case==p_FAST%NumSSCases .or. n_case==1) then
         call WrOver( ' Case '//trim(num2lstr(n_case))//' of '//trim(num2lstr(p_FAST%NumSSCases)) )
      end if
      
   
      if (p_FAST%WindSpeedOrTSR==1) then
         caseData%windSpeed = p_FAST%WS_TSR(n_case)
         caseData%tsr       = p_FAST%RotSpeed(n_case) * AD%p%rotors(1)%BEMT%rTipFixMax / caseData%windSpeed
      else
         caseData%tsr = p_FAST%WS_TSR(n_case)
         caseData%windSpeed = p_FAST%RotSpeed(n_case) * AD%p%rotors(1)%BEMT%rTipFixMax  / caseData%tsr
      end if
      caseData%pitch = p_FAST%Pitch(n_case)
      caseData%RotSpeed = p_FAST%RotSpeed(n_case)
      
      ! Call steady-state solve for this pitch and rotor speed
      call SolveSteadyState(caseData, Jmat, p_FAST, y_FAST, m_FAST, ED, BD, AD, MeshMapData, ErrStat2, ErrMsg2)
      
      if (ErrStat2 >= ErrID_Severe) then
         ! we didn't converge; let's try a different operating point and see if that helps:
         caseData_try2%RotSpeed  = caseData%RotSpeed
         caseData_try2%Pitch     = caseData%Pitch     * 0.5_ReKi
         caseData_try2%TSR       = caseData%TSR       * 0.5_ReKi
         caseData_try2%WindSpeed = caseData%WindSpeed * 0.5_ReKi
         
         call WrScr('Retrying case '//trim(num2lstr(n_case))//', first trying to get a better initial guess. Average error is '// &
                                      trim(num2lstr(y_FAST%DriverWriteOutput(SS_Indx_Err)))//'.')
         call SolveSteadyState(caseData_try2, Jmat, p_FAST, y_FAST, m_FAST, ED, BD, AD, MeshMapData, ErrStat2, ErrMsg2)
         
         ! if that worked, try the real case again:
         if (ErrStat2 < AbortErrLev) then
            call SolveSteadyState(caseData, Jmat, p_FAST, y_FAST, m_FAST, ED, BD, AD, MeshMapData, ErrStat2, ErrMsg2)
            call WrOver(BlankLine)
         end if
         
      end if
      
      if (ErrStat2 > ErrID_None) then
         ErrMsg2 = trim(ErrMsg2)//" case "//trim(num2lstr(n_case))//&
                                    ' (tsr='//trim(num2lstr(caseData%tsr))//&
                                    ', wind speed='//trim(num2lstr(caseData%windSpeed))//' m/s'//&
                                    ', pitch='//trim(num2lstr(caseData%pitch*R2D))//' deg'//&
                                    ', rotor speed='//trim(num2lstr(caseData%RotSpeed*RPS2RPM))//' rpm)'
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      end if
         
      !----------------------------------------------------------------------------------------
      ! Write results to file
      !----------------------------------------------------------------------------------------
      n_global = real(n_case, DbKi) ! n_global is double-precision so that we can reuse existing code.
      
      CALL WrOutputLine( n_global, p_FAST, y_FAST, UnusedAry, UnusedAry, ED%y%WriteOutput, UnusedAry, &
            AD%y, UnusedAry, UnusedAry, UnusedAry, UnusedAry, UnusedAry, UnusedAry, UnusedAry, &
            UnusedAry, UnusedAry, UnusedAry, UnusedAry, y_IceD, BD%y, ErrStat2, ErrMsg2 )

         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if
         
      
      ! in case we have a lot of error messages, let's print the non fatal ones here:
      if (ErrStat > ErrID_None) then
         call WrScr(trim(ErrMsg))
         call WrScr("")
         ErrStat = ErrID_None
         ErrMsg = ""
      end if

   END DO
   
CONTAINS
   SUBROUTINE Cleanup()
      if (allocated(Jmat)) deallocate(Jmat)
   END SUBROUTINE Cleanup
   
   
END SUBROUTINE FAST_SteadyState
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE InitSSVariables(p_FAST, y_FAST, m_FAST, ED, BD, AD, MeshMapData, JMat, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         !< Data for mapping between modules
   REAL(R8Ki), ALLOCATABLE , INTENT(INOUT) :: Jmat(:,:)           !< Matrix for storing Jacobian

   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None
   
   INTEGER(IntKi)                          :: NumBlades           !< number of blades

   
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   
   CHARACTER(*), PARAMETER                 :: RoutineName = 'SS_InitVariables'
   
   ErrStat = ErrID_None
   ErrMsg = ""

   NumBlades = size(AD%y%rotors(1)%BladeLoad)
   
   
   call AllocAry(Jmat, p_FAST%SizeJac_Opt1(1), p_FAST%SizeJac_Opt1(1), 'Jmat', ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
   CALL AllocAry( MeshMapData%Jacobian_pivot, p_FAST%SizeJac_Opt1(1), 'Pivot array for Jacobian LU decomposition', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   !CALL AllocAry( MeshMapData%HubOrient, 3, 3, NumBlades, 'Hub orientation matrix', ErrStat2, ErrMsg2 )
   !   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   if (ErrStat >= AbortErrLev) return
      
      
   CALL CopyStatesInputs( p_FAST, ED, BD, AD, ErrStat2, ErrMsg2, MESH_NEWCOPY )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )


END SUBROUTINE InitSSVariables
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE FAST_SS_Subs
!----------------------------------------------------------------------------------------------------------------------------------
